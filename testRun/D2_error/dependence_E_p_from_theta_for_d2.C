// RELATIVISTIC KINEMATICS (ROOT program)
// Program for relativistic kinematic calculations (elastic and inelastic)
// Programmed by R.E. Philips and S.T. Thornton
// Oak Ridge national lab, july 1967.
// Reload by Andreyanov A.

///////////////////////////////////////////////////////////////////////
// Reaction 1 + 2 --> 3 + 4

// c = 1

// INPUT DATA
// 1) m1, m2, m3, m4 -- particles masses of rest (Mev); 
// 2) E1_k -- kinetic energy of particle 1 in laboratory system;
// 3) theta -- laboratory scattering angle of particle 3.

// OUTPUT DATA
// 1) E1_cm; E2_cm; E3_cm; E4_cm -- energy in center of mass system 
// 2) Theta3_cm; Theta4_cm -- angles of outgoing particles in CM system 
// 3) Theta4 -- angle of particle 4 in lab system
// 3) E3; E4 -- energies of outgoing particles in lab system

#include <TMath.h>
#include <TLatex.h>
using namespace std;

void dependence_E_p_from_theta_for_d2(){
    double m1, m2, m3, m4, Q, E1_k, theta, EN, Q_gs, T, E1, E_th, ET, p1,\
    T_A, sin_th3, cos_th3, check, E_cmT, E1_cm, E2_cm, E3_cm,  E4_cm, \
    p_cm, alpha, beta, gamma, E1_k_min;
    //double E3, E4, p4, p3, sin_th3_cm, cos_th3_cm, th3_cm, th4_cm,\
    //sin_th4, cos_th4, theta4, dsigma_cm__dsigma_lab;
    TString g_name; 

    // Input data (Mev and degrees)
    double m1_arr[2] = {1875.61, 1875.61}; // D2
    m2 = 1875.61;
    double m4_arr[2] = {2808.92, 2808.39};
    double m3_arr[2] = {938.272, 939.565};
    double E1_k_arr[2] = {0.025, 0.05};

    int n3 = 180, n3_2 = 360, n4=180;

    double E_d_arr[n3_2], E_p_arr[n3], theta_p_arr[n3], dE[n4], dTheta[n4], E_T_arr[n3], theta_T_arr[n3];

    auto c1 = new TCanvas("c1","Multigraph",700,500);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    //c1->Divide(2);

    auto mg = new TMultiGraph();
    auto mg2 = new TMultiGraph();


    for (int ii=0; ii<2; ii++){
        E1_k = E1_k_arr[ii];
        m1 = m1_arr[0];
        
        for (int i4=0; i4<2; i4++){
            m4 = m4_arr[i4];
            m3 = m3_arr[i4];

            for (int num=0; num<n3; num++){
                theta = num;
                theta_p_arr[num] = theta;

                Q_gs = m1 + m2 - m3 - m4; // Q-value groundstate
                E1_k_min = (-Q_gs*(m1+m2+m3+m4)/(2*m2)); // Threshold energy
                // Check it is enough energy for this reaction
                if (E1_k <= E1_k_min){
                    cout << "Not enough kinetic energy for this reaction" << "\n";
                    E_p_arr[num] = 0;
                    return;
                }
                E1 = E1_k + m1; // Full energy of particle 1
                E_th = E1_k_min + m1; // Full threshold energy
                ET = E1 + m2; // Total energy (ET = E1+E2 = E3+E4)
                p1 = sqrt(E1*E1 - m1*m1); //  Momentum of particle 1 in lab
                
                //Change angles to radians
                theta = theta*0.017453292;
                //printf("theta  %f\n", theta);
                sin_th3 = sin(theta);
                cos_th3 = cos(theta);

                // Check it is enough energy for this angle
                check = pow(m2*E1 + (m1*m1 + m2*m2 - m3*m3 -m4*m4)/2,2) \
                - pow(m3*m4,2) - pow(p1*m3*sin_th3,2);
                //printf("check %f\n", check);
                if (check < 0){
                    cout << "Not enough energy for this angle " << "\n";
                    E_p_arr[num] = 0;
                    return;
                }

                // Calculate CENTER OF MASS energies
                E_cmT = sqrt(m1*m1 + m2*m2 + 2*m2*E1);
                E1_cm = (m1*m1 + m2*E1)/E_cmT;
                E2_cm = (m2*m2 + m2*E1)/E_cmT;
                E3_cm = (E_cmT*E_cmT + m3*m3 - m4*m4)/(2*E_cmT);
                E4_cm = (E_cmT*E_cmT + m4*m4 - m3*m3)/(2*E_cmT);

                p_cm = sqrt(E4_cm*E4_cm - m4*m4); // Center of mass momentum 
                                                    // of outgoing particles 

                // Alpha determines if there two or one solutions for E3
                alpha = (p1*(1+(m3*m3 - m4*m4)/(E_cmT*E_cmT)))/\
                (ET*sqrt((1-pow((m3 + m4)/E_cmT,2))*(1-pow((m3-m4)/E_cmT,2))));
                //printf("alpha %f\n", alpha);
                int n;
                if (alpha<=0){
                    cout << "One root for E3" << "\n";
                    n = 1;
                } else{
                    cout << "Two roots for E3" << "\n";
                    n = 2;
                }

                double E4[n], E3[n], p3[n], p4[n], sin_th3_cm[n], cos_th3_cm[n], th3_cm[n], th4_cm[n],\
                sin_th4[n], cos_th4[n], theta4[n], dsigma_cm__dsigma_lab[n]; 
                beta = p1/ET;
                //printf("beta %f\n", beta);
                gamma = ET/E_cmT;
                //printf("gamma %f\n", gamma);
                // Make loop for finding roots: one step for alpha < 0 and 
                // two steps for alpha > 0
                for (int i=0; i<n; i++){
                    E3[i] = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2)\
                    + (1 -2*i)*p1*cos_th3*sqrt(check))/(ET*ET - pow(p1*cos_th3,2));
                    E4[i] = ET - E3[i];
                    p3[i] = sqrt(E3[i]*E3[i] - m3*m3);
                    p4[i] = sqrt(E4[i]*E4[i] - m4*m4);
                    // Calculate center of mass angles of particles 3 and 4
                    sin_th3_cm[i] = p3[i]*sin_th3/p_cm;
                    cos_th3_cm[i] = gamma*(p3[i]*cos_th3 - beta*E3[i])/p_cm;
                    th3_cm[i] = TMath::ATan(sin_th3_cm[i]/cos_th3_cm[i]); // Theta3_cm
                    if (cos_th3_cm[i] < 0){
                        th3_cm[i] = 3.14159265 + th3_cm[i];
                    }
                    th4_cm[i] = 3.14159265 - th3_cm[i]; // In CM system theta3_cm + theta4_cm =180
                    // Calculate lab scattering angle of particle 4
                    sin_th4[i] = p3[i]*sin_th3/p4[i];
                    cos_th4[i] = (p1 - p3[i]*cos_th3)/p4[i];
                    theta4[i] = TMath::ATan(sin_th4[i]/cos_th4[i]);
                    if (cos_th4[i] < 0){
                        theta4[i] = 3.14159265 + theta4[i];
                    }

                    // Ratio differential cross section in CM to LAB
                    dsigma_cm__dsigma_lab[i] = ET/E_cmT*(1 + (1 -2*i)*p1/ET*(ET*ET - p1*p1 + m3*m3 -\
                    m4*m4)/sqrt(pow(ET*ET - p1*p1 + m3*m3 - m4*m4, 2) - 4*E_cmT*E_cmT*m3*m3)*\
                    cos_th3_cm[i])*pow(sin_th3,3)/pow(sin_th3_cm[i],3);

                    // Translate angles to degrees
                    th4_cm[i] = th4_cm[i]/0.017453292;
                    th3_cm[i] = th3_cm[i]/0.017453292;
                    theta4[i] = theta4[i]/0.017453292;
                    

                    // Change energies to kinetic (MeV)
                    E3[i] = E3[i] - m3;
                    E4[i] = E4[i] - m4;
                    
                }
                theta = theta/0.017453292;

                E_p_arr[num] = E3[0];
                E_T_arr[num] = E4[0];
                theta_T_arr[num] = theta4[0];

                //printf("E3_1 %f, E3_2 %f\n", E3[0], E3[1]);
                //printf("E4_1 %f, E4_2 %f\n", E4[0], E4[1]);

                // for(int i2=0; i2<n; i2++){
                //     if (i2!=0){
                //         E_p_arr[num+n3] = E3[i2];
                //         E_T_arr[num+n3] = E4[i2];
                //         theta_p_arr[num+n3] = theta;
                //         theta_T_arr[num+n3] = theta4[i2];
                //     } 
                //     else{
                //         E_p_arr[num] = E3[i2];
                //         E_T_arr[num] = E4[i2];
                //         theta_T_arr[num] = theta4[i2];
                //     }
                //     //printf("SOLUTION %d\n", i2+1);
                //     //printf("E3 %f\n", E3[i2]);

                // }

            }

            if (ii==0){

                if(i4==0){
                    TGraph* gr = new TGraph(n3, theta_p_arr, E_p_arr);
                    gr->GetXaxis()->SetTitle("#theta(degres)");
                    gr->GetYaxis()->SetTitle("E_{p kinetic}(Mev)");
                    gr->SetMarkerColor(1);
                    gr->SetMarkerStyle(20);
                    gr->SetLineColor(1);
                    g_name.Form("proton_%f",E1_k);
                    gr->SetTitle(g_name);
                    mg->Add(gr);
                    //gr->SetMinimum(1.9);
                    //gr->SetMaximum(2.23);

                }
                
                
                TGraph* gr2 = new TGraph(n3, theta_T_arr, E_T_arr);
                gr2->GetXaxis()->SetTitle("#theta(degres)");
                gr2->GetYaxis()->SetTitle("#E_{T kinetic}(Mev)");
                
                if(i4==0){
                    g_name.Form("Triton_%f",E1_k);
                    gr2->SetMarkerColor(2);
                    gr2->SetMarkerStyle(21);
                    gr2->SetLineColor(2);
                }else{
                    g_name.Form("Helium_%f",E1_k);
                    gr2->SetMarkerColor(3);
                    gr2->SetMarkerStyle(22);
                    gr2->SetLineColor(3);
                }
                
                gr2->SetTitle(g_name);
                mg->Add(gr2);

            }else{
                if(i4==0){
                    TGraph* gr3 = new TGraph(n3, theta_p_arr, E_p_arr);
                    gr3->GetXaxis()->SetTitle("#theta(degres)");
                    gr3->GetYaxis()->SetTitle("E_{p kinetic}(Mev)");
                    gr3->SetMarkerColor(4);
                    gr3->SetMarkerStyle(23);
                    gr3->SetLineColor(4);
                    g_name.Form("proton_%f",E1_k);
                    gr3->SetTitle(g_name);
                    mg->Add(gr3);
                    //gr->SetMinimum(1.9);
                    //gr->SetMaximum(2.23);
                }
                
                
                TGraph* gr4 = new TGraph(n3, theta_T_arr, E_T_arr);
                gr4->GetXaxis()->SetTitle("#theta(degres)");
                gr4->GetYaxis()->SetTitle("#E_{T kinetic}(Mev)");
                
                if(i4==0){
                    g_name.Form("Triton_%f",E1_k);
                    gr4->SetMarkerColor(7);
                    gr4->SetMarkerStyle(33);
                    gr4->SetLineColor(7);
                }else{
                    g_name.Form("Helium_%f",E1_k);
                    gr4->SetMarkerColor(9);
                    gr4->SetMarkerStyle(34);
                    gr4->SetLineColor(9);
                }
                
                gr4->SetTitle(g_name);
                mg->Add(gr4);

            }

        }

        

        
        
    }

    //c1->cd(1);
    mg->SetTitle("Dependence of E_{kin} p, 3H and 3He from angle (#theta) in lab system two energies = 25 kev and 50 kev");
    mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
    mg->GetYaxis()->SetTitle("E_{kinetic}(Mev)");
    gPad->Modified();
    mg->Draw("ap");
    c1->BuildLegend(0.35,0.55,0.65,0.89);

    //printf("AA %d", mg->GetListOfGraphs()->GetEntries());
    //c1->BuildLegend(0.35,0.55,0.65,0.89);
    // c1->cd(2);
    // mg2->SetTitle("Dependence of E_{kin} from angle (#theta) in lab system for D nucleus");
    // mg2->GetXaxis()->SetTitle("#it{#theta} (degres)");
    // mg2->GetYaxis()->SetTitle("#E_{kinetic}(Mev)");
    // gPad->Modified();
    // mg2->Draw("ap");
    // c1->cd(2)->BuildLegend(0.35,0.55,0.65,0.89);

    //printf("AA %d", mg->GetListOfGraphs()->GetEntries());
    
    

    // // Change energies to kinetic (MeV)
    // E1_cm = E1_cm - m1;
    // E2_cm = E2_cm - m2;
    // E3_cm = E3_cm - m3;
    // E4_cm = E4_cm - m4;
    // theta = theta/0.017453292;

    // // PRINT OUTPUT
    // cout << "OUTPUT \nINPUT DATA" << "\n";
    // printf("m1 %f; m2 %f; m3 %f; m4 %f; E_kin %f; theta %f \n", m1, m2, m3, m4, E1_k, theta);
    // cout << "OUTPUT DATA" << "\n";
    // printf("E1_cm %f\n", E1_cm);
    // printf("E2_cm %f\n", E2_cm);
    // printf("E3_cm %f\n", E3_cm);
    // printf("E4_cm %f\n", E4_cm);

    // for(int i2=0; i2<n; i2++){
        
    //     printf("SOLUTION %d\n", i2+1);
    //     printf("Theta3_cm %f\n", th3_cm[i2]);
    //     printf("Theta4_cm %f\n", th4_cm[i2]);
    //     printf("Theta4 %f\n", theta4[i2]);

    //     printf("E3 %f\n", E3[i2]);
    //     printf("E4 %f\n", E4[i2]);
    // }
   
}
