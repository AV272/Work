#include <TLatex.h>
#include <TMath.h>

void difcs() 
{
   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   //c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);

   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Angular distribution in ^{2}H(d,n)^{3}He");

   TLatex latex;

   const Int_t n = 12;
   Double_t Ed[n]  = {19.5, 26.9, 32.0, 45.1, 71.0, 96.6, 122.0, 147.5, 194.0, 248.3, 298.5, 348.7};
   Double_t An[n]  = {0.26, 0.51, 0.46, 0.57, 0.79, 0.87, 0.97, 1.09, 1.11, 1.32, 1.23, 1.29};
   Double_t eAn[n] = {0.22, 0.12, 0.06, 0.03, 0.03, 0.02, 0.02, 0.02, 0.05, 0.04, 0.04, 0.04};
   Double_t Bn[n]  = {0.0, 0.0, 0.0, 0.04, 0.01, 0.12, 0.17, 0.22, 0.41, 0.43, 0.75, 0.76};
   Double_t eBn[n] = {0.0, 0.0, 0.0, 0.03, 0.03, 0.03 ,0.02, 0.02, 0.06, 0.05, 0.05, 0.05};
   Double_t Ap[n]  = {0.20, 0.27, 0.20, 0.34, 0.45, 0.53, 0.58, 0.65, 0.67, 0.76, 0.80, 0.65};
   Double_t eAp[n] = {0.15, 0.08, 0.04, 0.02, 0.02, 0.01, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03};
   Double_t Bp[n]  = {0.0, 0.0, 0.0, -0.03, 0.01, 0.03, 0.06, 0.11, 0.24, 0.32, 0.41, 0.66};
   Double_t eBp[n] = {0.0, 0.0, 0.0, 0.02, 0.02, 0.02 ,0.02, 0.02, 0.04, 0.03, 0.03, 0.03};
   //auto dn = new TGraphErrors(n,Ed,An,0,eAn);
   //auto dp = new TGraphErrors(n,Ed,Ap,0,eAp);
   //dn->SetTitle("Energy-dependent assymetry coefficient A");
   //dn->SetMarkerColor(4);
   //dn->SetMarkerStyle(21);
   //mg->Add(dn);
   //mg->Add(dp);
   //mg->Draw();
   //dn->Draw();

   const Int_t nn=15;
   Double_t dcsn[nn];
   Double_t dcsp[nn];
   Double_t theta[nn];
   Double_t edcsn[nn];
   Double_t edcsp[nn];
   Double_t theta_cm[nn];
   Double_t dth = (3.1415926*0.5)/180;

   TString g_name;



   for (int i = 0; i < n; i++)
   {  
      
      for (int j = 0; j < nn; j++)
      {  
         
         
         Double_t theta0 = j*10 +20;
         printf("Theta00_%f\n", theta0);
         theta[j] = theta0;
         theta0 = (3.1415926*theta0)/180;
         
         // Translation angles to the center mass system
         // Input parameters: m1, m2, m3, m4 -- particles masses; Ed[] -- kinetic energy of deutron
         // in lab system; theta[] -- laboratory scattering angle of particle 3. 
         Double_t th_cm, cos_th_cm, beta_cm, gamma_cm, E_cm, p_cm, m1, m2, m3, m4, chi, s;
         
         m_1 = m_2 = 1875.61*pow(10,3); // deuterium MeV
         m_3 = 2808.92*pow(10,3); // He3
         m_4 = 939.565*pow(10,3); // n
         



         E1 = Ed[i] + m1;
         E_cmT = sqrt(m1*m1 + m2*m2 + 2*m2*E1);

         // alfa determines if there two or one solutions for E3
         alfa = (p1*(1+(m3*m3 - m4*m4)/(E_cmT*E_cmT)))/\
         (ET*sqrt((1-pow((m3 + m4)/E_cmT,2))(1-pow((m3-m4)/E_cmT,2))));
         if (alfa>1){

         } else {

         }
         E3 = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2))
         
         E_cm4 = (E_cmT*E_cmT + m4*m4 - m3*m3)/(2*E_cmT);
         p_cm = sqrt(E_cm4*E_cm4 - m4*m4 );
         p3 = sqrt(E3*E3 - m3*m3);
         sin_cm = p3*sin_3/p_cm;
         ////////////////////////////////////////////////////////

         //Double_t th = (3.1415926*th_cm)/180;
         Double_t th = th_cm;
         dcsn[j]  = 1 + An[i]*pow(cos(th),2) + Bn[i]*pow(cos(th),4)+0.1*i - 0.2;
         dcsp[j]  = 1 + Ap[i]*pow(cos(th),2) + Bp[i]*pow(cos(th),4)+0.1*i;
         edcsn[j] = pow(cos(th),2)*eAn[i] + pow(cos(th),4)*eBn[i] +\
         An[i]*2*cos(th)*sin(th)*dth + Bn[i]*4*pow(cos(th),3)*sin(th)*dth;
         edcsp[j] = pow(cos(th),2)*eAp[i] + pow(cos(th),4)*eBp[i] +\
         Ap[i]*2*cos(th)*sin(th)*dth + Bp[i]*4*pow(cos(th),3)*sin(th)*dth;
         
      }

      if ((i==0)||(i==2)){continue;}
      auto dcsng = new TGraphErrors(nn,theta_cm,dcsn,0,edcsn);
      //dcsng->SetTitle("Differential cross-section / angle");
      dcsng->SetMarkerColor(4);
      dcsng->SetMarkerStyle(21);
      dcsng->SetLineColor(i+1);
      g_name.Form("Energy_%f",Ed[i]);
      dcsng->SetTitle(g_name);
      
      mg->Add(dcsng);
      mg->Draw("ACP");
      mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
      mg->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
      gPad->Modified();
      //mg->GetXaxis()->SetLimits(1.5,7.5);
      mg->SetMinimum(0.8);
      mg->SetMaximum(4.);
      
      //if (i==4)
      //   {
      //      dcsng->Draw();;
      //   }
      

   }
   c1->BuildLegend(0.35,0.55,0.65,0.89);
   
   //mg->Draw();
   //auto c2 = new TCanvas("c2","Differential cross-section / angle",200,10,700,500);
   //c2->SetFillColor(42);
   //c2->SetGrid();
   //c2->GetFrame()->SetFillColor(21);
   //c2->GetFrame()->SetBorderSize(12);
   //auto d_th = new TGraphErrors(n,theta,dcsn,0,edcsn);
   //dn->SetTitle("Differential cross-section / angle");
   //dn->SetMarkerColor(4);
   //dn->SetMarkerStyle(21);
   //d_th->Draw();


}
