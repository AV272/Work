#include <TLatex.h>
#include <TMath.h>
#include "TH2F.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void difcs() 
{
   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   //c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);

   // Read DATA file
   ifstream myfile ("Th66_Data_dsigm[theta_cm]_DdnHe.txt");

   // Put DATA in new file
   //ofstream myfile2;

   // Create MultiGraph structure in order to put fits of all energies
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Angular distribution in ^{2}H(d,n)^{3}He");

   TLatex latex;

   // Introduce data and energies from article Th66
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

   //Introduce arrays for data from picture (take from file)
   const Int_t num = 15;
   Double_t theta_picture[num];
   Double_t dcsn_picture[num];
   Double_t theta_error_picture[num];
   Double_t dcsn_error_picture[num];
   string line;
   string s;
   Double_t buffer;
   
   // Introduce arrays for angles (one step for each angle degree)
   const Int_t nn=180;
   Double_t dcsn[nn];
   Double_t dcsn2[nn];
   Double_t dcsp[nn];
   Double_t theta[nn];
   Double_t edcsn[nn];
   Double_t edcsn2[nn];
   Double_t edcsp[nn];
   Double_t theta_cm[nn], theta_cm2[nn];
   Double_t dth = (3.1415926*0.5)/180;

   const Int_t nnn =nn;
   Double_t dcsn_cm_all[nnn], theta_cm_all[nnn], edcsn_cm_all[nnn];

   TString g_name;
   TString myfile2_name;



   for (int i = 0; i < n; i++)
   {  

      for (int j = 0; j < nn; j++)
      {  
         Double_t theta0 = j;
         theta[j] = theta0;
         theta0 = (3.1415926*theta0)/180;

 ///<<<///////////////////Translation angles to the center mass system//////////////////////////
         // Input parameters: m1, m2, m3, m4 -- particles masses; Ed[] -- kinetic energy of deutron
         // in lab system; theta[] -- laboratory scattering angle of particle 3. 
         Double_t th, m1, m2, m3, m4, E1, E_cmT, p1, check, E_cm4, alpha, E3, E3_2, ET,\
         p_cm, p3, sin_cm, p3_2, sin_cm_2, th_2, cos_cm;
         Double_t th_cm =0, th_cm2 =0;

         // Introduce rest mass for reaction particles        
         m1 = m2 = 1875.61*pow(10,3); // deuterium MeV
         m3 = 2808.92*pow(10,3); // He3
         m4 = 939.565*pow(10,3); // n
         th = theta0;

         E1 = Ed[i] + m1;
         ET = E1 + m2;
         p1 = sqrt(E1*E1 - m1*m1);
         E_cmT = sqrt(m1*m1 + m2*m2 + 2*m2*E1);
         //Check to see if particle enough energy for this angle
         check = pow(m2*E1 + (m1*m1 + m2*m2 - m3*m3 -m4*m4)/2,2) - pow(m3*m4,2) - pow(p1*m3*sin(th),2);
         E_cm4 = (E_cmT*E_cmT + m4*m4 - m3*m3)/(2*E_cmT);
         p_cm = sqrt(E_cm4*E_cm4 - m4*m4 );

         // alpha determines if there two or one solutions for E3 (in this case only one solution)
         alpha = (p1*(1+(m3*m3 - m4*m4)/(E_cmT*E_cmT)))/\
         (ET*sqrt((1-pow((m3 + m4)/E_cmT,2))*(1-pow((m3-m4)/E_cmT,2))));
         if (alpha>1){
            E3 = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2) + p1*cos(th)*sqrt(check))/\
            (ET*ET - pow(p1*cos(th),2));
            E3_2 = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2) - p1*cos(th)*sqrt(check))/\
            (ET*ET - pow(p1*cos(th),2));
            p3 = sqrt(E3*E3 - m3*m3);
            sin_cm = p3*sin(th)/p_cm;

            p3_2 = sqrt(E3_2*E3_2 - m3*m3);
            sin_cm_2 = p3_2*sin(th)/p_cm;

            th_cm = TMath::ASin(sin_cm);
            printf("Theta_cm_%f\n",th_cm);
            th_cm2 = TMath::ASin(sin_cm_2);
            printf("Theta_cm_222%f\n",th_cm2);
            theta_cm[j] = 180*th_cm/3.1415926;
            theta_cm2[j] = 180*th_cm2/3.1415926;

            dcsn[j]  = 1 + An[i]*pow(cos(th_cm),2) + Bn[i]*pow(cos(th_cm),4)+0.1*i - 0.2;
            edcsn[j] = pow(cos(th_cm),2)*eAn[i] + pow(cos(th_cm),4)*eBn[i] +\
            An[i]*2*cos(th_cm)*sin(th_cm)*dth + Bn[i]*4*pow(cos(th_cm),3)*sin(th_cm)*dth;

            dcsn2[j]  = 1 + An[i]*pow(cos(th_cm2),2) + Bn[i]*pow(cos(th_cm2),4)+0.1*i - 0.2;
            edcsn2[j] = pow(cos(th_cm2),2)*eAn[i] + pow(cos(th_cm2),4)*eBn[i] +\
            An[i]*2*cos(th_cm2)*sin(th_cm2)*dth + Bn[i]*4*pow(cos(th_cm2),3)*sin(th_cm2)*dth;

         } else {
            E3 = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2) + p1*cos(th)*sqrt(check))/\
            (ET*ET - pow(p1*cos(th),2));
            p3 = sqrt(E3*E3 - m3*m3);
            cos_cm = (ET*(p3*cos(th) - (p1*E3/ET)))/(E_cmT*p_cm);
            //sin_cm = p3*sin(th)/p_cm;
            th_cm = TMath::ACos(cos_cm);
            theta_cm[j] = 180*th_cm/3.1415926;

            dcsn[j]  = 1 + An[i]*pow(cos(th_cm),2) + Bn[i]*pow(cos(th_cm),4)+0.1*(Ed[i]/25.) - 0.1;
            dcsp[j]  = 1 + Ap[i]*pow(cos(th),2) + Bp[i]*pow(cos(th),4)+0.1*i;
            edcsn[j] = pow(cos(th_cm),2)*eAn[i] + pow(cos(th_cm),4)*eBn[i] +\
            An[i]*2*cos(th_cm)*sin(th_cm)*dth + Bn[i]*4*pow(cos(th_cm),3)*sin(th_cm)*dth;
            edcsp[j] = pow(cos(th),2)*eAp[i] + pow(cos(th),4)*eBp[i] +\
            Ap[i]*2*cos(th)*sin(th)*dth + Bp[i]*4*pow(cos(th),3)*sin(th)*dth;
         }
         
//>>>//////////////////////////////////////////////////////////////////////////////////
      }

      // Delete energy 19.6 and 32.0 because they abscent in article picture
      if ((i==0)||(i==2)){continue;}

      for (int k =0; k< nn; k++){
            dcsn_cm_all[k] = dcsn[k];
            edcsn_cm_all[k] = edcsn[k];
            theta_cm_all[k] = theta_cm[k];
            //dcsn_cm_all[k+nn] = dcsn2[k];
            //edcsn_cm_all[k+nn] = edcsn2[k];
            //theta_cm_all[k+nn] = theta_cm2[k];
         
      }

//<<<//////////////////////////GET DATA FROM FILE////////////////////////////////////////////
      if (myfile.is_open())
      {
         myfile >>buffer;
         //printf("Ed_%f\n",Ed[i]);
         //printf("buffer_%f\n",buffer);
         if (Ed[i]==buffer){
            for (int i2 =0; i2<num; i2++){
               // Read theta from file
               myfile >> theta_picture[i2];
               //cout << theta_picture[i2] << endl;
               // Theta error
               theta_error_picture[i2] = 0.5;
               Double_t th_cm = theta_picture[i2];
               // Cross section error (quadratic error)
               dcsn_error_picture[i2] = sqrt(pow(pow(cos(th_cm),2)*eAn[i],2) +\
               pow(pow(cos(th_cm),2)*eAn[i],2) + pow(An[i]*2*cos(th_cm)*sin(th_cm)*dth,2)+\
               pow(Bn[i]*4*pow(cos(th_cm),3)*sin(th_cm)*dth,2));
            }
            for (int i3 =0; i3<num; i3++){
               // Read CS from file
               myfile >> dcsn_picture[i3];
               //cout << dcsn_picture[i3] << endl; 
               dcsn_picture[i3] = dcsn_picture[i3] +  0.1*(Ed[i]/25.) - 0.1;                
            }
         }else{
            // if energy not equal then jump on three lines
            getline(myfile,s);
            getline(myfile,s);
            getline(myfile,s);
         }
         
      }
      else {cout << "Unable to open file";} 

      //////////////////PUT DATA IN FILE///////////////////////////////////

      //myfile2_name.Form("difsig_%f_n_unpol_Th66.txt",Ed[i]);
      //myfile2.open(myfile2_name);
      //for(int i4=0; i4<num; i4++){
      //   myfile2 << theta_picture[i4]<< " "<< theta_error_picture[i4] <<\
      //   " "<< dcsn_picture[i4] - 0.1*(Ed[i]/25.) + 0.1<< " "<< dcsn_error_picture[i4] << "\n";
      //} 
      //myfile2.close();
      

//>>>///////////////////////////////////////////////////////////////////////

      auto dcsng = new TGraphErrors(nnn,theta_cm_all,dcsn_cm_all,0,0);
      dcsng->SetMarkerColor(1);
      dcsng->SetMarkerStyle(1);
      if (i == 9){dcsng->SetLineColor(3);}
      else{dcsng->SetLineColor(i+1);}      
      g_name.Form("Energy_%f",Ed[i]);
      dcsng->SetTitle(g_name);

      auto dcsng_picture = new TGraphErrors(num,theta_picture,dcsn_picture,\
      theta_error_picture,dcsn_error_picture);
      dcsng_picture->SetMarkerColor(1);
      dcsng_picture->SetMarkerStyle(2);
      if (i == 9){dcsng_picture->SetLineColor(3);}
      else{dcsng_picture->SetLineColor(i+1);}
      
      // Adding graphics
      mg->Add(dcsng, "CP");
      mg->Add(dcsng_picture, "AP");

      mg->Draw("A");
      mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
      mg->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
      gPad->Modified();
      //mg->GetXaxis()->SetLimits(1.5,7.5);
      mg->SetMinimum(0.8);
      mg->SetMaximum(4.5);
   }
   c1->BuildLegend(0.35,0.55,0.65,0.89);
   myfile.close();
}
