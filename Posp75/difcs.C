#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void difcs() 
{
   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Divide(2,1);

   // auto c2 = new TCanvas("c2","A Simple Graph with error bars",200,10,700,500);
   // c2->SetGrid();
   // c2->GetFrame()->SetFillColor(21);
   // c2->GetFrame()->SetBorderSize(12);

   // Read DATA file
   //ifstream myfile ("Th66_Data_dsigm[theta_cm]_DdnHe.txt");

   // Put DATA in new file
   //ofstream myfile2;

   // Create MultiGraph structure in order to put fits of all energies
   TMultiGraph *mg1 = new TMultiGraph();
   mg1->SetTitle("Angular distribution in ^{2}H(d,n)^{3}He");

   TMultiGraph *mg2 = new TMultiGraph();
   mg2->SetTitle("Angular distribution in ^{2}H(d,p)^{3}H");

   TLatex latex;

   // Introduce data and energies from article Posp75
   const Int_t n = 18;
   Double_t Ed[n]  = {70, 75, 79, 84, 85, 89, 94, 94, 103, 104, 108, 110,\
   110, 114, 120, 124, 134, 144};
   Double_t B1_B0_p[n]  = {0.263, 0.252, 0.3, 0.344, 0.304, 0.272, 0.29, \
   0.33, 0.34, 0.299, 0.333, 0.367, 0.349, 0.328, 0.38, 0.366, 0.412, 0.401};
   Double_t error_B1_B0_p[n] = {0.032, 0.013, 0.011, 0.022, 0.02, 0.01, \
   0.018, 0.028, 0.015, 0.014, 0.02, 0.032, 0.023, 0.014, 0.016, 0.012, 0.016, 0.021};
   Double_t B1_B0_n[n]  = {0.442, 0.392, 0.43, 0.372, 0.368, 0.514, 0.447, \
   0.382, 0.532, 0.594, 0.536, 0.527, 0.558, 0.577, 0.514, 0.561, 0.605, 0.582};
   Double_t error_B1_B0_n[n] = {0.032, 0.035, 0.013, 0.01, 0.044, 0.024, \
   0.033, 0.032, 0.043, 0.03, 0.023, 0.046, 0.033, 0.036, 0.018, 0.031, 0.014, 0.045};
   
   
   //Introduce arrays for data from picture (take from file)
   // const Int_t num = 15;
   // Double_t theta_picture[num];
   // Double_t dcsn_picture[num];
   // Double_t theta_error_picture[num];
   // Double_t dcsn_error_picture[num];
   // string line;
   // string s;
   // Double_t buffer;
   
   // Introduce arrays for angles (one step for each angle degree)
   const Int_t nn=180;
   Double_t dcsn[nn];
   Double_t dcsp[nn];
   Double_t theta[nn];
   Double_t edcsn[nn];
   Double_t edcsp[nn];
   Double_t dth = (3.1415926*0.5)/180;

   TString g_name;
   TString myfile2_name;

   for (int i = 0; i < n; i++)
   {  
      for (int j = 0; j < nn; j++)
      {  
         Double_t theta0 = j;
         theta[j] = theta0;
         theta0 = (3.1415926*theta0)/180;
         Double_t pol_Legandre = 0.5*(3*pow(cos(theta0),2)-1);
         dcsn[j]  = 1 + B1_B0_n[i]*pol_Legandre;
         edcsn[j] = sqrt(pow(pol_Legandre*error_B1_B0_n[i],2) +\
         pow(B1_B0_n[i]*3*cos(theta0)*sin(theta0)*dth,2));

         dcsp[j]  = 1 + B1_B0_p[i]*pol_Legandre;
         edcsp[j] = sqrt(pow(pol_Legandre*error_B1_B0_p[i],2) +\
         pow(B1_B0_p[i]*3*cos(theta0)*sin(theta0)*dth,2));
         
      }

//<<<//////////////////////////GET DATA FROM FILE////////////////////////////////////////////
      // if (myfile.is_open())
      // {
      //    myfile >>buffer;
      //    //printf("Ed_%f\n",Ed[i]);
      //    //printf("buffer_%f\n",buffer);
      //    if (Ed[i]==buffer){
      //       for (int i2 =0; i2<num; i2++){
      //          // Read theta from file
      //          myfile >> theta_picture[i2];
      //          //cout << theta_picture[i2] << endl;
      //          // Theta error
      //          theta_error_picture[i2] = 0.5;
      //          Double_t th_cm = theta_picture[i2];
      //          // Cross section error (quadratic error)
      //          dcsn_error_picture[i2] = sqrt(pow(pow(cos(th_cm),2)*eAn[i],2) +\
      //          pow(pow(cos(th_cm),2)*eAn[i],2) + pow(An[i]*2*cos(th_cm)*sin(th_cm)*dth,2)+\
      //          pow(Bn[i]*4*pow(cos(th_cm),3)*sin(th_cm)*dth,2));
      //       }
      //       for (int i3 =0; i3<num; i3++){
      //          // Read CS from file
      //          myfile >> dcsn_picture[i3];
      //          //cout << dcsn_picture[i3] << endl; 
      //          dcsn_picture[i3] = dcsn_picture[i3] +  0.1*(Ed[i]/25.) - 0.1;                
      //       }
      //    }else{
      //       // if energy not equal then jump on three lines
      //       getline(myfile,s);
      //       getline(myfile,s);
      //       getline(myfile,s);
      //    }
         
      // }
      // else {cout << "Unable to open file";} 

      //////////////////PUT DATA IN FILE///////////////////////////////////

      //myfile2_name.Form("difsig_%f_n_unpol_Th66.txt",Ed[i]);
      //myfile2.open(myfile2_name);
      //for(int i4=0; i4<num; i4++){
      //   myfile2 << theta_picture[i4]<< " "<< theta_error_picture[i4] <<\
      //   " "<< dcsn_picture[i4] - 0.1*(Ed[i]/25.) + 0.1<< " "<< dcsn_error_picture[i4] << "\n";
      //} 
      //myfile2.close();
      

//>>>///////////////////////////////////////////////////////////////////////

      // Choose only energies on the puicture
      //if ((i==0)||(i==2)||(i==3)||(i==5)||(i==6)||(i==8)||(i==10)\
      ||(i==11)||(i==12)||(i==14)){continue;}

      // Shift grafics up with 10 kev interval
      // int it;
      // it = it +1;
      // for(int i2=0; i2<nn;i2++){
      //    dcsn[i2] = dcsn[i2] + 0.1*it;
      //    dcsp[i2] = dcsp[i2] + 0.1*it;
      // }
   
      auto dcsng = new TGraphErrors(nn,theta,dcsn,0,0);
      dcsng->SetMarkerColor(1);
      dcsng->SetMarkerStyle(1);
      dcsng->SetLineColor(i+1);
      g_name.Form("Energy_%f",Ed[i]);
      dcsng->SetTitle(g_name);

      // auto dcsng_picture = new TGraphErrors(num,theta_picture,dcsn_picture,\
      // theta_error_picture,dcsn_error_picture);
      // dcsng_picture->SetMarkerColor(1);
      // dcsng_picture->SetMarkerStyle(2);
      // if (i == 9){dcsng_picture->SetLineColor(3);}
      // else{dcsng_picture->SetLineColor(i+1);}
      
      // Adding graphics
      
      mg1->Add(dcsng, "ACP"); 
     
      
      auto dcspg = new TGraphErrors(nn,theta,dcsp,0,0);
      dcspg->SetMarkerColor(1);
      dcspg->SetMarkerStyle(1);
      dcspg->SetLineColor(i+1);
      g_name.Form("Energy_%f",Ed[i]);
      dcspg->SetTitle(g_name);
     
      mg2->Add(dcspg, "ACP");
      
   }

   c1->cd(2);
   mg1->Draw("A");
   mg1->GetXaxis()->SetTitle("#it{#theta} (degres)");
   mg1->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
   gPad->Modified();
   mg1->SetMinimum(0.8);
   mg1->SetMaximum(3.5);
   //c1->BuildLegend(0.35,0.55,0.65,0.89);


   c1->cd(1);
   mg2->Draw("A");
   mg2->GetXaxis()->SetTitle("#it{#theta} (degres)");
   mg2->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
   gPad->Modified();
   mg2->SetMinimum(0.8);
   mg2->SetMaximum(3.5);
   //myfile.close();
}
