#include <TLatex.h>
#include <TMath.h>

void cross_sect_from_articles_25_50_kev() {
   auto c1 = new TCanvas("c1","Multigraph",700,500);
   //c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Divide(2);

   auto mg = new TMultiGraph();
   auto mg2 = new TMultiGraph();
     

   TLatex latex;
   int nn=180;
   double theta[nn], Kr_arr[nn], Br_arr[nn], Kr_n_arr[nn], Br_n_arr[nn];
   double dth = (3.1415926*0.5)/180;
   TString g_name;
   double norm_angle = 3.1415926*0.5;

   ////////////Th66
   // double dcsp_th66_26k[nn], edcsp_th66_26k[nn], CS_norm_th66_26k,\
   // Ap_th66_26k = 0.27, Bp_th66_26k = 0.0, eAp_th66_26k = 0.08, eBp_th66_26k = 0.;

   // double dcsp_th66_45k[nn], edcsp_th66_45k[nn], CS_norm_th66_45k,\
   // Ap_th66_45k = 0.34, Bp_th66_45k = -0.03, eAp_th66_45k = 0.02, eBp_th66_45k = 0.02;


   // double dcsn_th66_26k[nn], edcsn_th66_26k[nn], CSn_norm_th66_26k,\
   // An_th66_26k = 0.51, Bn_th66_26k = 0.0, eAn_th66_26k = 0.12, eBn_th66_26k = 0.;

   // double dcsn_th66_45k[nn], edcsn_th66_45k[nn], CSn_norm_th66_45k,\
   // An_th66_45k = 0.57, Bn_th66_45k = 0.04, eAn_th66_45k = 0.03, eBn_th66_45k = 0.;

   ////////////Kr86
   double dcsp_kr86_20k[nn], edcsp_kr86_20k[nn], CS_norm_kr86_20k,\
   Ap_kr86_20k = 0.14, Bp_kr86_20k = 0.04, eAp_kr86_20k = 0.01, eBp_kr86_20k = 0.03;

   double dcsp_kr86_30k[nn], edcsp_kr86_30k[nn], CS_norm_kr86_30k,\
   Ap_kr86_30k = 0.21, Bp_kr86_30k = 0.04, eAp_kr86_30k = 0.01, eBp_kr86_30k = 0.03;

   double dcsp_kr86_50k[nn], edcsp_kr86_50k[nn], CS_norm_kr86_50k, q2 = 0.999, q4 = 0.996,\
   Ap_kr86_50k = 0.23, Bp_kr86_50k = 0.03, eAp_kr86_50k = 0.01, eBp_kr86_50k = 0.03, leg_2, leg_4;


   double dcsn_kr86_20k[nn], edcsn_kr86_20k[nn], CSn_norm_kr86_20k,\
   An_kr86_20k = 0.19, Bn_kr86_20k = -0.04, eAn_kr86_20k = 0.04, eBn_kr86_20k = 0.05;

   double dcsn_kr86_30k[nn], edcsn_kr86_30k[nn], CSn_norm_kr86_30k,\
   An_kr86_30k = 0.35, Bn_kr86_30k = 0.02, eAn_kr86_30k = 0.04, eBn_kr86_30k = 0.05;

   double dcsn_kr86_50k[nn], edcsn_kr86_50k[nn], CSn_norm_kr86_50k,\
   An_kr86_50k = 0.39, Bn_kr86_50k = 0.02, eAn_kr86_50k = 0.03, eBn_kr86_50k = 0.05;


   ////////////Br90
   double dcsp_br90_20k[nn], edcsp_br90_20k[nn], CS_norm_br90_20k,\
   Ap_br90_20k = 0.0208, Bp_br90_20k = 0.0023, eAp_br90_20k = 0.0010, eBp_br90_20k = 0.0022;

   double dcsp_br90_30k[nn], edcsp_br90_30k[nn], CS_norm_br90_30k,\
   Ap_br90_30k = 0.0886, Bp_br90_30k = 0.0184, eAp_br90_30k = 0.0020, eBp_br90_30k = 0.0044;

   double dcsp_br90_50k[nn], edcsp_br90_50k[nn], CS_norm_br90_50k,\
   Ap_br90_50k = 0.3215, Bp_br90_50k = 0.076, eAp_br90_50k = 0.0048, eBp_br90_50k = 0.027;


   double dcsn_br90_20k[nn], edcsn_br90_20k[nn], CSn_norm_br90_20k,\
   An_br90_20k = 0.0181, Bn_br90_20k = 0.0108, eAn_br90_20k = 0.0014, eBn_br90_20k = 0.0031;

   double dcsn_br90_30k[nn], edcsn_br90_30k[nn], CSn_norm_br90_30k,\
   An_br90_30k = 0.0782, Bn_br90_30k = 0.0425, eAn_br90_30k = 0.0028, eBn_br90_30k = 0.0066;

   double dcsn_br90_50k[nn], edcsn_br90_50k[nn], CSn_norm_br90_50k,\
   An_br90_50k = 0.2994, Bn_br90_50k = 0.212, eAn_br90_50k = 0.0059, eBn_br90_50k = 0.014;

   ////////////We52
   // double dcsp_we52_35k[nn], edcsp_we52_35k[nn], CS_norm_we52_35k,\
   // Ap_we52_35k = 0.155, Bp_we52_35k = 0.02, coef = 1.94/(4*3.1415926);

   // double dcsp_we52_50k[nn], edcsp_we52_50k[nn], CS_norm_we52_50k,\
   // Ap_we52_50k = 0.195, Bp_we52_50k = 0.02, coef50 = 4.81/(4*3.1415926);
   
   
   
   for (int j = 0; j < nn; j++){  
      double theta0 = 1.0*j;
      theta[j] = theta0;
      double th = (3.1415926*theta0)/180;

      // Th66
      // CS_norm_th66_26k = 1 + Ap_th66_26k*pow(cos(norm_angle),2) +\
      // Bp_th66_26k*pow(cos(norm_angle),4);
      // dcsp_th66_26k[j]  = (1 + Ap_th66_26k*pow(cos(th),2) +\
      // Bp_th66_26k*pow(cos(th),4))/(CS_norm_th66_26k);
      // edcsp_th66_26k[j] = sqrt(pow(pow(cos(th),2)*eAp_th66_26k,2) +\
      // pow(pow(cos(th),4)*eBp_th66_26k,2) + pow(Ap_th66_26k*2*cos(th)*(-sin(th))*dth +\
      // Bp_th66_26k*4*pow(cos(th),3)*(-sin(th))*dth,2));

      // CS_norm_th66_45k = 1 + Ap_th66_45k*pow(cos(norm_angle),2) +\
      // Bp_th66_45k*pow(cos(norm_angle),4);
      // dcsp_th66_45k[j]  = (1 + Ap_th66_45k*pow(cos(th),2) +\
      // Bp_th66_45k*pow(cos(th),4))/(CS_norm_th66_45k);
      // edcsp_th66_45k[j] = sqrt(pow(pow(cos(th),2)*eAp_th66_45k,2) +\
      // pow(pow(cos(th),4)*eBp_th66_45k,2) + pow(Ap_th66_45k*2*cos(th)*(-sin(th))*dth +\
      // Bp_th66_45k*4*pow(cos(th),3)*(-sin(th))*dth,2));


      // CSn_norm_th66_26k = 1 + An_th66_26k*pow(cos(norm_angle),2) +\
      // Bn_th66_26k*pow(cos(norm_angle),4);
      // dcsn_th66_26k[j]  = (1 + An_th66_26k*pow(cos(th),2) +\
      // Bn_th66_26k*pow(cos(th),4))/(CSn_norm_th66_26k);
      // edcsn_th66_26k[j] = sqrt(pow(pow(cos(th),2)*eAn_th66_26k,2) +\
      // pow(pow(cos(th),4)*eBn_th66_26k,2) + pow(An_th66_26k*2*cos(th)*(-sin(th))*dth +\
      // Bn_th66_26k*4*pow(cos(th),3)*(-sin(th))*dth,2));

      // CSn_norm_th66_45k = 1 + An_th66_45k*pow(cos(norm_angle),2) +\
      // Bn_th66_45k*pow(cos(norm_angle),4);
      // dcsn_th66_45k[j]  = (1 + An_th66_45k*pow(cos(th),2) +\
      // Bn_th66_45k*pow(cos(th),4))/(CSn_norm_th66_45k);
      // edcsn_th66_45k[j] = sqrt(pow(pow(cos(th),2)*eAn_th66_45k,2) +\
      // pow(pow(cos(th),4)*eBn_th66_45k,2) + pow(An_th66_45k*2*cos(th)*(-sin(th))*dth +\
      // Bn_th66_45k*4*pow(cos(th),3)*(-sin(th))*dth,2));

      //Kr86
      leg_2 = 0.5*(3.0*pow(cos(th),2)-1);
      leg_4 = 4.0628*pow(cos(th),4) - 2.375*pow(cos(th),2) + 0.3125;
      double leg_2_N = 0.5*(3.0*pow(cos(norm_angle),2)-1);
      double leg_4_N = 4.0628*pow(cos(norm_angle),4) - 2.375*pow(cos(norm_angle),2) + 0.3125;

      CS_norm_kr86_20k = 1 + Ap_kr86_20k*q2*leg_2_N + Bp_kr86_20k*q4*leg_4_N;
      dcsp_kr86_20k[j]  = (1 + Ap_kr86_20k*q2*leg_2 + Bp_kr86_20k*q4*leg_4)/(CS_norm_kr86_20k);
      edcsp_kr86_20k[j] = sqrt(pow(q2*leg_2*eAp_kr86_20k,2) + pow(q4*leg_4*eBp_kr86_20k,2) +\
      pow(Ap_kr86_20k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bp_kr86_20k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));

      CS_norm_kr86_30k = 1 + Ap_kr86_30k*q2*leg_2_N + Bp_kr86_30k*q4*leg_4_N;
      dcsp_kr86_30k[j]  = (1 + Ap_kr86_30k*q2*leg_2 + Bp_kr86_30k*q4*leg_4)/(CS_norm_kr86_30k);
      edcsp_kr86_30k[j] = sqrt(pow(q2*leg_2*eAp_kr86_30k,2) + pow(q4*leg_4*eBp_kr86_30k,2) +\
      pow(Ap_kr86_30k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bp_kr86_30k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));

      CS_norm_kr86_50k = 1 + Ap_kr86_50k*q2*leg_2_N + Bp_kr86_50k*q4*leg_4_N;
      dcsp_kr86_50k[j]  = (1 + Ap_kr86_50k*q2*leg_2 + Bp_kr86_50k*q4*leg_4)/(CS_norm_kr86_50k);
      edcsp_kr86_50k[j] = sqrt(pow(q2*leg_2*eAp_kr86_50k,2) + pow(q4*leg_4*eBp_kr86_50k,2) +\
      pow(Ap_kr86_50k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bp_kr86_50k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));


      CSn_norm_kr86_20k = 1 + An_kr86_20k*q2*leg_2_N + Bn_kr86_20k*q4*leg_4_N;
      dcsn_kr86_20k[j]  = (1 + An_kr86_20k*q2*leg_2 + Bn_kr86_20k*q4*leg_4)/(CSn_norm_kr86_20k);
      edcsn_kr86_20k[j] = sqrt(pow(q2*leg_2*eAn_kr86_20k,2) + pow(q4*leg_4*eBn_kr86_20k,2) +\
      pow(An_kr86_20k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bn_kr86_20k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));

      CSn_norm_kr86_30k = 1 + An_kr86_30k*q2*leg_2_N + Bn_kr86_30k*q4*leg_4_N;
      dcsn_kr86_30k[j]  = (1 + An_kr86_30k*q2*leg_2 + Bn_kr86_30k*q4*leg_4)/(CSn_norm_kr86_30k);
      edcsn_kr86_30k[j] = sqrt(pow(q2*leg_2*eAn_kr86_30k,2) + pow(q4*leg_4*eBn_kr86_30k,2) +\
      pow(An_kr86_30k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bn_kr86_30k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));

      CSn_norm_kr86_50k = 1 + An_kr86_50k*q2*leg_2_N + Bn_kr86_50k*q4*leg_4_N;
      dcsn_kr86_50k[j]  = (1 + An_kr86_50k*q2*leg_2 + Bn_kr86_50k*q4*leg_4)/(CSn_norm_kr86_50k);
      edcsn_kr86_50k[j] = sqrt(pow(q2*leg_2*eAn_kr86_50k,2) + pow(q4*leg_4*eBn_kr86_50k,2) +\
      pow(An_kr86_50k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bn_kr86_50k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));


      //Br90
      CS_norm_br90_20k = Ap_br90_20k + Bp_br90_20k*pow(cos(norm_angle),2);
      dcsp_br90_20k[j]  = (Ap_br90_20k + Bp_br90_20k*pow(cos(th),2))/(CS_norm_br90_20k);
      edcsp_br90_20k[j] = sqrt(pow(eAp_br90_20k,2) + pow(pow(cos(th),2)*eBp_br90_20k,2) +\
      pow(Bp_br90_20k*2*cos(th)*(-sin(th))*dth,2)); 

      CS_norm_br90_30k = Ap_br90_30k + Bp_br90_30k*pow(cos(norm_angle),2);
      dcsp_br90_30k[j]  = (Ap_br90_30k + Bp_br90_30k*pow(cos(th),2))/(CS_norm_br90_30k);
      edcsp_br90_30k[j] = sqrt(pow(eAp_br90_30k,2) + pow(pow(cos(th),2)*eBp_br90_30k,2) +\
      pow(Bp_br90_30k*2*cos(th)*(-sin(th))*dth,2)); 

      CS_norm_br90_50k = Ap_br90_50k + Bp_br90_50k*pow(cos(norm_angle),2);
      dcsp_br90_50k[j]  = (Ap_br90_50k + Bp_br90_50k*pow(cos(th),2))/(CS_norm_br90_50k);
      edcsp_br90_50k[j] = sqrt(pow(eAp_br90_50k,2) + pow(pow(cos(th),2)*eBp_br90_50k,2) +\
      pow(Bp_br90_50k*2*cos(th)*(-sin(th))*dth,2));


      CSn_norm_br90_20k = An_br90_20k + Bn_br90_20k*pow(cos(norm_angle),2);
      dcsn_br90_20k[j]  = (An_br90_20k + Bn_br90_20k*pow(cos(th),2))/(CSn_norm_br90_20k);
      edcsn_br90_20k[j] = sqrt(pow(eAn_br90_20k,2) + pow(pow(cos(th),2)*eBn_br90_20k,2) +\
      pow(Bn_br90_20k*2*cos(th)*(-sin(th))*dth,2)); 

      CSn_norm_br90_30k = An_br90_30k + Bn_br90_30k*pow(cos(norm_angle),2);
      dcsn_br90_30k[j]  = (An_br90_30k + Bn_br90_30k*pow(cos(th),2))/(CSn_norm_br90_30k);
      edcsn_br90_30k[j] = sqrt(pow(eAn_br90_30k,2) + pow(pow(cos(th),2)*eBn_br90_30k,2) +\
      pow(Bn_br90_30k*2*cos(th)*(-sin(th))*dth,2)); 

      CSn_norm_br90_50k = An_br90_50k + Bn_br90_50k*pow(cos(norm_angle),2);
      dcsn_br90_50k[j]  = (An_br90_50k + Bn_br90_50k*pow(cos(th),2))/(CSn_norm_br90_50k);
      edcsn_br90_50k[j] = sqrt(pow(eAn_br90_50k,2) + pow(pow(cos(th),2)*eBn_br90_50k,2) +\
      pow(Bn_br90_50k*2*cos(th)*(-sin(th))*dth,2)); 


      //We52
      // CS_norm_we52_35k = coef*(1 + Ap_we52_35k*leg_2_N + Bp_we52_35k*leg_4_N);
      // dcsp_we52_35k[j]  = (coef*(1 + Ap_we52_35k*leg_2 + Bp_we52_35k*leg_4))/(CS_norm_we52_35k);

      // CS_norm_we52_50k = coef50*(1 + Ap_we52_50k*leg_2_N + Bp_we52_50k*leg_4_N);
      // dcsp_we52_50k[j]  = (coef50*(1 + Ap_we52_50k*leg_2 + Bp_we52_50k*leg_4))/(CS_norm_we52_50k);


   }
   
   for (int i2=0; i2<nn; i2++){
      Kr_arr[i2] = dcsp_kr86_50k[i2]/dcsp_kr86_20k[i2];
      Br_arr[i2] = dcsp_br90_50k[i2]/dcsp_br90_20k[i2];
      Kr_n_arr[i2] = dcsn_kr86_50k[i2]/dcsn_kr86_20k[i2];
      Br_n_arr[i2] = dcsn_br90_50k[i2]/dcsn_br90_20k[i2];
      //printf("Kr_arr_%f\n", Kr_arr[i2]);
      //printf("Br_arr_%f\n", Br_arr[i2]);
   }

   TGraph* gr2 = new TGraph(nn, theta, Kr_arr);
   gr2->SetMarkerColor(4);
   gr2->SetMarkerStyle(23);
   gr2->SetLineColor(4);
   gr2->SetTitle("Krauss86 50/20");
   TGraph* gr3 = new TGraph(nn, theta, Br_arr);
   gr3->SetMarkerColor(1);
   gr3->SetMarkerStyle(20);
   gr3->SetLineColor(1);
   gr3->SetTitle("Brown90 50/20");

   TGraph* gr4 = new TGraph(nn, theta, Kr_n_arr);
   gr4->SetMarkerColor(2);
   gr4->SetMarkerStyle(21);
   gr4->SetLineColor(2);
   gr4->SetTitle("Krauss86 n 50/20");
   TGraph* gr5 = new TGraph(nn, theta, Br_n_arr);
   gr5->SetMarkerColor(3);
   gr5->SetMarkerStyle(22);
   gr5->SetLineColor(3);
   gr5->SetTitle("Brown90 n 50/20");

   // auto gr_th66_45k = new TGraphErrors(nn,theta,dcsp_th66_45k,0,0);
   // gr_th66_45k->SetMarkerColor(1);
   // gr_th66_45k->SetMarkerStyle(20);
   // gr_th66_45k->SetLineColor(1);
   // gr_th66_45k->SetTitle("Th66 45.1 keV");

   // auto gr_th66_26k = new TGraphErrors(nn,theta,dcsp_th66_26k,0,0);
   // gr_th66_26k->SetMarkerColor(2);
   // gr_th66_26k->SetMarkerStyle(21);
   // gr_th66_26k->SetLineColor(2);
   // gr_th66_26k->SetTitle("Th66 26.9 keV");

   
   auto gr_kr86_20k = new TGraphErrors(nn,theta,dcsp_kr86_20k,0,edcsp_kr86_20k);
   gr_kr86_20k->SetMarkerColor(4);
   gr_kr86_20k->SetMarkerStyle(23);
   gr_kr86_20k->SetLineColor(4);
   gr_kr86_20k->SetTitle("Krauss86 20 keV");

   auto gr_kr86_30k = new TGraphErrors(nn,theta,dcsp_kr86_30k,0,edcsp_kr86_30k);
   gr_kr86_30k->SetMarkerColor(6);
   gr_kr86_30k->SetMarkerStyle(24);
   gr_kr86_30k->SetLineColor(6);
   gr_kr86_30k->SetTitle("Krauss86 30 keV");

   auto gr_kr86_50k = new TGraphErrors(nn,theta,dcsp_kr86_50k,0,edcsp_kr86_50k);
   gr_kr86_50k->SetMarkerColor(3);
   gr_kr86_50k->SetMarkerStyle(22);
   gr_kr86_50k->SetLineColor(3);
   gr_kr86_50k->SetTitle("Krauss86 50 keV");

   auto gr_br90_20k = new TGraphErrors(nn,theta,dcsp_br90_20k,0,edcsp_br90_20k);
   gr_br90_20k->SetMarkerColor(7);
   gr_br90_20k->SetMarkerStyle(25);
   gr_br90_20k->SetLineColor(7);
   gr_br90_20k->SetTitle("Brown90 20 keV");

   auto gr_br90_30k = new TGraphErrors(nn,theta,dcsp_br90_30k,0,edcsp_br90_30k);
   gr_br90_30k->SetMarkerColor(8);
   gr_br90_30k->SetMarkerStyle(26);
   gr_br90_30k->SetLineColor(8);
   gr_br90_30k->SetTitle("Brown90 30 keV");

   auto gr_br90_50k = new TGraphErrors(nn,theta,dcsp_br90_50k,0,edcsp_br90_50k);
   gr_br90_50k->SetMarkerColor(5);
   gr_br90_50k->SetMarkerStyle(27);
   gr_br90_50k->SetLineColor(5);
   gr_br90_50k->SetTitle("Brown90 50 keV");

   // auto gr_we52_35k = new TGraphErrors(nn,theta,dcsp_we52_35k,0,0);
   // gr_we52_35k->SetMarkerColor(9);
   // gr_we52_35k->SetMarkerStyle(27);
   // gr_we52_35k->SetLineColor(9);
   // gr_we52_35k->SetTitle("Wenzel 35 keV");

   // auto gr_we52_50k = new TGraphErrors(nn,theta,dcsp_we52_50k,0,0);
   // gr_we52_50k->SetMarkerColor(12);
   // gr_we52_50k->SetMarkerStyle(28);
   // gr_we52_50k->SetLineColor(12);
   // gr_we52_50k->SetTitle("Wenzel 50 keV");
   
   // mg->Add(gr_th66_45k);
   // mg->Add(gr_th66_26k);
   // mg->Add(gr_kr86_50k);
   // mg->Add(gr_kr86_20k);
   // mg->Add(gr_kr86_30k);
   // mg->Add(gr_br90_20k);
   // mg->Add(gr_br90_30k);
   // mg->Add(gr_br90_50k);
   // mg->Add(gr_we52_35k);
   // mg->Add(gr_we52_50k);

   mg->Add(gr2);
   mg->Add(gr3);

   c1->cd(1);
   //mg->SetTitle("Cross section ^{2}H(d,p)^{3}H for energies 25, 50 kev normalize to CS(90 degres)");
   mg->SetTitle("Relation for cross sections ^{2}H(d,p)^{3}H for energies 25, 50 kev normalize to CS(90 degres)");
   mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
   mg->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
   gPad->Modified();
   //mg->SetMinimum(0.78);
   //mg->SetMaximum(1.2);
   mg->Draw("ACLP");

   //printf("AA %d", mg->GetListOfGraphs()->GetEntries());
   c1->cd(1)->BuildLegend(0.35,0.55,0.65,0.89);





   // auto gr_th66_n45k = new TGraphErrors(nn,theta,dcsn_th66_45k,0,0);
   // gr_th66_n45k->SetMarkerColor(1);
   // gr_th66_n45k->SetMarkerStyle(20);
   // gr_th66_n45k->SetLineColor(1);
   // gr_th66_n45k->SetTitle("Th66 n 45.1 keV");

   // auto gr_th66_n26k = new TGraphErrors(nn,theta,dcsn_th66_26k,0,0);
   // gr_th66_n26k->SetMarkerColor(2);
   // gr_th66_n26k->SetMarkerStyle(21);
   // gr_th66_n26k->SetLineColor(2);
   // gr_th66_n26k->SetTitle("Th66 n 26.9 keV");

   auto gr_kr86_n20k = new TGraphErrors(nn,theta,dcsn_kr86_20k,0,edcsn_kr86_20k);
   gr_kr86_n20k->SetMarkerColor(4);
   gr_kr86_n20k->SetMarkerStyle(23);
   gr_kr86_n20k->SetLineColor(4);
   gr_kr86_n20k->SetTitle("Krauss86 n 20 keV");

   auto gr_kr86_n30k = new TGraphErrors(nn,theta,dcsn_kr86_30k,0,edcsn_kr86_30k);
   gr_kr86_n30k->SetMarkerColor(6);
   gr_kr86_n30k->SetMarkerStyle(24);
   gr_kr86_n30k->SetLineColor(6);
   gr_kr86_n30k->SetTitle("Krauss86 n 30 keV");

   auto gr_kr86_n50k = new TGraphErrors(nn,theta,dcsn_kr86_50k,0,edcsn_kr86_50k);
   gr_kr86_n50k->SetMarkerColor(3);
   gr_kr86_n50k->SetMarkerStyle(22);
   gr_kr86_n50k->SetLineColor(3);
   gr_kr86_n50k->SetTitle("Krauss86 n 50 keV");

   auto gr_br90_n20k = new TGraphErrors(nn,theta,dcsn_br90_20k,0,edcsn_br90_20k);
   gr_br90_n20k->SetMarkerColor(7);
   gr_br90_n20k->SetMarkerStyle(25);
   gr_br90_n20k->SetLineColor(7);
   gr_br90_n20k->SetTitle("Brown90 n 20 keV");

   auto gr_br90_n30k = new TGraphErrors(nn,theta,dcsn_br90_30k,0,edcsn_br90_30k);
   gr_br90_n30k->SetMarkerColor(8);
   gr_br90_n30k->SetMarkerStyle(26);
   gr_br90_n30k->SetLineColor(8);
   gr_br90_n30k->SetTitle("Brown90 n 30 keV");

   auto gr_br90_n50k = new TGraphErrors(nn,theta,dcsn_br90_50k,0,edcsn_br90_50k);
   gr_br90_n50k->SetMarkerColor(5);
   gr_br90_n50k->SetMarkerStyle(27);
   gr_br90_n50k->SetLineColor(5);
   gr_br90_n50k->SetTitle("Brown90 n 50 keV");
   
   // mg2->Add(gr_th66_n45k);
   // mg2->Add(gr_th66_n26k);
   // mg2->Add(gr_kr86_n50k);
   // mg2->Add(gr_kr86_n20k);
   // mg2->Add(gr_kr86_n30k);
   // mg2->Add(gr_br90_n20k);
   // mg2->Add(gr_br90_n30k);
   // mg2->Add(gr_br90_n50k);

   mg2->Add(gr4);
   mg2->Add(gr5);


   c1->cd(2);
   //mg2->SetTitle("Cross section ^{2}H(d,n)^{3}He for energies 25, 50 kev normalize to CS(90 degres)");
   mg2->SetTitle("Relation of cross sections ^{2}H(d,n)^{3}He for energies 25, 50 kev normalize to CS(90 degres)");
   mg2->GetXaxis()->SetTitle("#it{#theta} (degres)");
   mg2->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
   gPad->Modified();
   //mg->SetMinimum(0.78);
   //mg->SetMaximum(1.2);
   mg2->Draw("ACLP");

   //printf("AA %d", mg->GetListOfGraphs()->GetEntries());
   c1->cd(2)->BuildLegend(0.35,0.55,0.65,0.89);



}
