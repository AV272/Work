#include <TLatex.h>
#include <TMath.h>

void cross_sect_normalized() {
   auto c1 = new TCanvas("c1","Multigraph",700,500);
   //c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);

   auto mg = new TMultiGraph();
   

   TLatex latex;
   const Int_t nn=180;
   Double_t theta[nn];
   Double_t dth = (3.1415926*0.5)/180;
   TString g_name;
   double norm_angle = 2.8;

   Double_t dcsp_th66_19k[nn], edcsp_th66_19k[nn], CS_norm_th66_19k,\
   Ap_th66_19k = 0.2, Bp_th66_19k = 0.0, eAp_th66_19k = 0.15, eBp_th66_19k = 0.;

   Double_t dcsp_th66_26k[nn], edcsp_th66_26k[nn], CS_norm_th66_26k,\
   Ap_th66_26k = 0.27, Bp_th66_26k = 0.0, eAp_th66_26k = 0.08, eBp_th66_26k = 0.;

   Double_t dcsp_kr86_15k[nn], edcsp_kr86_15k[nn], CS_norm_kr86_15k, q2 = 0.999, q4 = 0.996,\
   Ap_kr86_15k = 0.05, Bp_kr86_15k = 0.03, eAp_kr86_15k = 0.04, eBp_kr86_15k = 0., leg_2, leg_4;

   Double_t dcsp_kr86_20k[nn], edcsp_kr86_20k[nn], CS_norm_kr86_20k,\
   Ap_kr86_20k = 0.14, Bp_kr86_20k = 0.04, eAp_kr86_20k = 0.01, eBp_kr86_20k = 0.03;

   Double_t dcsp_kr86_30k[nn], edcsp_kr86_30k[nn], CS_norm_kr86_30k,\
   Ap_kr86_30k = 0.21, Bp_kr86_30k = 0.04, eAp_kr86_30k = 0.01, eBp_kr86_30k = 0.03;
   
   
   
   for (int j = 0; j < nn; j++){  
      Double_t theta0 = j;
      theta[j] = theta0;
      Double_t th = (3.1415926*theta0)/180;

      // Th66
      CS_norm_th66_19k = 1 + Ap_th66_19k*pow(cos(norm_angle),2) +\
      Bp_th66_19k*pow(cos(norm_angle),4);
      dcsp_th66_19k[j]  = (1 + Ap_th66_19k*pow(cos(th),2) +\
      Bp_th66_19k*pow(cos(th),4))/(CS_norm_th66_19k);
      edcsp_th66_19k[j] = sqrt(pow(pow(cos(th),2)*eAp_th66_19k,2) +\
      pow(pow(cos(th),4)*eBp_th66_19k,2) + pow(Ap_th66_19k*2*cos(th)*(-sin(th))*dth +\
      Bp_th66_19k*4*pow(cos(th),3)*(-sin(th))*dth,2));

      CS_norm_th66_26k = 1 + Ap_th66_26k*pow(cos(norm_angle),2) +\
      Bp_th66_26k*pow(cos(norm_angle),4);
      dcsp_th66_26k[j]  = (1 + Ap_th66_26k*pow(cos(th),2) +\
      Bp_th66_26k*pow(cos(th),4))/(CS_norm_th66_26k);
      edcsp_th66_26k[j] = sqrt(pow(pow(cos(th),2)*eAp_th66_26k,2) +\
      pow(pow(cos(th),4)*eBp_th66_26k,2) + pow(Ap_th66_26k*2*cos(th)*(-sin(th))*dth +\
      Bp_th66_26k*4*pow(cos(th),3)*(-sin(th))*dth,2));

      //Kr86
      leg_2 = 0.5*(3.0*pow(cos(th),2)-1);
      leg_4 = 4.0628*pow(cos(th),4) - 2.375*pow(cos(th),2) + 0.3125;
      double leg_2_N = 0.5*(3.0*pow(cos(norm_angle),2)-1);
      double leg_4_N = 4.0628*pow(cos(norm_angle),4) - 2.375*pow(cos(norm_angle),2) + 0.3125;

      CS_norm_kr86_15k = 1 + Ap_kr86_15k*q2*leg_2_N + Bp_kr86_15k*q4*leg_4_N;
      dcsp_kr86_15k[j]  = (1 + Ap_kr86_15k*q2*leg_2 + Bp_kr86_15k*q4*leg_4)/(CS_norm_kr86_15k);
      edcsp_kr86_15k[j] = sqrt(pow(q2*leg_2*eAp_kr86_15k,2) + pow(q4*leg_4*eBp_kr86_15k,2) +\
      pow(Ap_kr86_15k*q2*(0.5*(3*2*cos(th)*(-sin(th))*dth)-1) +\
      Bp_kr86_15k*q4*(4.0628*4*pow(cos(th),3)*(-sin(th))*dth - 2.375*2*cos(th)*(-sin(th))*dth +0.3125),2));

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

         

   }
   
   //auto gr_th66_19k = new TGraphErrors(nn,theta,dcsp_th66_19k,0,edcsp_th66_19k);
   auto gr_th66_19k = new TGraphErrors(nn,theta,dcsp_th66_19k,0,0);
   gr_th66_19k->SetMarkerColor(1);
   gr_th66_19k->SetMarkerStyle(20);
   gr_th66_19k->SetLineColor(1);
   gr_th66_19k->SetTitle("Th66 19.5 keV");

   auto gr_th66_26k = new TGraphErrors(nn,theta,dcsp_th66_26k,0,0);
   gr_th66_26k->SetMarkerColor(2);
   gr_th66_26k->SetMarkerStyle(21);
   gr_th66_26k->SetLineColor(2);
   gr_th66_26k->SetTitle("Th66 26.9 keV");

   auto gr_kr86_15k = new TGraphErrors(nn,theta,dcsp_kr86_15k,0,0);
   gr_kr86_15k->SetMarkerColor(3);
   gr_kr86_15k->SetMarkerStyle(22);
   gr_kr86_15k->SetLineColor(3);
   gr_kr86_15k->SetTitle("Krauss86 15 keV");

   auto gr_kr86_20k = new TGraphErrors(nn,theta,dcsp_kr86_20k,0,0);
   gr_kr86_20k->SetMarkerColor(4);
   gr_kr86_20k->SetMarkerStyle(23);
   gr_kr86_20k->SetLineColor(4);
   gr_kr86_20k->SetTitle("Krauss86 20 keV");

   auto gr_kr86_30k = new TGraphErrors(nn,theta,dcsp_kr86_30k,0,0);
   gr_kr86_30k->SetMarkerColor(6);
   gr_kr86_30k->SetMarkerStyle(24);
   gr_kr86_30k->SetLineColor(6);
   gr_kr86_30k->SetTitle("Krauss86 30 keV");
   
   mg->Add(gr_th66_19k);
   mg->Add(gr_th66_26k);
   mg->Add(gr_kr86_15k);
   mg->Add(gr_kr86_20k);
   mg->Add(gr_kr86_30k);

   mg->SetTitle("Cross section ^{2}H(d,p)^{3}H normalize to CS(160(2.8) degres)");
   mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
   mg->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (arbitrary units)");
   gPad->Modified();
   //mg->SetMinimum(0.78);
   //mg->SetMaximum(1.2);
   mg->Draw("ACLP");

   //printf("AA %d", mg->GetListOfGraphs()->GetEntries());
   c1->BuildLegend(0.35,0.55,0.65,0.89);





}
