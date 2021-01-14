#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void Th66_dcs_real_data(){

    auto c1 = new TCanvas("c1","Multigraph",700,500);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);

    auto mg = new TMultiGraph();
    mg->SetTitle("Absolute differential cross sections for ^{3}He,\
    normalaze to the absolute cross section taken by Ruby and Crawford");

    // Read DATA file
    ifstream myfile ("Th66_Data_DA[theta_cm]_EXP_He3.txt");

    int num = 20;
    double buffer;
    double Ed, Ed0;
    double ang_cm[num], DA[num], DA_err[num], Ed_arr[num];
    TString g_name;

    if (myfile.is_open()){
        
        int n = 12;
        int nn;
        myfile >> buffer;
        Ed = buffer;
        for (int i=0; i<n; i++){
            
            Ed0 = Ed;
            Ed_arr[i] = Ed0;
            printf("Ed0_%f\n",Ed0);
            printf("Ed_%f\n", Ed);
            int i1 = 0;
            while (Ed == Ed0){
                myfile >> ang_cm[i1];
                myfile >> DA[i1];
                myfile >> DA_err[i1];
                if (myfile.eof()){break;}
                myfile >> Ed;
                i1 = i1 +1;
            }

            // printf("Energy_%f\n", Ed0);
            // for (int i2=0; i2<num; i2++){
            //     if (DA0[i2]==0){
            //         nn=i2;
            //         break;
            //     }
                
            // }
            // double DA[nn], ang_cm[nn], DA_err[nn];
            // for (int i2=0; i2<nn; i2++){
            //     ang_cm[i2] = ang_cm0[i2];
            //     DA[i2] = DA0[i2];
            //     DA_err[i2] = DA_err0[i2]; 
            // }
            // cout << "\n";
            // for (int i2=0; i2<num; i2++){
            //     printf("%f_", DA_err[i2]);
            // }
            // cout << "\n";
            
            

            auto dcsng = new TGraphErrors(num, ang_cm, DA, 0, DA_err);
            dcsng->SetMarkerColor(i+1);
            dcsng->SetMarkerStyle(20+i);
            if (i == 9){dcsng->SetLineColor(3);}
            else{dcsng->SetLineColor(i+1);}
            g_name.Form("Energy_%f",Ed_arr[i]);
            dcsng->SetTitle(g_name);
      
            // Adding graphics
            mg->Add(dcsng, "AP");
            mg->Draw("AP");
            mg->GetXaxis()->SetTitle("#it{#theta} (degres)");
            mg->GetYaxis()->SetTitle("#it{#frac{d#sigma}{d#omega} (#theta)} (mb/sr)");
            gPad->Modified();
            //mg->GetXaxis()->SetLimits(1.5,7.5);
            mg->SetMinimum(0.0);
            mg->SetMaximum(100);

            // Clear arrays (fill by zeros)
            std::fill_n(ang_cm, num, 0);
            std::fill_n(DA, num, 0);
            std::fill_n(DA_err, num, 0);
            if (myfile.eof()){break;}
        }
        
        

         
    }else {cout << "Unable to open file";} 

    c1->BuildLegend(0.35,0.55,0.65,0.89);
    myfile.close();


}