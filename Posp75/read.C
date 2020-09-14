#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int read() {
  string line;
  string s;
  ifstream myfile ("Th66_Data_dsigm[theta_cm]_DdnHe.txt");

  const Int_t n = 12;
  Double_t Ed[n]  = {19.5, 26.9, 32.0, 45.1, 71.0, 96.6, 122.0, 147.5, 194.0, 248.3, 298.5, 348.7};
  Double_t dd;
  const Int_t nn = 15;
  Double_t theta[nn];
  Double_t dcsn[nn];
  

  if (myfile.is_open())
  {
    //myfile.getline(s,nn,'\n');
    //cout << s << endl;
    for (int i =0; i<n; i++){
        if ((i==0)||(i==2)){continue;}
        myfile >>dd;
        printf("Ed_%f\n",Ed[i]);
        printf("dd_%f\n",dd);
        if (Ed[i]==dd){

            for (int i =0; i<nn; i++){
                //myfile.getline(buffer,nn,' ');
                myfile >> theta[i];
                cout << theta[i] << endl;  
            }
            for (int i =0; i<nn; i++){
                myfile >> dcsn[i];
                cout << dcsn[i] << endl;  

            }
        }else{
            getline(myfile,s);
            getline(myfile,s);
            getline(myfile,s);
        }
    }
    

    //while ( getline (myfile,line) )
    //{
        //cout << typeid(line).name() << endl;
        //for (int i = 0; i < n; i++){
        //    if (Ed[i]== std::stod(line)){cout << line << '\n';}
        //}
    //}
    myfile.close();
  }

  else cout << "Unable to open file"; 

  return 0;
}
