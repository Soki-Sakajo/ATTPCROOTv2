#include "TMath.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
static double sq(double x){ return x*x; }

bool Readfile(const char* fname, Double_t *n, Double_t *m1, Double_t *m2, Double_t *m3, Double_t *m4, Double_t *Ex3, Double_t *Ex4, Double_t *K1){
  ifstream infi(fname);
  if(!infi){
    cout << "Cannot open input file: " << fname << endl;
    gSystem->Exit(0);
  }

  string line;
  Int_t nline=0;
  Double_t val;

  while(getline(infi,line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    nline++;
    istringstream iss(line);
    if(nline==1){
      iss>> *m1 >> *m2 >> *m3 >> *m4;
    }
    else if(nline==2){
      iss>> *Ex3 >> *Ex4;
    }
    else if(nline==3){
      iss>> *K1;
    }
    else if(nline==4){
      iss>> *n;
    }
  }
  return true;
}

void kinema(const char* filename=""){

  TString indi ="./in_file/";
  TString outdi="./out_file/";
  TString infi= indi+filename;
  TString oufi=outdi+filename;
  Double_t n=1800;
  Double_t m1,m2,m3,m4,Ex3,Ex4,K1_lab;
  Double_t dtheta= TMath::Pi()/n;
  //  cout<<"dtheta:"<<dtheta<<endl;
  
  if(!Readfile(infi, &n, &m1, &m2, &m3, &m4, &Ex3, &Ex4, &K1_lab)){
    return;
    gSystem->Exit(0);
  }
  Readfile(infi, &n, &m1, &m2, &m3, &m4, &Ex3, &Ex4, &K1_lab);

  //calcurate mass include Ex
  Double_t m3p = m3 + Ex3;
  Double_t m4p = m4 + Ex4;

  //starting condition
  Double_t E1_lab = K1_lab + m1;
  Double_t p1_lab = sqrt(max(0.0, E1_lab*E1_lab - m1*m1));
  Double_t s = m1*m1 + m2*m2 + 2.0*m2*E1_lab; //Lorentz invariant
  Double_t sqrt_s = sqrt(s);

  if(s < sq(m3p+m4p)){
    cout << "Reaction forbidden: sqrt(s) < m3'+m4'" << endl;
    return;
  }

  Double_t term1 = s - sq(m3p + m4p);
  Double_t term2 = s - sq(m3p - m4p);
  Double_t p_cm = 0.5 * sqrt(max(0.0, term1*term2)) / sqrt_s;

  Double_t E3_cm = sqrt(sq(p_cm) + sq(m3p));
  Double_t E4_cm = sqrt(sq(p_cm) + sq(m4p));

  //CM to Lab
  Double_t P_lab = p1_lab;
  Double_t E_lab_total = E1_lab + m2;
  Double_t beta = P_lab / E_lab_total;
  Double_t gamma = 1.0 / sqrt(1.0 - beta*beta);

  //set output
  ofstream fout(oufi);
  if(!fout){
    cout << "Cannot open output file." << endl;
    return;
  }

  //calculation
  fout<<"#theta_cm, theta3_lab, Ek3_lab, theta4_lab, Ek4_lab"<< "\n";
  for(int i=0;i<n+1;i++){
    if(i%100==0){
      cout<<"Finish "<<100*i/n<<" % !      \r"<<flush;
    }

    Double_t theta_cm_rad = i * dtheta;
    Double_t theta_cm_deg = theta_cm_rad * 180.0 / TMath::Pi();
    Double_t ct = cos(theta_cm_rad);
    Double_t st = sin(theta_cm_rad);
    //    cout<<"i:"<<i<<", theta:"<<theta_cm/TMath::Pi()<<endl;

    //particle 3
    Double_t E3_lab = gamma*(E3_cm + beta*p_cm*ct);
    Double_t p3_para = gamma*(p_cm*ct + beta*E3_cm);
    Double_t p3_perp = p_cm*st;
    //    Double_t theta3_lab = atan2(p3_perp, p3_para);
    Double_t theta3_lab = atan2(p3_perp, p3_para);
    Double_t theta3_deg = theta3_lab * 180.0 / TMath::Pi();
    //    if(theta3_deg < 0) theta3_deg += 180.0;

    Double_t Ek3_lab = E3_lab - m3p;
    if(Ek3_lab < 0) Ek3_lab = 0;

    //particle 4  (cm angle to pi-theta_cm)
    Double_t ct4 = -ct;
    Double_t st4 =  st;
    Double_t E4_lab = gamma*(E4_cm + beta*p_cm*ct4);
    Double_t p4_para = gamma*(p_cm*ct4 + beta*E4_cm);
    Double_t p4_perp = p_cm*st4;
    //    Double_t theta4_lab = atan2(p4_perp,p4_para);
    Double_t theta4_lab = atan2(p4_perp, p4_para);
    Double_t theta4_deg = theta4_lab * 180.0 / TMath::Pi();
    //    if(theta4_deg < 0) theta4_deg += 180.0;

    Double_t Ek4_lab = E4_lab - m4p;
    if(Ek4_lab < 0) Ek4_lab = 0;

    //output (cm angle, theta3, Ek3, theta4, Ek4)
    fout<<theta_cm_deg<<" "<<theta3_deg<<" "<<Ek3_lab<<" "<<theta4_deg<<" "<<Ek4_lab<< "\n";
  }

  fout.close();
  cout << "Finish 100 % !"<< endl;
  cout << "Output saved to: " << oufi << endl;
  gSystem->Exit(0);

}
