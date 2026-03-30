#include <iostream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TEventList.h>
#include <TH1F.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include"TString.h"
#include"TFile.h"
#include"TTree.h"
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TClonesArray.h"
#include"TCanvas.h"
#include"TMath.h"
#include<fstream>
#include<vector>

int test_eventlist(){
  
  Int_t z_s=6;          // Z of scatterd particl in simfile
  Int_t z_r1=6;         // Z of an recoil particl in simfile
  Int_t z_r2=1;         // Z of another recoil particl in simfile
  Int_t dz=1;           // mm; bin of depth
  Double_t da=1;        // deg; bin of angle
  Double_t dr=1;        // mm; bin of reaction rate
  Int_t nbin=300/dr;    // bin of rmax-reaction_rate hist
  Int_t nbin_zn=400/dz; // bin of z-reaction_rate hist
  Int_t theta_min=3;    // deg; rutherford scattar angle
  Double_t dtheta=1;    // dig; bin of angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  int ctom=10;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t da_r=da*deg_to_rad;
  
  //set variable
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Int_t eve;
  Int_t n_data;
  Int_t n_p1;
  Double_t rmax1;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t theta1_cms_r;
  Double_t Y1_rmax=0;
  Double_t l11[div_l];
  Double_t l12[div_l];
  Double_t ang11[div_a];
  Double_t ang_r11[div_a];
  Double_t ang12[div_a];
  Double_t ang_r12[div_a];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  Double_t n_react1[div_l];
  //  Double_t n_beam1[div_l];
  Int_t n_rmax1=0;

  Double_t z[div_r];
  Double_t z1_max;
  //  vector<vector<Double_t>> elist(div_l, vector<Double_t>(div_a,0));
  
  //input
  TString infname = "data/sim_macro_c12.root";
  TFile *file1 =new TFile(infname); 
  TTree *tree1 = (TTree*)file1->Get("anadata");

  TEventList *elist[div_l][div_a];
  
  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);
  Double_t dep1;
  Double_t dir1;
  Double_t theta1[Entry1];

  vector<vector<vector<Int_t>>> n_ve(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Int_t>>> np_ve(div_l, vector<vector<Int_t>>(div_a,vector<Int_t>(0)));

  bool eol=true;
  int z_n=0;
  int a_n=0;
  
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    if(z1_max<z1[0]){z1_max=z1[0];}
    if(z1_max<z1[1]){z1_max=z1[1];}
    if(theta1_cms[0]<=theta1_cms[1]){theta1[j]=theta1_cms[0];}
    else if(theta1_cms[0]>theta1_cms[1]){theta1[j]=theta1_cms[1];}

    if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      z_n=int(z1[0]/dz);
      a_n=int(theta1[j]/da);
      n_ve[z_n][a_n].push_back(j);
      E_ve[z_n][a_n].push_back(E1);
      r_ve[z_n][a_n].push_back(rmax1);
      np_ve[z_n][a_n].push_back(n_p1);

      cout<<j<<"  "<<n_ve.at(z_n).at(a_n).at(0)<<"   "<<E_ve.at(z_n).at(a_n).at(0)<<endl;
      //	cout<<n_ve[int(z1[0]/dz)][int(theta1[j]/da)]<<endl;
      
    }
  }
  
  return 0;
}
