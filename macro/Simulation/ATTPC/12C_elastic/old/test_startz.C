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

//int ruther_macro_c12_each_angle_0126_24(){
int test_startz(){ 
  
  //This file is only 12C12C and 12C1H elastic scattering
  
  Double_t startz[2];
  Double_t z[2];

  
  //input
  //  TString infname = "data/sim_macro_c12.root";
  TString infname = "data/sim_macro_c12_0.root";
  //  TString infname = "data1/sim_macro_c12_" + to_string(subnum) + ".root";
  TFile *file1 =new TFile(infname); 
  TTree *tree1 = (TTree*)file1->Get("anadata");

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("startz",&startz);
  Double_t dz[Entry1];
  
  for(int j=0; j< Entry1; j++){
    tree1->SetBranchStatus("*",0);
    tree1->SetBranchStatus("startz",1);
    z[0]=startz[0]*10;
    z[1]=startz[1]*10;

    dz[j]=fabs(z[0]-z[1]);
      cout<<j<<"    "<<dz[j]<<endl;    
    if(j%2==1 && 1<fabs(z[0]-z[1])){

    }      
  }
  return 0;
}
