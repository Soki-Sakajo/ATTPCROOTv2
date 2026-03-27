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
#include"TChain.h"

int comp_hists(){ 

  Int_t nbin=300;
  Int_t nbin_zn=400;
  
  Int_t  eve1;
  Double_t l1;
  Double_t ang1;
  Double_t rmax1;
  Double_t Y1;
  Double_t E1;

  Int_t  eve2;
  Double_t l2;
  Double_t ang2;
  Double_t rmax2;
  Double_t Y2;
  Double_t E2;

  //input
  TString inf1 = "data1/rutdata_c12_10.root";
  TString inf2 = "data1/rutdata_c12_30.root";  
  TFile *file1 =new TFile(inf1);
  TFile *file2 =new TFile(inf2); 
  TTree *tree1 = (TTree*)file1->Get("rutdata");
  TTree *tree2 = (TTree*)file2->Get("rutdata");

  cout<<"end reading files"<<endl<<endl;

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("reaction_position",&l1);
  tree1->SetBranchAddress("angle",&ang1);
  tree1->SetBranchAddress("rmax_c12",&rmax1);  
  tree1->SetBranchAddress("reaction_rate",&Y1);
  tree1->SetBranchAddress("Eave",&E1);

  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("reaction_position",&l2);
  tree2->SetBranchAddress("angle",&ang2);
  tree2->SetBranchAddress("rmax_c12",&rmax2);  
  tree2->SetBranchAddress("reaction_rate",&Y2);
  tree2->SetBranchAddress("Eave",&E2);

  //  cout<<"get entry"<<endl;
  
  //set histgram
  TH1F *h_rmax1=new TH1F("h_rmax1","rmax:reaction_rate 100",nbin,0,300);
  TH1F *h_rmax2=new TH1F("h_rmax2","rmax:reaction_rate 300",nbin,0,300);
  TH1F *h_z1=new TH1F("h_z1","reaction position:reaction_rate 100",nbin_zn,0,400);
  TH1F *h_z2=new TH1F("h_z2","reaction position:reaction_rate 300",nbin_zn,0,400);  

  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    h_z1->Fill(l1,Y1);
    h_rmax1->Fill(rmax1,Y1);
  }
  cout<<"filled entry1"<<endl;
  for(int j=0; j<Entry2; j++){
    tree2->GetEntry(j);
    h_z2->Fill(l2,Y2);
    h_rmax2->Fill(rmax2,Y2);
  }

  cout<<"filled data"<<endl;
  
  //make error bar
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax1->GetBinContent(m));
    h_rmax1->SetBinError(m,err);
  }

  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax2->GetBinContent(m));
    h_rmax2->SetBinError(m,err);
  }

  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z1->GetBinContent(m));
    h_z1->SetBinError(m,err_z);
  }
  cout<<"fill entry1"<<endl;
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z2->GetBinContent(m));
    h_z2->SetBinError(m,err_z);
  }

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(2,2);
  
  // z-reaction rate
  //  h_z1->SetStats(0);
  //  h_z2->SetStats(0);
  h_z1->SetStats(1);

  can->cd(1);
  gPad->SetLogy(1);
  h_z1->SetMinimum(1);
  h_z1->SetMaximum(10000);  
  h_z1->Draw("hist");

  h_z2->SetStats(1);    
  can->cd(2);
  gPad->SetLogy(1);
  h_z2->SetMinimum(1);
  h_z2->SetMaximum(10000);  
  h_z2->Draw("hist");
  
  // rmax-reaction rate
  //  h_rmax1->SetStats(0);
  //  h_rmax2->SetStats(0);
  h_rmax1->SetStats(1);
  can->cd(3);
  gPad->SetLogy(1);
  h_rmax1->SetMinimum(1e-5);
  h_rmax1->SetMaximum(1e5);  
  h_rmax1->Draw("hist");

  h_rmax2->SetStats(1);
  can->cd(4);
  gPad->SetLogy(1);
  h_rmax2->SetMinimum(1e-5);
  h_rmax2->SetMaximum(1e5);  
  h_rmax2->Draw("hist");

  //output for canvas
  //  TString outcan = "data/ruther_12c_can_vector.C";
  //  TString outcan = "data/ruther_12c_can_" + to_string(subnum) + ".C";
  TString outcan = "data1/comp_hists_can_10_30.C";
  can->SaveAs(outcan);
  /*
  h_z1->Write();
  h_z2->Write();
  h_rmax1->Write();
  h_rmax2->Write();
  */
  
  cout << endl;
  cout << "output: "<< outcan << endl<<endl;

  return 0;
}
