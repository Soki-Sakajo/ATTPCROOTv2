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
#include"TChain.h"
#include <fstream>
#include <vector>
#include <algorithm>

int check_dEdz(){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t dz=1;           // mm; bin of depth
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  Int_t ctom=10;

  //set variable
  // common
  Int_t div_l=(int)(l_attpc/dz+0.9);
  
  // for 12C12C
  Int_t eve1;
  Int_t n_p1;
  Int_t n_fi1;
  Int_t n_ev1;
  Int_t n_fi1_temp=0;
  Int_t eve1_temp=0;
  Int_t z_n1=0;
  Int_t z_n2=0;
  Double_t E1;
  Double_t z1_max=0;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1;
  Double_t Eave1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t l1[div_l];
  Double_t n1[div_l];
  Double_t Esum1[div_l];

  vector<Double_t> ev_ve1(div_l);
  vector<Double_t> fi_ve1(div_l);

  TString refi="data1/sim_macro_12c_nd_0.root";
  TFile* file =new TFile(refi);
  TTree *tree1 = (TTree*)file->Get("anadata");
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("n",&n_p1);

  vector<vector<Double_t>> E_ve1(Entry1, vector<Double_t>(2,0));
  vector<vector<Double_t>> E_ve2(Entry1, vector<Double_t>(2,0));

  //set histgram
  TH2F *h_ener=new TH2F("h_ener","reaction position:reaction energy",nbin_zn,0,400,500,0,50);
  TH2F *h_dedz=new TH2F("h_dEdz","Energy:dEdz",500,0,50,50,0,50);
  //  TH2F *h_dedz=new TH2F("h_dEdz","Energy:dEdz",500,0,50,10000,0,100);
  
  //set output for making rutdata
  TString outdata = "data1/check_dEdz.root";
  TFile *rutfile=new TFile(outdata, "recreate");
  TTree *rut12c = new TTree("rut12c","rut12c");
  rut12c->Branch("reaction_position",&dep1,"dep1/D");
  rut12c->Branch("energy",&Eave1,"Eave1/D");

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    if(eve1_temp>eve1){n_fi1_temp++;}
    if(z1_max<z1[0]){z1_max=z1[0];}
    if(z1_max<z1[1]){z1_max=z1[1];}
    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0 && theta1<=theta_max){
    if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      z_n1=int(z1[0]/dz);
      z_n2=int(z1[0]);
      E_ve1[j][0]=z_n1;
      E_ve1[j][1]=E1;
      E_ve2[j][0]=E1;
      E_ve2[j][1]=E1/z_n1;
    }
  }
  
  //sort
  sort(E_ve1.begin(),E_ve1.end(),[](const vector<double> &alp, const vector<double> &bet){return alp[1] < bet[1];});
  sort(E_ve2.begin(),E_ve2.end(),[](const vector<double> &gan, const vector<double> &del){return gan[1] < del[1];});

  
  for(int j=0; j< Entry1; j++){
    h_ener->Fill(E_ve1[j][0],E_ve1[j][1]);
    h_dedz->Fill(E_ve2[j][0],E_ve2[j][1]);
    dep1 =E_ve2[j][0];
    Eave1=E_ve2[j][1];
    rut12c->Fill();
    //need
    cout<<"12C12C finish "<<100*j/Entry1<<"% \r"<<flush;
  }
  /*  
  //make error bar
  // for 12c12c
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z->GetBinContent(m));
    h_z->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax->GetBinContent(m));
    h_rmax->SetBinError(m,err);
  }
  */

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(1,2);
  
  //  h_ener->SetStats(0);
  h_ener->SetStats(1);
  can->cd(1);
  h_ener->Draw("colz");

  h_dedz->SetStats(0);
  //  h_dedz->SetStats(1);
  can->cd(2);
  //  gPad->SetLogy(1);
  h_dedz->Draw("colz");
  //  h_dedz->Draw("hist");

  //output for canvas
  TString outcan = "data1/check_can_dEdz.C";
  can->SaveAs(outcan);
  h_ener->Write();
  h_dedz->Write();

  //save branch
  rut12c ->AutoSave();
  rutfile->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "output: " << outdata << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
