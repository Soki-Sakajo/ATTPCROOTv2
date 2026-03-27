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

int comp_rmax(Int_t subnum=1){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set variable
  // want to get
  Int_t eve;
  Int_t eve2;
  Int_t z1;
  Int_t z2;
  Double_t rmax;
  Double_t rmax_nd;
  Double_t startz[2];
  Double_t stopz[2];
  Double_t drmax;
  //  Double_t dth;
  //  Double_t bo_dth=1.8;
  //  Double_t Ea;
  
  //input
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
    TString dir="data1";
    TString refi="sim_macro_12c_" + to_string(i) + ".root";
    TString direfi=dir+"/"+refi;
    tree1->Add(direfi);
    if(i==0){cout<<"   data read from  "<<refi<<flush;}
    if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
  }
  TChain *tree2=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    TString dir="data1";
    TString refi="sim_macro_12c_nd_" + to_string(i) + ".root";
    TString direfi=dir+"/"+refi;
    tree2->Add(direfi);
    if(i==0){cout<<"nd data read from  "<<refi<<flush;}
    if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
  }
  cout<<"end reading files"<<endl<<endl;
  if(subnum==1){
    cout<<endl;
    cout<<"ATTENTION!! read ONLY one file!!"<<endl<<endl<<endl;
  }
  
  Int_t nEvents1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve);
  tree1->SetBranchAddress("rmax_nd",&rmax_nd);
  tree1->SetBranchAddress("rmax",&rmax);
  //  tree1->SetBranchAddress("TotalE",&Ea);
  tree1->SetBranchAddress("startz",&startz);
  tree1->SetBranchAddress("stopz",&stopz);

  Int_t nEvents2=tree2->GetEntries();
  //  tree2->SetBranchAddress("event_number",&eve2);
  //  tree2->SetBranchAddress("rmax",&rmax_nd);
  
  //  TH2F *h_theta=new TH2F("h_theta","theta:dth",200,0,200,300,1.5,2.1);
  TH2F *h_rmax=new TH2F("h_rmax","rmax_nd:rmax",400,0,200,400,0,200);
  TH2F *h_rmax1=new TH2F("h_rmax1","rmax_nd:rmax",400,0,200,400,0,200);
  TH2F *h_rmax2=new TH2F("h_rmax2","rmax_nd:rmax",400,0,200,400,0,200);
  TH2F *h_rmax3=new TH2F("h_rmax3","rmax_nd:rmax",400,0,200,400,0,200);
  TH2F *h_z=new TH2F("h_z","z:delta-rmax",400,0,400,400,0,200);
  TH2F *h_rs=new TH2F("h_rs","stop position:rmax",410,-10,400,400,0,200);
  TF1 *f= new TF1("f","[0]*x+[1]",0,170);
  
  //  cout<<"finished set. reading event"<<endl;
  //  cout<<"border of dth "<<bo_dth<<endl<<endl;

  if(nEvents1==nEvents2){cout<<"enable to comp"<<endl;}
  else {cout<<"not enable to comp"<<endl;}

  /*
  for(int i=0;i < nEvents1; i++){
    tree1->GetEntry(i);
    if(fabs(startz[0]-startz[1])<0.1){
      z1=int(startz[1]*10);
    }
    else if(startz[1]==0){
      z1=0;
      //	cout<<i<<"  beam event "<<endl;
    }
    else{cout<<i<<"  not startz equal"<<endl;}
    
    if(stopz[0]>=stopz[1]){z2=int(stopz[0]*10);}
    else if(stopz[1]>stopz[0]){z2=int(stopz[1]*10);}
    else{cout<<i<<"  "<<stopz[0]<<"   "<<stopz[1]<<"   something wrong in stop point"<<endl;}

    drmax=rmax-rmax_nd;
    h_z->Fill(z2,drmax);
    h_rs->Fill(z2,rmax);
    
    if(z1!=0){
	h_rmax->Fill(rmax_nd,rmax);
      if(z2==100){
	h_rmax1->Fill(rmax_nd,rmax);
      }
      if(z2==125){
	h_rmax2->Fill(rmax_nd,rmax);
      }
      if(z2==150){
	h_rmax3->Fill(rmax_nd,rmax);
      }
    }
    else if(startz[0]==-100){continue;}
    else{
      cout<<i<<"   "<<z1<<"   "<<startz[0]<<"   "<<z2<<"   "<<stopz[0]<<"   "<<stopz[1]<<endl;
      continue;
    }
  }
  f->SetParameters(1,0);
  //  h_theta->SetStats(0);

  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,2);

  can->cd(1);
  h_rmax->Draw("colz");
  f->Draw("same");
  //  h_rmax->Draw();

  can->cd(2);
  h_z->Draw("cont");

  can->cd(3);
  h_rs->Draw("colz");

  can->cd(4);
  h_rmax1->Draw("colz");
  f->Draw("same");

  can->cd(5);
  h_rmax2->Draw("colz");
  f->Draw("same");

  can->cd(6);
  h_rmax3->Draw("colz");
  f->Draw("same");

  */
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
