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

int check_z_x(Int_t subnum=50){
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  const Int_t ctom=10;
  const Int_t dz=1;
  Int_t eve1;
  Int_t n_p1;
  Double_t startx[2];
  Double_t starty[2];
  Double_t startz[2];
  Double_t x[2];
  Double_t y[2];
  Double_t z[2];
  Double_t E1;
  
  //input multi-file by using TChain
  //common set for input
  TString dir="data1";
  
  // Tchain of 12c12c
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
    //      TString refi="sim_macro_c12_" + to_string(i) + ".root";
    TString refi="sim_macro_12c_" + to_string(i) + ".root";
    TString direfi=dir+"/"+refi;
    tree1->Add(direfi);
    if(i==0){cout<<"read from  "<<refi<<flush;}
    if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
  }
  cout<<"end reading files"<<endl<<endl;
  if(subnum==1){
    cout<<endl;
    cout<<"ATTENTION!! read ONLY one 12C12C file!!"<<endl<<endl<<endl;
  }
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("startx",&startx);
  tree1->SetBranchAddress("starty",&starty);
  tree1->SetBranchAddress("startz",&startz);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("n",&n_p1);

  //set histgram
  TH1F *h_x=new TH1F("h_x","reaction position x:number of event",200,-100,100);
  //  TH1F *h_x=new TH1F("h_x","reaction position x:number of event",600,-300,300);
  TH1F *h_y=new TH1F("h_y","reaction position y:number of event",200,-100,100);
  //  TH1F *h_y=new TH1F("h_y","reaction position y:number of event",600,-300,300);
  TH1F *h_z=new TH1F("h_z","reaction position z:number of event",400,0,400);
  TH2F *h_E=new TH2F("h_E","reaction position z:beam energy",400,0,400,50,0,50);
  TH2F *h_z_x=new TH2F("h_z_x","reaction position z:x ",400,0,400,200,-100,100);
  TH2F *h_z_y=new TH2F("h_z_y","reaction position z:y ",400,0,400,200,-100,100);
  
  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    for(int i=0; i<2; i++){
      x[i]=startx[i]*ctom;
      y[i]=starty[i]*ctom;
      z[i]=startz[i]*ctom;
    }
    if(n_p1==2 && fabs(z[0]-z[1])<dz && fabs(y[0]-y[1])<dz && fabs(x[0]-x[1])<dz){
      h_x->Fill(x[0]);
      h_y->Fill(y[0]);
      h_z->Fill(z[0]);
      h_E->Fill(z[0],E1);
      h_z_x->Fill(z[0],x[0]);
      h_z_y->Fill(z[0],y[0]);
    }
    else if(fabs(z[0]-z[1])<dz && fabs(y[0]-y[1])<dz && fabs(x[0]-x[1])<dz){
      //      cout<<"eve  "<<eve1<<endl;
    }
    else if(n_p1==2 && fabs(z[0]-z[1])<dz && fabs(y[0]-y[1])<dz){
      cout<<"eve "<<eve1<<"  x[0] "<<x[0]<<"     x[1] "<<x[1]<<endl;
    }
    else if(n_p1==2 && fabs(z[0]-z[1])<dz && fabs(x[0]-x[1])<dz){
      cout<<"eve "<<eve1<<"  y[0] "<<y[0]<<"     y[1] "<<y[1]<<endl;
    }
    else if(n_p1==2 && fabs(y[0]-y[1])<dz && fabs(x[0]-x[1])<dz){
      cout<<"eve "<<eve1<<"  z[0] "<<z[0]<<"     z[1] "<<z[1]<<endl;
    }
    else{
      cout<<"else "<<endl;
      cout<<"eve "<<eve1<<" number of particle "<<n_p1<<"  x[0] "<<x[0]<<"     x[1] "<<x[1]
	  <<"  y[0] "<<y[0]<<"     y[1] "<<y[1]<<"  z[0] "<<z[0]<<"     z[1] "<<z[1]<<endl;
    }
    if(j%1000==0){
      //need
      cout<<"filling 12c12c data "<<100*j/Entry1 <<" % !    \r"<<flush;
    }
  }
  
  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,2);
  
  //x
  //  h_x->SetStats(0);
  h_x->SetStats(1);
  can->cd(1);
  gPad->SetLogy(1);
  h_x->Draw("hist");

  //y
  //  h_y->SetStats(0);
  h_y->SetStats(1);
  can->cd(2);
  gPad->SetLogy(1);
  h_y->Draw("hist");

  //z
  //  h_z->SetStats(0);
  h_z->SetStats(1);
  can->cd(3);
  h_z->Draw("hist");

  //z-energy
  //  h_z_x->SetStats(0);
  h_E->SetStats(1);
  can->cd(4);
  h_E->Draw("colz");

  //z-x
  //  h_z_x->SetStats(0);
  h_z_x->SetStats(1);
  can->cd(5);
  h_z_x->Draw("colz");

  //z-y
  //  h_z_y->SetStats(0);
  h_z_y->SetStats(1);
  can->cd(6);
  h_z_y->Draw("colz");

  //output for canvas
  TString outcan = "data1/zx_can_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  /*
  h_x->Write();
  h_y->Write();
  h_z->Write();
  h_E->Write();
  h_z_x->Write();
  h_z_y->Write();
  */
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "output: " << outcan << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
