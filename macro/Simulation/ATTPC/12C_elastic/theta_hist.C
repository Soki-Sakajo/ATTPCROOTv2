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

int theta_hist(Int_t subnum=1){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set variable
  // want to get
  Int_t eve;
  Double_t theta_lab[2];
  Double_t theta_cms[2];
  Double_t dth;
  Double_t bo_dth=1.8;
  Double_t Ea;
  
  //input
  TChain *tree1=new TChain("anadata","anadata");
    for(int i=0; i<subnum; i++){
      //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
      TString dir="data1";
      TString refi="sim_macro_c12_" + to_string(i) + ".root";
      TString direfi=dir+"/"+refi;
      tree1->Add(direfi);
      if(i==0){cout<<"read from  "<<refi<<flush;}
      if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
    }
  cout<<"end reading files"<<endl<<endl;
  if(subnum==1){
    cout<<endl;
    cout<<"ATTENTION!! read ONLY one file!!"<<endl<<endl<<endl;
  }

  Int_t nEvents=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve);
  tree1->SetBranchAddress("theta_cms",&theta_cms);
  tree1->SetBranchAddress("theta_lab",&theta_lab);
  tree1->SetBranchAddress("TotalE",&Ea);
  //  tree1->SetBranchAddress("startz",&startz);

  //set output
  /*
  TString outfname = "data1/angle_hist_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  //  TTree *PdgCode=new TTree("PdgCode","PdgCode");
  TTree *angdata = new TTree("angdata","angdata");

  // new branch for simmacro.root
  angdata->Branch("event_number",&eve,"eve/I");
  angdata->Branch("theta1_cms",&theta1,"theta1/D");
  angdata->Branch("theta2_cms",&theta2,"theta2/D");
  */

  //  TH2F *h_theta=new TH2F("h_theta","theta:dth",200,0,200,300,1.5,2.1);
  TH2F *h_theta=new TH2F("h_theta","theta:dth",100,0,50,300,1.5,2.1);
  //  TH1F *h_sin=new TH1F("h_sin","ang:number_of_event",200,0,200);
  //  TF1 *f= new TF1("f","[0]*sin([1]*x)",0,180);
  
  //  cout<<"finished set. reading event"<<endl;
  cout<<"border of dth "<<bo_dth<<endl<<endl;
  for(int i=0;i < nEvents; i++){
    tree1->GetEntry(i);
    if(eve%2==1){
      for(int j=0; j<2; j++){
	dth=theta_cms[j]/theta_lab[j];
	//	h_theta->Fill(theta_cms[j],dth);
	h_theta->Fill(Ea,dth);
	if(dth<bo_dth){
	  cout<<"eve "<<eve<<"  par "<<j<<"  dth "<<dth<<endl;
      }
    }
    }
    else{continue;}
  }
  /*
  f->SetParameters(650,deg_to_rad);
  f->SetParLimits(1,deg_to_rad,deg_to_rad);
  h_ang->Fit("f","","",0,180);
  */


  //  h_theta->SetStats(0);
  h_theta->Draw("colz");

  
  /*
  double pa=f->GetParameter(0);
  double da=f->GetParError(0);
  double pb=f->GetParameter(1);
  double db=f->GetParError(1);

  cout<<endl;
  cout<<"set depth "<<z<<" mm"<<endl;
  cout<<"a  "<<pa<<"  da  "<<da<<"  b  "<<pb<<"  db  "<<db<<endl;
  */  
  //  angdata->AutoSave();
  //  simmacro->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
