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

int position_angle_hist(Int_t subnum=1){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set variable
  // want to get
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Int_t eve;
  Double_t lab_theta[2];
  Double_t theta_cms[2];
  Double_t theta1;
  Double_t theta2;
  Double_t w;
  Double_t phi;
  Double_t startz[2];
  Double_t z;
  Int_t dz=1;
  Double_t Ea;
  Double_t E[2];
  
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
  tree1->SetBranchAddress("startz",&startz);
  tree1->SetBranchAddress("TotalE",&Ea);
  tree1->SetBranchAddress("Energy",&E);

  /*
  //set output
  TString outfname = "data1/angle_hist_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  //  TTree *PdgCode=new TTree("PdgCode","PdgCode");
  TTree *angdata = new TTree("angdata","angdata");

  // new branch for simmacro.root
  angdata->Branch("event_number",&eve,"eve/I");
  angdata->Branch("theta1_cms",&theta1,"theta1/D");
  angdata->Branch("theta2_cms",&theta2,"theta2/D");
  */
  
  TH2F *h_ang=new TH2F("h_ang","reaction_position:angle",500,0,50,100,0,100);
  
  //  cout<<"finished set. reading event"<<endl;
  for(int i=0;i < nEvents; i++){
    tree1->GetEntry(i);
    //    cout<<i<<endl;
    if (theta_cms[0]>theta_cms[1]) {theta1=theta_cms[1];}
    else if(theta_cms[0]<theta_cms[1]){theta1=theta_cms[0];}
    //    if(fabs(startz[0]-startz[0])<dz && startz[0]>=z && z+1>startz[0]){
    if(fabs(startz[0]-startz[1])<dz){
	eve;
	//	theta1=theta_cms[0];
	//	theta2=theta_cms[1];
	//      cout<<theta1<<"   "<<theta2<<endl;

	h_ang->Fill(Ea,theta1);
	//	h_ang->Fill(E[0],theta1);

	//	angdata->Fill();      
    }
  }
  /*
  f->SetParameters(650,deg_to_rad);
  f->SetParLimits(1,deg_to_rad,deg_to_rad);
  h_ang->Fit("f","","",0,180);
  */
  h_ang->Draw("colz");

  /*
  double pa=f->GetParameter(0);
  double da=f->GetParError(0);
  double pb=f->GetParameter(1);
  double db=f->GetParError(1);
  */

  cout<<endl;
  cout<<"set depth "<<z<<" mm"<<endl;
  //  cout<<"a  "<<pa<<"  da  "<<da<<"  b  "<<pb<<"  db  "<<db<<endl;
  
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
