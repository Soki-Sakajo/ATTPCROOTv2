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

//int comp_detection_0129_25(Int_t pre1=0){
int comp_detection_0129_25(Int_t pre1=80){

  if(pre1==0){
    cout<<"need to choose gass pressure!!!!!"<<endl;
    gSystem->Exit(0);
  }
  
  Int_t nbin=300;
  //  Int_t nbin_zn=450;
  Int_t nbin_zn=350;
  Int_t ctom=10;
  //  Double_t cutrmax=270;
  Double_t cutrmax=250;
  
  Int_t  n1;
  Int_t  n_r1=0;
  Int_t  n_eve1=0;
  Int_t  eve1;
  Double_t startz1[6];
  Double_t stopz1[6];
  Double_t l1;
  Double_t rmax1;
  Double_t lowlimt1=-1;
  Double_t E1[7];
  Double_t detef1;

  //input
  //  TString inf1 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre1)+"matm_0.root";
  //  TString inf1 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre1)+"matm_1.root";
  TString inf1 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre1)+"matm_11.root";
  TFile *file1 =new TFile(inf1);
  TTree *tree1 = (TTree*)file1->Get("anadata");
  
  cout<<"end reading files"<<endl;

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("n",&n1);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("stopz",&stopz1);
  tree1->SetBranchAddress("rmax",&rmax1);  
  tree1->SetBranchAddress("energy",&E1);

  cout<<"get entry"<<endl;
  
  //set histgram
  TH2F *h_rmax1=new TH2F("h_rmax1","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_cut1=new TH2F("h_cut1","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);

  TH1F *h_z1=new TH1F("h_z1","reaction position:reaction_rate 100",nbin_zn,0,400);
  TH1F *h_z2=new TH1F("h_z2","reaction position:reaction_rate 300",nbin_zn,0,400);  

  /*
  if(pre1==70){lowlimt1=92;}
  else if(pre1==80){lowlimt1=77;}
  else if(pre1==90){lowlimt1=67;}
  else if(pre1==100){lowlimt1=63;}
  else if(pre1==110){lowlimt1=59;}
  //  else if(pre1==120){lowlimt1=51;}
  else if(pre1==120){lowlimt1=56;}
  else if(pre1==130){lowlimt1=53;}
  */
  //adapt results of fit
  if(pre1==70){lowlimt1=84;}
  else if(pre1==80){lowlimt1=74;}
  else if(pre1==90){lowlimt1=66;}
  else if(pre1==100){lowlimt1=60;}
  else if(pre1==110){lowlimt1=55;}
  else if(pre1==120){lowlimt1=49;}
  else if(pre1==130){lowlimt1=48;}
  else{
    cout<<"pre1 value is not 70, 80, 90,100,110,120,130. end the macro."<<endl;
    gSystem->Exit(0);
  }

  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    if(n1==6){
      n_eve1++;
      l1=startz1[0];
      h_rmax1->Fill(l1,rmax1);
      if(lowlimt1<rmax1 && rmax1<cutrmax){
	n_r1++;
	h_cut1->Fill(l1,rmax1);
	/*
	if(n_r1<20){
	  cout<<eve1<<endl;
	}
	*/
      }
    }
    else if(n1==1){
      l1=stopz1[0];
      h_rmax1->Fill(l1,rmax1);
    }
  }

  //calculate detection efficiency
  detef1=100*n_r1/n_eve1;
  cout<<endl;
  cout<<"information about detection efficiency (rmax is under cutmax)"<<endl;
  cout<<" cutrmax is "<<cutrmax<<endl;
  cout<<"  "<<pre1<<" matm  all event "<<Entry1<<", decay event "<<n_eve1<<", under cutmax event "<<n_r1<<endl;
  cout<<"            detection efficiency is "<<detef1<<" %"<<endl;

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(2,1);

  h_rmax1->SetStats(0);
  //  h_rmax1->SetStats(1);
  can->cd(1);
  h_rmax1->Draw("colz");

  h_cut1->SetStats(1);
  can->cd(5);
  h_cut1->Draw("colz");

  //output for canvas
  //  TString outcan ="data1/comp_detection_0129_25_"+to_string(pre1)+"matm.C";
  TString outcan ="data1/comp_detection_0129_25_"+to_string(pre1)+"matm_1.C";
  can->SaveAs(outcan);
  
  cout << endl;
  cout << "output: "<< outcan << endl<<endl;

  return 0;
}
