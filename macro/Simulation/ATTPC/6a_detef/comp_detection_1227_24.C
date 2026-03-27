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

int comp_detection_1227(Int_t pre1=0,Int_t pre2=0,Int_t pre3=0){

  if(pre1==0||pre2==0||pre3==0){
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
  Int_t  eve1;
  Double_t startz1[6];
  Double_t stopz1[6];
  Double_t l1;
  Double_t rmax1;
  Double_t lowlimt1;
  Double_t E1[7];

  Int_t  n2;
  Int_t  n_r2=0;
  Int_t  eve2;
  Double_t startz2[6];
  Double_t stopz2[6];
  Double_t l2;
  Double_t rmax2;
  Double_t lowlimt2;
  Double_t E2[6];

  Int_t  n3;
  Int_t  n_r3=0;
  Int_t  eve3;
  Double_t startz3[6];
  Double_t stopz3[6];
  Double_t l3;
  Double_t rmax3;
  Double_t lowlimt3;
  Double_t E3[6];

  //input
  TString inf1 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre1)+"matm_0.root";
  TString inf2 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre2)+"matm_0.root";
  TString inf3 = "data1/sim_macro_sequential_hoylehoyle_"+to_string(pre3)+"matm_0.root";
  TFile *file1 =new TFile(inf1);
  TFile *file2 =new TFile(inf2);
  TFile *file3 =new TFile(inf3);
  TTree *tree1 = (TTree*)file1->Get("anadata");
  TTree *tree2 = (TTree*)file2->Get("anadata");
  TTree *tree3 = (TTree*)file3->Get("anadata");
  
  cout<<"end reading files"<<endl;

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("n",&n1);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("stopz",&stopz1);
  tree1->SetBranchAddress("rmax",&rmax1);  
  tree1->SetBranchAddress("energy",&E1);

  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("n",&n2);
  tree2->SetBranchAddress("startz",&startz2);
  tree2->SetBranchAddress("stopz",&stopz2);
  tree2->SetBranchAddress("rmax",&rmax2);  
  tree2->SetBranchAddress("energy",&E2);

  Int_t Entry3=tree3->GetEntries();
  tree3->SetBranchAddress("event_number",&eve3);
  tree3->SetBranchAddress("n",&n3);
  tree3->SetBranchAddress("startz",&startz3);
  tree3->SetBranchAddress("stopz",&stopz3);
  tree3->SetBranchAddress("rmax",&rmax3);  
  tree3->SetBranchAddress("energy",&E3);

  cout<<"get entry"<<endl;
  
  //set histgram
  TH2F *h_rmax1=new TH2F("h_rmax1","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_rmax2=new TH2F("h_rmax2","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_rmax3=new TH2F("h_rmax3","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);

  //  TH2F *h_rmax3=new TH2F("h_rmax3","reaction position: rmax",250,0,250,nbin,0,300);

  TH2F *h_cut1=new TH2F("h_cut1","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_cut2=new TH2F("h_cut2","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_cut3=new TH2F("h_cut3","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);

  TH1F *h_z1=new TH1F("h_z1","reaction position:reaction_rate 100",nbin_zn,0,400);
  TH1F *h_z2=new TH1F("h_z2","reaction position:reaction_rate 300",nbin_zn,0,400);  

  if(pre1==70){lowlimt1=92;}
  else if(pre1==80){lowlimt1=77;}
  else if(pre1==90){lowlimt1=67;}
  else if(pre1==100){lowlimt1=63;}
  else if(pre1==110){lowlimt1=59;}
  else if(pre1==120){lowlimt1=51;}
  else if(pre1==130){lowlimt1=55;}
  else{
    cout<<"pre1 value is not 70, 80, 90,100,110,120,130. end the macro."<<endl;
    gSystem->Exit(0);
  }

  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    if(n1==6){
      l1=startz1[0];
      h_rmax1->Fill(l1,rmax1);
      if(lowlimt1<rmax1 && rmax1<cutrmax){
	n_r1++;
	h_cut1->Fill(l1,rmax1);
      }
    }
    else if(n1==1){
      l1=stopz1[0];
      h_rmax1->Fill(l1,rmax1);
    }
  }

  if(pre2==70){lowlimt2=92;}
  else if(pre2==80){lowlimt2=77;}
  else if(pre2==90){lowlimt2=67;}
  else if(pre2==100){lowlimt2=63;}
  else if(pre2==110){lowlimt2=59;}
  else if(pre2==120){lowlimt2=51;}
  else if(pre2==130){lowlimt2=55;}
  else{
    cout<<"pre2 value is not 70, 80, 90,100,110,120,130. end the macro."<<endl;
    gSystem->Exit(0);
  }

  for(int i=0; i<Entry2; i++){
    tree2->GetEntry(i);
    if(n2==6){
      l2=startz2[0];
      h_rmax2->Fill(l2,rmax2);
      if(lowlimt2<rmax2 && rmax2<cutrmax){
	n_r2++;
	h_cut2->Fill(l2,rmax2);
      }
    }
    else if(n2==1){
      l2=stopz2[0];
      h_rmax2->Fill(l2,rmax2);
    }
  }

  if(pre3==70){lowlimt3=92;}
  else if(pre3==80){lowlimt3=77;}
  else if(pre3==90){lowlimt3=67;}
  else if(pre3==100){lowlimt3=63;}
  else if(pre3==110){lowlimt3=59;}
  else if(pre3==120){lowlimt3=51;}
  else if(pre3==130){lowlimt3=55;}
  else{
    cout<<"pre3 value is not 70, 80, 90,100,110,120,130. end the macro."<<endl;
    gSystem->Exit(0);
  }

  for(int i=0; i<Entry3; i++){
    tree3->GetEntry(i);
    if(n3==6){
      l3=startz3[0];
      h_rmax3->Fill(l3,rmax3);
      if(lowlimt3<rmax3 && rmax3<cutrmax){
	n_r3++;
	h_cut3->Fill(l3,rmax3);
      }
    }
    else if(n3==1){
      l3=stopz3[0];
      h_rmax3->Fill(l3,rmax3);
    }
  }

  //calculate detection efficiency
  cout<<endl;
  cout<<"information about detection efficiency (rmax is under cutmax)"<<endl;
  cout<<" cutrmax is "<<cutrmax<<endl;
  cout<<"  "<<pre1<<" matm  all event "<<Entry1<<" under cutmax event "<<n_r1<<endl;
  cout<<"            detection efficiency is "<<100*n_r1/Entry1<<" %"<<endl;
  cout<<"  "<<pre2<<" matm  all event "<<Entry2<<" under cutmax event "<<n_r2<<endl;
  cout<<"            detection efficiency is "<<100*n_r2/Entry2<<" %"<<endl;
  cout<<"  "<<pre3<<" matm  all event "<<Entry3<<" under cutmax event "<<n_r3<<endl;
  cout<<"            detection efficiency is "<<100*n_r3/Entry3<<" %"<<endl;

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,2);

  h_rmax1->SetStats(0);
  //  h_rmax1->SetStats(1);
  can->cd(1);
  h_rmax1->Draw("colz");

  h_rmax2->SetStats(0);
  //  h_rmax2->SetStats(1);
  can->cd(2);
  h_rmax2->Draw("colz");

  h_rmax3->SetStats(0);
  //  h_rmax3->SetStats(1);
  can->cd(3);
  h_rmax3->Draw("colz");

  h_cut1->SetStats(1);
  can->cd(4);
  h_cut1->Draw("colz");

  h_cut2->SetStats(1);
  can->cd(5);
  h_cut2->Draw("colz");

  h_cut3->SetStats(1);
  can->cd(6);
  h_cut3->Draw("colz");

  
  //output for canvas
  TString outcan =
    "data1/comp_detection_1227_"+to_string(pre1)+"matm_"+to_string(pre2)+"matm_"+to_string(pre3)+"matm.C";
  can->SaveAs(outcan);
  
  cout << endl;
  cout << "output: "<< outcan << endl<<endl;

  return 0;
}
