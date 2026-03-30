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

int comp_detection_0825_24(){ 

  Int_t nbin=300;
  //  Int_t nbin_zn=450;
  Int_t nbin_zn=1000;
  Int_t ctom=10;
  //  Double_t cutrmax=270;
  Double_t cutrmax=250;
  
  Int_t  n1;
  Int_t  n_r1=0;
  Int_t  eve1;
  Double_t l1;
  Double_t rmax1;
  Double_t rmax1_temp=0;
  Double_t rmax1_max=0;
  Double_t E1[7];

  Int_t  n_r2=0;
  Int_t  eve2;
  Double_t l2;
  Double_t rmax2;
  Double_t E2[7];

  Int_t  n_r3=0;
  Int_t  eve3;
  Double_t l3;
  Double_t rmax3;
  Double_t E3[7];

  Int_t  n_r4=0;
  Int_t  eve4;
  Double_t l4;
  Double_t rmax4;
  Double_t rmax4_temp=0;
  Double_t rmax4_max=0;
  Double_t E4[7];

  //input
  TString inf1 = "data1/sim_macro_10matm_0.root";
  TString inf2 = "data1/sim_macro_20matm_0.root";
  TString inf3 = "data1/sim_macro_30matm_0.root";
  //  TString inf4 = "data1/sim_macro_100matm_0.root";
  //  TString inf4 = "data1/sim_macro_130matm_0.root";
  TFile *file1 =new TFile(inf1);
  TFile *file2 =new TFile(inf2);
  TFile *file3 =new TFile(inf3);
  //  TFile *file4 =new TFile(inf4);
  TTree *tree1 = (TTree*)file1->Get("anadata");
  TTree *tree2 = (TTree*)file2->Get("anadata");
  TTree *tree3 = (TTree*)file3->Get("anadata");
  //inf4
  //  TTree *tree4 = (TTree*)file4->Get("anadata");

  cout<<"end reading files"<<endl;

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("n",&n1);
  tree1->SetBranchAddress("stopz",&l1);
  tree1->SetBranchAddress("rmax",&rmax1);  
  tree1->SetBranchAddress("energy",&E1);

  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("stopz",&l2);
  tree2->SetBranchAddress("rmax",&rmax2);  
  tree2->SetBranchAddress("energy",&E2);

  Int_t Entry3=tree3->GetEntries();
  tree3->SetBranchAddress("event_number",&eve3);
  tree3->SetBranchAddress("stopz",&l3);
  tree3->SetBranchAddress("rmax",&rmax3);  
  tree3->SetBranchAddress("energy",&E3);

  //inf4
  /*
  Int_t Entry4=tree4->GetEntries();
  tree4->SetBranchAddress("event_number",&eve4);
  tree4->SetBranchAddress("stopz",&l4);
  tree4->SetBranchAddress("rmax",&rmax4);  
  tree4->SetBranchAddress("energy",&E4);
  */
  //inf4
  cout<<"get entry"<<endl;
  
  //set histgram
  TH2F *h_rmax1=new TH2F("h_rmax1","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_rmax2=new TH2F("h_rmax2","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_rmax3=new TH2F("h_rmax3","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  //inf4
  //  TH2F *h_rmax4=new TH2F("h_rmax4","reaction position: rmax",nbin_zn,0,nbin_zn,nbin,0,300);
  //  TH2F *h_rmax3=new TH2F("h_rmax3","reaction position: rmax",250,0,250,nbin,0,300);

  TH2F *h_cut1=new TH2F("h_cut1","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_cut2=new TH2F("h_cut2","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);
  TH2F *h_cut3=new TH2F("h_cut3","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);
  //inf4
  //  TH2F *h_cut4=new TH2F("h_cut4","reaction position: rmax (use eff)",nbin_zn,0,nbin_zn,nbin,0,300);

  TH1F *h_z1=new TH1F("h_z1","reaction position:reaction_rate 100",nbin_zn,0,400);
  TH1F *h_z2=new TH1F("h_z2","reaction position:reaction_rate 300",nbin_zn,0,400);  

  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    //    if(n1==7){
    h_rmax1->Fill(l1*ctom,rmax1);
    if(rmax1_temp<rmax1){rmax1_temp=rmax1;}
    if(rmax1<cutrmax){
      n_r1++;
      h_cut1->Fill(l1*ctom,rmax1);
    }
    rmax1_max=rmax1_temp;
    //    }
  }
  for(int i=0; i<Entry2; i++){
    tree2->GetEntry(i);
    h_rmax2->Fill(l2*ctom,rmax2);
    if(rmax2<cutrmax){
      n_r2++;
      h_cut2->Fill(l2*ctom,rmax2);
    }
  }
  for(int i=0; i<Entry3; i++){
    tree3->GetEntry(i);
    h_rmax3->Fill(l3*ctom,rmax3);
    if(rmax3<cutrmax){
      n_r3++;
      h_cut3->Fill(l3*ctom,rmax3);
    }
  }
  //inf4
  /*
  for(int i=0; i<Entry4; i++){
    tree4->GetEntry(i);
    h_rmax4->Fill(l4*ctom,rmax4);
    if(rmax4_temp<rmax4){rmax4_temp=rmax4;}
    if(rmax4<cutrmax){
      n_r4++;
      h_cut4->Fill(l4*ctom,rmax4);
    }
    rmax4_max=rmax4_temp;
  }
  */
  //inf4
  //calculate detection efficiency
  cout<<endl;
  cout<<"information about detection efficiency (rmax is under cutmax)"<<endl;
  cout<<" cutrmax is "<<cutrmax<<endl;
  cout<<"  10 matm  all event "<<Entry1<<" under cutmax event "<<n_r1<<endl;
  cout<<"            detection efficiency is "<<100*n_r1/Entry1<<" %"<<endl;
  //  cout<<"            rmax1_rmax is "<<rmax1_max<<" mm"<<endl;  
  cout<<"  20 matm  all event "<<Entry2<<" under cutmax event "<<n_r2<<endl;
  cout<<"            detection efficiency is "<<100*n_r2/Entry2<<" %"<<endl;
  cout<<"  30 matm  all event "<<Entry3<<" under cutmax event "<<n_r3<<endl;
  cout<<"            detection efficiency is "<<100*n_r3/Entry3<<" %"<<endl;
  //inf4
  /*
  cout<<"  130 matm  all event "<<Entry4<<" under cutmax event "<<n_r4<<endl;
  cout<<"            detection efficiency is "<<100*n_r4/Entry4<<" %"<<endl;
  //  cout<<"            rmax4_rmax is "<<rmax4_max<<" mm"<<endl;  
  */
  //inf4
  /*
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
  */

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,1);
  //inf4
  //  can->Divide(4,1);

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
  //inf4
  /*
  h_rmax4->SetStats(0);
  //  h_rmax3->SetStats(1);
  can->cd(4);
  h_rmax4->Draw("colz");
  */
  //inf4
  /*
  h_cut1->SetStats(1);
  can->cd(4);
  h_cut1->Draw("colz");

  h_cut2->SetStats(1);
  can->cd(5);
  h_cut2->Draw("colz");

  h_cut3->SetStats(1);
  can->cd(6);
  h_cut3->Draw("colz");
  */
  
  //output for canvas
  TString outcan = "data1/comp_detection_0825.C";
  can->SaveAs(outcan);
  
  cout << endl;
  cout << "output: "<< outcan << endl<<endl;

  return 0;
}
