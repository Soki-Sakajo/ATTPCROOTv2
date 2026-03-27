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

int comp_detection_direct_0202_25(){

  Int_t nbin=300;
  //  Int_t nbin_zn=450;
  Int_t nbin_zn=350;
  Int_t ctom=10;
  Double_t lowlimt=0;
  Double_t cutrmax=250;
  
  Int_t pre1=10;
  Int_t pre2=20;
  Int_t pre3=30;
  Int_t pre4=40;
  Int_t pre5=50;
  Int_t pre6=60;
  Int_t pre7=70;
  Int_t pre8=80;
  Int_t pre9=90;
  Int_t pre10=100;
  Int_t pre11=110;
  Int_t pre12=120;
  Int_t pre13=130;
  Int_t  n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13;
  Int_t  n_r1,n_r2,n_r3,n_r4,n_r5,n_r6,n_r7,n_r8,n_r9,n_r10,n_r11,n_r12,n_r13;
  Int_t  n_eve1,n_eve2,n_eve3,n_eve4,n_eve5,n_eve6,n_eve7,n_eve8,n_eve9,n_eve10,n_eve11,n_eve12,n_eve13;
  Int_t  eve1,eve2,eve3,eve4,eve5,eve6,eve7,eve8,eve9,eve10,eve11,eve12,eve13;
  Double_t rmax1,rmax2,rmax3,rmax4,rmax5,rmax6,rmax7,rmax8,rmax9,rmax10,rmax11,rmax12,rmax13;
  Double_t detef1,detef2,detef3,detef4,detef5,detef6,detef7,detef8,detef9,detef10,detef11,detef12,detef13;
  //  Double_t E1[7];
  //  Double_t startz1[6];
  //  Double_t stopz1[6];
  //  Double_t l1;

  //input
  TString inf1 = "data1/sim_macro_"+to_string(pre1)+"matm_0.root";
  TFile *file1 =new TFile(inf1);
  TTree *tree1 = (TTree*)file1->Get("anadata");
  
  TString inf2 = "data1/sim_macro_"+to_string(pre2)+"matm_0.root";
  TFile *file2 =new TFile(inf2);
  TTree *tree2 = (TTree*)file2->Get("anadata");
  
  TString inf3 = "data1/sim_macro_"+to_string(pre3)+"matm_0.root";
  TFile *file3 =new TFile(inf3);
  TTree *tree3 = (TTree*)file3->Get("anadata");
  
  TString inf4 = "data1/sim_macro_"+to_string(pre4)+"matm_0.root";
  TFile *file4 =new TFile(inf4);
  TTree *tree4 = (TTree*)file4->Get("anadata");
  
  TString inf5 = "data1/sim_macro_"+to_string(pre5)+"matm_0.root";
  TFile *file5 =new TFile(inf5);
  TTree *tree5 = (TTree*)file5->Get("anadata");
  
  TString inf6 = "data1/sim_macro_"+to_string(pre6)+"matm_0.root";
  TFile *file6 =new TFile(inf6);
  TTree *tree6 = (TTree*)file6->Get("anadata");
  
  TString inf7 = "data1/sim_macro_"+to_string(pre7)+"matm_0.root";
  TFile *file7 =new TFile(inf7);
  TTree *tree7 = (TTree*)file7->Get("anadata");
  
  TString inf8 = "data1/sim_macro_"+to_string(pre8)+"matm_0.root";
  TFile *file8 =new TFile(inf8);
  TTree *tree8 = (TTree*)file8->Get("anadata");
  
  TString inf9 = "data1/sim_macro_"+to_string(pre9)+"matm_0.root";
  TFile *file9 =new TFile(inf9);
  TTree *tree9 = (TTree*)file9->Get("anadata");
  
  TString inf10 = "data1/sim_macro_"+to_string(pre10)+"matm_0.root";
  TFile *file10 =new TFile(inf10);
  TTree *tree10 = (TTree*)file10->Get("anadata");
  
  TString inf11 = "data1/sim_macro_"+to_string(pre11)+"matm_0.root";
  TFile *file11 =new TFile(inf11);
  TTree *tree11 = (TTree*)file11->Get("anadata");
  
  TString inf12 = "data1/sim_macro_"+to_string(pre12)+"matm_0.root";
  TFile *file12 =new TFile(inf12);
  TTree *tree12 = (TTree*)file12->Get("anadata");
  
  TString inf13 = "data1/sim_macro_"+to_string(pre13)+"matm_0.root";
  TFile *file13 =new TFile(inf13);
  TTree *tree13 = (TTree*)file13->Get("anadata");
  cout<<"end setting files"<<endl;

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("n",&n1);
  tree1->SetBranchAddress("rmax",&rmax1);
  //  tree1->SetBranchAddress("startz",&startz1);
  //  tree1->SetBranchAddress("stopz",&stopz1);

  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("n",&n2);
  tree2->SetBranchAddress("rmax",&rmax2);

  Int_t Entry3=tree3->GetEntries();
  tree3->SetBranchAddress("event_number",&eve3);
  tree3->SetBranchAddress("n",&n3);
  tree3->SetBranchAddress("rmax",&rmax3);

  Int_t Entry4=tree4->GetEntries();
  tree4->SetBranchAddress("event_number",&eve4);
  tree4->SetBranchAddress("n",&n4);
  tree4->SetBranchAddress("rmax",&rmax4);

  Int_t Entry5=tree5->GetEntries();
  tree5->SetBranchAddress("event_number",&eve5);
  tree5->SetBranchAddress("n",&n5);
  tree5->SetBranchAddress("rmax",&rmax5);

  Int_t Entry6=tree6->GetEntries();
  tree6->SetBranchAddress("event_number",&eve6);
  tree6->SetBranchAddress("n",&n6);
  tree6->SetBranchAddress("rmax",&rmax6);

  Int_t Entry7=tree7->GetEntries();
  tree7->SetBranchAddress("event_number",&eve7);
  tree7->SetBranchAddress("n",&n7);
  tree7->SetBranchAddress("rmax",&rmax7);

  Int_t Entry8=tree8->GetEntries();
  tree8->SetBranchAddress("event_number",&eve8);
  tree8->SetBranchAddress("n",&n8);
  tree8->SetBranchAddress("rmax",&rmax8);

  Int_t Entry9=tree9->GetEntries();
  tree9->SetBranchAddress("event_number",&eve9);
  tree9->SetBranchAddress("n",&n9);
  tree9->SetBranchAddress("rmax",&rmax9);

  Int_t Entry10=tree10->GetEntries();
  tree10->SetBranchAddress("event_number",&eve10);
  tree10->SetBranchAddress("n",&n10);
  tree10->SetBranchAddress("rmax",&rmax10);

  Int_t Entry11=tree11->GetEntries();
  tree11->SetBranchAddress("event_number",&eve11);
  tree11->SetBranchAddress("n",&n11);
  tree11->SetBranchAddress("rmax",&rmax11);

  Int_t Entry12=tree12->GetEntries();
  tree12->SetBranchAddress("event_number",&eve12);
  tree12->SetBranchAddress("n",&n12);
  tree12->SetBranchAddress("rmax",&rmax12);

  Int_t Entry13=tree13->GetEntries();
  tree13->SetBranchAddress("event_number",&eve13);
  tree13->SetBranchAddress("n",&n13);
  tree13->SetBranchAddress("rmax",&rmax13);
  cout<<"get entry"<<endl;
  
  /*
  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    if(n1==6){
      n_eve1++;
      //      l1=startz1[0];
      //      h_rmax1->Fill(l1,rmax1);
      if(lowlimt<rmax1 && rmax1<cutrmax){
	n_r1++;
	//	h_cut1->Fill(l1,rmax1);
      }
    }
    else if(n1==1){
      //      l1=stopz1[0];
      //      h_rmax1->Fill(l1,rmax1);
    }
  }
  */

  n_r1=0;
  n_eve1=0;
  for(int i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    if(n1==6){
      n_eve1++;
      if(lowlimt<rmax1 && rmax1<cutrmax){
	n_r1++;
      }
    }
  }

  n_r2=0;
  n_eve2=0;
  for(int i=0; i<Entry2; i++){
    tree2->GetEntry(i);
    if(n2==6){
      n_eve2++;
      if(lowlimt<rmax2 && rmax2<cutrmax){
	n_r2++;
      }
    }
  }

  n_r3=0;
  n_eve3=0;
  for(int i=0; i<Entry3; i++){
    tree3->GetEntry(i);
    if(n3==6){
      n_eve3++;
      if(lowlimt<rmax3 && rmax3<cutrmax){
	n_r3++;
      }
    }
  }

  n_r4=0;
  n_eve4=0;
  for(int i=0; i<Entry4; i++){
    tree4->GetEntry(i);
    if(n4==6){
      n_eve4++;
      if(lowlimt<rmax4 && rmax4<cutrmax){
	n_r4++;
      }
    }
  }

  n_r5=0;
  n_eve5=0;
  for(int i=0; i<Entry5; i++){
    tree5->GetEntry(i);
    if(n5==6){
      n_eve5++;
      if(lowlimt<rmax5 && rmax5<cutrmax){
	n_r5++;
      }
    }
  }

  n_r6=0;
  n_eve6=0;
  for(int i=0; i<Entry6; i++){
    tree6->GetEntry(i);
    if(n6==6){
      n_eve6++;
      if(lowlimt<rmax6 && rmax6<cutrmax){
	n_r6++;
      }
    }
  }

  n_r7=0;
  n_eve7=0;
  for(int i=0; i<Entry7; i++){
    tree7->GetEntry(i);
    if(n7==6){
      n_eve7++;
      if(lowlimt<rmax7 && rmax7<cutrmax){
	n_r7++;
      }
    }
  }

  n_r8=0;
  n_eve8=0;
  for(int i=0; i<Entry8; i++){
    tree8->GetEntry(i);
    if(n8==6){
      n_eve8++;
      if(lowlimt<rmax8 && rmax8<cutrmax){
	n_r8++;
      }
    }
  }

  n_r9=0;
  n_eve9=0;
  for(int i=0; i<Entry9; i++){
    tree9->GetEntry(i);
    if(n9==6){
      n_eve9++;
      if(lowlimt<rmax9 && rmax9<cutrmax){
	n_r9++;
      }
    }
  }

  n_r10=0;
  n_eve10=0;
  for(int i=0; i<Entry10; i++){
    tree10->GetEntry(i);
    if(n10==6){
      n_eve10++;
      if(lowlimt<rmax10 && rmax10<cutrmax){
	n_r10++;
      }
    }
  }

  n_r11=0;
  n_eve11=0;
  for(int i=0; i<Entry11; i++){
    tree11->GetEntry(i);
    if(n11==6){
      n_eve11++;
      if(lowlimt<rmax11 && rmax11<cutrmax){
	n_r11++;
      }
    }
  }

  n_r12=0;
  n_eve12=0;
  for(int i=0; i<Entry12; i++){
    tree12->GetEntry(i);
    if(n12==6){
      n_eve12++;
      if(lowlimt<rmax12 && rmax12<cutrmax){
	n_r12++;
      }
    }
  }

  n_r13=0;
  n_eve13=0;
  for(int i=0; i<Entry13; i++){
    tree13->GetEntry(i);
    if(n13==6){
      n_eve13++;
      if(lowlimt<rmax13 && rmax13<cutrmax){
	n_r13++;
      }
    }
  }

  //calculate detection efficiency
  detef1=100*n_r1/n_eve1;
  detef2=100*n_r2/n_eve2;
  detef3=100*n_r3/n_eve3;
  detef4=100*n_r4/n_eve4;
  detef5=100*n_r5/n_eve5;
  detef6=100*n_r6/n_eve6;
  detef7=100*n_r7/n_eve7;
  detef8=100*n_r8/n_eve8;
  detef9=100*n_r9/n_eve9;
  detef10=100*n_r10/n_eve10;
  detef11=100*n_r11/n_eve11;
  detef12=100*n_r12/n_eve12;
  detef13=100*n_r13/n_eve13;
  cout<<endl;
  cout<<"information about detection efficiency (rmax is under cutmax)"<<endl;
  cout<<" cutrmax is "<<cutrmax<<endl;
  cout<<"  "<<pre1<<" matm  all event "<<Entry1<<", decay event "<<n_eve1<<", under cutmax event "<<n_r1<<endl;
  cout<<"            detection efficiency is "<<detef1<<" %"<<endl<<endl;
  cout<<"  "<<pre2<<" matm  all event "<<Entry2<<", decay event "<<n_eve2<<", under cutmax event "<<n_r2<<endl;
  cout<<"            detection efficiency is "<<detef2<<" %"<<endl<<endl;
  cout<<"  "<<pre3<<" matm  all event "<<Entry3<<", decay event "<<n_eve3<<", under cutmax event "<<n_r3<<endl;
  cout<<"            detection efficiency is "<<detef3<<" %"<<endl<<endl;
  cout<<"  "<<pre4<<" matm  all event "<<Entry4<<", decay event "<<n_eve4<<", under cutmax event "<<n_r4<<endl;
  cout<<"            detection efficiency is "<<detef4<<" %"<<endl<<endl;
  cout<<"  "<<pre5<<" matm  all event "<<Entry5<<", decay event "<<n_eve5<<", under cutmax event "<<n_r5<<endl;
  cout<<"            detection efficiency is "<<detef5<<" %"<<endl<<endl;
  cout<<"  "<<pre6<<" matm  all event "<<Entry6<<", decay event "<<n_eve6<<", under cutmax event "<<n_r6<<endl;
  cout<<"            detection efficiency is "<<detef6<<" %"<<endl<<endl;
  cout<<"  "<<pre7<<" matm  all event "<<Entry7<<", decay event "<<n_eve7<<", under cutmax event "<<n_r7<<endl;
  cout<<"            detection efficiency is "<<detef7<<" %"<<endl<<endl;
  cout<<"  "<<pre8<<" matm  all event "<<Entry8<<", decay event "<<n_eve8<<", under cutmax event "<<n_r8<<endl;
  cout<<"            detection efficiency is "<<detef8<<" %"<<endl<<endl;
  cout<<"  "<<pre9<<" matm  all event "<<Entry9<<", decay event "<<n_eve9<<", under cutmax event "<<n_r9<<endl;
  cout<<"            detection efficiency is "<<detef9<<" %"<<endl<<endl;
  cout<<"  "<<pre10<<" matm  all event "<<Entry10<<", decay event "<<n_eve10<<", under cutmax event "<<n_r10<<endl;
  cout<<"            detection efficiency is "<<detef10<<" %"<<endl<<endl;
  cout<<"  "<<pre11<<" matm  all event "<<Entry11<<", decay event "<<n_eve11<<", under cutmax event "<<n_r11<<endl;
  cout<<"            detection efficiency is "<<detef11<<" %"<<endl<<endl;
  cout<<"  "<<pre12<<" matm  all event "<<Entry12<<", decay event "<<n_eve12<<", under cutmax event "<<n_r12<<endl;
  cout<<"            detection efficiency is "<<detef12<<" %"<<endl<<endl;
  cout<<"  "<<pre13<<" matm  all event "<<Entry13<<", decay event "<<n_eve13<<", under cutmax event "<<n_r13<<endl;
  cout<<"            detection efficiency is "<<detef13<<" %"<<endl<<endl;

  return 0;
}
