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

int check_depth_ener(Int_t subnum=0){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t dz=5;           // mm; bin of depth
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  Int_t ctom=10;
  Int_t div_l=(int)(l_attpc/dz+0.9);
  
  // for 12C12C
  Int_t eve1;
  Int_t n_p;
  Int_t massnum[7];
  Int_t atomnum[7];
  Int_t eve1_temp=0;
  Int_t z_n1=0;
  Double_t E[7];
  Double_t Emax;
  Double_t Etotal;
  Double_t rmax;
  Double_t cmax;
  Double_t cmin;
  Double_t rmax_nd;
  Double_t z1_max=0;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1;
  Double_t stopx[7];
  Double_t stopy[7];
  Double_t stopz[7];
  Double_t startx[7];
  Double_t starty[7];
  Double_t startz[7];
  Double_t theta_lab[7];
  Double_t theta_cms[7];
  Double_t x[7];
  Double_t y[7];
  Double_t z[7];

  vector<Double_t> ev_ve1;
  vector<Double_t> E_ve1;
  vector<Double_t> Emax_ve1;

  bool flag_dx7=false;
  bool flag_dy7=false;
  bool flag_dz7=false;

  //input
  //  TString refi="data/sim_macro_12c_nd_"+ to_string(subnum)+".root";
  //  TString refi="data1/sim_macro_hoylehoyle_nd_"+ to_string(subnum)+".root";
  TString refi="data1/sim_macro_hoylehoyle_"+ to_string(subnum)+".root";
  TFile* file =new TFile(refi);
  TTree *tree1 = (TTree*)file->Get("anadata");
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("startz",&startz);
  tree1->SetBranchAddress("beam_energy",&Etotal);
  tree1->SetBranchAddress("n",&n_p);

  //set histgram
  //  TH2F *h_dep=new TH2F("h_dep","reaction position: total energy of 6 alpha",nbin_zn,0,400,52,0,52);
  TH2F *h_dep=new TH2F("h_dep","reaction position: total energy of 6 alpha",400,0,400,52,0,52);

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    flag_dx7=true;
    flag_dy7=true;
    flag_dz7=true;
    if(n_p==7){
      for(int i=0; i<7; i++){
	z[i]=startz[i]*ctom;
      }
      for(int i=0; i<7; i++){
	if(i<6){
	  if(fabs(z[i]-z[i+1])>dz){
	    flag_dz7=false;
	    cout<<j<<"  "<<i<<"  "<<fabs(z[i]-z[i+1])<<endl;
	  }
	}
	if(i==6){
	  if(fabs(z[i]-z[0])>dz){
	    flag_dz7=false;
	  }
	}
      }
      if(flag_dz7){
	h_dep->Fill(z[0],Etotal);
      }
      if(flag_dz7==false){
	cout<<"this event is starting z position is not same  "<<j<<endl;
      }
    }
    else if(n_p==1){
    }
    else {
      cout<<"number of particle is not 1 or 7  "<<j<<"   "<<n_p<<"   "<<z[0]<<endl;
    }
    //    cout<<100*j/Entry1<<" % finished   \r "<<flush;
  }

  //  h_rmax->SetStats(0);
  h_dep->SetStats(1);
  //  gPad->SetLogy(1);
  h_dep->Draw("colz");
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
