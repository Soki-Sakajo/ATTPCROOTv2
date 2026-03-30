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

int check_sim_macro(Int_t subnum=0){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t dz=1;           // mm; bin of depth
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
  //  TString outdir="./data/";
  TString outdir="./data1/";
  //  TString refi="data/sim_macro_12c_nd_" + to_string(subnum) + ".root";
  //  TString refi="data1/sim_macro_hoylehoyle_nd_" + to_string(subnum) + ".root";
  //  TString refi="data1/demo_sim_macro_hoyle_nd_" + to_string(subnum) + ".root";
  //  TString refi="data/sim_macro_sequential_hoylehoyle_nd_" + to_string(subnum) + ".root";
  //  TString refi="data1/demo_sim_macro_hoyle_nd_" + to_string(subnum) + ".root";
  TString refi=outdir + "sim_macro_sequential_hoylehoyle_nd_" + to_string(subnum) + ".root";

  TFile* file =new TFile(refi);
  TTree *tree1 = (TTree*)file->Get("anadata");
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("startx",&startx);
  tree1->SetBranchAddress("starty",&starty);
  tree1->SetBranchAddress("startz",&startz);
  tree1->SetBranchAddress("stopx",&stopx);
  tree1->SetBranchAddress("stopy",&stopy);
  tree1->SetBranchAddress("stopz",&stopz);
  tree1->SetBranchAddress("massnum",&massnum);
  tree1->SetBranchAddress("atomnum",&atomnum);
  tree1->SetBranchAddress("energy",&E);
  tree1->SetBranchAddress("energy_max",&Emax);
  tree1->SetBranchAddress("beam_energy",&Etotal);
  tree1->SetBranchAddress("n",&n_p);
  tree1->SetBranchAddress("rmax_nd",&rmax_nd);
  tree1->SetBranchAddress("rmax",&rmax);
  tree1->SetBranchAddress("cmax",&cmax);
  tree1->SetBranchAddress("cmin",&cmin);

  //set histgram
  TH1F *h_x=new TH1F("h_x","reaction position x",100,-5,5);
  TH1F *h_y=new TH1F("h_y","reaction position y",100,-5,5);
  TH1F *h_z=new TH1F("h_z","reaction position z",400,0,400);
  TH1F *h_n=new TH1F("h_n","number of particle",20,0,10);
  TH1F *h_mass=new TH1F("h_mass","mass",30,0,15);
  TH1F *h_atom=new TH1F("h_atom","atom",20,0,10);
  TH2F *h_posi=new TH2F("h_posi","reaction position:energy of alpha",nbin_zn,0,400,350,0,35);
  TH2F *h_cmax=new TH2F("h_cmax","reaction position: max of charge",nbin_zn,0,400,250,0,2.5e6);
  TH2F *h_cmin=new TH2F("h_cmin","reaction position: min of charge in scat",nbin_zn,0,400,150,0,1.5e5);
  TH2F *h_cminb=new TH2F("h_cminb","reaction position: charge in beam",nbin_zn,0,400,300,0,3e6);
  TH2F *h_rmax=new TH2F("h_rmax","reaction position: rmax",nbin_zn,0,400,160,0,160);
  TH2F *h_rmaxnd=new TH2F("h_rmaxnd","reaction position: rmax_nd",nbin_zn,0,400,160,0,160);

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    flag_dx7=true;
    flag_dy7=true;
    flag_dz7=true;
    h_n->Fill(n_p);
    //    if(n_p==7){
    //      for(int i=0; i<7; i++){
    //7->6
    if(n_p==6){
      for(int i=0; i<6; i++){
	/*	
	x[i]=startx[i]*ctom;
	y[i]=starty[i]*ctom;
	z[i]=startz[i]*ctom;
	*/
	x[i]=startx[i];
	y[i]=starty[i];
	z[i]=startz[i];

	h_x->Fill(x[i]);
	h_y->Fill(y[i]);
	h_z->Fill(z[i]);
	h_mass->Fill(massnum[i]);
	h_atom->Fill(atomnum[i]);
      }
      //      for(int i=0; i<7; i++){
	//	if(i<6){
      //7->6
      for(int i=0; i<6; i++){
	if(i<5){
	  if(fabs(x[i]-x[i+1])>1){
	    flag_dx7=false;
	    cout<<j<<"  "<<i<<"  "<<fabs(x[i]-x[i+1])<<endl;
	  }
	  if(fabs(y[i]-y[i+1])>1){
	    flag_dy7=false;
	    cout<<j<<"  "<<i<<"  "<<fabs(y[i]-y[i+1])<<endl;
	  }
	  if(fabs(z[i]-z[i+1])>dz){
	    flag_dz7=false;
	    cout<<j<<"  "<<i<<"  "<<fabs(z[i]-z[i+1])<<endl;
	  }
	}
	//	if(i==6){
	//7->6
	if(i==5){
	  if(fabs(x[i]-x[0])>1){
	    flag_dx7=false;
	  }
	  if(fabs(y[i]-y[0])>1){
	    flag_dy7=false;
	  }
	  if(fabs(z[i]-z[0])>dz){
	    flag_dz7=false;
	  }
	}
      }
      if(flag_dz7){
	h_cmax->Fill(z[0],cmax);
	h_cmin->Fill(z[0],cmin);
	h_rmax->Fill(z[0],rmax);
	h_rmaxnd->Fill(z[0],rmax_nd);
	//	for(int i=0; i<7; i++){
	//7->6
	for(int i=0; i<6; i++){
	  h_posi->Fill(z[i],E[i]);
	}
      }
      if(flag_dx7==false || flag_dy7==false){
	cout<<"this event is starting x or y position is not same  "<<j<<endl;
      }
      if(flag_dz7==false){
	cout<<"this event is starting z position is not same  "<<j<<endl;
      }
    }
    else if(n_p==1){
      //      x[0]=startx[0]*ctom;
      //      y[0]=starty[0]*ctom;
      //      z[0]=stopz[0]*ctom;
      z[0]=stopz[0];
      /*
      h_x->Fill(x[0]);
      h_y->Fill(y[0]);
      h_z->Fill(z[0]);
      h_mass->Fill(massnum[0]);
      h_atom->Fill(atomnum[0]);
      */
      h_cminb->Fill(z[0],cmin);
      h_cminb->Fill(z[0],cmax);
    }
    else {
      cout<<"number of particle is not 1 or 7  "<<j<<"   "<<n_p<<"   "<<z[0]<<endl;
    }
    cout<<100*j/Entry1<<" % finished   \r "<<flush;
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
  //need

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,4);
  
  //  h_z->SetStats(0);
  h_x->SetStats(1);
  can->cd(1);
  //  gPad->SetLogy(1);
  h_x->Draw("hist");

  //  h_rmax->SetStats(0);
  h_y->SetStats(1);
  can->cd(2);
  //  gPad->SetLogy(1);
  h_y->Draw("hist");

  //  h_rmax->SetStats(0);
  h_z->SetStats(1);
  can->cd(3);
  //  gPad->SetLogy(1);
  h_z->Draw("hist");

  //  h_rmax->SetStats(0);
  h_n->SetStats(1);
  can->cd(4);
  //  gPad->SetLogy(1);
  h_n->Draw("hist");

  //  h_rmax->SetStats(0);
  h_mass->SetStats(1);
  can->cd(5);
  //  gPad->SetLogy(1);
  h_mass->Draw("hist");

  //  h_rmax->SetStats(0);
  h_atom->SetStats(1);
  can->cd(6);
  //  gPad->SetLogy(1);
  h_atom->Draw("hist");

  //  h_rmax->SetStats(0);
  h_posi->SetStats(1);
  can->cd(7);
  //  gPad->SetLogy(1);
  h_posi->Draw("colz");

  //  h_rmax->SetStats(0);
  h_cmax->SetStats(1);
  can->cd(8);
  //  gPad->SetLogy(1);
  h_cmax->Draw("colz");

  //  h_rmax->SetStats(0);
  h_cmin->SetStats(1);
  can->cd(9);
  //  gPad->SetLogy(1);
  h_cmin->Draw("colz");

  //  h_rmax->SetStats(0);
  h_rmax->SetStats(1);
  can->cd(10);
  //  gPad->SetLogy(1);
  h_rmax->Draw("colz");

  //  h_rmax->SetStats(0);
  h_rmaxnd->SetStats(1);
  can->cd(11);
  //  gPad->SetLogy(1);
  h_rmaxnd->Draw("colz");

  //  h_rmax->SetStats(0);
  h_cminb->SetStats(1);
  can->cd(12);
  //  gPad->SetLogy(1);
  h_cminb->Draw("colz");

  //output for canvas
  //  TString outcan = "data/check_sim_macro_"+ to_string(subnum) +".C";
  TString outcan = "data/check_sim_macro_sequential_"+ to_string(subnum) +".C";
  //  TString outcan = "data1/check_sim_macro_hoylehoyle_nd_"+ to_string(subnum) +".C";
  can->SaveAs(outcan);
  
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
