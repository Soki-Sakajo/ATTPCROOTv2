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

int check_hoyle_state_0810(Int_t subnum=0){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  const Int_t dz=1;           // mm; bin of depth
  const Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  const Int_t r_attpc=275;    // mm; radius of attpc
  const Int_t l_attpc=1000;   // mm; depth of attpc
  const Int_t div_l=(int)(l_attpc/dz+0.9);

  //set constants
  const Int_t ctom=10;
  const Int_t GtoM =1e3;
  const Int_t GtoM2=1e6;
  const Double_t amu=931.49432; //atomic mass unit in MeV
  const Double_t c12_mass=12*amu;     //MeV
  const Double_t mg24_mass=24*amu-13.933578;    //MeV
  const Double_t h1_mass=amu+7.28897;    //MeV
  const Double_t alpha_mass=4*amu+2.4249;  //MeV
  const Double_t pi=3.1415926535;
  const Double_t rad_to_deg=180/pi;

  // for 12C12C
  Int_t eve1;
  Int_t n_p;
  Int_t massnum[7];
  Int_t atomnum[7];
  Int_t eve1_temp=0;
  Int_t z_n1=0;
  Double_t E[7];
  Double_t Esum6;
  Double_t Emax;
  Double_t Ealpha[7];
  Double_t Kgalpha[7];
  Double_t Esum12[2];
  Double_t Ex12c[2];
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
  Double_t stoppx[7];
  Double_t stoppy[7];
  Double_t stoppz[7];
  Double_t startx[7];
  Double_t starty[7];
  Double_t startz[7];
  Double_t x[7];
  Double_t y[7];
  Double_t z[7];
  Double_t px[7];
  Double_t py[7];
  Double_t pz[7];
  Double_t pxg[7];
  Double_t pyg[7];
  Double_t pzg[7];
  Double_t p[7];
  Double_t palpha[7];
  Double_t pbx;
  Double_t pby;
  Double_t pbz;
  Double_t pba;
  Double_t p_gx;
  Double_t p_gy;
  Double_t p_gz;
  Double_t p_ga;
  Double_t psum12x[2];
  Double_t psum12y[2];
  Double_t psum12z[2];
  Double_t psum12[2];
  Double_t px12[2];
  Double_t py12[2];
  Double_t pz12[2];
  Double_t pxc[2];
  Double_t pyc[2];
  Double_t pzc[2];
  Double_t pac[2];
  Double_t psp[2];
  Double_t k12c[2];
  Double_t cms_theta[2];

  vector<Double_t> ev_ve1;
  vector<Double_t> E_ve1;
  vector<Double_t> Emax_ve1;

  bool flag_dx7=false;
  bool flag_dy7=false;
  bool flag_dz7=false;
  
  //input
  //  TString refi="data/sim_macro_12c_nd_"+ to_string(subnum) +".root";
  TString refi="data1/sim_macro_hoylehoyle_nd_"+ to_string(subnum) +".root";
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
  tree1->SetBranchAddress("startpx",&px);
  tree1->SetBranchAddress("startpy",&py);
  tree1->SetBranchAddress("startpz",&pz);
  tree1->SetBranchAddress("stoppx",&stoppx);
  tree1->SetBranchAddress("stoppy",&stoppy);
  tree1->SetBranchAddress("stoppz",&stoppz);
  tree1->SetBranchAddress("energy",&E);
  tree1->SetBranchAddress("energy_max",&Emax);
  tree1->SetBranchAddress("beam_energy",&Etotal);
  tree1->SetBranchAddress("n",&n_p);
  tree1->SetBranchAddress("rmax_nd",&rmax_nd);
  tree1->SetBranchAddress("rmax",&rmax);
  tree1->SetBranchAddress("cmax",&cmax);
  tree1->SetBranchAddress("cmin",&cmin);

  //set histgram
  TH2F *h_posi=new TH2F("h_posi","reaction position:energy of alpha",nbin_zn,0,400,150,0,15);
  TH2F *h_dep=new TH2F("h_dep","reaction position: total energy of 6 alpha",200,0,200,40,0,40);
  TH2F *h_beam=new TH2F("h_beam","reaction position:beam energy",200,0,200,520,0,52);
  TH2F *h_hoy=new TH2F("h_hoy","reaction position: exenergy of hoyle",200,0,200,50,7,8);
  TH2F *h_hang=new TH2F("h_hang","reaction angle: energy of hoyle",180,0,180,30,5,10);
  TH2F *h_h1=new TH2F("h_h1","reaction position: reaction angle of hoyle",400,0,400,180,0,180);
  TH1F *h_ex=new TH1F("h_ex","excitation energy",100,7.6,7.8);
  TH2F *h_ang=new TH2F("h_ang","reaction angle: kienergy of hoyle",180,0,180,350,0,35);
  TH2F *h_ang1=new TH2F("h_ang1","reaction angle: energy of hoyle",180,0,180,350,0,35);
  TH2F *h_ex_dep=new TH2F("h_ex_dep","reaction position: excitation energy in dep range",20,100,110,50,5,10);
  TH2F *h_ex_ang=new TH2F("h_ex_ang","reaction angle: excitation energy in dep range",180,0,180,20,5,10);


  TH2F *h_rmax=new TH2F("h_rmax","reaction position: rmax",nbin_zn,0,400,160,0,160);
  TH2F *h_rmaxnd=new TH2F("h_rmaxnd","reaction position: rmax_nd",nbin_zn,0,400,160,0,160);

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    flag_dx7=true;
    flag_dy7=true;
    flag_dz7=true;

    if(n_p==7){
      for(int i=0; i<7; i++){
	x[i]=startx[i]*ctom;
	y[i]=starty[i]*ctom;
	z[i]=startz[i]*ctom;
	/*
	h_x->Fill(x[i]);
	h_y->Fill(y[i]);
	h_z->Fill(z[i]);
	h_mass->Fill(massnum[i]);
	h_atom->Fill(atomnum[i]);
      */
      }
      for(int i=1; i<7; i++){
	if(i<6){
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
	if(i==6){
	  if(fabs(x[i]-x[1])>1){
	    flag_dx7=false;
	  }
	  if(fabs(y[i]-y[1])>1){
	    flag_dy7=false;
	  }
	  if(fabs(z[i]-z[1])>dz){
	    flag_dz7=false;
	  }
	}
      }
      if(flag_dz7){
	h_beam->Fill(z[1],Etotal);
	p_gx=(px[1]+px[2]+px[3]+px[4]+px[5]+px[6])*GtoM/6;
	p_gy=(py[1]+py[2]+py[3]+py[4]+py[5]+py[6])*GtoM/6;
	p_gz=(pz[1]+pz[2]+pz[3]+pz[4]+pz[5]+pz[6])*GtoM/6;
	p_ga=sqrt(p_gx*p_gx+p_gy*p_gy+p_gz*p_gz); //abs vector in cms
	for(int i=1; i<7; i++){
	  h_posi->Fill(z[i],E[i]);
	  pxg[i]=px[i]*GtoM-p_gx;
	  pyg[i]=py[i]*GtoM-p_gy;
	  pzg[i]=pz[i]*GtoM-p_gz;
	  p[i]=sqrt(pxg[i]*pxg[i]+pyg[i]*pyg[i]+pzg[i]*pzg[i]);
	  Ealpha[i]=E[i]+alpha_mass;
	  palpha[i]=sqrt(Ealpha[i]*Ealpha[i]-alpha_mass*alpha_mass);
	}

	Esum6=E[1]+E[2]+E[3]+E[4]+E[5]+E[6];
	//need to check the delta E of hoyle and 3alpha
	//hoyle 1
	Esum12[0]=Ealpha[1]+Ealpha[2]+Ealpha[3];
	psum12x[0]=(px[1]+px[2]+px[3])*GtoM;
	psum12y[0]=(py[1]+py[2]+py[3])*GtoM;
	psum12z[0]=(pz[1]+pz[2]+pz[3])*GtoM;
	psum12[0]=sqrt(psum12x[0]*psum12x[0]+psum12y[0]*psum12y[0]+psum12z[0]*psum12z[0]);
	Ex12c[0]=sqrt(Esum12[0]*Esum12[0]-psum12[0]*psum12[0])-c12_mass;
	k12c[0]=Esum12[0]-sqrt(Esum12[0]*Esum12[0]-psum12[0]*psum12[0]);
	//hoyle 2
	Esum12[1]=Ealpha[4]+Ealpha[5]+Ealpha[6];
	psum12x[1]=(px[4]+px[5]+px[6])*GtoM;
	psum12y[1]=(py[4]+py[5]+py[6])*GtoM;
	psum12z[1]=(pz[4]+pz[5]+pz[6])*GtoM;
	psum12[1]=sqrt(psum12x[1]*psum12x[1]+psum12y[1]*psum12y[1]+psum12z[1]*psum12z[1]);
	Ex12c[1]=sqrt(Esum12[1]*Esum12[1]-psum12[1]*psum12[1])-c12_mass;
	k12c[1]=Esum12[1]-sqrt(Esum12[1]*Esum12[1]-psum12[1]*psum12[1]);

	px12[0]=(px[1]+px[2]+px[3])*GtoM/3/alpha_mass*c12_mass;
	py12[0]=(py[1]+py[2]+py[3])*GtoM/3/alpha_mass*c12_mass;
	pz12[0]=(pz[1]+pz[2]+pz[3])*GtoM/3/alpha_mass*c12_mass;
	px12[1]=(px[4]+px[5]+px[6])*GtoM/3/alpha_mass*c12_mass;
	py12[1]=(py[4]+py[5]+py[6])*GtoM/3/alpha_mass*c12_mass;
	pz12[1]=(pz[4]+pz[5]+pz[6])*GtoM/3/alpha_mass*c12_mass;

	/*
	if((Ex12c[0]<7.65 || Ex12c[0]>7.66)||(Ex12c[1]<7.65 || Ex12c[1]>7.66)){
	  //	if((Ex12c[0]>7.65 && Ex12c[0]<7.66)&&(Ex12c[1]>7.65 && Ex12c[1]<7.66)){
	  cout<<j<<"  excitation 0 "<<Ex12c[0]<<"  excitation 1 "<<Ex12c[1]<<endl;
	}
	*/
	//	cout<<j<<"  "<<k12c[0]<<"     "<<Esum12[0]<<"     "<<psum12[0]<<endl;
	//about angle in cms
	for(int k=0; k<2; k++){
	  pxc[k]=px12[k]-p_gx/alpha_mass*c12_mass;
	  pyc[k]=py12[k]-p_gy/alpha_mass*c12_mass;
	  pzc[k]=pz12[k]-p_gz/alpha_mass*c12_mass;
	  pac[k]=sqrt(pxc[k]*pxc[k]+pyc[k]*pyc[k]+pzc[k]*pzc[k]); //abs vector in cms
	  psp[k]=p_gx*pxc[k]+p_gy*pyc[k]+p_gz*pzc[k]; // scalar product
	  cms_theta[k]=rad_to_deg*acos(psp[k]/(p_ga*pac[k]));
	  h_ang->Fill(cms_theta[k],k12c[k]);
	  h_ex->Fill(Ex12c[k]);
	  h_hoy->Fill(z[1],Ex12c[k]);
	  if(Etotal>=45 && Etotal<50){
	    h_ang1->Fill(cms_theta[k],k12c[k]);
	  }
	}
	h_dep->Fill(z[1],Esum6);
	h_h1->Fill(z[1],cms_theta[0]);
	h_h1->Fill(z[1],cms_theta[1]);
	h_hang->Fill(cms_theta[0],Ex12c[0]);
	h_hang->Fill(cms_theta[1],Ex12c[1]);
	if(z[1]>100 && z[1]<=110){
	  h_ex_dep->Fill(z[1],Ex12c[0]);
	  h_ex_dep->Fill(z[1],Ex12c[1]);
	  h_ex_ang->Fill(cms_theta[0],Ex12c[0]);
	  h_ex_ang->Fill(cms_theta[1],Ex12c[1]);

	}
      }
      if(flag_dx7==false || flag_dy7==false){
	cout<<"this event is starting x or y position is not same  "<<j<<endl;
      }
      if(flag_dz7==false){
	cout<<"this event is starting z position is not same  "<<j<<endl;
      }
    }
    else if(n_p!=1){
      cout<<"number of particle is not 1 or 7  "<<j<<"   "<<n_p<<"   "<<z[0]<<endl;
    }
    //need
    //    cout<<100*j/Entry1<<" % finished   \r "<<flush;
  }
  
  //need

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(4,3);
  
  //  h_z->SetStats(0);
  h_posi->SetStats(1);
  can->cd(1);
  //  gPad->SetLogy(1);
  h_posi->Draw("colz");

  h_dep->SetStats(1);
  can->cd(2);
  //  gPad->SetLogy(1);
  h_dep->Draw("colz");

  h_beam->SetStats(1);
  can->cd(3);
  //  gPad->SetLogy(1);
  h_beam->Draw("colz");

  h_h1->SetStats(1);
  can->cd(4);
  //  gPad->SetLogy(1);
  h_h1->Draw("colz");

  h_ex->SetStats(1);
  can->cd(5);
  //  gPad->SetLogy(1);
  h_ex->Draw("colz");

  h_hoy->SetStats(1);
  can->cd(6);
  //  gPad->SetLogy(1);
  h_hoy->Draw("colz");

  h_hang->SetStats(1);
  can->cd(7);
  //  gPad->SetLogy(1);
  h_hang->Draw("colz");

  h_ang->SetStats(1);
  can->cd(8);
  //  gPad->SetLogy(1);
  h_ang->Draw("colz");

  h_ex_dep->SetStats(1);
  can->cd(9);
  //  gPad->SetLogy(1);
  h_ex_dep->Draw("colz");

  h_ex_ang->SetStats(1);
  can->cd(10);
  //  gPad->SetLogy(1);
  h_ex_ang->Draw("colz");

  h_rmax->SetStats(1);
  can->cd(11);
  //  gPad->SetLogy(1);
  h_rmax->Draw("colz");

  h_rmaxnd->SetStats(1);
  can->cd(12);
  //  gPad->SetLogy(1);
  h_rmaxnd->Draw("colz");
  
  //output for canvas
  //  TString outcan = "data/check_hoyle_state_"+ to_string(subnum) +".C";
  TString outcan = "data1/check_hoyle_state_"+ to_string(subnum) +".C";
  can->SaveAs(outcan);

  //need to here
  /*
  //  h_rmax->SetStats(0);
  h_ang1->SetStats(1);
  //  gPad->SetLogy(1);
  h_ang1->Draw("colz");
  */

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outcan << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
