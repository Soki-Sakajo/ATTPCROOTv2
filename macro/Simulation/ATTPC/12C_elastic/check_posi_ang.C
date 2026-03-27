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

int check_posi_ang(Int_t subnum=1){
  
  //set timer
  TStopwatch timer;
  timer.Start();

  //set situation
  Int_t z_s=6;          // Z of scatterd particl in simfile
  Int_t z_r1=6;         // Z of an recoil particl in simfile
  Int_t z_r2=1;         // Z of another recoil particl in simfile
  Int_t dz=2;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  Double_t P=0.13;      // atm; gas pressure
  Double_t da=1;        // deg; bin of angle
  Double_t dr=1;        // mm; bin of reaction rate
  Int_t nbin=int(275/dr);    // bin of rmax-reaction_rate hist
  Int_t zn_max=520; //hist max val
  Int_t nbin_zn=int(zn_max/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t N_rate2=10;     // (number of recoli partcle)/molecular  
  Int_t theta_min=3;    // deg; rutherford scattar angle
  Int_t theta_max=177;    // deg; rutherford scattar angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  Int_t ctom=10;
  //  Int_t pdg12c=1000060120; // Pdg Code of 12C
  //  Int_t pdgp=2212; // Pdg Code of p
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; // Pa
  Double_t alfa=7.297352569e-3; // fine structure constant
  Double_t hc=197.463587*1e-15; // MeV•m
  Double_t k_b=1.38065e-23; // J/K ; Bolzman constant
  Double_t rho=P*atmtoPa/(k_b*T);
  Double_t m12=11177.73600; // MeV/c^2 ; 12C mass
  
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  Double_t b=alfa*hc*z_s*z_r2/4;  
  
  //set variable
  // common
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Double_t da_r=da*deg_to_rad;
  
  // for 12C12C
  bool setbeam1=true;
  bool fin_beam1=false;
  bool fin_rea1=false;
  bool checkang=true;
  Int_t eve1, n_p1, n_fi1, n_ev1;
  Int_t n_ef1p=0;
  Int_t n_ee1p=0;
  Int_t n_ef1l=0;
  Int_t n_ee1l=0;
  Int_t n_ef1a=0;
  Int_t n_ee1a=0;
  Int_t n_ep1a=0;
  Int_t n_ef1v=0;
  Int_t n_ee1v=0;
  Int_t n_ep1v=0;
  Int_t n_ef1b=0;
  Int_t n_ee1b=0;
  Int_t n_ef1o=0;
  Int_t n_ee1o=0;
  Int_t n_fi1_temp=0;
  Int_t eve1_temp=0;
  Int_t z_n1=0;
  Int_t a_n1=0;
  Int_t a_n1_an=0;
  Int_t n_rmax1=0;
  Int_t n_rmax1_nd=0;
  Double_t z1_max=0;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1,dir1, dir1_an, Y1_eve, theta1, theta1_an;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t l1[div_l];
  Double_t n1[div_l];
  Double_t ang1[div_a];
  Double_t ang_r1[div_a];

  vector<vector<vector<Int_t>>> fi_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ang_an1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  
  //for 12Cp
  bool setbeam2=true;
  bool fin_beam2=false;
  bool fin_rea2=false;
  Int_t eve2, n_p2, n_fi2, n_ev2;
  Int_t n_ef2p=0;
  Int_t n_ee2p=0;
  Int_t n_ef2l=0;
  Int_t n_ee2l=0;
  Int_t n_ef2a=0;
  Int_t n_ee2a=0;
  Int_t n_ep2a=0;
  Int_t n_ef2v=0;
  Int_t n_ee2v=0;
  Int_t n_ep2v=0;
  Int_t n_ef2b=0;
  Int_t n_ee2b=0;
  Int_t n_ef2o=0;
  Int_t n_ee2o=0;
  Int_t n_fi2_temp=0;
  Int_t eve2_temp=0;
  Int_t z_n2=0;
  Int_t a_n2=0;
  Int_t n_rmax2=0;
  Int_t n_rmax2_nd=0;

  Double_t z2_max;
  Double_t i2_max;
  Double_t j2_max;
  Double_t dep2;
  Double_t dir2;
  Double_t theta2;
  Double_t theta2_an;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];
  Double_t l2[div_l];
  Double_t n2[div_l];
  Double_t ang2[div_a];
  Double_t ang_r2[div_a];

  vector<vector<vector<Int_t>>> fi_ve2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ang_an2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  
  //input multi-file by using TChain
  //common set for input
  TString dir="data1";
  
  // Tchain of 12c12c
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
    //      TString refi="sim_macro_c12_" + to_string(i) + ".root";
    //    TString refi="check_posi_ang_12c_" + to_string(i) + ".root";
    TString refi="check_posi_ang_12c_130matm_" + to_string(i) + ".root";
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
  
  // Tchain of 12Cp
  TChain *tree2=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    //    TString refi="check_posi_ang_p_" + to_string(i) + ".root";
    TString refi="check_posi_ang_p_130matm_" + to_string(i) + ".root";
    TString direfi=dir+"/"+refi;
    tree2->Add(direfi);
    if(i==0){cout<<"read from  "<<refi<<flush;}
    if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
  }
  cout<<"end reading files"<<endl<<endl;
  if(subnum==1){
    cout<<endl;
    cout<<"ATTENTION!! read ONLY one 12Cp file!!"<<endl<<endl<<endl;
  }
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve1);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);
  /*
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("rmax_nd",&rmax1_nd);
  tree1->SetBranchAddress("TotalE",&E1);
  */
  //  for 12C12p
  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("startz",&startz2);
  tree2->SetBranchAddress("theta_cms",&theta2_cms);
  tree2->SetBranchAddress("n",&n_p2);
  /*
  tree2->SetBranchAddress("rmax",&rmax2);
  tree2->SetBranchAddress("rmax_nd",&rmax2_nd);
  tree2->SetBranchAddress("TotalE",&E2);
  */

  //set histgram
  TH2F *h_ztheta_12c=new TH2F("h_ztheta_12c","reaction position:angle (12c12c)",nbin_zn,0,zn_max,180,0,180);
  TH2F *h_ztheta_p=new TH2F("h_ztheta_p","reaction position:angle (12cp)",nbin_zn,0,zn_max,180,0,180);
  
  //set output for making rutdata by each subnum
  //  TString outdata = "data/rutdata_c12_vector.root";
  //  TString outdata = "data/rutdata_c12_" + to_string(subnum) + ".root";
  TString outdata = "data1/check_" + to_string(subnum) + ".root";
  TFile *rutfile=new TFile(outdata, "recreate");
  TTree *rut12c = new TTree("rut12c","rut12c");
  TTree *exp12cp = new TTree("exp12cp","exp12cp");
  TTree *exp12cl = new TTree("exp12cl","exp12cl");
  TTree *exp12ca = new TTree("exp12ca","exp12ca");
  TTree *exp12cv = new TTree("exp12cv","exp12cv");
  TTree *exp12cb = new TTree("exp12cb","exp12cb");
  TTree *exp12co = new TTree("exp12co","exp12co");
  rut12c->Branch("file_number",&n_fi1,"n_fi1/I");
  rut12c->Branch("event_number",&n_ev1,"n_ev1/I");
  rut12c->Branch("reaction_position",&dep1,"dep1/D");
  rut12c->Branch("angle",&dir1,"dir1/D");
  /*
  rut12c->Branch("rmax_c12",&rmax1,"rmax1/D");
  rut12c->Branch("rmax_nd_c12",&rmax1_nd,"rmax1/D");
  rut12c->Branch("reaction_rate",&Y1_eve,"Y1_eve/D");
  rut12c->Branch("Eave",&Eave1,"Eave1/D");
  */
  exp12cp->Branch("file_particle",&n_ef1p,"n_ef1p/I");
  exp12cp->Branch("event_number_particle",&n_ee1p,"n_ee1p/I");
  exp12cl->Branch("file_startz",&n_ef1l,"n_ef1l/I");
  exp12cl->Branch("event_number_startz",&n_ee1l,"n_ee1l/I");
  exp12ca->Branch("file_angle",&n_ef1a,"n_fi1a/I");
  exp12ca->Branch("event_number_angle",&n_ee1a,"n_ee1a/I");
  exp12ca->Branch("number_of_particle_angle",&n_ep1a,"n_ep1a/I");
  exp12cv->Branch("file_angle_over90",&n_ef1v,"n_fi1v/I");
  exp12cv->Branch("event_number_angle_over90",&n_ee1v,"n_ee1v/I");
  exp12cv->Branch("number_of_particle_angle_over90",&n_ep1v,"n_ep1v/I");
  exp12cb->Branch("file_angle",&n_ef1b,"n_fi1b/I");
  exp12cb->Branch("event_number_angle",&n_ee1b,"n_ee1b/I");
  exp12co->Branch("file_other",&n_ef1o,"n_fi1o/I");
  exp12co->Branch("event_number_other",&n_ee1o,"n_ee1o/I");

  //for 12cp
  TTree *rutp = new TTree("rutp","rutp");
  TTree *exppp = new TTree("exppp","exppp");
  TTree *exppl = new TTree("exppl","exppl");
  TTree *exppa = new TTree("exppa","exppa");
  TTree *exppv = new TTree("exppv","exppv");
  TTree *exppb = new TTree("exppb","exppb");
  TTree *exppo = new TTree("exppo","exppo");
  rutp->Branch("file_number",&n_fi2,"n_fi2/I");
  rutp->Branch("event_number",&n_ev2,"n_ev2/I");
  rutp->Branch("reaction_position",&dep2,"dep2/D");
  rutp->Branch("angle",&dir2,"dir2/D");
  /*
  rutp->Branch("rmax_p",&rmax2,"rmax2/D");
  rutp->Branch("rmax_nd_p",&rmax2_nd,"rmax2/D");
  rutp->Branch("reaction_rate",&Y2_eve,"Y2_eve/D");
  rutp->Branch("Eave",&Eave2,"Eave2/D");
  */
  exppp->Branch("file_particle",&n_ef2p,"n_ef2p/I");
  exppp->Branch("event_number_particle",&n_ee2p,"n_ee2p/I");
  exppl->Branch("file_startz",&n_ef2l,"n_ef2l/I");
  exppl->Branch("event_number_startz",&n_ee2l,"n_ee2l/I");
  exppa->Branch("file_angle",&n_ef2a,"n_fi2a/I");
  exppa->Branch("event_number_angle",&n_ee2a,"n_ee2a/I");
  exppa->Branch("number_of_particle_angle",&n_ep2a,"n_ep2a/I");
  exppv->Branch("file_angle_over90",&n_ef2v,"n_fi2v/I");
  exppv->Branch("event_number_angle_over90",&n_ee2v,"n_ee2v/I");
  exppv->Branch("number_of_particle_angle_over90",&n_ep2v,"n_ep2v/I");
  exppb->Branch("file_angle",&n_ef2b,"n_fi2b/I");
  exppb->Branch("event_number_angle",&n_ee2b,"n_ee2b/I");
  exppo->Branch("file_other",&n_ef2o,"n_fi2o/I");
  exppo->Branch("event_number_other",&n_ee2o,"n_ee2o/I");

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    if(eve1_temp>eve1){n_fi1_temp++;}
    if(z1_max<z1[0]){z1_max=z1[0];}
    if(z1_max<z1[1]){z1_max=z1[1];}
    theta1=theta1_cms[0];
    theta1_an=theta1_cms[1];

    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0){
    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0 && theta1+theta_min<181){
    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0 && theta1<=theta_max){
      z_n1=int(z1[0]/dz);
      a_n1=int((theta1-theta_min)/da);
      a_n1_an=int(theta1_an+0.99999);
      ev_ve1[z_n1][a_n1].push_back(eve1);
      ang_an1[z_n1][a_n1].push_back(a_n1_an);
      fi_ve1[z_n1][a_n1].push_back(n_fi1_temp);
      /*
      E_ve1[z_n1][a_n1].push_back(E1);
      r_ve1[z_n1][a_n1].push_back(rmax1);
      r_nd_ve1[z_n1][a_n1].push_back(rmax1_nd);
      */
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;

    }
    else if(fabs(z1[0]-z1[1])<dz && theta1-theta_min>0){
      //    else if(fabs(z1[0]-z1[1])<dz && theta1-theta_min>0 && theta1<=90){
      n_ef1p=n_fi1_temp;
      n_ee1p=eve1;
      exp12cp->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p1==2 && theta1-theta_min>0){
      //    else if(n_p1==2 && theta1-theta_min>0 && theta1<=90){
      n_ef1l=n_fi1_temp;
      n_ee1l=eve1;
      exp12cl->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      //    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1<=90){
      n_ef1a=n_fi1_temp;
      n_ee1a=eve1;
      n_ep1a=n_p1;
      exp12ca->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      n_ef1v=n_fi1_temp;
      n_ee1v=eve1;
      n_ep1v=n_p1;
      exp12cv->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve1%2==0){
      n_ef1b=n_fi1_temp;
      n_ee1b=eve1;
      exp12cb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      n_ef1o=n_fi1_temp;
      n_ee1o=eve1;
      exp12co->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve1_temp=eve1;
  }
  i1_max=int(z1_max/dz);
  
  //12c12c calculate
  // devide depth
  
  //  for(int i=0;i < div_l ; i++){
  for(int i=0;i < ev_ve1.size(); i++){
    l1[i]=i*dz;
    n1[i]=0;

    if(n1[i]<1){continue;}
    
    //divide angle
    for(int j=0; j< ev_ve1.at(i).size(); j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;

      //need
      cout<<"12c calculating reaction rate at "<<l1[i]<<" mm "<<ang1[j]<<" deg  \r"<<flush;
    }
    if (int(l1[i])%100==0){
      //need
      cout<<"12C12C finish "<<l1[i]<<" mm !!                              "<<endl;
    }
  }
  
  //make conections between rmax and reaction rate
  for(int i=0; i< ev_ve1.size(); i++){
    l1[i]=i*dz;
    for(int j=0;j < ev_ve1.at(i).size() ; j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;
      for(int k=0; k< ev_ve1.at(i).at(j).size() ; k++){
	n_fi1=fi_ve1.at(i).at(j).at(k);
	n_ev1=ev_ve1.at(i).at(j).at(k);
	dep1=l1[i];
	dir1=ang1[j];
	dir1_an=ang_an1.at(i).at(j).at(k);
	h_ztheta_12c->Fill(l1[i],ang1[j]);
	rut12c->Fill();
      }
    }
    //need
    cout<<"Filling 12c data "<<100*i/ev_ve1.size() <<" % !    \r"<<flush;
  }
  
  //make event lists  for 12cp
  for(int j=0; j<Entry2; j++){
    tree2->GetEntry(j);
    z2[0]=startz2[0]*ctom;
    z2[1]=startz2[1]*ctom;
    if(eve2_temp>eve2){n_fi2_temp++;}
    if(z2_max<z2[0]){z2_max=z2[0];}
    if(z2_max<z2[1]){z2_max=z2[1];}
    theta2=theta2_cms[0];
    theta2_an=theta2_cms[1];

    //    if(n_p2==2 && fabs(z2[0]-z2[1])<dz && theta2-theta_min>=0){
    //    if(n_p2==2 && fabs(z2[0]-z2[1])<dz && theta2-theta_min>0 && theta2<=90){
    if(n_p2==2 && fabs(z2[0]-z2[1])<dz && theta2-theta_min>=0 && theta2<=theta_max){
      z_n2=int(z2[0]/dz);
      a_n2=int((theta2-theta_min)/da);
      ev_ve2[z_n2][a_n2].push_back(eve2);
      fi_ve2[z_n2][a_n2].push_back(n_fi2_temp);
      /*
      E_ve2[z_n2][a_n2].push_back(E2);
      r_ve2[z_n2][a_n2].push_back(rmax2);
      r_nd_ve2[z_n2][a_n2].push_back(rmax2_nd);
      */
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;
    }
    else if(fabs(z2[0]-z2[1])<dz && theta2-theta_min>0){
      //    else if(fabs(z1[0]-z1[1])<dz && theta1-theta_min>0 && theta1<=90){
      n_ef2p=n_fi2_temp;
      n_ee2p=eve2;
      exppp->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p2==2 && theta2-theta_min>0){
      //    else if(n_p1==2 && theta1-theta_min>0 && theta1<=90){
      n_ef2l=n_fi2_temp;
      n_ee2l=eve2;
      exppl->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p2==2 && fabs(z2[0]-z2[1])<dz){
      //    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1<=90){
      n_ef2a=n_fi2_temp;
      n_ee2a=eve2;
      n_ep2a=n_p2;
      exppa->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(n_p2==2 && fabs(z2[0]-z2[1])<dz){
      n_ef2v=n_fi2_temp;
      n_ee2v=eve2;
      n_ep2v=n_p2;
      exppv->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve2%2==0){
      n_ef2b=n_fi2_temp;
      n_ee2b=eve2;
      exppb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      n_ef2o=n_fi2_temp;
      n_ee2o=eve2;
      exppo->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve2_temp=eve2;
  }
  i2_max=int(z2_max/dz);
  
  //12cp calculate
  // devide depth
  
  //  for(int i=0;i < div_l ; i++){
  for(int i=0;i < ev_ve2.size(); i++){
    l2[i]=i*dz;
    n2[i]=0;
    if(n2[i]<1){continue;}
    
    //divide angle
    for(int j=0; j< ev_ve2.at(i).size(); j++){
      ang2[j]=j*da+theta_min;
      ang_r2[j]=ang2[j]*deg_to_rad;
      
      //need
      cout<<"p calculating reaction rate at "<<l2[i]<<" mm "<<ang2[j]<<" deg  \r"<<flush;
    }
    if (int(l2[i])%100==0){
      //need
      cout<<"12Cp finish "<<l2[i]<<" mm !!                              "<<endl;
    }
  }
  
  //make conections between rmax and reaction rate
  for(int i=0; i< ev_ve2.size(); i++){
    l2[i]=i*dz;
    for(int j=0;j < ev_ve2.at(i).size() ; j++){
      ang2[j]=j*da+theta_min;
      ang_r2[j]=ang2[j]*deg_to_rad;
      for(int k=0; k< ev_ve2.at(i).at(j).size() ; k++){
	n_fi2=fi_ve2.at(i).at(j).at(k);
	n_ev2=ev_ve2.at(i).at(j).at(k);
	dep2=l2[i];
	dir2=ang2[j];
	h_ztheta_p->Fill(dep2,dir2);
	rutp->Fill();

	//for check 12cp
	/*	
	//	if(l2[i]==l_ext){
	if(l2[i]==356 || l2[i]==30 ){
	  //	  h_thetaY1->Fill(ang1[j],Y1_eve);
	  //	  h_thetaY1->Fill(ang1[j],1);
	  //	  if(checkang &&(ang1[j]<=20 || ang1[j]>160)){
	  if(ang2[j]<=8){
	  //	  if(checkang){
	    cout<<l2[i]<<"   "<<ang2[j]<<"                              "<<endl;
	  //	    checkang=false;
	  }
	}	
	*/
      }
    }
    //need
    cout<<"Filling p data "<<100*i/ev_ve2.size() <<" % !    \r"<<flush;
  }
  
  //make error bar
  // for 12c12c
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta_12c->GetBinContent(m));
    h_ztheta_12c->SetBinError(m,err);
  }
  // for 12cp
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta_p->GetBinContent(m));
    h_ztheta_p->SetBinError(m,err);
  }

  //need

  //set margin and ofset
  Double_t ma_l=0.1;
  Double_t ma_r=0.15;
  Double_t ma_t=0.08;
  Double_t ma_b=0.13;
  Double_t of_y=0.8;
  Double_t of_z=0.8;

  //make and save the hist 
  //  TCanvas *can=new TCanvas("hists","hists");
  TCanvas *can=new TCanvas("hists","hists",2000,1000);
  //  can->Divide(1,2);
  // z-angle only 12c
  can->cd(1);
  gPad->SetLeftMargin(ma_l);
  gPad->SetRightMargin(ma_r);
  gPad->SetTopMargin(ma_t);
  gPad->SetBottomMargin(ma_b);
  h_ztheta_12c->SetStats(0);
  //  h_ztheta_12c->SetStats(1);
  h_ztheta_12c->GetXaxis()->SetTitle("Depth of AT-TPC [mm]");
  h_ztheta_12c->GetYaxis()->SetTitle("scattering angle [deg]");
  h_ztheta_12c->GetZaxis()->SetTitle("count");
  h_ztheta_12c->Draw("colz");
  h_ztheta_12c->GetXaxis()->SetLabelSize(0.06);
  h_ztheta_12c->GetYaxis()->SetLabelSize(0.06);
  h_ztheta_12c->GetZaxis()->SetLabelSize(0.06);
  h_ztheta_12c->GetXaxis()->SetTitleSize(0.06);
  h_ztheta_12c->GetYaxis()->SetTitleSize(0.06);
  h_ztheta_12c->GetZaxis()->SetTitleSize(0.06);
  h_ztheta_12c->GetYaxis()->SetTitleOffset(of_y);
  h_ztheta_12c->GetZaxis()->SetTitleOffset(of_z);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex1;
  latex1.SetTextSize(0.05);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(370,160,"^{12}C beam @ 50 MeV");
  latex1.DrawLatex(375,145,"iso-C_{4}H_{10} 0.13 atm.");
  latex1.DrawLatex(390,128,"^{12}C +^{12}C elastic");
  can->Update();

  TCanvas *can2=new TCanvas("hists2","hists2",2000,1000);
  // z-angle only p
  can2->cd(1);
  gPad->SetLeftMargin(ma_l);
  gPad->SetRightMargin(ma_r);
  gPad->SetTopMargin(ma_t);
  gPad->SetBottomMargin(ma_b);
  h_ztheta_p->SetStats(0);
  //  h_ztheta_p->SetStats(1);
  h_ztheta_p->GetXaxis()->SetTitle("Depth of AT-TPC [mm]");
  h_ztheta_p->GetYaxis()->SetTitle("scattering angle [deg]");
  h_ztheta_p->GetZaxis()->SetTitle("count");
  h_ztheta_p->Draw("colz");
  h_ztheta_p->GetXaxis()->SetLabelSize(0.06);
  h_ztheta_p->GetYaxis()->SetLabelSize(0.06);
  h_ztheta_p->GetZaxis()->SetLabelSize(0.06);
  h_ztheta_p->GetXaxis()->SetTitleSize(0.06);
  h_ztheta_p->GetYaxis()->SetTitleSize(0.06);
  h_ztheta_p->GetZaxis()->SetTitleSize(0.06);
  h_ztheta_p->GetYaxis()->SetTitleOffset(of_y);
  h_ztheta_p->GetZaxis()->SetTitleOffset(of_z);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.05);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(370,160,"^{12}C beam @ 50 MeV");
  latex2.DrawLatex(375,145,"iso-C_{4}H_{10} 0.13 atm.");
  latex2.DrawLatex(390,128,"^{12}C + p elastic");
  can->Update();

  //output for canvas
  /*
  TString outcan = "data1/check_posiang_hist_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  h_ztheta_12c->Write();
  h_ztheta_p->Write();
  */

  //  TString outcan = "test.pdf";
  //  TString outcan2 = "test_rut.pdf";
  //  TString outcan = "pdf_output/hist_posiang_130matm_" + to_string(subnum) + ".pdf";
  TString outcan = "pdf_output/hist_posiang_12c_130matm_" + to_string(subnum) + ".pdf";
  TString outcan2 = "pdf_output/hist_posiang_p_130matm_" + to_string(subnum) + ".pdf";
  can->SaveAs(outcan);
  can2->SaveAs(outcan2);
  
  //save branch
  rut12c ->AutoSave();
  exp12cp->AutoSave();
  exp12cl->AutoSave();
  exp12ca->AutoSave();
  exp12cv->AutoSave();
  exp12cb->AutoSave();
  exp12co->AutoSave();
  rutp ->AutoSave();
  exppp->AutoSave();
  exppl->AutoSave();
  exppa->AutoSave();
  exppv->AutoSave();
  exppb->AutoSave();
  exppo->AutoSave();
  rutfile->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "output: " << outdata << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
