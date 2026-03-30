//set result
#define RES_ON
//set canvas
#define CAN_ON

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

int ruther_macro_each_angle_nd_0925_24(Int_t subnum=1){ 
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  const Int_t z_s=6;               // Z of scatterd particl in simfile
  const Int_t z_r1=6;              // Z of an recoil particl in simfile
  const Int_t z_r2=1;              // Z of another recoil particl in simfile
  const Int_t dz=2;                // mm; bin of depth
  const Int_t n_b=1e5;             // /s; beam intensity
  const Double_t T=300;            // K; thermodynamic temparature
  const Double_t P=0.13;           // atm; gas pressure
  const Double_t da=1;             // deg; bin of angle
  const Double_t dr=1;             // mm; bin of reaction rate
  const Int_t nbin=int(275/dr);    // bin of rmax-reaction_rate hist
  const Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  const Int_t N_rate1=4;           // (number of scatterd particle)/molecular
  const Int_t N_rate2=10;          // (number of recoli partcle)/molecular  
  const Int_t theta_min=3;         // deg; rutherford scattar angle
  const Int_t theta_max=177;       // deg; rutherford scattar angle
  const Int_t r_attpc=275;         // mm; radius of attpc
  const Int_t l_attpc=1000;        // mm; depth of attpc
  const Int_t l_ext=100;           // mm; extract information
  Double_t ra_fit=-100;
  Double_t y_fit=-100;
  
  //set constants
  const Int_t ctom=10;
  const Double_t pi=3.1415926535;
  const Double_t rad_to_deg=180/pi;
  const Double_t deg_to_rad=pi/180;
  const Double_t atmtoPa=101325;       // Pa
  const Double_t alfa=7.297352569e-3;  // fine structure constant
  const Double_t hc=197.463587*1e-15;  // MeV•m
  const Double_t k_b=1.38065e-23;      // J/K ; Bolzman constant
  const Double_t rho=P*atmtoPa/(k_b*T);

  //constants in ruthford scatter
  const Double_t a=alfa*hc*z_s*z_r1/4;  
  const Double_t b=alfa*hc*z_s*z_r2/4;  
  
  //set variable
  // common
  const Int_t div_l=(int)(l_attpc/dz+0.9);
  const Int_t div_r=(int)(r_attpc/dr+0.9);
  const Int_t div_a=(int)((180-theta_min)/da+0.9);
  const Double_t da_r=da*deg_to_rad;
  
  // for 12C12C
  bool setbeam1=true;
  bool fin_beam1=false;
  bool fin_rea1=false;
  bool checkang=true;
  Int_t eve1;
  Int_t n_p1;
  Int_t n_fi1;
  Int_t n_ev1;
  Int_t n_en1p=0;
  Int_t n_ef1p=0;
  Int_t n_ee1p=0;
  Int_t n_ep1p=0;
  Int_t n_en1l=0;
  Int_t n_ef1l=0;
  Int_t n_ee1l=0;
  Double_t n_ez1l[2]={0,0};
  Int_t n_en1a=0;
  Int_t n_ef1a=0;
  Int_t n_ee1a=0;
  Double_t n_ea1a[2]={0,0};
  Int_t n_en1b=0;
  Int_t n_ef1b=0;
  Int_t n_ee1b=0;
  Int_t n_en1o=0;
  Int_t n_ef1o=0;
  Int_t n_ee1o=0;
  Int_t n_ep1o=0;
  Double_t n_ez1o[2]={0,0};
  Double_t n_ea1o[2]={0,0};
  Int_t n_fi1_temp=0;
  Int_t eve1_temp=0;
  Int_t z_n1=0;
  Int_t a_n1=0;
  Int_t a_n1_an=0;
  Int_t n_rmax1=0;
  Int_t n_rmax1_nd=0;
  Double_t rmax1;
  Double_t rmax1_nd;
  Double_t E1;
  Double_t z1_max=0;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1;
  Double_t dir1;
  Double_t dir1_an;
  Double_t Y1_eve;
  Double_t theta1;
  Double_t theta1_an;
  Double_t Y1_sum=0;
  Double_t Y1_rmax=0;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t l1[div_l];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  Double_t Y1[div_l];
  Double_t ang1[div_a];
  Double_t ang_r1[div_a];

  vector<vector<Double_t>> I1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> n_beam1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_da(div_l, vector<Double_t>(div_a,0));  
  vector<vector<vector<Int_t>>> fi_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ang_an1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve1(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve1(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_nd_ve1(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  
  //for 12Cp
  bool setbeam2=true;
  bool fin_beam2=false;
  bool fin_rea2=false;
  Int_t eve2;
  Int_t n_p2;
  Int_t n_fi2;
  Int_t n_ev2;
  Int_t n_en2p=0;
  Int_t n_ef2p=0;
  Int_t n_ee2p=0;
  Int_t n_ep2p=0;
  Int_t n_en2l=0;
  Int_t n_ef2l=0;
  Int_t n_ee2l=0;
  Double_t n_ez2l[2]={0,0};
  Int_t n_en2a=0;
  Int_t n_ef2a=0;
  Int_t n_ee2a=0;
  Double_t n_ea2a[2]={0,0};
  Int_t n_en2b=0;
  Int_t n_ef2b=0;
  Int_t n_ee2b=0;
  Int_t n_en2o=0;
  Int_t n_ef2o=0;
  Int_t n_ee2o=0;
  Int_t n_ep2o=0;
  Double_t n_ez2o[2]={0,0};
  Double_t n_ea2o[2]={0,0};
  Int_t n_fi2_temp=0;
  Int_t eve2_temp=0;
  Int_t z_n2=0;
  Int_t a_n2=0;
  Int_t n_rmax2=0;
  Int_t n_rmax2_nd=0;
  Double_t rmax2;
  Double_t rmax2_nd;
  Double_t E2;
  Double_t z2_max;
  Double_t i2_max;
  Double_t j2_max;
  Double_t dep2;
  Double_t dir2;
  Double_t Y2_eve;
  Double_t theta2;
  Double_t theta2_an;
  Double_t Y2_sum=0;
  Double_t Y2_rmax=0;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];
  Double_t l2[div_l];
  Double_t n2[div_l];
  Double_t Esum2[div_l];
  Double_t Eave2[div_l];
  Double_t Y2[div_l];
  Double_t ang2[div_a];
  Double_t ang_r2[div_a];

  vector<vector<Double_t>> I2(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> n_beam2(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y2_da(div_l, vector<Double_t>(div_a,0));  
  vector<vector<vector<Int_t>>> fi_ve2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ang_an2(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve2(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve2(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_nd_ve2(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  
  //for elastic (12c12c, 12cp)
  bool setbeam=true;
  bool fin_beam=false;
  bool fin_rea=false;
  string scatter="";
  Int_t n_fi;
  Int_t n_ev;
  string scattap="";
  Int_t n_enp=0;
  Int_t n_efp=0;
  Int_t n_eep=0;
  Int_t n_epp=0;
  string scattal="";
  Int_t n_enl=0;
  Int_t n_efl=0;
  Int_t n_eel=0;
  Double_t n_ezl[2]={0,0};
  string scattaa="";
  Int_t n_ena=0;
  Int_t n_efa=0;
  Int_t n_eea=0;
  Double_t n_eaa[2]={0,0};
  string scattab="";
  Int_t n_enb=0;
  Int_t n_efb=0;
  Int_t n_eeb=0;
  string scattao="";
  Int_t n_eno=0;
  Int_t n_efo=0;
  Int_t n_eeo=0;
  Int_t n_epo=0;
  Double_t n_ezo[2]={0,0};
  Double_t n_eao[2]={0,0};
  Int_t n_fi_temp=0;
  Int_t eve_temp=0;
  Int_t z_n=0;
  Int_t a_n=0;
  Int_t a_n_an=0;
  Int_t n_rmax=0;
  Int_t n_rmax_nd=0;
  Double_t z_max;
  Double_t i_max;
  Double_t j_max;
  Double_t rmax;
  Double_t rmax_nd;
  Double_t dep;
  Double_t dir;
  Double_t dir_an;
  Double_t Y_eve;
  Double_t theta;
  Double_t theta_an;
  Double_t Y_sum=0;
  Double_t Y_rmax=0;
  Double_t z[2];
  Double_t l[div_l];
  Double_t n[div_l][2];
  Double_t Esum[div_l][2];
  Double_t Eave[div_l][2];
  Double_t Y[div_l];
  Double_t ang[div_a];
  Double_t ang_r[div_a];
  
  vector<vector<Int_t>> evsig(div_l, vector<Int_t>(div_a,0));
  vector<vector<Double_t>> n_beam(div_l, vector<Double_t>(div_a,0));
  vector<vector<vector<Double_t>>> I(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(2,0)));
  vector<vector<vector<Double_t>>> Y_da(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(2,0)));

  vector<vector<vector<Int_t>>> fi_ve(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<string>>> evche(div_l,vector<vector<string>>(div_a,vector<string>(0)));
  vector<vector<vector<Int_t>>> ang_an(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_nd_ve(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  
  //input multi-file by using TChain
  //common set for input
  TString direc="data1";
  
  // Tchain of 12c12c
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    TString refi="sim_macro_12c_nd_" + to_string(i) + ".root";
    TString direfi=direc+"/"+refi;
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
    TString refi="sim_macro_p_nd_" + to_string(i) + ".root";
    TString direfi=direc+"/"+refi;
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
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("rmax_nd",&rmax1_nd);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);
  //  for 12C12p
  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("startz",&startz2);
  tree2->SetBranchAddress("rmax",&rmax2);
  tree2->SetBranchAddress("rmax_nd",&rmax2_nd);
  tree2->SetBranchAddress("TotalE",&E2);
  tree2->SetBranchAddress("theta_cms",&theta2_cms);
  tree2->SetBranchAddress("n",&n_p2);

  //set histgram
  TH1F *h_z=new TH1F("h_z","reaction position:reaction_rate",nbin_zn,0,400);
  TH1F *h_rmax=new TH1F("h_rmax","rmax:reaction_rate",nbin,0,275);
  TH1F *h_rmax_nd=new TH1F("h_rmax_nd","rmax_nd:reaction_rate",nbin,0,275);
  TH1F *h_rmax_a_nd=new TH1F("h_rmax_a_nd","rmax_nd:reaction_rate (all)",nbin,0,275);
  
  TH2F *h_ztheta_12c=new TH2F("h_ztheta_12c","reaction position:angle (12c12c)",nbin_zn,0,400,100,0,200);
  TH1F *h_z_12c=new TH1F("h_z_12c","reaction position:reaction_rate (12c12c)",nbin_zn,0,400);
  TH1F *h_rmax_12c=new TH1F("h_rmax_12c","rmax:reaction_rate (12c12c)",nbin,0,275);
  TH2F *h_ztheta_p=new TH2F("h_ztheta_p","reaction position:angle (12cp)",nbin_zn,0,400,100,0,200);
  TH1F *h_z_p=new TH1F("h_z_p","reaction position:reaction_rate (12cp)",nbin_zn,0,400);
  TH1F *h_rmax_p=new TH1F("h_rmax_p","rmax:reaction_rate (12cp)",nbin,0,275);
  TH2F *h_ztheta_all=new TH2F("h_ztheta_all","reaction position:angle (all)",nbin_zn,0,400,100,0,200);
  TH1F *h_z_all=new TH1F("h_z_all","reaction position:reaction_rate (all)",nbin_zn,0,400);
  TH1F *h_rmax_all=new TH1F("h_rmax_all","rmax:reaction_rate (all)",nbin,0,275);
  
  //  TH1F *h_thetaY1=new TH1F("h_thetaY1","angle :reaction_rate",180,0,180);
  TH2F *h_thetaY1=new TH2F("h_thetaY1","reaction position :rmax",nbin_zn,0,400,nbin,0,275);
  TH1F *h_thetaY2=new TH1F("h_thetaY2","angle :reaction_rate",180,0,180);
  TH1F *h_thetaY3=new TH1F("h_thetaY3","angle :reaction_rate",180,0,180);
  TH1F *h_fit=new TH1F("h_fit","angle :ratio of number of reaction and fit",180,0,180);
  
  TF1 *f1=new TF1("f1","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2))+1/(cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)))",0,180);
  TF1 *f2=new TF1("f2","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)))",0,180);
  TF1 *f3=new TF1("f3","log10([0])+log10(sin([1]*x))-4*(log10(sin([1]*x/2)))",theta_min,theta_max);
  
  //set output for making rutdata by each subnum
  TString outdata = "data1/rutdata_nd_0925_" + to_string(subnum) + ".root";
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
  rut12c->Branch("rmax_c12",&rmax1,"rmax1/D");
  rut12c->Branch("rmax_nd_c12",&rmax1_nd,"rmax1/D");
  rut12c->Branch("reaction_rate",&Y1_eve,"Y1_eve/D");
  rut12c->Branch("Eave",&Eave1,"Eave1/D");
  exp12cp->Branch("file_particle",&n_ef1p,"n_ef1p/I");
  exp12cp->Branch("event_number_particle",&n_ee1p,"n_ee1p/I");
  exp12cp->Branch("number_of_particle",&n_ep1p,"n_ee1p/I");
  exp12cl->Branch("file_startz",&n_ef1l,"n_ef1l/I");
  exp12cl->Branch("event_number_startz",&n_ee1l,"n_ee1l/I");
  exp12cl->Branch("startz_not_same_startz",&n_ez1l,"n_ez1l/D");
  exp12ca->Branch("file_angle",&n_ef1a,"n_fi1a/I");
  exp12ca->Branch("event_number_angle",&n_ee1a,"n_ee1a/I");
  exp12ca->Branch("angle_not_range_angle",&n_ea1a,"n_ea1a/D");
  exp12cb->Branch("file_beam",&n_ef1b,"n_fi1b/I");
  exp12cb->Branch("event_number_beam",&n_ee1b,"n_ee1b/I");
  exp12co->Branch("file_other",&n_ef1o,"n_fi1o/I");
  exp12co->Branch("event_number_other",&n_ee1o,"n_ee1o/I");
  exp12co->Branch("startz_other",&n_ez1o,"n_ez1o/D");
  exp12co->Branch("angle_other",&n_ea1o,"n_ea1o/D");
  
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
  rutp->Branch("rmax_p",&rmax2,"rmax2/D");
  rutp->Branch("rmax_nd_p",&rmax2_nd,"rmax2/D");
  rutp->Branch("reaction_rate",&Y2_eve,"Y2_eve/D");
  rutp->Branch("Eave",&Eave2,"Eave2/D");
  exppp->Branch("file_particle",&n_ef2p,"n_ef2p/I");
  exppp->Branch("event_number_particle",&n_ee2p,"n_ee2p/I");
  exppp->Branch("number_of_particle",&n_ep2p,"n_ep2p/I");
  exppl->Branch("file_startz",&n_ef2l,"n_ef2l/I");
  exppl->Branch("event_number_startz",&n_ee2l,"n_ee2l/I");
  exppl->Branch("startz_not_same_startz",&n_ez2l,"n_ez2l/D");
  exppa->Branch("file_angle",&n_ef2a,"n_fi2a/I");
  exppa->Branch("event_number_angle",&n_ee2a,"n_ee2a/I");
  exppa->Branch("angle_not_range_angle",&n_ea2a,"n_ep2a/I");
  exppb->Branch("file_beam",&n_ef2b,"n_fi2b/I");
  exppb->Branch("event_number_beam",&n_ee2b,"n_ee2b/I");
  exppo->Branch("file_other",&n_ef2o,"n_fi2o/I");
  exppo->Branch("event_number_other",&n_ee2o,"n_ee2o/I");
  exppo->Branch("startz_other",&n_ez2o,"n_ez2o/D");
  exppo->Branch("angle_other",&n_ea2o,"n_ea2o/D");
  
  //for all (12c12c 12cp)
  TTree *rutal = new TTree("rutal","rutal");
  TTree *expap = new TTree("expap","expap");
  TTree *expal = new TTree("expal","expal");
  TTree *expaa = new TTree("expaa","expaa");
  TTree *expav = new TTree("expav","expav");
  TTree *expab = new TTree("expab","expab");
  TTree *expao = new TTree("expao","expao");
  rutal->Branch("file_number",&n_fi,"n_fi/I");
  rutal->Branch("event_number",&n_ev,"n_ev/I");
  rutal->Branch("event_check",&scatter,"scatter/C");
  rutal->Branch("reaction_position",&dep,"dep/D");
  rutal->Branch("angle",&dir,"dir/D");
  rutal->Branch("rmax_p",&rmax,"rmax/D");
  rutal->Branch("rmax_nd_p",&rmax_nd,"rmax/D");
  rutal->Branch("reaction_rate",&Y_eve,"Y_eve/D");
  rutal->Branch("Eave",&Eave,"Eave/D");
  expap->Branch("event",&scattap,"scattap/C");
  expap->Branch("file_particle",&n_efp,"n_efp/I");
  expap->Branch("event_number_particle",&n_eep,"n_eep/I");
  expap->Branch("number_of_particle",&n_epp,"n_epp/I");
  expal->Branch("event",&scattal,"scattal/C");
  expal->Branch("file_startz",&n_efl,"n_efl/I");
  expal->Branch("event_number_startz",&n_eel,"n_eel/I");
  expal->Branch("startz_not_same_startz",&n_ezl,"n_ezl/I");
  expaa->Branch("event",&scattaa,"scattaa/C");
  expaa->Branch("file_angle",&n_efa,"n_fia/I");
  expaa->Branch("event_number_angle",&n_eea,"n_eea/I");
  expaa->Branch("angle_not_range_angle",&n_eaa,"n_eaa/I");
  expab->Branch("event",&scattab,"scattab/C");
  expab->Branch("file_angle",&n_efb,"n_fib/I");
  expab->Branch("event_number_angle",&n_eeb,"n_eeb/I");
  expao->Branch("event",&scattao,"scattao/C");
  expao->Branch("file_other",&n_efo,"n_fio/I");
  expao->Branch("event_number_other",&n_eeo,"n_eeo/I");
  exppo->Branch("startz_other",&n_ezo,"n_ezo/D");
  exppo->Branch("angle_other",&n_eao,"n_eao/D");


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
    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1>=theta_min && theta1<=theta_max){
      z_n1=int(z1[0]/dz);
      a_n1=int((theta1-theta_min)/da);
      a_n1_an=int(theta1_an+0.99999);
      ev_ve1[z_n1][a_n1].push_back(eve1);
      ang_an1[z_n1][a_n1].push_back(a_n1_an);
      E_ve1[z_n1][a_n1].push_back(E1);
      r_ve1[z_n1][a_n1].push_back(rmax1);
      r_nd_ve1[z_n1][a_n1].push_back(rmax1_nd);
      fi_ve1[z_n1][a_n1].push_back(n_fi1_temp);
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;
    }
    else if(fabs(z1[0]-z1[1])<dz && theta1>=theta_min && theta1<=theta_max){
      n_en1p++;
      n_ef1p=n_fi1_temp;
      n_ee1p=eve1;
      n_ep1p=n_p1;
      exp12cp->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p1==2 && theta1>=theta_min && theta1<=theta_max){
      n_en1l++;
      n_ef1l=n_fi1_temp;
      n_ee1l=eve1;
      n_ez1l[0]=z1[0];
      n_ez1l[1]=z1[1];
      exp12cl->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      n_en1a++;
      n_ef1a=n_fi1_temp;
      n_ee1a=eve1;
      n_ea1a[0]=theta1;
      n_ea1a[1]=theta1_an;
      exp12ca->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve1%2==0){
      n_en1b++;
      n_ef1b=n_fi1_temp;
      n_ee1b=eve1;
      exp12cb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      n_en1o++;
      n_ef1o=n_fi1_temp;
      n_ee1o=eve1;
      n_ep1o=n_p1;
      n_ez1o[0]=z1[0];
      n_ez1o[1]=z1[1];
      n_ea1o[0]=theta1;
      n_ea1o[1]=theta1_an;
      exp12co->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve1_temp=eve1;
    if(j%1000==0){
      //need
      cout<<"Sorting 12c12c data "<<100*j/Entry1 <<" % !    \r"<<flush;
    }
  }
  i1_max=int(z1_max/dz);
  cout<<"about eliminated 12c12c events"<<endl;
  cout<<" number of particles not 2 "<<n_en1p<<"  events, not same start position "<<n_en1l<<"  events,"<<endl;
  cout<<" angle not in the range "<<n_en1a<<"  events, complex reason "<<n_en1o<<"  events"<<endl;
  
  //12c12c calculate
  //devide depth
  for(int i=0;i < ev_ve1.size(); i++){
    l1[i]=i*dz;
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;
    Y1[i]=0;

    //calculate Eave
    for(int j=0; j< ev_ve1.at(i).size(); j++){
      for(int k=0; k< ev_ve1.at(i).at(j).size(); k++){
      	n1[i]++;
      	Esum1[i]+=E_ve1.at(i).at(j).at(k);
      	if(i==i1_max){
      	  if(ev_ve1.at(i).at(j).size()>0){j1_max=j;}
      	}
      }
    }
    if(n1[i]<1){continue;}
    else{Eave1[i]=Esum1[i]/(2*n1[i]);}

    //calculate reaction rate
    //divide angle
    for(int j=0; j< ev_ve1.at(i).size(); j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;

      //need
      cout<<"12c calculating reaction rate at "<<l1[i]<<" mm "<<ang1[j]<<" deg  \r"<<flush;

      if(ev_ve1.at(i).at(j).size()<1){continue;}
      I1.at(i).at(j)=2*pi*a*a*da_r*sin(ang_r1[j])
        /(Eave1[i]*Eave1[i]
          *sin(ang_r1[j]/2)*sin(ang_r1[j]/2)
          *sin(ang_r1[j]/2)*sin(ang_r1[j]/2));      
      if(setbeam1){
        n_beam1.at(i).at(j)=n_b;
      	setbeam1=false;
      }
      else{n_beam1.at(i).at(j)=n_b-Y1_sum;}
      if(n_beam1.at(i).at(j)<0){
        cout<<"                                                                 "<<endl;
      	cout<<"end the loop at"<<l1[i]<<" mm~"<<l1[i]+dz<<" mm, "<<ang1[j]<<"deg~"<<ang1[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
      	fin_beam1=true;
      	break;
      }
      else{
        Y1_da.at(i).at(j)=I1.at(i).at(j)*n_beam1.at(i).at(j)*rho*dz*1e-3*N_rate1;
        Y1[i]+=Y1_da.at(i).at(j);
      	if(l1[i]==l_ext && theta_min<=ang1[j] && theta_max>=ang1[j]){
            Double_t lY1=log10(Y1_da.at(i).at(j));
            //	  h_thetaY3->Fill(ang1[j],log10(Y1_eve));
            h_thetaY3->Fill(ang1[j],lY1);
            h_fit->Fill(ang1[j],Y1_da.at(i).at(j));
            //	  cout<<l1[i]<<"   "<<ang1[j]<<"    "<<lY1<<"                          "<<endl;
        	}
      }
      //      cout<<i<<"       "<<j<<"      "<<ev_ve.at(i).at(j).size()<<endl;
      if(i==i1_max && j==j1_max){
        cout<<"                                                                 "<<endl;
        cout<<"end the loop at "<<l1[i]<<" mm~"<<l1[i]+dz<<" mm "
        <<ang1[j]<<" deg~"<<ang1[j]+da<<" deg, because of no reaction in over "
        <<l1[i]+dz<<" mm "<<ang1[j]+da<<" deg"<<endl<<endl;
        //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
        fin_rea1=true;
        break;
      }
    }
    Y1_sum+=Y1[i];
    
    //z-reaction rate hist filling data
    if(n1[i]>0){
      h_z->Fill(l1[i],Y1[i]);
      h_z_12c->Fill(l1[i],Y1[i]);
    }
    if (int(l1[i])%100==0){
      //need
      cout<<"12C12C finish "<<l1[i]<<" mm !!                              "<<endl;
    }
    if(fin_beam1){break;}
    if(fin_rea1){break;}
    //    if(z1_max<l1[i]+dz){}
  }

  //make conections between rmax and reaction rate
  for(int i=0; i< ev_ve1.size(); i++){
    l1[i]=i*dz;
    for(int j=0;j < ev_ve1.at(i).size() ; j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;
      checkang=true;
      for(int k=0; k< ev_ve1.at(i).at(j).size() ; k++){
        n_fi1=fi_ve1.at(i).at(j).at(k);
        n_ev1=ev_ve1.at(i).at(j).at(k);
        dep1=l1[i];
        dir1=ang1[j];
        dir1_an=ang_an1.at(i).at(j).at(k);
        Y1_eve=Y1_da.at(i).at(j)/ev_ve1.at(i).at(j).size();
        rmax1=r_ve1.at(i).at(j).at(k);
        rmax1_nd=r_nd_ve1.at(i).at(j).at(k);
        h_rmax->Fill(rmax1,Y1_eve);
        h_rmax_12c->Fill(rmax1,Y1_eve);
        h_rmax_nd->Fill(rmax1_nd,Y1_eve);
        n_rmax1++;
        n_rmax1_nd++;
        Y1_rmax+=Y1_eve;
        h_ztheta_12c->Fill(l1[i],ang1[j]);
        rut12c->Fill();
        /*
        //rmax test
        h_thetaY1->Fill(dep1,rmax1);
        if(rmax1>150){
          cout<<"strange event?  "<<n_fi1<<"   "<<n_ev1<<"   "<<dep1<<"   "<<dir1<<"   "<<dir1_an<<"   "<<rmax1<<endl;
        }
        //	if(n_fi1==0 && rmax1>150){
        if(rmax1>150){
          cout<<"file "<<n_fi1<<"   event "<<n_ev1<<"  depth "<<dep1<<" mm  angle "<<dir1<<" deg   event is rmax "<<rmax1<<" mm "<<endl;
        }

        //for check 12c12c
        if(l1[i]==l_ext){
        //	  h_thetaY1->Fill(ang1[j],Y1_eve);
        //	  h_thetaY1->Fill(ang1[j],1);
        //	  if(checkang &&(ang1[j]<=20 || ang1[j]>160)){
        //	  if(ang1[j]<=8){
        //	  if(checkang){
        //	    cout<<l1[i]<<"   "<<ang1[j]<<"                              "<<endl;
        //	    checkang=false;
        //    }
        h_thetaY1->Fill(dir1,Y1_eve);
        h_thetaY1->Fill(dir1_an,Y1_eve);
        h_thetaY2->Fill(ang1[j],Y1_eve);
        }
        if(l1[i]==142){
          h_thetaY2->Fill(ang1[j],Y1_eve);
          //	  h_thetaY2->Fill(180-ang1[j],Y1_eve);
        }
        if(l1[i]==5){
          h_thetaY3->Fill(ang1[j],Y1_eve);
          h_thetaY3->Fill(180-ang1[j],Y1_eve);
        }}
      	*/
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
    if(n_p2==2 && fabs(z2[0]-z2[1])<dz && theta2>=theta_min && theta2<=theta_max){
      z_n2=int(z2[0]/dz);
      a_n2=int((theta2-theta_min)/da);
      ev_ve2[z_n2][a_n2].push_back(eve2);
      E_ve2[z_n2][a_n2].push_back(E2);
      r_ve2[z_n2][a_n2].push_back(rmax2);
      r_nd_ve2[z_n2][a_n2].push_back(rmax2_nd);
      fi_ve2[z_n2][a_n2].push_back(n_fi2_temp);
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;
    }
    else if(fabs(z2[0]-z2[1])<dz && theta2>=theta_min && theta2<=theta_max){
      n_en2p++;
      n_ef2p=n_fi2_temp;
      n_ee2p=eve2;
      n_ep2p=n_p2;
      exppp->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p2==2 && theta2>=theta_min && theta2<=theta_max){
      n_en2l++;
      n_ef2l=n_fi2_temp;
      n_ee2l=eve2;
      n_ez2l[0]=z2[0];
      n_ez2l[1]=z1[1];
      exppl->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p2==2 && fabs(z2[0]-z2[1])<dz){
      n_en2a++;
      n_ef2a=n_fi2_temp;
      n_ee2a=eve2;
      n_ea2a[0]=theta2;
      n_ea2a[1]=theta2_an;
      exppa->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve2%2==0){
      n_en2b++;
      n_ef2b=n_fi2_temp;
      n_ee2b=eve2;
      exppb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      n_en2o++;
      n_ef2o=n_fi2_temp;
      n_ee2o=eve2;
      n_ep2o=n_p2;
      n_ez2o[0]=z2[0];
      n_ez2o[1]=z2[1];
      n_ea2o[0]=theta2;
      n_ea2o[1]=theta2_an;
      exppo->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve2_temp=eve2;
    if(j%1000==0){
      //need
      cout<<"Sorting 12cp data "<<100*j/Entry2 <<" % !    \r"<<flush;
    }
  }
  i2_max=int(z2_max/dz);
  cout<<"about eliminated 12cp events"<<endl;
  cout<<" number of particles not 2 "<<n_en2p<<"  events, not same start position "<<n_en2l<<"  events,"<<endl;
  cout<<" angle not in the range "<<n_en2a<<"  events, complex reason "<<n_en2o<<"  events"<<endl;
  
  //12cp calculate
  // devide depth
  for(int i=0;i < ev_ve2.size(); i++){
    l2[i]=i*dz;
    n2[i]=0;
    Esum2[i]=0;    
    Eave2[i]=0;
    Y2[i]=0;
    
    //calculate Eave
    for(int j=0; j< ev_ve2.at(i).size(); j++){
      for(int k=0; k< ev_ve2.at(i).at(j).size(); k++){
        n2[i]++;
        Esum2[i]+=E_ve2.at(i).at(j).at(k);
        if(i==i2_max){
          if(ev_ve2.at(i).at(j).size()>0){j2_max=j;}
        }
      }
    }
    if(n2[i]<1){continue;}
    else{Eave2[i]=Esum2[i]/(2*n2[i]);}
    
    //calculate reaction rate
    //divide angle
    for(int j=0; j< ev_ve2.at(i).size(); j++){
      ang2[j]=j*da+theta_min;
      ang_r2[j]=ang2[j]*deg_to_rad;
      
      //need
      cout<<"p calculating reaction rate at "<<l2[i]<<" mm "<<ang2[j]<<" deg  \r"<<flush;
      
      if(ev_ve2.at(i).at(j).size()<1){continue;}
      I2.at(i).at(j)=2*pi*b*b*da_r*sin(ang_r2[j])
        /(Eave2[i]*Eave2[i]
          *sin(ang_r2[j]/2)*sin(ang_r2[j]/2)
          *sin(ang_r2[j]/2)*sin(ang_r2[j]/2));
      if(setbeam2){
        n_beam2.at(i).at(j)=n_b;
        setbeam2=false;
      }
      else{n_beam2.at(i).at(j)=n_b-Y2_sum;}
      if(n_beam2.at(i).at(j)<0){
        cout<<"                                                                 "<<endl;
        cout<<"end the loop at"<<l2[i]<<" mm~"<<l2[i]+dz<<" mm, "<<ang2[j]<<"deg~"<<ang2[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
        fin_beam2=true;
        break;
      }
      else{
        Y2_da.at(i).at(j)=I2.at(i).at(j)*n_beam2.at(i).at(j)*rho*dz*1e-3*N_rate2;
        Y2[i]+=Y2_da.at(i).at(j);
      }
      //      cout<<i<<"       "<<j<<"      "<<ev_ve.at(i).at(j).size()<<endl;
      if(i==i2_max && j==j2_max){
        cout<<"                                                                 "<<endl;
        cout<<"end the loop at "<<l2[i]<<" mm~"<<l2[i]+dz<<" mm "
          <<ang2[j]<<" deg~"<<ang2[j]+da<<" deg, because of no reaction in over "
          <<l2[i]+dz<<" mm "<<ang2[j]+da<<" deg"<<endl<<endl;
          //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
          fin_rea2=true;
          break;
      }
    }
    Y2_sum+=Y2[i];
    
    //z-reaction rate hist filling data
    if(n2[i]>0){
      h_z->Fill(l2[i],Y2[i]);
      h_z_p->Fill(l2[i],Y2[i]);
    }
    if (int(l2[i])%100==0){
      //need
      cout<<"12Cp finish "<<l2[i]<<" mm !!                              "<<endl;
    }
    if(fin_beam2){break;}
    if(fin_rea2){break;}
    //    if(z2_max<l2[i]+dz){}
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
        Y2_eve=Y2_da.at(i).at(j)/ev_ve2.at(i).at(j).size();
        rmax2=r_ve2.at(i).at(j).at(k);
        rmax2_nd=r_nd_ve2.at(i).at(j).at(k);
        h_rmax->Fill(rmax2,Y2_eve);
        h_rmax_p->Fill(rmax2,Y2_eve);
        h_rmax_nd->Fill(rmax2_nd,Y2_eve);
        n_rmax2++;
        n_rmax2_nd++;
        Y2_rmax+=Y2_eve;
        h_ztheta_p->Fill(dep2,dir2);
        rutp->Fill();
        /*
      	//for check 12cp
      	//	if(l2[i]==l_ext){
      	if(l2[i]==356){
        	//	  h_thetaY1->Fill(ang1[j],Y1_eve);
        	//	  h_thetaY1->Fill(ang1[j],1);
        	//	  if(checkang &&(ang1[j]<=20 || ang1[j]>160)){
          if(ang2[j]<=8){
        	  //	  if(checkang){
          	//	    cout<<l1[i]<<"   "<<ang1[j]<<"                              "<<endl;
          	//	    checkang=false;
          	//	  }}
          	h_thetaY1->Fill(dir1,Y1_eve);
          	h_thetaY1->Fill(dir1_an,Y1_eve);
          	h_thetaY2->Fill(ang1[j],Y1_eve);
          }
        }
      	*/
      }
    }
    //need
    cout<<"Filling p data "<<100*i/ev_ve2.size() <<" % !    \r"<<flush;
  }
  
  //make event lists  for all 12c12c, 12cp
    //fill event 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    z[0]=startz1[0]*ctom;
    z[1]=startz1[1]*ctom;
    scatter="12c12c";
    if(eve_temp>eve1){n_fi_temp++;}
    if(z_max<z[0]){z_max=z[0];}
    if(z_max<z[1]){z_max=z[1];}
    theta=theta1_cms[0];
    theta_an=theta1_cms[1];
    if(n_p1==2 && fabs(z[0]-z[1])<dz && theta>=theta_min && theta<=theta_max){
      z_n=int(z[0]/dz);
      a_n=int((theta-theta_min)/da);
      a_n_an=int(theta_an+0.99999);
      evsig[z_n][a_n]=1;
      ev_ve[z_n][a_n].push_back(eve1);
      evche[z_n][a_n].push_back(scatter);
      ang_an[z_n][a_n].push_back(a_n_an);
      E_ve[z_n][a_n].push_back(E1);
      r_ve[z_n][a_n].push_back(rmax1);
      r_nd_ve[z_n][a_n].push_back(rmax1_nd);
      fi_ve[z_n][a_n].push_back(n_fi_temp);
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;
    }
    else if(fabs(z[0]-z[1])<dz && theta>=theta_min && theta<=theta_max){
      scattap=scatter;
      n_enp++;
      n_efp=n_fi_temp;
      n_eep=eve1;
      n_epp=n_p1;
      expap->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p1==2 && theta>=theta_min && theta<=theta_max){
      scattal=scatter;
      n_enl++;
      n_efl=n_fi_temp;
      n_eel=eve1;
      n_ezl[0]=z[0];
      n_ezl[1]=z[1];
      expal->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p1==2 && fabs(z[0]-z[1])<dz){
      scattaa=scatter;
      n_ena++;
      n_efa=n_fi_temp;
      n_eea=eve1;
      n_eaa[0]=theta;
      n_eaa[1]=theta_an;
      expaa->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve1%2==0){
      scattab=scatter;
      n_enb++;
      n_efb=n_fi_temp;
      n_eeb=eve1;
      expab->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      scattao=scatter;
      n_eno++;
      n_efo=n_fi_temp;
      n_eeo=eve1;
      n_epo=n_p1;
      n_ezo[0]=z[0];
      n_ezo[1]=z[1];
      n_eao[0]=theta;
      n_eao[1]=theta_an;
      expao->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve_temp=eve1;
    if(j%1000==0){
      //need
      cout<<"Sorting elastic data(12c12c) "<<100*j/(Entry1+Entry2) <<" % !    \r"<<flush;
    }
  }
  i_max=int(z_max/dz);
  cout<<"about eliminated elastic(12c12c only) events"<<endl;
  cout<<" number of particles not 2 "<<n_enp<<"  events, not same start position "<<n_enl<<"  events,"<<endl;
  cout<<" angle not in the range "<<n_ena<<"  events, complex reason "<<n_eno<<"  events"<<endl<<endl;

  //make event lists  for all 12c12c, 12cp
    // fill event 12cp
  n_fi_temp=0;
  z_max=0;
  for(int j=0; j<Entry2; j++){
    tree2->GetEntry(j);
    z[0]=startz2[0]*ctom;
    z[1]=startz2[1]*ctom;
    scatter="12cp";
    if(eve_temp>eve2){n_fi_temp++;}
    if(z_max<z[0]){z_max=z[0];}
    if(z_max<z[1]){z_max=z[1];}
    theta=theta2_cms[0];
    theta_an=theta2_cms[1];
    if(n_p2==2 && fabs(z[0]-z[1])<dz && theta>=theta_min && theta<=theta_max){
      if(evsig[z_n][a_n]==0){
        //exist in first time (but only 12cp)
        evsig[z_n][a_n]=3;
      }
      else if(evsig[z_n][a_n]==1){
        //exist in first time (both events exist)
        evsig[z_n][a_n]=2;
      }
      else if(evsig[z_n][a_n]==2){
        //already exist (both events exist)
        //        cout<<"both 12c12c and 12cp events exist in "<<z[0]<<" mm "<<theta<<" deg "<<endl;
      }
      else if(evsig[z_n][a_n]==3){
        //already exist (but only 12cp)
        //        cout<<"12cp event only exist in "<<z[0]<<" mm "<<theta<<" deg "<<endl;
      }
      else {
        cout<<"something wrong with the event signature dep "<<z[0]<<" mm "<<theta<<" deg event signature "<<evsig[z_n][a_n]<<endl;
      }
      z_n=int(z[0]/dz);
      a_n=int((theta-theta_min)/da);
      a_n_an=int(theta_an+0.99999);
      ev_ve[z_n][a_n].push_back(eve2);
      evche[z_n][a_n].push_back(scatter);
      ang_an[z_n][a_n].push_back(a_n_an);
      E_ve[z_n][a_n].push_back(E2);
      r_ve[z_n][a_n].push_back(rmax2);
      r_nd_ve[z_n][a_n].push_back(rmax2_nd);
      fi_ve[z_n][a_n].push_back(n_fi_temp);
      //      cout<<j<<"  "<<n_ve.at(z_n1).at(a_n1).at(0)<<"   "<<E_ve1.at(z_n1).at(a_n1).at(0)<<endl;
    }
    else if(fabs(z[0]-z[1])<dz && theta>=theta_min && theta<=theta_max){
      scattap=scatter;
      n_enp++;
      n_efp=n_fi_temp;
      n_eep=eve2;
      n_epp=n_p2;
      expap->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p2==2 && theta>=theta_min && theta<=theta_max){
      scattal=scatter;
      n_efl=n_fi_temp;
      n_eel=eve2;
      n_ezl[0]=z[0];
      n_ezl[1]=z[1];
      expal->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p2==2 && fabs(z[0]-z[1])<dz){
      scattaa=scatter;
      n_efa=n_fi_temp;
      n_eea=eve2;
      n_eaa[0]=theta;
      n_eaa[1]=theta_an;
      expaa->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve2%2==0){
      scattab=scatter;
      n_enb++;
      n_efb=n_fi_temp;
      n_eeb=eve2;
      exppb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      scattao=scatter;
      n_eno++;
      n_efo=n_fi_temp;
      n_eeo=eve1;
      n_epo=n_p1;
      n_ezo[0]=z[0];
      n_ezo[1]=z[1];
      n_eao[0]=theta;
      n_eao[1]=theta_an;
      exppo->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    eve_temp=eve2;
    if(j%1000==0){
      //need
      cout<<"Sorting elastic data(12cp) "<<100*j/(Entry1+Entry2) <<" % !    \r"<<flush;
    }
  }
  if(i_max<int(z_max/dz)){
    i_max=int(z_max/dz);
  }
  cout<<"about eliminated elastic events"<<endl;
  cout<<" number of particles not 2 "<<n_enp<<"  events, not same start position "<<n_enl<<"  events,"<<endl;
  cout<<" angle not in the range "<<n_ena<<"  events, complex reason "<<n_eno<<"  events"<<endl;
  
  //elastic calculate
  // devide depth
  for(int i=0;i < ev_ve.size(); i++){
    l[i]=i*dz;
    n[i][0]=0;
    n[i][1]=0;
    Esum[i][0]=0;    
    Esum[i][1]=0;    
    Eave[i][0]=0;
    Eave[i][1]=0;
    Y[i]=0;
    
    //calculate Eave
    for(int j=0; j< ev_ve.at(i).size(); j++){
      if(i==i_max && ev_ve.at(i).at(j).size()>0){j_max=j;}
      for(int k=0; k< ev_ve.at(i).at(j).size(); k++){
        if(evche.at(i).at(j).at(k)=="12c12c"){
          n[i][0]++;
          Esum[i][0]+=E_ve.at(i).at(j).at(k);
        }
        else if(evche.at(i).at(j).at(k)=="12cp"){
          n[i][1]++;
          Esum[i][0]+=E_ve.at(i).at(j).at(k);
        }
        else{
          cout<<"                                                                 "<<endl;
          cout<<"something wrong with the vector of event checker "
            <<i<<"  "<<j<<" "<<k<<"  "<<evche.at(i).at(j).at(k)<<endl;
        }
      }
    }
    if(n[i][0]<1){continue;}
    else{Eave[i][0]=Esum[i][0]/(2*n[i][0]);}
    if(n[i][1]<1){continue;}
    else{Eave[i][1]=Esum[i][1]/(2*n[i][1]);}
    
    //calculate reaction rate
    // divide angle
    for(int j=0; j< ev_ve.at(i).size(); j++){
      ang[j]=j*da+theta_min;
      ang_r[j]=ang[j]*deg_to_rad;
      
      //need
      cout<<"elastic calculating reaction rate at "<<l[i]<<" mm "<<ang[j]<<" deg  \r"<<flush;
      
      if(ev_ve.at(i).at(j).size()<1){continue;}
      I.at(i).at(j).at(0)=2*pi*a*a*da_r*sin(ang_r[j])
        /(Eave[i][0]*Eave[i][0]
          *sin(ang_r[j]/2)*sin(ang_r[j]/2)
          *sin(ang_r[j]/2)*sin(ang_r[j]/2));

      I.at(i).at(j).at(1)=2*pi*b*b*da_r*sin(ang_r[j])
        /(Eave[i][1]*Eave[i][1]
          *sin(ang_r[j]/2)*sin(ang_r[j]/2)
          *sin(ang_r[j]/2)*sin(ang_r[j]/2));

      /*
      if(l1[i]<=5 && (ang1[j]<=20 || ang1[j]>=170)){
        cout<<"about I1  "<<l1[i]<<"   "<<ang1[j]<<"   "<<pi<<"     "<<a<<"     "<<da_r<<"     "<<Eave1[i]<<endl;
      }

    	if(l1[i]==10){
        Double_t tr=sin(ang_r1[j])
          /(sin(ang_r1[j]/2)*sin(ang_r1[j]/2)
            *sin(ang_r1[j]/2)*sin(ang_r1[j]/2));
          +sin(ang_r1[j])
          /(cos(ang_r1[j]/2)*cos(ang_r1[j]/2)
            *cos(ang_r1[j]/2)*cos(ang_r1[j]/2));
	
	      //	h_thetaY1->Fill(ang1[j],tr);
        h_thetaY1->Fill(ang1[j],I1.at(i).at(j)*1e30);
        //	h_thetaY2->Fill(180-ang1[j],I1.at(i).at(j));
        //	h_thetaY3->Fill(I1.at(i).at(j));
      }
      */
      
      if(setbeam){
        n_beam.at(i).at(j)=n_b;
        setbeam=false;
      }
      else{n_beam.at(i).at(j)=n_b-Y1_sum;}
      //      cout<<l1[i]<<"       "<<ang1[j]<<"      "<<n_beam1.at(i).at(j)<<"          "<<Y1_sum<<"       "<<endl;
      if(n_beam.at(i).at(j)<=0){
        cout<<"                                                                 "<<endl;
        cout<<"end the loop at"<<l[i]<<" mm~"<<l[i]+dz<<" mm, "<<ang[j]<<"deg~"<<ang[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
        fin_beam=true;
        break;
      }
      else{
        //check the existences of 12c12c and 12cp events
        if(evsig[z_n][a_n]==0){
          //no events
          Y_da.at(i).at(j).at(0)=0;
          Y_da.at(i).at(j).at(1)=0;
          //        cout<<"no event in "<<l[i]<<" mm "<<ang[j]<<" deg "<<endl;
        }
        else if(evsig[z_n][a_n]==1){
          //12c12c events only
          Y_da.at(i).at(j).at(0)=I.at(i).at(j).at(0)*n_beam.at(i).at(j)*rho*dz*1e-3*N_rate1;
          Y_da.at(i).at(j).at(1)=0;
          //        cout<<"12c12c event only exist in "<<l[i]<<" mm "<<ang[j]<<" deg "<<endl;
        }
        else if(evsig[z_n][a_n]==2){
          //both 12c12c and 12cp events
          Y_da.at(i).at(j).at(0)=I.at(i).at(j).at(0)*n_beam.at(i).at(j)*rho*dz*1e-3*N_rate1;
          Y_da.at(i).at(j).at(1)=I.at(i).at(j).at(1)*n_beam.at(i).at(j)*rho*dz*1e-3*N_rate2;
          //        cout<<"both 12c12c and 12cp events exist in "<<l[i]<<" mm "<<ang[j]<<" deg "<<endl;
        }
        else if(evsig[z_n][a_n]==3){
          //12cp events only
          Y_da.at(i).at(j).at(0)=0;
          Y_da.at(i).at(j).at(1)=I.at(i).at(j).at(1)*n_beam.at(i).at(j)*rho*dz*1e-3*N_rate2;
          //        cout<<"12cp event only exist in "<<l[i]<<" mm "<<ang[j]<<" deg "<<endl;
        }
        else {
          //no suggestion
          cout<<"something wrong with the event signature dep "<<l[i]<<" mm "<<ang[j]<<" deg event signature "<<evsig[z_n][a_n]<<endl;
        }

        Y[i]+=Y_da.at(i).at(j).at(0)+Y_da.at(i).at(j).at(1);
        //	Y1_sum+=Y1_da.at(i).at(j);
        //	cout<<i<<"   "<<k<<"          "<<Y1_da.at(i).at(k)<<"          "<<Y1[i]<<"          "<<Y1_sum<<endl;
        /*
        if(l1[i]<=5 && (ang1[j]<=20 || ang1[j]>=170)){
          cout<<"about Y1  "<<l1[i]<<"   "<<ang1[j]<<"   "<<rho<<"     "<<dz<<"     "<<N_rate1<<endl;
        }
        if(l[i]==l_ext && theta_min<=ang[j] && theta_max>=ang[j]){
          Double_t lY=log10(Y_da.at(i).at(j));
          //	  h_thetaY3->Fill(ang1[j],log10(Y1_eve));
          h_thetaY3->Fill(ang[j],lY);
          h_fit->Fill(ang[j],Y_da.at(i).at(j));
          //	  cout<<l1[i]<<"   "<<ang1[j]<<"    "<<lY1<<"                          "<<endl;
        }
        */
      }
      //      cout<<i<<"       "<<j<<"      "<<ev_ve.at(i).at(j).size()<<endl;
      if(i==i_max && j==j_max){
        cout<<"                                                                 "<<endl;
        cout<<"end the loop at "<<l[i]<<" mm~"<<l[i]+dz<<" mm "
          <<ang[j]<<" deg~"<<ang[j]+da<<" deg, because of no reaction in over "
          <<l[i]+dz<<" mm "<<ang[j]+da<<" deg"<<endl<<endl;
        //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
        fin_rea=true;
        break;
      }
    }
    Y_sum+=Y[i];
    
    //z-reaction rate hist filling data
    if(n[i][0]+n[i][1]>0){
      h_z->Fill(l[i],Y[i]);
      h_z_12c->Fill(l[i],Y[i]);
    }
    if (int(l[i])%100==0){
      //need
      cout<<"elastic calculation finish "<<l[i]<<" mm !!                              "<<endl;
    }
    if(fin_beam){break;}
    if(fin_rea){break;}
  }
  
  //make conections between rmax and reaction rate
  for(int i=0; i< ev_ve.size(); i++){
    l[i]=i*dz;
    for(int j=0;j < ev_ve.at(i).size() ; j++){
      ang[j]=j*da+theta_min;
      ang_r[j]=ang[j]*deg_to_rad;
      checkang=true;
      for(int k=0; k< ev_ve.at(i).at(j).size() ; k++){
        if(evche.at(i).at(j).at(k)=="12c12c"){
          scatter=evche.at(i).at(j).at(k);
          n_fi=fi_ve.at(i).at(j).at(k);
          n_ev=ev_ve.at(i).at(j).at(k);
          dep=l[i];
          dir=ang[j];
          dir_an=ang_an.at(i).at(j).at(k);
          Y_eve=Y_da.at(i).at(j).at(0)/ev_ve.at(i).at(j).size();
          rmax=r_ve.at(i).at(j).at(k);
          rmax_nd=r_nd_ve.at(i).at(j).at(k);
          h_rmax_all->Fill(rmax,Y_eve);
          h_rmax_a_nd->Fill(rmax_nd,Y_eve);
          n_rmax++;
          n_rmax_nd++;
          Y_rmax+=Y_eve;
          h_ztheta_all->Fill(dep,dir);
          rutal->Fill();
          /*
          //rmax test
          h_thetaY1->Fill(dep1,rmax1);
          if(rmax1>150){
            cout<<"strange event?  "<<n_fi1<<"   "<<n_ev1<<"   "<<dep1<<"   "<<dir1<<"   "<<dir1_an<<"   "<<rmax1<<endl;
          }

          //	if(n_fi1==0 && rmax1>150){
          if(rmax1>150){
            cout<<"file "<<n_fi1<<"   event "<<n_ev1<<"  depth "<<dep1<<" mm  angle "<<dir1<<" deg   event is rmax "<<rmax1<<" mm "<<endl;
          }

          //for check 12c12c
          if(l1[i]==l_ext){
            //	  h_thetaY1->Fill(ang1[j],Y1_eve);
            //	  h_thetaY1->Fill(ang1[j],1);
            //	  if(checkang &&(ang1[j]<=20 || ang1[j]>160)){
            //	  if(ang1[j]<=8){
            //	  if(checkang){
            //	    cout<<l1[i]<<"   "<<ang1[j]<<"                              "<<endl;
            //	    checkang=false;
            //	  }
            h_thetaY1->Fill(dir1,Y1_eve);
            h_thetaY1->Fill(dir1_an,Y1_eve);
            h_thetaY2->Fill(ang1[j],Y1_eve);
          }
          if(l1[i]==142){
            h_thetaY2->Fill(ang1[j],Y1_eve);
            //	  h_thetaY2->Fill(180-ang1[j],Y1_eve);
          }
          if(l1[i]==5){
            h_thetaY3->Fill(ang1[j],Y1_eve);
            h_thetaY3->Fill(180-ang1[j],Y1_eve);
          }
          */
        }
        else if(evche.at(i).at(j).at(k)=="12cp"){
          scatter=evche.at(i).at(j).at(k);
          n_fi=fi_ve.at(i).at(j).at(k);
          n_ev=ev_ve.at(i).at(j).at(k);
          dep=l[i];
          dir=ang[j];
          Y_eve=Y_da.at(i).at(j).at(1)/ev_ve.at(i).at(j).size();
          rmax=r_ve.at(i).at(j).at(k);
          rmax_nd=r_nd_ve.at(i).at(j).at(k);
          h_rmax_all->Fill(rmax,Y_eve);
          h_rmax_a_nd->Fill(rmax_nd,Y_eve);
          n_rmax++;
          n_rmax_nd++;
          Y_rmax+=Y_eve;
          h_ztheta_all->Fill(dep,dir);
          rutal->Fill();

          //for check 12cp
          /*
          //	if(l2[i]==l_ext){
          if(l2[i]==356){
            //	  h_thetaY1->Fill(ang1[j],Y1_eve);
            //	  h_thetaY1->Fill(ang1[j],1);
            //	  if(checkang &&(ang1[j]<=20 || ang1[j]>160)){
            if(ang2[j]<=8){
            	//	  if(checkang){
              //	    cout<<l1[i]<<"   "<<ang1[j]<<"                              "<<endl;
              //	    checkang=false;
            }}
              h_thetaY1->Fill(dir1,Y1_eve);
              h_thetaY1->Fill(dir1_an,Y1_eve);
              h_thetaY2->Fill(ang1[j],Y1_eve);
          }
          */
        }
        else{
          cout<<"something wrong with the vector of event checker "
            <<i<<"  "<<j<<" "<<k<<"  "<<evche.at(i).at(j).at(k)<<endl;
        }
      }
    }
    //need
    cout<<"Filling elastic data "<<100*i/ev_ve.size() <<" % !    \r"<<flush;
  }

#ifdef RES_ON
  cout<<"                                                                 "<<endl;
  cout<<"12c12c n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax
    <<"  n_rmax1_nd  "<<n_rmax1_nd<<"                   "<<endl;
  cout<<"12cp   n_rmax2 "<<n_rmax2<<"    Y2sum "<<Y2_rmax
    <<"  n_rmax2_nd  "<<n_rmax2_nd<<"                   "<<endl;
  cout<<"elastic   n_rmax "<<n_rmax<<"    Ysum "<<Y_rmax
    <<"  n_rmax_nd  "<<n_rmax_nd<<"                   "<<endl;
  cout<<endl;
  
  //search rmax when over 100Hz
  double hist_cnt = 0;
  int flag_100hz=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt += h_rmax->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt > 100 && flag_100hz==0 ){
      printf("sum of 12c12c and 12cp over 100Hz %d mm \n\n", ihist);
      flag_100hz=1;
    }
  }
  double hist_cnt_nd = 0;
  int flag_nd_100hz=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt_nd += h_rmax_nd->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt_nd > 100 && flag_nd_100hz==0 ){
      printf("nd  over 100Hz %d mm \n\n", ihist);
      flag_nd_100hz=1;
    }
  }
  double hist_cnt_12c = 0;
  int flag_100hz_12c=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt_12c += h_rmax_12c->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax_12c->GetBinContent(nbin));
    if(hist_cnt_12c > 100 && flag_100hz_12c==0 ){
      printf("12c over 100Hz %d mm \n\n", ihist);
      flag_100hz_12c=1;
    }
  }
  double hist_cnt_p = 0;
  int flag_100hz_p=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt_p += h_rmax_p->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax_12c->GetBinContent(nbin));
    if(hist_cnt_p > 100 && flag_100hz_p==0 ){
      printf("p over 100Hz %d mm \n\n", ihist);
      flag_100hz_p=1;
    }
  }
  double hist_cnt_all = 0;
  int flag_100hz_all=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt_all += h_rmax_all->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt_all > 100 && flag_100hz_all==0 ){
      printf("elastic (12c12c and 12cp) over 100Hz %d mm \n\n", ihist);
      flag_100hz_all=1;
    }
  }
#endif
  //need
#ifdef CAN_ON
  //make error bar
  // for 12c12c and 12cp
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z->GetBinContent(m));
    h_z->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax->GetBinContent(m));
    h_rmax->SetBinError(m,err);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax_nd->GetBinContent(m));
    h_rmax_nd->SetBinError(m,err);
  }
  // for 12c12c
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta_12c->GetBinContent(m));
    h_ztheta_12c->SetBinError(m,err);
  }
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z_12c->GetBinContent(m));
    h_z_12c->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax_12c->GetBinContent(m));
    h_rmax_12c->SetBinError(m,err);
  }
  // for 12cp
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta_p->GetBinContent(m));
    h_ztheta_p->SetBinError(m,err);
  }
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z_p->GetBinContent(m));
    h_z_p->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax_p->GetBinContent(m));
    h_rmax_p->SetBinError(m,err);
  }

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,4);
  
  // z-reaction rate all
  h_z->SetStats(0);
  //  h_z->SetStats(1);
  can->cd(1);
  gPad->SetLogy(1);
  h_z->Draw("hist");

  // rmax-reaction rate all
  //  h_rmax->SetStats(0);
  h_rmax->SetStats(1);
  can->cd(2);
  gPad->SetLogy(1);
  h_rmax->Draw("hist");

  // rmax-reaction rate all
  //  h_rmax->SetStats(0);
  h_rmax_nd->SetStats(1);
  can->cd(3);
  gPad->SetLogy(1);
  h_rmax_nd->Draw("hist");

  // z-angle only 12c
  //  h_ztheta_12c->SetStats(0);
  h_ztheta_12c->SetStats(1);
  can->cd(4);
  h_ztheta_12c->Draw("colz");

  // z-reaction rate only 12c
  h_z_12c->SetStats(0);
  //  h_z_12c->SetStats(1);
  can->cd(5);
  gPad->SetLogy(1);
  h_z_12c->Draw("hist");

  // rmax-reaction rate only 12c
  //  h_rmax_12c->SetStats(0);
  h_rmax_12c->SetStats(1);
  can->cd(6);
  gPad->SetLogy(1);
  h_rmax_12c->Draw("hist");

  // z-angle only p
  //  h_ztheta_p->SetStats(0);
  h_ztheta_p->SetStats(1);
  can->cd(7);
  h_ztheta_p->Draw("colz");

  // z-reaction rate only p
  h_z_p->SetStats(0);
  //  h_z_p->SetStats(1);
  can->cd(8);
  gPad->SetLogy(1);
  h_z_p->Draw("hist");

  // rmax-reaction rate only p
  //  h_rmax_p->SetStats(0);
  h_rmax_p->SetStats(1);
  can->cd(9);
  gPad->SetLogy(1);
  h_rmax_p->Draw("hist");

  //output for canvas
  TString outcan = "data1/ruther_can_nd_0925_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  h_z->Write();
  h_rmax->Write();
  h_rmax_nd->Write();
  h_ztheta_12c->Write();
  h_z_12c->Write();
  h_rmax_12c->Write();
  h_ztheta_p->Write();
  h_z_p->Write();
  h_rmax_p->Write();
#endif
  //need to here
  /*koko
  // for check 12c12c
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_thetaY1->GetBinContent(m));
    h_thetaY1->SetBinError(m,err);
  }kokomade

  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_thetaY2->GetBinContent(m));
    h_thetaY2->SetBinError(m,err);
  }
  for(int m=0; m<=nbin;m++){
    //    double err=sqrt(h_thetaY3->GetBinContent(m));
    double err=1e-5;
    h_thetaY3->SetBinError(m,err);
  }
  */
  /*koko
  TCanvas *can1=new TCanvas("hists","hists");
  //  can1->Divide(2,2);
  kokomade
  // theta-reaction rate
  //  can1->cd(1);
  gPad->SetLogy(0);
  h_thetaY1->SetStats(0);

  f1->SetParameters(1,deg_to_rad);
  f1->SetParLimits(0,0,2);
  f1->SetParLimits(1,deg_to_rad,deg_to_rad);
  //  h_thetaY1->Fit("f1","","",0,180);
  //  h_thetaY1->Fit("f1","","",theta_min,180);
  //  h_thetaY1->Fit("f1","","",theta_min,181-theta_min);
  //  h_thetaY1->Fit("f1","","",theta_min,theta_max);
  h_thetaY1->Fit("f1","","",10,170);
  */
  /*koko
  h_thetaY1->Draw("colz");
  */
  /*
  f1->Draw("same");
  Double_t chi1=f1->GetChisquare();
  Double_t ndf1=f1->GetNDF();
  cout<<chi1/ndf1<<endl;

  can1->cd(2);
  gPad->SetLogy(1);
  h_thetaY2->SetStats(0);
  f2->SetParameters(1,deg_to_rad);
  f2->SetParLimits(0,0,2);
  f2->SetParLimits(1,deg_to_rad,deg_to_rad);
  //  h_thetaY2->Fit("f2","","",0,180);
  //  h_thetaY2->Fit("f2","","",theta_min,180);
  //  h_thetaY2->Fit("f2","","",theta_min,181-theta_min);
  h_thetaY2->Fit("f2","","",theta_min,theta_max);
  //  h_thetaY2->Fit("f2","","",10,170);
  h_thetaY2->Draw("hist");
  f2->Draw("same");
  Double_t chi2=f2->GetChisquare();
  Double_t ndf2=f2->GetNDF();
  cout<<chi2/ndf2<<endl;

  can1->cd(3);
  gPad->SetLogy(0);
  h_thetaY3->SetStats(0);
  f3->SetParameters(1,deg_to_rad);
  f3->SetParLimits(0,0,2);
  f3->SetParLimits(1,deg_to_rad,deg_to_rad);
  //  h_thetaY3->Fit("f3","","",theta_min,180);
  h_thetaY3->Fit("f3","","",theta_min,theta_max);
  //  h_thetaY3->Fit("f3","","",10,170);
  h_thetaY3->Draw("hist");
  f3->Draw("same");
  Double_t chi3=f3->GetChisquare();
  Double_t ndf3=f3->GetNDF();
  cout<<chi3/ndf3<<endl;
  */
  /*
  for(int i=0; i< ev_ve1.size(); i++){
    l1[i]=i*dz;
    if(l1[i]==l_ext){
      for(int j=0;j < ev_ve1.at(i).size() ; j++){
	ang1[j]=j*da+theta_min;
	y_fit=f3->Eval(ang1[j]);
	ra_fit=log10(Y1_da.at(i).at(j))/y_fit;
	//	ra_fit=Y1_da.at(i).at(j);
	h_fit->Fill(ang1[j],ra_fit);
      }
    }
  }  
  can1->cd(4);
  gPad->SetLogy(1);
  h_fit->SetStats(0);
  h_fit->Draw("hist");
  */

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
  rutal->AutoSave();
  expap->AutoSave();
  expal->AutoSave();
  expaa->AutoSave();
  expav->AutoSave();
  expab->AutoSave();
  expao->AutoSave();
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
