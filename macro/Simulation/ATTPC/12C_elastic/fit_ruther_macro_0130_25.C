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
#include<cmath>

//int fit_ruther_macro_0130_25(Int_t gas=0){
int fit_ruther_macro_0130_25(Int_t gas=130){

  if(gas==70 ||gas==80 ||gas==90 ||gas==100 ||gas==110 ||gas==120 ||gas==130 ){
    cout<<endl;
    cout<<"this is calculation of gas pressure "<<gas<<" matm"<<endl;
  }
  else{
    cout<<"chose "<<gas<<" atm. it is no supported"<<endl;
    cout<<"need to choose from 70--130."<<endl;
    gSystem->Exit(0);
  }
  Int_t subnum=30;
  
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
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t N_rate2=10;     // (number of recoli partcle)/molecular  
  Int_t theta_min=3;    // deg; rutherford scattar angle
  Int_t theta_max=177;    // deg; rutherford scattar angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc

  Int_t l_ext=100;   // mm; extract information
  Double_t ra_fit=-100;
  Double_t y_fit=-100;
  Double_t fit_st=50;


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
  Int_t eve1;
  Int_t n_p1;
  Int_t n_fi1;
  Int_t n_ev1;
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
  //  Int_t pdg2[2];

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
  
  //input multi-file by using TChain
  //common set for input
  //  TString dir="data1";
  TString dir="data2";
  
  // Tchain of 12c12c
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    TString refi="sim_macro_12c_"+to_string(gas)+"matm_" + to_string(i) + ".root";
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
    TString refi="sim_macro_p_"+to_string(gas)+"matm_" + to_string(i) + ".root";
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
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("rmax_nd",&rmax1_nd);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);
  //  for 12C12p
  Int_t Entry2=tree2->GetEntries();
  //  tree2->SetBranchAddress("PdgCode",&pdg2);
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

  TH2F *h_ztheta_12c=new TH2F("h_ztheta_12c","reaction position:angle (12c12c)",nbin_zn,0,400,100,0,200);
  TH1F *h_z_12c=new TH1F("h_z_12c","reaction position:reaction_rate (12c12c)",nbin_zn,0,400);
  TH1F *h_rmax_12c=new TH1F("h_rmax_12c","rmax:reaction_rate (12c12c)",nbin,0,275);
  TH2F *h_ztheta_p=new TH2F("h_ztheta_p","reaction position:angle (12cp)",nbin_zn,0,400,100,0,200);
  TH1F *h_z_p=new TH1F("h_z_p","reaction position:reaction_rate (12cp)",nbin_zn,0,400);
  TH1F *h_rmax_p=new TH1F("h_rmax_p","rmax:reaction_rate (12cp)",nbin,0,275);

  //  TH1F *h_thetaY1=new TH1F("h_thetaY1","angle :reaction_rate",180,0,180);
  TH2F *h_thetaY1=new TH2F("h_thetaY1","reaction position :rmax",nbin_zn,0,400,nbin,0,275);
  TH1F *h_thetaY2=new TH1F("h_thetaY2","angle :reaction_rate",180,0,180);
  TH1F *h_thetaY3=new TH1F("h_thetaY3","angle :reaction_rate",180,0,180);
  TH1F *h_fit=new TH1F("h_fit","angle :ratio of number of reaction and fit",180,0,180);

  //  TF1 *f1=new TF1("f1","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2))+1/(cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)))",0,180);
  //  TF1 *f2=new TF1("f2","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)))",0,180);
  TF1 *f1=new TF1("f1","[0]*x^[1]",0,300);
  TF1 *f2=new TF1("f2","[0]*x^[1]+[2]*x^[3]",0,300);
  //  TF1 *f2=new TF1("f2","[2]*x^[3]+[4]*x^[5]",0,300);
  //  TF1 *f1=new TF1("f1","[0]*x^[1]+[2]*x^[3]+[4]*x^[5]",0,300);

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
    eve1_temp=eve1;
    if(j%1000==0){
      //need
      cout<<"Sorting 12c12c data "<<100*j/Entry1 <<" % !    \r"<<flush;
    }
  }
  i1_max=int(z1_max/dz);
  
  //12c12c calculate
  // devide depth
  for(int i=0;i < ev_ve1.size(); i++){
    l1[i]=i*dz;
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;
    
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
    Y1[i]=0;
    
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
      //      cout<<l1[i]<<"       "<<ang1[j]<<"      "<<n_beam1.at(i).at(j)<<"          "<<Y1_sum<<"       "<<endl;
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
    eve2_temp=eve2;
    if(j%1000==0){
      //need
      cout<<"Sorting 12cp data "<<100*j/Entry2 <<" % !    \r"<<flush;
    }
  }
  i2_max=int(z2_max/dz);
  
  //12cp calculate
  // devide depth
  for(int i=0;i < ev_ve2.size(); i++){
    l2[i]=i*dz;
    n2[i]=0;
    Esum2[i]=0;    
    Eave2[i]=0;
    
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
    Y2[i]=0;
    
    //divide angle
    for(int j=0; j< ev_ve2.at(i).size(); j++){
      ang2[j]=j*da+theta_min;
      ang_r2[j]=ang2[j]*deg_to_rad;
      
      //need
      cout<<"p calculating reaction rate at "<<l2[i]<<" mm "<<ang2[j]<<" deg  \r"<<flush;
      
      if(ev_ve2.at(i).at(j).size()<1){continue;}
      
      I2.at(i).at(j)=
	2*pi*b*b*da_r*sin(ang_r2[j])
	/(Eave2[i]*Eave2[i]
	  *sin(ang_r2[j]/2)*sin(ang_r2[j]/2)
	  *sin(ang_r2[j]/2)*sin(ang_r2[j]/2));
      if(setbeam2){
	n_beam2.at(i).at(j)=n_b;
	setbeam2=false;
      }
      else{n_beam2.at(i).at(j)=n_b-Y2_sum;}
      //      cout<<i<<"       "<<j<<"      "<<n_beam1.at(i).at(j)<<endl;
      if(n_beam2.at(i).at(j)<0){
	cout<<"                                                                 "<<endl;
	cout<<"end the loop at"<<l2[i]<<" mm~"<<l2[i]+dz<<" mm, "<<ang2[j]<<"deg~"<<ang2[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
	fin_beam2=true;
	break;
      }
      else{
	Y2_da.at(i).at(j)=I2.at(i).at(j)*n_beam2.at(i).at(j)*rho*dz*1e-3*N_rate2;
	Y2[i]+=Y2_da.at(i).at(j);
	//	cout<<i<<"   "<<k<<"          "<<Y1_da.at(i).at(k)<<"          "<<Y1[i]<<"          "<<Y1_sum<<endl;
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
      }
    }
    //need
    cout<<"Filling p data "<<100*i/ev_ve2.size() <<" % !    \r"<<flush;
  }
  
  cout<<"                              "<<endl;
  cout<<"12c12c n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax<<"  n_rmax1_nd  "<<n_rmax1_nd<<"                   "<<endl;
  
  cout<<"12cp   n_rmax2 "<<n_rmax2<<"    Y2sum "<<Y2_rmax<<"  n_rmax2_nd  "<<n_rmax2_nd<<"                   "<<endl;
  cout<<endl;
  
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
  can->Divide(3,1);
  /*

  // z-reaction rate all
  h_z->SetStats(0);
  //  h_z->SetStats(1);
  can->cd(1);
  gPad->SetLogy(1);
  h_z->Draw("hist");
  */
  // rmax-reaction rate all
  h_rmax->SetStats(0);
  //  h_rmax->SetStats(1);
  can->cd(1);
  gPad->SetLogy(1);
  //  h_rmax->GetXaxis()->SetRangeUser(0,275);
  h_rmax->Draw("hist");

  //  can->cd(3);
  can->cd(2);
  gPad->SetLogy(1);
  //  gPad->SetLogx(1);
  //  f1->SetParameters(1e13,-10,-1e13,-10);
  //  f1->SetParLimits(0,9.96e12,9.97e12);
  //  f1->SetParLimits(1,-7.11,-7.10);
  h_rmax->SetStats(0);
  h_rmax->Fit("f1","","",fit_st,146);
  //  h_rmax->Fit("f1","","",26,146);
  //  h_rmax->GetXaxis()->SetRangeUser(50,146);
  h_rmax->Draw("hist");
  f1->Draw("same");
  Double_t chi1=f1->GetChisquare();
  Double_t ndf1=f1->GetNDF();
  cout<<endl;
  cout<<"Chi/NDF: "<<chi1/ndf1<<endl;

  Double_t p1=f1->GetParameter(0);
  Double_t p2=f1->GetParameter(1);

  //  can->cd(4);
  can->cd(3);
  gPad->SetLogy(1);
  //  gPad->SetLogx(1);
  f2->SetParameter(0,p1);
  f2->SetParameter(1,p2);
  f2->FixParameter(0,p1);
  f2->FixParameter(1,p2);
  //  f2->SetParLimits(0,p1-1e-6,p1+1e-6);
  //  f2->SetParLimits(1,p2-1e-6,p2+1e-6);
  h_rmax->SetStats(0);
  h_rmax->Fit("f2","","",fit_st,146);
  //  h_rmax->Fit("f2","","",26,146);
  //  h_rmax->GetXaxis()->SetRangeUser(50,146);
  //  gPad->SetLogx(1);
  h_rmax->Draw("hist");
  f2->Draw("same");

  TCanvas *can2=new TCanvas("hists2","hists2",1500,1000);
  can2->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.13);
  h_rmax->SetStats(0);
  h_rmax->SetTitle("r_{max} distribution (Rutherford scattering 12c12c and 12cp)");
  h_rmax->GetXaxis()->SetTitle("r_{max} [mm]");
  h_rmax->GetYaxis()->SetTitle("reaction rate [cps]");
  h_rmax->Draw("hist");
  f2->Draw("same");
  h_rmax->SetTitleSize(0.08);
  h_rmax->GetXaxis()->SetLabelSize(0.06);
  h_rmax->GetYaxis()->SetLabelSize(0.06);
  h_rmax->GetXaxis()->SetTitleSize(0.06);
  h_rmax->GetYaxis()->SetTitleSize(0.06);
  h_rmax->GetYaxis()->SetTitleOffset(1.1);
  //  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.06);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(155,3.5e3,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  
  //output for canvas
  //  TString outcan = dir+"/fit_ruther_can_0130_130matm_" + to_string(subnum) + ".C";
  TString outcan = "test_rut.pdf";
  //  TString outcan2 = "test.pdf";
  //  TString outcan2 = "./pdf_output_fit_ruther_can_0130_"+to_string(gas)+"matm_" + to_string(subnum) + ".pdf";
  TString outcan2 = "pdf_output/fit_ruther_can_0130_"+to_string(gas)+"matm_" + to_string(subnum) + ".pdf";
  can->SaveAs(outcan);
  can2->SaveAs(outcan2);
  cout<<endl;


  Double_t p3=f2->GetParameter(2);
  Double_t p4=f2->GetParameter(3);

  Int_t flag_100hz_fit=0;
  Double_t hist_cnt_fit = 0;
  Double_t r_nbin=0;
  Double_t fx=0;
  Double_t val_bin=0;
  for(int ihist=nbin; ihist>0; ihist--){
    r_nbin=dr*ihist;
    
    if(r_nbin<fit_st || 146<r_nbin){
      hist_cnt_fit += h_rmax->GetBinContent(ihist);
    }
    else if(fit_st<=r_nbin && r_nbin<=146){
      //      fx=p1*pow(r_nbin,p2);
      fx=p1*pow(r_nbin,p2)+p3*pow(r_nbin,p4);
      hist_cnt_fit += fx;
      val_bin= h_rmax->GetBinContent(ihist);
      //      cout<<"ihist: "<<ihist<<" val of bin "<<val_bin<<" fx "<<fx<<endl;
    }
    else {
      cout<<endl;
      cout<<"Something worng with count over 100Hz!! end this macro."<<endl;
      cout<<"ihist: "<<ihist<<" r_nbin: "<<r_nbin<<endl;

      gSystem->Exit(0);
    }
    if(hist_cnt_fit > 100 && flag_100hz_fit==0 ){
      cout<<endl;
      cout<<"this is result of "<<gas<<" matm data"<<endl;
      printf("all over 100Hz from fit %d mm \n\n", ihist);
      flag_100hz_fit=1;
    }
  }

  double hist_cnt = 0;
  int flag_100hz=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt += h_rmax->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt > 100 && flag_100hz==0 ){
      printf("all over 100Hz %d mm \n\n", ihist);
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
      printf("p   over 100Hz %d mm \n\n", ihist);
      flag_100hz_p=1;
    }
  }

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
