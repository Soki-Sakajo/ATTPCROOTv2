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

int ruther_macro_c12_each_angle(){ 
  
  //This file is only 12C12C and 12C1H elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t z_s=6; // this is Z of scatterd particl in simfile
  Int_t z_r1=6; // this is Z of an recoil particl in simfile
  Int_t z_r2=1; // this is Z of another recoil particl in simfile
  Int_t dz=1; //mm; this is thickness you think
  Int_t n_b=1e5; // /s; beam intensity
  Double_t T=300; //K; thermodynamic temparature
  Double_t P=0.13; //atm; gas pressure
  Double_t dr=1;//mm; bin of reaction rate
  Int_t nbin=300;//bin of hist
  Int_t N_rate1=4;//(number of scatterd particle)/molecular
  Int_t N_rate2=10;//(number of recoli partcle)/molecular
  
  Int_t theta_min=3; //deg; rutherford scattar angle
  Double_t dtheta=1;//dig; bin of angle
  Int_t r_attpc=275; //mm
  Int_t l_attpc=1000; //drift length
  
  //set constants
  int ctom=10;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; //MPa
  Double_t alfa=7.297352569e-3;//fine structure constant
  Double_t hc=197.463587*1e-15;//MeV•m
  //  Double_t hc=197.463587;//MeV•fm
  Double_t k_b=1.38065e-23; //J/K ; Bolzman constant
  Double_t dtheta_r=dtheta*deg_to_rad;
  Double_t rho=P*atmtoPa/(k_b*T);
  Double_t I=-2/(sin((180-theta_min)*deg_to_rad/2)*sin((180-theta_min)*deg_to_rad/2))+2/(sin(theta_min*deg_to_rad/2)*sin(theta_min*deg_to_rad/2));
  //  cout<<I<<endl;
  
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  Double_t b=alfa*hc*z_s*z_r2/4;  
  
  //set variable
  Double_t rmax1;
  vector<Double_t> *r1=0;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t theta1_cms_r[2];
  Double_t rmax2;
  vector<Double_t> *r2=0;
  Double_t E2;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];
  Double_t theta2_cms_r[2];
  vector<Double_t> r;

  Int_t n_p1;
  Int_t n_p2;
  Int_t n_b1=0;
  Int_t n_b2=0;
  Int_t evec;
  Double_t theta_s;
  Double_t rmax_s;
  Double_t Y_s;

  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/dtheta+0.9);
  Int_t n_rmax1[div_l]; //
  Int_t n_rmax2[div_l]; //  
  Int_t n_rmax[div_l]; //
  Double_t l1[div_l];
  Double_t l2[div_l];
  Double_t theta_r[div_a];
  Double_t theta_r1[div_a];
  Double_t theta_r2[div_a];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  //  Double_t Y1[div_a];
  Double_t Y1[div_l];
  //  Double_t Y1_ang[div_a];
  Double_t Y1_ang;
  Double_t n_react[div_l];
  Double_t n_beam[div_l];
  
  Double_t n2[div_l];
  Double_t Esum2[div_l];
  Double_t Eave2[div_l];
  Double_t Y2[div_a];

  Double_t sigma[div_l];
  Double_t Y1_sum;
  
  //input
  TFile *file1 =new TFile("data/sim_macro_c12.root"); 
  TTree *tree1 = (TTree*)file1->Get("anadata");
  
  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("startz",startz1);
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);
  //  tree1->SetBranchAddress("r",&r1);
  
  //set output
  /*
  TFile *rutmac=new TFile("data/ruther_macro_c12.root","recreate");
  TTree *rutdat=new TTree("rutdat","ratdat");
  rutdat->Branch("eventcode",&evec,"evec/I");
  rutdat->Branch("theta",&theta_s,"theta_s/D");
  rutdat->Branch("rmax",&rmax_s,"rmax_s/D");
  rutdat->Branch("r",&r);
  rutdat->Branch("reaction rate",&Y_s,"Y_s/D");
  */
  
  TH1F *h_rmax=new TH1F("h_rmax","r:reaction_rate",nbin,0,300);
  //  TH1F *h_r=new TH1F("h_r","r:reaction_rate",nbin,0,150);

  /*  
  //check beam event
  for(int m=0;m<Entry1;m++){
    tree1->GetEntry(m);
    tree1->SetBranchStatus("*",0);
    tree1->SetBranchStatus("n",1);

    if(m%2==0 && n_p1==1){
      n_b1++;
    }
    //      if(m%100==0){
    cout<<"check beam event "<<100*m/Entry1<<" %\r"<<flush;
    //      } 
  }
  //  cout<<n_b1<<" beam events (12c12c) "<<endl;
  n_b=n_b1;
  cout<<n_b<<" beam events     "<<endl;
  */
  
  //12c12c calculate
  for(int i=0;i < div_l ; i++){
  //  for(int i=0;i <100; i++){
    l1[i]=i*dz;
    l2[i]=(i+1)*dz;

    // initialization
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;

    for(int j=0; j< Entry1; j++){
      cout<<"Reading events "<<100*j/Entry1 <<" % !\r"<<flush;

      tree1->GetEntry(j);
      tree1->SetBranchStatus("*",0);
      tree1->SetBranchStatus("startz",1);
      tree1->SetBranchStatus("TotalE",1);
      z1[0]=startz1[0]*ctom;
      z1[1]=startz1[1]*ctom;

      if(j%2==1 && fabs(z1[0]-z1[1])<dz && z1[0]>=l1[i] && z1[0]<l2[i]){
	    n1[i]++;
	    Esum1[i]+=E1;    
      }
      else if(j%2==0){
	//	cout<<i<<" "<<j<<" This is beam event"<<endl;
	continue;
      }
      else {
	continue;
      }
    }
    if(n1[i]<1){
      //      Eave1[i]=0;
      continue;
    }
    else{
      Eave1[i]=Esum1[i]/n1[i];
      //      cout<<i<<"     "<<Eave1[i]<<"            "<<endl;
    }

    Y1[i]=0;
    n_react[i]=0;

    if(i==0){n_beam[i]=n_b;}
    else{n_beam[i]=n_b-Y1_sum;}

    if(n1[i]<2){
      continue;
    }
    else{
      
      sigma[i]=2*pi*a*a*I/(Eave1[i]*Eave1[i]);
      //      sigma[i]=2*pi*a*a*dtheta_r*sin(theta_r[k])/(Eave1[i]*Eave1[i]*sin(theta_r[k]/2)*sin(theta_r[k]/2)*sin(theta_r[k]/2)*sin(theta_r[k]/2));
      
      //	Y1[i]=2*pi*a*a*(n_b-n_react[i-1])*dtheta_r*I/(Eave1[i]*Eave1[i]);
      Y1[i]=sigma[i]*n_beam[i]*rho*dz*1e-3*N_rate1;
      //            Y1[i]=sigma[i]*n_beam[i]*rho*dz*1e-3*N_rate1*(n1[i]/Entry1);

      Y1_sum+=Y1[i];
    }
    cout<<i<<"     "<<n_beam[i]<<"         "<<sigma[i]<<"      "<<Y1[i]<<"         "<<Y1_sum<<endl;

    for(int k=0; k<div_a; k++){
      theta_r1[k]=(k*dtheta+theta_min)*deg_to_rad;
      theta_r[k]=theta_r1[k]+dtheta_r/2;
      theta_r2[k]=theta_r1[k]+dtheta_r;
      
      for(int l=0; l<Entry1; l++){
	tree1->GetEntry(l);
	tree1->SetBranchStatus("*",0);
	tree1->SetBranchStatus("startz",1);
	tree1->SetBranchStatus("theta_cms",1);
	tree1->SetBranchStatus("rmax",1);

	z1[0]=startz1[0]*ctom;
	z1[1]=startz1[1]*ctom;
	theta1_cms_r[0]=theta1_cms[0]*deg_to_rad;
	theta1_cms_r[1]=theta1_cms[1]*deg_to_rad;

	Y1_ang=Y1[i]*dtheta_r*sin(theta_r[k])
	  /(I*Eave1[i]*Eave1[i]
	    *sin(theta_r[k]/2)*sin(theta_r[k]/2)*sin(theta_r[k]/2)*sin(theta_r[k]/2));

	if(l%2==1 && fabs(z1[0]-z1[1])<dz && z1[0]>=l1[i] && z1[0]<l2[i] &&
	   fabs(theta1_cms_r[0]-theta_r[k])<dtheta_r/2){
	  h_rmax->Fill(rmax1,Y1_ang);
	  //	  h_rmax->Fill(rmax1/n_rmax[i],Y1[k]);
	  //	    h_r->Fill(r1,Y2[k]); //if use, need to modify to vector
	}
	else if (l%2==0){
	  h_rmax -> Fill(rmax1, 1.0e5/n_b);
	}
	else{continue;}
	//	if(l%1000==0) {
	cout<<"Entry "<<100*l/Entry1 <<"  %  angle "<<100*k/div_a<<"%  Filling data finished "<<100*l/Entry1 <<" % !\r"<<flush;
	//	}
      }
	
      //    h_rmax->Fill(l1[i],Y1_ang[k]);
    }
    if (i%100==0){
      cout<<"12C12C finish "<<100*i/div_l<<"% !!                        "<<endl;
      //      cout<<"finish "<<100*(i+div_l)/(2*div_l)<<"% !!"<<endl;
    }
  }
  
  for(int m=0; m<=nbin; m++){
    double err = sqrt(h_rmax->GetBinContent(m));
    h_rmax->SetBinError(m, err);  
  }
  h_rmax->Draw("hist");
  //  h_rmax->SetLogy();
  //  h_rmax->GetYaxis()->SetMoreLogLabels();
  h_rmax->SaveAs("data/ruther_macro_result_c12_per_angle.C");
  //  rutdat->AutoSave();
  //  rutmac->Close();
  
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
