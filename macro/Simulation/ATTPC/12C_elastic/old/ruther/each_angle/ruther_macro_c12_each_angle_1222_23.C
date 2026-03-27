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

int ruther_macro_c12_each_angle_1222(){ 
  
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
  Int_t nbin=300/dr;//bin of rmax-reaction_rate hist
  Int_t nbin_zn=400/dz;//bin of z-reaction_rate hist
  Int_t N_rate1=4;//(number of scatterd particle)/molecular
  Int_t N_rate2=10;//(number of recoli partcle)/molecular  
  Int_t theta_min=3; //deg; rutherford scattar angle
  Double_t dtheta=1;//dig; bin of angle

  Int_t r_attpc=275; //mm
  Int_t l_attpc=1000;
  
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
  /*  Double_t I=-2/(sin((180-theta_min)*deg_to_rad/2)
		 *sin((180-theta_min)*deg_to_rad/2))
    +2/(sin(theta_min*deg_to_rad/2)*sin(theta_min*deg_to_rad/2));
  //  cout<<I<<endl;
  */
    
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  Double_t b=alfa*hc*z_s*z_r2/4;  
  
  //set variable
  Int_t div_l=(int)(l_attpc/dz+0.9);
  //  Int_t div_ll;
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/dtheta+0.9);

  Double_t rmax1;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t theta1_cms_r;
  
  Double_t rmax2;
  Double_t E2;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];

  Int_t n_p1;
  Int_t n_p2;
  Double_t Y1_rmax=0;

  Int_t n_rmax2[div_l]; 
  Int_t n_rmax[div_l];
  Double_t l1[div_l];
  Double_t l2[div_l];

  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  Double_t n_react1[div_l];
  //  Double_t n_beam1[div_l];
  Int_t n_rmax1=0;
  Double_t Y1[div_l];
  
  Double_t n2[div_l];
  Double_t Esum2[div_l];
  Double_t Eave2[div_l];
  Double_t Y2[div_l];

  //  Double_t sigma[div_l];
  Double_t Y1_sum;
  Double_t z[div_r];

  bool fin1=false;
  Double_t z1_max;
  
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
  Double_t I[Entry1];  
  vector<vector<Double_t>> n_beam1(div_l, vector<Double_t>(Entry1,0));
  vector<vector<Double_t>> sigma1(div_l, vector<Double_t>(Entry1,0));
  vector<vector<Double_t>> Y1_da(div_l, vector<Double_t>(Entry1,0));

  
  //set histgram
  TH1F *h_rmax=new TH1F("h_rmax","rmax:reaction_rate",nbin,0,300);
  TH1F *h_z=new TH1F("h_z","depth:reaction_rate",nbin_zn,0,400);
  
  //12c12c calculate
  for(int i=0;i < div_l ; i++){
    //  for(int i=0;i <100; i++){
    l1[i]=i*dz;
    l2[i]=(i+1)*dz;
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;

    //calculate Eave
    for(int j=0; j< Entry1; j++){
      cout<<"at "<<l1[i]<<" mm Reading events "<<100*j/Entry1 <<" % !               \r"<<flush;
      tree1->GetEntry(j);
      tree1->SetBranchStatus("*",0);
      tree1->SetBranchStatus("startz",1);
      tree1->SetBranchStatus("TotalE",1);
      z1[0]=startz1[0]*ctom;
      z1[1]=startz1[1]*ctom;

      if(z1_max<z1[0]){z1_max=z1[0];}
      else if(z1_max<z1[1]){z1_max=z1[1];}
      
      if(j%2==1 && fabs(z1[0]-z1[1])<dz && z1[0]>=l1[i] && z1[0]<l2[i]){
	n1[i]++;
	Esum1[i]+=E1;    
	//	    if(i>=350){cout<<i<<"   "<<E1<<"               "<<endl}
      }
      else if(j%2==0){//cout<<i<<" "<<j<<" This is beam event"<<endl;
	continue;
      }
      else {continue;}
    }
    //    cout<<"z1_max "<<z1_max<<endl;
    
    if(n1[i]<1){//Eave1[i]=0;
      continue;
    }
    else{
      Eave1[i]=Esum1[i]/n1[i];
      //      cout<<i<<"     "<<Eave1[i]<<"            "<<endl;
    }
    
    //calculate reaction rate
    // need to loop about angle like z and div_l (l1, l2 etc)
    //  each 1 deg in the loop
    Y1[i]=0;
    n_react1[i]=0;

    //need to modify
    /*
    for(int j=0; j< Entry1; j++){
      cout<<"at "<<l1[i]<<" mm calculate reaction rate "<<100*j/Entry1 <<" % !\r"<<flush;
      tree1->GetEntry(j);
      tree1->SetBranchStatus("*",0);
      tree1->SetBranchStatus("startz",1);
      tree1->SetBranchStatus("theta_cms",1);
      z1[0]=startz1[0]*ctom;
      z1[1]=startz1[1]*ctom;
      
      I[j]=0;
      
      if(j%2==1 && fabs(z1[0]-z1[1])<dz && z1[0]>=l1[i] && z1[0]<l2[i]){
	if(theta1_cms[0]<theta1_cms[1]){
	  theta1_cms_r=theta1_cms[1]*deg_to_rad;
	}
	else if(theta1_cms[1]<theta1_cms[0]){
	  theta1_cms_r=theta1_cms[0]*deg_to_rad;
	}
      }
      else if(j%2==0){//cout<<i<<" "<<j<<" This is beam event"<<endl;
	continue;
      }
      else {//cout<<i<<" "<<j<<endl;
	continue;
      }
      
      I[j]=2*pi*dtheta_r*sin(theta1_cms_r)
	/(sin(theta1_cms_r/2)*sin(theta1_cms_r/2)
	  *sin(theta1_cms_r/2)*sin(theta1_cms_r/2));
      //      printf(" %f \n", I[j]);      
      
      //      if(i==0 && j==0){n_beam1[i]=n_b;}
      if(i==0 && j==0){n_beam1.at(i).at(j)=n_b;}
      //      else{n_beam1[i]=n_b-Y1_sum;}
      else{n_beam1.at(i).at(j)=n_b-Y1_sum;}
      
      //      if(n_beam1[i]<0){
      if(n_beam1.at(i).at(j)<0){
	//	div_ll=l1[i];
	cout<<"end the loop at"<<l1[i]<<" mm~"<<l2[i]<<" mm, angle  "<<theta1_cms_r*rad_to_deg<<" deg "<<endl;

	fin1=true;
	break;
      }
      if(n1[i]<1){
	continue;
      }
      else{      
	//	sigma[i]=2*pi*a*a*I/(Eave1[i]*Eave1[i]);
	sigma1.at(i).at(j)=a*a*I[j]/(Eave1[i]*Eave1[i]);
	//	Y1[i]=sigma[i]*n_beam1[i]*rho*dz*1e-3*N_rate1;
	Y1_da.at(i).at(j)=sigma1.at(i).at(j)*n_beam1.at(i).at(j)*rho*dz*1e-3*N_rate1;
	//	Y1_sum+=Y1[i];
	Y1[i]+=Y1_da.at(i).at(j);
	Y1_sum+=Y1_da.at(i).at(j);
	
	//	cout<<i<<" "<<j<<"  Y1_da "<<Y1_da.at(i).at(j)<<"        Y1[i] "<<Y1[i]<<"     Y1_sum "<<Y1_sum<<endl;
      }
      //      if(i>=350){cout<<i<<"     "<<n_beam1[i]<<"         "<<sigma[i]<<"      "<<Y1[i]<<"         "<<Y1_sum<<"     "<<Eave1[i]<<endl}
      //      if(i>330){cout<<i<<" "<<j<<"         Y1_sum  "<<Y1_sum<<"                      "<<endl;}
      //      if(i==200){cout<<i<<" "<<j<<"         Y1_sum  "<<Y1_sum<<"                      "<<endl;}
    }
    if(fin1){break;}

    */
    
    //z-reaction rate hist filling data
    if(n1[i]>0){
      h_z->Fill(l1[i],Y1[i]);
    }
    if (i%100==0){
      cout<<"12C12C finish "<<l1[i]<<" mm !!                                "<<endl;
    }

    //need to modify
    /*
    if(z1_max<l2[i]){
      //      cout<<"end the loop at"<<l1[i]<<" mm~"<<l2[i]<<" mm, angle  "<<theta1_cms_r*rad_to_deg<<" deg "<<endl;
      cout<<"                                                  "<<endl;
      cout<<"break the loop at "<<l1[i]<<" mm~"<<l2[i]<<" mm, because of no reaction in over "<<l1[i]<<" mm"<<endl;
      //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
      break;
    }
    */
  }
  
  //make conections between rmax and reaction rate
  for(int j=0; j< Entry1; j++){
    tree1->GetEntry(j);
    tree1->SetBranchStatus("*",0);
    tree1->SetBranchStatus("startz",1);
    tree1->SetBranchStatus("rmax",1);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    
    for(int i=0;i < div_l ; i++){
      l1[i]=i*dz;
      l2[i]=(i+1)*dz;
      
      if(j%2==1 && fabs(z1[0]-z1[1])<dz && z1[0]>=l1[i] && z1[0]<l2[i]){
	if(n1[i]>0){
	  h_rmax->Fill(rmax1,Y1[i]/n1[i]);
	  n_rmax1++;
	  Y1_rmax+=Y1[i]/n1[i];
	}
      }
      /*
      else{//	cout<<j<<"   "<<l1[i]<<"   "<<z1[0]<<"       "<<z2[1]<<endl;
	continue;
      }
      */
    }
    cout<<"Filling data "<<100*j/Entry1 <<" % !  \r"<<flush;
  }
  cout<<"n_rmax1"<<n_rmax1<<"    Y1sum"<<Y1_rmax<<"                            "<<endl;
  
  //search rmax when over 100Hz
  double hist_cnt = 0;
  int flag_100hz=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt += h_rmax->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt > 100 && flag_100hz==0 ){
      printf("over 100Hz %d mm \n\n", ihist);
      flag_100hz=1;
    }
  }

  //make error bar
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax->GetBinContent(m));
    h_rmax->SetBinError(m,err);
  }
  for(int l=0; l<=nbin_zn;l++){
    double err_z=sqrt(h_z->GetBinContent(l));
    h_z->SetBinError(l,err_z);
  }

  //make and save the hist 
  // z-reaction rate
  h_z->SetStats(0);
  //  h_z->Draw("hist");
  //  h_z->SaveAs("data/ruther_macro_result_c12_each_angle_depth.C");

  // rmax-reaction rate
  h_rmax->SetStats(0);
  h_rmax->Draw("hist");
  //  h_rmax->SaveAs("data/ruther_macro_result_c12_each_angle_rmax.C");

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

  return 0;
}
