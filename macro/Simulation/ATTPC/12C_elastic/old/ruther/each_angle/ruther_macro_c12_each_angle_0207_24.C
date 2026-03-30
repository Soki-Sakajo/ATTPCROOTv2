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

int ruther_macro_c12_each_angle_0207_24(){
  //int ruther_macro_c12_each_angle_0207_24(Int_t subnum=0){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t z_s=6;          // Z of scatterd particl in simfile
  Int_t z_r1=6;         // Z of an recoil particl in simfile
  Int_t z_r2=1;         // Z of another recoil particl in simfile
  Int_t dz=1;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  Double_t P=0.13;      // atm; gas pressure
  Double_t da=1;        // deg; bin of angle
  Double_t dr=1;        // mm; bin of reaction rate
  Int_t nbin=300/dr;    // bin of rmax-reaction_rate hist
  Int_t nbin_zn=400/dz; // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t N_rate2=10;     // (number of recoli partcle)/molecular  
  Int_t theta_min=3;    // deg; rutherford scattar angle
  //  Double_t dtheta=1;    // dig; bin of angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  int ctom=10;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; // MPa
  Double_t alfa=7.297352569e-3; // fine structure constant
  Double_t hc=197.463587*1e-15; // MeV•m
  Double_t k_b=1.38065e-23; // J/K ; Bolzman constant
  //  Double_t dtheta_r=da*deg_to_rad;
  Double_t da_r=da*deg_to_rad;
  Double_t rho=P*atmtoPa/(k_b*T);
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  
  //set variable
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Int_t eve;
  Int_t n_data;
  Int_t n_p1;
  Double_t rmax1;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t theta1_cms_r;
  Double_t Y1_rmax=0;

  //  Int_t n_rmax[div_l];
  Double_t l11[div_l];
  Double_t l12[div_l];
  Double_t ang11[div_a];
  Double_t ang_r11[div_a];
  Double_t ang12[div_a];
  Double_t ang_r12[div_a];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];

  //  Double_t n_a1[div_a];
  vector<vector<Double_t>> n_a1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> I(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> n_beam1(div_l, vector<Double_t>(div_a,0));
  //  vector<vector<vector<Double_t>>> n_beam1(div_l, vector<vector<Double_t>>(div_a, vector<Double_t>(Entry1,0)));
  vector<vector<Double_t>> sigma1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_da(div_l, vector<Double_t>(div_a,0));
  //  Double_t Y1_nom[div_l];
  vector<vector<Double_t>> Y1_nom(div_l, vector<Double_t>(div_a,0));

  Double_t n_react1[div_l];
  //  Double_t n_beam1[div_l];
  Int_t n_rmax1=0;
  Double_t Y1[div_l];

  Double_t Y1_sum=0;
  Double_t z[div_r];
  Double_t z1_max;
  Double_t dep1;
  Double_t dir1;
  Double_t Y1_eve;

  bool setbeam=true;
  bool fin_beam1=false;
  //  bool fin_depth1=false;

  //for 12Cp
  Double_t b=alfa*hc*z_s*z_r2/4;  
  Double_t rmax2;
  Double_t E2;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];
  Int_t n_p2;
  Int_t n_rmax2[div_l];   
  Double_t n2[div_l];
  Double_t n_a2[div_a];
  Double_t Esum2[div_l];
  Double_t Eave2[div_l];
  Double_t Y2[div_l];
  
  //input
  TString infname = "data/sim_macro_c12.root";
  //  TString infname = "data/sim_macro_c12_" + to_string(subnum) + ".root";
  //  TString infname = "data1/sim_macro_c12_" + to_string(subnum) + ".root";
  TFile *file1 =new TFile(infname); 
  TTree *tree1 = (TTree*)file1->Get("anadata");

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);

  Double_t theta1[Entry1];
 
  //set histgram
  TH1F *h_rmax=new TH1F("h_rmax","rmax:reaction_rate",nbin,0,300);
  TH1F *h_z=new TH1F("h_z","reaction position:reaction_rate",nbin_zn,0,400);
  TH2F *h_ztheta=new TH2F("h_ztheta","reaction position:angle",400,0,400,100,0,100);

  //set output for making ratdata by each subnum
  TString outdata = "data/rutdata_c12.root";
  //  TString outdata = "data/rutdata_c12_" + to_string(subnum) + ".root";
  //  TString outdata = "data1/rutdata_c12_" + to_string(subnum) + ".root";
  TFile *rutfile=new TFile(outdata, "recreate");
  TTree *rutdata = new TTree("rutdata","rutdata");
  
  rutdata->Branch("event_number",&eve,"eve/I");
  rutdata->Branch("reaction_position",&dep1,"dep1/D");
  rutdata->Branch("angle",&dir1,"dir1/D");
  rutdata->Branch("rmax_c12",&rmax1,"rmax1/D");
  rutdata->Branch("reaction_rate",&Y1_eve,"Y1_eve/D");
  rutdata->Branch("depth",&l11,"l11/D");
  rutdata->Branch("Eave",&Eave1,"Eave1/D");
  
  //make event lists
  /*
  for(int i=0;i < div_l ; i++){
    l11[i]=i*dz;
    l12[i]=(i+1)*dz;
    for(int k=0; k< div_a; k++){
      ang11[k]=k*da+theta_min;
      ang12[k]=(k+1)*da+theta_min;
      ang_r11[k]=ang11[k]*deg_to_rad;
      ang_r12[k]=ang12[k]*deg_to_rad;


      tree->Draw();


   */
  
  //12c12c calculate
  // devide depth
  for(int i=0;i < div_l ; i++){
    l11[i]=i*dz;
    l12[i]=(i+1)*dz;
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;

    //calculate Eave
    for(int j=0; j< Entry1; j++){
      tree1->GetEntry(j);
      z1[0]=startz1[0]*ctom;
      z1[1]=startz1[1]*ctom;
      if(z1_max<z1[0]){z1_max=z1[0];}
      if(z1_max<z1[1]){z1_max=z1[1];}
      if(n_p1==2 && fabs(z1[0]-z1[1])<dz && z1[0]>=l11[i] && z1[0]<l12[i]){
	n1[i]++;
	Esum1[i]+=E1;    
      }
      else if(n_p1==1){continue;}
      else {continue;}
    }
    if(n1[i]<1){continue;}
    else{Eave1[i]=Esum1[i]/n1[i];}

    //calculate reaction rate
    Y1[i]=0;
    n_react1[i]=0;

    //divide angle
    for(int k=0; k< div_a; k++){
      ang11[k]=k*da+theta_min;
      ang12[k]=(k+1)*da+theta_min;
      ang_r11[k]=ang11[k]*deg_to_rad;
      ang_r12[k]=ang12[k]*deg_to_rad;
      n_a1.at(i).at(k)=0;
      
      for(int j=0; j< Entry1; j++){
	//need
	cout<<"at "<<l11[i]<<" mm "<<ang11[k]<<" deg, Reading event "<<100*j/Entry1 <<" % !               \r"<<flush;

	tree1->GetEntry(j);
	z1[0]=startz1[0]*ctom;
	z1[1]=startz1[1]*ctom;

	if(theta1_cms[0]<=theta1_cms[1]){theta1[j]=theta1_cms[0];}
	else if(theta1_cms[0]>theta1_cms[1]){theta1[j]=theta1_cms[1];}
	if(n_p1==2 && fabs(z1[0]-z1[1])<dz && z1[0]>=l11[i] && z1[0]<l12[i]){
	  if(theta1[j]>=ang11[k] && theta1[j]<ang12[k]){
	    n_a1.at(i).at(k)++;
	    //	    cout<<j<<"  "<<i<<"   "<<k<<"  "<<n_a1.at(i).at(k)<<endl;

	  }
	}
	else {continue;}
      }
      if(n_a1.at(i).at(k)<1){continue;}

      I.at(i).at(k)=2*pi*da_r*sin(ang_r11[k])
	  /(sin(ang_r11[k]/2)*sin(ang_r11[k]/2)
	    *sin(ang_r11[k]/2)*sin(ang_r11[k]/2))
	+2*pi*da_r*sin(ang_r11[k])
	/(sin((pi-ang_r11[k])/2)*sin((pi-ang_r11[k])/2)
	  *sin((pi-ang_r11[k])/2)*sin((pi-ang_r11[k])/2));      

      if(setbeam){
	n_beam1.at(i).at(k)=n_b;
	setbeam=false;
      }
      else{n_beam1.at(i).at(k)=n_b-Y1_sum;}
      //      cout<<i<<"       "<<k<<"      "<<n_beam1.at(i).at(k)<<endl;
      if(n_beam1.at(i).at(k)<0){cout<<"end the loop at"<<l11[i]<<" mm~"<<l12[i]<<" mm, "<<ang11[k]<<" deg~"<<ang12[k]<<" deg because beam particle is 0."<<endl;
	fin_beam1=true;
	break;
      }
      else{

	sigma1.at(i).at(k)=2*pi*a*a*I.at(i).at(k)/(Eave1[i]*Eave1[i]);
	Y1_da.at(i).at(k)=sigma1.at(i).at(k)*n_beam1.at(i).at(k)*rho*dz*1e-3*N_rate1;
	Y1_nom.at(i).at(k)+=Y1_da.at(i).at(k)/n_a1.at(i).at(k);
	Y1[i]+=Y1_da.at(i).at(k);
	Y1_sum+=Y1_da.at(i).at(k);
	//	cout<<i<<"   "<<k<<"          "<<Y1_da.at(i).at(k)<<"          "<<Y1[i]<<"          "<<Y1_sum<<endl;
	
      }
    }
    
    //z-reaction rate hist filling data
    if(n1[i]>0){
      //      h_z->Fill(l11[i],Y1_nom[i]);
      h_z->Fill(l11[i],Y1[i]);
    }
    if (i%100==0){
      //need
      cout<<"12C12C finish "<<l11[i]<<" mm !!                                "<<endl;
    }
    if(fin_beam1){break;}    
    if(z1_max<l12[i]){
      cout<<"                                                                 "<<endl;
      cout<<"end the loop at "<<l11[i]<<" mm~"<<l12[i]<<" mm, because of no reaction in over "<<l12[i]<<" mm"<<endl;
      //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
      break;
    }
  }

  //make conections between rmax and reaction rate
  for(int j=0; j< Entry1; j++){
    tree1->GetEntry(j);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    
    for(int i=0;i < div_l ; i++){
      l11[i]=i*dz;
      l12[i]=(i+1)*dz;
      if(n_p1==2 && fabs(z1[0]-z1[1])<dz && z1[0]>=l11[i] && z1[0]<l12[i]){
	
	for(int k=0; k< div_a; k++){
	  ang11[k]=k*da+theta_min;
	  ang12[k]=(k+1)*da+theta_min;
	  ang_r11[k]=ang11[k]*deg_to_rad;
	  ang_r12[k]=ang12[k]*deg_to_rad;	  
	  if(theta1_cms[0]<=theta1_cms[1]){theta1[j]=theta1_cms[0];}
	  else if(theta1_cms[0]>theta1_cms[1]){theta1[j]=theta1_cms[1];}

	  if(theta1[j]>=ang11[k] && theta1[j]<ang12[k] && n_a1.at(i).at(k)>0 ){
	    dep1=l11[i];
	    dir1=ang11[k];

	    Y1_eve=Y1_nom.at(i).at(k);
	    //	    h_rmax->Fill(rmax1,Y1_nom.at(i).at(k));
	    h_rmax->Fill(rmax1,Y1_eve);
	    n_rmax1++;
	    Y1_rmax+=Y1_nom.at(i).at(k);	    
	    h_ztheta->Fill(l11[i],ang11[k]);
	    //	    cout<<i<<"  "<<k<<"  "<<j<<"  "<<l11[i]<<"  "<<dep1<<"  "<<ang11[k]<<"  "<<dir1<<endl;
	    rutdata->Fill();

	  }
	  else{continue;}
	}
      }      
      else{continue;}
    }
    //need
    cout<<"Filling data "<<100*j/Entry1 <<" % !  \r"<<flush;
  }

  cout<<"n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax<<"                          "<<endl;
  
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
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z->GetBinContent(m));
    h_z->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta->GetBinContent(m));
    h_ztheta->SetBinError(m,err);
  }

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(2,2);
  
  // z-angle
  h_ztheta->SetStats(0);
  can->cd(1);
  h_ztheta->Draw("colz");

  // z-reaction rate
  h_z->SetStats(0);
  can->cd(2);
  gPad->SetLogy(1);
  h_z->Draw("hist");

  // rmax-reaction rate
  h_rmax->SetStats(0);
  can->cd(3);
  gPad->SetLogy(1);
  h_rmax->Draw("hist");

  //output for canvas
  TString outcan = "data/ruther_12c_can.C";
  //  TString outcan = "data/ruther_12c_can_" + to_string(subnum) + ".C";
  //  TString outcan = "data1/ruther_12c_can_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  h_ztheta->Write();
  h_z->Write();
  h_rmax->Write();
  rutdata->AutoSave();
  rutfile->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
