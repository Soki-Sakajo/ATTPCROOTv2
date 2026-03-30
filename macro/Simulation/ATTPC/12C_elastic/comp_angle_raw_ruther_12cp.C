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

int comp_angle_raw_ruther_12cp(Int_t subnum=30){

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

  Int_t l_ext=10;   // mm; extract information
  Double_t ra_fit=-100;
  Double_t y_fit=-100;
  
  //set constants
  Int_t ctom=10;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; // Pa
  Double_t alfa=7.297352569e-3; // fine structure constant
  Double_t hc=197.463587*1e-15; // MeV•m
  Double_t k_b=1.38065e-23; // J/K ; Bolzman constant
  Double_t rho=P*atmtoPa/(k_b*T);
  //  Double_t m12=11177.73600; // MeV/c^2 ; 12C mass
  
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  Double_t b=alfa*hc*z_s*z_r2/4;  
  
  //set variable
  // common
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Double_t da_r=da*deg_to_rad;
  
  //for 12Cp
  bool setbeam2=true;
  bool fin_beam2=false;
  bool fin_rea2=false;
  bool checkang=true;
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

  Double_t rmax2;
  Double_t rmax2_nd;
  Double_t E2;
  Double_t z2_max=0;
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
  TString dir="data2";
  
  // Tchain of 12Cp
  TChain *tree2=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    TString refi="sim_macro_p_130matm_" + to_string(i) + ".root";
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
  //  for 12C12p
  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve2);
  tree2->SetBranchAddress("startz",&startz2);
  tree2->SetBranchAddress("TotalE",&E2);
  tree2->SetBranchAddress("theta_cms",&theta2_cms);
  tree2->SetBranchAddress("n",&n_p2);

  //set histgram
  //  TH1F *h_ang=new TH1F("h_ang","ang:number_of_event",200,0,200);
  TH1F *h_ang=new TH1F("h_ang","ang:number_of_event",180,0,180);
  //  TH1F *h_sin=new TH1F("h_sin","ang:number_of_event",200,0,200);
  TH1F *h_sin=new TH1F("h_sin","ang:number_of_event",180,0,180);
  TF1 *f= new TF1("f","[0]*sin([1]*x)",0,180);
  
  TH1F *h_thetaY1=new TH1F("h_thetaY1","angle :reaction_rate",180,0,180);
  TH1F *h_thetaY2=new TH1F("h_thetaY2","angle :reaction_rate",180,0,180);
  TH1F *h_thetaY3=new TH1F("h_thetaY3","angle :reaction_rate",180,0,180);

  TF1 *f1=new TF1("f1","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2))+1/(cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)*cos([1]*x/2)))",0,180);
  TF1 *f2=new TF1("f2","[0]*sin([1]*x)*(1/(sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)*sin([1]*x/2)))",0,180);
  
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
      if(z2[0]>=l_ext && l_ext+dz>z2[0]){
	h_ang->Fill(theta2);
      }
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

      if(l2[i]==l_ext){
	//	h_thetaY1->Fill(ang2[j],I2.at(i).at(j)*1e30);
	h_thetaY2->Fill(ang2[j],I2.at(i).at(j)*1e30);
	//	h_thetaY2->Fill(180-ang2[j],I2.at(i).at(j)*1e30);
      }
    }
    if (int(l2[i])%100==0){
      //need
      cout<<"12Cp finish "<<l2[i]<<" mm !!                              "<<endl;
    }
  }
  
  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists",1500,1000);
  can->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.15);
  h_ang->SetStats(0);
  h_ang->SetTitle("Scattering angle distribution (isotropic) (12cp)");
  h_ang->GetXaxis()->SetTitle("angle in cm [deg]");
  h_ang->GetYaxis()->SetTitle("number of event  N(#theta) #propto #frac{d#sigma}{d#Omega}(#theta) #Delta#Omega");
  h_ang->Draw("hist");
  h_ang->SetTitleSize(0.08);
  h_ang->GetXaxis()->SetLabelSize(0.06);
  h_ang->GetYaxis()->SetLabelSize(0.06);
  h_ang->GetXaxis()->SetTitleSize(0.06);
  h_ang->GetXaxis()->SetTitleOffset(0.8);
  h_ang->GetYaxis()->SetTitleSize(0.06);
  h_ang->GetYaxis()->SetTitleOffset(1.1);
  //  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex1;
  latex1.SetTextSize(0.06);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(150,1e2,"#propto sin#theta");
  can->Update();

  TCanvas *can2=new TCanvas("hist2","hist2",1500,1000);
  can2->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.18);
  h_thetaY2->SetStats(0);
  h_thetaY2->SetTitle("Rutherford scattering angle distribution (12cp)");
  h_thetaY2->GetXaxis()->SetTitle("angle in cm [deg]");
  h_thetaY2->GetYaxis()->SetTitle("reaction rate [cps]");
  h_thetaY2->Draw("hist");
  h_ang->SetTitleSize(0.08);
  h_thetaY2->GetXaxis()->SetLabelSize(0.08);
  h_thetaY2->GetYaxis()->SetLabelSize(0.08);
  h_thetaY2->GetXaxis()->SetTitleSize(0.08);
  h_thetaY2->GetYaxis()->SetTitleSize(0.08);
  //  h_thetaY2->GetYaxis()->SetTitleOffset(1);
  h_thetaY2->GetYaxis()->SetTitleOffset(0.95);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.09);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(40,14,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  TLatex latex3;
  latex3.SetTextSize(0.09);
  latex3.SetTextColor(kBlack);
  latex3.DrawLatex(65,2e-1,"#propto #frac{sin#theta}{sin^{4}#frac{#theta}{2}}");
  can->Update();

  /*
  TString outcan1="test.pdf";
  TString outcan2="test_rut.pdf";

  TString outcan1="test_rut.pdf";
  TString outcan2="test.pdf";
  */
  TString outcan1="./pdf_output/can_12cp_angle_raw.pdf";
  TString outcan2="./pdf_output/can_12cp_angle_rut.pdf";

  can->SaveAs(outcan1);
  can2->SaveAs(outcan2);

  /*
  TCanvas *can=new TCanvas("hists","hists");
  //  can->Divide(3,2);
  //  can->Divide(3,1);
  //  can->Divide(2,1);
  can->cd(2);
  gPad->SetLogy(1);
  h_thetaY1->SetStats(0);
  h_thetaY1->Draw("hist");
  h_thetaY1->GetXaxis()->SetLabelSize(0.06);
  h_thetaY1->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();
  //  can->cd(3);
  //  can->cd(2);
  can->cd(4);
  gPad->SetLogy(1);
  h_ang->SetStats(0);
  f->SetParameters(650,deg_to_rad);
  f->SetParLimits(1,deg_to_rad,deg_to_rad);
  h_ang->Fit("f","","",theta_min,theta_max);
  h_ang->Draw("hist");
  f->Draw("same");
  double pa=f->GetParameter(0);
  double ea=f->GetParError(0);
  double pb=f->GetParameter(1);
  double eb=f->GetParError(1);
  cout<<endl;
  cout<<"set depth "<<l_ext<<" mm"<<endl;
  cout<<"a  "<<pa<<"  ea  "<<ea<<"  b  "<<pb<<"  eb  "<<eb<<endl;
  can->cd(5);
  gPad->SetLogy(1);
  h_thetaY1->SetStats(0);
  f2->SetParameters(1,deg_to_rad);
  f2->SetParLimits(0,0,2);
  f2->SetParLimits(1,deg_to_rad,deg_to_rad);
  h_thetaY1->Fit("f2","","",theta_min,theta_max);
  h_thetaY1->Draw("hist");
  f2->Draw("same");
  Double_t chi2=f2->GetChisquare();
  Double_t ndf2=f2->GetNDF();
  cout<<chi2/ndf2<<endl;
  can->cd(6);
  gPad->SetLogy(1);
  h_thetaY2->SetStats(0);
  f1->SetParameters(1,deg_to_rad);
  f1->SetParLimits(0,0,2);
  f1->SetParLimits(1,deg_to_rad,deg_to_rad);
  h_thetaY2->Fit("f1","","",theta_min,theta_max);
  h_thetaY2->Draw("hist");
  f1->Draw("same");
  Double_t chi3=f1->GetChisquare();
  Double_t ndf3=f1->GetNDF();
  cout<<chi3/ndf3<<endl;
  */

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
