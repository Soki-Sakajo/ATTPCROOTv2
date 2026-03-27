//#define CAN_FLL

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

int comp_depth_raw_ruther_12c12c(Int_t subnum=30){

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
  
  //input multi-file by using TChain
  //common set for input
  //  TString dir="data1";
  TString dir="data2";
  
  // Tchain of 12c12c
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<subnum; i++){
    //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
    //      TString refi="sim_macro_c12_" + to_string(i) + ".root";
    TString refi="sim_macro_12c_130matm_" + to_string(i) + ".root";
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

  //set histgram
  TH1F *h_gene_12c=new TH1F("h_gene_12c","reaction position:generate (12c12c)",nbin_zn,0,400);
  TH1F *h_z_12c=new TH1F("h_z_12c","reaction position:reaction_rate (12c12c)",nbin_zn,0,400);

  TH1F *h_geen_12c=new TH1F("h_geen_12c","beam energy:generate (12c12c)",500,0,50);
  TH1F *h_en_12c=new TH1F("h_en_12c","beam energy:reaction_rate (12c12c)",50,0,50);

  TH2F *h_12c=new TH2F("h_12c","reaction position:beam energy (12c12c)",nbin_zn,0,400,104,0,52);
  TH1F *h_zen_12c=new TH1F("h_zen_12c","reaction position:ave beam energy (12c12c)",nbin_zn,0,400);

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
    /*
    if(theta1_cms[0]<=theta1_cms[1]){theta1=theta1_cms[0];}
    else if(theta1_cms[0]>theta1_cms[1]){theta1=theta1_cms[1];}
    */
    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0){
    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0 && theta1+theta_min<181){
    //    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>=0 && theta1<=theta_max){
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
      h_gene_12c->Fill(z1[0]);
      h_geen_12c->Fill(E1);
      h_12c->Fill(z1[0],E1);
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
  
  //  for(int i=0;i < div_l ; i++){
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
    h_zen_12c->Fill(l1[i],Eave1[i]*2);
    
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
	cout<<"end the loop at"<<l1[i]<<" mm~"<<l1[i]+dz<<" mm, "
	    <<ang1[j]<<"deg~"<<ang1[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
	cout<<"max point at "<<l1[i]<<" mm~"<<l1[i]+dz<<" mm "
	    <<ang1[j]<<" deg~"<<ang1[j]+da<<" deg, because of no reaction in over "
	    <<l1[i]+dz<<" mm "<<ang1[j]+da<<" deg"<<endl<<endl;

	fin_beam1=true;
	break;
      }
      else{
	Y1_da.at(i).at(j)=I1.at(i).at(j)*n_beam1.at(i).at(j)*rho*dz*1e-3*N_rate1;
	Y1[i]+=Y1_da.at(i).at(j);
      }
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
      h_z_12c->Fill(l1[i],Y1[i]);
      h_en_12c->Fill(Eave1[i]*2,Y1[i]);
    }
    if (int(l1[i])%100==0){
      //need
      cout<<"12C12C finish "<<l1[i]<<" mm !!                              "<<endl;
    }
    if(fin_beam1){break;}
    if(fin_rea1){break;}
    //    if(z1_max<l1[i]+dz){
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
	n_rmax1++;
	n_rmax1_nd++;
	Y1_rmax+=Y1_eve;
      }
    }
    //need
    cout<<"Filling 12c data "<<100*i/ev_ve1.size() <<" % !    \r"<<flush;
  }
  cout<<"                              "<<endl;
  cout<<"12c12c n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax<<"  n_rmax1_nd  "<<n_rmax1_nd<<"                   "<<endl;
  cout<<endl;

  #ifdef CAN_FLL
  //make and save the hist 
  TCanvas *can=new TCanvas("can","can");
  can->Divide(2,3);
  
  can->cd(1);
  gPad->SetLogy(1);
  h_gene_12c->SetStats(0);
  h_gene_12c->Draw("hist");
  h_gene_12c->GetXaxis()->SetLabelSize(0.06);
  h_gene_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(2);
  gPad->SetLogy(1);
  h_z_12c->SetStats(0);
  h_z_12c->Draw("hist");
  h_z_12c->GetXaxis()->SetLabelSize(0.06);
  h_z_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(3);
  gPad->SetLogy(1);
  h_geen_12c->SetStats(0);
  h_geen_12c->Draw("hist");
  h_geen_12c->GetXaxis()->SetLabelSize(0.06);
  h_geen_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(4);
  gPad->SetLogy(1);
  h_en_12c->SetStats(0);
  h_en_12c->Draw("hist");
  h_en_12c->GetXaxis()->SetLabelSize(0.06);
  h_en_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(5);
  h_12c->SetStats(0);
  h_12c->Draw("colz");
  h_12c->GetXaxis()->SetLabelSize(0.06);
  h_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(6);
  h_12c->SetStats(0);
  h_zen_12c->Draw("hist");
  h_zen_12c->GetXaxis()->SetLabelSize(0.06);
  h_zen_12c->GetYaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();
  #endif

  TCanvas *can1=new TCanvas("hists1","hists1",1500,1000);
  can1->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.16);
  h_geen_12c->SetStats(0);
  //  h_geen_12c->SetMinimum(1e2);
  h_geen_12c->SetMinimum(1e-1);
  h_geen_12c->SetMaximum(2e4);
  //  h_geen_12c->SetTitle("Energy distribution (uniform energy loss) (12c12c)");
  //  h_geen_12c->SetTitle("Energy distribution (uniform energy loss)");
  h_geen_12c->GetXaxis()->SetTitle("Energy [MeV]");
  h_geen_12c->GetYaxis()->SetTitle("Generated rate");
  h_geen_12c->Draw("hist");
  h_geen_12c->SetTitleSize(0.08);
  h_geen_12c->GetXaxis()->SetLabelSize(0.08);
  h_geen_12c->GetYaxis()->SetLabelSize(0.08);
  h_geen_12c->GetXaxis()->SetTitleSize(0.08);
  //  h_geen_12c->GetXaxis()->SetTitleOffset(0.8);
  h_geen_12c->GetYaxis()->SetTitleSize(0.08);
  h_geen_12c->GetXaxis()->SetTitleOffset(0.9);
  h_geen_12c->GetYaxis()->SetTitleOffset(0.9);
  gStyle->SetOptTitle(0);
  gPad->Update();

  TCanvas *can2=new TCanvas("hists2","hists2",1500,1000);
  can2->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.16);
  h_gene_12c->SetStats(0);
  h_gene_12c->SetMinimum(1e-1);
  h_gene_12c->SetMaximum(2e4);
  h_gene_12c->SetTitle("Energy distribution (uniform energy loss) (12c12c)");
  //  h_gene_12c->SetTitle("Energy distribution (uniform energy loss)");
  h_gene_12c->GetXaxis()->SetTitle("Depth in AT-TPC [mm]");
  h_gene_12c->GetYaxis()->SetTitle("Generated rate");
  h_gene_12c->Draw("hist");
  h_gene_12c->SetTitleSize(0.08);
  h_gene_12c->GetXaxis()->SetLabelSize(0.08);
  h_gene_12c->GetYaxis()->SetLabelSize(0.08);
  h_gene_12c->GetXaxis()->SetTitleSize(0.08);
  //  h_gene_12c->GetXaxis()->SetTitleOffset(0.8);
  h_gene_12c->GetYaxis()->SetTitleSize(0.08);
  h_gene_12c->GetXaxis()->SetTitleOffset(0.9);
  h_gene_12c->GetYaxis()->SetTitleOffset(0.9);
  gStyle->SetOptTitle(0);
  gPad->Update();
  /*
  TLatex latex1;
  latex1.SetTextSize(0.1);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(100,5e1,"N(z) #propto E_{beam}");
  can1->Update();
  */

  TCanvas *can3=new TCanvas("hist3","hist3",1500,1000);
  can3->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.16);
  h_z_12c->SetStats(0);
  h_z_12c->SetMinimum(1e-1);
  h_z_12c->SetMaximum(2e4);
  h_z_12c->SetTitle("Rutherford scattering energy distribution (12c12c)");
  h_z_12c->GetXaxis()->SetTitle("Depth in AT-TPC [mm]");
  h_z_12c->GetYaxis()->SetTitle("reaction rate [cps]");
  h_z_12c->Draw("hist");
  h_z_12c->SetTitleSize(0.08);
  h_z_12c->GetXaxis()->SetLabelSize(0.08);
  h_z_12c->GetYaxis()->SetLabelSize(0.08);
  h_z_12c->GetXaxis()->SetTitleSize(0.08);
  h_z_12c->GetYaxis()->SetTitleSize(0.08);
  h_z_12c->GetXaxis()->SetTitleOffset(0.9);
  h_z_12c->GetYaxis()->SetTitleOffset(0.9);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.1);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(20,2e3,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  TLatex latex3;
  latex3.SetTextSize(0.1);
  latex3.SetTextColor(kBlack);
  //  latex3.DrawLatex(80,2e3,"#propto E^{-2}_{z}");
  latex3.DrawLatex(50,1e2,"#propto E^{ -2}");
  can2->Update();
  /*
  TString outcan1="test.pdf";
  TString outcan2="test1.pdf";
  TString outcan3="test2.pdf";

  TString outcan1="test1.pdf";
  TString outcan2="test.pdf";
  TString outcan3="test2.pdf";

  TString outcan1="test1.pdf";
  TString outcan2="test2.pdf";
  TString outcan3="test.pdf";
  */
  TString outcan1="./pdf_output/can_12c12c_energy.pdf";
  TString outcan2="./pdf_output/can_12c12c_depth_raw.pdf";
  TString outcan3="./pdf_output/can_12c12c_depth_rut.pdf";

  //  TString outcan1="./pdf_output/can_depth.pdf";
  can1->SaveAs(outcan1);
  can2->SaveAs(outcan2);
  can3->SaveAs(outcan3);

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
