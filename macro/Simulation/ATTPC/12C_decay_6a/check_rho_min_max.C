#include <iostream>
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

int check_rho_min_max(Int_t subnum=0){

  Double_t rho_max;
  Double_t psi_max;
  Double_t theta_max;
  Double_t phi_max;
  Double_t x_max;
  Double_t y_max;
  Double_t z_max;
  Double_t rho_min;
  Double_t psi_min;
  Double_t theta_min;
  Double_t phi_min;
  Double_t x_min;
  Double_t y_min;
  Double_t z_min;

  //for hist
  const Int_t z_bin=10;
  const Int_t z_und=-10;
  const Int_t z_upp=1010;
  const Int_t xy_bin=10;
  const Int_t xy_und=-300;
  const Int_t xy_upp=300;
  const Int_t ang_bin=30;
  const Int_t psi_und=0;
  const Int_t psi_upp=360;
  const Int_t theta_und=0;
  //  const Int_t theta_und=60;
  const Int_t theta_upp=90;
  const Int_t phi_und=0;
  const Int_t phi_upp=360;
  const Int_t rho_bin=100;
  const Double_t rho_max_cen=1075.45;
  const Double_t rho_max_und=0.085;
  const Double_t rho_max_upp=0.1;
  const Double_t rho_min_cen=-1075.45;
  const Double_t rho_min_upp=-0.085;
  const Double_t rho_min_und=-0.1;
  
  //in and out dire set
  TString dir= "./data1/";

  //input
  TString infname= dir+ "cal_rho_"+ to_string(subnum) +".root";
  TFile *file = new TFile(infname);
  TTree *tree1 = (TTree*)file->Get("calmax");
  TTree *tree2 = (TTree*)file->Get("calmin");

  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("max_rho",&rho_max);
  tree1->SetBranchAddress("max_psi",&psi_max);
  tree1->SetBranchAddress("max_theta",&theta_max);
  tree1->SetBranchAddress("max_phi",&phi_max);
  tree1->SetBranchAddress("max_x",&x_max);
  tree1->SetBranchAddress("max_y",&y_max);
  tree1->SetBranchAddress("max_z",&z_max);

  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("min_rho",&rho_min);
  tree2->SetBranchAddress("min_psi",&psi_min);
  tree2->SetBranchAddress("min_theta",&theta_min);
  tree2->SetBranchAddress("min_phi",&phi_min);
  tree2->SetBranchAddress("min_x",&x_min);
  tree2->SetBranchAddress("min_y",&y_min);
  tree2->SetBranchAddress("min_z",&z_min);

  //set histgram
  // for data about max
  TH1F *h_x_max= new TH1F("h_x_max","x_max",xy_bin,xy_und,xy_upp);
  TH1F *h_y_max= new TH1F("h_y_max","y_max",xy_bin,xy_und,xy_upp);
  TH1F *h_z_max= new TH1F("h_z_max","z_max",z_bin,z_und,z_upp);
  TH1F *h_psi_max= new TH1F("h_psi_max","psi_max",ang_bin,psi_und,psi_upp);
  TH1F *h_theta_max= new TH1F("h_theta_max","theta_max",ang_bin,theta_und,theta_upp);
  TH1F *h_phi_max= new TH1F("h_phi_max","phi_max",ang_bin,phi_und,phi_upp);
  TH2F *h_xy_max= new TH2F("h_xy_max","x-y max",xy_bin,xy_und,xy_upp,xy_bin,xy_und,xy_upp);
  TH3D *h_ang_max= new TH3D("h_ang_max","psi-theta-phi max",
			    ang_bin,psi_und,psi_upp,ang_bin,theta_und,theta_upp,ang_bin,phi_und,phi_upp);
  TH1F *h_rho_max= new TH1F("h_rho_max","rho_max",rho_bin,rho_max_und,rho_max_upp);

  // for data about min
  TH1F *h_x_min= new TH1F("h_x_min","x_min",xy_bin,xy_und,xy_upp);
  TH1F *h_y_min= new TH1F("h_y_min","y_min",xy_bin,xy_und,xy_upp);
  TH1F *h_z_min= new TH1F("h_z_min","z_min",z_bin,z_und,z_upp);
  TH1F *h_psi_min= new TH1F("h_psi_min","psi_min",ang_bin,psi_und,psi_upp);
  TH1F *h_theta_min= new TH1F("h_theta_min","theta_min",ang_bin,theta_und,theta_upp);
  TH1F *h_phi_min= new TH1F("h_phi_min","phi_min",ang_bin,phi_und,phi_upp);
  TH2F *h_xy_min= new TH2F("h_xy_min","x-y min",xy_bin,xy_und,xy_upp,xy_bin,xy_und,xy_upp);
  TH3D *h_ang_min= new TH3D("h_ang_min","psi-theta-phi min",
			    ang_bin,psi_und,psi_upp,ang_bin,theta_und,theta_upp,ang_bin,phi_und,phi_upp);
  TH1F *h_rho_min= new TH1F("h_rho_min","rho_min",rho_bin,rho_min_und,rho_min_upp);
  
  //fill
  for(Int_t i=0; i<Entry1; i++){
    tree1->GetEntry(i);
    h_x_max->Fill(x_max);
    h_y_max->Fill(y_max);
    h_z_max->Fill(z_max);
    h_psi_max->Fill(psi_max);
    h_theta_max->Fill(theta_max);
    h_phi_max->Fill(phi_max);
    h_rho_max->Fill(rho_max-rho_max_cen);
    h_xy_max->Fill(x_max,y_max);
    h_ang_max->Fill(psi_max,theta_max,phi_max);
  }

  for(Int_t i=0; i<Entry2; i++){
    tree2->GetEntry(i);
    h_x_min->Fill(x_min);
    h_y_min->Fill(y_min);
    h_z_min->Fill(z_min);
    h_psi_min->Fill(psi_min);
    h_theta_min->Fill(theta_min);
    h_phi_min->Fill(phi_min);
    h_rho_min->Fill(rho_min-rho_min_cen);
    h_xy_min->Fill(x_min,y_min);
    h_ang_min->Fill(psi_min,theta_min,phi_min);
  }

  //make and save the hist
  TCanvas *can=new TCanvas("hists","hists");
  //  can->Divide(3,2);
  can->Divide(4,2);

  /*
  h_xy_max->SetStats(1);
  can->cd(1);
  h_xy_max->Draw("colz");

  h_z_max->SetStats(1);
  can->cd(2);
  h_z_max->Draw("hist");

  h_ang_max->SetStats(1);
  can->cd(3);
  h_ang_max->Draw("colz");

  h_xy_min->SetStats(1);
  can->cd(4);
  h_xy_min->Draw("colz");

  h_z_min->SetStats(1);
  can->cd(5);
  h_z_min->Draw("hist");

  h_ang_min->SetStats(1);
  can->cd(6);
  h_ang_min->Draw("colz");
  */
  
  h_psi_max->SetStats(1);
  can->cd(1);
  h_psi_max->Draw("hist");
  
  h_theta_max->SetStats(1);
  can->cd(2);
  h_theta_max->Draw("hist");
  
  h_phi_max->SetStats(1);
  can->cd(3);
  h_phi_max->Draw("hist");
  
  h_rho_max->SetStats(1);
  can->cd(4);
  h_rho_max->Draw("hist");

  /*
  h_x_max->SetStats(0);
  //  h_x_max->SetStats(1);
  can->cd(4);
  h_x_max->Draw("hist");
  
  h_y_max->SetStats(0);
  //  h_y_max->SetStats(1);
  can->cd(5);
  h_y_max->Draw("hist");
  
  h_z_max->SetStats(0);
  //  h_z_max->SetStats(1);
  can->cd(6);
  h_z_max->Draw("hist");
  */
  
  h_psi_min->SetStats(1);
  can->cd(5);
  h_psi_min->Draw("hist");
  
  h_theta_min->SetStats(1);
  can->cd(6);
  h_theta_min->Draw("hist");
  
  h_phi_min->SetStats(1);
  can->cd(7);
  h_phi_min->Draw("hist");
  
  h_rho_min->SetStats(1);
  can->cd(8);
  h_rho_min->Draw("hist");

  /*
  h_x_min->SetStats(0);
  //  h_x_min->SetStats(1);
  can->cd(10);
  h_x_min->Draw("hist");
  
  h_y_min->SetStats(0);
  //  h_y_min->SetStats(1);
  can->cd(11);
  h_y_min->Draw("hist");
  
  h_z_min->SetStats(0);
  //  h_z_min->SetStats(1);
  can->cd(12);
  h_z_min->Draw("hist");
  */
  
  //output for canvas
  TString outcan="data1/hists_check_rho_"+ to_string(subnum) +".C";
  can->SaveAs(outcan);

  cout<<"output: "<<outcan<<endl;
  return 0;
}
