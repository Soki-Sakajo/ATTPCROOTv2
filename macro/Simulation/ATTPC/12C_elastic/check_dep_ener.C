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

int check_dep_ener(Int_t subnum=30){
  
  //set timer
  TStopwatch timer;
  timer.Start();

  //set situation
  Int_t dz=1;           // mm; bin of depth
  Double_t da=1;        // deg; bin of angle
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t l_attpc=1000;   // mm; depth of attpc

  Double_t dz2=1;           // mm; bin of depth
  Int_t nbin_zn2=int(140/dz2+0.9); // bin of z-reaction_rate hist
  Double_t de=0.1;
  Int_t nbin_e=int(8/de+0.9); // bin of z-reaction_rate hist
  
  //set constants
  Int_t ctom=10;
  Double_t pi=3.1415926535;
  Double_t deg_to_rad=pi/180;
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_a=(int)(180/da+0.9);
  
  // for 12C12C
  Int_t n1;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];

  /*
  //  TString refi="data1/sim_macro_12c_nd_0.root";
  TString refi="data2/sim_macro_12c_130matm_0.root";
  TFile* file =new TFile(refi);
  TTree *tree1 = (TTree*)file->Get("anadata");
  */

  // Tchain of 12c12c
  TString dir="data2";
  TChain *tree1=new TChain("anadata","anadata");
  for(int i=0; i<30; i++){
    TString refi="sim_macro_12c_130matm_" + to_string(i) + ".root";
    //    TString refi="sim_macro_12c_70matm_" + to_string(i) + ".root";
    TString direfi=dir+"/"+refi;
    tree1->Add(direfi);
    if(i==0){cout<<"read from  "<<refi<<flush;}
    if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
  }
  cout<<"end reading files"<<endl<<endl;
  
  // set reading data
  //  for 12C12C
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("n",&n1);
  tree1->SetBranchAddress("TotalE",&E1);

  //set histgram
  TH2F *h_dep1=new TH2F("h_dep1","reaction position: beam energy",nbin_zn,0,400, 104,0,52);
  //  TH2F *h_dep1=new TH2F("h_dep1","reaction position: beam energy",1000,0,1000,104,0,52);
  //  TH2F *h_dep2=new TH2F("h_dep2","reaction position: reaction energy of cm",nbin_zn,0,400,52,0,26);
  TH2F *h_dep2=new TH2F("h_dep2","reaction position: beam energy",nbin_zn,0,400, 1040,0,52);  
  TH2F *h_dep3=new TH2F("h_dep3","reaction position: reaction energy of cm",nbin_zn2,0,140,nbin_e,17,25);

  //make event lists  for 12c12c
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    if(n1==2){
      z1[0]=startz1[0]*ctom;
      z1[1]=startz1[1]*ctom;
      //      z1[0]=startz1[0];
      //      z1[1]=startz1[1];
      h_dep1->Fill(z1[0],E1);
      //      h_dep2->Fill(z1[0],E1/2);
      h_dep2->Fill(z1[0],E1);
      h_dep3->Fill(z1[0],E1/2);
    }
    if(j%1000==0){
      //need
      cout<<"Filling 12c data "<<100*j/Entry1 <<" % !    \r"<<flush;
    }
  }

  /*
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,1);
  
  //  h_ztheta_12c->SetStats(0);
  h_dep1->SetStats(1);
  can->cd(1);
  h_dep1->Draw("colz");

  //  h_ztheta_12c->SetStats(0);
  h_dep2->SetStats(1);
  can->cd(2);
  h_dep2->Draw("colz");

  h_dep3->SetStats(1);
  can->cd(3);
  h_dep3->Draw("colz");
  */

  TCanvas *can=new TCanvas("hists","hists",1500,900);
  //  TCanvas *can=new TCanvas("hists","hists");
  can->cd(1);
  h_dep1->SetStats(0);
  can->SetLeftMargin(0.1);
  can->SetRightMargin(0.17);
  can->SetTopMargin(0.01);
  can->SetBottomMargin(0.13);
  h_dep1->GetXaxis()->SetTitle("Depth in AT-TPC [mm]");
  h_dep1->GetYaxis()->SetTitle("Energy of ^{12}C [MeV]");
  h_dep1->GetZaxis()->SetTitle("number of events");
  h_dep1->Draw("colz");

  h_dep1->GetXaxis()->SetLabelSize(0.06);
  h_dep1->GetYaxis()->SetLabelSize(0.06);
  h_dep1->GetZaxis()->SetLabelSize(0.06);
  //set title
  h_dep1->GetXaxis()->SetTitleSize(0.06);
  h_dep1->GetYaxis()->SetTitleSize(0.06);
  h_dep1->GetYaxis()->SetTitleOffset(0.65);
  h_dep1->GetZaxis()->SetTitleSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex;
  latex.SetTextSize(0.06);
  latex.SetTextColor(kBlack);
  latex.DrawLatex(220,45.5,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  can->Update();

  //  TString outcan="test.pdf";
  TString outcan="./pdf_output/can_check_dep_ener.pdf";
  can->SaveAs(outcan);

  TCanvas *can2=new TCanvas("hists2","hists2",1500,900);
  //  TCanvas *can=new TCanvas("hists","hists");
  can2->cd(1);
  h_dep2->SetStats(0);
  can2->SetLeftMargin(0.1);
  can2->SetRightMargin(0.17);
  can2->SetTopMargin(0.01);
  can2->SetBottomMargin(0.13);
  h_dep2->GetXaxis()->SetTitle("Depth in AT-TPC [mm]");
  h_dep2->GetYaxis()->SetTitle("Beam Energy [MeV]");
  h_dep2->GetZaxis()->SetTitle("number of events");
  h_dep2->Draw("colz");
  h_dep2->GetXaxis()->SetLabelSize(0.06);
  h_dep2->GetYaxis()->SetLabelSize(0.06);
  h_dep2->GetZaxis()->SetLabelSize(0.06);
  //set title
  h_dep2->GetXaxis()->SetTitleSize(0.06);
  h_dep2->GetYaxis()->SetTitleSize(0.06);
  h_dep2->GetYaxis()->SetTitleOffset(0.65);
  h_dep2->GetZaxis()->SetTitleSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.06);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(220,45.5,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  can->Update();

  /*
  h_dep3->SetStats(1);
  h_dep3->Draw("colz");
  */
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

  TFile *file_out = new TFile("dep_ene.root", "RECREATE");
  h_dep2->Write();
  file_out->Close();

  return 0;
}
