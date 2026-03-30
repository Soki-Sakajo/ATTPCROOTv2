#include <iostream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
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

int pressure_sum(){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set val
  vector<vector<Double_t>> padund(0,vector<Double_t>(0,0));
  vector<vector<Double_t>> elasti(0,vector<Double_t>(0,0));
  vector<vector<Double_t>> detect(0,vector<Double_t>(0,0));

  //set output for making rutdata by each subnum
  TString outdata = "data1/pressure_sum.root";
  TFile *file=new TFile(outdata, "recreate");
  TTree *data = new TTree("data","data");
  data->Branch("under_250mm",&padund);
  data->Branch("over_100hz_rmax",&elasti);
  data->Branch("detection_efficiency",&detect);

  //output
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
  TString outcan1="./pdf_output/can_gas_rmax.pdf";
  TString outcan2="./pdf_output/can_gas_elastic.pdf";
  TString outcan3="./pdf_output/can_gas_detection.pdf";


  //set data
  /*
  vector<pair<Double_t,Int_t>> paddata={
    {0.01,88},{0.02,83},{0.03,84},{0.04,90},{0.05,96},{0.06,98},{0.07,99},
    {0.08,100},{0.09,100},{0.10,100},{0.11,100},{0.12,100},{0.13,100}
  };
  vector<pair<Double_t,Int_t>> eladata={
    {0.07,92},{0.08,77},{0.09,67},{0.10,63},{0.11,59},{0.12,51},{0.13,53}
  };
  vector<pair<Double_t,Double_t>> detdata={
    {0.07,55.9},{0.08,60.1},{0.09,61.6},{0.10,60.4},{0.11,59.7},{0.12,63.1},{0.13,57.4}
  };
  */
  // adapt results of fit
  vector<pair<Double_t,Int_t>> paddata={
    {0.01,47},{0.02,62},{0.03,68},{0.04,80},{0.05,92},{0.06,97},{0.07,99},
    {0.08,100},{0.09,100},{0.10,100},{0.11,100},{0.12,100},{0.13,100}
  };
  vector<pair<Double_t,Int_t>> eladata={
    {0.07,84},{0.08,74},{0.09,66},{0.10,60},{0.11,55},{0.12,49},{0.13,48}
  };

  vector<pair<Double_t,Double_t>> detdata={
    {0.07,60.74},{0.08,62.43},{0.09,62.65},{0.10,63.05},{0.11,63.47},{0.12,65.83},{0.13,63.54}
  };
  /*
  vector<pair<Double_t,Double_t>> detdata={
    {0.07,60.74},{0.08,62.43},{0.09,62.65},{0.10,63.05},{0.11,63.47}
  };
  */
  Int_t n_p=paddata.size();
  Int_t n_e=eladata.size();
  Int_t n_d=detdata.size();

  //like fill for TGraphErrors
  padund.resize(2,vector<Double_t>(n_p));
  for(Int_t i=0; i<n_p; i++){
    padund.at(0).at(i)=paddata[i].first;
    padund.at(1).at(i)=paddata[i].second;
  }
  elasti.resize(2,vector<Double_t>(n_e));
  for(Int_t i=0; i<n_e; i++){
    elasti.at(0).at(i)=eladata[i].first;
    elasti.at(1).at(i)=eladata[i].second;
  }
  detect.resize(2,vector<Double_t>(n_d));
  for(Int_t i=0; i<n_d; i++){
    detect.at(0).at(i)=detdata[i].first;
    detect.at(1).at(i)=detdata[i].second;
  }
  cout<<"Filled data"<<endl;

  TGraph *h_pad= new TGraph(n_p,padund.at(0).data(),padund.at(1).data());
  TGraph *h_ela= new TGraph(n_e,elasti.at(0).data(),elasti.at(1).data());
  TGraph *h_det= new TGraph(n_d,detect.at(0).data(),detect.at(1).data());

  //make and save the hist
  TCanvas *can1=new TCanvas("hist1","hist1",1500,1000);
  can1->cd(1);
  can1->SetLeftMargin(0.15);
  can1->SetRightMargin(0.05);
  can1->SetTopMargin(0.03);
  can1->SetBottomMargin(0.16);
  h_pad->SetTitle("");
  h_pad->SetMarkerStyle(20);
  h_pad->SetMarkerSize(3);
  h_pad->SetMarkerColor(kBlue);
  h_pad->GetXaxis()->SetTitle("Gas pressure [atm]");
  h_pad->GetYaxis()->SetTitle("P [%]");
  h_pad->GetXaxis()->SetLabelSize(0.08);
  h_pad->GetYaxis()->SetLabelSize(0.08);
  //set title
  h_pad->GetXaxis()->SetTitleSize(0.08);
  h_pad->GetYaxis()->SetTitleSize(0.08);
  h_pad->GetXaxis()->SetTitleOffset(0.9);
  h_pad->GetYaxis()->SetTitleOffset(0.9);
  gPad->Update();
  h_pad->Draw("AP");

  TCanvas *can2=new TCanvas("hist2","hist2",1500,1000);
  can2->cd(1);
  can2->SetLeftMargin(0.15);
  can2->SetRightMargin(0.05);
  can2->SetTopMargin(0.03);
  can2->SetBottomMargin(0.16);
  h_ela->SetTitle("");
  h_ela->SetMarkerStyle(20);
  h_ela->SetMarkerSize(3);
  h_ela->SetMarkerColor(kBlue);
  h_ela->GetXaxis()->SetTitle("Gas pressure [atm]");
  h_ela->GetYaxis()->SetTitle("r_{trg} [mm]");
  h_ela->GetXaxis()->SetLabelSize(0.08);
  h_ela->GetYaxis()->SetLabelSize(0.08);
  //set title
  h_ela->GetXaxis()->SetTitleSize(0.08);
  h_ela->GetYaxis()->SetTitleSize(0.08);
  h_ela->GetXaxis()->SetTitleOffset(0.9);
  h_ela->GetYaxis()->SetTitleOffset(0.8);
  gPad->Update();
  h_ela->Draw("AP");
  TLatex latex1;
  latex1.SetTextSize(0.08);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(0.093,83,"^{12}C beam @ 50 MeV");
  can2->Update();

  TCanvas *can3=new TCanvas("hist3","hist3",1500,1000);
  can3->cd(1);
  can3->SetLeftMargin(0.15);
  can3->SetRightMargin(0.05);
  can3->SetTopMargin(0.05);
  can3->SetBottomMargin(0.16);
  h_det->SetTitle("");
  h_det->SetMarkerStyle(20);
  //  h_det->SetMarkerSize(3);
  h_det->SetMarkerSize(2);
  h_det->SetMarkerColor(kBlue);
  h_det->GetXaxis()->SetTitle("Gas pressure [atm]");
  h_det->GetYaxis()->SetTitle("P_{6#alpha} [%]");
  h_det->GetXaxis()->SetLabelSize(0.08);
  h_det->GetYaxis()->SetLabelSize(0.08);
  //set X axis val
  //  h_det->GetXaxis()->SetNdivisions(5);
  //set title
  h_det->GetXaxis()->SetTitleSize(0.08);
  h_det->GetYaxis()->SetTitleSize(0.08);
  h_det->GetXaxis()->SetTitleOffset(0.9);
  h_det->GetYaxis()->SetTitleOffset(0.8);
  h_det->GetYaxis()->SetRangeUser(0,100);
  gPad->Update();
  h_det->Draw("AP");
  /*
  TLatex latex2;
  latex2.SetTextSize(0.08);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(0.067,65.7,"^{12}C beam @ 50 MeV");
  */
  can3->Update();

  can1->SaveAs(outcan1);
  can2->SaveAs(outcan2);
  can3->SaveAs(outcan3);

  //save branch
  data->AutoSave();
  file->Close();
  
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
