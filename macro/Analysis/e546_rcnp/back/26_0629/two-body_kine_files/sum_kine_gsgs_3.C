#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TFile.h>
#include <vector>
#include"TString.h"
#include"TFile.h"
#include"TTree.h"
#include"TH1.h"
#include"TH2.h"
#include"TCanvas.h"

void sum_kine_gsgs(){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  std::cout << std::endl;

  //set para
  char comma;
  Int_t n1, n2, n3;
  Double_t E, sigma;
  vector<vector<Double_t>> plot(0,vector<Double_t>(0,0));
  std::vector<std::pair<Double_t, Double_t>> data1;
  std::vector<std::pair<Double_t, Double_t>> data2;
  std::vector<std::pair<Double_t, Double_t>> data3;

  //input
  ifstream infi1("kine_12c12c_gsgs_55.0_verz_100.txt");
  ifstream infi2("kine_12c12c_gsgs_41.5_verz_300.txt");
  ifstream infi3("kine_12c12c_gsgs_33.8_verz_400.txt");
  if(!infi1.is_open() || !infi2.is_open() || !infi3.is_open()){
    std::cerr << "failed to open the file" << std::endl;
    //    return 1;
  }

  //set output for making rutdata by each subnum
  TString outdata = "sum_kine_gsgs_10.root";
  TFile *file=new TFile(outdata, "recreate");

  //set energy, range, cross section
  string line;
  while(getline(infi1, line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    istringstream iss(line);
    if(iss >> E >> sigma){
      data1.emplace_back(E, sigma);
    }
    else{
      std::cerr << "failed to parse line: " << line << std::endl;
    }
  }
  while(getline(infi2, line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    istringstream iss(line);
    if(iss >> E >> sigma){
      data2.emplace_back(E, sigma);
    }
    else{
      std::cerr << "failed to parse line: " << line << std::endl;
    }
  }
  while(getline(infi3, line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    istringstream iss(line);
    if(iss >> E >> sigma){
      data3.emplace_back(E, sigma);
    }
    else{
      std::cerr << "failed to parse line: " << line << std::endl;
    }
  }
  std::cout << "finished reading data from txt files." << std::endl << std::endl;

  n1=data1.size();
  n2=data2.size();
  n3=data3.size();
  sort(data1.begin(), data1.end(),[](const auto &a, auto &b){return a.first > b.first;});
  sort(data2.begin(), data2.end(),[](const auto &a, auto &b){return a.first > b.first;});
  sort(data3.begin(), data3.end(),[](const auto &a, auto &b){return a.first > b.first;});
  std::cout << "check vector size of data : " << n1 << ", " << n2 << ", " << n3 << std::endl;
  if (n1 == 0 || n2 == 0 || n3 == 0){
    std::cout << "data size is 0. end this macro!!" << std::endl;
    gSystem->Exit(0);
  }
  if (n1 != n2 || n1 != n3){
    std::cout << "data size is different. end this macro!!" << std::endl;
    gSystem->Exit(0);
  }

  plot.resize(6,vector<Double_t> (0));
  //like fill for TGraph
  for(Int_t i=0; i<n1; i++){
    plot[0].push_back(data1.at(i).first);
    plot[1].push_back(data1.at(i).second);
    plot[2].push_back(data2.at(i).first);
    plot[3].push_back(data2.at(i).second);
    plot[4].push_back(data3.at(i).first);
    plot[5].push_back(data3.at(i).second);
  }

  if(n1 != plot.at(0).size()){
    std::cout << "something wrong about filling data!!! finish this macro"<<std::endl;
    gSystem->Exit(0);
  }

  TGraph *kine_gsgs_100 = new TGraph(n1, plot.at(0).data(), plot.at(1).data());
  TGraph *kine_gsgs_300 = new TGraph(n2, plot.at(2).data(), plot.at(3).data());
  TGraph *kine_gsgs_400 = new TGraph(n3, plot.at(4).data(), plot.at(5).data());

  /*
  TCanvas *c1 = new TCanvas("c1","c1", 1200, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  kine_gsgs_100->SetMarkerStyle(20);
  kine_gsgs_100->SetMarkerSize(1);
  kine_gsgs_100->SetMarkerColor(kBlue);
  kine_gsgs_100->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  kine_gsgs_100->GetYaxis()->SetTitle("E_{beam} [MeV]");
  kine_gsgs_100->SetTitle("#theta_{lab} vs E_{beam} (vertex z = 100 mm)");
  kine_gsgs_100->Draw("AP");
  c1->cd(2);
  kine_gsgs_300->SetMarkerStyle(20);
  kine_gsgs_300->SetMarkerSize(1);
  kine_gsgs_300->SetMarkerColor(kBlue);
  kine_gsgs_300->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  kine_gsgs_300->GetYaxis()->SetTitle("E_{beam} [MeV]");
  kine_gsgs_300->SetTitle("#theta_{lab} vs E_{beam} (vertex z = 300 mm)");
  kine_gsgs_300->Draw("AP");
  c1->cd(3);
  kine_gsgs_400->SetMarkerStyle(20);
  kine_gsgs_400->SetMarkerSize(1);
  kine_gsgs_400->SetMarkerColor(kBlue);
  kine_gsgs_400->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  kine_gsgs_400->GetYaxis()->SetTitle("E_{beam} [MeV]");
  kine_gsgs_400->SetTitle("#theta_{lab} vs E_{beam} (vertex z = 400 mm)");
  kine_gsgs_400->Draw("AP");
  */

  // save lines
  kine_gsgs_100->Write("kine_12c_gsgs_z100");
  kine_gsgs_300->Write("kine_12c_gsgs_z300");
  kine_gsgs_400->Write("kine_12c_gsgs_z400");
  file->Close();

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl;
  std::cout << "output: " << outdata << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << std::endl << std::endl;

}
