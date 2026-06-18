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

TGraph* kine_read(TString kineFile);

void sum_kine_gsgs(){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  std::cout << std::endl;

  TGraph *kine_gsgs_0 =   kine_read("kine_12c12c_gsgs_60.7_verz_0.txt"  );
  TGraph *kine_gsgs_50 =  kine_read("kine_12c12c_gsgs_57.9_verz_50.txt" );
  TGraph *kine_gsgs_100 = kine_read("kine_12c12c_gsgs_55.0_verz_100.txt");
  TGraph *kine_gsgs_150 = kine_read("kine_12c12c_gsgs_51.8_verz_150.txt");
  TGraph *kine_gsgs_200 = kine_read("kine_12c12c_gsgs_48.5_verz_200.txt");
  TGraph *kine_gsgs_300 = kine_read("kine_12c12c_gsgs_41.5_verz_300.txt");
  TGraph *kine_gsgs_350 = kine_read("kine_12c12c_gsgs_37.8_verz_350.txt");
  TGraph *kine_gsgs_400 = kine_read("kine_12c12c_gsgs_33.8_verz_400.txt");
  TGraph *kine_gsgs_450 = kine_read("kine_12c12c_gsgs_29.5_verz_450.txt");
  TGraph *kine_gsgs_500 = kine_read("kine_12c12c_gsgs_24.6_verz_500.txt");

  //set output for making rutdata by each subnum
  TString outdata = "sum_kine_gsgs_10.root";
  TFile *file=new TFile(outdata, "recreate");

  TCanvas *c1 = new TCanvas("c1","c1");
  /*
  kine_gsgs_0->SetMarkerStyle(20);
  kine_gsgs_0->SetMarkerSize(1);
  kine_gsgs_0->SetMarkerColor(kBlue);
  */
  kine_gsgs_0->SetLineWidth(2);
  kine_gsgs_50->SetLineWidth(2);
  kine_gsgs_100->SetLineWidth(2);
  kine_gsgs_150->SetLineWidth(2);
  kine_gsgs_200->SetLineWidth(2);
  kine_gsgs_300->SetLineWidth(2);
  kine_gsgs_350->SetLineWidth(2);
  kine_gsgs_400->SetLineWidth(2);
  kine_gsgs_450->SetLineWidth(2);
  kine_gsgs_500->SetLineWidth(2);
  kine_gsgs_0->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  kine_gsgs_0->GetYaxis()->SetTitle("E_{beam} [MeV]");
  kine_gsgs_0->SetTitle("#theta_{lab} vs E_{beam} (vertex z = 100 mm)");
  kine_gsgs_0->  Draw("AL");
  kine_gsgs_50-> Draw("AL same");
  kine_gsgs_100->Draw("L same");
  kine_gsgs_150->Draw("L same");
  kine_gsgs_200->Draw("L same");
  kine_gsgs_300->Draw("L same");
  kine_gsgs_350->Draw("L same");
  kine_gsgs_400->Draw("L same");
  kine_gsgs_450->Draw("L same");
  kine_gsgs_500->Draw("L same");

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
  kine_gsgs_0   -> Write("kine_gsgs_z0"  );
  kine_gsgs_50  -> Write("kine_gsgs_z50" );
  kine_gsgs_100 -> Write("kine_gsgs_z100"); 
  kine_gsgs_150 -> Write("kine_gsgs_z150"); 
  kine_gsgs_200 -> Write("kine_gsgs_z200"); 
  kine_gsgs_300 -> Write("kine_gsgs_z300"); 
  kine_gsgs_350 -> Write("kine_gsgs_z350"); 
  kine_gsgs_400 -> Write("kine_gsgs_z400"); 
  kine_gsgs_450 -> Write("kine_gsgs_z450"); 
  kine_gsgs_500 -> Write("kine_gsgs_z500"); 
  file->Close();

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl;
  std::cout << "output: " << outdata << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << std::endl << std::endl;

}

TGraph* kine_read(TString kineFile){
  //set para
  char comma;
  Int_t n1, n2, n3;
  Double_t E, sigma;
  vector<vector<Double_t>> plot(0,vector<Double_t>(0,0));
  std::vector<std::pair<Double_t, Double_t>> data1;
  std::vector<std::pair<Double_t, Double_t>> data2;
  std::vector<std::pair<Double_t, Double_t>> data3;

  //input
  ifstream infi1(kineFile);
  if(!infi1.is_open()){
    std::cerr << "failed to open the file" << std::endl;
    //    return 1;
  }

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
  //  std::cout << "finished reading data from txt files." << std::endl << std::endl;

  n1=data1.size();
  sort(data1.begin(), data1.end(),[](const auto &a, auto &b){return a.first > b.first;});
  std::cout << "check vector size of data : " << n1 << std::endl;
  if (n1 == 0){
    std::cout << "data size is 0. end this macro!!" << std::endl;
    gSystem->Exit(0);
  }

  plot.resize(6,vector<Double_t> (0));
  //like fill for TGraph
  for(Int_t i=0; i<n1; i++){
    plot[0].push_back(data1.at(i).first);
    plot[1].push_back(data1.at(i).second);
  }

  if(n1 != plot.at(0).size()){
    std::cout << "something wrong about filling data!!! finish this macro"<<std::endl;
    gSystem->Exit(0);
  }

  TGraph *kine = new TGraph(n1, plot.at(0).data(), plot.at(1).data());

  return kine;
}
