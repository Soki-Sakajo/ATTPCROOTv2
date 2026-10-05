#include <iostream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TGraphErrors.h>
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

TGraph* read_crosssection(TString crossFile);

void comp_cross_12c12c_gsgs_cm90(){
   //read file (theta_cm, theta3_lab, E3. theta4_lab, E4)

   TString infi1 = "c12c12_gsgs_theta_90_Ecm_mb_sr_10-38MeV.txt"; //H.Emling et.al NPA 211 (1973) 600-616 fig.3
   TString infi2 = "c12c12_gsgs_theta_90_Ecm_mb_sr_20-30MeV.txt"; //R.Wieland et.al PRC vol.8 num.1 (1973) 37-45 fig.8
   TString infi3 = "c12c12_gsgs_theta_90_Ecm_mb_sr_13-38MeV.txt"; //W.Reilly et.al Nuovo Cimento A 13 (1973) 913-922 fig.2
   TString infi4 = "c12c12_gsgs_theta_90_Ecm_mb_sr_5-36MeV.txt";  //Y.Kucuk et.al NPA 764 (2006) 160-180 fig.6

   //hist difinition
   TGraph *emlin = read_crosssection(infi1);
   TGraph *wiela = read_crosssection(infi2);
   TGraph *reill = read_crosssection(infi3);
   TGraph *kucuk = read_crosssection(infi4);

   //write
   /*
   emlin->SetLineColor(kRed);
   emlin->SetLineWidth(2);
   emlin->SetMarkerStyle(21);
   wiela->SetLineColor(kBlue);
   wiela->SetLineWidth(2);
   wiela->SetMarkerStyle(4);
   reill->SetLineColor(kViolet);
   reill->SetLineWidth(2);
   reill->SetMarkerStyle(20);
   kucuk->SetLineColor(kBlack);
   kucuk->SetLineWidth(2);
   kucuk->SetMarkerStyle(3);

   TCanvas* c1 = new TCanvas("c1","c1");
   c1->cd();
   kucuk->GetXaxis()->SetRangeUser(3, 39);
   kucuk->GetYaxis()->SetRangeUser(5e-3, 3e+3);
   kucuk->SetTitle("E_{c.m.} and cross section of ^{12}C + ^{12}C at #theta_{c.m.} = 90 deg");
   kucuk->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
   kucuk->GetYaxis()->SetTitle("Cross section (d#sigma/d#Omega) [mb/sr]");
   kucuk->Draw("AL");
   emlin->Draw("L same");
   wiela->Draw("L same");
   reill->Draw("L same");
   gPad->SetLogy();

   TLegend *leg1 = new TLegend(0.68, 0.68, 0.90, 0.90);
   leg1->SetTextFont(42);
   leg1->SetTextSize(0.035);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(1);
   leg1->AddEntry(emlin, "H.Emling et.al", "L");
   leg1->AddEntry(wiela, "R.Wieland et.al", "L");
   leg1->AddEntry(reill, "W.Reilly et.al", "L");
   leg1->AddEntry(kucuk, "Y.Kucuk et.al", "L");
   leg1->Draw("same");
   c1->Update();
   */

   emlin->SetLineWidth(2);
   emlin->SetMarkerStyle(20);
   emlin->SetLineColor(kRed);
   emlin->SetMarkerColor(kRed);
   wiela->SetLineWidth(2);
   wiela->SetMarkerStyle(20);
   wiela->SetLineColor(kBlue);
   wiela->SetMarkerColor(kBlue);
   reill->SetLineWidth(2);
   reill->SetMarkerStyle(20);
   reill->SetMarkerSize(0.8);
   reill->SetLineColor(kBlack);
   reill->SetMarkerColor(kBlack);

   TCanvas* c2 = new TCanvas("c2","c2");
   c2->cd();
   emlin->GetXaxis()->SetRangeUser(10, 39);
   emlin->GetYaxis()->SetRangeUser(5e-3, 3e+3);
   emlin->SetTitle("E_{c.m.} and cross section of ^{12}C + ^{12}C at #theta_{c.m.} = 90 deg");
   emlin->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
   emlin->GetYaxis()->SetTitle("Cross section (d#sigma/d#Omega) [mb/sr]");
   /*
   emlin->Draw("AL");
   wiela->Draw("L same");
   reill->Draw("L same");
   */
   emlin->Draw("APL");
   wiela->Draw("PL same");
   reill->Draw("PL same");
   gPad->SetLogy();

   TLegend *leg2 = new TLegend(0.68, 0.70, 0.90, 0.90);
   leg2->SetTextFont(42);
   leg2->SetTextSize(0.035);
   leg2->SetFillStyle(0);
   leg2->SetBorderSize(1);
   leg2->AddEntry(emlin, "H.Emling et.al", "L");
   leg2->AddEntry(wiela, "R.Wieland et.al", "L");
   leg2->AddEntry(reill, "W.Reilly et.al", "L");
   leg2->Draw("same");
   c2->Update();

}

TGraph* read_crosssection(TString crossFile){

   Int_t n_cs = 0;
   Double_t E, cs;
   std::vector<std::vector<Double_t>> plot(0, std::vector<Double_t>(0,0));
   std::vector<std::pair<Double_t, Double_t>> data;

   #ifdef debug_mode
   // Debug: Check if file exists and current working directory
   std::cout << "DEBUG ReadKinematics: pwd = " << gSystem->pwd() << std::endl;
   std::cout << "DEBUG ReadKinematics: Attempting to read: " << crossFile << std::endl;
   if (gSystem->AccessPathName(crossFile.Data(), kFileExists)) {
      std::cout << "DEBUG ReadKinematics: FILE NOT FOUND!" << std::endl;
   } else {
      std::cout << "DEBUG ReadKinematics: File found" << std::endl;
   }
#endif
   std::ifstream *crossStr = new std::ifstream(crossFile.Data());

   if (!crossStr -> is_open() || !crossStr -> good()){
      std::cout << "Warning : No kinematics file found or cannot read: " << crossFile << std::endl;
      return new TGraph();
   }
   string line;
   while (getline(*crossStr, line)){
      if(line.empty() || line[0] == '#'){
         continue;
      }
      istringstream iss(line);
      if (iss >> E >> cs){
         data.emplace_back(E, cs);
      }
      else{
      std::cerr << "failed to parse line: " << line << std::endl;
      }
   }
#ifdef debug_mode
      std::cout << "DEBUG ReadKinematics: Successfully read " << data.size() << " points" << std::endl;
#endif
   if (data.size() == 0){
      std::cout << " Warning : No data read from kinematics file: " << crossFile << std::endl;
   }
   n_cs = data.size();
   sort(data.begin(), data.end(),[](const auto &a, auto &b){return a.first > b.first;});
#ifdef debug_mode
   std::cout << "DEBUG ReadKinematics: check vector size of Ecm : " << n_cs << std::endl;
#endif

   plot.resize(2, vector<Double_t> (0));
   for (Int_t i = 0; i < n_cs; i++){
      plot.at(0).push_back(data[i].first);
      plot.at(1).push_back(data[i].second);
   }
   if(n_cs != plot.at(0).size()){
      std::cout << "something wrong about filling data in making TGraph of cross section!!!"<<std::endl;
      std::cout << "  skipping drawing of cross section!!!"<<std::endl;
   }
   TGraph *gcs = new TGraph(n_cs, plot.at(0).data(), plot.at(1).data());
   return gcs;
}
