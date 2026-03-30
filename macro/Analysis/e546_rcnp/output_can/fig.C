#include <iostream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TH1.h>
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
#include"TCanvas.h"

void fig(Int_t start=52, Int_t stop=56){ 
  
  //input
  TString infi="can_raw_eve_run" + to_string(start) + "-"+to_string(stop)+".root";
  TFile *file= new TFile(infi, "READ");
  if(!file || file->IsZombie()){
    cerr<<"Error: file not found or cannot be opened"<<endl;
    return;
  }
  TH1D *h =(TH1D*)file->Get("histADC");
  h->SetDirectory(0);
  TH1F *h2 =(TH1F*)file->Get("histADCvT");
  h2->SetDirectory(0);
  file->Close();

  /*
  Int_t xnbin=h2->GetNbinsX();
  TH1F *h3 = new TH1F("h3","h3",xnbin,0,h->GetXaxis()->GetXmax()-10);
  for(int i=1; i< xnbin+1; i++){
    Double_t content = h->GetBinContent(i);
    Double_t x = h->GetBinCenter(i)-10;
    h_new->Fill(x, content)
  }
  */

  //make and save the hist 
  /*
  TCanvas *can=new TCanvas("can","can",800,600);
  can->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.16);
  h->SetStats(0);
  //  h->SetTitle("ADC hist");
  h->GetXaxis()->SetTitle("ADC");
  h->GetYaxis()->SetTitle("");
  h->Draw("hist");
  h->SetTitleSize(0.08);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetXaxis()->SetRangeUser(-50,450);
  //  h->GetXaxis()->SetNdivisions(505);
  //  h->SetMinimum(1e5);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex1;
  latex1.SetTextSize(0.05);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(240,1.3e9,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.08 atm}");
  //  latex1.DrawLatex(240,8e8,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.13 atm}");
  can->Update();
  */

  TCanvas *c2=new TCanvas("c2","c2",800,600);
  c2->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.16);
  h2->SetStats(0);
  //  h2->SetTitle("ADC hist");
  h2->GetXaxis()->SetTitle("T [#mus]");
  h2->GetYaxis()->SetTitle("#SigmaADC");
  h2->Draw();
  h2->SetTitleSize(0.08);
  h2->GetXaxis()->SetLabelSize(0.06);
  h2->GetYaxis()->SetLabelSize(0.06);
  h2->GetXaxis()->SetTitleSize(0.06);
  h2->GetYaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleOffset(0.9);
  h2->GetYaxis()->SetTitleOffset(1.15);
  h2->GetXaxis()->SetRangeUser(-50,450);
  //  h2->GetXaxis()->SetNdivisions(505);
  h2->SetMinimum(3e6);
  h2->SetMaximum(7e9);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.05);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(2,3e9,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.08 atm}");
  //  latex2.DrawLatex(2,3e9,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.13 atm}");
  c2->Update();

  /*
  TCanvas *c3=new TCanvas("c3","c3",800,600);
  c3->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.16);
  h3->SetStats(0);
  //  h3->SetTitle("ADC hist");
  h3->GetXaxis()->SetTitle("T [#mus]");
  h3->GetYaxis()->SetTitle("#SigmaADC");
  h3->Draw("hist");
  h3->SetTitleSize(0.08);
  h3->GetXaxis()->SetLabelSize(0.06);
  h3->GetYaxis()->SetLabelSize(0.06);
  h3->GetXaxis()->SetTitleSize(0.06);
  h3->GetYaxis()->SetTitleSize(0.06);
  h3->GetXaxis()->SetTitleOffset(0.9);
  h3->GetYaxis()->SetTitleOffset(0.8);
  h3->GetXaxis()->SetRangeUser(0,450);
  //  h3->GetXaxis()->SetNdivisions(505);
  h3->SetMinimum(1e6);
  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex2;
  latex2.SetTextSize(0.05);
  latex2.SetTextColor(kBlack);
  latex2.DrawLatex(240,1.3e9,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.08 atm}");
  //  latex2.DrawLatex(240,8e8,"#splitline{^{12}C beam @ 60.7 MeV}{   iso-C_{4}H_{10}: 0.13 atm}");
  c3->Update();
  */
}
