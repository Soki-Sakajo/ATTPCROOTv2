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

int yeild_calculation_130matm_gs_hoyle(){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  //  Int_t dz=2;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  Double_t P=0.13;      // atm; gas pressure
  //  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc

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
  cout<<"rho: "<<rho*N_rate1<<" /m^3"<<endl;

  Int_t n=0;
  Int_t n_p=0;
  Double_t z_min,z_max, E,sigma,stat,syst,dif_z;
  Double_t z_tem[5];
  Double_t Ysum=0;
  Double_t Yesum=0;
  Double_t Ydasum=0;
  Double_t Yrunsum=0;
  Double_t Yrunesum=0;
  vector<vector<Double_t>> range_attpc(2,vector<Double_t>(0,0));
  vector<Double_t> dz(0,0);
  vector<Double_t> sig(0,0);
  vector<Double_t> sier(0,0);
  vector<Double_t> Ecm(0,0);
  vector<Double_t> Y1(0,0);
  vector<Double_t> Ye(0,0);
  vector<Double_t> Yda(0,0);
  vector<Double_t> Yrun(0,0);
  vector<Double_t> Yrunerr(0,0);
  vector<vector<Double_t>> plot(0,vector<Double_t>(0,0));
  
  //set histgram
  TH1F *h_sig=new TH1F("h_sig","Ecm :cross-section in fujikawa experiment",85,17,25.5);
  TH1F *h_dif=new TH1F("h_dif","Ecm :thick of gass (dz)",85,17,25.5);
  //  TH1F *h_yei=new TH1F("h_yei","Ecm :yeild of attpc experiment",85,17,25.5);
  TH1F *h_yei=new TH1F("h_yei","Ecm :yeild of attpc experiment",90,16.5,25.5);
  TH2F *h_ran=new TH2F("h_ran","Ecm :range of attpc",85,17,25.5,150,0,150);
  TH1F *h_yeda=new TH1F("h_yeda","Ecm :count of 1day attpc experiment",85,17,25.5);
  TH1F *h_yeru=new TH1F("h_yeru","Ecm :count attpc experiment",85,17,25.5);
  //  TH1F *h_yeru=new TH1F("h_yeru","Ecm :count attpc experiment",70,17,25.5);

  //input
  /*
  ifstream infi("constant_yeild_calculation.dat");
  if(!infi.is_open()){
    cerr<<"failed to open the file"<<endl;
    return 1;
  }
  */
  ifstream infi("constant_yeild_calculation_gs_hoyle_dz.dat");
  if(!infi.is_open()){
    cerr<<"failed to open the file"<<endl;
    return 1;
  }

  //set output for making rutdata by each subnum
  TString outdata = "data1/yeild_calculation_130matm_gs_hoyle.root";
  TFile *file=new TFile(outdata, "recreate");
  TTree *data = new TTree("data","data");
  data->Branch("number_of_data",&n,"n/I");
  data->Branch("reaction_energy",&Ecm);
  data->Branch("cross_section",&sig);
  data->Branch("range",&range_attpc);
  data->Branch("dz",&dz);
  data->Branch("yeild",&Y1);
  data->Branch("error_of_yeild",&Ye);
  data->Branch("count",&Yrun);
  data->Branch("error_of_count",&Yrunerr);

  cout<<endl;

  //set energy, range, cross section
  string line;
  while(getline(infi, line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    istringstream iss(line);
    //    if(iss >> E >> z_min >> z_max >> sigma >> stat >> syst ){
    if(iss >> E >> z_min >> sigma >> stat >> syst ){
      //      cout<<E<<" "<<z_min<<" "<<z_max<<" "<<sigma<<" "<<stat<<" "<<syst<<endl;
      Ecm.push_back(E);
      range_attpc[0].push_back(z_min);
      range_attpc[1].push_back(z_max);
      dz.push_back(z_min);
      sig.push_back(sigma*1e-3);
      sier.push_back((stat+syst)*1e-3);
    }
    else{
      cerr<<"faild to parse line: "<<line<<endl;
    }
  }

  cout<<"check vector size of Ecm : "<<Ecm.size()<<endl;
  n=Ecm.size();
  if(n==0){
    cout<<"Ecm size is 0. end this macro!!"<<endl;
    gSystem->Exit(0);
  }
  Y1.resize(n);
  Ye.resize(n);
  Yda.resize(n);
  Yrun.resize(n);
  Yrunerr.resize(n);
  
  for(Int_t i=0; i<n; i++){
    //    cout<<i<<" "<<Ecm.at(i)<<" "<<dz.at(i)<<" "<<rho*dz.at(i)*1e-3*N_rate1/6e23*12*1e3*1e-4<<endl;
    Y1.at(i)=sig.at(i)*n_b*rho*dz.at(i)*1e-34*N_rate1;
    Ye.at(i)=sier.at(i)*n_b*rho*dz.at(i)*1e-34*N_rate1;
    Yda.at(i)=Y1.at(i)*24*3600;
    Yrun.at(i)=Y1.at(i)*24*3600*5;
    Yrunerr.at(i)=Ye.at(i)*24*3600*5;
    //    h_sig->Fill(Ecm.at(i),sig.at(i));
    h_sig->Fill(Ecm.at(i),sig.at(i)*1e3);
    h_ran->Fill(Ecm.at(i),range_attpc.at(0).at(i));
    h_ran->Fill(Ecm.at(i),range_attpc.at(1).at(i));
    h_dif->Fill(Ecm.at(i),dz.at(i));
    h_yei->Fill(Ecm.at(i),Y1.at(i));
    h_yeru->Fill(Ecm.at(i),Yrun.at(i));
    h_yeda->Fill(Ecm.at(i),Yda.at(i));
    if(Y1.at(i)>0){
      n_p++;
    }
    Ysum+=Y1.at(i);
    Yesum+=Ye.at(i);
    Ydasum+=Yda.at(i);
    Yrunsum+=Yrun.at(i);
    Yrunesum+=Yrunerr.at(i);
    /*
    if(Ecm.at(i)==22.2){
      cout<<i<<" energy "<<Ecm.at(i)<<" sig "<<sig.at(i)<<" n_b "<<n_b<<"  rho "<<rho*N_rate1
	  <<" range_attpc "<<range_attpc.at(0).at(i)<<" "<<range_attpc.at(1).at(i)<<" dz "<<dz.at(i)
	  <<" yeild "<<Y1.at(i)<<endl;
      cout<<endl<<endl<<endl;
      cout<<"rho dz: "<<rho*N_rate1*dz.at(i)*1e-3*1e-4<<endl;
    }
    */
    dif_z=dz.at(i)/5;
    for(Int_t j=0; j<5; j++){
      z_tem[i]=range_attpc.at(0).at(i)+dif_z*j;
      h_ran->Fill(Ecm.at(i),z_tem[i]);
    }
  }

  plot.resize(5,vector<Double_t> (0));
  //like fill for TGraphErrors
  for(Int_t i=0; i<n; i++){
    //    if(Y1.at(i)>1e-6){
    if(Y1.at(i)>0){
      plot[0].push_back(Ecm.at(i));
      plot[1].push_back(0);
      plot[2].push_back(Y1.at(i));
      plot[3].push_back(Ye.at(i));
      plot[4].push_back(Ye.at(i));
      //      cout<<Ecm.at(i)<<endl;
    }
  }

  if(n_p!=plot.at(0).size()){
    cout<<endl;
    cout<<"something wrong about filling data!!! finish this macro"<<endl;
    data->AutoSave();
    file->Close();
    gSystem->Exit(0);
  }

  cout<<endl;
  //  cout<<"sum of yeild "<<Ysum<<", sum of day count "<<Ydasum<<", sum of all count "<<Yrunsum<<endl;
  cout<<"sum of yeild "<<Ysum<<" error "<<Yesum<<", all count "<<Yrunsum<<" error "<<Yrunesum<<endl;

  for(Int_t i=0; i<plot.at(0).size(); i++){
    if(plot.at(2).at(i)-plot.at(3).at(i)<0){
      cout<<"error bar is under 0 ! energy: "<<plot.at(0).at(i)<<", Yeild: "<<plot.at(2).at(i)
	  <<", Error: "<<plot.at(3).at(i)<<endl;
      plot.at(3).at(i)=plot.at(2).at(i);
    }
  }

  //  TGraphErrors *h_yp=new TGraphErrors(n_p,plot.at(0).data(),plot.at(2).data(),plot.at(1).data(),plot.at(3).data());

  TGraphAsymmErrors *h_yp=
    new TGraphAsymmErrors(n_p,plot.at(0).data(),plot.at(2).data(),
			  plot.at(1).data(),plot.at(1).data(),plot.at(3).data(),plot.at(4).data());

  /*
  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,2);
  //  can->Divide(2,1);
  
  can->cd(1);
  h_sig->SetStats(0);
  //  h_sig->SetStats(1);
  h_sig->Draw("hist");
  h_sig->GetXaxis()->SetLabelSize(0.06);
  h_sig->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(2);
  h_ran->SetStats(0);
  //  h_ran->SetStats(1);
  h_ran->Draw("colz");
  h_ran->GetXaxis()->SetLabelSize(0.06);
  h_ran->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(3);
  h_dif->SetStats(0);
  //  h_dif->SetStats(1);
  h_dif->Draw("hist");
  h_dif->GetXaxis()->SetLabelSize(0.06);
  h_dif->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  TCanvas *can2=new TCanvas("hist2","hist2");
  can2->cd(1);
  //  can->cd(4);
  //  can->cd(2);
  h_yei->SetStats(0);
  //  h_yei->SetStats(1);
  h_yei->Draw("hist");
  h_yei->GetXaxis()->SetLabelSize(0.06);
  h_yei->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(5);
  //  can->cd(3);
  h_yeda->SetStats(0);
  //  h_yeda->SetStats(1);
  h_yeda->Draw("hist");
  h_yeda->GetXaxis()->SetLabelSize(0.06);
  h_yeda->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can->cd(6);
  //  can->cd(4);
  h_yeru->SetStats(0);
  //  h_yeru->SetStats(1);
  h_yeru->Draw("hist");
  h_yeru->GetXaxis()->SetLabelSize(0.06);
  h_yeru->GetYaxis()->SetLabelSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  //output for canvas
  TString outcan = "data1/yeild_calculation_can_130matm.C";
  can->SaveAs(outcan);
  h_sig->Write();
  h_ran->Write();
  h_dif->Write();
  h_yei->Write();
  */
  //  TCanvas *can1=new TCanvas("hists","hists");
  //  TCanvas *can1=new TCanvas("hists","hists",1200,900);
  TCanvas *can1=new TCanvas("hists","hists",1500,900);
  //  can1->Divide(2,1);
  can1->cd(1);
  /*
  h_yei->SetStats(0);
  can1->SetLeftMargin(0.13);
  can1->SetRightMargin(0.01);
  can1->SetTopMargin(0.1);
  can1->SetBottomMargin(0.13);
  h_yei->GetXaxis()->SetTitle("Reaction Energy (Ecm) [MeV]");
  h_yei->GetYaxis()->SetTitle("Yeild [MeV]");
  //  h_yei->GetZaxis()->SetTitle("number of events");
  h_yei->Draw("hist");

  h_yei->GetXaxis()->SetLabelSize(0.06);
  h_yei->GetYaxis()->SetLabelSize(0.06);
  //  h_yei->GetZaxis()->SetLabelSize(0.06);
  //set title
  h_yei->GetXaxis()->SetTitleSize(0.06);
  h_yei->GetYaxis()->SetTitleSize(0.06);
  h_yei->GetYaxis()->SetTitleOffset(1);
  //  h_yei->GetZaxis()->SetTitleSize(0.06);
  //  gStyle->SetOptTitle(0);
  gPad->Update();
  TLatex latex1;
  latex1.SetTextSize(0.06);
  latex1.SetTextColor(kBlack);
  latex1.DrawLatex(17,3.5e-5,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  can1->Update();

  can1->cd(2);

  h_yp->SetTitle("yeild;E [MeV];Y(E) [cps]");
  h_yp->SetMarkerStyle(20);
  h_yp->SetMarkerSize(1.0);
  h_yp->SetMarkerColor(kBlue);
  h_yp->Draw("AP");
  can1->Update();
  */
  can1->SetLeftMargin(0.12);
  can1->SetRightMargin(0.05);
  can1->SetTopMargin(0.1);
  can1->SetBottomMargin(0.13);
  h_yp->SetTitle("yeild;E [MeV];Y(E) [cps]");
  h_yp->SetTitle("");
  h_yp->SetMarkerStyle(20);
  h_yp->SetMarkerSize(1.5);
  h_yp->SetMarkerColor(kBlue);

  h_yp->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  h_yp->GetYaxis()->SetTitle("Y(E) [cps]");
  h_yp->GetXaxis()->SetLabelSize(0.06);
  h_yp->GetYaxis()->SetLabelSize(0.06);
  //set title
  h_yp->GetXaxis()->SetTitleSize(0.06);
  h_yp->GetYaxis()->SetTitleSize(0.06);
  h_yp->GetYaxis()->SetTitleOffset(0.9);
  gPad->Update();
  h_yp->Draw("AP");
  TLatex latex;
  latex.SetTextSize(0.06);
  latex.SetTextColor(kBlack);
  latex.DrawLatex(17.2,3e-4,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  can1->Update();


  //  TString outcan1="test.pdf";
  TString outcan1="./pdf_output/can_yeild_gs_hoyle.pdf";
  can1->SaveAs(outcan1);


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
