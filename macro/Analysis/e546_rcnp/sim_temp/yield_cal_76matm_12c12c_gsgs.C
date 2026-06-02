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

std::vector<std::vector<Double_t>> cal_Ebeam_depth();
Double_t est_depth(std::vector<std::vector<Double_t>> &Ebeam_depth, Double_t Ecm);

void yield_cal_76matm_12c12c_gsgs(){
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t ang = 90;       // deg; angle of 12c12c scattering for cross section
  Int_t da = 3;         // deg; width of angle for cross section
  //  Int_t dz=2;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  //  Double_t P=0.13;      // atm; gas pressure
  Double_t P=0.076;      // atm; gas pressure
  //  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  Double_t Ebeamm = 60.8;

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

  char comma;
  Int_t n=0;
  Int_t n_p=0;
  Double_t z_min, z_max, E, sigma, stat, syst, dif_z;
  Double_t z_tem[5];
  Double_t Ecmm = 0;
  Double_t domega = 0;
  Double_t Ysum = 0;
  Double_t Yesum = 0;
  Double_t Ydasum = 0;
  Double_t Yrunsum = 0;
  Double_t Yrunesum = 0;
  //  Double_t dz = 1.0;
  vector<Double_t> dz(0,0);
  vector<Double_t> z(0,0);
  //  vector<Double_t> sig(0,0);
  //  vector<Double_t> sier(0,0);
  //  vector<Double_t> Ecm(0,0);
  vector<Double_t> Y1(0,0);
  //  vector<Double_t> Ye(0,0);
  vector<Double_t> Yda(0,0);
  vector<Double_t> Yrun(0,0);
  vector<Double_t> Ya(0,0);
  //  vector<Double_t> Yrunerr(0,0);
  vector<vector<Double_t>> plot(0,vector<Double_t>(0,0));
  //  vector<vector<Double_t>> range_attpc(2,vector<Double_t>(0,0));
  //  vector<vector<Double_t>> data(2,vector<Double_t>(0,0));
  vector<vector<Double_t>> Ebeam_depth(2,vector<Double_t>(0,0));
  std::vector<std::vector<Double_t>> v_data(1, vector<Double_t>(0,0));
  std::vector<std::pair<Double_t, Double_t>> data;
  
  //set histgram
  //  TH1F *h_sig=new TH1F("h_sig","Ecm :cross-section in thesis", 370, 0, 37);
  //  TH2F *h_ran=new TH2F("h_ran","Ecm :range of attpc",85,17,25.5,150,0,150);

  Ecmm = Ebeamm/2;
  domega = 2*pi * (cos((ang/2 - da) * deg_to_rad) - cos((ang/2 + da) * deg_to_rad));
  Ebeam_depth = cal_Ebeam_depth();

  //input
  ifstream infi(Form("c12c12_gsgs_theta_%d_Ecm_mb_sr.txt", ang));
  if(!infi.is_open()){
    cerr<<"failed to open the file"<<endl;
    return 1;
  }
  /*
  //set output for making rutdata by each subnum
  // if you want to set the branch, you need to change tree name.
  TString outdata = "yeild_calculation_76matm.root";
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
  */

  cout<<endl;
  //set energy, range, cross section
  string line;
  while(getline(infi, line)){
    if(line.empty()||line[0]=='#'){
      continue;
    }
    istringstream iss(line);
    //    if(iss >> E >> z_min >> z_max >> sigma >> stat >> syst ){
    //    if(iss >> E >> z_min >> sigma >> stat >> syst ){
    if(iss >> E >> comma >> sigma){
      //      cout<<E<<" "<<z_min<<" "<<z_max<<" "<<sigma<<" "<<stat<<" "<<syst<<endl;
      /*
      Ecm.push_back(E);
      range_attpc[0].push_back(z_min);
      range_attpc[1].push_back(z_max);
      dz.push_back(z_min);
      sig.push_back(sigma);
      sier.push_back((stat+syst));
      data[0].push_back(E);
      data[1].push_back(sigma);
      */
      data.emplace_back(E, sigma);
    }
    else{
      cerr<<"faild to parse line: "<<line<<endl;
    }
  }

  n=data.size();
  sort(data.begin(), data.end(),[](const auto &a, auto &b){return a.first > b.first;});
  std::cout << "check vector size of Ecm : " << n << std::endl;
  if (n == 0){
    cout<<"Ecm size is 0. end this macro!!"<<endl;
    gSystem->Exit(0);
  }
  z.resize(n+1);
  dz.resize(n);
  Y1.resize(n);
  //  Ye.resize(n);
  Yda.resize(n);
  Yrun.resize(n);
  Ya.resize(n);
  //  Yrunerr.resize(n);
  
  for(Int_t i = 0; i < n; i++){
    z.at(i) = est_depth(Ebeam_depth, data.at(i).first);
  }
  z.at(n) = 680;
  std::cout << "end the calculation of z" << endl;
  for(Int_t i=0; i<n; i++){
    if(data.at(i).first > Ecmm){
      dz.at(i) = 0;
    }
    else {
      dz.at(i) = abs(z.at(i+1) - z.at(i));
    }
    //    cout<<i<<" "<<Ecm.at(i)<<" "<<dz.at(i)<<" "<<rho*dz.at(i)*1e-3*N_rate1/6e23*12*1e3*1e-4<<endl;
    //    Y1.at(i)=sig.at(i)*n_b*rho*dz.at(i)*1e-34*N_rate1;
    Y1.at(i) = data.at(i).second * domega * n_b * rho * dz.at(i) * 1e-34 * N_rate1;
    //  Ye.at(i)=sier.at(i)*n_b*rho*dz.at(i)*1e-34*N_rate1;
    Yda.at(i)=Y1.at(i)*24*3600;
    Yrun.at(i)=Y1.at(i)*24*3600*5;
    Ya.at(i)=Y1.at(i)*24*3600*5;
    //    Yrunerr.at(i)=Ye.at(i)*24*3600*5;
    //    cout << setprecision(3) << "test: " << i << ", dz: " << dz.at(i) << endl;
    if(Y1.at(i)>0){
      n_p++;
    }
    Ysum += Y1.at(i);
    //    Yesum+=Ye.at(i);
    Ydasum += Yda.at(i);
    Yrunsum += Yrun.at(i);
    //    Yrunesum+=Yrunerr.at(i);
    /*
    if(Ecm.at(i)==22.2){
      cout<<i<<" energy "<<Ecm.at(i)<<" sig "<<sig.at(i)<<" n_b "<<n_b<<"  rho "<<rho*N_rate1
	  <<" range_attpc "<<range_attpc.at(0).at(i)<<" "<<range_attpc.at(1).at(i)<<" dz "<<dz.at(i)
	  <<" yeild "<<Y1.at(i)<<endl;
      cout<<endl<<endl<<endl;
      cout<<"rho dz: "<<rho*N_rate1*dz.at(i)*1e-3*1e-4<<endl;
    }
    dif_z=dz.at(i)/5;
    for(Int_t j=0; j<5; j++){
      z_tem[i]=range_attpc.at(0).at(i)+dif_z*j;
    }
    */
  }

  plot.resize(6,vector<Double_t> (0));
  //like fill for TGraph
  for(Int_t i=0; i<n; i++){
    //    if(Y1.at(i)>1e-6){
    if(Y1.at(i)>0){
      plot[0].push_back(data.at(i).first);
      plot[1].push_back(data.at(i).second);
      plot[2].push_back(dz.at(i));
      plot[3].push_back(Y1.at(i));
      plot[4].push_back(z.at(i));
      plot[5].push_back(Ya.at(i));
      //      cout<<Ecm.at(i)<<endl;
    }
  }

  if(n_p!=plot.at(0).size()){
    cout<<endl;
    cout<<"something wrong about filling data!!! finish this macro"<<endl;
    /*
    data->AutoSave();
    file->Close();
    */
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

  /*
  TGraphAsymmErrors *h_yp=
    new TGraphAsymmErrors(n_p,plot.at(0).data(),plot.at(2).data(),
			  plot.at(1).data(),plot.at(1).data(),plot.at(3).data(),plot.at(4).data());
  */

  TGraph *h_sig = new TGraph(n_p, plot.at(0).data(), plot.at(1).data());
  TGraph *h_edz = new TGraph(n_p, plot.at(0).data(), plot.at(2).data());
  TGraph *h_yei = new TGraph(n_p, plot.at(0).data(), plot.at(3).data());
  TGraph *h_Evz = new TGraph(n_p, plot.at(0).data(), plot.at(4).data());
  TGraph *h_vzy = new TGraph(n_p, plot.at(4).data(), plot.at(3).data());

  TCanvas *c1 = new TCanvas("c1","c1", 1200, 1000);
  c1->Divide(3,2);
  c1->cd(1);
  h_sig->SetMarkerStyle(20);
  h_sig->SetMarkerSize(1);
  h_sig->SetMarkerColor(kBlue);
  h_sig->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  h_sig->GetYaxis()->SetTitle("d#sigma/d#Omega");
  h_sig->SetTitle(Form("E_{c.m.} d#sigma/d#Omega in %d", ang));
  h_sig->GetXaxis()->SetLimits(0, Ecmm);
  gPad->SetLogy();
  h_sig->Draw("AP");
  c1->cd(2);
  h_edz->SetMarkerStyle(20);
  h_edz->SetMarkerSize(1);
  h_edz->SetMarkerColor(kBlue);
  h_edz->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  h_edz->GetYaxis()->SetTitle("dz [mm]");
  h_edz->SetTitle(Form("E_{c.m.} dz"));
  h_edz->GetXaxis()->SetLimits(0, Ecmm);
  gPad->SetLogy();
  h_edz->Draw("AP");
  c1->cd(3);
  h_yei->SetMarkerStyle(20);
  h_yei->SetMarkerSize(1);
  h_yei->SetMarkerColor(kBlue);
  h_yei->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  h_yei->GetYaxis()->SetTitle("Y(E) [cps]");
  h_yei->SetTitle(Form("E_{c.m.} Y(E) in %d (n_b constant, da = %d)", ang, da));
  h_yei->GetXaxis()->SetLimits(0, Ecmm);
  gPad->SetLogy();
  h_yei->Draw("AP");
  c1->cd(4);
  h_Evz->SetMarkerStyle(20);
  h_Evz->SetMarkerSize(1);
  h_Evz->SetMarkerColor(kBlue);
  h_Evz->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  h_Evz->GetYaxis()->SetTitle("Estimated z [mm]");
  h_Evz->SetTitle(Form("E_{c.m.} estimated z"));
  h_Evz->GetXaxis()->SetLimits(0, Ecmm);
  h_Evz->Draw("AP");
  c1->cd(5);
  h_vzy->SetMarkerStyle(20);
  h_vzy->SetMarkerSize(1);
  h_vzy->SetMarkerColor(kBlue);
  h_vzy->GetXaxis()->SetTitle("Estimated z [mm]");
  h_vzy->GetYaxis()->SetTitle("Y(E) [cps]");
  h_vzy->SetTitle(Form("estimated z Y(E)"));
  h_vzy->GetXaxis()->SetLimits(0, 680);
  gPad->SetLogy();
  h_vzy->Draw("AP");

  /*
  TLatex latex;
  latex.SetTextSize(0.06);
  latex.SetTextColor(kBlack);
  latex.DrawLatex(18.8,4.5e-5,"#splitline{^{12}C beam @ 50 MeV}{iso-C_{4}H_{10} 0.13 atm.}");
  c1->Update();

  //  TString outcan1="test.pdf";
  TString outcan1="./pdf_output/can_yeild.pdf";
  c1->SaveAs(outcan1);

  //save branch
  data->AutoSave();
  file->Close();
  */
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  //  cout << "output: " << outdata << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

}

std::vector<std::vector<Double_t>> cal_Ebeam_depth(){

  // set data
  vector<pair<Double_t, Double_t>> lise_data={
      {0, 60.724}, {50, 57.862}, {100, 54.873}, {150, 51.775}, {200, 48.529}, {250, 45.143}, {300, 41.562},
      {350, 37.779},{400, 33.736},{450, 29.363}, {500, 24.589}, {550, 19.089},{600, 13.089}, {650, 5.513},
      {680, 0.360}
  };

  Int_t n_data = lise_data.size();
  std::vector<std::vector<Double_t>> Ebeam(2, std::vector<Double_t>(n_data, 0));
  std::vector<std::vector<Double_t>> Vertz(2, std::vector<Double_t>(n_data, 0));
  for (Int_t i = 0; i < n_data; i++){
    Ebeam.at(0).at(i) = lise_data[i].first;
    Ebeam.at(1).at(i) = lise_data[i].second;
    Vertz.at(0).at(i) = lise_data[i].second;
    Vertz.at(1).at(i) = lise_data[i].first;
  }
  TGraph *h_Ebeam = new TGraph(n_data, Ebeam.at(0).data(), Ebeam.at(1).data());
  TGraph *h_vertz = new TGraph(n_data, Vertz.at(0).data(), Vertz.at(1).data());
  TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0, 1000);
  TF1 *f2 = new TF1("f2", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0, 65);

  h_Ebeam->Fit(f1, "QRN", "", 0, 700);
  h_vertz->Fit(f2, "QRN", "", 0, 60.8);

  TCanvas *c90 = new TCanvas("c90", "c90");
  h_Ebeam->SetMarkerStyle(20);
  h_Ebeam->SetMarkerSize(1.2);
  h_Ebeam->SetMarkerColor(kBlue);
  //   h_Ebeam->SetLineColor(kRed);
  h_Ebeam->GetXaxis()->SetTitle("Depth [mm]");
  h_Ebeam->GetYaxis()->SetTitle("beam energy [MeV]");
  h_Ebeam->SetTitle("Beam energy vs Depth");
  h_Ebeam->GetXaxis()->SetLimits(0, 700);
  h_Ebeam->Draw("AP");
  f1->SetNpx(1000);
  f1->SetLineColor(kBlack);
  f1->Draw("P same");

  TCanvas *c91 = new TCanvas("c91", "c91");
  h_vertz->SetMarkerStyle(20);
  h_vertz->SetMarkerSize(1.2);
  h_vertz->SetMarkerColor(kBlue);
  //  h_vertz->SetLineColor(kRed);
  h_vertz->GetXaxis()->SetTitle("beam energy [MeV]");
  h_vertz->GetYaxis()->SetTitle("Depth [mm]");
  h_vertz->SetTitle("Beam energy vs Depth");
  h_vertz->GetXaxis()->SetLimits(0, 65);
  h_vertz->Draw("AP");
  f2->SetNpx(1000);
  f2->SetLineColor(kBlack);
  f2->Draw("P same");

  Double_t az = f1->GetParameter(0);
  Double_t bz = f1->GetParameter(1);
  Double_t cz = f1->GetParameter(2);
  Double_t dz = f1->GetParameter(3);
  Double_t ez = f1->GetParameter(4);
  Double_t ae = f2->GetParameter(0);
  Double_t be = f2->GetParameter(1);
  Double_t ce = f2->GetParameter(2);
  Double_t de = f2->GetParameter(3);
  Double_t ee = f2->GetParameter(4);
  std::vector<std::vector<Double_t>> Ebeam_depth(2, std::vector<Double_t>(0,0));
  Ebeam_depth[0].push_back(az);
  Ebeam_depth[0].push_back(bz);
  Ebeam_depth[0].push_back(cz);
  Ebeam_depth[0].push_back(dz);
  Ebeam_depth[0].push_back(ez);
  Ebeam_depth[1].push_back(ae);
  Ebeam_depth[1].push_back(be);
  Ebeam_depth[1].push_back(ce);
  Ebeam_depth[1].push_back(de);
  Ebeam_depth[1].push_back(ee);

  std::cout << std::setprecision(4) << std::endl;
  std::cout << "Set estimation of beam energy." << std::endl;
  std::cout << "  E_{beam}(z) = a + b * z + c * z^2 + d * z^3 + e * z^4" << std::endl;
  std::cout << "  a:" << az << ", b:" << bz << ", c:" << cz << ", d:" << dz << ", e:" << ez << std::endl;
  std::cout << "Set estimation of depth." << std::endl;
  std::cout << "  z(E) = a + b * E + c * E^2 + d * E^3 + e * E^4" << std::endl;
  std::cout << "  a:" << ae << ", b:" << be << ", c:" << ce << ", d:" << de << ", e:" << ee << std::endl;

  return Ebeam_depth;
}

Double_t est_depth(std::vector<std::vector<Double_t>> &Ebeam_depth, Double_t Ecm){

  Int_t n = Ebeam_depth.at(0).size();
  Double_t Depth = -1;
  Double_t E = Ecm * 2;
  if (n != 5){
    std::cout << "Error: Ebeam_para should have 5 parameters. Current size: " << n << std::endl;
    return -1;
  }
  else {
    Double_t a = Ebeam_depth[1][0];
    Double_t b = Ebeam_depth[1][1];
    Double_t c = Ebeam_depth[1][2];
    Double_t d = Ebeam_depth[1][3];
    Double_t e = Ebeam_depth[1][4];
    Depth = a + b * E + c * pow(E, 2) + d * pow(E, 3) + e * pow(E, 4);
  }

  return Depth;
}