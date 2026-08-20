//#define debug_mode
#include <string>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TObject.h"
#include "TCanvas.h"

TGraph* ReadKinematics(TString kineFile);
std::vector<Double_t> cal_Ebeam_para();
Double_t est_Ebeam(std::vector<Double_t> &Ebeam_para, Double_t vertz);
TGraph* read_crosssection(TString crossFile);

void comp_kineE_beam_cross(){
   //copy from check_vd_76matm.C 2026/07/08 16:00

   //set timer
   TStopwatch timer;
   timer.Start();

   // directory
   //   Double_t vd_val = 4.05;
   Double_t vd_val = 4.07;

   Double_t shiftMeV = 3.0; // shift of beam energy in MeV

   Int_t run_start = 28;
   Int_t run_end = 112;
   std::vector<Double_t> Ebeam_para(0);

   // read file
   //   TFile * file = new TFile(Form("vd_data/vd_results_all_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val),"read");
   TFile * file = new TFile(Form("vd_data/vd_results_all_hists_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val),"read");

   // Kinematic lines.
   TGraph *kine_gsex_60_7 = (TGraph*) file->Get("kine_12c_gsex");
   TGraph *kine_exex_60_7 = ReadKinematics("./two-body_kine_files/kine_12c12c_exex_4.44_60.7.txt");

   TGraph *kine_gsgs_0   = (TGraph*) file->Get("kine_12c_gsgs_E60_7_z0");
   TGraph *kine_gsgs_50  = (TGraph*) file->Get("kine_12c_gsgs_E57_9_z50");
   TGraph *kine_gsgs_100 = (TGraph*) file->Get("kine_12c_gsgs_E55_0_z100");
   TGraph *kine_gsgs_150 = (TGraph*) file->Get("kine_12c_gsgs_E51_8_z150");
   TGraph *kine_gsgs_200 = (TGraph*) file->Get("kine_12c_gsgs_E48_5_z200");

   TGraph *kine_gsgs_300 = (TGraph*) file->Get("kine_12c_gsgs_E41_5_z300");
   TGraph *kine_gsgs_350 = (TGraph*) file->Get("kine_12c_gsgs_E37_8_z350");
   TGraph *kine_gsgs_400 = (TGraph*) file->Get("kine_12c_gsgs_E33_8_z400");
   TGraph *kine_gsgs_450 = (TGraph*) file->Get("kine_12c_gsgs_E29_5_z450");
   TGraph *kine_gsgs_500 = (TGraph*) file->Get("kine_12c_gsgs_E24_6_z500");

   TF1 *angle_gsgs_60_7    = (TF1*) file->Get("angle_12c_gsgs");
   TGraph *angle_gsex_60_7 = (TGraph*) file->Get("angle_12c_gsex");
   TGraph *angle_exex_60_7 = (TGraph*) file->Get("angle_12c_exex");

   TGraph *cross_gsgs_cm90 = read_crosssection("./sim_temp/c12c12_gsgs_theta_90_Ecm_mb_sr_10-38MeV.txt");

   Ebeam_para = cal_Ebeam_para();

   // Histogram definitions.
   // ... vertex of tracks
   TH1D *h_verz_gsgs = (TH1D*) file->Get("h_verz_gsgs");
   TH1D *h_verz_gsgs_cm50 = (TH1D*) file->Get("h_verz_gsgs_cm50");
   TH1D *h_verz_gsgs_cm60 = (TH1D*) file->Get("h_verz_gsgs_cm60");
   TH1D *h_verz_gsgs_cm70 = (TH1D*) file->Get("h_verz_gsgs_cm70");
   TH1D *h_verz_gsgs_cm80 = (TH1D*) file->Get("h_verz_gsgs_cm80");
   TH1D *h_verz_gsgs_cm90 = (TH1D*) file->Get("h_verz_gsgs_cm90");


   // ... Beam energy estimation
   TH1D *h_Ebeam = (TH1D*) file->Get("h_Ebeam");
   TH1D *h_Ebeam_gsgs = (TH1D*) file->Get("h_Ebeam_gsgs");
   TH1D *h_Ebeam_gsgs_cm50 = (TH1D*) file->Get("h_Ebeam_gsgs_cm50");
   TH1D *h_Ebeam_gsgs_cm60 = (TH1D*) file->Get("h_Ebeam_gsgs_cm60");
   TH1D *h_Ebeam_gsgs_cm70 = (TH1D*) file->Get("h_Ebeam_gsgs_cm70");
   TH1D *h_Ebeam_gsgs_cm80 = (TH1D*) file->Get("h_Ebeam_gsgs_cm80");
   TH1D *h_Ebeam_gsgs_cm90 = (TH1D*) file->Get("h_Ebeam_gsgs_cm90");

   TH1D *h_Ebcm = (TH1D*) file->Get("h_Ebcm");
   TH1D *h_Ebcm_gsgs = (TH1D*) file->Get("h_Ebcm_gsgs");
   TH1D *h_Ebcm_gsgs_cm50 = (TH1D*) file->Get("h_Ebcm_gsgs_cm50");
   TH1D *h_Ebcm_gsgs_cm60 = (TH1D*) file->Get("h_Ebcm_gsgs_cm60");
   TH1D *h_Ebcm_gsgs_cm70 = (TH1D*) file->Get("h_Ebcm_gsgs_cm70");
   TH1D *h_Ebcm_gsgs_cm80 = (TH1D*) file->Get("h_Ebcm_gsgs_cm80");
   TH1D *h_Ebcm_gsgs_cm90 = (TH1D*) file->Get("h_Ebcm_gsgs_cm90");

   TH2F *h_kineE_thetalab_gsgs_index_0 = (TH2F*) file->Get("h_kineE_thetalab_gsgs_index_0");
   TH2F *h_kineE_thetalab_gsgs_index_1 = (TH2F*) file->Get("h_kineE_thetalab_gsgs_index_1");
   TH2F *h_kineE_thetalab_gsgs_index_2 = (TH2F*) file->Get("h_kineE_thetalab_gsgs_index_2");
   TH2F *h_kineE_thetalab_gsgs_index_3 = (TH2F*) file->Get("h_kineE_thetalab_gsgs_index_3");
   TH2F *h_kineE_thetalab_gsgs_index_4 = (TH2F*) file->Get("h_kineE_thetalab_gsgs_index_4");

   //  set color
   kine_exex_60_7  -> SetLineColor(kRed);
   angle_gsgs_60_7 -> SetLineColor(kRed);
   angle_gsgs_60_7 -> SetLineWidth(1);
   angle_exex_60_7 -> SetLineColor(kRed);

   // scale hist of cross section
   Double_t bw = h_Ebcm_gsgs_cm90->GetBinWidth(1);
   Int_t shift_bin = TMath::Nint(shiftMeV / bw);

   Double_t xmin = h_Ebcm_gsgs_cm90->GetXaxis()->GetXmin();
   Double_t xmax = h_Ebcm_gsgs_cm90->GetXaxis()->GetXmax();
   TH1D *h_Ebcm_gsgs_cm90_shift =
      new TH1D("h_Ebcm_gsgs_cm90_shift", "shifted histgram;E{re}  [MeV];", h_Ebcm_gsgs_cm90->GetNbinsX() + shift_bin, xmin, xmax + shift_bin * bw);
   for (Int_t i = 1; i <= h_Ebcm_gsgs_cm90->GetNbinsX(); i++){
      h_Ebcm_gsgs_cm90_shift->SetBinContent(i + shift_bin, h_Ebcm_gsgs_cm90->GetBinContent(i));
      h_Ebcm_gsgs_cm90_shift->SetBinError(i + shift_bin, h_Ebcm_gsgs_cm90->GetBinError(i));
   }

   //   Double_t scale = 1e-1; // scale factor for cross section
   TGraph *scale_gsgs_cm90 = (TGraph*) cross_gsgs_cm90 -> Clone("scale_gsgs_cm90");
   Double_t scale = h_Ebcm_gsgs_cm90 -> GetMinimum(0.0) / TMath::MinElement(scale_gsgs_cm90->GetN(), scale_gsgs_cm90->GetY());
   Double_t gmin_ori = TMath::MinElement(scale_gsgs_cm90->GetN(), scale_gsgs_cm90->GetY()) * scale;
   Double_t gmax_ori = TMath::MaxElement(scale_gsgs_cm90->GetN(), scale_gsgs_cm90->GetY()) * scale;
   for (Int_t i = 0; i < scale_gsgs_cm90 -> GetN(); i++){
      Double_t x, y;
      scale_gsgs_cm90 -> GetPoint(i, x, y);
      scale_gsgs_cm90 -> SetPoint(i, x, y * scale);
   }
   scale_gsgs_cm90 -> SetMarkerColor(kBlack);
   scale_gsgs_cm90 -> SetMarkerStyle(20);
   scale_gsgs_cm90 -> SetMarkerSize(0.8);
   //   scale_gsgs_cm90 -> SetMinimum(gmin_ori);
   //   scale_gsgs_cm90 -> SetMaximum(gmax_ori);

   TCanvas *c1 = new TCanvas("c1", "c1");
   c1->cd();
   h_Ebcm_gsgs_cm90->SetDirectory(0);
   h_Ebcm_gsgs_cm90->SetStats(0);
   h_Ebcm_gsgs_cm90->SetLineColor(kBlue);
   h_Ebcm_gsgs_cm90->SetLineWidth(2);
   h_Ebcm_gsgs_cm90->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm90->GetXaxis()->SetRangeUser(0, 35);
   h_Ebcm_gsgs_cm90->SetMinimum(0.5);
   h_Ebcm_gsgs_cm90->SetMaximum(50);
   h_Ebcm_gsgs_cm90->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm90->Draw("HIST");
   scale_gsgs_cm90->SetMinimum(5);
   scale_gsgs_cm90->SetMaximum(50);
   scale_gsgs_cm90->Draw("PL same");
   gPad->Update();
   // set second y-axis for cross section
   Double_t xr = gPad->GetUxmax();
   Double_t ymin = h_Ebcm_gsgs_cm90->GetMinimum();
   Double_t ymax = h_Ebcm_gsgs_cm90->GetMaximum();
   TGaxis *r_axis = new TGaxis(xr, ymin, xr, ymax, 5, 500, 510, "+LG");
   r_axis->SetTitle("cross section [mb/sr]");
   r_axis->SetLineColor(kBlack);
   r_axis->SetLabelColor(kBlack);
   r_axis->SetTitleColor(kBlack);
   r_axis->Draw();
   c1->Update();

   TCanvas *c2 = new TCanvas("c2", "c2");
   c2->cd();
   h_Ebcm_gsgs_cm90_shift->SetDirectory(0);
   h_Ebcm_gsgs_cm90_shift->SetStats(0);
   h_Ebcm_gsgs_cm90_shift->SetLineColor(kBlue);
   h_Ebcm_gsgs_cm90_shift->SetLineWidth(2);
   h_Ebcm_gsgs_cm90_shift->GetXaxis()->SetTitle("shifted E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm90_shift->GetXaxis()->SetRangeUser(0, 35);
   h_Ebcm_gsgs_cm90_shift->SetMinimum(0.5);
   h_Ebcm_gsgs_cm90_shift->SetMaximum(50);
   h_Ebcm_gsgs_cm90_shift->SetTitle(Form("%.1f MeV Shifted E_{beam.cm} (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)", (Double_t)shiftMeV));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm90_shift->Draw("HIST");
   scale_gsgs_cm90->SetMinimum(5);
   scale_gsgs_cm90->SetMaximum(50);
   scale_gsgs_cm90->Draw("PL same");
   gPad->Update();
   // set second y-axis for cross section
   Double_t xr2 = gPad->GetUxmax();
   Double_t ymin2 = h_Ebcm_gsgs_cm90_shift->GetMinimum();
   Double_t ymax2 = h_Ebcm_gsgs_cm90_shift->GetMaximum();
   TGaxis *r_axis2 = new TGaxis(xr2, ymin2, xr2, ymax2, 5, 500, 510, "+LG");
   r_axis2->SetTitle("cross section [mb/sr]");
   r_axis2->SetLineColor(kBlack);
   r_axis2->SetLabelColor(kBlack);
   r_axis2->SetTitleColor(kBlack);
   r_axis2->Draw();
   c2->Update();

   TCanvas *c100 = new TCanvas("c100","c100");
   c100->Divide(2,2);
   c100->cd(1);
   h_kineE_thetalab_gsgs_index_0->SetStats(0);
   gPad->SetLogz();
   h_kineE_thetalab_gsgs_index_0->Draw("colz");
   kine_gsgs_0->Draw("PL same");
   //  kine_gsgs_50->Draw("same");
   kine_gsgs_100->Draw("PL same");

   c100->cd(2);
   h_kineE_thetalab_gsgs_index_1->SetStats(0);
   gPad->SetLogz();
   h_kineE_thetalab_gsgs_index_1->Draw("colz");
   kine_gsgs_100->Draw("PL same");
   //  kine_gsgs_150->Draw("same");
   kine_gsgs_200->Draw("PL same");

   c100->cd(3);
   h_kineE_thetalab_gsgs_index_3->SetStats(0);
   gPad->SetLogz();
   h_kineE_thetalab_gsgs_index_3->Draw("colz");
   kine_gsgs_300->Draw("PL same");
   //  kine_gsgs_350->Draw("same");
   kine_gsgs_400->Draw("PL same");

   c100->cd(4);
   h_kineE_thetalab_gsgs_index_4->SetStats(0);
   gPad->SetLogz();
   h_kineE_thetalab_gsgs_index_4->Draw("colz");
   kine_gsgs_400->Draw("PL same");
   //  kine_gsgs_450->Draw("same");
   kine_gsgs_500->Draw("PL same");

   // stop timer
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl;
   cout << "Real time: " << rtime << " s, CPU time: " << ctime << " s" << endl << endl;

}

TGraph* ReadKinematics(TString kineFile){

   std::vector<Double_t> ThetacmS, ThetalabR, ElabR, ThetalabS, ElabS, MomlabR;
#ifdef debug_mode
   // Debug: Check if file exists and current working directory
   std::cout << "DEBUG ReadKinematics: CWD = " << gSystem->pwd() << std::endl;
   std::cout << "DEBUG ReadKinematics: Attempting to read: " << kineFile << std::endl;
   if (gSystem->AccessPathName(kineFile.Data(), kFileExists)) {
      std::cout << "DEBUG ReadKinematics: FILE NOT FOUND!" << std::endl;
   } else {
      std::cout << "DEBUG ReadKinematics: File found" << std::endl;
   }
#endif
   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr -> is_open() || !kineStr -> good()){
      std::cout << "Warning : No kinematics file found or cannot read: " << kineFile << std::endl;
      return new TGraph();
   }
   double tcm, tlabr, elabr, tlabs, elabs;
   string line;
   while (getline(*kineStr, line)){
      if(line.empty() || line[0] == '#'){
         continue;
      }
      istringstream iss(line);
      if (iss >> tlabr >> elabr){
         ThetalabR.push_back(tlabr);
         ElabR.push_back(elabr);
      }
   }
#ifdef debug_mode
      std::cout << "DEBUG ReadKinematics: Successfully read " << ThetalabR.size() << " points" << std::endl;
#endif
   if (ThetalabR.size() == 0){
      std::cout << " Warning : No data read from kinematics file: " << kineFile << std::endl;
   }
   //   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   TGraph *kine = new TGraph(ThetalabR.size(), ThetalabR.data(), ElabR.data());
   //   kine->Sort();
   return kine;
}

std::vector<Double_t> cal_Ebeam_para(){

   // set data
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 60.724}, {50, 57.862}, {100, 54.873}, {150, 51.775}, {200, 48.529}, {250, 45.143}, {300, 41.562},
         {350, 37.779},{400, 33.736},{450, 29.363}, {500, 24.589}, {550, 19.089},{600, 13.089}, {650, 5.513},
         {680, 0.360}
   };

   Int_t n_data = lise_data.size();
   std::vector<std::vector<Double_t>> Ebeam(2, std::vector<Double_t>(n_data, 0));
   for (Int_t i = 0; i < n_data; i++){
      Ebeam.at(0).at(i) = lise_data[i].first;
      Ebeam.at(1).at(i) = lise_data[i].second;
   }
   TGraph *h_Ebeam_est = new TGraph(n_data, Ebeam.at(0).data(), Ebeam.at(1).data());
   TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0, 1000);

   h_Ebeam_est->Fit(f1, "QRN", "", 0, 700);

   TCanvas *c0 = new TCanvas("c0", "c0");
   h_Ebeam_est->SetMarkerStyle(20);
   h_Ebeam_est->SetMarkerSize(1.2);
   h_Ebeam_est->SetMarkerColor(kBlue);
   h_Ebeam_est->SetLineColor(kRed);
   h_Ebeam_est->GetXaxis()->SetTitle("Depth [mm]");
   h_Ebeam_est->GetYaxis()->SetTitle("beam energy [MeV]");
   h_Ebeam_est->SetTitle("Beam energy vs Depth");
   h_Ebeam_est->GetXaxis()->SetLimits(0, 700);
   h_Ebeam_est->Draw("AP");
   f1->SetNpx(1000);
   f1->SetLineColor(kBlack);
   f1->Draw("P same");

   Double_t a = f1->GetParameter(0);
   Double_t b = f1->GetParameter(1);
   Double_t c = f1->GetParameter(2);
   Double_t d = f1->GetParameter(3);
   Double_t e = f1->GetParameter(4);
   std::vector<Double_t> Ebeam_para;
   Ebeam_para.push_back(a);
   Ebeam_para.push_back(b);
   Ebeam_para.push_back(c);
   Ebeam_para.push_back(d);
   Ebeam_para.push_back(e);

   std::cout << std::setprecision(4) << std::endl;
   std::cout << "Set estimation of beam energy." << std::endl;
   std::cout << "  function: a + b * x + c * x^2 + d * x^3 + e * x^4" << std::endl;
   std::cout << "  a:" << a << ", b:" << b << ", c:" << c << ", d:" << d << ", e:" << e << std::endl;

   return Ebeam_para;
}

Double_t est_Ebeam(std::vector<Double_t> &Ebeam_para, Double_t vertz){

   Double_t Ebeam = -1;
   if (Ebeam_para.size() != 5){
      std::cout << "Error: Ebeam_para should have 5 parameters. Current size: " << Ebeam_para.size() << std::endl;
      return -1;
   }
   else {
      Double_t a = Ebeam_para[0];
      Double_t b = Ebeam_para[1];
      Double_t c = Ebeam_para[2];
      Double_t d = Ebeam_para[3];
      Double_t e = Ebeam_para[4];
      Ebeam = a + b * vertz + c * pow(vertz, 2) + d * pow(vertz, 3) + e * pow(vertz, 4);
   }

   return Ebeam;
}

TGraph* read_crosssection(TString crossFile){

   Int_t n_cs = 0;
   Double_t E, cs;
   std::vector<std::vector<Double_t>> plot(0, std::vector<Double_t>(0,0));
   std::vector<std::pair<Double_t, Double_t>> data;

   #ifdef debug_mode
   // Debug: Check if file exists and current working directory
   std::cout << "DEBUG ReadKinematics: CWD = " << gSystem->pwd() << std::endl;
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
   std::cout << "check vector size of Ecm : " << n_cs << std::endl;

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
