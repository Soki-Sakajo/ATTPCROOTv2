#define cal_Ebeam
#include <fstream>
#include "TFile.h"
#include "TObject.h"
#include "TCanvas.h"

void cal_Ebeam_para(){

   // set para
   Double_t verz = -50;

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

   TFile *outfile = new TFile("data/Ebeam_para.root", "RECREATE");

   TGraph *h_Ebeam_est = new TGraph(n_data, Ebeam.at(0).data(), Ebeam.at(1).data());
   TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0, 1000);

   h_Ebeam_est->Fit(f1, "QRN", "", 0, 700);

   TCanvas *c0 = new TCanvas("c0", "c0");
   h_Ebeam_est->SetMarkerStyle(20);
   h_Ebeam_est->SetMarkerSize(1.2);
   h_Ebeam_est->SetMarkerColor(kBlue);
   //   h_Ebeam_est->SetLineColor(kRed);
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

   h_Ebeam_est->Write("h_Ebeam_depth");
   f1->Write("f_Ebeam_depth");
   outfile->Close();

   std::cout << std::setprecision(4) << std::endl;
   std::cout << "Set estimation of beam energy." << std::endl;
   std::cout << "  function: a + b * x + c * x^2 + d * x^3 + e * x^4" << std::endl;
   std::cout << "  a:" << a << ", b:" << b << ", c:" << c << ", d:" << d << ", e:" << e << std::endl << std::endl;

#ifdef cal_Ebeam
   std::cout << "E_beam estimation at z = " << verz << " mm; E_beam = " 
   << a + b * verz + c * pow(verz, 2) + d * pow(verz, 3) + e * pow(verz, 4) << " MeV" << std::endl << std::endl;

#endif

}
