//#define debug_mode
//#define catima_check
//#define nom_check
//#define peak_check
//#define c12_check
//#define vertex_check
//#define vertex_index
//#define vertex_depth
//#define gsgs_check
//#define kine_comp
//#define check_verz0
#include <string>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TString.h"
#include "TObject.h"
#include "TCanvas.h"
#include "Math/Point3D.h"

TGraph* ReadKinematics(TString kineFile, TString ReturnType);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);
void draw_ind(TString cname, TString states, Int_t n_group, Int_t n_h_z, TH1D* h_ver, TH2F* h_verxy, TH2F* h_E_theta, TH2F* h_theta_theta, 
         std::vector<TH1D*> &h_verz_i, std::vector<TH2F*> &h_verxy_i, std::vector<TH2F*> &h_E_theta_i, std::vector<TH2F*> &h_theta_theta_i,
         TGraph *kine_gsgs, TGraph *kine_gsex, TGraph *kine_exex, TGraph *ang_gsgs, TGraph *ang_gsex, TGraph *ang_exex);
void draw_dep(TString cname, TString states, TString LineType, Int_t n_group, Int_t n_div,
         std::vector<TH2F*> &h_E_theta_i, std::vector<TGraph*> &kine_i);

std::vector<Double_t> cal_Ebeam_para();
Double_t est_Ebeam(std::vector<Double_t> &Ebeam_para, Double_t vertz);
TGraph* read_crosssection(TString crossFile);

void c12_12c12c_ana_76matm_rough(){
   //copy from kine.C 2026/05/06 12:20

   //set timer
   TStopwatch timer;
   timer.Start();

   // set parameters
   Double_t del_phi = 10; // cut value; phi1 - phi2 - 180 deg < del_phi

   // directory (drift velocity)
   Double_t vd_val = 4.07;

   // files.
   //   std::vector runNums = {52};
   //   std::vector runNums = {50,51,52,53,54,55,56,57,58};

   std::vector runNums = {
      28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,50,
      51,52,53,54,55,56,57,58,62,63,64,66,67,68,69,70,71,75,76,77,
      78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,95,96,97,98,99,
      100,101,102,103,104,105,106,107,108,109,110,111,112
   };
   /*
   */

   const bool vertex_corr = true; // true: use the values calculated from vertex z; false: use rough values
   const bool kine_comp_b = false;
   const Int_t n_group = 9; // group number of vertex z: 0-100, ... 700-800, 800-1000.
   //   const Int_t verz_h = 600;
   const Int_t verz_h = (n_group - 1) * 100;
   const Int_t n_div = 3; // number of divided canvas for each line
   const Int_t n_dep = 17; // number of hists by vertex z

   Int_t run_start = runNums.front();
   Int_t run_end = runNums.back();
   Double_t r = 0;
   std::vector<Double_t> Ebeam_para(0);
#ifdef kine_comp
   kine_comp_b = true;
   //   TString f_Re = TString::Format("data2/c12_12c12c/ana_results_all_hists_kine_comp_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val);
   TString f_Re = TString::Format("data2/c12_12c12c/ana_results_kine_comp_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val);
#else
   kine_comp_b = false;
   //   TString f_Re = TString::Format("data2/c12_12c12c/ana_results_all_hists_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val);
   TString f_Re = TString::Format("data2/c12_12c12c/ana_results_run%d-run%d_vd%.2f.root", run_start, run_end, vd_val);
#endif

   TFile * Results = new TFile(f_Re,"recreate");

   FairRunAna *run = new FairRunAna(); // Forcing a dummy run
#ifdef debug_mode
   // Debug: Check working directory and file existence
   std::cout << "DEBUG: Current working directory: " << gSystem->pwd() << std::endl;
   std::cout << "DEBUG: Checking kinematic files..." << std::endl;
   std::string testFile = "./two-body_kine_files/kine_12c12c_gsex_4.44_60.7.txt";
   if (gSystem->AccessPathName(testFile.c_str(), kFileExists)) {
      std::cout << "ERROR: File NOT found: " << testFile << std::endl;
   } else {
      std::cout << "OK: File found: " << testFile << std::endl;
   }
   kine_comp_b = true;
#endif

   // AtMap to check if a hit belong to a big pad or small pad.
   TString scriptfile = "rcnp_map_size.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   AtTpcMap *map = new AtTpcMap();
   map->ParseXMLMap(mapDir.Data());
   map->GeneratePadPlane();

   // Punch through filter.
   double punchThroughThreshold = 20;
   AtTools::AtPunchThroughChecker punchThroughChecker = AtTools::AtPunchThroughChecker();
   punchThroughChecker.SetDistanceThreshold(punchThroughThreshold);

   // ELoss model for kinetic energy estimations.
   // ... ATTPC
   double density = 1.7884e-4; // 76 matm
   //   double density = 3.0406e-4; // 129 matm
   std::vector<std::tuple<int, int, int>> materialComponents;
   materialComponents.push_back(std::make_tuple(12, 6, 4));
   materialComponents.push_back(std::make_tuple(1, 1, 10));

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC4H10_p = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C4H10_p");
   eLossModelC4H10_p->SetMaterial(materialComponents);
   eLossModelC4H10_p->SetProjectile(1, 1, 1.007825031898);
   //   eLossModels.push_back(std::move(eLossModelC4H10_p));
   //   eLossModelC4H10_p->SetPDGCode("1000010010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC4H10_alpha = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C4H10_alpha");
   eLossModelC4H10_alpha->SetMaterial(materialComponents);
   eLossModelC4H10_alpha->SetProjectile(4, 2, 4.00260325413);
   //   eLossModels.push_back(std::move(eLossModelC4H10_alpha));
   //   eLossModelC4H10_alpha->SetPDGCode("1000010010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC4H10_12C = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C4H10_12C");
   eLossModelC4H10_12C->SetMaterial(materialComponents);
   eLossModelC4H10_12C->SetProjectile(12, 6, 12);
   //   eLossModels.push_back(std::move(eLossModelC4H10_12C));
   //   eLossModelC4H10_12C->SetPDGCode("1000120060");

   //For warming up CATIMA calculations.
   const Double_t warmupRangeP = eLossModelC4H10_p->GetRange(0.1);
   const Double_t warmupRangeAlpha = eLossModelC4H10_alpha->GetRange(0.1);
   const Double_t warmupRange12C = eLossModelC4H10_12C->GetRange(0.1);

   if (!std::isfinite(warmupRangeP) || !std::isfinite(warmupRangeAlpha) || !std::isfinite(warmupRange12C)) {
      throw std::runtime_error( "Failed to initialize CATIMA energy-loss tables");
   }

#ifdef catima_check
   std::cout << "Catima check... After definition. " << std::flush;
   r = eLossModelC4H10_12C->GetRange(0.1);
   std::cout << " range = " << r << ", density = " << eLossModelC4H10_12C->GetDensity()
            << ", A = " << eLossModelC4H10_12C->GetAtomicMassNumber() << ", Z = " << eLossModelC4H10_12C->GetChargeNumber()
            << ", mass = " << eLossModelC4H10_12C->GetMassAmu() << ", dEdx = " << eLossModelC4H10_12C->GetdEdx(0.1) << std::endl;
   /*
   std::cout << std::endl << std::endl;
   std::cout << "DEBUG: Checking catima calculations for C4H10_12C." << std::endl;
   for (double e : {0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0}) {
      const double r = eLossModelC4H10_12C->GetRange(e);
      std::cout << "E = " << e << " range = " << r << " finite = " << std::isfinite(r) << std::endl;
   }
   */
#endif

   // Cut files.
   gROOT->ProcessLine(".x ./cut_files/charge_range_12c.C");
   TCutG *cut12c = (TCutG*) gROOT->FindObject("charge_range_12c");

   gROOT->ProcessLine(".x ./cut_files/charge_range_alpha.C");
   TCutG *cutalpha = (TCutG*) gROOT->FindObject("charge_range_alpha");

   gROOT->ProcessLine(".x ./cut_files/charge_range_p.C");
   TCutG *cutproton = (TCutG*) gROOT->FindObject("charge_range_p");

   gROOT->ProcessLine(".x ./cut_files/range_theta_ext.C");
   TCutG *cutext = (TCutG*) gROOT->FindObject("range_theta_ext");

   gROOT->ProcessLine(".x ./cut_files/charge_range_mult_1.C");
   TCutG *cutmul1 = (TCutG*) gROOT->FindObject("charge_range_mult_1");

   gROOT->ProcessLine(".x ./cut_files/charge_range_mult_2.C");
   TCutG *cutmul2 = (TCutG*) gROOT->FindObject("charge_range_mult_2");

   gROOT->ProcessLine(".x ./cut_files/theta_theta_12c_gs-gs.C");
   TCutG *theta_gsgs = (TCutG*) gROOT->FindObject("theta_theta_gs-gs");

   gROOT->ProcessLine(".x ./cut_files/theta_theta_12c_gs-ex.C");
   TCutG *theta_gsex = (TCutG*) gROOT->FindObject("theta_theta_gs-ex");

   gROOT->ProcessLine(".x ./cut_files/theta_theta_12c_ex-ex.C");
   TCutG *theta_exex = (TCutG*) gROOT->FindObject("theta_theta_ex-ex");

   // Kinematic lines.
   TGraph *kine_gsgs_0   = ReadKinematics("gsgs_69.5_verz_0",   "kine");
   TGraph *kine_gsgs_50  = ReadKinematics("gsgs_66.9_verz_50",  "kine");
   TGraph *kine_gsgs_100 = ReadKinematics("gsgs_64.2_verz_100", "kine");
   TGraph *kine_gsgs_150 = ReadKinematics("gsgs_61.5_verz_150", "kine");
   TGraph *kine_gsgs_200 = ReadKinematics("gsgs_58.6_verz_200", "kine");
   TGraph *kine_gsgs_250 = ReadKinematics("gsgs_55.7_verz_250", "kine");
   TGraph *kine_gsgs_300 = ReadKinematics("gsgs_52.6_verz_300", "kine");
   TGraph *kine_gsgs_350 = ReadKinematics("gsgs_49.4_verz_350", "kine");
   TGraph *kine_gsgs_400 = ReadKinematics("gsgs_46.0_verz_400", "kine");
   TGraph *kine_gsgs_450 = ReadKinematics("gsgs_42.5_verz_450", "kine");
   TGraph *kine_gsgs_500 = ReadKinematics("gsgs_38.8_verz_500", "kine");
   TGraph *kine_gsgs_550 = ReadKinematics("gsgs_34.8_verz_550", "kine");
   TGraph *kine_gsgs_600 = ReadKinematics("gsgs_30.5_verz_600", "kine");
   TGraph *kine_gsgs_650 = ReadKinematics("gsgs_25.9_verz_650", "kine");
   TGraph *kine_gsgs_700 = ReadKinematics("gsgs_20.7_verz_700", "kine");
   TGraph *kine_gsgs_750 = ReadKinematics("gsgs_14.8_verz_750", "kine");
   TGraph *kine_gsgs_800 = ReadKinematics( "gsgs_7.7_verz_800", "kine");

   TGraph *kine_gsex_0   = ReadKinematics("gsex_69.5_verz_0",   "kine");
   TGraph *kine_gsex_50  = ReadKinematics("gsex_66.9_verz_50",  "kine");
   TGraph *kine_gsex_100 = ReadKinematics("gsex_64.2_verz_100", "kine");
   TGraph *kine_gsex_150 = ReadKinematics("gsex_61.5_verz_150", "kine");
   TGraph *kine_gsex_200 = ReadKinematics("gsex_58.6_verz_200", "kine");
   TGraph *kine_gsex_250 = ReadKinematics("gsex_55.7_verz_250", "kine");
   TGraph *kine_gsex_300 = ReadKinematics("gsex_52.6_verz_300", "kine");
   TGraph *kine_gsex_350 = ReadKinematics("gsex_49.4_verz_350", "kine");
   TGraph *kine_gsex_400 = ReadKinematics("gsex_46.0_verz_400", "kine");
   TGraph *kine_gsex_450 = ReadKinematics("gsex_42.5_verz_450", "kine");
   TGraph *kine_gsex_500 = ReadKinematics("gsex_38.8_verz_500", "kine");
   TGraph *kine_gsex_550 = ReadKinematics("gsex_34.8_verz_550", "kine");
   TGraph *kine_gsex_600 = ReadKinematics("gsex_30.5_verz_600", "kine");
   TGraph *kine_gsex_650 = ReadKinematics("gsex_25.9_verz_650", "kine");
   TGraph *kine_gsex_700 = ReadKinematics("gsex_20.7_verz_700", "kine");
   TGraph *kine_gsex_750 = ReadKinematics("gsex_14.8_verz_750", "kine");

   TGraph *kine_exex_0   = ReadKinematics("exex_69.5_verz_0",   "kine");
   TGraph *kine_exex_50  = ReadKinematics("exex_66.9_verz_50",  "kine");
   TGraph *kine_exex_100 = ReadKinematics("exex_64.2_verz_100", "kine");
   TGraph *kine_exex_150 = ReadKinematics("exex_61.5_verz_150", "kine");
   TGraph *kine_exex_200 = ReadKinematics("exex_58.6_verz_200", "kine");
   TGraph *kine_exex_250 = ReadKinematics("exex_55.7_verz_250", "kine");
   TGraph *kine_exex_300 = ReadKinematics("exex_52.6_verz_300", "kine");
   TGraph *kine_exex_350 = ReadKinematics("exex_49.4_verz_350", "kine");
   TGraph *kine_exex_400 = ReadKinematics("exex_46.0_verz_400", "kine");
   TGraph *kine_exex_450 = ReadKinematics("exex_42.5_verz_450", "kine");
   TGraph *kine_exex_500 = ReadKinematics("exex_38.8_verz_500", "kine");
   TGraph *kine_exex_550 = ReadKinematics("exex_34.8_verz_550", "kine");
   TGraph *kine_exex_600 = ReadKinematics("exex_30.5_verz_600", "kine");
   TGraph *kine_exex_650 = ReadKinematics("exex_25.9_verz_650", "kine");
   TGraph *kine_exex_700 = ReadKinematics("exex_20.7_verz_700", "kine");

   TGraph *angle_gsgs_0 = ReadKinematics("gsgs_69.5_verz_0", "angle");
   TGraph *angle_gsex_0 = ReadKinematics("gsex_69.5_verz_0", "angle");
   TGraph *angle_exex_0 = ReadKinematics("exex_69.5_verz_0", "angle");

   TF1 *sum_kine_beam = new TF1("sum_kine_beam", "x", 0, 70);
   TF1 *angle_xy = new TF1("angle_xy", "x", -200, 200);

   TGraph *cross_gsgs_cm90 = read_crosssection("./sim_temp/c12c12_gsgs_theta_90_Ecm_mb_sr_10-38MeV.txt");

   std::vector<TGraph*> kine_gsgs_i = {
      kine_gsgs_0, kine_gsgs_50, kine_gsgs_100, kine_gsgs_150, kine_gsgs_200, kine_gsgs_250,
      kine_gsgs_300, kine_gsgs_350, kine_gsgs_400, kine_gsgs_450, kine_gsgs_500, kine_gsgs_550,
      kine_gsgs_600, kine_gsgs_650, kine_gsgs_700, kine_gsgs_750, kine_gsgs_800
   };
   std::vector<TGraph*> kine_gsex_i = {
      kine_gsex_0, kine_gsex_50, kine_gsex_100, kine_gsex_150, kine_gsex_200, kine_gsex_250,
      kine_gsex_300, kine_gsex_350, kine_gsex_400, kine_gsex_450, kine_gsex_500, kine_gsex_550,
      kine_gsex_600, kine_gsex_650, kine_gsex_700, kine_gsex_750
   };
   std::vector<TGraph*> kine_exex_i = {
      kine_exex_0, kine_exex_50, kine_exex_100, kine_exex_150, kine_exex_200, kine_exex_250,
      kine_exex_300, kine_exex_350, kine_exex_400, kine_exex_450, kine_exex_500, kine_exex_550,
      kine_exex_600, kine_exex_650, kine_exex_700
   };

   if (vertex_corr){
      std::cout << std::endl << std::endl;
      std::cout << "Vertex correction is ON." << std::endl << std::endl;
   }
   else if (!vertex_corr){
      std::cout << std::endl << std::endl;
      std::cout << "Vertex correction is OFF." << std::endl << std::endl;
   }
   else{
      std::cout << std::endl << std::endl;
      std::cerr << "Please set the correct vertex_corr: true or false" << std::endl << std::endl;
      gSystem->Exit(1);
   }

   Ebeam_para = cal_Ebeam_para();
   TF1 *f_Ebeam = new TF1("f_Ebeam", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[5]*x^5", 0, 1000);
   f_Ebeam->SetParameters(Ebeam_para[0], Ebeam_para[1], Ebeam_para[2], Ebeam_para[3], Ebeam_para[4], Ebeam_para[5]);

   // Characteristic definitions
   bool tracks_vertex = false;
   bool alpha_tracks = false;
   bool proton_tracks = false;
   const Int_t n_h_z = n_group +1;
   Int_t nruns = 1;
   Int_t nevents = 0;
   Int_t narray = 10;
   Int_t ntrack = 0;
   Int_t itrack = 0;
   Int_t nalpha = 0;
   Int_t nproton = 0;
   Int_t nbragg = 0;
   Int_t nbrain = 0;
   Int_t nbrano = 0;
   Int_t vindex = 0;
   Int_t n_bragg_true = 0;
   Int_t n_bragg_false = 0;
   Int_t n_12c12c = 0;
   Int_t n_no_vertex = 0;
   Int_t n_ver_z_01 = 0;
   Double_t rad = -100;
   Double_t r_tem = 0;
   Double_t r_max = 0;
   Double_t r_tri = 300;
   Double_t r_tri_12c = 300;
   Double_t max_r_max = 0;
   Double_t max_r_max_12c = 0;
   Double_t track_lastx = 0;
   Double_t track_lasty = 0;
   Double_t track_lastz = 0;
   Double_t vtx = 0;
   Double_t vty = 0;
   Double_t vtz = 0;
   Double_t Ebeam = 0;
   Double_t Ebeam_cm = 0;
   Double_t track_theta[narray];
   Double_t track_phi[narray];
   Double_t track_range[narray];
   Double_t track_charge[narray];
   Double_t track_r[narray];
   Double_t track_dedx[narray];
   Double_t track_KinE[narray];
   Double_t vertx[narray];
   Double_t verty[narray];
   Double_t vertz[narray];
   Double_t vertex_theta[narray];
   Double_t vertex_phi[narray];
   Double_t vertex_range[narray];
   Double_t vertex_KinE[narray];
   ROOT::Math::XYZPoint vertexPoint(0, 0, 0);
   ROOT::Math::XYZPoint vertexPoint2(0, 0, 0);
   std::vector<Int_t> track6(0);
   std::vector<Int_t> peak1(0);
   std::vector<Int_t> peak2(0);
   std::vector<Int_t> peak3(0);
   std::vector<Double_t> v_distance(0);
   std::vector<ROOT::Math::XYZPoint> v_lastPoint{};
   std::vector<std::vector<Int_t>> track2_ver0(2, std::vector<Int_t>(0));
   std::vector<std::vector<Int_t>> vertex_z_01(2, std::vector<Int_t>(0));

   // Histogram definitions.
   // ... TH1 hist for checking something.
   TH1D *h_ntra = new TH1D("h_ntra", "h_ntra;NTracks", 11, -0.5, 10.5);
   TH1D *h_rmax = new TH1D("h_rmax", "h_rmax;Rmax [mm]", 150, 0, 300);
   TH1D *h_rmax_12c = new TH1D("h_rmax_12c", "h_rmax_12c;Rmax [mm]", 100, 0, 300);
   //   TH1D *h_sum_theta_cut12c_ela = new TH1D("h_sum_theta_cut12c_ela", "h_sum_theta_cut12c_ela;theta", 100, 40, 140);
   TH1D *h_sum_theta_cut12c_ela = new TH1D("h_sum_theta_cut12c_ela", "h_sum_theta_cut12c_ela;theta", 200, 40, 140);
   TH1D *h_sum_theta_cut12c_run52 = new TH1D("h_sum_theta_cut12c_run52", "h_sum_theta_cut12c_run52;theta", 200, 40, 140);
   TH1D *h_sum_theta_gsgs = new TH1D("h_sum_theta_gsgs", "h_sum_theta_gsgs;theta", 60, 70, 100);
   TH1D *h_sum_theta_gsex = new TH1D("h_sum_theta_gsex", "h_sum_theta_gsex;theta", 60, 70, 100);
   TH1D *h_sum_theta_exex = new TH1D("h_sum_theta_exex", "h_sum_theta_exex;theta", 60, 70, 100);

   // ... ATTPC PID
   TH2F *h_charge_range = new TH2F("h_charge_range", "h_charge_range;roughRange [mm];Charge [ADC]", 600, 0, 1200, 600, 0, 6e5);
   TH2F *h_charge_range_cutphi = new TH2F("h_charge_range_cutphi", "h_charge_range_cutphi;roughRange [mm];Charge [ADC]", 300, 0, 1200, 300, 0, 6e5);
   TH2F *h_charge_range_cutphi_2tra = new TH2F("h_charge_range_cutphi_2tra", "h_charge_range_cutphi_2tra;roughRange [mm];Charge [ADC]", 300, 0, 1200, 300, 0, 6e5);
   TH2F *h_charge_range_cut12c_ela = new TH2F("h_charge_range_cut12c_ela", "h_charge_range_cut12c_ela;roughRange [mm];Charge [ADC]", 300, 0, 1200, 300, 0, 6e5);
   
   // ... .. using dEdx
   TH2F *h_dEdx_range = new TH2F("h_dEdx_range", "h_dEdx_range;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_backwards = new TH2F("h_dEdx_range_backwards", "h_dEdx_range_backwards;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_cutphi = new TH2F("h_dEdx_range_cutphi", "h_dEdx_range_cutphi;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_cutphi_12c = new TH2F("h_dEdx_range_cutphi_12c", "h_dEdx_range_cutphi_12c;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);

   // ... kinematics 
   TH2F *h_kineE_thetalab = new TH2F("h_kineE_thetalab", "h_kineE_thetalab;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 600, 0, 60);

   TH2F *h_kineE_thetalab_carbon = new TH2F("h_kineE_thetalab_carbon", "h_kineE_thetalab_carbon;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 600, 0, 60);
   TH2F *h_kineE_thetalab_gsgs = new TH2F("h_kineE_thetalab_gsgs", "h_kineE_thetalab_gsgs;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 60, 0, 60);
   TH2F *h_kineE_thetalab_gsex = new TH2F("h_kineE_thetalab_gsex", "h_kineE_thetalab_gsex;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 60, 0, 60);
   TH2F *h_kineE_thetalab_exex = new TH2F("h_kineE_thetalab_exex", "h_kineE_thetalab_exex;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 60, 0, 60);
   /*
   TH2F *h_kineE_thetalab_carbon = new TH2F("h_kineE_thetalab_carbon", "h_kineE_thetalab_carbon;#theta_{LAB} [deg];roughKinE [MeV]", 50, 0, 100, 300, 0, 60);
   TH2F *h_kineE_thetalab_gsgs = new TH2F("h_kineE_thetalab_gsgs", "h_kineE_thetalab_gsgs;#theta_{LAB} [deg];roughKinE [MeV]", 50, 0, 100, 30, 0, 60);
   TH2F *h_kineE_thetalab_gsex = new TH2F("h_kineE_thetalab_gsex", "h_kineE_thetalab_gsex;#theta_{LAB} [deg];roughKinE [MeV]", 50, 0, 100, 30, 0, 60);
   TH2F *h_kineE_thetalab_exex = new TH2F("h_kineE_thetalab_exex", "h_kineE_thetalab_exex;#theta_{LAB} [deg];roughKinE [MeV]", 50, 0, 100, 30, 0, 60);
   */
   TH2F *h_sumkine_kineE_gsgs = new TH2F("h_sumkine_kineE_gsgs", "h_sumkine_kineE_gsgs;E_{kine} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);
   TH2F *h_sumkine_kineE_ver_gsgs = new TH2F("h_sumkine_kineE_ver_gsgs", "h_sumkine_kineE_gsgs;E_{kine} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);
   TH2F *h_sumkine_kineE_gsgs_cm90 = new TH2F("h_sumkine_kineE_gsgs_cm90", "h_sumkine_kineE_gsgs_cm90;E_{kine} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);
   TH2F *h_sumkine_verz_cut12c_ela = new TH2F("h_sumkine_verz_cut12c_ela", "h_sumkine_verz_gsgs;E_{beam} [MeV];sum_kineE [MeV]", 121, -10, 1200, 160, 0, 80);
   TH2F *h_sumkine_verz_gsgs = new TH2F("h_sumkine_verz_gsgs", "h_sumkine_verz_gsgs;E_{beam} [MeV];sum_kineE [MeV]", 121, -10, 1200, 160, 0, 80);

   // ... angle correlations
   TH2F *h_range_thetalab = new TH2F("h_range_thetalab", "h_range_thetalab", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutphi = new TH2F("h_range_thetalab_cutphi", "h_range_thetalab_cutphi", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutphi_2tra = new TH2F("h_range_thetalab_cutphi_2tra", "h_range_thetalab_cutphi_2tra", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cut12c_ela = new TH2F("h_range_thetalab_cut12c_ela", "h_range_thetalab_cut12c_ela", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_gsgs = new TH2F("h_range_thetalab_gsgs", "h_range_thetalab_gsgs", 50, 0, 100, 300, 0, 600);
   TH2F *h_range_thetalab_gsex = new TH2F("h_range_thetalab_gsex", "h_range_thetalab_gsex", 50, 0, 100, 300, 0, 600);
   TH2F *h_range_thetalab_exex = new TH2F("h_range_thetalab_exex", "h_range_thetalab_exex", 50, 0, 100, 300, 0, 600);
   // ... .. theta vs theta
   TH2F *h_thetalab_thetalab_cutphi = new TH2F("h_thetalab_thetalab_cutphi", "h_thetalab_thetalab_cutphi", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_cutphi_2tra = new TH2F("h_thetalab_thetalab_cutphi_2tra", "h_thetalab_thetalab_cutphi_2tra", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_cut12c_ela = new TH2F("h_thetalab_thetalab_cut12c_ela", "h_thetalab_thetalab_cut12c_ela", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs = new TH2F("h_thetalab_thetalab_gsgs", "h_thetalab_thetalab_gsgs", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_ver_gsgs = new TH2F("h_thetalab_thetalab_ver_gsgs", "h_thetalab_thetalab_ver_gsgs", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs_cm50 = new TH2F("h_thetalab_thetalab_gsgs_cm50", "h_thetalab_thetalab_gsgs_cm50", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs_cm60 = new TH2F("h_thetalab_thetalab_gsgs_cm60", "h_thetalab_thetalab_gsgs_cm60", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs_cm70 = new TH2F("h_thetalab_thetalab_gsgs_cm70", "h_thetalab_thetalab_gsgs_cm70", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs_cm80 = new TH2F("h_thetalab_thetalab_gsgs_cm80", "h_thetalab_thetalab_gsgs_cm80", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_gsgs_cm90 = new TH2F("h_thetalab_thetalab_gsgs_cm90", "h_thetalab_thetalab_gsgs_cm90", 200, 0, 100, 200, 0, 100);

   // ... .. phi vs phi
   TH2F *h_philab_philab = new TH2F("h_philab_philab", "h_philab_philab", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi = new TH2F("h_philab_philab_cutphi", "h_philab_philab_cutphi", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi_12c = new TH2F("h_philab_philab_cutphi_12c", "h_philab_philab_cutphi_12c", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_gsgs = new TH2F("h_philab_philab_gsgs", "h_philab_philab_gsgs", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_ver_gsgs = new TH2F("h_philab_philab_ver_gsgs", "h_philab_philab_ver_gsgs", 360, -180, 180, 360, -180, 180);

   // ... vertex of tracks
   TH1D *h_verz = new TH1D("h_verz", "h_verz;Vertex Z [mm]", 1020, -10, 1010);
   TH1D *h_verz01 = new TH1D("h_verz01", "h_verz;Vertex Z [mm]", 51, -0.1, 5);
   TH1D *h_verz_cut12c_ela = new TH1D("h_verz_cut12c_ela", "h_verz_cut12c_ela;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs = new TH1D("h_verz_gsgs", "h_verz_gsgs;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsex = new TH1D("h_verz_gsex", "h_verz_gsex;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_exex = new TH1D("h_verz_exex", "h_verz_exex;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs_cm50 = new TH1D("h_verz_gsgs_cm50", "h_verz_gsgs_cm50;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs_cm60 = new TH1D("h_verz_gsgs_cm60", "h_verz_gsgs_cm60;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs_cm70 = new TH1D("h_verz_gsgs_cm70", "h_verz_gsgs_cm70;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs_cm80 = new TH1D("h_verz_gsgs_cm80", "h_verz_gsgs_cm80;Vertex Z [mm]", 102, -10, 1010);
   TH1D *h_verz_gsgs_cm90 = new TH1D("h_verz_gsgs_cm90", "h_verz_gsgs_cm90;Vertex Z [mm]", 102, -10, 1010);

   TH2F *h_ntraver_ntra = new TH2F("h_ntraver_ntra", "h_ntraver_ntra", 11, -0.5, 10.5, 11, -0.5, 10.5);
   TH2F *h_ntraver_ntra_cutphi = new TH2F("h_ntraver_ntra_cutphi", "h_ntraver_ntra_cutphi", 11, -0.5, 10.5, 11, -0.5, 10.5);
   TH2F *h_ntraver_ntra_cutphi_2tra = new TH2F("h_ntraver_ntra_cutphi_2tra", "h_ntraver_ntra_cutphi_2tra", 11, -0.5, 10.5, 11, -0.5, 10.5);
   TH2F *h_ntraver_ntra_cut12c_ela = new TH2F("h_ntraver_ntra_cut12c_ela", "h_ntraver_ntra_cut12c_ela", 11, -0.5, 10.5, 11, -0.5, 10.5);
   //   TH2F *h_ntra_verz = new TH2F("h_ntra_verz", "h_ntra_verz", 112, -10, 1010, 11, -0.5, 10.5);
   TH2F *h_ntra_verz = new TH2F("h_ntra_verz", "h_ntra_verz", 1020, -10, 1010, 11, -0.5, 10.5);
   TH2F *h_verxy = new TH2F("h_verxy", "h_verxy", 50, -50, 50,  50, -50, 50);
   TH2F *h_verxz = new TH2F("h_verxz", "h_verxz", 102, -10, 1010, 50, -50, 50);
   TH2F *h_veryz = new TH2F("h_veryz", "h_veryz", 102, -10, 1010, 50, -50, 50);
   TH2F *h_verxy_cut12c_ela = new TH2F("h_verxy_cut12c_ela", "h_verxy_cut12c_ela", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_gsgs = new TH2F("h_verxy_gsgs", "h_verxy_gsgs", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_gsex = new TH2F("h_verxy_gsex", "h_verxy_gsex", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_exex = new TH2F("h_verxy_exex", "h_verxy_exex", 50, -50, 50, 50, -50, 50);

   // ... Beam energy estimation
   TH1D *h_Ebeam = new TH1D("h_Ebeam", "h_Ebeam;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs = new TH1D("h_Ebeam_gsgs", "h_Ebeam_gsgs;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs_cm50 = new TH1D("h_Ebeam_gsgs_cm50", "h_Ebeam_gsgs_cm50;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs_cm60 = new TH1D("h_Ebeam_gsgs_cm60", "h_Ebeam_gsgs_cm60;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs_cm70 = new TH1D("h_Ebeam_gsgs_cm70", "h_Ebeam_gsgs_cm70;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs_cm80 = new TH1D("h_Ebeam_gsgs_cm80", "h_Ebeam_gsgs_cm80;E_{beam} [MeV]", 70, 0, 70);
   TH1D *h_Ebeam_gsgs_cm90 = new TH1D("h_Ebeam_gsgs_cm90", "h_Ebeam_gsgs_cm90;E_{beam} [MeV]", 70, 0, 70);
   TH2F *h_sumkine_ver_Ebeam_cut12c_ela = new TH2F("h_sumkine_ver_Ebeam_cut12c_ela", "h_sumkine_Ebeam_gsgs;E_{beam} [MeV];sum_kineE [MeV]", 400, 0, 90, 400, 0, 90);
   TH2F *h_sumkine_Ebeam_gsgs = new TH2F("h_sumkine_Ebeam_gsgs", "h_sumkine_Ebeam_gsgs;E_{beam} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);
   TH2F *h_sumkine_Ebeam_gsgs_cm90 = new TH2F("h_sumkine_Ebeam_gsgs_cm90", "h_sumkine_Ebeam_gsgs_cm90;E_{beam} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);
   TH2F *h_sumkine_ver_Ebeam_gsgs = new TH2F("h_sumkine_ver_Ebeam_gsgs", "h_sumkine_Ebeam_gsgs;E_{beam} [MeV];sum_kineE [MeV]", 160, 0, 80, 160, 0, 80);

   TH1D *h_Ebcm = new TH1D("h_Ebcm", "h_Ebcm;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs = new TH1D("h_Ebcm_gsgs", "h_Ebcm_gsgs;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs_cm50 = new TH1D("h_Ebcm_gsgs_cm50", "h_Ebcm_gsgs_cm50;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs_cm60 = new TH1D("h_Ebcm_gsgs_cm60", "h_Ebcm_gsgs_cm60;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs_cm70 = new TH1D("h_Ebcm_gsgs_cm70", "h_Ebcm_gsgs_cm70;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs_cm80 = new TH1D("h_Ebcm_gsgs_cm80", "h_Ebcm_gsgs_cm80;E_{bc} [MeV]", 80, 0, 40);
   TH1D *h_Ebcm_gsgs_cm90 = new TH1D("h_Ebcm_gsgs_cm90", "h_Ebcm_gsgs_cm90;E_{bc} [MeV]", 80, 0, 40);

   /*
   // ... Excitation energy 
   TH1F *h_Exdp = new TH1F("h_Exdp", "h_Exdp;Ex [MeV]", 80, -5, 15);
   TH1F *h_Exdp_extended = new TH1F("h_Exdp_extended", "h_Exdp_extended;Ex [MeV]", 80, -5, 15);
   TH1F *h_Exdp_CarbonSi = new TH1F("h_Exdp_CarbonSi", "h_Exdp_CarbonSi;Ex [MeV]", 80, -5, 15);
   TH1F *h_Exdd = new TH1F("h_Exdd", "h_Exdd;Ex [MeV]",80 , -5, 15);

   // ... center-of-mass angular distributions
   TH1F *h_AngDist_elastic = new TH1F("h_AngDist_elastic", "h_AngDist_elastic;#theta_{c.m.} [deg]", 180, 0, 180);
   TH1F *h_AngDist_dp = new TH1F("h_AngDist_dp", "h_AngDist_dp;#theta_{c.m.} [deg];Counts / deg", 180, 0, 180);
   TH1F *h_AngDist_dp_CarbonSi = new TH1F("h_AngDist_dp_CarbonSi", "h_AngDist_dp_CarbonSi;#theta_{c.m.} [deg];Counts / deg", 180, 0, 180);
   */
   // ... others ...
   TH2F *h_Esmall_range = new TH2F("h_ESmallVTotalRange", "h_ESmallVTotalRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *h_Ebig_rangebig = new TH2F("h_EBigVBigRange", "h_EBigVBigRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *h_kineE_thetalab_2H  = new TH2F("h_kineE_thetalab_2H", "h_kineE_thetalab_2H", 180, 0, 180, 250, 0, 20);
   TH2F *h_kineE_thetalab_1H  = new TH2F("h_kineE_thetalab_1H", "h_kineE_thetalab_1H", 180, 0, 180, 250, 0, 20);
   TH2F *h_lastpoint_comp = new TH2F("h_lastpoint_comp", "h_lastpoint_comp;first_data [mm];second_data [mm]", 150, 0, 1500, 150, 0, 1500);
   TH2F *h_vertex_comp = new TH2F("h_vertex_comp", "h_vertex_comp;roughRange [mm];vertexRange [mm]", 150, 0, 1500, 150, 0, 1500);
   TH2F *h_distance_comp = new TH2F("h_distance_comp", "h_distance_comp;roughRange [mm];vertexRange [mm]", 150, 0, 1500, 150, 0, 1500);
   TH2F *h_range_comp = new TH2F("h_range_comp", "h_range_comp;roughRange [mm];vertexRange [mm]", 151, -10, 1500, 151, -10, 1500);
   TH2F *h_thetalab_comp = new TH2F("h_thetalab_comp", "h_thetalab_comp", 200, 0, 200, 200, 0, 200);
   TH2F *h_philab_comp = new TH2F("h_philab_comp", "h_philab_comp", 200, -200, 200, 200, -200, 200);
   TH2F *h_range_comp_gsgs = new TH2F("h_range_comp_gsgs", "h_range_comp_gsgs;roughRange [mm];vertexRange [mm]", 151, -10, 1500, 151, -10, 1500);
   TH2F *h_kineE_comp_gsgs = new TH2F("h_kineE_comp_gsgs", "h_kineE_comp_gsgs;first_data [MeV];second_data [MeV]", 70, 0, 70, 70, 0, 70);
   TH2F *h_thetalab_comp_gsgs = new TH2F("h_thetalab_comp_gsgs", "h_thetalab_comp_gsgs", 200, 0, 200, 200, 0, 200);
   TH2F *h_philab_comp_gsgs = new TH2F("h_philab_comp_gsgs", "h_philab_comp_gsgs", 200, -200, 200, 200, -200, 200);

   // ... index hist.
#if defined (vertex_index) || defined (vertex_depth)
   std::vector<TH1D*> h_verz_gsgs_index(n_group);
   std::vector<TH1D*> h_verz_gsex_index(n_group);
   std::vector<TH1D*> h_verz_exex_index(n_group);
   std::vector<TH2F*> h_verxy_gsgs_index(n_group);
   std::vector<TH2F*> h_verxy_gsex_index(n_group);
   std::vector<TH2F*> h_verxy_exex_index(n_group);
   std::vector<TH2F*> h_kineE_thetalab_gsgs_index(n_group);
   std::vector<TH2F*> h_kineE_thetalab_gsex_index(n_group);
   std::vector<TH2F*> h_kineE_thetalab_exex_index(n_group);
   std::vector<TH2F*> h_thetalab_thetalab_gsgs_index(n_group);
   std::vector<TH2F*> h_thetalab_thetalab_gsex_index(n_group);
   std::vector<TH2F*> h_thetalab_thetalab_exex_index(n_group);

   for (int i = 0; i < n_group; i++){
      if (i == n_group - 1){
         h_verz_gsgs_index[i] 
            = new TH1D("h_verz_gsgs_index_other", "h_verz_gsgs_z_others;Vertex Z [mm]", 112, -10, 1010);
         h_verz_gsex_index[i] 
            = new TH1D("h_verz_gsex_index_other", "h_verz_gsex_z_others;Vertex Z [mm]", 112, -10, 1010);
         h_verz_exex_index[i] 
            = new TH1D("h_verz_exex_index_other", "h_verz_exex_z_others;Vertex Z [mm]", 112, -10, 1010);
         h_verxy_gsgs_index[i] 
            = new TH2F("h_verxy_gsgs_index_other", "h_verxy_gsgs_z_others;Vertex X [mm];Vertex Y [mm]", 50, -50, 50, 50, -50, 50);
         h_verxy_gsex_index[i] 
            = new TH2F("h_verxy_gsex_index_other", "h_verxy_gsex_z_others;Vertex X [mm];Vertex Y [mm]", 50, -50, 50, 50, -50, 50);
         h_verxy_exex_index[i] 
            = new TH2F("h_verxy_exex_index_other", "h_verxy_exex_z_others;Vertex X [mm];Vertex Y [mm]", 50, -50, 50, 50, -50, 50);
         h_kineE_thetalab_gsgs_index[i] 
            = new TH2F("h_kineE_thetalab_gsgs_index_other","h_kineE_gsgs_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 700, 0, 70);
         h_kineE_thetalab_gsex_index[i] 
            = new TH2F("h_kineE_thetalab_gsex_index_other","h_kineE_gsex_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 700, 0, 70);
         h_kineE_thetalab_exex_index[i] 
            = new TH2F("h_kineE_thetalab_exex_index_other","h_kineE_exex_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 100, 0, 100, 700, 0, 70);
         h_thetalab_thetalab_gsgs_index[i]
            = new TH2F("h_thetalab_thetalab_gsgs_index_other","h_thetalab_thetalab_gsgs_z_others;#theta_{LAB} [deg];#theta_{LAB} [deg]", 200, 0, 100, 200, 0, 100);
         h_thetalab_thetalab_gsex_index[i]
            = new TH2F("h_thetalab_thetalab_gsex_index_other","h_thetalab_thetalab_gsex_z_others;#theta_{LAB} [deg];#theta_{LAB} [deg]", 200, 0, 100, 200, 0, 100);
         h_thetalab_thetalab_exex_index[i]
            = new TH2F("h_thetalab_thetalab_exex_index_other","h_thetalab_thetalab_exex_z_others;#theta_{LAB} [deg];#theta_{LAB} [deg]", 200, 0, 100, 200, 0, 100);
      }
      else {
         h_verz_gsgs_index[i] 
            = new TH1D(Form("h_verz_gsgs_index_%d", i),
               Form("h_verz_gsgs_z_%d--%d;Vertex Z [mm]", i*100, (i+1)*100), 112, -10, 1010);
         h_verz_gsex_index[i] 
            = new TH1D(Form("h_verz_gsex_index_%d", i),
               Form("h_verz_gsex_z_%d--%d;Vertex Z [mm]", i*100, (i+1)*100), 112, -10, 1010);
         h_verz_exex_index[i] 
            = new TH1D(Form("h_verz_exex_index_%d", i), 
               Form("h_verz_exex_z_%d--%d;Vertex Z [mm]", i*100, (i+1)*100), 112, -10, 1010);
         h_verxy_gsgs_index[i] 
            = new TH2F(Form("h_verxy_gsgs_index_%d", i),
               Form("h_verxy_gsgs_z_%d--%d;Vertex X [mm];Vertex Y [mm]", i*100, (i+1)*100), 50, -50, 50, 50, -50, 50);
         h_verxy_gsex_index[i] 
            = new TH2F(Form("h_verxy_gsex_index_%d", i),
               Form("h_verxy_gsex_z_%d--%d;Vertex X [mm];Vertex Y [mm]", i*100, (i+1)*100), 50, -50, 50, 50, -50, 50);
         h_verxy_exex_index[i] 
            = new TH2F(Form("h_verxy_exex_index_%d", i),
               Form("h_verxy_exex_z_%d--%d;Vertex X [mm];Vertex Y [mm]", i*100, (i+1)*100), 50, -50, 50, 50, -50, 50);
         h_kineE_thetalab_gsgs_index[i] 
            = new TH2F(Form("h_kineE_thetalab_gsgs_index_%d", i),
               Form("h_kineE_thetalab_gsgs_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 100, 0, 100, 700, 0,70);
         h_kineE_thetalab_gsex_index[i] 
            = new TH2F(Form("h_kineE_thetalab_gsex_index_%d", i),
               Form("h_kineE_thetalab_gsex_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 100, 0, 100, 700, 0, 70);
         h_kineE_thetalab_exex_index[i] 
            = new TH2F(Form("h_kineE_thetalab_exex_index_%d", i),
               Form("h_kineE_thetalab_exex_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 100, 0, 100, 700, 0, 70);
         h_thetalab_thetalab_gsgs_index[i]
            = new TH2F(Form("h_thetalab_thetalab_gsgs_index_%d", i),
               Form("h_thetalab_thetalab_gsgs_z_%d--%d;#theta_{LAB} [deg];#theta_{LAB} [deg]", i*100, (i+1)*100), 200, 0, 100, 200, 0, 100);
         h_thetalab_thetalab_gsex_index[i]
            = new TH2F(Form("h_thetalab_thetalab_gsex_index_%d", i),
               Form("h_thetalab_thetalab_gsex_z_%d--%d;#theta_{LAB} [deg];#theta_{LAB} [deg]", i*100, (i+1)*100), 200, 0, 100, 200, 0, 100);
         h_thetalab_thetalab_exex_index[i]
            = new TH2F(Form("h_thetalab_thetalab_exex_index_%d", i),
               Form("h_thetalab_thetalab_exex_z_%d--%d;#theta_{LAB} [deg];#theta_{LAB} [deg]", i*100, (i+1)*100), 200, 0, 100, 200, 0, 100);
      }
   }

#endif

#ifdef catima_check
   std::cout << "Catima check... Before loop of run. " << std::flush;
   r = eLossModelC4H10_12C->GetRange(0.1);
   std::cout << " range = " << r << ", density = " << eLossModelC4H10_12C->GetDensity()
            << ", A = " << eLossModelC4H10_12C->GetAtomicMassNumber() << ", Z = " << eLossModelC4H10_12C->GetChargeNumber()
            << ", mass = " << eLossModelC4H10_12C->GetMassAmu() << ", dEdx = " << eLossModelC4H10_12C->GetdEdx(0.1) << std::endl;
#endif

   std::cout << std::endl;
   for (int runNum: runNums) {
      // Open the digitalization file and get the TTree.
      //      TString unpackFileName = TString::Format("./vd_check_data/run_%04d.root", runNum);
      TString unpackFileName = TString::Format("./vd_check_data/vd%.2f_files/run_%04d.root", vd_val, runNum);
      TFile *unpackFile = new TFile(unpackFileName, "READ");
      TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
      int nUnpackEvents = unpackTree->GetEntries();
      //      std::cout << "Number of unpacked events in run " << runNum << ": " << nUnpackEvents << std::endl;
      std::cout << "Reading files: " << nruns << "/" << runNums.size() << "; run: " << runNum
                  << ", events: " << nUnpackEvents << std::string(50, ' ') << "\r" << std::flush;
      int nEventsWith2Tracks = 0;
      // Creare the TTreeReader to read the AtTrackingEvents and simulation.
      TTreeReader unpackReader("cbmsim", unpackFile);
      //      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");
      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEventModified");

      // Loop over events.
#ifdef catima_check
      std::cout << "Catima check... Before loop of event. " << std::flush;
      r = eLossModelC4H10_12C->GetRange(0.1);
      std::cout << " range = " << r << ", density = " << eLossModelC4H10_12C->GetDensity()
               << ", A = " << eLossModelC4H10_12C->GetAtomicMassNumber() << ", Z = " << eLossModelC4H10_12C->GetChargeNumber()
               << ", mass = " << eLossModelC4H10_12C->GetMassAmu() << ", dEdx = " << eLossModelC4H10_12C->GetdEdx(0.1) << std::endl;
#endif

      for (int i = 0; i < nUnpackEvents; i++) {
         unpackReader.Next();
         // First, we obtain some rough kinematics just by using the AtPatternEvent.
         AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
         if (!patternEvent){continue;}

         auto &tracks = patternEvent->GetTrackCand();
         ntrack = tracks.size();
         h_ntra->Fill(ntrack);
         std::vector<bool> track_braggd(ntrack, false);
         std::vector<bool> track_12c(ntrack, false);
         std::vector<bool> track_alpha(ntrack, false);
         std::vector<bool> track_proton(ntrack, false);
         tracks_vertex = false;
         alpha_tracks = false;
         proton_tracks = false;
         rad = 0;
         r_max = 0;
         r_tem = 0;
         itrack = 0;
         nbragg = 0;
         nbrain = 0;
         nbrano = 0;
         vindex = 0;
         n_bragg_true = 0;
         n_bragg_false = 0;
         nalpha = 0;
         nproton = 0;
         vtx = 0;
         vty = 0;
         vtz = 0;
         for (Int_t i = 0; i < narray; i ++){
            vertex_theta[i] = 0;
            vertex_phi[i] = 0;
            vertex_range[i] = 0;
            vertex_KinE[i] = 0;
         }
         vertexPoint.SetXYZ(0, 0, 0);
         vertexPoint2.SetXYZ(0, 0, 0);
         v_distance.clear();
         v_lastPoint.clear();

         if (ntrack == 6){
            track6.push_back(i);
         }
         if (ntrack > narray){
            std::cout << "Warning: Too many tracks. (narray = " << narray << "), run: " << runNum << ", event: " << i 
                        << ", ntrack: " << ntrack << std::endl;
         }
         //         int maxTrackNum{4};
         //         if (tracks.size() > maxTrackNum) continue;
         //         if (tracks.size() == 2) nEventsWith2Tracks++;
         int trackIndex = 0;
         double thetaLABArray[2] = {0.0,0.0};
         double estimatedKineEArray[2] = {0.0,0.0};

#ifdef catima_check
         std::cout << "Catima check... Before loop of track. " << std::flush;
         r = eLossModelC4H10_12C->GetRange(0.1);
         std::cout << " range = " << r << ", density = " << eLossModelC4H10_12C->GetDensity()
                  << ", A = " << eLossModelC4H10_12C->GetAtomicMassNumber() << ", Z = " << eLossModelC4H10_12C->GetChargeNumber()
                  << ", mass = " << eLossModelC4H10_12C->GetMassAmu() << ", dEdx = " << eLossModelC4H10_12C->GetdEdx(0.1) << std::endl;
#endif

         // Iterate over AtTracks and extract their kinematics.
         for (auto &track: tracks) {
            bool isPunchThrough = punchThroughChecker.IsPunchThrough(&track);
            // if (isPunchThrough) continue;
            double smallPadCharge{};
            double bigPadCharge{};
            double rangeInSmallPads{};
            auto *pattern = track.GetPattern();
            auto firstPoint = track.GetFirstPoint();
            auto lastPoint = track.GetLastPoint();
            auto pseudoVertex = pattern->ClosestPointOnPattern(firstPoint);
            auto braggCurvePairs = track.GetBraggCurveValues();
            auto &hits = track.GetHitArray();
            auto braggCurve = track.GetBraggCurve();
            if(braggCurve.RangeValues.size() > 0){
               track_braggd[itrack] = true;
               //               std::cout << "test! run:" << runNum << ", event:"<< i << ", track:" << itrack << ", check_bragg == true !!"<< std::endl;
            }

            v_lastPoint.push_back(lastPoint);
            track_range[itrack] = pattern->DistanceAlongPattern(lastPoint, firstPoint);
            track_charge[itrack] = track.GetGeoQEnergy();
            //            track_theta[itrack] = 180 - track.GetGeoTheta() * 180 / TMath::Pi();
            track_theta[itrack] = track.GetGeoTheta() * 180 / TMath::Pi();
            track_phi[itrack] = track.GetGeoPhi() * 180 / TMath::Pi();

            track_lastx = lastPoint.X();
            track_lasty = lastPoint.Y();
            track_lastz = lastPoint.Z();
            rad = TMath::Sqrt(track_lastx * track_lastx + track_lasty * track_lasty);

            if(track_braggd[itrack]){
               vertx[itrack] = braggCurve.vertexX;
               verty[itrack] = braggCurve.vertexY;
               vertz[itrack] = 1000.0 - braggCurve.vertexZ;
               if(abs(vertx[itrack]) < 1e-6 && abs(verty[itrack]) < 1e-6 && abs(vertz[itrack] + 999 ) < 1e-6 ){
                  nbrain ++;
                  std::cout << " Something wrong with vertex! run:" << runNum << ", event;" << i << ", itracks:" << itrack 
                        << ", vetex: (" << vertx[itrack] << ", " << verty[itrack] << ", " << vertz[itrack] << ")"  << std::endl;
               }
               if(abs(vtx) < 1e-6 && abs(vty) < 1e-6 && abs(vtz) < 1e-6){
                  nbragg ++;
                  vtx = vertx[itrack];
                  vty = verty[itrack];
                  vtz = vertz[itrack];
                  vertexPoint.SetXYZ(vtx, vty, 1000.0 - vtz);
               }
               else if (abs(vertx[itrack] - vtx) < 1e-6 && abs(verty[itrack] - vty) < 1e-6 && abs(vertz[itrack] - vtz) < 1e-6){
                  nbragg ++;
               }
               else {
                  std::cout << std::endl;
                  std::cout << " Multiple vertices in one run! run:" << runNum << ", event;" << i << ", track:" << itrack << std::endl;
               }
               /*
               if(ntrack > 3){
                  std::cout << "test! run:" << runNum << ", event:"<< i << ", track:" << itrack 
                  << ", vertex: (" << vertx[itrack] <<",  " << verty[itrack] << ", " << vertz[itrack] << " )"<< std::endl;
               }
               */
            }
            else {
               nbrano ++;
               //               std::cout << "test! check_bragg == false !!"<< std::endl;
            }

            for (auto &hit: hits) {
               auto pos = hit->GetPosition();
               int padNum = hit->GetPadNum();
               int sizeID = map->GetPadSize(padNum);
               if (sizeID == 1) {
                  bigPadCharge += hit->GetCharge();
                  continue;
               }
               smallPadCharge += hit->GetCharge();
               double currentRangeInSmallPads = pattern->DistanceAlongPattern(hit->GetPosition(), firstPoint);
               if (currentRangeInSmallPads > rangeInSmallPads){
                  rangeInSmallPads = currentRangeInSmallPads;
               }
            }
            if (rad > r_tem){
               r_tem = rad;
            }
            track_r[itrack] = r_tem;
            if(track_r[itrack] > r_max){
               r_max = track_r[itrack];
            }

            track_dedx[itrack] = smallPadCharge / rangeInSmallPads;
            double rangeInBigPads = track_range[itrack] - rangeInSmallPads;
            bool reachedBigPads = true;
            if (rangeInBigPads / track_range[itrack] < 0.05){
               reachedBigPads = false;
            }

            h_charge_range->Fill(track_range[itrack],track_charge[itrack]);
            h_range_thetalab->Fill(track_theta[itrack], track_range[itrack]);
            h_dEdx_range->Fill(track_range[itrack], track_dedx[itrack]);
            /*
            if (!reachedBigPads){
               h_Esmall_range->Fill(track_range[itrack], smallPadCharge);
            }
            else {
               h_Ebig_rangebig->Fill(rangeInBigPads, bigPadCharge);
            }
            if (track_theta[itrack] > 100){
               h_dEdx_range_backwards->Fill(track_range[itrack], track_dedx[itrack]);
            }
            */

            //            double estimatedKinE{0.1};
            track_KinE[itrack] = 0.1;

#ifdef catima_check
            std::cout << "Catima check... Before calculation. " << std::flush;
            r = eLossModelC4H10_12C->GetRange(0.1);
            std::cout << " range = " << r << ", density = " << eLossModelC4H10_12C->GetDensity()
                     << ", A = " << eLossModelC4H10_12C->GetAtomicMassNumber() << ", Z = " << eLossModelC4H10_12C->GetChargeNumber()
                     << ", mass = " << eLossModelC4H10_12C->GetMassAmu() << ", dEdx = " << eLossModelC4H10_12C->GetdEdx(0.1) << std::endl;
#endif

            if (cut12c->IsInside(track_range[itrack], track_charge[itrack])){
#if defined (debug_mode) || defined (catima_check)
               const double measuredRange = track_range[itrack];
               const double modelRange = eLossModelC4H10_12C->GetRange(track_KinE[itrack]);
               std::cout <<  "Catima check. run:" << runNum << ", event:" << i << ", ntrack = " << ntrack << " itrack = " << itrack
                           << " measuredRange = " << measuredRange << " E = " << track_KinE[itrack] << " modelRange=" << modelRange
                           << " punchThrough=" << isPunchThrough << " first=(" << firstPoint.X() << "," << firstPoint.Y() << "," << firstPoint.Z()
                           << ") last=(" << lastPoint.X() << "," << lastPoint.Y() << "," << lastPoint.Z() << ")" << std::endl;
#endif
               while (eLossModelC4H10_12C->GetRange(track_KinE[itrack]) < track_range[itrack]){
                  track_KinE[itrack] += 0.01;
               }
            }
            else if (cutalpha->IsInside(track_range[itrack], track_charge[itrack])){
               while (eLossModelC4H10_alpha->GetRange(track_KinE[itrack]) < track_range[itrack]){
                  track_KinE[itrack] += 0.01;
               }
            }
            else if (cutproton->IsInside(track_range[itrack], track_charge[itrack])){
               while (eLossModelC4H10_p->GetRange(track_KinE[itrack]) < track_range[itrack]){
                  track_KinE[itrack] += 0.01;
               }
            }
            else {
               track_KinE[itrack] = -1;
            }

            // Kinematics plots
            h_kineE_thetalab->Fill(track_theta[itrack], track_KinE[itrack]);
            /*
            // .... 12c in ATTPC
            if (cutPIDdeuteron->IsInside(track_range[itrack], track_dedx[itrack])) {
               histSiPIDADCMax_DeuteronATTPC->Fill(maxADCFront2, maxADCFront1);
               histEstimatedKinEVThetaLAB_DeuteronATTPC->Fill(trackThetaLAB, track_KinE[itrack]);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_d, m_17C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), track_KinE[itrack]);
               histExdd->Fill(Ex);

               if ((Ex >= -4.75) && (Ex <= 4.73)) {
                  histAngDist_elastic->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));
               }
            }
            // .... alpha in ATTPC
            if (cutPIDdeuteron->IsInside(track_range[itrack], track_dedx[itrack])) {
               histSiPIDADCMax_DeuteronATTPC->Fill(maxADCFront2, maxADCFront1);
               histEstimatedKinEVThetaLAB_DeuteronATTPC->Fill(trackThetaLAB, track_KinE[itrack]);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_d, m_17C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), track_KinE[itrack]);
               histExdd->Fill(Ex);

               if ((Ex >= -4.75) && (Ex <= 4.73)) {
                  histAngDist_elastic->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));
               }
            }
            // .... protons in ATTPC
            if (cutPIDproton->IsInside(track_range[itrack], track_dedx[itrack])) {
               histSiPIDADCMax_ProtonATTPC->Fill(maxADCFront2, maxADCFront1);
               histEstimatedKinEVThetaLAB_ProtonATTPC->Fill(trackThetaLAB, track_KinE[itrack]);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_p, m_18C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), track_KinE[itrack]);
               histExdp->Fill(Ex);            
            }
            if (cutPIDproton->IsInside(track_range[itrack], track_dedx[itrack]) or cutPIDproton_extension->IsInside(track_range[itrack], track_dedx[itrack])) {
               histEstimatedKinEVThetaLAB_ProtonATTPC_extended->Fill(trackThetaLAB, track_KinE[itrack]);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_p, m_18C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), track_KinE[itrack]);
               histExdp_extended->Fill(Ex);

               histAngDist_dp->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));

               // Also add condition of Carbon in Si
               if (cutSiC->IsInside(maxADCFront2, maxADCFront1)) {
                  histEstimatedKinEVThetaLAB_ProtonATTPC_CarbonSi->Fill(trackThetaLAB, track_KinE[itrack]);
                  histExdp_CarbonSi->Fill(Ex);
                  histAngDist_dp_CarbonSi->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));
               }
            }
            */
            itrack ++;
         }

         // check bragg curve
         n_bragg_true = std::count(track_braggd.begin(), track_braggd.end(), true);
         n_bragg_false = std::count(track_braggd.begin(), track_braggd.end(), false);

         if(ntrack != nbragg + nbrain + nbrano){
            std::cout << " Something wrong with vertex! run:" << runNum << ", event;" << i << ", ntracks:" << ntrack 
                     << ", nbragg:" << nbragg << ", nbrain:" << nbrain << ", nbrano:" << nbrano << std::endl;
         }
         if(nbragg != n_bragg_true){
            std::cout << " Something wrong with vertex! run:" << runNum << ", event;" << i << ", ntracks:" << ntrack 
                     << ", nbragg:" << nbragg << ", nbrain:" << nbrain << ", nbrano:" << nbrano << std::endl;
         }

         if(ntrack == n_bragg_true){
            // all true
            tracks_vertex = true;
         }
         else if(ntrack == n_bragg_false){
            // all false
         }
         else{
            // mixed
            tracks_vertex = true;
            /*
            std::cout << "bragg check run:" << runNum << ", event;" << i 
               << ", ntracks:" << ntrack << ", nbragg:" << nbragg << ", nbrain:" << nbrain << ", nbrano:" << nbrano << std::endl;
            */
         }
         if(vtz >= 0 && vtz < verz_h){
            vindex = int (vtz/100.0);
         }
         else {
            vindex = n_group - 1;
         }

         h_rmax->Fill(r_max);
         if(r_tri > r_max && r_max > 0){
            r_tri = r_max;
         }
         if(r_max > max_r_max){
            max_r_max = r_max;
         }

         h_philab_philab -> Fill(track_phi[0], track_phi[1]);
         h_ntraver_ntra -> Fill(ntrack, n_bragg_true);
         if(tracks_vertex){
            h_verxy -> Fill(vtx, vty);
            h_verz -> Fill(vtz);
            h_verxz -> Fill(vtz, vtx);
            h_veryz -> Fill(vtz, vty);
            h_ntra_verz -> Fill(vtz, ntrack);
            itrack = 0;
            for(auto &track: tracks){
               auto *pattern = track.GetPattern();
               auto lastPoint = track.GetLastPoint();
               auto direvertra = lastPoint - vertexPoint;
               vertex_range[itrack] = pattern->DistanceAlongPattern(lastPoint, vertexPoint);
               vertex_theta[itrack] = 180 - direvertra.Theta() * 180 / TMath::Pi();
               vertex_phi[itrack] = direvertra.Phi() * 180 / TMath::Pi();
               if (kine_comp_b){
                  h_lastpoint_comp->Fill(v_lastPoint[itrack].R(), lastPoint.R());
                  h_distance_comp -> Fill((vertexPoint - v_lastPoint[itrack]).R(), (vertexPoint - lastPoint).R());
                  h_range_comp->Fill(track_range[itrack], vertex_range[itrack]);
                  h_thetalab_comp->Fill(track_theta[itrack], vertex_theta[itrack]);
                  h_philab_comp->Fill(track_phi[itrack], vertex_phi[itrack]);
                  auto braggCurve = track.GetBraggCurve();
                  if (braggCurve.RangeValues.size() > 0){
                     vertexPoint2.SetXYZ(braggCurve.vertexX, braggCurve.vertexY, braggCurve.vertexZ);
                     h_vertex_comp->Fill(vertexPoint.R(), vertexPoint2.R());
                  }
                  /*
                  if (abs(vertexPoint2.R() - 1000.0) <= 1){
                     std::cout << " vertexPoint2 check. run: " << runNum << ", event: " << i << ", itrack: " << itrack << std::endl;
                     std::cout << "    vertexPoint  = (" << vertexPoint.X() << ", " << vertexPoint.Y() << ", " << vertexPoint.Z() << ")" << std::endl;
                     std::cout << "    vertexPoint2 = (" << vertexPoint2.X() << ", " << vertexPoint2.Y() << ", " << vertexPoint2.Z() << ")" << std::endl;
                  }
                  */
               }
               vertex_KinE[itrack] = 0.1;
               if (cut12c->IsInside(track_range[itrack], track_charge[itrack])){
                  const Double_t warmupRange12C_2 = eLossModelC4H10_12C->GetRange(0.1);
                  if (abs(warmupRange12C - warmupRange12C_2) > 1e-6){
                     std::cout << "Something worng with warmupRange12C. run: " << runNum << ", event: " << i
                              << ", itrack: " << itrack << ", warmupRange12C: " << warmupRange12C
                              << ", warmupRange12C_2: " << warmupRange12C_2 << ", distance: "
                              << abs(warmupRange12C - warmupRange12C_2) << std::endl;
                  }
                  while (eLossModelC4H10_12C->GetRange(vertex_KinE[itrack]) < vertex_range[itrack]){
                     vertex_KinE[itrack] += 0.01;
                  }
               }
               else {
                  vertex_KinE[itrack] = -1.0;
               }
               itrack ++;
            }
            if (ntrack != itrack){
               std::cout << " Something wrong with kineE calculation by using vertex!" << std::endl;
               std::cout << "  run:" << runNum << ", event;" << i << ", ntracks:" << ntrack << ", itrack:" << itrack << std::endl;
            }
         }
         if(abs(abs(track_phi[0] - track_phi[1]) - 180) < del_phi){
            // pid
            for (Int_t k = 0; k < ntrack; k++){
               if (cut12c->IsInside(track_range[k], track_charge[k])) {
                  //                  std::cout << "event num:" << i << ", tracks: " << ntrack << ", track: "<< k << ", particle: 12c " << std::endl;
                  track_12c[k]=true;
                  h_philab_philab_cutphi_12c -> Fill(track_phi[0], track_phi[1]);                  
               }
            }
            // fill to histograms
            h_charge_range_cutphi -> Fill(track_range[0], track_charge[0]);
            h_charge_range_cutphi -> Fill(track_range[1], track_charge[1]);
            h_dEdx_range_cutphi->Fill(track_range[0], track_dedx[0]);
            h_dEdx_range_cutphi->Fill(track_range[1], track_dedx[1]);
            h_range_thetalab_cutphi -> Fill(track_theta[0], track_range[0]);
            h_range_thetalab_cutphi -> Fill(track_theta[1], track_range[1]);
            h_thetalab_thetalab_cutphi -> Fill(track_theta[0], track_theta[1]);
            h_philab_philab_cutphi -> Fill(track_phi[0], track_phi[1]);
            h_ntraver_ntra_cutphi -> Fill(ntrack, n_bragg_true);
            if (ntrack == 2){
               h_charge_range_cutphi_2tra -> Fill(track_range[0], track_charge[0]);
               h_charge_range_cutphi_2tra -> Fill(track_range[1], track_charge[1]);
               h_range_thetalab_cutphi_2tra -> Fill(track_theta[0], track_range[0]);
               h_range_thetalab_cutphi_2tra -> Fill(track_theta[1], track_range[1]);
               h_thetalab_thetalab_cutphi_2tra -> Fill(track_theta[0], track_theta[1]);
               h_ntraver_ntra_cutphi_2tra -> Fill(ntrack, n_bragg_true);
               if (track_12c[0] && track_12c[1]){
                  Double_t sum_theta = track_theta[0] + track_theta[1];
                  h_charge_range_cut12c_ela -> Fill(track_range[0], track_charge[0]);
                  h_charge_range_cut12c_ela -> Fill(track_range[1], track_charge[1]);
                  h_range_thetalab_cut12c_ela -> Fill(track_theta[0], track_range[0]);
                  h_range_thetalab_cut12c_ela -> Fill(track_theta[1], track_range[1]);
                  h_kineE_thetalab_carbon -> Fill(track_theta[0], track_KinE[0]);
                  h_kineE_thetalab_carbon -> Fill(track_theta[1], track_KinE[1]);
                  h_thetalab_thetalab_cut12c_ela -> Fill(track_theta[0], track_theta[1]);
                  h_sum_theta_cut12c_ela -> Fill(sum_theta);
                  h_ntraver_ntra_cut12c_ela -> Fill(ntrack, n_bragg_true);
                  h_rmax_12c->Fill(r_max);
                  if(r_tri_12c > r_max && r_max > 0){
                     r_tri_12c = r_max;
                  }
                  if(r_max > max_r_max_12c){
                     max_r_max_12c = r_max;
                  }
                  n_12c12c ++;
                  if(tracks_vertex){
                     h_verxy_cut12c_ela -> Fill(vtx, vty);
                     h_verz_cut12c_ela -> Fill(vtz);
                     Ebeam = est_Ebeam(Ebeam_para, vtz);
                     Ebeam_cm = Ebeam/2.0;
                     h_Ebeam->Fill(Ebeam);
                     h_Ebcm->Fill(Ebeam_cm);
                     h_sumkine_verz_cut12c_ela->Fill(vtz, vertex_KinE[0] + vertex_KinE[1]);
                     h_sumkine_ver_Ebeam_cut12c_ela->Fill(Ebeam, vertex_KinE[0] + vertex_KinE[1]);
#ifdef check_verz0
                     if(vtz > 0 && vtz <= 1){
                        std::cout << "  Vertex z ~ 0 (12c12c); run: " << runNum << ", event: " << i << ", vertexz = " << vtz << std::string(10, ' ') << std::endl;
                        vertex_z_01[0].push_back(runNum);
                        vertex_z_01[1].push_back(i);
                        h_verz01 -> Fill(vtz);
                        n_ver_z_01 ++;
                     }
#endif
                  }
                  else {
#ifdef vertex_check
                     n_no_vertex ++;
                     std::cout << "  No vertex (12c12c); run: " << runNum << ", event:" << i << ", count:" << n_no_vertex << std::string(30, ' ') << std::endl;
                     track2_ver0[0].push_back(runNum);
                     track2_ver0[1].push_back(i);
#endif
                  }
                  if (runNum == 52){
                     h_sum_theta_cut12c_run52 -> Fill(sum_theta);
                  }
                  //                  if (sum_theta > 76.0 && sum_theta < 84.0){
                  if (theta_gsgs->IsInside(track_theta[0],track_theta[1])){
                     if (runNum == 52){
                        peak1.push_back(i);
                     }
                     h_sum_theta_gsgs -> Fill(sum_theta);
                     if (vertex_corr){
                        // Using new kinematics energy by calculation of new range from vertex to last point of track
                        h_range_thetalab_gsgs -> Fill(track_theta[0], vertex_range[0]);
                        h_range_thetalab_gsgs -> Fill(track_theta[1], vertex_range[1]);
                        h_kineE_thetalab_gsgs -> Fill(track_theta[0], vertex_KinE[0]);
                        h_kineE_thetalab_gsgs -> Fill(track_theta[1], vertex_KinE[1]);
                     }
                     else if (!vertex_corr){
                        // Using rough kinematics energy by calculation of rough range from first to last point of track
                        h_range_thetalab_gsgs -> Fill(track_theta[0], track_range[0]);
                        h_range_thetalab_gsgs -> Fill(track_theta[1], track_range[1]);
                        h_kineE_thetalab_gsgs -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_gsgs -> Fill(track_theta[1], track_KinE[1]);
                     }
                     if (!kine_comp_b){
                        if (vertex_corr){
                           h_sumkine_kineE_gsgs->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[0]);
                           h_sumkine_kineE_gsgs->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[1]);
                           h_thetalab_thetalab_gsgs -> Fill(track_theta[0], track_theta[1]);
                           h_philab_philab_gsgs -> Fill(track_phi[0], track_phi[1]);
                        }
                        else if (!vertex_corr){
                           h_sumkine_kineE_gsgs->Fill(track_KinE[0] + track_KinE[1], track_KinE[0]);
                           h_sumkine_kineE_gsgs->Fill(track_KinE[0] + track_KinE[1], track_KinE[1]);
                           h_thetalab_thetalab_gsgs -> Fill(track_theta[0], track_theta[1]);
                           h_philab_philab_gsgs -> Fill(track_phi[0], track_phi[1]);
                        }
                     }
                     if(tracks_vertex){
                        h_verz_gsgs->Fill(vtz);
                        h_verxy_gsgs->Fill(vtx, vty);
                        h_Ebeam_gsgs->Fill(Ebeam);
                        h_Ebcm_gsgs->Fill(Ebeam_cm);
                        if(!kine_comp_b){
                           if (vertex_corr){
                              h_sumkine_verz_gsgs->Fill(vtz, vertex_KinE[0] + vertex_KinE[1]);
                              h_sumkine_Ebeam_gsgs->Fill(Ebeam, vertex_KinE[0] + vertex_KinE[1]);
                           }
                           else if (!vertex_corr){
                              h_sumkine_verz_gsgs->Fill(vtz, track_KinE[0] + track_KinE[1]);
                              h_sumkine_Ebeam_gsgs->Fill(Ebeam, track_KinE[0] + track_KinE[1]);
                           }
                        }
                        else if (kine_comp_b){
                           h_thetalab_thetalab_gsgs -> Fill(track_theta[0], track_theta[1]);
                           h_philab_philab_gsgs -> Fill(track_phi[0], track_phi[1]);
                           h_thetalab_thetalab_ver_gsgs -> Fill(vertex_theta[0], vertex_theta[1]);
                           h_philab_philab_ver_gsgs -> Fill(vertex_phi[0], vertex_phi[1]);
                           h_range_comp_gsgs->Fill(track_range[0], vertex_range[0]);
                           h_range_comp_gsgs->Fill(track_range[1], vertex_range[1]);
                           h_kineE_comp_gsgs->Fill(track_KinE[0], vertex_KinE[0]);
                           h_kineE_comp_gsgs->Fill(track_KinE[1], vertex_KinE[1]);
                           h_thetalab_comp_gsgs->Fill(track_theta[0], vertex_theta[0]);
                           h_thetalab_comp_gsgs->Fill(track_theta[1], vertex_theta[1]);
                           h_philab_comp_gsgs->Fill(track_phi[0], vertex_phi[0]);
                           h_philab_comp_gsgs->Fill(track_phi[1], vertex_phi[1]);
                           h_sumkine_verz_gsgs->Fill(vtz, vertex_KinE[0] + vertex_KinE[1]);
                           h_sumkine_Ebeam_gsgs->Fill(Ebeam, track_KinE[0] + track_KinE[1]);
                           h_sumkine_ver_Ebeam_gsgs->Fill(Ebeam, vertex_KinE[0] + vertex_KinE[1]);
                           h_sumkine_kineE_gsgs->Fill(track_KinE[0] + track_KinE[1], track_KinE[0]);
                           h_sumkine_kineE_gsgs->Fill(track_KinE[0] + track_KinE[1], track_KinE[1]);
                           h_sumkine_kineE_ver_gsgs->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[0]);
                           h_sumkine_kineE_ver_gsgs->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[1]);
                        }
                        if ((track_theta[0] > 23 && track_theta[0] < 27) || (track_theta[1] > 23 && track_theta[1] < 27)){
                           h_thetalab_thetalab_gsgs_cm50 -> Fill(track_theta[0], track_theta[1]);
                           h_verz_gsgs_cm50 -> Fill(vtz);
                           h_Ebeam_gsgs_cm50->Fill(Ebeam);
                           h_Ebcm_gsgs_cm50->Fill(Ebeam_cm);
                        }
                        if ((track_theta[0] > 28 && track_theta[0] < 32) || (track_theta[1] > 28 && track_theta[1] < 32)){
                           h_thetalab_thetalab_gsgs_cm60 -> Fill(track_theta[0], track_theta[1]);
                           h_verz_gsgs_cm60 -> Fill(vtz);
                           h_Ebeam_gsgs_cm60->Fill(Ebeam);
                           h_Ebcm_gsgs_cm60->Fill(Ebeam_cm);
                        }
                        if ((track_theta[0] > 33 && track_theta[0] < 37) || (track_theta[1] > 33 && track_theta[1] < 37)){
                           h_thetalab_thetalab_gsgs_cm70 -> Fill(track_theta[0], track_theta[1]);
                           h_verz_gsgs_cm70 -> Fill(vtz);
                           h_Ebeam_gsgs_cm70->Fill(Ebeam);
                           h_Ebcm_gsgs_cm70->Fill(Ebeam_cm);
                        }
                        if ((track_theta[0] > 38 && track_theta[0] < 42) || (track_theta[1] > 38 && track_theta[1] < 42)){
                           h_thetalab_thetalab_gsgs_cm80 -> Fill(track_theta[0], track_theta[1]);
                           h_verz_gsgs_cm80 -> Fill(vtz);
                           h_Ebeam_gsgs_cm80->Fill(Ebeam);
                           h_Ebcm_gsgs_cm80->Fill(Ebeam_cm);
                        }
                        if ((track_theta[0] > 43 && track_theta[0] < 47) || (track_theta[1] > 43 && track_theta[1] < 47)){
                           h_verz_gsgs_cm90 -> Fill(vtz);
                           h_Ebeam_gsgs_cm90->Fill(Ebeam);
                           h_Ebcm_gsgs_cm90->Fill(Ebeam_cm);
                           h_thetalab_thetalab_gsgs_cm90 -> Fill(track_theta[0], track_theta[1]);
                           if (vertex_corr){
                              h_sumkine_kineE_gsgs_cm90->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[0]);
                              h_sumkine_kineE_gsgs_cm90->Fill(vertex_KinE[0] + vertex_KinE[1], vertex_KinE[1]);
                              h_sumkine_Ebeam_gsgs_cm90->Fill(Ebeam, vertex_KinE[0] + vertex_KinE[1]);
                           }
                           else if (!vertex_corr){
                              h_sumkine_kineE_gsgs_cm90->Fill(track_KinE[0] + track_KinE[1], track_KinE[0]);
                              h_sumkine_kineE_gsgs_cm90->Fill(track_KinE[0] + track_KinE[1], track_KinE[1]);
                              h_sumkine_Ebeam_gsgs_cm90->Fill(Ebeam, track_KinE[0] + track_KinE[1]);
                           }
                        }
#if defined (vertex_index) || defined (vertex_depth)
                        h_verz_gsgs_index[vindex] -> Fill(vtz);
                        h_verxy_gsgs_index[vindex] -> Fill(vtx, vty);
                           h_thetalab_thetalab_gsgs_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                        if (vertex_corr){
                           // For using new kinematics energy by calculation of new range from vertex to last point of track
                           h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[0], vertex_KinE[0]);
                           h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[1], vertex_KinE[1]);
                        }
                        else if (!vertex_corr){
                           // For using rough kinematics energy by calculation of rough range from first to last point of track
                           h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                           h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
                        }
#endif
                     }
                  } 
                  else if (theta_gsex->IsInside(track_theta[0],track_theta[1])){
                     if (runNum == 52){
                        peak2.push_back(i);
                     }
                     h_sum_theta_gsex -> Fill(sum_theta);
                     if (vertex_corr){
                        h_range_thetalab_gsex -> Fill(track_theta[0], vertex_range[0]);
                        h_range_thetalab_gsex -> Fill(track_theta[1], vertex_range[1]);
                        h_kineE_thetalab_gsex -> Fill(track_theta[0], vertex_KinE[0]);
                        h_kineE_thetalab_gsex -> Fill(track_theta[1], vertex_KinE[1]);
                     }
                     else if (!vertex_corr){
                        h_range_thetalab_gsex -> Fill(track_theta[0], track_range[0]);
                        h_range_thetalab_gsex -> Fill(track_theta[1], track_range[1]);
                        h_kineE_thetalab_gsex -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_gsex -> Fill(track_theta[1], track_KinE[1]);
                     }
                     if(tracks_vertex){
                        h_verxy_gsex -> Fill(vtx, vty);
                        h_verz_gsex -> Fill(vtz);
#if defined (vertex_index) || defined (vertex_depth)
                        h_verz_gsex_index[vindex] -> Fill(vtz);
                        h_verxy_gsex_index[vindex] -> Fill(vtx, vty);
                        h_thetalab_thetalab_gsex_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                        if (vertex_corr){
                           // For using new kinematics energy by calculation of new range from vertex to last point of track
                           h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[0], vertex_KinE[0]);
                           h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[1], vertex_KinE[1]);
                        }
                        else if (!vertex_corr){
                           // For using rough kinematics energy by calculation of rough range from first to last point of track
                           h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                           h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
                        }
#endif
                     }
                  }
                  else if (theta_exex->IsInside(track_theta[0],track_theta[1])){
                     if (runNum == 52){
                        peak3.push_back(i);
                     }
                     h_sum_theta_exex -> Fill(sum_theta);
                     if (vertex_corr){
                        h_range_thetalab_exex -> Fill(track_theta[0], vertex_range[0]);
                        h_range_thetalab_exex -> Fill(track_theta[1], vertex_range[1]);
                        h_kineE_thetalab_exex -> Fill(track_theta[0], vertex_KinE[0]);
                        h_kineE_thetalab_exex -> Fill(track_theta[1], vertex_KinE[1]);            
                     }
                     else if (!vertex_corr){
                        h_range_thetalab_exex -> Fill(track_theta[0], track_range[0]);
                        h_range_thetalab_exex -> Fill(track_theta[1], track_range[1]);
                        h_kineE_thetalab_exex -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_exex -> Fill(track_theta[1], track_KinE[1]);            
                     }
                     if(tracks_vertex){
                        h_verxy_exex -> Fill(vtx, vty);
                        h_verz_exex -> Fill(vtz);
#if defined (vertex_index) || defined (vertex_depth)
                        h_verz_exex_index[vindex] -> Fill(vtz);
                        h_verxy_exex_index[vindex] -> Fill(vtx, vty);
                        h_thetalab_thetalab_exex_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                     if (vertex_corr){
                        // For using new kinematics energy by calculation of new range from vertex to last point of track
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[0], vertex_KinE[0]);
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[1], vertex_KinE[1]);
                     }
                     else if (!vertex_corr){
                        // For using rough kinematics energy by calculation of rough range from first to last point of track
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
                     }
#endif
                     }
                  }
               }
            }
         }
         if(i%500==0){
            std::cout << "Reading files: " << nruns << "/" << runNums.size() << "; run: " << runNum << ", events: " << nUnpackEvents 
                        << ". Filling data: " << 100*i/nUnpackEvents << " %!" << std::string(10, ' ') << "\r" << std::flush;
         }
         nevents ++;
      }
      //      std::cout << "  Filled data: 100 %!, " << nUnpackEvents << " events" << std::endl;
      //      std::cout << "Number of 2 tracks events in run" << runNum << ":" << nEventsWith2Tracks << std::endl;
      // Close files.
      unpackFile->Close();
      std::cout << "Reading files: " << nruns << "/" << runNums.size() << "; run: " << runNum << ", events: " << nUnpackEvents 
                  << ". Filling data: 100 %!" << std::string(10, ' ') << "\r" << std::flush;
      nruns ++;
   }

   std::cout << "Finished read data. Starting drawing histograms." << std::string(30, ' ') << std::endl;
   // Draw histograms in TCanvas.
   //  Reset canvas c1
   if (gROOT->FindObject("c1")){
      delete gROOT->FindObject("c1");
   }
   //  set color
   kine_exex_0  -> SetLineColor(kRed);
   angle_gsgs_0 -> SetLineColor(kRed);
   angle_gsgs_0 -> SetLineWidth(1);
   angle_exex_0 -> SetLineColor(kRed);

   f_Ebeam->SetMinimum(0);
   f_Ebeam->SetLineWidth(2);
   f_Ebeam->SetLineColor(kRed);
   sum_kine_beam->SetLineWidth(2);
   sum_kine_beam->SetLineColor(kRed);
   angle_xy->SetLineWidth(1);
   angle_xy->SetLineColor(kRed);

   // scale hist of cross section

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

#ifdef nom_check
   TCanvas *c0 = new TCanvas("c0", "c0");
   c0->cd();
   h_rmax->SetDirectory(0);
   h_rmax->GetXaxis()->SetTitle("Rmax [mm]");
   h_rmax->Draw();

   TCanvas *c1 = new TCanvas("c1", "c1");
   c1->Divide(2,1);
   c1->cd(1);
   h_rmax->SetDirectory(0);
   h_rmax->GetXaxis()->SetTitle("Rmax [mm]");
   h_rmax->Draw();
   c1->cd(2);
   h_rmax_12c->SetDirectory(0);
   h_rmax_12c->GetXaxis()->SetTitle("Rmax [mm]");
   h_rmax_12c->Draw();

   TCanvas *c2 = new TCanvas("c2", "c2");
   c2->cd();
   h_ntra->SetDirectory(0);
   h_ntra->GetXaxis()->SetTitle("Number of tracks");
   gPad->SetLogy();
   h_ntra->Draw();

   TCanvas *c3 = new TCanvas("c3", "c3");
   c3->cd();
   h_charge_range->SetDirectory(0);
   h_charge_range->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range->GetYaxis()->SetTitle("Charge [ADC]");
   h_charge_range->Draw("colz");
   //   cut12c->Draw("same");
   //   cutalpha->Draw("same");
   //   cutp->Draw("same");

   TCanvas *c4 = new TCanvas("c4", "c4");
   c4->cd();
   h_dEdx_range->SetDirectory(0);
   h_dEdx_range->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range->GetYaxis()->SetTitle("dE/dx [ADC/mm]");
   h_dEdx_range->Draw("colz");

   TCanvas *c5 = new TCanvas("c5", "c5");
   c5->cd();
   h_range_thetalab->SetDirectory(0);
   h_range_thetalab->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab->Draw("colz");

   TCanvas *c6 = new TCanvas("c6", "c6");
   c6->cd();
   h_kineE_thetalab->SetDirectory(0);
   h_kineE_thetalab->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab->GetYaxis()->SetTitle("roughKinE [MeV]");
   h_kineE_thetalab->Draw("colz");
   kine_gsgs_0->Draw("same");

   TCanvas *c7 = new TCanvas("c7", "c7");
   c7->cd();
   h_philab_philab->SetDirectory(0);
   h_philab_philab->GetXaxis()->SetTitle("#phi_{LAB} [deg]");
   h_philab_philab->GetYaxis()->SetTitle("#phi_{LAB} [deg]");
   h_philab_philab->Draw("colz");

   TCanvas *c8 = new TCanvas("c8", "c8");
   c8->cd();
   h_charge_range_cutphi->SetDirectory(0);
   h_charge_range_cutphi->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range_cutphi->GetYaxis()->SetTitle("Charge [ADC]");
   h_charge_range_cutphi->SetTitle(Form("Charge Range (phi1-phi2-180 < %d )", (int)del_phi));
   h_charge_range_cutphi->Draw("colz");
   cut12c->Draw("same");
   cutalpha->Draw("same");
   cutproton->Draw("same");

   TCanvas *c9 = new TCanvas("c9", "c9");
   c9->cd();
   h_dEdx_range_cutphi->SetDirectory(0);
   h_dEdx_range_cutphi->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range_cutphi->GetYaxis()->SetTitle("dE/dx [ADC/mm]");
   h_dEdx_range_cutphi->SetTitle(Form("dE/dx Range (phi1-phi2-180 < %d )", (int)del_phi));
   h_dEdx_range_cutphi->Draw("colz");
   //   cut12c_dedx->Draw("same");
   //   cutalpha_dedx->Draw("same");
   //   cutproton_dedx->Draw("same");

   TCanvas *c10 = new TCanvas("c10", "c10");
   c10->cd();
   h_range_thetalab_cutphi->SetDirectory(0);
   h_range_thetalab_cutphi->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cutphi->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cutphi->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d )", (int)del_phi));
   h_range_thetalab_cutphi->Draw("colz");
   //kine_d3He_tt->Draw("same");

   TCanvas *c11 = new TCanvas("c11", "c11");
   c11->cd();
   h_thetalab_thetalab_cutphi->SetDirectory(0);
   h_thetalab_thetalab_cutphi->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d )", (int)del_phi));
   h_thetalab_thetalab_cutphi->Draw("colz");
   angle_gsgs_0->Draw("same");

   TCanvas *c12 = new TCanvas("c12", "c12");
   c12->cd();
   h_philab_philab_cutphi->SetDirectory(0);
   h_philab_philab_cutphi->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_cutphi->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_cutphi->SetTitle(Form("Phi_LAB Phi_LAB (phi1-phi2-180 < %d )", (int)del_phi));
   h_philab_philab_cutphi->Draw("colz");
   //kine_d3He_tt->Draw("same");

   TCanvas *c13 = new TCanvas("c13", "c13");
   c13->cd();
   h_charge_range_cutphi_2tra->SetDirectory(0);
   h_charge_range_cutphi_2tra->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range_cutphi_2tra->GetYaxis()->SetTitle("Charge [ADC]");
   h_charge_range_cutphi_2tra->SetTitle(Form("Charge Range (phi1-phi2-180 < %d, track == 2 )", (int)del_phi));
   h_charge_range_cutphi_2tra->Draw("colz");
   //   cut12c->Draw("same");

   TCanvas *c14 = new TCanvas("c14", "c14");
   c14->cd();
   h_range_thetalab_cutphi_2tra->SetDirectory(0);
   h_range_thetalab_cutphi_2tra->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cutphi_2tra->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cutphi_2tra->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d, track == 2 )", (int)del_phi));
   h_range_thetalab_cutphi_2tra->Draw("colz");

   TCanvas *c15 = new TCanvas("c15", "c15");
   c15->cd();
   h_thetalab_thetalab_cutphi_2tra->SetDirectory(0);
   h_thetalab_thetalab_cutphi_2tra->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi_2tra->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi_2tra->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2 )", (int)del_phi));
   h_thetalab_thetalab_cutphi_2tra->Draw("colz");
   angle_gsgs_0->Draw("same");

   TCanvas *c16 = new TCanvas("c16", "c16");
   c16->cd();
   h_charge_range_cut12c_ela->SetDirectory(0);
   h_charge_range_cut12c_ela->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range_cut12c_ela->GetYaxis()->SetTitle("Charge [ADC]");
   h_charge_range_cut12c_ela->SetTitle(Form("Charge Range (phi1-phi2-180 < %d, track == 2, 12c12c )", (int)del_phi));
   h_charge_range_cut12c_ela->Draw("colz");
   //   cut12c->Draw("same");

   TCanvas *c17 = new TCanvas("c17", "c17");
   c17->cd();
   h_range_thetalab_cut12c_ela->SetDirectory(0);
   h_range_thetalab_cut12c_ela->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cut12c_ela->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cut12c_ela->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c)", (int)del_phi));
   h_range_thetalab_cut12c_ela->Draw("colz");

   TCanvas *c18 = new TCanvas("c18", "c18");
   c18->cd();
   h_kineE_thetalab_carbon->SetDirectory(0);
   h_kineE_thetalab_carbon->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_carbon->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_carbon->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c)", (int)del_phi));
   h_kineE_thetalab_carbon->Draw("colz");
   kine_gsgs_0->Draw("same");

   TCanvas *c19 = new TCanvas("c19", "c19");
   c19->cd();
   h_thetalab_thetalab_cut12c_ela->SetDirectory(0);
   h_thetalab_thetalab_cut12c_ela->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c )", (int)del_phi));
   h_thetalab_thetalab_cut12c_ela->Draw("colz");
   angle_gsgs_0->Draw("same");
   angle_gsex_0->Draw("same");
   angle_exex_0->Draw("same");
   theta_gsgs->Draw("same");
   theta_gsex->Draw("same");
   theta_exex->Draw("same");

   TCanvas *c20 = new TCanvas("c20", "c20");
   c20->cd();
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();

   TCanvas *c21 = new TCanvas("c21", "c21", 600, 1000);
   c21->Divide(2,4);
   c21->cd(1);
   h_ntraver_ntra->SetDirectory(0);
   h_ntraver_ntra->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra->SetTitle(Form("tracks with vertex"));
   gPad->SetLogz();
   h_ntraver_ntra->Draw("colz");
   c21->cd(2);
   h_ntra_verz->SetDirectory(0);
   h_ntra_verz->GetXaxis()->SetTitle("vertex z [mm]");
   h_ntra_verz->GetYaxis()->SetTitle("number of tracks");
   h_ntra_verz->SetTitle(Form("vertex z vs number of tracks"));
   gPad->SetLogz();
   h_ntra_verz->Draw("colz");
   c21->cd(3);
   h_verxy->SetDirectory(0);
   h_verxy->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy->SetTitle(Form("vertex xy"));
   gPad->SetLogz();
   h_verxy->Draw("colz");
   c21->cd(4);
   h_verz->SetDirectory(0);
   h_verz->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz->SetTitle(Form("vertex z"));
   gPad->SetLogy();
   h_verz->Draw();
   c21->cd(5);
   h_verxz->SetDirectory(0);
   h_verxz->GetXaxis()->SetTitle("vertex z [mm]");
   h_verxz->GetYaxis()->SetTitle("vertex x [mm]");
   h_verxz->SetTitle(Form("vertex zx"));
   gPad->SetLogz();
   h_verxz->Draw("colz");
   c21->cd(6);
   h_veryz->SetDirectory(0);
   h_veryz->GetXaxis()->SetTitle("vertex z [mm]");
   h_veryz->GetYaxis()->SetTitle("vertex y [mm]");
   h_veryz->SetTitle(Form("vertex zy"));
   gPad->SetLogz();
   h_veryz->Draw("colz");
   c21->cd(7);
   h_verxy_cut12c_ela->SetDirectory(0);
   h_verxy_cut12c_ela->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_cut12c_ela->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_cut12c_ela->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verxy_cut12c_ela->Draw("colz");
   c21->cd(8);
   h_verz_cut12c_ela->SetDirectory(0);
   h_verz_cut12c_ela->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_cut12c_ela->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verz_cut12c_ela->Draw();
   //   c21->SaveAs("can_output/check_vd_76matm_c21_vertex.pdf");

#endif

#ifdef peak_check
   TCanvas *c30 = new TCanvas("c30", "c30");
   c30->Divide(2,1);
   c30->cd(1);
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->SetTitle(Form("Sum Theta_LAB (track == 2, 12c12c)"));
   h_sum_theta_cut12c_ela->Draw();
   c30->cd(2);
   h_sum_theta_cut12c_run52->SetDirectory(0);
   h_sum_theta_cut12c_run52->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_run52->SetTitle(Form("Sum Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c, run 52)", (int)del_phi));
   h_sum_theta_cut12c_run52->Draw();

   TCanvas *c31 = new TCanvas("c31", "c31", 1200,600);
   c31->Divide(4,1);
   c31->cd(1);
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();
   c31->cd(2);
   h_sum_theta_gsgs->SetDirectory(0);
   h_sum_theta_gsgs->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c gsgs");
   h_sum_theta_gsgs->Draw();
   c31->cd(3);
   h_sum_theta_gsex->SetDirectory(0);
   h_sum_theta_gsex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (gsex)");
   h_sum_theta_gsex->Draw();
   c31->cd(4);
   h_sum_theta_exex->SetDirectory(0);
   h_sum_theta_exex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (exex)");
   h_sum_theta_exex->Draw();

   TCanvas *c32 = new TCanvas("c32", "c32", 600, 600);
   c32->cd();
   h_sum_theta_gsgs->SetDirectory(0);
   h_sum_theta_gsgs->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c gsgs");
   h_sum_theta_gsgs->GetXaxis()->SetRangeUser(85, 100);
   h_sum_theta_gsgs->Draw();

#endif

#ifdef c12_check
   TCanvas *c40 = new TCanvas("c40", "c40",1200,1000);
   h_thetalab_thetalab_cut12c_ela->SetDirectory(0);
   h_thetalab_thetalab_cut12c_ela->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c )", (int)del_phi));
   gPad->SetLogz();
   h_thetalab_thetalab_cut12c_ela->SetMinimum(1);
   h_thetalab_thetalab_cut12c_ela->Draw("colz");
   angle_gsgs_0->Draw("same");
   angle_gsex_0->Draw("same");
   angle_exex_0->Draw("same");
   theta_gsgs->Draw("same");
   theta_gsex->Draw("same");
   theta_exex->Draw("same");

   TCanvas *c41 = new TCanvas("c41", "c41", 1200, 1000);
   c41->Divide(4,4);
   c41->cd(1);
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();
   c41->cd(2);
   h_sum_theta_gsgs->SetDirectory(0);
   h_sum_theta_gsgs->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c gsgs");
   h_sum_theta_gsgs->Draw();
   c41->cd(3);
   h_sum_theta_gsex->SetDirectory(0);
   h_sum_theta_gsex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (gsex)");
   h_sum_theta_gsex->Draw();
   c41->cd(4);
   h_sum_theta_exex->SetDirectory(0);
   h_sum_theta_exex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (exex)");
   h_sum_theta_exex->Draw();
   c41->cd(5);
   h_kineE_thetalab_carbon->SetDirectory(0);
   h_kineE_thetalab_carbon->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_carbon->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_carbon->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_carbon->SetMinimum(1);
   h_kineE_thetalab_carbon->Draw("colz");
   kine_gsgs_0->Draw("same");
   kine_gsex_0->Draw("same");
   kine_exex_0->Draw("same");
   c41->cd(6);
   h_kineE_thetalab_gsgs->SetDirectory(0);
   h_kineE_thetalab_gsgs->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_gsgs->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_gsgs->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, gsgs)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_gsgs->SetMinimum(1);
   h_kineE_thetalab_gsgs->Draw("colz");
   kine_gsgs_0->Draw("same");
   kine_gsex_0->Draw("same");
   kine_exex_0->Draw("same");
   c41->cd(7);
   h_kineE_thetalab_gsex->SetDirectory(0);
   h_kineE_thetalab_gsex->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_gsex->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_gsex->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, gsex)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_gsex->SetMinimum(1);
   h_kineE_thetalab_gsex->Draw("colz");
   kine_gsgs_0->Draw("same");
   kine_gsex_0->Draw("same");
   kine_exex_0->Draw("same");
   c41->cd(8);
   h_kineE_thetalab_exex->SetDirectory(0);
   h_kineE_thetalab_exex->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_exex->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_exex->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, exex)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_exex->SetMinimum(1);
   h_kineE_thetalab_exex->Draw("colz");
   kine_gsgs_0->Draw("same");
   kine_gsex_0->Draw("same");
   kine_exex_0->Draw("same");
   c41->cd(9);
   h_verxy_cut12c_ela->SetDirectory(0);
   h_verxy_cut12c_ela->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_cut12c_ela->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_cut12c_ela->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verxy_cut12c_ela->Draw("colz");
   c41->cd(10);
   h_verxy_gsgs->SetDirectory(0);
   h_verxy_gsgs->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_gsgs->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_gsgs->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   h_verxy_gsgs->Draw("colz");
   c41->cd(11);
   h_verxy_gsex->SetDirectory(0);
   h_verxy_gsex->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_gsex->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_gsex->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, gsex)", (int)del_phi));
   h_verxy_gsex->Draw("colz");
   c41->cd(12);
   h_verxy_exex->SetDirectory(0);
   h_verxy_exex->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_exex->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_exex->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, exex)", (int)del_phi));
   h_verxy_exex->Draw("colz");
   c41->cd(13);
   h_verz_cut12c_ela->SetDirectory(0);
   h_verz_cut12c_ela->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_cut12c_ela->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verz_cut12c_ela->Draw();
   c41->cd(14);
   h_verz_gsgs->SetDirectory(0);
   h_verz_gsgs->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   h_verz_gsgs->Draw();
   c41->cd(15);
   h_verz_gsex->SetDirectory(0);
   h_verz_gsex->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsex->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, gsex)", (int)del_phi));
   h_verz_gsex->Draw();
   c41->cd(16);
   h_verz_exex->SetDirectory(0);
   h_verz_exex->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_exex->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, exex)", (int)del_phi));
   h_verz_exex->Draw();
   //   c41->SaveAs("can_output/check_vd_76matm_c41_vertex_states.pdf");

   TCanvas *c42 = new TCanvas("c42", "c42", 600, 600);
   c42->Divide(2,2);
   c42->cd(1);
   h_ntraver_ntra->SetDirectory(0);
   h_ntraver_ntra->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra->SetTitle(Form("tracks with vertex"));
   gPad->SetLogz();
   h_ntraver_ntra->Draw("colz");
   c42->cd(2);
   h_ntraver_ntra_cutphi->SetDirectory(0);
   h_ntraver_ntra_cutphi->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra_cutphi->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra_cutphi->SetTitle(Form("tracks with vertex (phi1-phi2-180 < %d)", (int)del_phi));
   gPad->SetLogz();
   h_ntraver_ntra_cutphi->Draw("colz");
   c42->cd(3);
   h_ntraver_ntra_cutphi_2tra->SetDirectory(0);
   h_ntraver_ntra_cutphi_2tra->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra_cutphi_2tra->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra_cutphi_2tra->SetTitle(Form("tracks with vertex (phi1-phi2-180 < %d, track == 2)", (int)del_phi));
   gPad->SetLogz();
   h_ntraver_ntra_cutphi_2tra->Draw("colz");
   c42->cd(4);
   h_ntraver_ntra_cut12c_ela->SetDirectory(0);
   h_ntraver_ntra_cut12c_ela->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra_cut12c_ela->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra_cut12c_ela->SetTitle(Form("tracks with vertex (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   gPad->SetLogz();
   h_ntraver_ntra_cut12c_ela->Draw("colz");
   //   c50->SaveAs("can_output/check_vd_76matm_c50_vertex_tracks.pdf");

#endif

#ifdef gsgs_check
   TCanvas *c51 = new TCanvas("c51", "c51", 1600, 1000);
   c51->Divide(6, 3);
   c51->cd(1);
   h_thetalab_thetalab_gsgs->SetDirectory(0);
   h_thetalab_thetalab_gsgs->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs->SetMinimum(1);
   h_thetalab_thetalab_gsgs->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c51->cd(2);
   h_thetalab_thetalab_gsgs_cm50->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm50->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm50->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm50->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 23<#theta_1<27 or 23<#theta_2<27)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm50->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm50->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c51->cd(3);
   h_thetalab_thetalab_gsgs_cm60->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm60->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm60->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm60->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 28<#theta_1<32 or 28<#theta_2<32)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm60->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm60->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c51->cd(4);
   h_thetalab_thetalab_gsgs_cm70->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm70->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm70->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm70->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 33<#theta_1<37 or 33<#theta_2<37)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm70->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm70->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c51->cd(5);
   h_thetalab_thetalab_gsgs_cm80->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm80->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm80->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm80->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 38<#theta_1<42 or 38<#theta_2<42)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm80->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm80->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c51->cd(6);
   h_thetalab_thetalab_gsgs_cm90->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm90->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm90->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm90->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm90->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm90->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");

   c51->cd(7);
   h_verz_gsgs->SetDirectory(0);
   h_verz_gsgs->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs->SetTitle(Form("vertex z (12c12c, gsgs)"));
   gPad->SetLogy();
   h_verz_gsgs->Draw();
   c51->cd(8);
   h_verz_gsgs_cm50->SetDirectory(0);
   h_verz_gsgs_cm50->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs_cm50->SetTitle(Form("vertex z (12c12c, gsgs, 23<#theta_1<27 or 23<#theta_2<27)"));
   gPad->SetLogy();
   h_verz_gsgs_cm50->Draw();
   c51->cd(9);
   h_verz_gsgs_cm60->SetDirectory(0);
   h_verz_gsgs_cm60->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs_cm60->SetTitle(Form("vertex z (12c12c, gsgs, 28<#theta_1<32 or 28<#theta_2<32)"));
   gPad->SetLogy();
   h_verz_gsgs_cm60->Draw();
   c51->cd(10);
   h_verz_gsgs_cm70->SetDirectory(0);
   h_verz_gsgs_cm70->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs_cm70->SetTitle(Form("vertex z (12c12c, gsgs, 33<#theta_1<37 or 33<#theta_2<37)"));
   gPad->SetLogy();
   h_verz_gsgs_cm70->Draw();
   c51->cd(11);
   h_verz_gsgs_cm80->SetDirectory(0);
   h_verz_gsgs_cm80->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs_cm80->SetTitle(Form("vertex z (12c12c, gsgs, 38<#theta_1<42 or 38<#theta_2<42)"));
   gPad->SetLogy();
   h_verz_gsgs_cm80->Draw();
   c51->cd(12);
   h_verz_gsgs_cm90->SetDirectory(0);
   h_verz_gsgs_cm90->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs_cm90->SetTitle(Form("vertex z (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogy();
   h_verz_gsgs_cm90->Draw();

   c51->cd(13);
   h_Ebeam_gsgs->SetDirectory(0);
   h_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs)"));
   gPad->SetLogy();
   h_Ebeam_gsgs->Draw();
   c51->cd(14);
   h_Ebeam_gsgs_cm50->SetDirectory(0);
   h_Ebeam_gsgs_cm50->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm50->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 23<#theta_1<27 or 23<#theta_2<27)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm50->Draw();
   c51->cd(15);
   h_Ebeam_gsgs_cm60->SetDirectory(0);
   h_Ebeam_gsgs_cm60->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm60->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 28<#theta_1<32 or 28<#theta_2<32)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm60->Draw();
   c51->cd(16);
   h_Ebeam_gsgs_cm70->SetDirectory(0);
   h_Ebeam_gsgs_cm70->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm70->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 33<#theta_1<37 or 33<#theta_2<37)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm70->Draw();
   c51->cd(17);
   h_Ebeam_gsgs_cm80->SetDirectory(0);
   h_Ebeam_gsgs_cm80->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm80->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 38<#theta_1<42 or 38<#theta_2<42)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm80->Draw();
   c51->cd(18);
   h_Ebeam_gsgs_cm90->SetDirectory(0);
   h_Ebeam_gsgs_cm90->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm90->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm90->Draw();

   TCanvas *c52 = new TCanvas("c52", "c52", 1600, 600);
   c52->Divide(6, 2);
   c52->cd(1);
   h_thetalab_thetalab_gsgs->SetDirectory(0);
   h_thetalab_thetalab_gsgs->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs->SetMinimum(1);
   h_thetalab_thetalab_gsgs->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c52->cd(2);
   h_thetalab_thetalab_gsgs_cm50->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm50->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm50->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm50->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 23<#theta_1<27 or 23<#theta_2<27)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm50->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm50->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c52->cd(3);
   h_thetalab_thetalab_gsgs_cm60->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm60->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm60->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm60->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 28<#theta_1<32 or 28<#theta_2<32)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm60->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm60->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c52->cd(4);
   h_thetalab_thetalab_gsgs_cm70->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm70->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm70->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm70->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 33<#theta_1<37 or 33<#theta_2<37)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm70->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm70->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c52->cd(5);
   h_thetalab_thetalab_gsgs_cm80->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm80->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm80->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm80->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 38<#theta_1<42 or 38<#theta_2<42)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm80->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm80->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c52->cd(6);
   h_thetalab_thetalab_gsgs_cm90->SetDirectory(0);
   h_thetalab_thetalab_gsgs_cm90->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm90->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs_cm90->SetTitle(Form("Theta_LAB Theta_LAB (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs_cm90->SetMinimum(1);
   h_thetalab_thetalab_gsgs_cm90->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");

   c52->cd(7);
   h_Ebcm_gsgs->SetDirectory(0);
   h_Ebcm_gsgs->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs)"));
   gPad->SetLogy();
   h_Ebcm_gsgs->Draw();
   c52->cd(8);
   h_Ebcm_gsgs_cm50->SetDirectory(0);
   h_Ebcm_gsgs_cm50->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm50->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 23<#theta_1<27 or 23<#theta_2<27)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm50->Draw();
   c52->cd(9);
   h_Ebcm_gsgs_cm60->SetDirectory(0);
   h_Ebcm_gsgs_cm60->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm60->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 28<#theta_1<32 or 28<#theta_2<32)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm60->Draw();
   c52->cd(10);
   h_Ebcm_gsgs_cm70->SetDirectory(0);
   h_Ebcm_gsgs_cm70->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm70->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 33<#theta_1<37 or 33<#theta_2<37)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm70->Draw();
   c52->cd(11);
   h_Ebcm_gsgs_cm80->SetDirectory(0);
   h_Ebcm_gsgs_cm80->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm80->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 38<#theta_1<42 or 38<#theta_2<42)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm80->Draw();
   c52->cd(12);
   h_Ebcm_gsgs_cm90->SetDirectory(0);
   h_Ebcm_gsgs_cm90->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
   h_Ebcm_gsgs_cm90->SetTitle(Form("Estimated E_{beam.cm} (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogy();
   h_Ebcm_gsgs_cm90->Draw();

   TCanvas *c53 = new TCanvas("c53", "c53");
   c53->cd();
   h_Ebcm_gsgs_cm90->SetDirectory(0);
   h_Ebcm_gsgs_cm90->SetLineColor(kBlue);
   h_Ebcm_gsgs_cm90->SetLineWidth(2);
   h_Ebcm_gsgs_cm90->SetMinimum(0.5);
   h_Ebcm_gsgs_cm90->SetMaximum(50);
   h_Ebcm_gsgs_cm90->GetXaxis()->SetRangeUser(0, 35);
   h_Ebcm_gsgs_cm90->GetXaxis()->SetTitle("E_{beam.cm} [MeV]");
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
   c53->Update();

   TCanvas *c54 = new TCanvas("c54", "c54", 1600, 600);
   c54->Divide(3,1);
   c54->cd(1);
   h_Ebeam_gsgs->SetDirectory(0);
   h_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs)"));
   gPad->SetLogy();
   h_Ebeam_gsgs->Draw();
   c54->cd(2);
   h_sumkine_Ebeam_gsgs->SetDirectory(0);
   h_sumkine_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_Ebeam_gsgs->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_sumkine_Ebeam_gsgs->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (12c12c, gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_Ebeam_gsgs->SetMinimum(1);
   h_sumkine_Ebeam_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");
   c54->cd(3);
   h_sumkine_kineE_gsgs->SetDirectory(0);
   h_sumkine_kineE_gsgs->GetXaxis()->SetTitle("E_{kine} [MeV]");
   h_sumkine_kineE_gsgs->GetYaxis()->SetTitle("sum E_{kine} [MeV]");
   h_sumkine_kineE_gsgs->SetTitle(Form("Sum E_{kine} vs E_{kine} , (12c12c, gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_kineE_gsgs->SetMinimum(1);
   h_sumkine_kineE_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");

   TCanvas *c55 = new TCanvas("c55", "c55", 1000, 1000);
   c55->Divide(2,2);
   c55->cd(1);
   h_Ebeam_gsgs_cm90->SetDirectory(0);
   h_Ebeam_gsgs_cm90->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_Ebeam_gsgs_cm90->SetTitle(Form("Estimated E_{beam} (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogy();
   h_Ebeam_gsgs_cm90->Draw();
   c55->cd(2);
   h_sumkine_Ebeam_gsgs_cm90->SetDirectory(0);
   h_sumkine_Ebeam_gsgs_cm90->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_Ebeam_gsgs_cm90->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_sumkine_Ebeam_gsgs_cm90->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogz();
   //   h_sumkine_Ebeam_gsgs_cm90->SetMaximum(100);
   h_sumkine_Ebeam_gsgs_cm90->SetMinimum(1);
   h_sumkine_Ebeam_gsgs_cm90->Draw("colz");
   sum_kine_beam->Draw("same");
   c55->cd(3);
   h_sumkine_kineE_gsgs_cm90->SetDirectory(0);
   h_sumkine_kineE_gsgs_cm90->GetXaxis()->SetTitle("E_{kine} [MeV]");
   h_sumkine_kineE_gsgs_cm90->GetYaxis()->SetTitle("sum E_{kine} [MeV]");
   h_sumkine_kineE_gsgs_cm90->SetTitle(Form("Sum E_{kine} vs E_{kine} , (12c12c, gsgs, 43<#theta_1<47 or 43<#theta_2<47)"));
   gPad->SetLogz();
   //   h_sumkine_Ebeam_gsgs_cm90->SetMaximum(100);
   h_sumkine_kineE_gsgs_cm90->SetMinimum(1);
   h_sumkine_kineE_gsgs_cm90->Draw("colz");
   sum_kine_beam->Draw("same");

#endif

#ifdef kine_comp
   TCanvas *c60 = new TCanvas("c60", "c60", 1600, 1000);
   c60->Divide(3, 2);
   c60->cd(1);
   h_lastpoint_comp->SetDirectory(0);
   h_lastpoint_comp->GetXaxis()->SetTitle("first term [mm]");
   h_lastpoint_comp->GetYaxis()->SetTitle("second term [mm]");
   h_lastpoint_comp->SetTitle(Form("lastpoint [rough vs vertex]"));
   gPad->SetLogz();
   //   h_lastpoint_comp->SetMaximum(100);
   h_lastpoint_comp->SetMinimum(1);
   h_lastpoint_comp->Draw("colz");
   c60->cd(2);
   h_vertex_comp->SetDirectory(0);
   h_vertex_comp->GetXaxis()->SetTitle("first term [mm]");
   h_vertex_comp->GetYaxis()->SetTitle("second term [mm]");
   h_vertex_comp->SetTitle(Form("vertex [rough vs vertex]"));
   gPad->SetLogz();
   //   h_lastpoint_comp->SetMaximum(100);
   h_vertex_comp->SetMinimum(1);
   h_vertex_comp->Draw("colz");
   c60->cd(3);
   h_distance_comp->SetDirectory(0);
   h_distance_comp->GetXaxis()->SetTitle("rough distance [mm]");
   h_distance_comp->GetYaxis()->SetTitle("vertex distance [mm]");
   h_distance_comp->SetTitle(Form("distance between 0 and lastpoint [rough vs vertex]"));
   gPad->SetLogz();
   //   h_distance_comp->SetMaximum(100);
   h_distance_comp->SetMinimum(1);
   h_distance_comp->Draw("colz");
   c60->cd(4);
   h_range_comp->SetDirectory(0);
   h_range_comp->GetXaxis()->SetTitle("roughRange [mm]");
   h_range_comp->GetYaxis()->SetTitle("vertexRange [mm]");
   h_range_comp->SetTitle(Form("Range [rough vs vertex]"));
   gPad->SetLogz();
   //   h_range_comp->SetMaximum(100);
   h_range_comp->SetMinimum(1);
   h_range_comp->Draw("colz");
   c60->cd(5);
   h_thetalab_comp->SetDirectory(0);
   h_thetalab_comp->GetXaxis()->SetTitle("#theta_{rough} [deg]");
   h_thetalab_comp->GetYaxis()->SetTitle("#theta_{vertex} [deg]");
   h_thetalab_comp->SetTitle(Form("#theta angle [rough vs vertex]"));
   gPad->SetLogz();
   //   h_range_comp->SetMaximum(100);
   h_thetalab_comp->SetMinimum(1);
   h_thetalab_comp->Draw("colz");
   c60->cd(6);
   h_philab_comp->SetDirectory(0);
   h_philab_comp->GetXaxis()->SetTitle("#phi_{rough} [deg]");
   h_philab_comp->GetYaxis()->SetTitle("#phi_{vertex} [deg]");
   h_philab_comp->SetTitle(Form("#phi angle [rough vs vertex]"));
   gPad->SetLogz();
   //   h_kineE_comp->SetMaximum(100);
   h_philab_comp->SetMinimum(1);
   h_philab_comp->Draw("colz");

   TCanvas *c61 = new TCanvas("c61", "c61", 1000, 1000);
   c61->Divide(2, 2);
   c61->cd(1);
   h_range_comp_gsgs->SetDirectory(0);
   h_range_comp_gsgs->GetXaxis()->SetTitle("roughRange [mm]");
   h_range_comp_gsgs->GetYaxis()->SetTitle("vertexRange [mm]");
   h_range_comp_gsgs->SetTitle(Form("Range [rough vs vertex] (gsgs)"));
   gPad->SetLogz();
   //   h_range_comp->SetMaximum(100);
   h_range_comp_gsgs->SetMinimum(1);
   h_range_comp_gsgs->Draw("colz");
   c61->cd(2);
   h_kineE_comp_gsgs->SetDirectory(0);
   h_kineE_comp_gsgs->GetXaxis()->SetTitle("roughkineE [mm]");
   h_kineE_comp_gsgs->GetYaxis()->SetTitle("vertexkineE [mm]");
   h_kineE_comp_gsgs->SetTitle(Form("kineE [rough vs vertex] (gsgs)"));
   gPad->SetLogz();
   //   h_kineE_comp->SetMaximum(100);
   h_kineE_comp_gsgs->SetMinimum(1);
   h_kineE_comp_gsgs->Draw("colz");
   c61->cd(3);
   h_thetalab_comp_gsgs->SetDirectory(0);
   h_thetalab_comp_gsgs->GetXaxis()->SetTitle("#theta_{rough} [deg]");
   h_thetalab_comp_gsgs->GetYaxis()->SetTitle("#theta_{vertex} [deg]");
   h_thetalab_comp_gsgs->SetTitle(Form("#theta angle [rough vs vertex] (gsgs)"));
   gPad->SetLogz();
   //   h_range_comp->SetMaximum(100);
   h_thetalab_comp_gsgs->SetMinimum(1);
   h_thetalab_comp_gsgs->Draw("colz");
   angle_xy->Draw("same");
   c61->cd(4);
   h_philab_comp_gsgs->SetDirectory(0);
   h_philab_comp_gsgs->GetXaxis()->SetTitle("#phi_{rough} [deg]");
   h_philab_comp_gsgs->GetYaxis()->SetTitle("#phi_{vertex} [deg]");
   h_philab_comp_gsgs->SetTitle(Form("#phi angle [rough vs vertex] (gsgs)"));
   gPad->SetLogz();
   //   h_kineE_comp->SetMaximum(100);
   h_philab_comp_gsgs->SetMinimum(1);
   h_philab_comp_gsgs->Draw("colz");
   angle_xy->Draw("same");

   TCanvas *c62 = new TCanvas("c62", "c62", 1000, 1000);
   c62->Divide(2, 2);
   c62->cd(1);
   h_thetalab_thetalab_gsgs->SetDirectory(0);
   h_thetalab_thetalab_gsgs->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_gsgs->SetTitle(Form("Theta_LAB Theta_LAB (gsgs)"));
   gPad->SetLogz();
   h_thetalab_thetalab_gsgs->SetMinimum(1);
   h_thetalab_thetalab_gsgs->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c62->cd(2);
   h_thetalab_thetalab_ver_gsgs->SetDirectory(0);
   h_thetalab_thetalab_ver_gsgs->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_ver_gsgs->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_ver_gsgs->SetTitle(Form("Theta_LAB Theta_LAB recalculated (gsgs)"));
   gPad->SetLogz();
   h_thetalab_thetalab_ver_gsgs->SetMinimum(1);
   h_thetalab_thetalab_ver_gsgs->Draw("colz");
   angle_gsgs_0->Draw("same");
   theta_gsgs->Draw("same");
   c62->cd(3);
   h_philab_philab_gsgs->SetDirectory(0);
   h_philab_philab_gsgs->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_gsgs->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_gsgs->SetTitle(Form("Phi_LAB Phi_LAB (gsgs)"));
   h_philab_philab_gsgs->Draw("colz");
   c62->cd(4);
   h_philab_philab_ver_gsgs->SetDirectory(0);
   h_philab_philab_ver_gsgs->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_ver_gsgs->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_ver_gsgs->SetTitle(Form("Phi_LAB Phi_LAB recalculated (gsgs)"));
   h_philab_philab_ver_gsgs->Draw("colz");

   TCanvas *c63 = new TCanvas("c63", "c63", 1000, 1000);
   c63->Divide(2,2);
   c63->cd(1);
   h_sumkine_Ebeam_gsgs->SetDirectory(0);
   h_sumkine_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_Ebeam_gsgs->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_sumkine_Ebeam_gsgs->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_Ebeam_gsgs->SetMinimum(1);
   h_sumkine_Ebeam_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");
   c63->cd(2);
   h_sumkine_ver_Ebeam_gsgs->SetDirectory(0);
   h_sumkine_ver_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_ver_Ebeam_gsgs->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_ver_Ebeam_gsgs->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_ver_Ebeam_gsgs->SetMinimum(1);
   h_sumkine_ver_Ebeam_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");
   c63->cd(3);
   h_sumkine_kineE_gsgs->SetDirectory(0);
   h_sumkine_kineE_gsgs->GetXaxis()->SetTitle("E_{kine} [MeV]");
   h_sumkine_kineE_gsgs->GetYaxis()->SetTitle("sum E_{kine} [MeV]");
   h_sumkine_kineE_gsgs->SetTitle(Form("Sum E_{kine} vs E_{kine} , (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_kineE_gsgs->SetMaximum(100);
   h_sumkine_kineE_gsgs->SetMinimum(1);
   h_sumkine_kineE_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");
   c63->cd(4);
   h_sumkine_kineE_ver_gsgs->SetDirectory(0);
   h_sumkine_kineE_ver_gsgs->GetXaxis()->SetTitle("E_{kine} [MeV]");
   h_sumkine_kineE_ver_gsgs->GetYaxis()->SetTitle("sum E_{kine} [MeV]");
   h_sumkine_kineE_ver_gsgs->SetTitle(Form("Sum E_{kine} vs E_{kine} , (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_kineE_ver_gsgs->SetMaximum(100);
   h_sumkine_kineE_ver_gsgs->SetMinimum(1);
   h_sumkine_kineE_ver_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");

   TCanvas *c64 = new TCanvas("c64", "c64", 1200, 1200);
   c64->Divide(2,2);
   c64->cd(1);
   h_sumkine_verz_cut12c_ela->SetDirectory(0);
   h_sumkine_verz_cut12c_ela->GetXaxis()->SetTitle("vertex z [mm]");
   h_sumkine_verz_cut12c_ela->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_verz_cut12c_ela->SetTitle(Form("Sum KineE vs vertex z, (12c12c)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_verz_cut12c_ela->SetMinimum(1);
   h_sumkine_verz_cut12c_ela->Draw("colz");
   f_Ebeam->Draw("same");
   c64->cd(2);
   h_sumkine_ver_Ebeam_cut12c_ela->SetDirectory(0);
   h_sumkine_ver_Ebeam_cut12c_ela->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_ver_Ebeam_cut12c_ela->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_ver_Ebeam_cut12c_ela->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (12c12c)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_ver_Ebeam_cut12c_ela->SetMinimum(1);
   h_sumkine_ver_Ebeam_cut12c_ela->Draw("colz");
   sum_kine_beam->Draw("same");
   c64->cd(3);
   h_sumkine_verz_gsgs->SetDirectory(0);
   h_sumkine_verz_gsgs->GetXaxis()->SetTitle("vertex z [mm]");
   h_sumkine_verz_gsgs->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_verz_gsgs->SetTitle(Form("Sum KineE vs vertex z, (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_verz_gsgs->SetMinimum(1);
   h_sumkine_verz_gsgs->Draw("colz");
   f_Ebeam->Draw("same");
   c64->cd(4);
   h_sumkine_ver_Ebeam_gsgs->SetDirectory(0);
   h_sumkine_ver_Ebeam_gsgs->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_ver_Ebeam_gsgs->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_ver_Ebeam_gsgs->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (gsgs)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_ver_Ebeam_gsgs->SetMinimum(1);
   h_sumkine_ver_Ebeam_gsgs->Draw("colz");
   sum_kine_beam->Draw("same");

#endif

   TCanvas *c65 = new TCanvas("c65", "c65", 1000, 1000);
   c65->cd();
   h_sumkine_ver_Ebeam_cut12c_ela->SetDirectory(0);
   h_sumkine_ver_Ebeam_cut12c_ela->GetXaxis()->SetTitle("E_{beam} [MeV]");
   h_sumkine_ver_Ebeam_cut12c_ela->GetYaxis()->SetTitle("vertexKineE [MeV]");
   h_sumkine_ver_Ebeam_cut12c_ela->SetTitle(Form("Sum KineE vs Estimated E_{beam} , (12c12c)"));
   gPad->SetLogz();
   //   h_sumkine_ver_Ebeam_gsgs->SetMaximum(100);
   h_sumkine_ver_Ebeam_cut12c_ela->SetMinimum(1);
   h_sumkine_ver_Ebeam_cut12c_ela->Draw("colz");
   sum_kine_beam->Draw("same");

#ifdef vertex_index
   draw_ind("c70", "gsgs", n_group, n_h_z, h_verz_gsgs, h_verxy_gsgs, h_kineE_thetalab_gsgs, h_thetalab_thetalab_cut12c_ela,
      h_verz_gsgs_index, h_verxy_gsgs_index, h_kineE_thetalab_gsgs_index, h_thetalab_thetalab_gsgs_index,
      kine_gsgs_0, kine_gsex_0, kine_exex_0, angle_gsgs_0, angle_gsex_0, angle_exex_0
   );

   draw_ind("c71", "gsex", n_group, n_h_z, h_verz_gsex, h_verxy_gsex, h_kineE_thetalab_gsex, h_thetalab_thetalab_cut12c_ela,
      h_verz_gsex_index, h_verxy_gsex_index, h_kineE_thetalab_gsex_index, h_thetalab_thetalab_gsex_index,
      kine_gsgs_0, kine_gsex_0, kine_exex_0, angle_gsgs_0, angle_gsex_0, angle_exex_0
   );

   draw_ind("c72", "exex", n_group, n_h_z, h_verz_exex, h_verxy_exex, h_kineE_thetalab_exex, h_thetalab_thetalab_cut12c_ela,
      h_verz_exex_index, h_verxy_exex_index, h_kineE_thetalab_exex_index, h_thetalab_thetalab_exex_index,
      kine_gsgs_0, kine_gsex_0, kine_exex_0, angle_gsgs_0, angle_gsex_0, angle_exex_0
   );

#endif

#ifdef vertex_depth
   draw_dep("c80", "gsgs", "full", n_group, n_div, h_kineE_thetalab_gsgs_index, kine_gsgs_i);
   //   draw_dep("c81", "gsgs", "center", n_group, n_div, h_kineE_thetalab_gsgs_index, kine_gsgs_i);
   draw_dep("c82", "gsex", "full", n_group, n_div, h_kineE_thetalab_gsex_index, kine_gsex_i);
   //   draw_dep("c83", "gsex", "center", n_group, n_div, h_kineE_thetalab_gsex_index, kine_gsex_i);
   draw_dep("c84", "exex", "full", n_group, n_div, h_kineE_thetalab_exex_index, kine_exex_i);
   //   draw_dep("c85", "exex", "center", n_group, n_div, h_kineE_thetalab_exex_index, kine_exex_i);

#endif

#ifdef check_verz0
   TCanvas *c100 = new TCanvas("c100", "c100");
   h_verz01->SetDirectory(0);
   h_verz01->Draw();
   h_verz01->GetXaxis()->SetTitle("vertex z [mm]");
#endif

   // Saving histograms in a .root file ...
   Results->cd();
   h_rmax->Write();
   h_rmax_12c->Write();
   h_ntra->Write();
   h_sum_theta_cut12c_ela->Write();
   h_sum_theta_cut12c_run52->Write();
   h_sum_theta_gsgs->Write();
   h_sum_theta_gsex->Write();
   h_sum_theta_exex->Write();

   // charge range
   h_charge_range->Write();
   h_charge_range_cutphi->Write();
   h_charge_range_cutphi_2tra->Write();
   h_charge_range_cut12c_ela->Write();

   // dEdx Vs Total Range
   h_dEdx_range->Write();
   //   h_dEdx_range_backwards->Write();
   h_dEdx_range_cutphi->Write();
   h_dEdx_range_cutphi_12c->Write();

   // Kinematics
   h_kineE_thetalab->Write();
   h_kineE_thetalab_carbon->Write();
   h_kineE_thetalab_gsgs->Write();
   h_kineE_thetalab_gsex->Write();
   h_kineE_thetalab_exex->Write();
   h_sumkine_kineE_gsgs->Write();
   h_sumkine_kineE_ver_gsgs->Write();
   h_sumkine_kineE_gsgs_cm90->Write();
   h_sumkine_verz_gsgs->Write();

   // angle correlations
   //  range vs theta
   h_range_thetalab->Write();
   h_range_thetalab_cutphi->Write();
   h_range_thetalab_cutphi_2tra->Write();
   h_range_thetalab_cut12c_ela->Write();
   h_range_thetalab_gsgs->Write();
   h_range_thetalab_gsex->Write();
   h_range_thetalab_exex->Write();
   //  theta vs theta
   h_thetalab_thetalab_cutphi->Write();
   h_thetalab_thetalab_cutphi_2tra->Write();
   h_thetalab_thetalab_cut12c_ela->Write();
   h_thetalab_thetalab_gsgs->Write();
   h_thetalab_thetalab_ver_gsgs->Write();
   h_thetalab_thetalab_gsgs_cm50->Write();
   h_thetalab_thetalab_gsgs_cm60->Write();
   h_thetalab_thetalab_gsgs_cm70->Write();
   h_thetalab_thetalab_gsgs_cm80->Write();
   h_thetalab_thetalab_gsgs_cm90->Write();

   //  phi vs phi
   h_philab_philab->Write();
   h_philab_philab_cutphi->Write();
   h_philab_philab_cutphi_12c->Write();

   // vertex of tracks
   h_verz->Write();
   h_verz_cut12c_ela->Write();
   h_verz01->Write();
   h_verz_gsgs->Write();
   h_verz_gsgs_cm50->Write();
   h_verz_gsgs_cm60->Write();
   h_verz_gsgs_cm70->Write();
   h_verz_gsgs_cm80->Write();
   h_verz_gsgs_cm90->Write();

   h_verz_gsex->Write();
   h_verz_exex->Write();
   h_verxy->Write();
   h_verxy_cut12c_ela->Write();
   h_verxy_gsgs->Write();
   h_verxy_gsex->Write();
   h_verxy_exex->Write();

   h_verxz->Write();
   h_veryz->Write();
   h_ntra_verz->Write();
   h_ntraver_ntra->Write();
   h_ntraver_ntra_cutphi->Write();
   h_ntraver_ntra_cutphi_2tra->Write();
   h_ntraver_ntra_cut12c_ela->Write();

   // Beam energy estimation
   h_Ebeam->Write();
   h_Ebeam_gsgs->Write();
   h_Ebeam_gsgs_cm50->Write();
   h_Ebeam_gsgs_cm60->Write();
   h_Ebeam_gsgs_cm70->Write();
   h_Ebeam_gsgs_cm80->Write();
   h_Ebeam_gsgs_cm90->Write();
   h_Ebcm->Write();
   h_Ebcm_gsgs->Write();
   h_Ebcm_gsgs_cm50->Write();
   h_Ebcm_gsgs_cm60->Write();
   h_Ebcm_gsgs_cm70->Write();
   h_Ebcm_gsgs_cm80->Write();
   h_Ebcm_gsgs_cm90->Write();
   h_sumkine_Ebeam_gsgs->Write();
   h_sumkine_ver_Ebeam_gsgs->Write();
   h_sumkine_Ebeam_gsgs_cm90->Write();

   // Others ...
   h_lastpoint_comp->Write();
   h_vertex_comp->Write();
   h_distance_comp->Write();
   h_range_comp->Write();
   h_thetalab_comp->Write();
   h_philab_comp->Write();
   h_range_comp_gsgs->Write();
   h_kineE_comp_gsgs->Write();
   h_thetalab_comp_gsgs->Write();
   h_philab_comp_gsgs->Write();

#if defined (vertex_index) || defined (vertex_depth)
   for(Int_t i = 0; i < n_group; i++){
      h_verz_gsgs_index[i]->Write();
      h_verz_gsex_index[i]->Write();
      h_verz_exex_index[i]->Write();
      h_verxy_gsgs_index[i]->Write();
      h_verxy_gsex_index[i]->Write();
      h_verxy_exex_index[i]->Write();
      h_kineE_thetalab_gsgs_index[i]->Write();
      h_kineE_thetalab_gsex_index[i]->Write();
      h_kineE_thetalab_exex_index[i]->Write();
      h_thetalab_thetalab_gsgs_index[i]->Write();
      h_thetalab_thetalab_gsex_index[i]->Write();
      h_thetalab_thetalab_exex_index[i]->Write();
   }

#endif

   // save cuts
   cut12c->Write("charge_range_12c");
   theta_gsgs->Write("theta_theta_gsgs");
   theta_gsex->Write("theta_theta_gsex");
   theta_exex->Write("theta_theta_exex");

   // save lines
   kine_gsgs_0  ->Write("kine_12c_gsgs_E69_0_z0"  );
   kine_gsgs_50 ->Write("kine_12c_gsgs_E66_4_z50" );
   kine_gsgs_100->Write("kine_12c_gsgs_E63_7_z100");
   kine_gsgs_150->Write("kine_12c_gsgs_E60_9_z150");
   kine_gsgs_200->Write("kine_12c_gsgs_E58_1_z200");
   kine_gsgs_250->Write("kine_12c_gsgs_E55_1_z250");
   kine_gsgs_300->Write("kine_12c_gsgs_E52_1_z300");
   kine_gsgs_350->Write("kine_12c_gsgs_E48_8_z350");
   kine_gsgs_400->Write("kine_12c_gsgs_E45_4_z400");
   kine_gsgs_450->Write("kine_12c_gsgs_E41_8_z450");
   kine_gsgs_500->Write("kine_12c_gsgs_E38_1_z500");
   kine_gsgs_550->Write("kine_12c_gsgs_E34_1_z550");
   kine_gsgs_600->Write("kine_12c_gsgs_E29_7_z600");
   kine_gsgs_650->Write("kine_12c_gsgs_E25_0_z650");
   kine_gsgs_700->Write("kine_12c_gsgs_E19_7_z700");
   kine_gsgs_750->Write("kine_12c_gsgs_E13_6_z750");
   kine_gsgs_800->Write("kine_12c_gsgs_E06_2_z800");

   kine_gsex_0->Write( "kine_12c_gsex_E69_0_z0" );

   kine_exex_0->Write( "kine_12c_exex_E69_0_z0" );

   angle_gsgs_0->Write("angle_12c_gsgs_E69_0_z0");
   angle_gsex_0->Write("angle_12c_gsex_E69_0_z0");
   angle_exex_0->Write("angle_12c_exex_E69_0_z0");

   sum_kine_beam->Write("sum_kine_beam");
   Results->Close();

   // cout of information
   std::cout << "                                                                " << std::endl << std::endl;
   std::cout << "Finished analysis. Displaying analysis results..." << std::endl;
   std::cout << " Number of all events: " << nevents << std::endl;
   std::cout << " Maximum radius of hits: " << max_r_max << " mm, Trigger radius: " << r_tri << " mm" << std::endl;
   std::cout << " Maximum radius of hits in 12c12c events: " << max_r_max_12c << " mm, Trigger radius: " << r_tri_12c << " mm" << std::endl;
   std::cout << std::endl;
#ifdef peak_check
   std::cout << " peak1 events: " << peak1.size() << std::endl;
   std::cout << "   Event of peak1 " << std::endl;
   for (auto &eventIndex: peak1){
      std::cout << eventIndex << ", " << std::flush;
   }
   std::cout << std::endl << std::endl;
   std::cout << " peak2 events: " << peak2.size() << std::endl;
   std::cout << "   Event of peak2: " << std::endl;
   for (auto &eventIndex: peak2){
      std::cout << eventIndex << ", " << std::flush;
   }
   std::cout << std::endl << std::endl;
   std::cout << " peak3 events: " << peak3.size() << std::endl;
   std::cout << "   Event of peak3: " << std::endl;
   for (auto &eventIndex: peak3){
      std::cout << eventIndex << ", " << std::flush;
   }
   std::cout << std::endl << std::endl;

#endif

#ifdef vertex_check
   std::cout << " 12c12c elastic events: " << n_12c12c << std::flush;
   std::cout << " (no vertex: " << n_no_vertex << " events.)" << std::endl;
   if(n_no_vertex != track2_ver0.at(1).size()){
      std::cout << "  something is wrong with the vector. vector size: " << track2_ver0.at(1).size() << std::endl;
   }
   if (n_no_vertex != 0){
      for (Int_t i = 0; i < track2_ver0.at(1).size(); i++){
         if (i == 0){
            std::cout << "  runNum: " << track2_ver0[0].at(i)<< std::flush;
            std::cout << "; events: " << track2_ver0.at(1).at(i) << std::flush;
         }
         else if(track2_ver0.at(0).at(i) != track2_ver0.at(0).at(i - 1)){
            std::cout << std::endl;
            std::cout << "  runNum: " << track2_ver0[0].at(i)<< std::flush;
            std::cout << "; events: " << track2_ver0.at(1).at(i) << std::flush;
         }
         else{
            std::cout << ", " << track2_ver0.at(1).at(i) << std::flush;
         }
      }
      std::cout << std::endl << std::endl;
   }
   else {
      std::cout << std::endl;
   }

#endif

#ifdef check_verz0
   std::cout << " vertex 0 < z <= 1 events: " << n_ver_z_01 << " events (12c12c)." << std::endl;
   if(n_ver_z_01 != vertex_z_01.at(1).size()){
      std::cout << "  something is wrong with the vector. vector size: " << vertex_z_01.at(1).size() << std::endl;
   }
   for (Int_t i = 0; i < vertex_z_01.at(1).size(); i++){
      if (i == 0){
         std::cout << "  runNum: " << vertex_z_01[0].at(i)<< std::flush;
         std::cout << "; events: " << vertex_z_01.at(1).at(i) << std::flush;
      }
      else if(vertex_z_01.at(0).at(i) != vertex_z_01.at(0).at(i - 1)){
         std::cout << std::endl;
         std::cout << "  runNum: " << vertex_z_01[0].at(i)<< std::flush;
         std::cout << "; events: " << vertex_z_01.at(1).at(i) << std::flush;
      }
      else{
         std::cout << ", " << vertex_z_01.at(1).at(i) << std::flush;
      }
   }
   std::cout << std::endl << std::endl;

#endif

   // stop timer
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl;
   cout << " Real time: " << rtime << " s, CPU time: " << ctime << " s" << endl << endl;

}

TGraph* ReadKinematics(TString kineFileName, TString ReturnType){

   TString kineder = "./two-body_kine_files/out_kinema_69.5/";
   TString kineFile = kineder + "kine_12c12c_" + kineFileName + ".txt";
   //   std::vector<Double_t> ThetacmS, ThetalabR, ElabR, ThetalabS, ElabS, MomlabR;
   std::vector<Double_t> Thetacm, Theta3lab, E3lab, Theta4lab, E4lab;
#ifdef debug_mode
   // Debug: Check if file exists and current working directory
   std::cout << "DEBUG ReadKinematics: pwd = " << gSystem->pwd() << std::endl;
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
   /*
   while (!kineStr->eof()){
   //	      *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >> EnerLabSca[numKin];
	   *kineStr >> ThetaLabRec[numKin] >> EnerLabRec[numKin];
      numKin++;
      }
   */
   double tcm, t3lab, e3lab, t4lab, e4lab;
   string line;
   while (getline(*kineStr, line)){
      if(line.empty() || line[0] == '#'){
         continue;
      }
      istringstream iss(line);
      if (iss >> tcm >> t3lab >> e3lab >> t4lab >> e4lab){
         Thetacm.push_back(tcm);
         Theta3lab.push_back(t3lab);
         E3lab.push_back(e3lab);
         Theta4lab.push_back(t4lab);
         E4lab.push_back(e4lab);
      }
   }

   /*
   while ( *kineStr >> tlabr >> elabr){
      ThetalabR.push_back(tlabr);
      ElabR.push_back(elabr);
   }
   */
#ifdef debug_mode
      std::cout << "DEBUG ReadKinematics: Successfully read " << Theta3lab.size() << " points" << std::endl;
#endif
   if (Theta3lab.size() == 0){
      std::cout << " Warning : No data read from kinematics file: " << kineFile << std::endl;
   }
   //   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   TGraph *kine  = new TGraph(Theta3lab.size(), Theta3lab.data(), E3lab.data());
   TGraph *angle = new TGraph(Theta3lab.size(), Theta3lab.data(), Theta4lab.data());
   //   kine->Sort();

   if (ReturnType == "kine"){
      return kine;
   }
   else if (ReturnType == "angle"){
      return angle;
   }
   else {
      std::cout << "Warning: Invalid ReturnType specified. Returning empty TGraph." << std::endl;
      return new TGraph();
   }

}

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{

   // in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables

   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) -
                 (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) /
                   (2 * pow(m2, 2)) + s + u - pow(m2, 2));
   Ex = m4_ex - m4;

   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;

   // for inverse kinematics Note: this angle corresponds to the recoil
    theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   /*theta_cm = acos((pow(s, 2) + s * (2 * u - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m4_ex, 2) - pow(m3, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m4_ex, 2), pow(m3, 2))));*/

   theta_cm = theta_cm * TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

void draw_ind(TString cname, TString states, Int_t n_group, Int_t n_h_z,
   TH1D* h_verz, TH2F* h_verxy, TH2F* h_E_theta, TH2F* h_theta_theta,
   std::vector<TH1D*> &h_verz_i, std::vector<TH2F*> &h_verxy_i, std::vector<TH2F*> &h_E_theta_i, std::vector<TH2F*> &h_theta_theta_i,
   TGraph *kine_gsgs, TGraph *kine_gsex, TGraph *kine_exex, TGraph *ang_gsgs, TGraph *ang_gsex, TGraph *ang_exex)
{
   TCanvas *c = new TCanvas(cname, cname, 1600, 1000);
   c->Divide(n_h_z,4);
   c->cd(1);
   h_verz->SetDirectory(0);
   h_verz->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz->SetTitle(Form("vertex z (12c12c, %s)", states.Data()));
   gPad->SetLogy();
   h_verz->Draw();
   c->cd(n_h_z + 1);
   h_verxy->SetDirectory(0);
   h_verxy->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy->SetTitle(Form("vertex xy (12c12c, %s)", states.Data()));
   gPad->SetLogz();
   h_verxy->Draw("colz");
   c->cd(n_h_z * 2 +1);
   h_theta_theta->SetDirectory(0);
   h_theta_theta->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_theta_theta->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_theta_theta->SetTitle(Form("Theta_LAB Theta_LAB (12c12c)"));
   gPad->SetLogz();
   h_theta_theta->SetMinimum(1);
   h_theta_theta->Draw("colz");
   ang_gsgs->Draw("same");
   ang_gsex->Draw("same");
   ang_exex->Draw("same");
   c->cd(n_h_z * 3 +1);
   h_E_theta->SetDirectory(0);
   h_E_theta->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_E_theta->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_E_theta->SetTitle(Form("KinE Theta_LAB (12c12c, %s)", states.Data()));
   gPad->SetLogz();
   h_E_theta->SetMinimum(1);
   h_E_theta->Draw("colz");
   kine_gsgs->Draw("same");
   kine_gsex->Draw("same");
   kine_exex->Draw("same");
   for (Int_t i = 0; i < n_group; i++){
      c->cd(i + 2);
      h_verz_i[i]->SetDirectory(0);
      h_verz_i[i]->GetXaxis()->SetTitle("vertex z [mm]");
      //      h_verz_i[i]->SetTitle(Form("h_verz_%s_z_%d--%d;Vertex Z [mm]", states.Data(),i*100, (i+1)*100));
      gPad->SetLogy();
      h_verz_i[i]->Draw();
      c->cd(n_h_z + i + 2);
      h_verxy_i[i]->SetDirectory(0);
      h_verxy_i[i]->GetXaxis()->SetTitle("vertex x [mm]");
      h_verxy_i[i]->GetYaxis()->SetTitle("vertex y [mm]");
      //      h_verxy_i[i]->SetTitle(Form("h_verxy_%s_z_%d--%d;Vertex X [mm];Vertex Y [mm]", states.Data(), i*100, (i+1)*100));
      gPad->SetLogz();
      h_verxy_i[i]->Draw("colz");
      c->cd(2 * n_h_z + i + 2);
      h_theta_theta_i[i]->SetDirectory(0);
      h_theta_theta_i[i]->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
      h_theta_theta_i[i]->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
      //      h_theta_theta_i[i]->SetTitle(Form("h_theta_theta_%s_z_%d--%d;Theta_lab;Theta_lab", states.Data(), i*100, (i+1)*100));
      gPad->SetLogz();
      h_theta_theta_i[i]->SetMinimum(1);
      h_theta_theta_i[i]->Draw("colz");
      ang_gsgs->Draw("same");
      ang_gsex->Draw("same");
      ang_exex->Draw("same");
      c->cd(3 * n_h_z + i + 2);
      h_E_theta_i[i]->SetDirectory(0);
      h_E_theta_i[i]->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
      h_E_theta_i[i]->GetYaxis()->SetTitle("roughKineE [MeV]");
      //      h_E_theta_i[i]->SetTitle(Form("h_KineE_theta_%s_z_%d--%d;KinE;Theta_lab", states.Data(), i*100, (i+1)*100));
      gPad->SetLogz();
      h_E_theta_i[i]->SetMinimum(1);
      h_E_theta_i[i]->Draw("colz");
      kine_gsgs->Draw("same");
      kine_gsex->Draw("same");
      kine_exex->Draw("same");

   }
   //   c->SaveAs(Form("can_output/check_vd_76matm_%s_vertex_%s.pdf", cname.Data(), states.Data()));

}

void draw_dep(TString cname, TString states, TString LineType, Int_t n_group, Int_t n_div,
            std::vector<TH2F*> &h_E_theta_i, std::vector<TGraph*> &kine_i)
{

   Int_t n_l = 0;
   Int_t n_kine = kine_i.size();
   if (LineType == "full" || LineType == "center"){
      std::cout << std::endl;
      std::cout << "Canvas info; canvas: " << cname << ", states: " << states << ", type: " << LineType << std::endl;
   }
   else{
      std::cout << "Warning!! Invalid LineType specified. No lines will be drawn." << std::endl;
   } 

   TCanvas *c = new TCanvas(cname, cname, 1600, 1000);
   c->Divide(n_div, n_div);
   for (Int_t i = 0; i < n_group; i++){
      c->cd(i + 1);
      h_E_theta_i[i]->SetDirectory(0);
      //      h_E_theta_i[i]->SetStats(0);
      gPad->SetLogz(1);
      h_E_theta_i[i]->SetMinimum(1);
      h_E_theta_i[i]->Draw("colz");

      if (LineType == "full"){
         if (n_kine > 2 * i + 1){
         kine_i[2 * i]->Draw("PL same");
         }
         if (n_kine > 2 * i + 2){
            kine_i[2 * i + 1]->Draw("same");
         }
         if (n_kine > 2 * i + 3){
            kine_i[2 * i + 2]->Draw("PL same");
         }
      }
      else if (LineType == "center"){
         if (n_kine > 2 * i + 2){
            kine_i[2 * i + 1]->Draw("same");
         }
      }
      n_l ++;
   }
   c->Update();
   //   c->SaveAs(Form("can_output/c12_" + states + "_kine_" + LineType + "_E69.5_z0-%d.pdf", (n_kine - 1) * 50));

}

std::vector<Double_t> cal_Ebeam_para(){

   // set data
   /*
   //  60.7 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 60.724}, {50, 57.862}, {100, 54.873}, {150, 51.775}, {200, 48.529}, {250, 45.143}, {300, 41.562},
         {350, 37.779},{400, 33.736},{450, 29.363}, {500, 24.589}, {550, 19.089},{600, 13.089}, {650, 5.513},
         {680, 0.360}
   };

   //  69.0 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 69.027}, {50, 66.411}, {100, 63.727}, {150, 60.949}, {200, 58.083}, {250, 55.105}, {300, 52.061},
         {350, 48.784},{400, 45.411},{450, 41.844}, {500, 38.084}, {550, 34.062},{600, 29.719}, {650, 24.984},
         {700, 19.708}, {750, 13.621}, {800, 6.194}, {830, 0.832}, {835, 0.271}
   };
   */

   //  69.5 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 69.500}, {50, 66.899}, {100, 64.228}, {150, 61.468}, {200, 58.618}, {250, 55.662}, {300, 52.594},
         {350, 49.391},{400, 46.046},{450, 42.514}, {500, 38.797}, {550, 34.822},{600, 30.545}, {650, 25.892},
         {700, 20.728}, {750, 14.821}, {800, 7.705}, {840, 0.723}, {848, 0.028}
   };

   /*
   //  70.0 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 70.000}, {50, 67.414}, {100, 64.757}, {150, 62.016}, {200, 59.183}, {250, 56.251}, {300, 53.203},
         {350, 50.031},{400, 46.713},{450, 43.222}, {500, 39.545}, {550, 35.617},{600, 31.411}, {650, 26.837},
         {700, 21.788}, {750, 16.056}, {800, 9.230}, {850, 0.695}, {858, 0.019}
   };
   */

   Int_t n_data = lise_data.size();
   std::vector<std::vector<Double_t>> Ebeam(2, std::vector<Double_t>(n_data, 0));
   for (Int_t i = 0; i < n_data; i++){
      Ebeam.at(0).at(i) = lise_data[i].first;
      Ebeam.at(1).at(i) = lise_data[i].second;
   }
   TGraph *h_Ebeam_est = new TGraph(n_data, Ebeam.at(0).data(), Ebeam.at(1).data());
   //   TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0, 1000);
   TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[5]*x^5", 0, 1000);

   f1->FixParameter(0, Ebeam.at(1).at(0));
   h_Ebeam_est->Fit(f1, "QRN", "", 0, Ebeam.at(0).at(n_data - 1));

   TCanvas *c0 = new TCanvas("c0", "c0");
   h_Ebeam_est->SetMarkerStyle(20);
   h_Ebeam_est->SetMarkerSize(1.2);
   h_Ebeam_est->SetMarkerColor(kBlue);
   h_Ebeam_est->SetLineColor(kRed);
   h_Ebeam_est->GetXaxis()->SetTitle("Depth [mm]");
   h_Ebeam_est->GetYaxis()->SetTitle("beam energy [MeV]");
   h_Ebeam_est->SetTitle("Beam energy vs Depth");
   h_Ebeam_est->GetXaxis()->SetLimits(0, 1000);
   h_Ebeam_est->Draw("AP");
   f1->SetNpx(1000);
   f1->SetLineColor(kBlack);
   f1->Draw("P same");

   Double_t a = f1->GetParameter(0);
   Double_t b = f1->GetParameter(1);
   Double_t c = f1->GetParameter(2);
   Double_t d = f1->GetParameter(3);
   Double_t e = f1->GetParameter(4);
   Double_t f = f1->GetParameter(5);
   std::vector<Double_t> Ebeam_para;
   Ebeam_para.push_back(a);
   Ebeam_para.push_back(b);
   Ebeam_para.push_back(c);
   Ebeam_para.push_back(d);
   Ebeam_para.push_back(e);
   Ebeam_para.push_back(f);

   std::cout << std::setprecision(4) << std::endl;
   std::cout << "Set estimation of beam energy." << std::endl;
   std::cout << "  function: a + b * x + c * x^2 + d * x^3 + e * x^4 + f * x^5" << std::endl;
   std::cout << "  a:" << a << ", b:" << b << ", c:" << c << ", d:" << d << ", e:" << e << ", f:" << f << std::endl;

   return Ebeam_para;
}

Double_t est_Ebeam(std::vector<Double_t> &Ebeam_para, Double_t vertz){

   Double_t Ebeam = -1;
   if (Ebeam_para.size() != 6){
      std::cout << "Error: Ebeam_para should have 6 parameters. Current size: " << Ebeam_para.size() << std::endl;
      return -1;
   }
   else {
      Double_t a = Ebeam_para[0];
      Double_t b = Ebeam_para[1];
      Double_t c = Ebeam_para[2];
      Double_t d = Ebeam_para[3];
      Double_t e = Ebeam_para[4];
      Double_t f = Ebeam_para[5];
      Ebeam = a + b * vertz + c * pow(vertz, 2) + d * pow(vertz, 3) + e * pow(vertz, 4) + f * pow(vertz, 5);
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
