//#define nom_check
//#define peak_check
#define c12_check
#define states_vertex
#include <fstream>
#include "TFile.h"
#include "TObject.h"
#include "TCanvas.h"

TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);
void draw_ind(TString cname, TString states, Int_t n_group, Int_t n_h_z, TH1D* h_ver, TH2F* h_verxy, TH2F* h_E_theta, TH2F* h_theta_theta, 
         std::vector<TH1D*> &h_verz_i, std::vector<TH2F*> &h_verxy_i, std::vector<TH2F*> &h_E_theta_i, std::vector<TH2F*> &h_theta_theta_i,
         TGraph *kine_gsgs, TGraph *kine_gsex, TGraph *kine_exex, TF1 *ang_gsgs, TGraph *ang_gsex, TGraph *ang_exex);

void check_vd_76matm(){
   //copy from kine.C 2026/05/06 12:20

   //set timer
   TStopwatch timer;
   timer.Start();

   // set parameters
   Double_t del_phi = 10; // cut value; phi1 - phi2 - 180 deg < del_phi
   /*
   // Masses.
   double u_to_MeV = 931.49401;
   double m_17C = 17.022578650 * u_to_MeV;
   double m_18C = 18.026751930 * u_to_MeV;
   double m_p = 1.007825031898 * u_to_MeV;
   double m_d = 2.014101777844 * u_to_MeV;

   // Beam energy.
   double E_beam = 27.091 * 17.022578650;
   */

   // directory
   //   double vd_val = 3.06;
   //   double vd_val = 3.50;
   //   double vd_val = 3.89;
   //   double vd_val = 4.00;
   double vd_val = 4.05;
   //   double vd_val = 4.07;
   //   double vd_val = 4.08;
   //   double vd_val = 4.10;
   //   double vd_val = 4.20;
   //   double vd_val = 4.50;

   // files.
   //   std::vector runNums = {52};
   //   std::vector runNums = {50,51,52,53,54,55,56,57,58};

   std::vector runNums = {
      28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,50,
      51,52,53,54,55,56,57,58,62,63,64,66,67,68,69,70,71,75,76,77,
      78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,95,96,97,98,99,
      100,101,102,103,104,105,106,107,108,109,110,111,112
   };

#ifdef states_vertex
   const Int_t n_group = 7;
   const Int_t verz_h = 600;
#endif

   Int_t run_start = runNums.front();
   Int_t run_end = runNums.back();
   TFile * Results = new TFile(Form("data/check_vd_results_run%d-run%d.root", run_start, run_end),"recreate");
      
   FairRunAna *run = new FairRunAna(); // Forcing a dummy run
   //   TString outfname="./canvas_kine.root";
   //   TFile *outfile=new TFile(outfname,"recreate");

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
   TGraph *kine_12c12c_gsgs_60_7 = ReadKinematics("./two-body_kine_files/kine_12c12c_gsgs_60.7.txt");
   TGraph *kine_12c12c_gsex_60_7 = ReadKinematics("./two-body_kine_files/kine_12c12c_gsex_4.44_60.7.txt");
   TGraph *kine_12c12c_exex_60_7 = ReadKinematics("./two-body_kine_files/kine_12c12c_exex_4.44_60.7.txt");
   TGraph *angle_12c12c_gsex_60_7 = ReadKinematics("./two-body_kine_files/angle_12c12c_gsex_4.44_60.7.txt");
   TGraph *angle_12c12c_exex_60_7 = ReadKinematics("./two-body_kine_files/angle_12c12c_exex_4.44_60.7.txt");
   /*
   TGraph *kine_dp_gs = ReadKinematics("./kineFiles/kine17C_dp_gs.txt");
   TGraph *kine_dd_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dd_gs_25MeVu.txt");
   TGraph *kine_dp_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dp_gs_25MeVu.txt");
   */

   // Characteristic definitions
   bool tracks_vertex = false;
   bool alpha_tracks = false;
   bool proton_tracks = false;
   const Int_t n_h_z = n_group +1;
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
   Int_t n_no_vertex = 0;
   Double_t rad = -100;
   Double_t r_tem = 0;
   Double_t r_max = 0;
   Double_t r_tri = 300;
   Double_t max_r_max = 0;
   Double_t track_lastx = 0;
   Double_t track_lasty = 0;
   Double_t track_lastz = 0;
   Double_t vtx = 0;
   Double_t vty = 0;
   Double_t vtz = 0;
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
   std::vector<Int_t> track6(0);
   std::vector<Int_t> peak1(0);
   std::vector<Int_t> peak2(0);
   std::vector<Int_t> peak3(0);

   // Histogram definitions.
   // ... TH1 hist for checking something.
   TH1D *h_ntra = new TH1D("h_ntra", "h_ntra;NTracks", 11, -0.5, 10.5);
   TH1D *h_rmax = new TH1D("h_rmax", "h_rmax;Rmax [mm]", 150, 0, 300);
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
   TH2F *h_dEdx_range_cutphi_alpha = new TH2F("h_dEdx_range_cutphi_alpha", "h_dEdx_range_cutphi_alpha;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_cutphi_proton = new TH2F("h_dEdx_range_cutphi_proton", "h_dEdx_range_cutphi_proton;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_cutmul1 = new TH2F("h_dEdx_range_cutmul1", "h_dEdx_range_cutmul1;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_cutmul2 = new TH2F("h_dEdx_range_cutmul2", "h_dEdx_range_cutmul2;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);

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
   TH2F *h_kineE_thetalab_alpha = new TH2F("h_kineE_thetalab_alpha", "h_kineE_thetalab_alpha;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 500, 0, 10);
   TH2F *h_kineE_thetalab_proton = new TH2F("h_kineE_thetalab_proton", "h_kineE_thetalab_proton;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 250, 0, 5);

   // ... angle correlations
   TH2F *h_range_thetalab = new TH2F("h_range_thetalab", "h_range_thetalab", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutphi = new TH2F("h_range_thetalab_cutphi", "h_range_thetalab_cutphi", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutphi_2tra = new TH2F("h_range_thetalab_cutphi_2tra", "h_range_thetalab_cutphi_2tra", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cut12c_ela = new TH2F("h_range_thetalab_cut12c_ela", "h_range_thetalab_cut12c_ela", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutalpha = new TH2F("h_range_thetalab_cutalpha", "h_range_thetalab_cutalpha", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutproton = new TH2F("h_range_thetalab_cutproton", "h_range_thetalab_cutproton", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_gsgs = new TH2F("h_range_thetalab_gsgs", "h_range_thetalab_gsgs", 50, 0, 100, 300, 0, 600);
   TH2F *h_range_thetalab_gsex = new TH2F("h_range_thetalab_gsex", "h_range_thetalab_gsex", 50, 0, 100, 300, 0, 600);
   TH2F *h_range_thetalab_exex = new TH2F("h_range_thetalab_exex", "h_range_thetalab_exex", 50, 0, 100, 300, 0, 600);
   // ... .. theta vs theta
   TH2F *h_thetalab_thetalab_cutphi = new TH2F("h_thetalab_thetalab_cutphi", "h_thetalab_thetalab_cutphi", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_cutphi_2tra = new TH2F("h_thetalab_thetalab_cutphi_2tra", "h_thetalab_thetalab_cutphi_2tra", 200, 0, 100, 200, 0, 100);
   TH2F *h_thetalab_thetalab_cut12c_ela = new TH2F("h_thetalab_thetalab_cut12c_ela", "h_thetalab_thetalab_cut12c_ela", 200, 0, 100, 200, 0, 100);

   // ... .. phi vs phi
   TH2F *h_philab_philab = new TH2F("h_philab_philab", "h_philab_philab", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi = new TH2F("h_philab_philab_cutphi", "h_philab_philab_cutphi", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi_12c = new TH2F("h_philab_philab_cutphi_12c", "h_philab_philab_cutphi_12c", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi_alpha = new TH2F("h_philab_philab_cutphi_alpha", "h_philab_philab_cutphi_alpha", 360, -180, 180, 360, -180, 180);
   TH2F *h_philab_philab_cutphi_proton = new TH2F("h_philab_philab_cutphi_proton", "h_philab_philab_cutphi_proton", 360, -180, 180, 360, -180, 180);

   // ... vertex of tracks
   TH1D *h_verz = new TH1D("h_verz", "h_verz;Vertex Z [mm]", 112, -10, 1010);
   TH1D *h_verz_cut12c_ela = new TH1D("h_verz_cut12c_ela", "h_verz_cut12c_ela;Vertex Z [mm]", 112, -10, 1010);
   TH1D *h_verz_gsgs = new TH1D("h_verz_gsgs", "h_verz_gsgs;Vertex Z [mm]", 112, -10, 1010);
   TH1D *h_verz_gsex = new TH1D("h_verz_gsex", "h_verz_gsex;Vertex Z [mm]", 112, -10, 1010);
   TH1D *h_verz_exex = new TH1D("h_verz_exex", "h_verz_exex;Vertex Z [mm]", 112, -10, 1010);

   TH2F *h_ntraver_ntra = new TH2F("h_ntraver_ntra", "h_ntraver_ntra", 11, -0.5, 10.5, 11, -0.5, 10.5);
   TH2F *h_ntra_verz = new TH2F("h_ntra_verz", "h_ntra_verz", 112, -10, 1010, 11, -0.5, 10.5);
   TH2F *h_verxy = new TH2F("h_verxy", "h_verxy", 50, -50, 50,  50, -50, 50);
   TH2F *h_verxz = new TH2F("h_verxz", "h_verxz", 112, -10, 1010, 50, -50, 50);
   TH2F *h_veryz = new TH2F("h_veryz", "h_veryz", 112, -10, 1010, 50, -50, 50);
   TH2F *h_verxy_cut12c_ela = new TH2F("h_verxy_cut12c_ela", "h_verxy_cut12c_ela", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_gsgs = new TH2F("h_verxy_gsgs", "h_verxy_gsgs", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_gsex = new TH2F("h_verxy_gsex", "h_verxy_gsex", 50, -50, 50, 50, -50, 50);
   TH2F *h_verxy_exex = new TH2F("h_verxy_exex", "h_verxy_exex", 50, -50, 50, 50, -50, 50);

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
   TH2F *h_nalp_ntra = new TH2F("h_nalp_ntra", "h_nalp_ntra", 11, -0.5, 10.5, 11, -0.5, 10.5);
   TH2F *h_npro_ntra = new TH2F("h_npro_ntra", "h_npro_ntra", 11, -0.5, 10.5, 11, -0.5, 10.5);

   TH2F *h_Esmall_range = new TH2F("h_ESmallVTotalRange", "h_ESmallVTotalRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *h_Ebig_rangebig = new TH2F("h_EBigVBigRange", "h_EBigVBigRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *h_kineE_thetalab_2H  = new TH2F("h_kineE_thetalab_2H", "h_kineE_thetalab_2H", 180, 0, 180, 250, 0, 20);
   TH2F *h_kineE_thetalab_1H  = new TH2F("h_kineE_thetalab_1H", "h_kineE_thetalab_1H", 180, 0, 180, 250, 0, 20);

   // ... index hist.
#ifdef states_vertex
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
            = new TH2F("h_kineE_thetalab_gsgs_index_other","h_kineE_gsgs_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 600, 0, 60);
         h_kineE_thetalab_gsex_index[i] 
            = new TH2F("h_kineE_thetalab_gsex_index_other","h_kineE_gsex_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 600, 0, 60);
         h_kineE_thetalab_exex_index[i] 
            = new TH2F("h_kineE_thetalab_exex_index_other","h_kineE_exex_thetalab_z_others;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 600, 0, 60);
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
               Form("h_kineE_thetalab_gsgs_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 180, 0, 180, 600, 0, 60);
         h_kineE_thetalab_gsex_index[i] 
            = new TH2F(Form("h_kineE_thetalab_gsex_index_%d", i),
               Form("h_kineE_thetalab_gsex_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 180, 0, 180, 600, 0, 60);
         h_kineE_thetalab_exex_index[i] 
            = new TH2F(Form("h_kineE_thetalab_exex_index_%d", i),
               Form("h_kineE_thetalab_exex_z_%d--%d;#theta_{LAB} [deg];roughKinE [MeV]", i*100, (i+1)*100), 180, 0, 180, 600, 0, 60);
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

   // Line definitions.
   TF1 *xy90 = new TF1("xy90", "-x + 90", 0, 90);

   std::cout << std::endl;
   for (int runNum: runNums) {
      // Open the digitalization file and get the TTree.
      //      TString unpackFileName = TString::Format("./vd_check_data/run_%04d.root", runNum);
      TString unpackFileName = TString::Format("./vd_check_data/vd%.2f_files/run_%04d.root", vd_val, runNum);
      TFile *unpackFile = new TFile(unpackFileName, "READ");
      TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
      int nUnpackEvents = unpackTree->GetEntries();
      //      std::cout << "Number of unpacked events in run " << runNum << ": " << nUnpackEvents << std::endl;
      std::cout << "Reading run:" << runNum << ", events:" << nUnpackEvents << std::endl;
      int nEventsWith2Tracks = 0;
      // Creare the TTreeReader to read the AtTrackingEvents and simulation.
      TTreeReader unpackReader("cbmsim", unpackFile);
      //      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");
      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEventModified");

      // Loop over events.
      for (int i = 0; i < nUnpackEvents; i++) {
         unpackReader.Next();
         // First, we obtain some rough kinematics just by using the AtPatternEvent.
         AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
         if (!patternEvent){continue;}

         // We want to focus on events with 2 or less tracks for now.
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

         if (ntrack == 6){
            track6.push_back(i);
         }
         //         int maxTrackNum{4};
         //         if (tracks.size() > maxTrackNum) continue;
         //         if (tracks.size() == 2) nEventsWith2Tracks++;
         int trackIndex = 0;
         double thetaLABArray[2] = {0.0,0.0};
         double estimatedKineEArray[2] = {0.0,0.0};

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
            if (cut12c->IsInside(track_range[itrack], track_charge[itrack])){
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

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_p, m_18C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), track_KinE[itrack]   );
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
         }
         if(abs(abs(track_phi[0] - track_phi[1]) - 180) < del_phi){
            // pid
            for (Int_t k = 0; k < ntrack; k++){
               if (cut12c->IsInside(track_range[k], track_charge[k])) {
                  //                  std::cout << "event num:" << i << ", tracks: " << ntrack << ", track: "<< k << ", particle: 12c " << std::endl;
                  track_12c[k]=true;
                  h_philab_philab_cutphi_12c -> Fill(track_phi[0], track_phi[1]);                  
               }
               else if (cutalpha->IsInside(track_range[k], track_charge[k])) {
                  //                  std::cout << "event num:" << i << ", tracks: " << ntrack << ", track: " << k << ", particle: alpha " << std::endl;
                  track_alpha[k]=true;
                  h_range_thetalab_cutalpha->Fill(track_theta[k], track_range[k]);
                  h_kineE_thetalab_alpha->Fill(track_theta[k], track_KinE[k]);
                  h_philab_philab_cutphi_alpha -> Fill(track_phi[0], track_phi[1]);
               }
               else if (cutproton->IsInside(track_range[k], track_charge[k])) {
                  //                  std::cout << "event num:" << i << ", tracks: " << ntrack << ", track: " << k << ", particle: proton " << std::endl;
                  track_proton[k]=true;
                  h_range_thetalab_cutproton->Fill(track_theta[k], track_range[k]);
                  h_kineE_thetalab_proton->Fill(track_theta[k], track_KinE[k]);
                  h_philab_philab_cutphi_proton -> Fill(track_phi[0], track_phi[1]);
               }
               if (cutmul1->IsInside(track_range[k], track_charge[k])) {
                  h_dEdx_range_cutmul1 -> Fill(track_range[k], track_dedx[k]);                  
               }
               else if (cutmul2->IsInside(track_range[k], track_charge[k])) {
                  h_dEdx_range_cutmul2 -> Fill(track_range[k], track_dedx[k]);                  
               }
               // pid using dEdx-range
                  //                  h_dEdx_range_cutphi_12c->Fill(track_range[k], track_dedx[k]);
                  //                  h_dEdx_range_cutphi_alpha->Fill(track_range[k], track_dedx[k]);
                  //                  h_dEdx_range_cutphi_proton->Fill(track_range[k], track_dedx[k]);
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
            if (ntrack == 2){
               h_charge_range_cutphi_2tra -> Fill(track_range[0], track_charge[0]);
               h_charge_range_cutphi_2tra -> Fill(track_range[1], track_charge[1]);
               h_range_thetalab_cutphi_2tra -> Fill(track_theta[0], track_range[0]);
               h_range_thetalab_cutphi_2tra -> Fill(track_theta[1], track_range[1]);
               h_thetalab_thetalab_cutphi_2tra -> Fill(track_theta[0], track_theta[1]);
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
                  if(tracks_vertex){
                     h_verxy_cut12c_ela -> Fill(vtx, vty);
                     h_verz_cut12c_ela -> Fill(vtz);
                  }
                  else {
                     std::cout << "  No vertex event:" << i << " (12c12c elastic)" << std::endl;
                     n_no_vertex ++;
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
                     h_range_thetalab_gsgs -> Fill(track_theta[0], track_range[0]);
                     h_range_thetalab_gsgs -> Fill(track_theta[1], track_range[1]);
                     h_kineE_thetalab_gsgs -> Fill(track_theta[0], track_KinE[0]);
                     h_kineE_thetalab_gsgs -> Fill(track_theta[1], track_KinE[1]);            
                     if(tracks_vertex){
                        h_verxy_gsgs -> Fill(vtx, vty);
                        h_verz_gsgs -> Fill(vtz);
#ifdef states_vertex
                        h_verz_gsgs_index[vindex] -> Fill(vtz);
                        h_verxy_gsgs_index[vindex] -> Fill(vtx, vty);
                        h_thetalab_thetalab_gsgs_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                        h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_gsgs_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
#endif
                     }
                  } 
                  //                  else if (sum_theta >= 84.0 && sum_theta < 88.0){
                  else if (theta_gsex->IsInside(track_theta[0],track_theta[1])){
                     if (runNum == 52){
                        peak2.push_back(i);
                     }
                     h_sum_theta_gsex -> Fill(sum_theta);
                     h_range_thetalab_gsex -> Fill(track_theta[0], track_range[0]);
                     h_range_thetalab_gsex -> Fill(track_theta[1], track_range[1]);
                     h_kineE_thetalab_gsex -> Fill(track_theta[0], track_KinE[0]);
                     h_kineE_thetalab_gsex -> Fill(track_theta[1], track_KinE[1]);            
                     if(tracks_vertex){
                        h_verxy_gsex -> Fill(vtx, vty);
                        h_verz_gsex -> Fill(vtz);
#ifdef states_vertex
                        h_verz_gsex_index[vindex] -> Fill(vtz);
                        h_verxy_gsex_index[vindex] -> Fill(vtx, vty);
                        h_thetalab_thetalab_gsex_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                        h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_gsex_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
#endif
                     }
                  }
                  //                  else if (sum_theta > 88.0 && sum_theta < 94.0){
                  else if (theta_exex->IsInside(track_theta[0],track_theta[1])){
                     if (runNum == 52){
                        peak3.push_back(i);
                     }
                     h_sum_theta_exex -> Fill(sum_theta);
                     h_range_thetalab_exex -> Fill(track_theta[0], track_range[0]);
                     h_range_thetalab_exex -> Fill(track_theta[1], track_range[1]);
                     h_kineE_thetalab_exex -> Fill(track_theta[0], track_KinE[0]);
                     h_kineE_thetalab_exex -> Fill(track_theta[1], track_KinE[1]);            
                     if(tracks_vertex){
                        h_verxy_exex -> Fill(vtx, vty);
                        h_verz_exex -> Fill(vtz);
#ifdef states_vertex
                        h_verz_exex_index[vindex] -> Fill(vtz);
                        h_verxy_exex_index[vindex] -> Fill(vtx, vty);
                        h_thetalab_thetalab_exex_index[vindex] -> Fill(track_theta[0], track_theta[1]);
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[0], track_KinE[0]);
                        h_kineE_thetalab_exex_index[vindex] -> Fill(track_theta[1], track_KinE[1]);
#endif
                     }
                  }
               }
            }
         }
         for (Int_t k = 0; k < ntrack; k++){
            if (track_alpha[k]){
               nalpha ++;
               alpha_tracks = true;
               //                  std::cout << " double check. event num:" << i << ", tracks: " << ntrack << ", alpha track: " << k << std::endl;
            }
            if (track_proton[k]){
               nproton ++;
               proton_tracks = true;
               //                  std::cout << " double check. event num:" << i << ", tracks: " << ntrack << ", proton track: " << k << std::endl;
            }
         }
         h_nalp_ntra->Fill(ntrack, nalpha);
         h_npro_ntra->Fill(ntrack, nproton);
         if(i%500==0){
            std::cout << "  Filling data: " << 100*i/nUnpackEvents << " %!    \r" << std::flush;
         }
      }
      //      std::cout << "  Filled data: 100 %!, " << nUnpackEvents << " events" << std::endl;
      //      std::cout << "Number of 2 tracks events in run" << runNum << ":" << nEventsWith2Tracks << std::endl;
      // Close files.
      unpackFile->Close();
   }

   // Draw histograms in TCanvas.
   //  Reset canvas c1
   if (gROOT->FindObject("c1")){
      delete gROOT->FindObject("c1");
   }
   //  set color
   kine_12c12c_exex_60_7 -> SetLineColor(kRed);
   xy90->SetLineColor(kRed);
   xy90->SetLineWidth(1);
   angle_12c12c_exex_60_7 -> SetLineColor(kRed);

#ifdef nom_check
   TCanvas *c1 = new TCanvas("c1", "c1");
   c1->cd();
   h_rmax->SetDirectory(0);
   h_rmax->GetXaxis()->SetTitle("Rmax [mm]");
   h_rmax->Draw();

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
   kine_12c12c_gsgs_60_7->Draw("same");

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
   xy90->Draw("same");

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
   xy90->Draw("same");

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
   kine_12c12c_gsgs_60_7->Draw("same");

   TCanvas *c19 = new TCanvas("c19", "c19");
   c19->cd();
   h_thetalab_thetalab_cut12c_ela->SetDirectory(0);
   h_thetalab_thetalab_cut12c_ela->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c )", (int)del_phi));
   h_thetalab_thetalab_cut12c_ela->Draw("colz");
   xy90->Draw("same");
   angle_12c12c_gsex_60_7->Draw("same");
   angle_12c12c_exex_60_7->Draw("same");
   theta_gsgs->Draw("same");
   theta_gsex->Draw("same");
   theta_exex->Draw("same");

   TCanvas *c20 = new TCanvas("c20", "c20");
   c20->cd();
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();

   TCanvas *c21 = new TCanvas("c21", "c21");
   c21->cd();
   h_nalp_ntra->SetDirectory(0);
   h_nalp_ntra->GetXaxis()->SetTitle("number of tracks");
   h_nalp_ntra->GetYaxis()->SetTitle("number of alphas tracks");
   h_nalp_ntra->SetTitle(Form("n_tracks n_alpha_tracks"));
   gPad->SetLogz();
   h_nalp_ntra->Draw("colz");

   TCanvas *c22 = new TCanvas("c22", "c22");
   c22->cd();
   h_range_thetalab_cutalpha->SetDirectory(0);
   h_range_thetalab_cutalpha->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cutalpha->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cutalpha->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d, alpha)", (int)del_phi));
   h_range_thetalab_cutalpha->Draw("colz");

   TCanvas *c23 = new TCanvas("c23", "c23");
   c23->cd();
   h_kineE_thetalab_alpha->SetDirectory(0);
   h_kineE_thetalab_alpha->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_alpha->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_alpha->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, alpha)", (int)del_phi));
   h_kineE_thetalab_alpha->Draw("colz");

   TCanvas *c24 = new TCanvas("c24", "c24");
   c24->cd();
   h_npro_ntra->SetDirectory(0);
   h_npro_ntra->GetXaxis()->SetTitle("number of tracks");
   h_npro_ntra->GetYaxis()->SetTitle("number of protons tracks");
   h_npro_ntra->SetTitle(Form("n_tracks n_proton_tracks"));
   gPad->SetLogz();
   h_npro_ntra->Draw("colz");

   TCanvas *c25 = new TCanvas("c25", "c25");
   c25->cd();
   h_range_thetalab_cutproton->SetDirectory(0);
   h_range_thetalab_cutproton->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cutproton->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cutproton->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d, proton)", (int)del_phi));
   h_range_thetalab_cutproton->Draw("colz");

   TCanvas *c26 = new TCanvas("c26", "c26");
   c26->cd();
   h_kineE_thetalab_proton->SetDirectory(0);
   h_kineE_thetalab_proton->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_proton->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_proton->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, proton)", (int)del_phi));
   h_kineE_thetalab_proton->Draw("colz");

   TCanvas *c27 = new TCanvas("c27", "c27", 600, 1000);
   c27->Divide(2,4);
   c27->cd(1);
   h_ntraver_ntra->SetDirectory(0);
   h_ntraver_ntra->GetXaxis()->SetTitle("number of tracks");
   h_ntraver_ntra->GetYaxis()->SetTitle("number of tracks belong to vertex");
   h_ntraver_ntra->SetTitle(Form("tracks with vertex"));
   gPad->SetLogz();
   h_ntraver_ntra->Draw("colz");
   c27->cd(2);
   h_ntra_verz->SetDirectory(0);
   h_ntra_verz->GetXaxis()->SetTitle("vertex z [mm]");
   h_ntra_verz->GetYaxis()->SetTitle("number of tracks");
   h_ntra_verz->SetTitle(Form("vertex z"));
   gPad->SetLogz();
   h_ntra_verz->Draw("colz");
   c27->cd(3);
   h_verxy->SetDirectory(0);
   h_verxy->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy->SetTitle(Form("vertex xy"));
   gPad->SetLogz();
   h_verxy->Draw("colz");
   c27->cd(4);
   h_verz->SetDirectory(0);
   h_verz->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz->SetTitle(Form("vertex z"));
   gPad->SetLogy();
   h_verz->Draw();
   c27->cd(5);
   h_verxz->SetDirectory(0);
   h_verxz->GetXaxis()->SetTitle("vertex z [mm]");
   h_verxz->GetYaxis()->SetTitle("vertex x [mm]");
   h_verxz->SetTitle(Form("vertex zx"));
   gPad->SetLogz();
   h_verxz->Draw("colz");
   c27->cd(6);
   h_veryz->SetDirectory(0);
   h_veryz->GetXaxis()->SetTitle("vertex z [mm]");
   h_veryz->GetYaxis()->SetTitle("vertex y [mm]");
   h_veryz->SetTitle(Form("vertex zy"));
   gPad->SetLogz();
   h_veryz->Draw("colz");
   c27->cd(7);
   h_verxy_cut12c_ela->SetDirectory(0);
   h_verxy_cut12c_ela->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_cut12c_ela->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_cut12c_ela->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verxy_cut12c_ela->Draw("colz");
   c27->cd(8);
   h_verz_cut12c_ela->SetDirectory(0);
   h_verz_cut12c_ela->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_cut12c_ela->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verz_cut12c_ela->Draw();

#endif

#ifdef peak_check
   TCanvas *c80 = new TCanvas("c80", "c80");
   c80->Divide(2,1);
   c80->cd(1);
   h_sum_theta->SetDirectory(0);
   h_sum_theta->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta->SetTitle(Form("Sum Theta_LAB (track == 2, 12c12c, run 52)"));
   h_sum_theta->Draw();
   c80->cd(2);
   h_sum_theta_cut12c_run52->SetDirectory(0);
   h_sum_theta_cut12c_run52->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_run52->SetTitle(Form("Sum Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c, run 52)", (int)del_phi));
   h_sum_theta_cut12c_run52->Draw();

   TCanvas *c81 = new TCanvas("c81", "c81", 1200,600);
   c81->Divide(3,1);
   c81->cd(1);
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();
   c81->cd(2);
   h_sum_theta_gsgs->SetDirectory(0);
   h_sum_theta_gsgs->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c gsgs");
   h_sum_theta_gsgs->Draw();
   c81->cd(3);
   h_sum_theta_gsex->SetDirectory(0);
   h_sum_theta_gsex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (gsex)");
   h_sum_theta_gsex->Draw();
   c81->cd(4);
   h_sum_theta_exex->SetDirectory(0);
   h_sum_theta_exex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (exex)");
   h_sum_theta_exex->Draw();

#endif

#ifdef c12_check
   TCanvas *c90 = new TCanvas("c90", "c90",1200,1000);
   h_thetalab_thetalab_cut12c_ela->SetDirectory(0);
   h_thetalab_thetalab_cut12c_ela->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cut12c_ela->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d, track == 2, 12c12c )", (int)del_phi));
   gPad->SetLogz();
   h_thetalab_thetalab_cut12c_ela->SetMinimum(1);
   h_thetalab_thetalab_cut12c_ela->Draw("colz");
   xy90->Draw("same");
   angle_12c12c_gsex_60_7->Draw("same");
   angle_12c12c_exex_60_7->Draw("same");
   theta_gsgs->Draw("same");
   theta_gsex->Draw("same");
   theta_exex->Draw("same");

   TCanvas *c91 = new TCanvas("c91", "c91", 1200, 1000);
   c91->Divide(4,4);
   c91->cd(1);
   h_sum_theta_cut12c_ela->SetDirectory(0);
   h_sum_theta_cut12c_ela->GetXaxis()->SetTitle("Sum of tracks [deg]");
   h_sum_theta_cut12c_ela->Draw();
   c91->cd(2);
   h_sum_theta_gsgs->SetDirectory(0);
   h_sum_theta_gsgs->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c gsgs");
   h_sum_theta_gsgs->Draw();
   c91->cd(3);
   h_sum_theta_gsex->SetDirectory(0);
   h_sum_theta_gsex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (gsex)");
   h_sum_theta_gsex->Draw();
   c91->cd(4);
   h_sum_theta_exex->SetDirectory(0);
   h_sum_theta_exex->GetXaxis()->SetTitle("Sum of tracks [deg] cut 12c (exex)");
   h_sum_theta_exex->Draw();
   c91->cd(5);
   h_kineE_thetalab_carbon->SetDirectory(0);
   h_kineE_thetalab_carbon->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_carbon->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_carbon->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_carbon->SetMinimum(1);
   h_kineE_thetalab_carbon->Draw("colz");
   kine_12c12c_gsgs_60_7->Draw("same");
   kine_12c12c_gsex_60_7->Draw("same");
   kine_12c12c_exex_60_7->Draw("same");
   c91->cd(6);
   h_kineE_thetalab_gsgs->SetDirectory(0);
   h_kineE_thetalab_gsgs->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_gsgs->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_gsgs->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, gsgs)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_gsgs->SetMinimum(1);
   h_kineE_thetalab_gsgs->Draw("colz");
   kine_12c12c_gsgs_60_7->Draw("same");
   kine_12c12c_gsex_60_7->Draw("same");
   kine_12c12c_exex_60_7->Draw("same");
   c91->cd(7);
   h_kineE_thetalab_gsex->SetDirectory(0);
   h_kineE_thetalab_gsex->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_gsex->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_gsex->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, gsex)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_gsex->SetMinimum(1);
   h_kineE_thetalab_gsex->Draw("colz");
   kine_12c12c_gsgs_60_7->Draw("same");
   kine_12c12c_gsex_60_7->Draw("same");
   kine_12c12c_exex_60_7->Draw("same");
   c91->cd(8);
   h_kineE_thetalab_exex->SetDirectory(0);
   h_kineE_thetalab_exex->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_kineE_thetalab_exex->GetYaxis()->SetTitle("roughKineE [MeV]");
   h_kineE_thetalab_exex->SetTitle(Form("KinE Theta_LAB (phi1-phi2-180 < %d, track ==2, 12c12c, exex)", (int)del_phi));
   gPad->SetLogz();
   h_kineE_thetalab_exex->SetMinimum(1);
   h_kineE_thetalab_exex->Draw("colz");
   kine_12c12c_gsgs_60_7->Draw("same");
   kine_12c12c_gsex_60_7->Draw("same");
   kine_12c12c_exex_60_7->Draw("same");
   c91->cd(9);
   h_verxy_cut12c_ela->SetDirectory(0);
   h_verxy_cut12c_ela->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_cut12c_ela->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_cut12c_ela->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verxy_cut12c_ela->Draw("colz");
   c91->cd(10);
   h_verxy_gsgs->SetDirectory(0);
   h_verxy_gsgs->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_gsgs->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_gsgs->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   h_verxy_gsgs->Draw("colz");
   c91->cd(11);
   h_verxy_gsex->SetDirectory(0);
   h_verxy_gsex->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_gsex->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_gsex->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, gsex)", (int)del_phi));
   h_verxy_gsex->Draw("colz");
   c91->cd(12);
   h_verxy_exex->SetDirectory(0);
   h_verxy_exex->GetXaxis()->SetTitle("vertex x [mm]");
   h_verxy_exex->GetYaxis()->SetTitle("vertex y [mm]");
   h_verxy_exex->SetTitle(Form("vertex xy (phi1-phi2-180 < %d, track == 2, 12c12c, exex)", (int)del_phi));
   h_verxy_exex->Draw("colz");
   c91->cd(13);
   h_verz_cut12c_ela->SetDirectory(0);
   h_verz_cut12c_ela->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_cut12c_ela->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c)", (int)del_phi));
   h_verz_cut12c_ela->Draw();
   c91->cd(14);
   h_verz_gsgs->SetDirectory(0);
   h_verz_gsgs->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsgs->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, gsgs)", (int)del_phi));
   h_verz_gsgs->Draw();
   c91->cd(15);
   h_verz_gsex->SetDirectory(0);
   h_verz_gsex->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_gsex->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, gsex)", (int)del_phi));
   h_verz_gsex->Draw();
   c91->cd(16);
   h_verz_exex->SetDirectory(0);
   h_verz_exex->GetXaxis()->SetTitle("vertex z [mm]");
   h_verz_exex->SetTitle(Form("vertex z (phi1-phi2-180 < %d, track == 2, 12c12c, exex)", (int)del_phi));
   h_verz_exex->Draw();

#endif

#ifdef states_vertex
   draw_ind("c92", "gsgs", n_group, n_h_z, h_verz_gsgs, h_verxy_gsgs, h_kineE_thetalab_gsgs, h_thetalab_thetalab_cut12c_ela,
      h_verz_gsgs_index, h_verxy_gsgs_index, h_kineE_thetalab_gsgs_index, h_thetalab_thetalab_gsgs_index,
      kine_12c12c_gsgs_60_7, kine_12c12c_gsex_60_7, kine_12c12c_exex_60_7, xy90, angle_12c12c_gsex_60_7, angle_12c12c_exex_60_7
   );

   draw_ind("c93", "gsex", n_group, n_h_z, h_verz_gsex, h_verxy_gsex, h_kineE_thetalab_gsex, h_thetalab_thetalab_cut12c_ela,
      h_verz_gsex_index, h_verxy_gsex_index, h_kineE_thetalab_gsex_index, h_thetalab_thetalab_gsex_index,
      kine_12c12c_gsgs_60_7, kine_12c12c_gsex_60_7, kine_12c12c_exex_60_7, xy90, angle_12c12c_gsex_60_7, angle_12c12c_exex_60_7
   );

   draw_ind("c94", "exex", n_group, n_h_z, h_verz_exex, h_verxy_exex, h_kineE_thetalab_exex, h_thetalab_thetalab_cut12c_ela,
      h_verz_exex_index, h_verxy_exex_index, h_kineE_thetalab_exex_index, h_thetalab_thetalab_exex_index,
      kine_12c12c_gsgs_60_7, kine_12c12c_gsex_60_7, kine_12c12c_exex_60_7, xy90, angle_12c12c_gsex_60_7, angle_12c12c_exex_60_7
   );

#endif

   /*
   TCanvas *c29 = new TCanvas("c29", "c29");
   c29->cd();
   h_philab_philab_cutphi_12c->SetDirectory(0);
   h_philab_philab_cutphi_12c->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_12c->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_12c->SetTitle(Form("Phi_LAB Phi_LAB (phi1-phi2-180 < %d, 12c)", (int)del_phi));
   h_philab_philab_cutphi_12c->Draw("colz");

   TCanvas *c30 = new TCanvas("c30", "c30");
   c30->cd();
   h_philab_philab_cutphi_alpha->SetDirectory(0);
   h_philab_philab_cutphi_alpha->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_alpha->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_alpha->SetTitle(Form("Phi_LAB Phi_LAB (phi1-phi2-180 < %d, alpha)", (int)del_phi));
   h_philab_philab_cutphi_alpha->Draw("colz");

   TCanvas *c31 = new TCanvas("c31", "c31");
   c31->cd();
   h_philab_philab_cutphi_proton->SetDirectory(0);
   h_philab_philab_cutphi_proton->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_proton->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_cutphi_proton->SetTitle(Form("Phi_LAB Phi_LAB (phi1-phi2-180 < %d, proton)", (int)del_phi));
   h_philab_philab_cutphi_proton->Draw("colz");

   TCanvas *c30 = new TCanvas("c30", "c30");
   c30->Divide(2,1);
   c30->cd(1);
   h_dEdx_range_cutmul1->SetDirectory(0);
   h_dEdx_range_cutmul1->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_dEdx_range_cutmul1->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_dEdx_range_cutmul1->Draw("colz");
   c30->cd(2);
   h_dEdx_range_cutmul2->SetDirectory(0);
   h_dEdx_range_cutmul2->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_dEdx_range_cutmul2->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_dEdx_range_cutmul2->Draw("colz");

   */

   /*
   TCanvas *c26 = new TCanvas("c26", "c26");
   c26->cd();
   h_dEdx_range_cutphi_12c->SetDirectory(0);
   h_dEdx_range_cutphi_12c->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range_cutphi_12c->GetYaxis()->SetTitle("dE/dx [ADC/mm]");
   h_dEdx_range_cutphi_12c->SetTitle(Form("dE/dx Range (phi1-phi2-180 < %d, 12c)", (int)del_phi));
   h_dEdx_range_cutphi_12c->Draw("colz");

   TCanvas *c27 = new TCanvas("c27", "c27");
   c27->cd();
   h_dEdx_range_cutphi_alpha->SetDirectory(0);
   h_dEdx_range_cutphi_alpha->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range_cutphi_alpha->GetYaxis()->SetTitle("dE/dx [ADC/mm]");
   h_dEdx_range_cutphi_alpha->SetTitle(Form("dE/dx Range (phi1-phi2-180 < %d, alpha)", (int)del_phi));
   h_dEdx_range_cutphi_alpha->Draw("colz");

   TCanvas *c28 = new TCanvas("c28", "c28");
   c28->cd();
   h_dEdx_range_cutphi_proton->SetDirectory(0);
   h_dEdx_range_cutphi_proton->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range_cutphi_proton->GetYaxis()->SetTitle("dE/dx [ADC/mm]");
   h_dEdx_range_cutphi_proton->SetTitle(Form("dE/dx Range (phi1-phi2-180 < %d, alpha)", (int)del_phi));
   h_dEdx_range_cutphi_proton->Draw("colz");

   TCanvas *c8 = new TCanvas();
   histExdp->SetDirectory(0);
   histExdp->Draw();
   histExdp->GetXaxis()->SetTitle("Ex [MeV]");

   TCanvas *c9 = new TCanvas();
   histExdd->SetDirectory(0);
   histExdd->Draw();
   histExdd->GetXaxis()->SetTitle("Ex [MeV]");
   */

   // Saving histograms in a .root file ...
   Results->cd();
   h_rmax->Write();
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
   h_dEdx_range_cutphi_alpha->Write();
   h_dEdx_range_cutphi_proton->Write();
   h_dEdx_range_cutmul1->Write();
   h_dEdx_range_cutmul2->Write();
   //   cutPIDproton->Write("PIDCutProton");
   //   cutPIDproton_extension->Write("PIDCutProtonExtension");
   //   cutPIDdeuteron->Write("PIDCutDeuteron");

   // Kinematics
   h_kineE_thetalab->Write();
   h_kineE_thetalab_carbon->Write();
   h_kineE_thetalab_alpha->Write();
   h_kineE_thetalab_proton->Write();
   h_kineE_thetalab_gsgs->Write();
   h_kineE_thetalab_gsex->Write();
   h_kineE_thetalab_exex->Write();
   /*
   kine_dp_gs->Write("kin_dp_gs");
   kine_dd_gs->Write("kin_dd_gs");
   kine_dd_gs_25MeVu->Write("kin_dd_gs_25MeVu");
   kine_dp_gs_25MeVu->Write("kin_dp_gs_25MeVu");
   */

   // angle correlations
   //  range vs theta
   h_range_thetalab->Write();
   h_range_thetalab_cutphi->Write();
   h_range_thetalab_cutphi_2tra->Write();
   h_range_thetalab_cut12c_ela->Write();
   h_range_thetalab_cutalpha->Write();
   h_range_thetalab_cutproton->Write();
   h_range_thetalab_gsgs->Write();
   h_range_thetalab_gsex->Write();
   h_range_thetalab_exex->Write();
   //  theta vs theta
   h_thetalab_thetalab_cutphi->Write();
   h_thetalab_thetalab_cutphi_2tra->Write();
   h_thetalab_thetalab_cut12c_ela->Write();
   //  phi vs phi
   h_philab_philab->Write();
   h_philab_philab_cutphi->Write();
   h_philab_philab_cutphi_12c->Write();
   h_philab_philab_cutphi_alpha->Write();
   h_philab_philab_cutphi_proton->Write();

   // vertex of tracks
   h_verz->Write();
   h_verz_cut12c_ela->Write();
   h_verz_gsgs->Write();
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

   /*
   // Excitation energy spectra
   histExdp->Write();
   histExdp_extended->Write();
   histExdp_CarbonSi->Write();
   histExdd->Write();
   // Angular distributions 
   histAngDist_elastic->Write();
   histAngDist_dp->Write();
   histAngDist_dp_CarbonSi->Write();
   */

   // Others ...
   h_nalp_ntra->Write();
   h_npro_ntra->Write();
   /*
   histESmallVTotalRange->Write();
   histEBigVBigRange->Write();
   histEstimatedKinEVThetaLAB2H->Write();
   histEstimatedKinEVThetaLAB1H->Write();
   */
   // save cuts
   cut12c->Write("charge_range_12c");
   cutalpha->Write("charge_range_alpha");
   cutproton->Write("charge_range_proton");
   theta_gsgs->Write("theta_theta_gsgs");
   theta_gsex->Write("theta_theta_gsex");
   theta_exex->Write("theta_theta_exex");
   // save lines
   kine_12c12c_gsgs_60_7->Write("kine_12c_gsgs");
   kine_12c12c_gsex_60_7->Write("kine_12c_gsex");
   kine_12c12c_exex_60_7->Write("kine_12c_exex");
   xy90->Write("angle_12c_gsgs");
   angle_12c12c_gsex_60_7->Write("angle_12c_gsex");
   angle_12c12c_exex_60_7->Write("angle_12c_exex");
   Results->Close();

   // cout of information
   std::cout << "                                                                " << std::endl;
   std::cout << "12c12c elastic events with no vertex is " << n_no_vertex << " events." << std::endl;
   //   std::cout << "Maximum radius of hits: " << max_r_max << " mm, Trigger radius: " << r_tri << " mm" << std::endl;
#ifdef peak_check
   std::cout << "peak1 events: " << peak1.size() << std::endl;
   std::cout << "  Event of peak1 " << std::endl;
   for (auto &eventIndex: peak1){
      std::cout << eventIndex << ", " << std::flush;
   }
   std::cout << std::endl << std::endl;
   std::cout << "peak2 events: " << peak2.size() << std::endl;
   std::cout << "  Event of peak2: " << std::endl;
   for (auto &eventIndex: peak2){
      std::cout << eventIndex << ", " << std::flush;
   }
   std::cout << std::endl << std::endl;
   std::cout << "peak3 events: " << peak3.size() << std::endl;
   std::cout << "  Event of peak3: " << std::endl;
   for (auto &eventIndex: peak3){
      std::cout << eventIndex << ", " << std::flush;
   }

#endif
   std::cout << std::endl;

   // stop timer
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl;
   cout << "Real time: " << rtime << " s, CPU time: " << ctime << " s" << endl << endl;

}

TGraph* ReadKinematics(TString kineFile){
   /*
   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];
   */

   std::vector<Double_t> ThetacmS, ThetalabR, ElabR, ThetalabS, ElabS, MomlabR;
   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()){
      /*
      while (!kineStr->eof()){
      //	      *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >> EnerLabSca[numKin];
	      *kineStr >> ThetaLabRec[numKin] >> EnerLabRec[numKin];
         numKin++;
      }
      */
      double tcm, tlabr, elabr, tlabs, elabs;
      while ( *kineStr >> tlabr >> elabr){
         ThetalabR.push_back(tlabr);
         ElabR.push_back(elabr);
      }
   }
   else {
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;
   }
   //   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   TGraph *kine = new TGraph(ThetalabR.size(), ThetalabR.data(), ElabR.data());
   //   kine->Sort();
   return kine;
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
   TGraph *kine_gsgs, TGraph *kine_gsex, TGraph *kine_exex, TF1 *ang_gsgs, TGraph *ang_gsex, TGraph *ang_exex)
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
}
