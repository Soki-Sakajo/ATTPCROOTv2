#include <fstream>
#include "TFile.h"
#include "TObject.h"
#include "TCanvas.h"

TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);

void kine(){
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

   // files.
   std::vector runNums = {52};
   //   std::vector runNums = {52,53,54,55,56,57,58,59};
   Int_t run_start = runNums.front();
   Int_t run_end = runNums.back();
   TFile * Results = new TFile(Form("data/kine_results_run%d-run%d.root", run_start, run_end),"recreate");
   std::ofstream Results_c(Form("can_output/kine_canvases_run%d-run%d.C", run_start, run_end),std::ios::out | std::ios::trunc);
   
   FairRunAna *run = new FairRunAna(); // Forcing a dummy run
   //   TString outfname="./canvas_kine.root";
   //   TFile *outfile=new TFile(outfname,"recreate");
   AtFindVertex *fver = new AtFindVertex();

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
   /*
   TFile *cutKineFile = new TFile("./cutFiles/cutKine.root", "READ");
   TCutG *cutKineBProton = (TCutG *)cutKineFile->Get("cutKineBProton");
   cutKineFile->Close();

   TFile *cutATTPCPIDFile = new TFile("./cutFiles/ATTPC_PID.root", "READ"); 
   TCutG *cutPIDproton = (TCutG *)cutATTPCPIDFile->Get("cutATTPCPIDproton");
   TCutG *cutPIDproton_extension = (TCutG *)cutATTPCPIDFile->Get("cutATTPCPIDproton_extension");
   TCutG *cutPIDdeuteron = (TCutG *)cutATTPCPIDFile->Get("cutATTPCPIDdeuteron");
   cutATTPCPIDFile->Close();
   */

   // Kinematic curve
   /*TGraph* kinecurve_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");
   kinecurve_3HeGS->SetLineWidth(2);
   kinecurve_3HeGS->SetLineColor(kRed);
   TGraph* ang_lab_cm_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");// 5 deg pitch in theta_cm
   ang_lab_cm_3HeGS->SetMarkerStyle(8);
   ang_lab_cm_3HeGS->SetMarkerSize(1);
   TLegend *legend = new TLegend(0.6,0.2,0.85,0.5);
   legend->AddEntry(kinecurve_3HeGS,"3He G.S.","l");
   legend->AddEntry(ang_lab_cm_3HeGS,"5 deg pitch in #theta_{cm}","p");
   legend->SetFillColor(0);*/

   // Characteristic definitions
   bool check_tracks = false; 
   Int_t narray = 10;
   Int_t ntrack = 0;
   Int_t itrack = 0;
   Double_t rad = -100;
   Double_t r_tem = 0;
   Double_t r_max = 0;
   Double_t r_tri = 300;
   Double_t max_r_max = 0;
   Double_t track_verx = 0;
   Double_t track_very = 0;
   Double_t track_verz = 0;
   Double_t track_theta[narray];
   Double_t track_phi[narray];
   Double_t track_range[narray];
   Double_t track_charge[narray];
   Double_t track_r[narray];

   // Histogram definitions.
   // ... TH1 hist for checking something.
   TH1D *h_ntra = new TH1D("h_ntra", "h_ntra;NTracks", 10, 0, 10);
   TH1D *h_rmax = new TH1D("h_rmax", "h_rmax;Rmax [mm]", 100, 0, 200);

   // ... ATTPC PID
   TH2F *h_charge_range = new TH2F("h_charge_range", "h_charge_range;roughRange [mm];Charge [ADC]", 600, 0, 1200, 600, 0, 6e5);
   TH2F *h_charge_range_cutphi = new TH2F("h_charge_range_cutphi", "h_charge_range_cutphi;roughRange [mm];Charge [ADC]", 300, 0, 1200, 300, 0, 6e5);
   TH2F *h_dEdx_range = new TH2F("h_dEdx_range", "h_dEdx_range;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);
   TH2F *h_dEdx_range_backwards = new TH2F("h_dEdx_range_backwards", "h_dEdx_range_backwards;roughRange [mm];dEdx [ADC/mm]", 515, 0, 1030, 2000, 0, 4000);

   // ... kinematics 
   TH2F *h_kinE_thetalab = new TH2F("h_kinE_thetalab", "h_kinE_thetalab;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);
   TH2F *h_kinE_thetalab_carbon = new TH2F("h_kinE_thetalab_carbon", "h_kinE_thetalab_carbon;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);
   TH2F *h_kinE_thetalab_alpha = new TH2F("h_kinE_thetalab_alpha", "h_kinE_thetalab_alpha;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);

   // ... angle correlations
   TH2F *h_thetalab_thetalab_cutphi = new TH2F("h_thetalab_thetalab_cutphi", "h_thetalab_thetalab_cutphi", 200, 0, 100, 200, 0, 100);
   TH2F *h_philab_philab_cutphi = new TH2F("h_philab_philab_cutphi", "h_philab_philab_cutphi", 360, -180, 180, 360, -180, 180);
   TH2F *h_range_thetalab = new TH2F("h_range_thetalab", "h_range_thetalab", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutphi = new TH2F("h_range_thetalab_cutphi", "h_range_thetalab_cutphi", 180, 0, 180, 1030, 0, 1030);
   TH2F *h_range_thetalab_cutverz = new TH2F("h_range_thetalab_cutverz", "h_range_thetalab_cutverz", 180, 0, 180, 1030, 0, 1030);

   // ... track vertex
   TH1D *h_verz = new TH1D("h_verz", "h_verz;Vertex Z [mm]", 1500, -500, 1000);
   TH2F *h_verxy = new TH2F("h_verxy", "h_verxy", 400, -200, 200, 400, -200, 200);
   TH2F *h_verz_ntra = new TH2F("h_verz_ntra", "h_verz_ntra", 10, 0, 10, 750, -500, 1000);

   
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
   TH2F *h_kinE_thetalab_2H  = new TH2F("h_kinE_thetalab_2H", "h_kinE_thetalab_2H", 180, 0, 180, 250, 0, 20);
   TH2F *h_kinE_thetalab_1H  = new TH2F("h_kinE_thetalab_1H", "h_kinE_thetalab_1H", 180, 0, 180, 250, 0, 20);

   for (int runNum: runNums) {
      // Open the digitalization file and get the TTree.
      TString unpackFileName = TString::Format("./decode_data/run_%04d.root", runNum);
      TFile *unpackFile = new TFile(unpackFileName, "READ");
      TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
      int nUnpackEvents = unpackTree->GetEntries();
      std::cout << " Number of unpacked events in run " << runNum << ": " << nUnpackEvents << std::endl;
      int nEventsWith2Tracks = 0;
      // Creare the TTreeReader to read the AtTrackingEvents and simulation.
      TTreeReader unpackReader("cbmsim", unpackFile);
      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");

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
         rad = 0;
         r_max = 0;
         r_tem = 0;
         itrack = 0;
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

            track_range[itrack] = pattern->DistanceAlongPattern(lastPoint, firstPoint);
            track_charge[itrack] = track.GetGeoQEnergy();
            track_theta[itrack] = 180 - track.GetGeoTheta() * 180 / TMath::Pi();
            track_phi[itrack] = track.GetGeoPhi() * 180 / TMath::Pi();

            for (auto &hit: hits) {
               auto pos = hit->GetPosition();
               int padNum = hit->GetPadNum();
               int sizeID = map->GetPadSize(padNum);
               rad = TMath::Sqrt(pos.X() * pos.X() + pos.Y() * pos.Y());
               if (sizeID == 1) {
                  bigPadCharge += hit->GetCharge();
                  continue;
               }
               smallPadCharge += hit->GetCharge();
               double currentRangeInSmallPads = pattern->DistanceAlongPattern(hit->GetPosition(), firstPoint);
               if (currentRangeInSmallPads > rangeInSmallPads){
                  rangeInSmallPads = currentRangeInSmallPads;
               }
               if (rad > r_tem){
                  r_tem = rad;
               }
            }
            track_r[itrack] = r_tem;
            if(track_r[itrack] > r_max){
               r_max = track_r[itrack];
            }

            double dEdx = smallPadCharge / rangeInSmallPads;
            double rangeInBigPads = track_range[itrack] - rangeInSmallPads;
            bool reachedBigPads = true;
            if (rangeInBigPads / track_range[itrack] < 0.05){
               reachedBigPads = false;
            }

            h_charge_range->Fill(track_range[itrack],track_charge[itrack]);
            h_range_thetalab->Fill(track_theta[itrack], track_range[itrack]);
            h_dEdx_range->Fill(track_range[itrack], dEdx);
            if (!reachedBigPads){
               h_Esmall_range->Fill(track_range[itrack], smallPadCharge);
            }
            else {
               h_Ebig_rangebig->Fill(rangeInBigPads, bigPadCharge);
            }

            if (track_theta[itrack] > 100){
               h_dEdx_range_backwards->Fill(track_range[itrack], dEdx);
            }

            double estimatedKinE{0.1};
            /*
            if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx)) {
               while (eLossModelC3D8_d->GetRange(estimatedKinE) < roughRangeEstimation)
                  estimatedKinE += 0.01;
            } else if (cutPIDproton->IsInside(roughRangeEstimation, dEdx) or cutPIDproton_extension->IsInside(roughRangeEstimation, dEdx)) {
               while (eLossModelC3D8_p->GetRange(estimatedKinE) < roughRangeEstimation)
                  estimatedKinE += 0.01;
            } else {
               estimatedKinE = -1;
            }

            // While we don't have cut files defined, assume deuteron.
            //while (eLossModelC3D8_d->GetRange(estimatedKinE) < roughRangeEstimation)
               //estimatedKinE += 0.01;

            h_kinE_thetalab->Fill(trackThetaLAB, estimatedKinE);
            */

            /*
            // Kinematics plots
            // .... protons in ATTPC
            if (cutPIDproton->IsInside(roughRangeEstimation, dEdx)) {
               histSiPIDADCMax_ProtonATTPC->Fill(maxADCFront2, maxADCFront1);
               histEstimatedKinEVThetaLAB_ProtonATTPC->Fill(trackThetaLAB, estimatedKinE);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_p, m_18C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), estimatedKinE);
               histExdp->Fill(Ex);            
            }
            if (cutPIDproton->IsInside(roughRangeEstimation, dEdx) or cutPIDproton_extension->IsInside(roughRangeEstimation, dEdx)) {
               histEstimatedKinEVThetaLAB_ProtonATTPC_extended->Fill(trackThetaLAB, estimatedKinE);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_p, m_18C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), estimatedKinE);
               histExdp_extended->Fill(Ex);

               histAngDist_dp->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));

               // Also add condition of Carbon in Si
               if (cutSiC->IsInside(maxADCFront2, maxADCFront1)) {
                  histEstimatedKinEVThetaLAB_ProtonATTPC_CarbonSi->Fill(trackThetaLAB, estimatedKinE);
                  histExdp_CarbonSi->Fill(Ex);
                  histAngDist_dp_CarbonSi->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));
               }
            }
            // .... deuterons in ATTPC
            if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx)) {
               histSiPIDADCMax_DeuteronATTPC->Fill(maxADCFront2, maxADCFront1);
               histEstimatedKinEVThetaLAB_DeuteronATTPC->Fill(trackThetaLAB, estimatedKinE);

               auto [Ex, thetaCM] = kine_2b(m_17C, m_d, m_d, m_17C, eLossModelC3D8_17C->GetEnergy(E_beam, 1000 - pseudoVertex.z() * 1000), trackThetaLAB * TMath::DegToRad(), estimatedKinE);
               histExdd->Fill(Ex);

               if ((Ex >= -4.75) && (Ex <= 4.73)) {
                  histAngDist_elastic->Fill(thetaCM, 1 / TMath::Sin(thetaCM * TMath::DegToRad()));
               }
            }
            */
            itrack ++;
         }
         check_tracks = false;
         if (ntrack == 0){
            //            cout<<"No tracks found in this event. run: "<<runNum<<" , event: "<<i<<endl;
            continue;
         }
         else if (ntrack == 1){
            fver->FindvertexAingleLine(tracks);
            check_tracks = true;
         }
         else if (ntrack > 1 && ntrack < 7){
            fver->FindVertexMultipleLines(tracks, ntrack);
            check_tracks = true;
         }
         else {
            //            cout<<"More than 6 tracks found in this event. run: "<<runNum<<" , event: "<<i<<endl;
            continue;
         }
         auto vtxlist = fver->GetTracksVertex();
         for (auto &v: vtxlist){
            track_verx = v.vertex.X();
            track_very = v.vertex.Y();
            track_verz = v.vertex.Z();
         }
         h_verxy->Fill(track_verx, track_very);
         h_verz->Fill(track_verz);
         h_verz_ntra->Fill(ntrack, track_verz);

         h_rmax->Fill(r_max);
         if(r_tri > r_max && r_max > 0){
            r_tri = r_max;
         }
         if(r_max > max_r_max){
            max_r_max = r_max;
         }

         if(abs(abs(track_phi[0] - track_phi[1]) - 180) < del_phi){
            h_charge_range_cutphi -> Fill(track_range[0], track_charge[0]);
            h_thetalab_thetalab_cutphi -> Fill(track_theta[0], track_theta[1]);
            h_range_thetalab_cutphi -> Fill(track_theta[0], track_range[0]);
            h_philab_philab_cutphi -> Fill(track_phi[0], track_phi[1]);
         }
         if(i%100==0){
            std::cout << "  Filling data: " << 100*i/nUnpackEvents << " %!    \r" << std::flush;
         }
      }
      //      std::cout << "  Filled data: 100 %!, " << nUnpackEvents << " events" << std::endl;
      //      std::cout << "Number of 2 tracks events in run" << runNum << ":" << nEventsWith2Tracks << std::endl;
      // Close files.
      unpackFile->Close();
   }
   std::cout << "Maximum radius of hits: " << max_r_max << " mm, Trigger radius: " << r_tri << " mm" << std::endl;

   /*
   // Kinematic lines.
   TGraph *kine_dd_gs = ReadKinematics("./kineFiles/kine17C_dd_gs.txt");
   TGraph *kine_dp_gs = ReadKinematics("./kineFiles/kine17C_dp_gs.txt");

   TGraph *kine_dd_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dd_gs_25MeVu.txt");
   TGraph *kine_dp_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dp_gs_25MeVu.txt");
   */

   // Draw histograms in TCanvas.
   TCanvas *c1 = new TCanvas("c1", "c1");
   h_rmax->SetDirectory(0);
   h_rmax->Draw();
   h_rmax->GetXaxis()->SetTitle("Rmax [mm]");

   TCanvas *c2 = new TCanvas("c2", "c2");
   h_ntra->SetDirectory(0);
   h_ntra->Draw();
   h_ntra->GetXaxis()->SetTitle("Number of tracks");

   TCanvas *c3 = new TCanvas("c3", "c3");
   h_charge_range->SetDirectory(0);
   h_charge_range->Draw("colz");
   h_charge_range->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range->GetYaxis()->SetTitle("Charge [ADC]");

   TCanvas *c4 = new TCanvas("c4", "c4");
   h_range_thetalab->SetDirectory(0);
   h_range_thetalab->Draw("colz");
   h_range_thetalab->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab->GetYaxis()->SetTitle("roughRange [mm]");

   TCanvas *c5 = new TCanvas("c5", "c5");
   h_dEdx_range->SetDirectory(0);
   h_dEdx_range->Draw("colz");
   //   cutPIDproton->Draw("same");
   //   cutPIDdeuteron->Draw("same");
   h_dEdx_range->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range->GetYaxis()->SetTitle("dEdx [ADC/mm]");

   TCanvas *c6 = new TCanvas("c6", "c6");
   h_thetalab_thetalab_cutphi->SetDirectory(0);
   h_thetalab_thetalab_cutphi->Draw("colz");
   //kine_d3He_tt->Draw("same");
   h_thetalab_thetalab_cutphi->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");
   h_thetalab_thetalab_cutphi->SetTitle(Form("Theta_LAB Theta_LAB (phi1-phi2-180 < %d )", (int)del_phi));

   TCanvas *c7 = new TCanvas("c7", "c7");
   h_range_thetalab_cutphi->SetDirectory(0);
   h_range_thetalab_cutphi->Draw("colz");
   h_range_thetalab_cutphi->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   h_range_thetalab_cutphi->GetYaxis()->SetTitle("roughRange [mm]");
   h_range_thetalab_cutphi->SetTitle(Form("Range Theta_LAB (phi1-phi2-180 < %d )", (int)del_phi));

   TCanvas *c8 = new TCanvas("c8", "c8");
   h_philab_philab_cutphi->SetDirectory(0);
   h_philab_philab_cutphi->Draw("colz");
   //kine_d3He_tt->Draw("same");
   h_philab_philab_cutphi->GetXaxis()->SetTitle("track1_#phi_{LAB} [deg]");
   h_philab_philab_cutphi->GetYaxis()->SetTitle("track2_#phi_{LAB} [deg]");
   h_philab_philab_cutphi->SetTitle(Form("Phi_LAB Phi_LAB (phi1-phi2-180 < %d )", (int)del_phi));

   TCanvas *c9 = new TCanvas("c9", "c9");
   h_verxy->SetDirectory(0);
   h_verxy->Draw("colz");
   h_verxy->GetXaxis()->SetTitle("Vertex X [mm]");
   h_verxy->GetYaxis()->SetTitle("Vertex Y [mm]");

   TCanvas *c10 = new TCanvas("c10", "c10");
   h_verz->SetDirectory(0);
   h_verz->Draw();
   h_verz->GetXaxis()->SetTitle("Vertex Z [mm]");

   TCanvas *c11 = new TCanvas("c11", "c11");
   h_verz_ntra->SetDirectory(0);
   h_verz_ntra->Draw("colz");
   h_verz_ntra->GetXaxis()->SetTitle("Number of tracks");
   h_verz_ntra->GetYaxis()->SetTitle("Vertex Z [mm]");

   /*
   TCanvas *c6 = new TCanvas("c6", "c6");
   h_dEdx_range_backwards->SetDirectory(0);
   h_dEdx_range_backwards->Draw("colz");
   //   cutPIDproton->Draw("same");
   //   cutPIDdeuteron->Draw("same");
   h_dEdx_range_backwards->GetXaxis()->SetTitle("roughRange [mm]");
   h_dEdx_range_backwards->GetYaxis()->SetTitle("#frac{dE}{dx} [ADC/mm]");

   TCanvas *c7 = new TCanvas("c7", "c7");
   h_charge_range_cutphi->SetDirectory(0);
   h_charge_range_cutphi->Draw("colz");
   h_charge_range_cutphi->GetXaxis()->SetTitle("roughRange [mm]");
   h_charge_range_cutphi->GetYaxis()->SetTitle("Charge [ADC]");
   h_charge_range_cutphi->SetTitle(Form("Charge Range (phi1-phi2-180 < %d )", (int)del_phi));

   TCanvas *c5 = new TCanvas();
   histEstimatedKinEVThetaLABTotal->SetDirectory(0);
   histEstimatedKinEVThetaLABTotal->Draw("colz");
   //kine_d3He->Draw("same");
   //kine_d3HeEx2_2->Draw("same");
   //kine_d3HeEx2_7->Draw("same");
   histEstimatedKinEVThetaLABTotal->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   histEstimatedKinEVThetaLABTotal->GetYaxis()->SetTitle("roughKinE [MeV]");

   TCanvas *c6 = new TCanvas();
   histEstimatedKinEVThetaLABTotal->SetDirectory(0);
   histEstimatedKinEVThetaLABTotal->Draw("colz");
   //   kine_dp_gs->Draw("same");
   //   kine_dd_gs->Draw("same");
   histEstimatedKinEVThetaLABTotal->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   histEstimatedKinEVThetaLABTotal->GetYaxis()->SetTitle("roughKinE [MeV]");

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
   h_charge_range->Write();
   h_charge_range_cutphi->Write();
   // dEdx Vs Total Range
   //   h_dEdx_range->Write();
   //   h_dEdx_range_backwards->Write();
   //   cutPIDproton->Write("PIDCutProton");
   //   cutPIDproton_extension->Write("PIDCutProtonExtension");
   //   cutPIDdeuteron->Write("PIDCutDeuteron");

   // Kinematics
   /*
   h_kinE_thetalab->Write();
   histEstimatedKinEVThetaLAB_ProtonATTPC->Write();
   histEstimatedKinEVThetaLAB_ProtonATTPC_extended->Write();
   histEstimatedKinEVThetaLAB_DeuteronATTPC->Write();
   histEstimatedKinEVThetaLAB_CarbonSi->Write();
   histEstimatedKinEVThetaLAB_NitrogenSi->Write();
   kine_dp_gs->Write("kin_dp_gs");
   kine_dd_gs->Write("kin_dd_gs");
   kine_dd_gs_25MeVu->Write("kin_dd_gs_25MeVu");
   kine_dp_gs_25MeVu->Write("kin_dp_gs_25MeVu");
   */

   // angle correlations
   h_range_thetalab->Write();
   h_range_thetalab_cutphi->Write();
   h_thetalab_thetalab_cutphi->Write();
   h_philab_philab_cutphi->Write();

   // track vertex
   h_verxy->Write();
   h_verz->Write();
   h_verz_ntra->Write();

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
   /*
   histESmallVTotalRange->Write();
   histEBigVBigRange->Write();
   histEstimatedKinEVThetaLAB2H->Write();
   histEstimatedKinEVThetaLAB1H->Write();
   */
   Results->Close();

   // save canvases as .C macros
   Results_c << "void kine_canvases(){\n" << std::endl;
   TSeqCollection *canList = gROOT->GetListOfCanvases();
   TIter next(canList);
   TCanvas *c;
   Int_t can_num = 1;
   while ((c = (TCanvas *)next())) {
      Results_c << Form("TCanvas *c%d = new TCanvas(\"c%d\", \"c%d\");\n", can_num, can_num, can_num);
      Results_c << Form("c%d->cd();\n", can_num);
      c->SavePrimitive(Results_c, "");
      Results_c << "\n";
      can_num++;
   }
   Results_c << "}\n" << std::endl;
   Results_c.close();
}

TGraph* ReadKinematics(TString kineFile){
   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()){
      while (!kineStr->eof()){
	         *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >> EnerLabSca[numKin];
         numKin++;
      }
   } else if (kineStr->fail())
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
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
