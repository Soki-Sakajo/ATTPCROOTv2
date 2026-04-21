TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);

void kine(){
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

   // events.
   std::vector runNums = {52};
   TFile * Results = new TFile("data/kine_results.root","recreate");

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
   Double_t r_tri = 300;
   Double_t r_max = 0;

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
   
   // Histogram definitions.
   // ... ATTPC PID
   TH2F *histChargeTotalRange = new TH2F("histChargeTotalRange", "histChargeTotalRange;roughRange [mm];Charge [ADC]", 500, 0, 1030, 600, 0, 6e5);
   TH2F *histdEdxVTotalRange = new TH2F("histdEdxVTotalRange", "histdEdxVTotalRange;roughRange [mm];dEdx [ADC/mm]", 500, 0, 1030, 1600, 0, 4000);
   TH2F *histdEdxVTotalRangeBackwards = new TH2F("histdEdxVTotalRangeBackwards", "histdEdxVTotalRangeBackwards;roughRange [mm];dEdx [ADC/mm]", 500, 0, 1030, 1600, 0, 4000);
   // ... kinematics 
   TH2F *histEstimatedKinEVThetaLABTotal = new TH2F("histEstimatedKinEVThetaLABTotal", "histEstimatedKinEVThetaLABTotal;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB_Carbon = new TH2F("histEstimatedKinEVThetaLAB_Carbon", "histEstimatedKinEVThetaLAB_Carbon;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB_Alpha = new TH2F("histEstimatedKinEVThetaLAB_Alpha", "histEstimatedKinEVThetaLAB_Alpha;#theta_{LAB} [deg];roughKinE [MeV]", 180, 0, 180, 300, 0, 30);

   /*
   // ... Excitation energy 
   TH1F *histExdp = new TH1F("histExdp", "histExdp;Ex [MeV]", 80, -5, 15);
   TH1F *histExdp_extended = new TH1F("histExdp_extended", "histExdp_extended;Ex [MeV]", 80, -5, 15);
   TH1F *histExdp_CarbonSi = new TH1F("histExdp_CarbonSi", "histExdp_CarbonSi;Ex [MeV]", 80, -5, 15);
   TH1F *histExdd = new TH1F("histExdd", "histExdd;Ex [MeV]",80 , -5, 15);

   // ... center-of-mass angular distributions
   TH1F *histAngDist_elastic = new TH1F("histAngDist_elastic", "histAngDist_elastic;#theta_{c.m.} [deg]", 180, 0, 180);
   TH1F *histAngDist_dp = new TH1F("histAngDist_dp", "histAngDist_dp;#theta_{c.m.} [deg];Counts / deg", 180, 0, 180);
   TH1F *histAngDist_dp_CarbonSi = new TH1F("histAngDist_dp_CarbonSi", "histAngDist_dp_CarbonSi;#theta_{c.m.} [deg];Counts / deg", 180, 0, 180);
   */
   // ... others ...
   TH2F *histThetaLABThetaLAB = new TH2F("histThetaLABThetaLAB", "histThetaLABThetaLAB", 360, 0, 180, 360, 0, 180);
   TH2F *histESmallVTotalRange = new TH2F("histESmallVTotalRange", "histESmallVTotalRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *histEBigVBigRange = new TH2F("histEBigVBigRange", "histEBigVBigRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *histRangeVThetaLAB = new TH2F("histRangeVThetaLAB", "histRangeVThetaLAB", 180, 0, 180, 1030, 0, 1030);
   TH2F *histEstimatedKinEVThetaLAB2H  = new TH2F("histEstimatedKinEVThetaLAB2H", "histEstimatedKinEVThetaLAB2H", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB1H  = new TH2F("histEstimatedKinEVThetaLAB1H", "histEstimatedKinEVThetaLAB1H", 180, 0, 180, 250, 0, 20);

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
         if (!patternEvent) continue;

         // We want to focus on events with 2 or less tracks for now.
         auto &tracks = patternEvent->GetTrackCand();
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

            auto *pattern = track.GetPattern();
            auto firstPoint = track.GetFirstPoint();
            auto lastPoint = track.GetLastPoint();
            auto charge = track.GetGeoQEnergy();
            double roughRangeEstimation = pattern->DistanceAlongPattern(lastPoint, firstPoint);
            auto pseudoVertex = pattern->ClosestPointOnPattern(firstPoint);
            double trackThetaLAB = track.GetGeoTheta() * 180 / TMath::Pi();
            double trackPhi = track.GetGeoPhi() * 180 / TMath::Pi();
            double smallPadCharge{};
            double bigPadCharge{};
            auto braggCurvePairs = track.GetBraggCurveValues();
            auto &hits = track.GetHitArray();
            double rangeInSmallPads{};
            for (auto &hit: hits) {
               auto pos = hit->GetPosition();
               auto rad = TMath::Sqrt(pos.X() * pos.X() + pos.Y() * pos.Y());
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
               if(rad > r_max){
                  r_max = rad;
               }
            }
            double dEdx = smallPadCharge / rangeInSmallPads;

            double rangeInBigPads = roughRangeEstimation - rangeInSmallPads;
            bool reachedBigPads = true;
            if (rangeInBigPads / roughRangeEstimation < 0.05)
               reachedBigPads = false;

            histChargeTotalRange->Fill(roughRangeEstimation,charge);
            histRangeVThetaLAB->Fill(trackThetaLAB, roughRangeEstimation);
            histdEdxVTotalRange->Fill(roughRangeEstimation, dEdx);
            if (!reachedBigPads)
               histESmallVTotalRange->Fill(roughRangeEstimation, smallPadCharge);
            else
               histEBigVBigRange->Fill(rangeInBigPads, bigPadCharge);

            if (trackThetaLAB > 100)
               histdEdxVTotalRangeBackwards->Fill(roughRangeEstimation, dEdx);

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

            histEstimatedKinEVThetaLABTotal->Fill(trackThetaLAB, estimatedKinE);
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
         }
         if(r_tri > r_max){
            r_tri = r_max;
         }
         if(i%100==0){
            std::cout << "  Filling data: " << 100*i/nUnpackEvents << " %!    \r" << std::flush;
         }
      }
      std::cout << "  Filled data: 100 %!, " << nUnpackEvents << " events" << std::endl;
      //      std::cout << "Number of 2 tracks events in run" << runNum << ":" << nEventsWith2Tracks << std::endl;
      // Close files.
      unpackFile->Close();
   }
   std::cout << "Maximum radius of hits: " << r_max << " mm, Trigger radius: " << r_tri << " mm" << std::endl;

   /*
   // Kinematic lines.
   TGraph *kine_dd_gs = ReadKinematics("./kineFiles/kine17C_dd_gs.txt");
   TGraph *kine_dp_gs = ReadKinematics("./kineFiles/kine17C_dp_gs.txt");

   TGraph *kine_dd_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dd_gs_25MeVu.txt");
   TGraph *kine_dp_gs_25MeVu = ReadKinematics("./kineFiles/kine17C_dp_gs_25MeVu.txt");
   */

   // Draw histograms in TCanvas.
   TCanvas *c1 = new TCanvas();
   histChargeTotalRange->SetDirectory(0);
   histChargeTotalRange->Draw("colz");
   histChargeTotalRange->GetXaxis()->SetTitle("roughRange [mm]");
   histChargeTotalRange->GetYaxis()->SetTitle("Charge [ADC]");

   TCanvas *c2 = new TCanvas();
   histRangeVThetaLAB->SetDirectory(0);
   histRangeVThetaLAB->Draw("colz");
   histRangeVThetaLAB->GetXaxis()->SetTitle("#theta_{LAB} [deg]");
   histRangeVThetaLAB->GetYaxis()->SetTitle("roughRange [mm]");

   TCanvas *c3 = new TCanvas();
   histdEdxVTotalRange->SetDirectory(0);
   histdEdxVTotalRange->Draw("zcol");
   //   cutPIDproton->Draw("same");
   //   cutPIDdeuteron->Draw("same");
   histdEdxVTotalRange->GetXaxis()->SetTitle("roughRange [mm]");
   histdEdxVTotalRange->GetYaxis()->SetTitle("dEdx [ADC/mm]");

   TCanvas *c4 = new TCanvas();
   histdEdxVTotalRangeBackwards->SetDirectory(0);
   histdEdxVTotalRangeBackwards->Draw("zcol");
   //   cutPIDproton->Draw("same");
   //   cutPIDdeuteron->Draw("same");
   histdEdxVTotalRangeBackwards->GetXaxis()->SetTitle("roughRange [mm]");
   histdEdxVTotalRangeBackwards->GetYaxis()->SetTitle("#frac{dE}{dx} [ADC/mm]");

   /*
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

   TCanvas *c7 = new TCanvas();
   histThetaLABThetaLAB->SetDirectory(0);
   histThetaLABThetaLAB->Draw("zcol");
   //kine_d3He_tt->Draw("same");
   histThetaLABThetaLAB->GetXaxis()->SetTitle("track1_#theta_{LAB} [deg]");
   histThetaLABThetaLAB->GetYaxis()->SetTitle("track2_#theta_{LAB} [deg]");

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
   histRangeVThetaLAB->Write();
   // dE Vs Total Range
   histdEdxVTotalRange->Write();
   histdEdxVTotalRangeBackwards->Write();
   //   cutPIDproton->Write("PIDCutProton");
   //   cutPIDproton_extension->Write("PIDCutProtonExtension");
   //   cutPIDdeuteron->Write("PIDCutDeuteron");

   // Kinematics
   histEstimatedKinEVThetaLABTotal->Write();
   //   histEstimatedKinEVThetaLAB_ProtonATTPC->Write();
   //   histEstimatedKinEVThetaLAB_ProtonATTPC_extended->Write();
   //   histEstimatedKinEVThetaLAB_DeuteronATTPC->Write();
   //   histEstimatedKinEVThetaLAB_CarbonSi->Write();
   //   histEstimatedKinEVThetaLAB_NitrogenSi->Write();
   //   kine_dp_gs->Write("kin_dp_gs");
   //   kine_dd_gs->Write("kin_dd_gs");
   //   kine_dd_gs_25MeVu->Write("kin_dd_gs_25MeVu");
   //   kine_dp_gs_25MeVu->Write("kin_dp_gs_25MeVu");
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
   histThetaLABThetaLAB->Write();
   Results->Close();

   // save canvases
   c1->SaveAs("data/kine_charge_range.C");
   c2->SaveAs("data/kine_thetaLAB_range.C");
   c3->SaveAs("data/kine_dEdx_range.C");
   c4->SaveAs("data/kine_dEdx_range_backwards.C");
   //   c5->SaveAs("data/kine_kinE_thetaLAB.C");
   //   c6->SaveAs("data/kine_kinE_thetaLAB_kinematics.C");
   //   c7->SaveAs("data/kine_thetaLAB_thetaLAB.C");
   //   c8->SaveAs("data/kine_Ex_dp.C");
   //   c9->SaveAs("data/kine_Ex_dd.C");

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
