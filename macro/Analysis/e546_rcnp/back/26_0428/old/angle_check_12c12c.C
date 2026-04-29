void angle_check_12c12c(){

	Double_t del_phi = 10; // cut value; phi1 - phi2 - 180 deg < del_phi
	Double_t del_theta = 10; // cut value; theta1 + theta2 - 90 deg < del_theta
	std::vector run_Nums = {52};
	//	std::vector run_Nums = {52,53,54,56,57,58,59};
	std::ofstream histfile("can_output/angle_check_12c12c.C", std::ios::out | std::ios::trunc);

	FairRunAna *run = new FairRunAna();

	// histogram definitions.
	// theta
	TH2D *h_theta = new TH2D("h_theta", "h_theta", 200, 0, 180, 200, 0, 180);
	//	TH2D *h_theta_cutphi = new TH2D("h_theta_cutphi", "h_theta_cutphi", 200, 0, 180, 200, 0, 180);
	TH2D *h_theta_cutphi_12c = new TH2D("h_theta_cutphi_12c", "h_theta_cutphi_12c", 200, 0, 180, 200, 0, 180);

	// phi
	TH2D *h_phi = new TH2D("h_phi", "h_phi", 200, -180, 180, 200, -180, 180);
	TH2D *h_phi_cuttheta = new TH2D("h_phi_cuttheta", "h_phi_cuttheta", 200, -180, 180, 200, -180, 180);
	//	TH2D *h_phi_cutphi = new TH2D("h_phi_cutphi", "h_phi_cutphi", 200, -180, 180, 200, -180, 180);
	TH2D *h_phi_cutphi_12c = new TH2D("h_phi_cutphi_12c", "h_phi_cutphi_12c", 200, -180, 180, 200, -180, 180);

	// range
	TH2D *h_range = new TH2D("h_range", "h_range", 600, 0, 1200, 600, 0, 1200);
  	TH2D *h_charge_range = new TH2D("h_range_charge", "h_range_charge", 600, 0, 1200, 600, 0, 6e5);
  	TH2D *h_charge_range_cut12c = new TH2D("h_range_charge_cut12c", "h_range_charge_cut12c", 600, 0, 1200, 600, 0, 6e5);
	TH2D *h_charge_range_cutphi_12c = new TH2D("h_range_charge_cutphi_12c", "h_range_charge_cutphi_12c", 600, 0, 1200, 600, 0, 6e5);
	TH2D *h_range_theta = new TH2D("h_range_theta", "h_range_theta", 180, 0, 180, 1030, 0, 1030);
	TH2D *h_range_theta_cut12c = new TH2D("h_range_theta_cut12c", "h_range_theta_cut12c", 180, 0, 180, 1030, 0, 1030);

	// others
	//	TH2D *h_xy = new TH2D("h_xy", "h_xy", 100, -300, 300, 100, -300, 300);
	//	TH2D *h_xz = new TH2D("h_xz", "h_xz", 100, -50, 1000, 100, -300, 300);
	//	TH3D *h_xyz = new TH3D("h_xyz", "h_xyz", 100, 0, 1000, 100, -300, 300, 100,-300, 300);

	// get cut from file.
	TFile *cutfile_12c = new TFile("cut_files/charge_range_12c.root", "READ");
	TCutG *cut_charge_range_12c = (TCutG *)cutfile_12c->Get("charge_range_12c");

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

		for (Int_t i = 0; i < nUnpackEvents; i++) {
			unpackReader.Next();
			AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
			if (event && patternEvent) {
				auto &tracks = patternEvent->GetTrackCand();
				//std::cout << " Number of hits : " << hitArray.size() << std::endl;
				//std::cout << " Number of tracks : " << tracks.size() << std::endl;

				if(tracks.size() == 2){
					AtTrack track1 = tracks[0];
					AtTrack track2 = tracks[1];
					Double_t theta1 = 180.0 - track1.GetGeoTheta()*180.0/TMath::Pi();
					Double_t theta2 = 180.0 - track2.GetGeoTheta()*180.0/TMath::Pi();
					Double_t phi1 = track1.GetGeoPhi()*180.0/TMath::Pi();
					Double_t phi2 = track2.GetGeoPhi()*180.0/TMath::Pi();
					Double_t range1 = track1.GetLinearRange();
					Double_t range2 = track2.GetLinearRange();
					Double_t charge1 = track1.GetGeoQEnergy();
					Double_t charge2 = track2.GetGeoQEnergy();
					bool is12c12c = cut_charge_range_12c->IsInside(range1, charge1) && cut_charge_range_12c->IsInside(range2, charge2);

					h_charge_range->Fill(range1, charge1);
					h_charge_range->Fill(range2, charge2);
					h_range_theta->Fill(theta1, range1);
					h_range_theta->Fill(theta2, range2);
					h_range->Fill(range1, range2);
					h_theta->Fill(theta1, theta2);
					h_phi->Fill(phi1, phi2);

					/*
					if(abs(theta1 + theta2 - 90) < del_theta){
						h_phi_cuttheta->Fill(phi1, phi2);
					}
					*/

					if(abs(abs(phi1 - phi2) - 180.0) < del_phi){
						//						h_theta_cutphi->Fill(theta1, theta2);
						//						h_phi_cutphi->Fill(phi1, phi2);
						if(is12c12c){
							h_charge_range_cutphi_12c->Fill(range1, charge1);
							h_charge_range_cutphi_12c->Fill(range2, charge2);
							h_theta_cutphi_12c->Fill(theta1, theta2);
							h_phi_cutphi_12c->Fill(phi1, phi2);
						}
					}
					if(is12c12c){
						h_charge_range_cut12c->Fill(range1, charge1);
						h_charge_range_cut12c->Fill(range2, charge2);
						h_range_theta_cut12c->Fill(theta1, range1);
						h_range_theta_cut12c->Fill(theta2, range2);
					}
				}
			}
			if(i % 100 == 0) {
				std::cout << "  Processing data: " << 100*i/nUnpackEvents << " %!    \r" << std::flush;
			}
		}
	unpackfile->Close();
	}

	//file->Close();
	TCanvas *c1 = new TCanvas("c1", "c1");
	h_charge_range->SetDirectory(0);
	h_charge_range->Draw("colz");
	h_charge_range->GetXaxis()->SetTitle("roughRange [mm]");
	h_charge_range->GetYaxis()->SetTitle("Charge [ADC]");
	h_charge_range->SetTitle("charge range");

	TCanvas *c2 = new TCanvas("c2", "c2");
	h_range_theta->SetDirectory(0);
	h_range_theta->Draw("colz");
	h_range_theta->GetXaxis()->SetTitle(" #theta [deg]");
	h_range_theta->GetYaxis()->SetTitle("roughRange [mm]");
	h_range_theta->SetTitle("range #theta");

	TCanvas *c3 = new TCanvas("c3", "c3");
	h_range->SetDirectory(0);
	h_range->Draw("colz");
	h_range->GetXaxis()->SetTitle("track1 roughRange [mm]");
	h_range->GetYaxis()->SetTitle("track2 roughRange [mm]");
	h_range->SetTitle("range track1 track2");

	TCanvas *c4 = new TCanvas("c4", "c4");
	h_theta->SetDirectory(0);
	h_theta->Draw("colz");
	h_theta->GetXaxis()->SetTitle("track1 #theta_1 [deg]");
	h_theta->GetYaxis()->SetTitle("track2 #theta_2 [deg]");
	h_theta->SetTitle("#theta_1 #theta_2");

	TCanvas *c5 = new TCanvas("c5", "c5");
	h_phi->SetDirectory(0);
	h_phi->Draw("colz");
	h_phi->GetXaxis()->SetTitle("track1 #phi_1 [deg]");
	h_phi->GetYaxis()->SetTitle("track2 #phi_2 [deg]");
	h_phi->SetTitle("#phi_1 #phi_2");

	TCanvas *c6 = new TCanvas("c6", "c6");
	h_charge_range_cut12c->SetDirectory(0);
	h_charge_range_cut12c->Draw("colz");
	h_charge_range_cut12c->GetXaxis()->SetTitle("roughRange [mm]");
	h_charge_range_cut12c->GetYaxis()->SetTitle("Charge [ADC]");
	h_charge_range_cut12c->SetTitle("charge range (cut 12C-12C)");

	TCanvas *c7 = new TCanvas("c7", "c7");
	h_range_theta_cut12c->SetDirectory(0);
	h_range_theta_cut12c->Draw("colz");
	h_range_theta_cut12c->GetXaxis()->SetTitle(" #theta [deg]");
	h_range_theta_cut12c->GetYaxis()->SetTitle("roughRange [mm]");
	h_range_theta_cut12c->SetTitle("range #theta (cut 12C-12C)");

	TCanvas *c8 = new TCanvas("c8", "c8");
	h_charge_range_cutphi_12c->SetDirectory(0);
	h_charge_range_cutphi_12c->Draw("colz");
	h_charge_range_cutphi_12c->GetXaxis()->SetTitle("roughRange [mm]");
	h_charge_range_cutphi_12c->GetYaxis()->SetTitle("Charge [ADC]");
	h_charge_range_cutphi_12c->SetTitle(Form("charge range (phi1-phi2-180 < %d )", (int)del_phi));

	TCanvas *c9 = new TCanvas("c9", "c9");
	h_theta_cutphi_12c->SetDirectory(0);
	h_theta_cutphi_12c->Draw("colz");
	h_theta_cutphi_12c->GetXaxis()->SetTitle("track1 #theta_1 [deg]");
	h_theta_cutphi_12c->GetYaxis()->SetTitle("track2 #theta_2 [deg]");
	h_theta_cutphi_12c->SetTitle(Form("#theta_1 #theta_2 (phi1-phi2-180 < %d )", (int)del_phi));

	TCanvas *c10 = new TCanvas("c10", "c10");
	h_phi_cutphi_12c->SetDirectory(0);
	h_phi_cutphi_12c->Draw("colz");
	h_phi_cutphi_12c->GetXaxis()->SetTitle("track1 #phi_1 [deg]");
	h_phi_cutphi_12c->GetYaxis()->SetTitle("track2 #phi_2 [deg]");
	h_phi_cutphi_12c->SetTitle(Form("#phi_1 #phi_2 (phi1-phi2-180 < %d )", (int)del_phi));

	// save canvases as .C macros
	TSeqCollection *canlist = gROOT->GetListOfCanvases();
	TIter next(canlist);
	Int_t can_num = 1;
	Int_t n_can = canlist->GetEntries();
	histfile << "void angle_check_12c12c(){\n" << std::endl;
	histfile <<  Form("std::cout << \"load %d canvases... \" << std::endl;\n", n_can) << std::endl;
	TCanvas *c;
	while ((c = (TCanvas *)next())) {
		histfile << Form("TCanvas *c%d = new TCanvas(\"c%d\", \"c%d\");\n", can_num, can_num, can_num);
		histfile << Form("c%d->cd();\n", can_num);
		c->SavePrimitive(histfile, "");
		histfile << "\n";
		histfile << Form("std::cout << \"Loading canvases: %d / %d \\r\" << std::flush;\n", can_num, n_can);
		can_num++;
	}
	histfile <<  Form("std::cout << \"Drawing canvases... \" << std::endl;\n") << std::endl;
	histfile << "}\n" << std::endl;
	histfile.close();
}
