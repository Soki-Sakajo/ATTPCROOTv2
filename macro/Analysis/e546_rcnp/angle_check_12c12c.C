void angle_check_12c12c()
{
	bool sigle_event = true;
	FairRunAna *run = new FairRunAna();

	Int_t run_number = 52;

	TString inputpath = "~/fair_install/decode_data/";
	TString inputfile = inputpath + TString::Format("run_%04d.root",run_number);
	TFile *file = new TFile(inputfile, "READ");

	TTree *tree = (TTree *)file->Get("cbmsim");
	Int_t nEvents = tree->GetEntries();
	std::cout << " Number of events : " << nEvents << std::endl;

	TTreeReader Reader("cbmsim", file);
	TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
	TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");

	TH2D *h1 = new TH2D("h1", "h1", 100, 0, 500, 200, 0, 1000);
	TH2D *htheta = new TH2D("htheta", "htheta", 200, 0, 180, 200, 0, 180);
	TH2D *hphi = new TH2D("hphi", "hphi", 200, -180, 180, 200, -180, 180);
	TH2D *hrange = new TH2D("hrange", "hrange", 200, 0, 1000, 200, 0, 1000);
  	TH2D *h_range_charge = new TH2D("h_range_charge", "hrange", 100, 0, 500, 200, 0, 1000);
	TH2D *h_xy = new TH2D("h_xy", "h_xy", 100, -300, 300, 100, -300, 300);
	TH2D *h_xz = new TH2D("h_xz", "h_xz", 100, -50, 1000, 100, -300, 300);
	TH3D *h_xyz = new TH3D("h_xyz", "h_xyz", 100, 0, 1000, 100, -300, 300, 100,-300, 300);
	gROOT->Macro("theta_theta.C");
	TCutG *theta_theta = (TCutG *)gROOT->FindObject("theta_theta");

	for (Int_t i = 0; i < nEvents; i++) {
		h_xy->Reset();
		h_xz->Reset();
		h_xyz->Reset();
		Reader.Next();
		AtEvent *event = (AtEvent *)eventHArray->At(0);
		AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
		if (event && patternEvent) {
			auto &hitArray = event->GetHits();
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
				// h1->Fill(theta1, theta2);
				if(abs(theta1 + theta2 - 90) < 10)hphi->Fill(phi1, phi2);
				if(abs(abs(phi1 - phi2) - 180.0) < 10.0)htheta->Fill(theta1, theta2); 
				hrange->Fill(track1.GetLinearRange(), track2.GetLinearRange());
				//for (auto &track : tracks) {
				//	Double_t theta1 = track.GetGeoThera();
				//std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
	//}
	//}
				for(auto &hit : hitArray) {
					auto pos = hit->GetPosition();
					auto charge = hit->GetCharge();
					auto time = hit->GetTimeStamp();
				//hitsFile  << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
				}

				for (auto &track : tracks) {
					auto &points = track.GetHitArray();
					//	  std::cout << "    Track ID: " << track.GetTrackID() << std::endl;
					for (auto &point : points) {
						auto pos = point->GetPosition();
						auto charge = point->GetCharge();
						auto time = point->GetTimeStamp();
						//	    std::cout <<"posx = "<< pos.X() << ", posy =  " << pos.Y() << ",posz = " << pos.Z() << ", time =  " << time << ", charge =  " << charge << std::endl;
						h_xy->Fill(pos.X(),pos.Y());
						h_xz->Fill(pos.Z(),pos.Y());
						h_xyz->Fill(pos.Z(),pos.X(),pos.Y());
						h_range_charge->Fill(range1,charge);
					}
					//std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
				}

				if(sigle_event ==1 ){
					TCanvas *c3 = new TCanvas("c1","c1",1200,600);
					c3->Divide(2,1);
					c3->cd(1); h_xy->Draw("colz");
					c3->cd(2); h_xz->Draw("colz");
					c3->Draw();
					c3->Update();
					TCanvas *c5 = new TCanvas();
					h_xyz->Draw("BOX2");
					c5->Draw();
					c5->Update();

					cout<<"Press ENTER to continue to next event" <<endl;
					cin.ignore();
					}
			}
			if(tracks.size() == 3){
				AtTrack track1 = tracks[0];
				AtTrack track2 = tracks[1];
				AtTrack track3 = tracks[2];
				Double_t theta1 = 180.0 - track1.GetGeoTheta()*180.0/TMath::Pi();
				Double_t theta2 = 180.0 - track2.GetGeoTheta()*180.0/TMath::Pi();
				Double_t theta3 = 180.0 - track3.GetGeoTheta()*180.0/TMath::Pi();
				Double_t phi1 = track1.GetGeoPhi()*180.0/TMath::Pi();
				Double_t phi2 = track2.GetGeoPhi()*180.0/TMath::Pi();
				Double_t phi3 = track3.GetGeoPhi()*180.0/TMath::Pi();
				Double_t range1 = track1.GetLinearRange();
				Double_t range2 = track2.GetLinearRange();
				Double_t range3 = track3.GetLinearRange();
				//	     h1->Fill(theta1, theta2);
				//	if(abs(theta1 + theta2 - 90) < 10)hphi->Fill(phi1, phi2);
				//	if(abs(abs(phi1 - phi2) - 180.0) < 10.0)htheta->Fill(theta1, theta2); 
				//	hrange->Fill(track1.GetLinearRange(), track2.GetLinearRange());
				//for (auto &track : tracks) {
				//	Double_t theta1 = track.GetGeoThera();
				//std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
				//}
				//}
				for(auto &hit : hitArray) {
					auto pos = hit->GetPosition();
					auto charge = hit->GetCharge();
					auto time = hit->GetTimeStamp();
					//hitsFile  << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
				}
				for (auto &track : tracks) {
					auto &points = track.GetHitArray();
					std::cout << "    Track ID: " << track.GetTrackID() << std::endl;
					for (auto &point : points) {
						auto pos = point->GetPosition();
						auto charge = point->GetCharge();
						auto time = point->GetTimeStamp();
						std::cout <<"posx = "<< pos.X() << ", posy =  " << pos.Y() << ",posz = " << pos.Z() << ", time =  " << time << ", charge =  " << charge << std::endl;
						h_xy->Fill(pos.X(),pos.Y());
						h_xz->Fill(pos.Z(),pos.Y());
						h_xyz->Fill(pos.Z(),pos.X(),pos.Y());
						h_range_charge->Fill(range1,charge);
					}
					//std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
				}

				if(sigle_event ==1 ){
					TCanvas *c3 = new TCanvas("c1","c1",1200,600);
					c3->Divide(2,1);
					c3->cd(1); h_xy->Draw("colz");
					c3->cd(2); h_xz->Draw("colz");
					c3->Draw();
					c3->Update();
					TCanvas *c5 = new TCanvas();
					h_xyz->Draw("BOX2");
					c5->Draw();
					c5->Update();
					cout<<"Press ENTER to continue to next event" <<endl;
					cin.ignore();
				}
			}
		}
		if(i % 10000 == 0)cout<<"processing : "<<i<<" / "<<nEvents<<endl;
	}

	//file->Close();
	TCanvas *c1 = new TCanvas();
	htheta->Draw("colz");
	c1->Draw();

	TCanvas *c2 = new TCanvas();
	hphi->Draw("colz");
	c2->Draw();

	TCanvas *c4 = new TCanvas();
	h_range_charge->Draw("colz");
	c4->Draw();
}
