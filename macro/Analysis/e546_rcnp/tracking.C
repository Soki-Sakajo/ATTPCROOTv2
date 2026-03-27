void tracking(Int_t run_number = 63)
{

   FairRunAna *run = new FairRunAna();

   std::ofstream hitsFile;
   hitsFile.open("hits.txt");

   TString inputpath = "~/fair_install/decode_data/";
   TString inputfile = inputpath + TString::Format("run_%04d.root",run_number);
   TFile *file = new TFile(inputfile, "READ");

   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");
   
   Int_t single_track = 0;
   Int_t multi_track = 0;

   for (Int_t i = 0; i < nEvents; i++) {

      Reader.Next();

      AtEvent *event = (AtEvent *)eventHArray->At(0);
      AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

      if (event && patternEvent) {

         auto &hitArray = event->GetHits();
         auto &tracks = patternEvent->GetTrackCand();
         //std::cout << " Number of hits : " << hitArray.size() << std::endl;
         //std::cout << " Number of tracks : " << tracks.size() << std::endl;
	 
	 if(tracks.size() == 1)single_track += 1;
	 else if(tracks.size() > 1)multi_track += 1;
	
	 if(tracks.size() >= 6)cout<<"find more than 5 tracks ! "<<i<<" track_num = "<<tracks.size()<<endl;
         /*for(auto &hit : hitArray) {
             auto pos = hit->GetPosition();
             auto charge = hit->GetCharge();
             auto time = hit->GetTimeStamp();
             hitsFile  << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
         }*/

         /*for (auto &track : tracks) {
            auto &points = track.GetHitArray();
            std::cout << "    Track ID: " << track.GetTrackID() << std::endl;
            for (auto &point : points) {
               auto pos = point->GetPosition();
               auto charge = point->GetCharge();
               auto time = point->GetTimeStamp();
               std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
            }
	    std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
         }*/
      }

      if(i % 10000 == 0)cout<<"processing : "<<i<<" / "<<nEvents<<endl;
   }
   
   cout<<" single track event number : "<<single_track<<endl;
   cout<<" multi  track event number : "<<multi_track<<endl;
   cout<<"multi/single ratio : "<< double(multi_track *1.0 / single_track *1.0)<<endl;

   hitsFile.close();
   file->Close();
}
