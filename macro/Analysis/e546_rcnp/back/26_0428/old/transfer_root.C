void suppressCout(){
  static std::ofstream nullStream("/dev/null");
  std::cout.rdbuf(nullStream.rdbuf());
}

void restoreCout(){
  static std::streambuf* defaultBuf = std::cout.rdbuf();
  std::cout.rdbuf(defaultBuf);
}

void transfer_root(int run_number=42)
{

   FairRunAna *run = new FairRunAna();

   //Int_t run_number = 40;

   TString inputpath = "~/fair_install/decode_data/";
   TString inputfile = inputpath + TString::Format("run_%04d.root",run_number);
   TFile *file = new TFile(inputfile, "READ");
   
   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   AtFindVertex *fver=new AtFindVertex();
  
   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");
   
   TString outputfile ="~/fair_install/transfer_data/" + TString::Format("run_%04d.root", run_number);
   TFile *opf = new TFile(outputfile, "RECREATE");
   TTree *new_tree = new TTree("tree", "tree"); 

   opf->cd();
   Int_t run_num = 0;
   Int_t ntrack = 0;
   Double_t track_theta[50];
   Double_t track_phi[50];
   Double_t track_range[50];
   Double_t track_charge[50];
   //   Double_t track_vertex[50];
   Double_t track_first_z[50];
   Double_t track_last_x[50];
   Double_t track_last_y[50];
   Double_t track_last_z[50];
   Double_t track_leng_vx[50];
   Double_t track_leng_vy[50];
   Double_t track_leng_vz[50];
   Double_t track_leng[50];
   Double_t track_r[50];

   Double_t track_verx;
   Double_t track_very;
   Double_t track_verz;

   bool check_tracks=false;

   new_tree->Branch("run_num", &run_num, "run_num/I");
   new_tree->Branch("ntrack", &ntrack, "ntrack/I");
   new_tree->Branch("theta", track_theta, "theta[ntrack]/D");
   new_tree->Branch("phi", track_phi, "phi[ntrack]/D");
   new_tree->Branch("range", track_range, "range[ntrack]/D");
   new_tree->Branch("charge", track_charge, "charge[ntrack]/D");
   //   new_tree->Branch("vertex", track_vertex, "vertex[ntrack]/D");
   new_tree->Branch("vertex_x", &track_verx, "track_verx/D");
   new_tree->Branch("vertex_y", &track_very, "track_very/D");
   new_tree->Branch("vertex_z", &track_verz, "track_verz/D");
   new_tree->Branch("first_z", track_first_z, "first_z[ntrack]/D");
   new_tree->Branch("last_x", track_last_x, "last_z[ntrack]/D");
   new_tree->Branch("last_y", track_last_y, "last_z[ntrack]/D");
   new_tree->Branch("last_z", track_last_z, "last_z[ntrack]/D");

   for (Int_t i = 0; i < nEvents; i++) {
     run_num = run_number;
     ntrack = 0;
     Reader.Next();

     AtEvent *event = (AtEvent *)eventHArray->At(0);
     AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

     if (event && patternEvent) {

       auto &hitArray = event->GetHits();
       auto &tracks = patternEvent->GetTrackCand();
       //std::cout << " Number of hits : " << hitArray.size() << std::endl;
       //std::cout << " Number of tracks : " << tracks.size() << std::endl;

       ntrack = tracks.size();
       int itrack = 0;
	   
       for (auto &track : tracks) {
	 track_theta[itrack] =180.0 - track.GetGeoTheta()*180.0/TMath::Pi();
	 track_phi[itrack] = track.GetGeoPhi()*180.0/TMath::Pi();
	 track_range[itrack] = track.GetLinearRange();
	 track_charge[itrack] = track.GetGeoQEnergy();
	 //	   track_vertex[itrack] = track.GetVertexToZDist();
	 track_first_z[itrack] = track.GetFirstPoint().Z();
	 track_last_x[itrack] = track.GetLastPoint().X();
	 track_last_y[itrack] = track.GetLastPoint().Y();
	 track_last_z[itrack] = track.GetLastPoint().Z();
	 //cout<<track.GetVertexToZDist()<<endl;
	 //std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
	 itrack ++;
       }
       check_tracks=false;
       //       suppressCout();
       //       for(int k=1; k<7; k++){
       for(int k=1; k<4; k++){
	 if(tracks.size() != k && tracks.size() != 1){continue;}
	 else if(tracks.size() == 1){
	   fver->FindVertexSingleLine(tracks);
	   check_tracks=true;
	 }
	 else if(tracks.size() == k){
	   fver->FindVertexMultipleLines(tracks, ntrack);
	   check_tracks=true;
	 }
	 auto vtxlist = fver->GetTracksVertex();
	 for(auto &v:vtxlist){
	   //	 cout<<"run number: "<<run_num<<endl;
	   //	 cout<<"Found vertex at ("<<v.vertex.X()<<", "<<v.vertex.Y()<<", "<<v.vertex.Z()<<")"<<endl;
	   track_verx = v.vertex.X();
	   track_very = v.vertex.Y();
	   track_verz = v.vertex.Z();
	 }
       }
       //       restoreCout();
       if(!check_tracks){continue;}
       itrack = 0;
       for (auto &track : tracks) {
	 track_leng_vx[itrack] = track_last_x[itrack]-track_verx;
	 track_leng_vy[itrack] = track_last_y[itrack]-track_very;
	 track_leng_vz[itrack] = track_last_z[itrack]-track_verz;
	 track_leng[itrack] =
	   sqrt(track_leng_vx[itrack]*track_leng_vx[itrack]+
		track_leng_vy[itrack]*track_leng_vy[itrack]+
		track_leng_vz[itrack]*track_leng_vz[itrack]);
	 track_r[itrack] =
	   sqrt(track_leng_vx[itrack]*track_leng_vx[itrack]+
		track_leng_vy[itrack]*track_leng_vy[itrack]);
	 //cout<<track.GetVertexToZDist()<<endl;
	 //std::cout<< "    Track theta :"<< track.GetGeoTheta() << std::endl;
       }
       new_tree->Fill();
     }
     if(i % 1000 == 0)cout<<"processing : "<<i<<" / "<<nEvents<<" \r"<<flush;
   }

   cout<<endl;
   new_tree->Write();
   opf->Close();
   file->Close(); 
}
