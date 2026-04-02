void trans_all1(TTree* tree)
{

   FairRunAna *run = new FairRunAna();

   Int_t run_number = 40;

   //TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   //   AtFindVertex *fver=new AtFindVertex();
  
   //TTreeReader Reader("cbmsim", file);
   TTreeReader Reader(tree);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");
   
   TString outputfile ="~/fair_install/transfer_data/" + TString::Format("trans_100torr.root");
   TFile *opf = new TFile(outputfile, "RECREATE");
   TTree *new_tree = new TTree("tree", "tree"); 

   opf->cd();
   Int_t run_num = 0;
   Int_t ntrack = 50;
   Double_t track_theta[ntrack];
   Double_t track_phi[ntrack];
   Double_t track_range[ntrack];
   Double_t vertex_range[ntrack];
   Double_t track_charge[ntrack];
   //   Double_t track_vertex[ntrack];
   Double_t track_first_z[ntrack];
   Double_t track_last_z[ntrack];

   Int_t    n_vertex1=0;
   Double_t track_verx[ntrack];
   Double_t track_very[ntrack];
   Double_t track_verz[ntrack];

   new_tree->Branch("run_num", &run_num, "run_num/I");
   new_tree->Branch("ntrack", &ntrack, "ntrack/I");
   new_tree->Branch("theta", track_theta, "theta[ntrack]/D");
   new_tree->Branch("phi", track_phi, "phi[ntrack]/D");
   //   new_tree->Branch("vtxrange", vertex_range, "vtxrange[ntrack]/D");
   new_tree->Branch("range", track_range, "range[ntrack]/D");
   new_tree->Branch("charge", track_charge, "charge[ntrack]/D");
   //   new_tree->Branch("vertex", track_vertex, "vertex[ntrack]/D");
   new_tree->Branch("n_vertex", &n_vertex1, "n_vertex/I");
   new_tree->Branch("vertexx", track_verx, "vertex[10]/D");
   new_tree->Branch("vertexy", track_very, "vertey[10]/D");
   new_tree->Branch("vertexz", track_verz, "vertez[10]/D");
   new_tree->Branch("first_z", track_first_z, "first_z[ntrack]/D");
   new_tree->Branch("last_z", track_last_z, "last_z[ntrack]/D");
   Int_t nevt=0;
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

       if(tracks.size() >1){
	 for (Int_t j=0;j<10;j++){
	   track_theta[j]=-400;
	   track_phi[j]  =-400;
	   track_range[j]=-100;
	   vertex_range[j]=-100;
	   track_charge[j]=-100;
	   track_first_z[j]=-2000;
	   track_last_z[j]=-2000;
	   track_verx[j]=-2000;
	   track_very[j]=-2000;
	   track_verz[j]=-2000;
	 }
	 ntrack = tracks.size();	     

	 AtFindVertex *fver=new AtFindVertex();
	 fver->FindVertexMultipleLines(tracks, ntrack);
	 auto vtxlist = fver->GetTracksVertex();
	 n_vertex1=0;
	 for(auto &v:vtxlist){
	   //  cout<<"vertex number: "<<n_vertex1<<"  size "<<vtxlist.size()<<"  "<<itrack<<endl;
	   //	   cout<<"Found vertex at ("<<v.vertex.X()<<", "<<v.vertex.Y()<<", "<<v.vertex.Z()<<")"<<endl;

	   track_verx[n_vertex1] = v.vertex.X();
	   track_very[n_vertex1] = v.vertex.Y();
	   track_verz[n_vertex1] = v.vertex.Z();
	   n_vertex1++;
	   
	   if(n_vertex1> vtxlist.size()) break;
	   
	 }
	 int itrack = 0;
         for (auto &track : tracks) {
	   
           track_theta[itrack] =180.0 - track.GetGeoTheta()*180.0/TMath::Pi();
           track_phi[itrack] = track.GetGeoPhi()*180.0/TMath::Pi();
           //if (track_phi[itrack]<0) track_phi[itrack]+=180;
           track_range[itrack] = track.GetLinearRange();
           track_charge[itrack] = track.GetGeoQEnergy();
           //      track_vertex[itrack] = track.GetVertexToZDist();
           track_first_z[itrack] = track.GetFirstPoint().Z();
           track_last_z[itrack] = track.GetLastPoint().Z();

	   /*/cout<<" Track range  "<<track_range[itrack]<<endl;
	   if (n_vertex1>0) {
             XYZPoint vtx(vtxlist[0].vertex.X(),vtxlist[0].vertex.Y(),vtxlist[0].vertex.Z());
	     XYZPoint lp = track.GetLastPoint();
	     vertex_range[itrack]= track.GetLinearRange(vtx,lp);
	     lp = track.GetFirstPoint();
	     //	     std::cout<< "    Track range :"<<track_range[itrack]<<"  vtxrange  "<< vertex_range[itrack]<<"  "<<track.GetLinearRange(vtx,lp) << std::endl;
	     } */
	   itrack++; 
         }
	 // cout<<i<<" nvertex  "<<n_vertex1<<endl;
	 new_tree->Fill();
	 nevt++;
	 delete fver;
       }
     }

     if(i % 1000 == 0)cout<<"processing : "<<i<<" / "<<nEvents<<"  "<<nevt<<" \r"<<flush;
   }
   cout<<endl;
   new_tree->Write();
   opf->Close();
   
   //   file->Close(); 
}

void trans_all(){
  TChain *ch=new TChain("cbmsim");
  TString inputpath = "~/fair_install/decode_data/";
  Int_t nfile=0;
  for (int i=115;i<138;i++){
    TString inputfile = inputpath + TString::Format("run_%04d.root",i);
    if (ch->Add(inputfile.Data())){
      cout<<" "<<nfile<<"  Added file: "<<inputfile<<endl;
      nfile++;
    }
  }
  long long int nevt=ch->GetEntries();
  cout<<endl<<" There are "<<nfile<<" runs were added in with "<<nevt<<" events. "<<endl;
  trans_all1(ch);
}
