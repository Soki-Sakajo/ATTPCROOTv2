void extract_eventDisplay(Int_t subnum=0, Int_t nEvent=0){
  
   //-----User Settings:-----------------------------------------------
  //   TString InputDataFile = "./data/run_sim_0_0/attpcsim.root";
  //   TString ParFile = "./data/run_sim_0_0/attpcpar.root";
  //   TString OutputDataFile = "./data/run_sim_0_0/attpctest.root";

  TString InputDataFile = "./data1/extract_c12_" + to_string(subnum) + "_" + to_string(nEvent) + ".root";

  TString ParFile = "./data1/attpcpar_12c12c_" + to_string(subnum) + ".root";

  TString OutputDataFile = "./data1/extract_attpctest.root";
  
  // -----   Reconstruction run   -------------------------------------------
  
  FairRunAna *fRun = new FairRunAna();
  FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
  FairFileSource *source = new FairFileSource(InputDataFile);
  fRun->SetSource(source);
  fRun->SetSink(sink);
  // fRun->SetGeomFile(GeoDataPath);
  
  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *parInput1 = new FairParRootFileIo();
  parInput1->open(ParFile.Data());
  rtdb->setFirstInput(parInput1);
  
  FairEventManager *fMan = new FairEventManager();
  
  //----------------------Traks and points -------------------------------------
  // FairMCTracks *Track = new FairMCTracks("Monte-Carlo Tracks");
  FairMCPointDraw *AtTpcPoints = new FairMCPointDraw("AtTpcPoint", kBlue, kFullSquare);
  
  // fMan->AddTask(Track);
  fMan->AddTask(AtTpcPoints);
  
  fMan->Init();
}
