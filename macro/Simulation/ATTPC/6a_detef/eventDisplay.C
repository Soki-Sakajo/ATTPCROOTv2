void eventDisplay(Int_t subnum=10)
{
  //-----User Settings:-----------------------------------------------
  //  TString inOutDir ="./data/";
  TString inOutDir ="./data1/";

  /*
  TString InputDataFile = inOutDir + "attpcsim_50matm_"+ to_string(subnum) +".root";
  TString ParFile = inOutDir + "attpcpar_50matm_"+ to_string(subnum) +".root";
  TString OutputDataFile = inOutDir + "attpctest.root";
  */

  TString InputDataFile = inOutDir + "attpcsim_sequential_hoylehoyle_80matm_"+ to_string(subnum) +".root";
  TString ParFile = inOutDir + "attpcpar_sequential_hoylehoyle_80matm_"+ to_string(subnum) +".root";
  TString OutputDataFile = inOutDir + "attpctest.root";


  //  TString InputDataFile = inOutDir + "attpcsim_sequential_hoylehoyle_100matm_"+ to_string(subnum) +".root";
  //  TString ParFile = inOutDir + "attpcpar_sequential_hoylehoyle_100matm_"+ to_string(subnum) +".root";
  //  TString OutputDataFile = inOutDir + "attpctest.root";

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
