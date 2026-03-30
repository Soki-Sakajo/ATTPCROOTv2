void eventDisplay(Int_t subnum=0)
{
  //-----User Settings:-----------------------------------------------
  //  TString inOutDir ="./data/";
  TString inOutDir ="./data1/";

  //  TString InputDataFile = "./data/attpcsim_C12C12.root";
  //  TString ParFile = "./data/attpcpar_C12C12.root";
  //  TString OutputDataFile = "./data/attpctest_C12C12.root";
  //  TString InputDataFile = inOutDir + "attpcsim_12c12c_"+ to_string(subnum) +".root";
  //  TString ParFile = inOutDir + "attpcpar_12c12c_"+ to_string(subnum) +".root";
  //  TString OutputDataFile = inOutDir + "attpctest_12c12c.root";

  /*
  TString InputDataFile = inOutDir + "attpcsim_hoylehoyle_"+ to_string(subnum) +".root";
  TString ParFile = inOutDir + "attpcpar_hoylehoyle_"+ to_string(subnum) +".root";
  TString OutputDataFile = inOutDir + "attpctest_hoylehoyle.root";
  */

  TString InputDataFile = inOutDir + "attpcsim_sequential_hoylehoyle_"+ to_string(subnum) +".root";
  TString ParFile = inOutDir + "attpcpar_sequential_hoylehoyle_"+ to_string(subnum) +".root";
  TString OutputDataFile = inOutDir + "attpctest_sequential_hoylehoyle_"+ to_string(subnum) +".root";


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
