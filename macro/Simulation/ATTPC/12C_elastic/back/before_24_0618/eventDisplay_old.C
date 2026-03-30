void eventDisplay()
{
   //-----User Settings:-----------------------------------------------
  //   TString InputDataFile = "./data/run_sim_0_0/attpcsim.root";
  //   TString ParFile = "./data/run_sim_0_0/attpcpar.root";
  //   TString OutputDataFile = "./data/run_sim_0_0/attpctest.root";

  //    TString InputDataFile = "/home/attpc/opt/ATTPCROOTv2/macro/tests/AT-TPC/data/symFissionMC.root";
  TString InputDataFile = "/home/attpc/opt/ATTPCROOTv2/macro/Simulation/ATTPC/12C_/data/attpcsim.root";   
  //TString ParFile = "/home/attpc/opt/ATTPCROOTv2/macro/tests/AT-TPC/data/par_attpc.root";
  TString ParFile = "/home/attpc/opt/ATTPCROOTv2/macro/Simulation/ATTPC/12C_/data/attpcpar.root";   
  TString OutputDataFile = "data/attpctest.root";
  
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
