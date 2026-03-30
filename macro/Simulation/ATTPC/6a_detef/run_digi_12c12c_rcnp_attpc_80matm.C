bool reduceFunc(AtRawEvent *evt);

void run_digi_12c12c_rcnp_attpc_80matm(int subnum=0)
{
  //  TString inOutDir = "./data/";
  TString inOutDir = "./data1/";

  TString outputFile=
    //    inOutDir + "output_digi_rcnp_hoylehoyle_" + std::to_string(subnum) + ".root";
    inOutDir + "output_digi_rcnp_sequential_hoylehoyle_80matm_" + std::to_string(subnum) + ".root";
  
  TString scriptfile = "Lookup20150611.xml";
  
  //  TString paramFile = "ATTPC.e20009_sim.par";
  //  TString paramFile = "rcnp_attpc.par";
  TString paramFile = "rcnp_attpc_80matm.par";
  
  TString dir = getenv("VMCWORKDIR");
  
  //  TString mcFile = inOutDir + "attpcsim_hoylehoyle_"+ to_string(subnum) +".root";
  TString mcFile = inOutDir + "attpcsim_sequential_hoylehoyle_80matm_"+ to_string(subnum) +".root";
  
  // Create the full parameter file paths
  TString digiParFile = dir + "/parameters/" + paramFile;
  TString mapParFile = dir + "/scripts/" + scriptfile;
  
  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  
  // ------------------------------------------------------------------------
  
  // __ Run ____________________________________________
  FairRunAna *fRun = new FairRunAna();
  FairFileSource *source = new FairFileSource(mcFile);
  fRun->SetSource(source);
  fRun->SetOutputFile(outputFile);
  
  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
  parIo1->open(digiParFile.Data(), "in");
  rtdb->setFirstInput(parIo1);
  
  // Create the detector map to pass to the simulation
  auto mapping = std::make_shared<AtTpcMap>();
  mapping->ParseXMLMap(mapParFile.Data());
  mapping->GeneratePadPlane();
  mapping->ParseInhibitMap("./data/inhibit.txt", AtMap::InhibitType::kTotal);
  
  AtClusterizeTask *clusterizer = new AtClusterizeTask(std::make_shared<AtClusterizeLine>());
  clusterizer->SetPersistence(kFALSE);
  
  AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulseLine>(mapping));
  pulse->SetPersistence(kTRUE);
  
  auto psa = std::make_unique<AtPSAMax>();
  psa->SetThreshold(0);

  // Create PSA task
  AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
  psaTask->SetPersistence(kTRUE);
  /*
  AtPRAtask *praTask = new AtPRAtask();
  praTask->SetPersistence(kTRUE);
  */

  AtRansacTask *ransacTask = new AtRansacTask();
  ransacTask->SetPersistence(kTRUE);
  ransacTask->SetVerbose(kFALSE);
  ransacTask->SetDistanceThreshold(15.0); //12
  ransacTask->SetMinHitsLine(30); //10
  // in AtRansacTask pattern tyepe set to line: auto patternType = AtPatternType::kLine;
  //1=Homemade Ransac(default); 2=Homemade Mlesac; 3=Homemade Lmeds; //4
  ransacTask->SetAlgorithm(1);
  // SampleMethod{kUniform=0,kChargeWeighted=1,kGaussian=2, kWeightGaussian=3, kWeightedY=4};//2
  ransacTask->SetRanSamMode(5);
  ransacTask->SetChargeThreshold(0); //150
  ransacTask->SetNumItera(500);

  fRun->AddTask(clusterizer);
  fRun->AddTask(pulse);
  fRun->AddTask(psaTask);
  //  fRun->AddTask(praTask);
  fRun->AddTask(ransacTask);
  
  //  __ Init and run ___________________________________
  fRun->Init();
  
  timer.Start();
  fRun->Run(0, 10000000);
  timer.Stop();
  
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl << std::endl;
  // -----   Finish   -------------------------------------------------------
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------
}

bool reduceFunc(AtRawEvent *evt)
{
  return (evt->GetNumPads() > 0) && evt->IsGood();
}
