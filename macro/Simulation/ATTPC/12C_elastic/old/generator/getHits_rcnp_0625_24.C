// Script to pull the hits from a single event and output them in a CSV file

void getHits_rcnp(int eventNumber = 1)
{
  
  int tpcRun = 200;
  double fThreshold = 0; // Only save hits if they're above this threshold
  
  TChain tpc_tree("cbmsim");
  //   tpc_tree.Add(TString::Format("/mnt/analysis/e12014/TPC/unpackedReducedFiltered/run_%04d.root", tpcRun));
  //   tpc_tree.Add("~/fair_install/ATTPCROOTv2/macro/Simulation/ATTPC/12C/data/output_digi_rcnp_12c_0.root");
  //   tpc_tree.Add("~/fair_install/ATTPCROOTv2/macro/Simulation/ATTPC/12C/data1/output_digi_rcnp_p_0.root");
  tpc_tree.Add("./data1/output_digi_rcnp_p_999.root");
  
  TTreeReader reader(&tpc_tree);
  TTreeReaderValue<TClonesArray> event(reader, "AtEventH");
  
  // Open file and skip header
  std::ofstream traceFile("./out_hit_rcnp.dat");
  if (!traceFile.good())
    std::cout << "Failed to open file" << std::endl;
  
  // Create pad plane and load map
  //   TString mapFile = "e12014_pad_mapping.xml"; //"Lookup20150611.xml";
  TString mapFile = "Lookup20150611.xml";
  // Set directories
  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString mapDir = dir + "/scripts/" + mapFile;
  
  auto fAtMapPtr = new AtTpcMap();
  fAtMapPtr->ParseXMLMap(mapDir.Data());
  //   fAtMapPtr->GenerateAtTpc();  // old
  //   auto fPadPlane = fAtMapPtr->GetAtTpcPlane(); // old
  fAtMapPtr->GeneratePadPlane();    // Adam suggest
  auto fPadPlane = fAtMapPtr->GetPadPlane(); // Adam suggest   
  
  std::vector<std::vector<Short_t>> traces; // traces[event][time]
  
  // Get the event to read
  reader.SetEntry(eventNumber);
  // Get the event
  AtEvent *eventPtr = (AtEvent *)(event->At(0));
  auto numHits = eventPtr->GetNumHits();
  
  // Write header
  traceFile << "#x,y,z,charge" << std::endl;
  
  for (int i = 0; i < numHits; ++i) {
    // Get the padNum for the hit
    auto hit = eventPtr->GetHit(i);
    //      auto position = hit->GetPosition();
    //      auto charge = hit->GetCharge();
    auto position = hit.GetPosition();      
    auto charge = hit.GetCharge();
    
    traceFile << position.X() << " " << position.Y() << " " << position.Z() << " " << charge << std::endl;
  }
  traceFile.close();
  std::cout << "Done writing file" << std::endl;
}
