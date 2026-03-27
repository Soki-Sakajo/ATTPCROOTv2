void pad_plot(Int_t eve_num=81){

  double fThreshold =0; //Only save hits if they're above this threshold
  TChain tpc_tree("cbmsim");

  tpc_tree.Add("./data1/output_digi_rcnp_sequential_hoylehoyle_0.root");

  TString outcan="test.pdf";
  //  TString outcan="can_padplane_12c_130matm_0-"+ to_string(eve_num) +".pdf";
  //  TString outcan="can_padplane/can_padplane_sequential_hoylehoyle_130matm_0-"+ to_string(eve_num) +".pdf";

  TTreeReader reader(&tpc_tree);
  TTreeReaderValue<TClonesArray> event(reader,"AtEventH");

  //Create pad plane and load map
  TString mapFile = "Lookup20150611.xml";
  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString mapDir = dir + "/scripts/"+mapFile;
  auto fAtMapPtr = new AtTpcMap();
  fAtMapPtr->ParseXMLMap(mapDir.Data());
  fAtMapPtr->GeneratePadPlane();
  auto fPadPlane = fAtMapPtr->GetPadPlane();
  
  // generate the pad plane histogram
  gSystem->Load("/home/attpc/fair_install/ATTPCROOTv2/build/lib/libAtMap.so");
  AtMap *maptest = new AtTpcMap();
  maptest->GeneratePadPlane();
  TH2Poly *padplane = maptest->GetPadPlane();

  Double_t x,y,z,q;
  
  //get hit pattern data from digi
  std::vector<std::vector<Short_t>> traces; // traces[event][time]

  //get the event to read
  reader.SetEntry(eve_num);
  //Get the event
  AtEvent *eventPtr =(AtEvent *)(event->At(0));
  auto numHits = eventPtr->GetNumHits();
  vector<vector<Double_t>> hitpoint(numHits,vector<Double_t>(5,0));

  for(Int_t i=0; i<numHits; i++){
    auto hit=eventPtr->GetHit(i);
    auto position=hit.GetPosition();
    auto charge=hit.GetCharge();
    hitpoint.at(i).at(0)=eve_num;
    hitpoint.at(i).at(1)=position.X();
    hitpoint.at(i).at(2)=position.Y();
    hitpoint.at(i).at(3)=position.Z();
    hitpoint.at(i).at(4)=charge;
  }

  // fill in the histogram
  for(Int_t i = 0; i <numHits; i++){    
    x=hitpoint.at(i).at(1);
    y=hitpoint.at(i).at(2);
    q=hitpoint.at(i).at(4);
    padplane->Fill(x,y,q);
  }
  /*
  TCanvas *can=new TCanvas("hists","hists",1000,900);
  can->cd(1);
  can->SetLeftMargin(0.12);
  can->SetRightMargin(0.2);
  padplane->SetStats(0);
  padplane->Draw("colz");

  padplane->GetXaxis()->SetLabelSize(0.06);
  padplane->GetYaxis()->SetLabelSize(0.06);
  padplane->GetZaxis()->SetLabelSize(0.06);
  gStyle->SetOptTitle(0);
  gPad->Update();
  */
  TCanvas *can=new TCanvas("hists","hists",1250,1000);
  can->cd(1);
  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.23);
  can->SetTopMargin(0.08);
  can->SetBottomMargin(0.13);
  padplane->SetStats(0);
  padplane->GetXaxis()->SetTitle("x [mm]");
  padplane->GetYaxis()->SetTitle("y [mm]");
  padplane->GetZaxis()->SetTitle("charge");
  padplane->Draw("colz");
  padplane->GetXaxis()->SetLabelSize(0.06);
  padplane->GetYaxis()->SetLabelSize(0.06);
  padplane->GetZaxis()->SetLabelSize(0.06);
  padplane->GetXaxis()->SetTitleSize(0.06);
  padplane->GetYaxis()->SetTitleSize(0.06);
  padplane->GetYaxis()->SetTitleOffset(1.2);
  padplane->GetZaxis()->SetTitleSize(0.06);
  padplane->GetZaxis()->SetTitleOffset(1.6);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can->SaveAs(outcan);

}
