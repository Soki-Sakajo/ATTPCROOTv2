{
  //input file
  //  TFile *file = new TFile("data/attpcsim.root");
  
  //input tree of the event data
  //  TTree *tree = (TTree*)file->Get("cbmsim");

  TFile *_file0 = TFile::Open("data/attpcsim.root");

  cbmsim->SetAlias("ckineE","(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)*1000000/2/11177.736");
  cbmsim->SetAlias("hkineE","(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)*1000000/2/938.7669711");
  cbmsim->SetAlias("ang","180*acos(MCTrack.fPz/(sqrt(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)))/3.141592");
  //  cbmsim->SetAlias("ang","180*acos(MCTrack.fPz/(sqrt(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)))/3.141592");  

  TCut cut_12c = "MCTrack.fPdgCode>=1.000e+09";
  TCut scatt   = "Entry$%2==1";  
  TCut cut_1h  = "MCTrack.fPdgCode==2212";   // should be checked

  //  TH2F *h1 = new TH2F("h1","kineE:ang",180,0,180,50,0,50,"colz");

  //  cbmsim->Draw("ckineE:ang>>h1(100,0,100,50,0,50)", scatt && cut_12c ,"colz");
  cbmsim->Draw("hkineE:ang>>h1(100,0,100,50,0,15)", scatt ,"colz"); 
  //&& "MCTrack.fStartZ>=25" && "MCTrack.fStartZ<35"
  
  h1->SetTitle("kineE:ang");
  //{scatt&&cut_1h && 25<=StartZ<35}
  
  //  h1->GetXaxis()->SetTitle("ang");
  //  h1->GetYaxis()->SetTitle("kineE");
    

  
  //  tree->Print();

  //  TString px=Form("MCTrack.fPx");

    //  TTree *ener=new TTree("(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)*1000000/2/938.2");
  
  
  //  tree->Draw("(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)*1000000/2/11177.736:MCTrack.fStartZ>>h1","Entry$%2==1&&MCTrack.fPdgCode<10000", "colz");
  //  tree->Draw("(MCTrack.fPx*MCTrack.fPx + MCTrack.fPy*MCTrack.fPy + MCTrack.fPz*MCTrack.fPz)*1000000/2/11177.736:180-acos(MCTrack.fPz*11177.736/sqrt(MCTrack.fPx*MCTrack.fPx+MCTrack.fPy*MCTrack.fPy+MCTrack.fPz*MCTrack.fPz))*180/3.141592653589793>>h1","MCTrack.fPdgCode>10000", "colz");  
    //  tree->Draw("ener>>h1","Entry$%2==1&&MCTrack.fPdgCode<10000", "colz");

  //  tree->Scan("px");
  //  cout<< px <<endl;


  
  


}
