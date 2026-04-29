void check(){

  TFile *f1=new TFile("r_under250.root");
  TH1F* h1 =(TH1F*)f1->Get("h2");
  TFile *f2=new TFile("r_over250.root");
  TH1F* h2 =(TH1F*)f2->Get("h3");

  if(!h1||!h2){
    std::cerr<<"h1 could NOT take the data!!!"<<std::endl;
    return;
  }
  if(!h2){
    std::cerr<<"h1 could NOT take the data!!!"<<std::endl;
    return;
  }


  TH1F* h_add=(TH1F*)h1->Clone("h_add");
  h_add->Add(h2);

  TH1F* h_div=(TH1F*)h2->Clone("h_div");
  h_div->Divide(h_add);

  TCanvas *can=new TCanvas();
  can->Divide(1,2);
  can->cd(1);
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  can->cd(2);
  h_div->SetLineColor(1);
  h_div->Draw();
  
}
