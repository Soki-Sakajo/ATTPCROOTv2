void anadep(){
  TFile *infile = new TFile("dep_ene.root");

  TH2 *hist = (TH2*)infile->Get("h_dep2");

  double dx=1.0;
  double x=0;
  int i=0;

  TF1 *fitfunc = new TF1("fitfunc", "gaus");

  vector<float> vec_x;
  vector<float> vec_y;  
  
  
  while(x<150){
    Int_t bin1 = hist->GetXaxis()->FindBin(x);
    Int_t bin2 = hist->GetXaxis()->FindBin(x+dx);
    //    if(x>100) break;

    hist->ProjectionY("proy",bin1, bin2);
    TH1F *h1 = (TH1F*)gROOT->FindObject("proy");

    double y_mean = h1->GetMean();
    vec_x.push_back(x);
    vec_y.push_back(y_mean);    
    
    x+= dx;
    i++;
  }

  TGraph *gr = new TGraph(vec_x.size(), &vec_y[0], &vec_x[0]);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);
  gr->Draw("AP");

  //  vector<float> vec_e_cen = {19.0, 19.2, 19.4, 19.5, 19.6, 19.8, 20.0, 20.2,
  //    20.5, 20.8, 21.0, 21.2, 21.5, 21.8, 22.0, 22.2, 22.5, 22.8, 23.0, 23.2, 23.5,
  //    24.0, 24.5, 25.0};

  vector<float> vec_e_cen = {17.5, 18.0, 18.5, 19.0, 19.2, 19.4, 19.5, 19.6, 19.8, 20.0, 20.2,
    20.5, 20.8, 21.0, 21.2, 21.5, 21.8, 22.0, 22.2, 22.5, 22.8, 23.0, 23.2, 23.5, 24.0};

  vector<float> vec_dep;
  vector<float> vec_depup;
  vector<float> vec_depdown;
  
  printf("energy dz\n");
  for(int i=0; i< vec_e_cen.size(); i++){
    float e_first = 2.0*vec_e_cen[i] -0.1;
    float e_last  = 2.0*vec_e_cen[i] +0.1;    

    float d_first = gr->Eval(e_first);
    float d_last  = gr->Eval(e_last);
    float dep = fabs(d_last - d_first);

    vec_dep.push_back(dep);
  }

  TGraph *gr2 = new TGraph(vec_e_cen.size(), &vec_e_cen[0], &vec_dep[0]);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(4);
  TCanvas *can = new TCanvas();

  TF1 *fitfunc2 = new TF1("fitfunc2", "pol3");
  gr2->Fit("fitfunc2");

  gr2->Draw("AP");
  fitfunc2->Draw("same");

  cout<<endl;
  for(int i=0; i< vec_e_cen.size(); i++){
    printf("%f %f\n", vec_e_cen[i], fitfunc2->Eval(vec_e_cen[i]));
  }
  
}
