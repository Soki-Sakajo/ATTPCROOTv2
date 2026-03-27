int sequential_hoyle_0811(){

  //set situation
  const Double_t AMU = 0.93149432;      // atomic mass unit [GeV]
  const Double_t alpha_mass = 4*AMU+2.4249/1000;     // alpha (= 4He) mass [GeV]
  const Double_t c12_mass = 12*AMU;      // 12C mass [GeV]
  const Double_t be8_mass = 8*AMU+4.9416/1000; // 8Be mass [GeV]
  const Double_t Ex_hoyle = 0.00765407; // 12C(0+2) excitation energy [GeV]
  const Double_t GtoM=1e3;
  
  Int_t n_decay=0;
  Double_t be8alpha_mass[2] ={be8_mass, alpha_mass};
  Double_t alphaalpha_mass[2] ={alpha_mass, alpha_mass};
  Double_t px_alpha[3];
  Double_t py_alpha[3];
  Double_t pz_alpha[3];
  Double_t E_alpha[3];
  Double_t E_hoyle=c12_mass+Ex_hoyle;
  Double_t weight;
  TLorentzVector c12_hoyle(0.0,0.0,0.0,E_hoyle);
  TLorentzVector be8_tmp;
  TLorentzVector *p_alpha[3];
  TGenPhaseSpace event_8BeAlpha,event_2alpha;

  //event generation 12C(hoyle) -> 3 alpha
  // 12C(hoyle) -> 8Be + alpha
  if(event_8BeAlpha.SetDecay(c12_hoyle, 2, be8alpha_mass)){
    weight = event_8BeAlpha.Generate();
    be8_tmp=*event_8BeAlpha.GetDecay(0);
    p_alpha[0]=event_8BeAlpha.GetDecay(1);
    n_decay+=1;
    // 8Be -> 2 alpha    
    if(event_2alpha.SetDecay(be8_tmp,2,alphaalpha_mass)){
      weight = event_2alpha.Generate();
      p_alpha[1]=event_2alpha.GetDecay(0);
      p_alpha[2]=event_2alpha.GetDecay(1);
      n_decay+=2;
    }
  }
  for(int j=0;j<3;j++){
    px_alpha[j]=p_alpha[j]->Px()*GtoM;
    py_alpha[j]=p_alpha[j]->Py()*GtoM;
    pz_alpha[j]=p_alpha[j]->Pz()*GtoM;
    E_alpha[j] =p_alpha[j]->E()*GtoM;
  }
  
  cout<<"alpha 1 (px,py,pz,Energy)= ("
      <<px_alpha[0]<<" ,"<<py_alpha[0]<<" ,"<<py_alpha[0]<<" ,"<<E_alpha[0]<<")"<<endl;
  cout<<"alpha 2 (px,py,pz,Energy)= ("
      <<px_alpha[1]<<" ,"<<py_alpha[1]<<" ,"<<py_alpha[1]<<" ,"<<E_alpha[1]<<")"<<endl;
  cout<<"alpha 3 (px,py,pz,Energy)= ("
      <<px_alpha[2]<<" ,"<<py_alpha[2]<<" ,"<<py_alpha[2]<<" ,"<<E_alpha[2]<<")"<<endl;
  cout<<"unit is (MeV/c, MeV/c, MeV/c, MeV)"<<endl;

  return 0;
}

