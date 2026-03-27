#include <iostream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TEventList.h>
#include <TH1F.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include"TString.h"
#include"TFile.h"
#include"TTree.h"
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TClonesArray.h"
#include"TCanvas.h"
#include"TMath.h"
#include<fstream>
#include<vector>
#include"TChain.h"

int sequential_hoyle_0808(Int_t nEvents=1000, Int_t subnum=0){

  //set situation
  Double_t p_hoyle=0; //z momentum hoyle in MeV
  
  const Double_t AMU = 0.93149432;      // atomic mass unit [GeV]
  const Double_t C = 299792458;       // light speed [m/s]
  const Double_t alpha_mass = 4*AMU+2.4249/1000;     // alpha (= 4He) mass [GeV]
  const Double_t M_4He = alpha_mass;     // alpha (= 4He) mass [GeV]
  const Double_t c12_mass = 12*AMU;      // 12C mass [GeV]
  const Double_t be8_mass = 8*AMU+4.9416/1000; // 8Be mass [GeV]
  const Double_t weight_max = 3;
  const Double_t Ex_hoyle = 0.00765407; // 12C(0+2) excitation energy [GeV]
  const Double_t GtoM=1e3;
  const Double_t GtoM2=1e6;
  const Double_t MtoG=1e-3;
  const Double_t MtoG2=1e-6;
  
  Int_t entry=0;
  Int_t n_decay=0;
  Int_t n_p=0;
  Double_t be8alpha_mass[2] ={be8_mass, alpha_mass};
  Double_t alphaalpha_mass[2] ={alpha_mass, alpha_mass};
  Double_t px_be8, py_be8, pz_be8, E_be8;
  Double_t px_alpha[3];
  Double_t py_alpha[3];
  Double_t pz_alpha[3];
  Double_t nomp_alpha[3];
  Double_t E_alpha[3];
  Double_t E_hoyle=sqrt((c12_mass+Ex_hoyle)*(c12_mass+Ex_hoyle)+p_hoyle*p_hoyle*MtoG2);
  Double_t uniform_ran, p_gx, p_gy, p_gz, p_ga, Esum6, Esum12, Ex12c, k12c;
  Double_t psum12, psum12x, psum12y, psum12z;
  Double_t px[3];
  Double_t py[3];
  Double_t pz[3];
  Double_t E[3];
  Double_t pxg[3];
  Double_t pyg[3];
  Double_t pzg[3];
  Double_t p[3];
  Double_t palpha[3];
  Double_t Ealpha[3];
  Double_t px12[3];
  Double_t py12[3];
  Double_t pz12[3];
  Double_t weight;
  TLorentzVector c12_hoyle(0.0,0.0,p_hoyle,E_hoyle);
  TLorentzVector be8_tmp;
  TLorentzVector *p_alpha[3];
  TGenPhaseSpace event_8BeAlpha,event_2alpha;

  //set histgram
  TH1F *h_ndecay=new TH1F("h_ndecay","number of decay particle",10,0,5);
  TH1F *h_px=new TH1F("h_px","px of alpha",100,0,40);
  TH1F *h_py=new TH1F("h_py","py of alpha",100,0,40);
  TH1F *h_pz=new TH1F("h_pz","pz of alpha",100,0,40);
  //  TH1F *h_p=new TH1F("h_p","p of alpha",100,0,50);
  TH1F *h_p=new TH1F("h_p","energy of alpha",30,0,0.3);
  //  TH1F *h_Eal=new TH1F("h_Eal","Energy of alpha",100,0,10);
  //  TH1F *h_Ealcal=new TH1F("h_Ealcal","calculated Energy of alpha",100,0,10);
  TH1F *h_Ealcal=new TH1F("h_Ealcal","diff of calculated or not Energy of alpha",100,-5,5);
  TH1F *h_pgx=new TH1F("h_pgx","px of CM",100,-0.55,0.5);
  TH1F *h_pgy=new TH1F("h_pgy","py of CM",100,-0.55,0.5);
  TH1F *h_pgz=new TH1F("h_pgz","pz of CM",100,-0.55,0.5);
  TH1F *h_Esa=new TH1F("h_Esa","sum energy of alpha",100,0,1);
  TH1F *h_Ex=new TH1F("h_Ex","Excitation energy of hoyle",100,7,8);
  TH1F *h_k=new TH1F("h_k","kinetic energy of hoyle",100,-1.5,1);

  //output file
  TString outdata="data/output_sequential_"+ to_string(subnum) +".root";
  TFile *tfile=new TFile(outdata,"recreate");
  TTree *tree=new TTree("hoyle","hoyle");
  tree->Branch("momentum_of_hoyle",&p_hoyle,"p_hoyle/I");
  tree->Branch("energy_of_hoyle",&E_hoyle,"E_hoyle/I");
  tree->Branch("entry",&entry,"entry/I");
  tree->Branch("number_of_particles",&n_p,"n_p/I");
  tree->Branch("px_of_be8",&px_be8,"px_be8/D");
  tree->Branch("py_of_be8",&py_be8,"py_be8/D");
  tree->Branch("pz_of_be8",&pz_be8,"pz_be8/D");
  tree->Branch("E_of_be8",&E_be8,"E_be8/D");
  tree->Branch("px_of_alpha",&px_alpha,"px_alpha/D");
  tree->Branch("py_of_alpha",&py_alpha,"py_alpha/D");
  tree->Branch("pz_of_alpha",&pz_alpha,"pz_alpha/D");
  tree->Branch("p_of_alpha",&nomp_alpha,"nomp_alpha/D");
  tree->Branch("E_of_alpha",&E_alpha,"E_alpha/D");
  tree->Branch("E_of_alpha_cal",&E,"E/D");
  tree->Branch("px_center",&p_gx,"p_gx/D");
  tree->Branch("py_center",&p_gy,"p_gy/D");
  tree->Branch("pz_center",&p_gz,"p_gz/D");
  tree->Branch("px_of_CM",&pxg,"pxg/D");
  tree->Branch("py_of_CM",&pyg,"pyg/D");
  tree->Branch("pz_of_CM",&pzg,"pzg/D");
  tree->Branch("sum_E_of_alpha",&Esum6,"Esum6/D");
  tree->Branch("Ex_of_hoyle",&Ex12c,"Ex12c/D");
  tree->Branch("kinemaE",&k12c,"k12c/D");
  tree->Branch("px_of_12c",&px12,"px12/D");
  tree->Branch("py_of_12c",&py12,"py12/D");
  tree->Branch("pz_of_12c",&pz12,"pz12/D");
  
  for(int i=0;i<nEvents;i++){
    // initialization
    n_decay=0;
    n_p=0;
    
    //event generation 12C(hoyle) -> 3-alpha    
    entry++;
    if(event_8BeAlpha.SetDecay(c12_hoyle, 2, be8alpha_mass)){
      do {
	weight = event_8BeAlpha.Generate();
	if (weight > weight_max) {
	  std::cerr << "weight_max is too small!" << std::endl;
	  exit(0);
	}
	uniform_ran = gRandom->Uniform(0.,weight_max);
      } while (uniform_ran > weight);      
      be8_tmp=*event_8BeAlpha.GetDecay(0);
      p_alpha[0]=event_8BeAlpha.GetDecay(1);
      n_decay+=1;

      if(event_2alpha.SetDecay(be8_tmp,2,alphaalpha_mass)){
	do {
	  weight = event_2alpha.Generate();
	  if (weight > weight_max) {
	    std::cerr << "weight_max is too small!" << std::endl;
	    exit(0);
	  }
	  uniform_ran = gRandom->Uniform(0.,weight_max);
	} while (uniform_ran > weight);
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
      px_be8=be8_tmp.Px()*GtoM;
      py_be8=be8_tmp.Py()*GtoM;
      pz_be8=be8_tmp.Pz()*GtoM;
      E_be8 =be8_tmp.E()*GtoM;

      nomp_alpha[j]=sqrt(px_alpha[j]*px_alpha[j]+py_alpha[j]*py_alpha[j]+pz_alpha[j]*pz_alpha[j]);
      
      //     cout<<i<<" "<<j<<" px_alpha "<<px_alpha[j]<<"               py_alpha "<<py_alpha[j]<<"               pz_alpha "<<pz_alpha[j]<<endl;
      //     cout<<i<<" "<<j<<" E_alpha "<<E_alpha[j]<<endl;
      
      h_px->Fill(px_alpha[j]);
      h_py->Fill(py_alpha[j]);
      h_pz->Fill(pz_alpha[j]);
      //      h_p->Fill(nomp_alpha[j]);
      
    }
    
    //copy from check hoyle
    // reconstruct hoyle
    for(int j=0; j<3; j++){
      px[j]=px_alpha[j];
      py[j]=py_alpha[j];
      pz[j]=pz_alpha[j];
      E[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])/(2*alpha_mass*GtoM);
      h_Ealcal->Fill(fabs(E[j]-E_alpha[j]));
      h_p->Fill(E[j]);

      //      cout<<i<<" "<<j<<" px "<<px[j]<<"               py "<<py[j]<<"               pz "<<pz[j]<<endl;
      //      cout<<i<<" "<<j<<"     E_alpha "<<E_alpha[j]<<"       E_alpha "<<E[j]<<endl;
    }
    
    p_gx=(px[0]+px[1]+px[2])/3;
    p_gy=(py[0]+py[1]+py[2])/3;
    p_gz=(pz[0]+pz[1]+pz[2])/3;
    p_ga=sqrt(p_gx*p_gx+p_gy*p_gy+p_gz*p_gz); //abs vector in cms

    //    cout<<i<<" pgx "<<p_gx<<"    pgy "<<p_gy<<"   pgz "<<p_gz<<endl;

    h_pgx->Fill(p_gx);
    h_pgy->Fill(p_gy);
    h_pgz->Fill(p_gz);
    
    for(int j=0; j<3; j++){
      //      h_posi->Fill(z[i],E[i]);
      pxg[j]=px[j]-p_gx;
      pyg[j]=py[j]-p_gy;
      pzg[j]=pz[j]-p_gz;
      p[j]=sqrt(pxg[j]*pxg[j]+pyg[j]*pyg[j]+pzg[j]*pzg[j]);
      Ealpha[j]=E[j]+alpha_mass*GtoM;
      palpha[j]=sqrt(Ealpha[j]*Ealpha[j]-alpha_mass*alpha_mass*GtoM2);

      //      cout<<i<<" "<<j<<" px "<<px[j]<<"    pxg "<<pxg[j]<<endl;
      if(j==0){
	//      cout<<i<<" "<<j<<" p "<<p[j]<<"    kine E "<<E[j]<<endl;
      }
      //      cout<<i<<" "<<j<<" p "<<p[j]<<"    kine E "<<E[j]<<"     palpha from relativity "<<palpha[j]<<endl;
      //      cout<<i<<" "<<j<<"     E_alpha "<<E_alpha[j]<<"     Ealpha "<<Ealpha[j]<<"       E "<<E[j]<<endl;
    }
    
    Esum6=E[0]+E[1]+E[2];
    Esum12=Ealpha[0]+Ealpha[1]+Ealpha[2];
    psum12x=pxg[0]+pxg[1]+pxg[2];
    psum12y=pyg[0]+pyg[1]+pyg[2];
    psum12z=pzg[0]+pzg[1]+pzg[2];
    psum12=sqrt(psum12x*psum12x+psum12y*psum12y+psum12z*psum12z);
    Ex12c=sqrt(Esum12*Esum12-psum12*psum12)-c12_mass*GtoM;
    //    Ex12c=Esum12-c12_mass*GtoM;
    k12c=Esum12-sqrt(Esum12*Esum12-psum12*psum12);

    //    cout<<i<<" "<<" kine E[0] "<<E[0]<<"     Esum3 "<<Esum6<<endl;
    //    cout<<i<<" Esum "<<Esum6<<"    Esum12 "<<Esum12<<"     psum12 "<<psum12<<endl;
    
    px12[0]=(px[0]+px[1]+px[2])/3/alpha_mass*c12_mass;
    py12[0]=(py[0]+py[1]+py[2])/3/alpha_mass*c12_mass;
    pz12[0]=(pz[0]+pz[1]+pz[2])/3/alpha_mass*c12_mass;
    
    n_p=n_decay;
    h_ndecay->Fill(n_decay);
    h_Esa->Fill(Esum6);
    h_Ex->Fill(Ex12c);
    h_k->Fill(k12c);
    tree->Fill();

  }

  //need

  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(3,4);

  h_ndecay->SetStats(1);
  can->cd(1);
  h_ndecay->Draw("hist");

  h_Ealcal->SetStats(1);
  can->cd(2);
  h_Ealcal->Draw("hist");

  h_p->SetStats(1);
  can->cd(3);
  h_p->Draw("hist");

  h_px->SetStats(1);
  can->cd(4);
  h_px->Draw("hist");

  h_py->SetStats(1);
  can->cd(5);
  h_py->Draw("hist");

  h_pz->SetStats(1);
  can->cd(6);
  h_pz->Draw("hist");

  h_pgx->SetStats(1);
  can->cd(7);
  h_pgx->Draw("hist");

  h_pgy->SetStats(1);
  can->cd(8);
  h_pgy->Draw("hist");

  h_pgz->SetStats(1);
  can->cd(9);
  h_pgz->Draw("hist");

  h_Esa->SetStats(1);
  can->cd(10);
  h_Esa->Draw("hist");

  h_Ex->SetStats(1);
  can->cd(11);
  h_Ex->Draw("hist");

  h_k->SetStats(1);
  can->cd(12);
  h_k->Draw("hist");
  
  TString outcan ="data/output_hists_sequential_"+ to_string(subnum) +".C";
  can->SaveAs(outcan);

  //need to here
  
  tree->AutoSave();
  tfile->Close();

  cout<<"output file "<<outdata<<" generated"<<endl;
  
  return 0;
}
