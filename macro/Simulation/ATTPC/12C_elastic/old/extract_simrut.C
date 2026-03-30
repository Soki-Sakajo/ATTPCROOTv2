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

int extract_simrut(Int_t subnum=0){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();

  //from simmacro
  //set constants
  Double_t mass[2];
  Double_t c12_mass=11177.736;  
  Double_t h1_mass=938.7669711;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  int GtoM2=1e6;
  //set variable
  // want to get
  Int_t n;
  Int_t n_d;
  Int_t s;
  Int_t p;
  Int_t eve;
  Int_t pdg[2]={0,0};
  Double_t lab_theta[2];
  Double_t lab_phi[2];
  Double_t p_beam;
  Double_t v_g;
  Double_t pz_cms[2];
  Double_t cms_theta[2];
  Double_t cms_phi[2];
  
  // get from data
  //  start position
  Double_t startx[2];
  Double_t starty[2];
  Double_t startz[2];  
  //  momentum
  Double_t px[2]={0,0};
  Double_t py[2]={0,0};
  Double_t pz[2]={0,0};

  //change

  //set situation
  Int_t z_s=6;          // Z of scatterd particl in simfile
  Int_t z_r1=6;         // Z of an recoil particl in simfile
  Int_t z_r2=1;         // Z of another recoil particl in simfile
  Int_t dz=1;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  Double_t P=0.13;      // atm; gas pressure
  Double_t da=1;        // deg; bin of angle
  Double_t dr=1;        // mm; bin of reaction rate
  Int_t nbin=int(300/dr);    // bin of rmax-reaction_rate hist
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t N_rate2=10;     // (number of recoli partcle)/molecular  
  Int_t theta_min=3;    // deg; rutherford scattar angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  int ctom=10;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; // Pa
  Double_t alfa=7.297352569e-3; // fine structure constant
  Double_t hc=197.463587*1e-15; // MeV•m
  Double_t k_b=1.38065e-23; // J/K ; Bolzman constant
  Double_t rho=P*atmtoPa/(k_b*T);
  Double_t m12=11177.73600; // MeV/c^2 ; 12C mass

  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  
  //set variable
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Double_t rmax;
  Double_t E;
  Double_t theta_cms[2];
  Int_t n_p;

  Int_t n_data;
  Double_t z1[2];
  Double_t theta1_cms_r;
  Double_t Y1_rmax=0;
  Double_t da_r=da*deg_to_rad;
  Double_t l1[div_l];
  Double_t ang1[div_a];
  Double_t ang_r1[div_a];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  
  vector<vector<vector<Int_t>>> ev_ve(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Int_t>>> np_ve(div_l, vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  //  vector<vector<vector<Double_t>>> Y1_ve(div_l, vector<vector<Int_t>>(div_a,vector<Int_t>(0)));  
  
  Int_t z_n=0;
  Int_t a_n=0;
  //  Int_t nEntry=0;
  
  vector<vector<Double_t>> n_a1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> I(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> n_beam1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> sigma1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_da(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_nom(div_l, vector<Double_t>(div_a,0));
  
  Double_t n_react1[div_l];
  Int_t n_rmax1=0;
  Double_t Y1[div_l];
  Double_t Y1_sum=0;
  Double_t z[div_r];
  Double_t z1_max;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1;
  Double_t dir1;
  Double_t Y1_eve;
  Double_t theta1;
  
  bool setbeam=true;
  bool fin_beam1=false;
  bool fin_rea1=false;

  
  //input multi-file by using TChain
  /*
  TChain *tree1=new TChain("cbmsim","cbmsim");
  if(subnum>0){
    for(int i=0; i<subnum; i++){
      //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
      TString dir1="data1";
      TString refi1="attpcsim_12c12c_" + to_string(i) + ".root";
      TString direfi1=dir1+"/"+refi1;
      tree1->Add(direfi1);
      if(i==0){cout<<"read from "<<refi1<<flush;}
      if(i==subnum-1){cout<<"read to "<<refi1<<endl;}
    }
  }
  */
  TChain *tree2=new TChain("anadata","anadata");
  if(subnum>0){
    for(int i=0; i<subnum; i++){
      //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
      TString dir2="data1";
      TString refi2="sim_macro_c12_" + to_string(i) + ".root";
      TString direfi2=dir2+"/"+refi2;
      tree2->Add(direfi2);
      if(i==0){cout<<"read from "<<refi2<<flush;}
      if(i==subnum-1){cout<<"read to "<<refi2<<endl;}
    }
  }
  cout<<"end reading files"<<endl<<endl;
  
  // set reading data
  Int_t Entry2=tree2->GetEntries();
  tree2->SetBranchAddress("event_number",&eve);
  tree2->SetBranchAddress("startz",&startz);
  tree2->SetBranchAddress("rmax",&rmax);
  tree2->SetBranchAddress("TotalE",&E);
  tree2->SetBranchAddress("theta_cms",&theta_cms);
  tree2->SetBranchAddress("n",&n_p);

  //  cout<<"set branch adress"<<endl;
 
  /*
  //set histgram
  TH1F *h_rmax=new TH1F("h_rmax","rmax:reaction_rate",nbin,0,300);
  TH1F *h_z=new TH1F("h_z","reaction position:reaction_rate",nbin_zn,0,400);
  TH2F *h_ztheta=new TH2F("h_ztheta","reaction position:angle",400,0,400,100,0,100);
  TH1F *h_thetaY1=new TH1F("h_thetaY1","angle :reaction_rate",100,0,100);
  TH1F *h_thetaY2=new TH1F("h_thetaY2","angle :reaction_rate",100,0,100);
  TH1F *h_thetaY3=new TH1F("h_thetaY3","angle :reaction_rate",100,0,100);
  */
  
  //change

  //need tp modify
  //set output
  TString outfname = "data1/extract_simmacro_c12_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  TTree *anadata = new TTree("anadata","anadata");
  anadata->Branch("event_number",&eve,"eve/I");
  anadata->Branch("PdgCode",&pdg,"pdg[2]/I");
  anadata->Branch("startx",&startx,"startx[2]/D");
  anadata->Branch("starty",&starty,"starty[2]/D");
  anadata->Branch("startz",&startz,"startz[2]/D"); 
  anadata->Branch("theta_lab",&lab_theta,"theta_lab[2]/D");
  anadata->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  anadata->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");
  anadata->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");  
  anadata->Branch("n",&n,"n/I");

  //change

  //set output for making ratdata by each subnum
  TString outdata = "data1/extract_rutdata_c12_" + to_string(subnum) + "_1.root";
  TFile *rutfile=new TFile(outdata, "recreate");
  TTree *rutdata = new TTree("rutdata","rutdata");
  rutdata->Branch("event_number",&eve,"eve/I");
  rutdata->Branch("PdgCode",&pdg,"pdg2[2]/I");
  rutdata->Branch("startx",&startx,"startx[2]/D");
  rutdata->Branch("starty",&starty,"starty[2]/D");
  rutdata->Branch("startz",&startz,"startz[2]/D"); 
  rutdata->Branch("theta_lab",&lab_theta,"theta_lab[2]/D");
  rutdata->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  rutdata->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");
  rutdata->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");  
  rutdata->Branch("n",&n_p,"n/I");
  TString outdata2 = "data1/extract_rutdata_c12_" + to_string(subnum) + "_2.root";
  TFile *rutfile2=new TFile(outdata2, "recreate");
  TTree *rutdata2 = new TTree("rutdata","rutdata");
  rutdata2->Branch("event_number",&eve,"eve/I");
  rutdata2->Branch("PdgCode",&pdg,"pdg[2]/I");
  rutdata2->Branch("startx",&startx,"startx[2]/D");
  rutdata2->Branch("starty",&starty,"starty[2]/D");
  rutdata2->Branch("startz",&startz,"startz[2]/D"); 
  rutdata2->Branch("theta_lab",&lab_theta,"theta_lab[2]/D");
  rutdata2->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  rutdata2->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");
  rutdata2->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");  
  rutdata2->Branch("n",&n_p,"n/I");
  TString outdata3= "data1/extract_rutdata_c12_" + to_string(subnum) + "_3.root";
  TFile *rutfile3=new TFile(outdata3, "recreate");
  TTree *rutdata3 = new TTree("rutdata","rutdata");
  rutdata3->Branch("event_number",&eve,"eve/I");
  rutdata3->Branch("PdgCode",&pdg,"pdg[2]/I");
  rutdata3->Branch("startx",&startx,"startx[2]/D");
  rutdata3->Branch("starty",&starty,"starty[2]/D");
  rutdata3->Branch("startz",&startz,"startz[2]/D"); 
  rutdata3->Branch("theta_lab",&lab_theta,"theta_lab[2]/D");
  rutdata3->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  rutdata3->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");
  rutdata3->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");  
  rutdata3->Branch("n",&n_p,"n/I");

  /*
  rutdata->Branch("event_number",&eve,"eve/I");
  rutdata->Branch("reaction_position",&dep1,"dep1/D");
  rutdata->Branch("angle",&dir1,"dir1/D");
  rutdata->Branch("rmax_c12",&rmax1,"rmax1/D");
  rutdata->Branch("reaction_rate",&Y1_eve,"Y1_eve/D");
  rutdata->Branch("Eave",&Eave1,"Eave1/D");
  */

  //change

  for(int i=0; i<subnum; i++){
    if(subnum==0){
      cout<<"subnum=0. stop the macro"<<endl<<endl;
      break;
    }
    TString dir1="data1";
    TString refi1="attpcsim_12c12c_" + to_string(i) + ".root";
    TString direfi1=dir1+"/"+refi1;
    TFile *file1 =new TFile(direfi1);
    //    TTree *tree1 = (TTree*)file1->Get("cbmsim");
    if(i==0){cout<<"read from "<<refi1<<endl;}
    if(i==subnum-1){cout<<"read to "<<refi1<<endl;}
   
    // set reading data
    AtMCTrack* track=new AtMCTrack();
    TClonesArray *trackArray=0;
    TTree *tree1= (TTree*) file1->Get("cbmsim");
    cout<<"set before setbranchaddress"<<endl;
    tree1->SetBranchAddress("MCTrack",&trackArray);
    Int_t nEvents=tree1->GetEntriesFast();
    //    Int_t nEvents=2;
    cout<<"set reading data "<<nEvents<<endl;
    
    // analysis of attpcsim.root
    for(int q=0;q < nEvents/2; q++){
      //      s=2*q+1;
      s=2*q+0;
      cout<<q<<"  "<<s<<endl;
      tree1->GetEvent(s);
      p=trackArray->GetEntries();
      n_d=p;
      //    cout<<"looking for need beam events "<<q<<" "<<n_d<<endl;
    }
    for(int j=0;j < nEvents; j++){
      //      cout<<"Finish "<<100*(j-1)/nEvents <<" % !  \r"<<flush;
      
      eve = j;
      //      cout<<eve<<endl;
      // initialization
      pdg[0]=-10;
      pdg[1]=-10;
      
      for(int pernum=0; pernum<2; pernum++){
	startx[pernum]=-1000;
	starty[pernum]=-1000;
	startz[pernum]=-1000;      
      }   
      lab_theta[0]=-1000;
      lab_phi[1]=-1000;    
      cms_theta[0]=-1000;
      cms_phi[1]=-1000;
      p_beam=-1000;
      v_g=-1000;
      pz_cms[0]=-1000;
      pz_cms[1]=-1000;
      
      //get event
      //      tree1->GetEvent(j);
      tree1->GetEntry(j);
      n=trackArray->GetEntries();
      //    cout<<i <<" "<<n<<endl;
      
      //calculate energy, ang
      //exclude beam and strange event(n>2)
      if (j%2==1 && (n>2 ||n==1 )) {
	for(int k=0; k<n ; k++){	
	  track=(AtMCTrack*)trackArray->At(k);
	  pdg[k]=track->GetPdgCode();
	  startx[k]=track->GetStartX();
	  starty[k]=track->GetStartY();
	  startz[k]=track->GetStartZ();      
	  px[k]=track->GetPx();
	  py[k]=track->GetPy();
	  pz[k]=track->GetPz();
	}
      }
      
      //set for beam event
      else if (j%2==0 && n>1) {
	
	track=(AtMCTrack*)trackArray->At(0);
	pdg[0]=track->GetPdgCode();
	pdg[1]=0;
	startx[0]=track->GetStartX();
	startx[1]=0;
	starty[0]=track->GetStartY();
	starty[1]=0;
	startz[0]=track->GetStartZ();
	startz[1]=0;
	px[0]=0;
	py[0]=0;
	pz[0]=0;
	px[1]=0;
	py[1]=0;
	pz[1]=0;
	pz_cms[0]=0;
	pz_cms[1]=0;
	lab_theta[0]=0;
	lab_theta[1]=0;
	lab_phi[0]=0;
	lab_phi[1]=0;
	cms_theta[0]=0;
	cms_theta[1]=0;
	cms_phi[0]=0;
	cms_phi[1]=0;
	
	if (pdg[0]==1000060120) {}
	else {//cout<<i<<" something is wrong with pdg at beam event !!!     "<<endl<<endl;
	}
      }
      //display n>3 event    
      else {
	//      cout<<i <<" "<<n<<endl;
	continue;      
      }
    }
   file1->Close();
  }
  anadata->AutoSave();
  simmacro->Close();
  cout<<"sim_macro Finished"<<endl<<endl;

  //make event lists  
  for(int j=0; j<Entry2; j++){
    tree2->GetEntry(j);
    z1[0]=startz[0]*ctom;
    z1[1]=startz[1]*ctom;

    if(theta1-theta_min>=0){rutdata->Fill();}
    if(n_p!=2){rutdata2->Fill();}
    if(fabs(z1[0]-z1[1])>dz){rutdata3->Fill();}
  }

  /*
  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(2,3);

  // z-angle
  //  h_ztheta->SetStats(0);
  h_ztheta->SetStats(1);
  can->cd(1);
  h_ztheta->Draw("colz");

  // z-reaction rate
  //  h_z->SetStats(0);
  h_z->SetStats(1);
  can->cd(2);
  gPad->SetLogy(1);
  h_z->Draw("hist");

  // rmax-reaction rate
  //  h_rmax->SetStats(0);
  h_rmax->SetStats(1);
  can->cd(3);
  gPad->SetLogy(1);
  h_rmax->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY1->SetStats(1);
  can->cd(4);
  gPad->SetLogy(1);
  h_thetaY1->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY2->SetStats(1);
  can->cd(5);
  gPad->SetLogy(1);
  h_thetaY2->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY3->SetStats(1);
  can->cd(6);
  gPad->SetLogy(1);
  h_thetaY3->Draw("hist");

  //output for canvas
  //  TString outcan = "data/ruther_12c_can_vector.C";
  //  TString outcan = "data/ruther_12c_can_" + to_string(subnum) + ".C";
  TString outcan = "data1/ruther_12c_can_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  h_ztheta->Write();
  h_z->Write();
  h_rmax->Write();
  */
  rutdata->AutoSave();
  rutfile->Close();
  rutdata2->AutoSave();
  rutfile2->Close();
  rutdata3->AutoSave();
  rutfile3->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "output "<<endl;
  cout << outfname <<endl;
  cout<< outdata <<endl;
  cout<< outdata2 <<endl;
  cout<< outdata3 << endl<<endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
