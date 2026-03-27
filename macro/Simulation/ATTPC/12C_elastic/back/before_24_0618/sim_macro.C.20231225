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

int sim_macro(){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
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
  Double_t rmax;
  Double_t r_tmp;
  Double_t ener[2]={0,0};
  Double_t ener_max;
  Double_t tener;
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

  std::vector<Double_t> r;
  
  //input
  TFile *file =new TFile("data/attpcsim.root");  
  TTree *tree = (TTree*)file->Get("cbmsim");
  
  // digi data input
  TChain tpc_tree("cbmsim");
  tpc_tree.Add("./data/output_digi_rcnp.root");   
  TTreeReader reader(&tpc_tree);
  TTreeReaderValue<TClonesArray> event(reader, "AtEventH");
  
  // Create pad plane and load map
  TString mapFile = "Lookup20150611.xml";
  // Set directories
  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString mapDir = dir + "/scripts/" + mapFile;
  
  auto fAtMapPtr = new AtTpcMap();
  fAtMapPtr->ParseXMLMap(mapDir.Data());
  fAtMapPtr->GeneratePadPlane();    // Adam suggest
  auto fPadPlane = fAtMapPtr->GetPadPlane(); // Adam suggest
  
  //set output
  TFile *simmacro=new TFile("data/simmacro.root","recreate");
  //  TTree *PdgCode=new TTree("PdgCode","PdgCode");
  
  TTree *anadata = new TTree("anadata","anadata");

  // new branch for simmacro.root
  anadata->Branch("event_number",&eve,"eve/I");
  anadata->Branch("PdgCode",&pdg,"pdg[2]/I");
  anadata->Branch("startx",&startx,"startx[2]/D");
  anadata->Branch("starty",&starty,"starty[2]/D");
  anadata->Branch("startz",&startz,"startz[2]/D"); 
  anadata->Branch("Energy",&ener,"ener[2]/D");
  anadata->Branch("ener_max",&ener_max,"ener_max/D");
  anadata->Branch("TotalE",&tener,"TotalE/D");
  anadata->Branch("theta_lab",&lab_theta,"theta_lab[2]/D");
  anadata->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  anadata->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");  
  anadata->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");  
  anadata->Branch("n",&n,"n/I");
  
  // the branch from run_digi.root
  anadata->Branch("rmax", &rmax, "rmax/D");  
  //  anadata->Branch("r",&r);
  
  // set reading data
  AtMCTrack* track=new AtMCTrack();
  TClonesArray *trackArray=0;
  tree= (TTree*) file->Get("cbmsim");
  tree->SetBranchAddress("MCTrack",&trackArray);
  Int_t nEvents=tree->GetEntriesFast();
  
  //  cout << nEvents <<endl;

  // analysis of attpcsim.root
  for(int q=0;q < nEvents/2; q++){
    s=2*q+1;
    //    cout<<q<<"  "<<s<<endl;
    tree->GetEvent(s);
    p=trackArray->GetEntries();
    n_d=p;
    //    cout<<"looking for need beam events "<<q<<" "<<n_d<<endl;
  }
  for(int i=0;i < nEvents; i++){
    if((i-1)%1000==0) {
      cout<<"Finish "<<100*(i-1)/nEvents <<" % !"<<endl<<endl;
    }

    eve = i;    

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
    tree->GetEvent(i);
    n=trackArray->GetEntries();
    
    //    cout<<i <<" "<<n<<endl;

    //calculate energy, ang
    //exclude beam and strange event(n>2)
    if (i%2==1 && n==2) {
      for(int j=0;j<n ; j++){	
	track=(AtMCTrack*)trackArray->At(j);
	pdg[j]=track->GetPdgCode();
	startx[j]=track->GetStartX();
	starty[j]=track->GetStartY();
	startz[j]=track->GetStartZ();      
	px[j]=track->GetPx();
	py[j]=track->GetPy();
	pz[j]=track->GetPz();
	
	if (pdg[j]==1000060120) {
	  mass[j]=11177.736;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/mass[j]/2;
	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j] +pz[j]*pz[j])));
	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));
	}
	else if (pdg[j]==2212) {
	  mass[j]=938.7669711;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/h1_mass/2;
	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])));
	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));	  
	}
	else {cout<<i<<" "<<j<<" something is wrong with pdg !!!     "<<endl<<endl;
	}
      }

      //total energy
      tener=(ener[0]+ener[1]);
      if(ener[0] >  ener[1]) ener_max = ener[0];
      if(ener[0] <= ener[1]) ener_max = ener[1];      
      p_beam=sqrt(2*c12_mass*tener/GtoM2);
      v_g=p_beam/(mass[0]+mass[1]);

      //about cms
      for(int k=0; k<n; k++){
	pz_cms[k]=pz[k]-mass[k]*v_g;
	cms_theta[k]=rad_to_deg*acos(pz_cms[k]/(sqrt(px[k]*px[k]+py[k]*py[k]+pz_cms[k]*pz_cms[k])));
	cms_phi[k]=rad_to_deg*acos(px[k]/sqrt(px[k]*px[k]+py[k]*py[k]));	
      }
    }
    //exclude the beam event
    else if (i%2==1 && n==1) {
      continue;
    }
    //set for beam event
    else if (i%2==0 && n==1) {
      
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
      ener[0]=0;
      ener[1]=0;
      ener_max=0;
      tener=0;
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
      
      if (pdg[0]==1000060120) {
      }
      else {cout<<i<<" something is wrong with pdg at beam event !!!     "<<endl<<endl;
      }
    }
    //display n>3 event    
    else {
      cout<<i <<" "<<n<<endl;
      continue;      
    }
    
    // analysis of attpcsim.root to here
    
    // analysis of digi out from here, get hit pattern
    
    std::vector<std::vector<Short_t>> traces; // traces[event][time]
    
    // Get the event to read
    reader.SetEntry(i);
    // Get the event

    AtEvent *eventPtr = (AtEvent *)(event->At(0));
    auto numHits = eventPtr->GetNumHits();
    
    double rmax_tmp=-100;

    for(int m = 0; m<numHits; m++) {
      // Get the padNum for the hit
      auto hit = eventPtr->GetHit(m);
      auto position = hit.GetPosition();      
      auto charge = hit.GetCharge();
      
      //      r[m] = sqrt(position.X()*position.X() + position.Y()*position.Y());
      r_tmp=sqrt(position.X()*position.X() + position.Y()*position.Y());
      r.push_back(r_tmp);
      
      //      if(r[m]>rmax_tmp) rmax_tmp = r[m];
      if(r_tmp>rmax_tmp) rmax_tmp = r_tmp;
    }
    rmax = rmax_tmp;

    if(i%2==1 && n<3) {anadata->Fill();
    }
    //    if(n<3) {anadata->Fill();// need to check; i%2==1 or no need    
    else if(i%2==0 && n_d<3) {anadata->Fill();    
    }
    //  anadata->Fill();
    else{
      continue;
    }
  }
  anadata->AutoSave();
  simmacro->Close();
  cout<<"Finish 100 % !"<<endl<<endl;
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
