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

int sim_macro_12c12c(Int_t subnum=0){ 

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
  Int_t eve;
  Int_t eenb;
  Int_t eens;
  Int_t eeno;
  Int_t n;
  Int_t np;
  Int_t enpb;
  Int_t enps;
  Int_t enpo;
  Int_t s;
  Int_t p;
  Int_t pdg[2]={0,0};
  Int_t ne=0;
  Int_t ne_b=0;
  Int_t ne_s=0;
  Int_t ne_es=0;
  Int_t ne_eb=0;
  Int_t ne_eo=0;
  Double_t rmax;
  Double_t r_tmp;
  Double_t ener[2]={0,0};
  Double_t ener_max;
  Double_t tener;
  Double_t lab_theta[2];
  //  Double_t lab_phi[2];
  
  Double_t pbx;
  Double_t pby;
  Double_t pbz;
  Double_t pba;
  Double_t v_gx;
  Double_t v_gy;
  Double_t v_gz;
  Double_t v_ga;
  Double_t vxc[2];
  Double_t vyc[2];
  Double_t vzc[2];
  Double_t vac[2];
  Double_t vsp[2];
  Double_t cms_theta[2];
  //  Double_t cms_phi[2];
  
  // get from data
  //  start position
  Double_t startx[2];
  Double_t starty[2];
  Double_t startz[2];  
  //  momentum
  Double_t px[2]={0,0};
  Double_t py[2]={0,0};
  Double_t pz[2]={0,0};
  Double_t pa[2]={0,0};
  Double_t sp[2]={0,0};

  std::vector<Double_t> r;
  
  //input
  TString infname = "data1/attpcsim_12c12c_" + to_string(subnum) + ".root";
  TFile *file =new TFile(infname);  
  TTree *tree = (TTree*)file->Get("cbmsim");
  
  // digi data input
  TChain tpc_tree("cbmsim");
  TString digifname = "data1/output_digi_rcnp_12c_" + to_string(subnum) + ".root";
  tpc_tree.Add(digifname);   
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
  TString outfname = "data1/sim_macro_c12_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  //  TTree *PdgCode=new TTree("PdgCode","PdgCode");
  TTree *anadata = new TTree("anadata","anadata");
  TTree *expbeam = new TTree("expbeam","expbeam");
  TTree *expscat = new TTree("expscat","expscat");
  TTree *expothe = new TTree("expothe","expothe");

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
  anadata->Branch("theta_cms",&cms_theta,"theta_cms[2]/D");
  anadata->Branch("n",&n,"n/I");
  //  anadata->Branch("phi_lab",&lab_phi,"phi_lab[2]/D");
  //  anadata->Branch("phi_cms",&cms_phi,"phi_cms[2]/D");

  expbeam->Branch("event_number_beam",&eenb,"eenb/I");
  expbeam->Branch("number_of_particle_beam",&enpb,"enpb/I");
  expscat->Branch("event_number_scat",&eens,"eens/I");
  expscat->Branch("number_of_particle_scat",&enps,"enps/I");
  expothe->Branch("event_number_other",&eeno,"eeno/I");
  expothe->Branch("number_of_particle_other",&enpo,"enpo/I");
  
  // the branch from run_digi.root
  anadata->Branch("rmax", &rmax, "rmax/D");  
  //  anadata->Branch("r",&r);
  
  // set reading data
  AtMCTrack* track=new AtMCTrack();
  TClonesArray *trackArray=0;
  tree= (TTree*) file->Get("cbmsim");
  tree->SetBranchAddress("MCTrack",&trackArray);

  //get AtMCPoint data
  FairMCPoint* point = new FairMCPoint();    
  TClonesArray *pointArray=0;
  tree->SetBranchAddress("AtTpcPoint",&pointArray);

  Int_t nEvents=tree->GetEntriesFast();
  Int_t n_d[nEvents];

  //  cout << nEvents <<endl;
  
  // analysis of attpcsim.root
  for(int i=0;i < nEvents; i++){
    //    s=2*i+1;
    //    cout<<i<<"  "<<s<<endl;
    //    tree->GetEvent(s);
    tree->GetEvent(i);
    p=trackArray->GetEntries();
    n_d[i]=p;
    //    cout<<"looking for need beam events "<<i<<" "<<n_d<<endl;
  }
  cout<<endl;
  for(int i=0;i < nEvents; i++){
  //  for(int i=0;i < 100; i++){
    if((i-1)%1000==0) {
      //need
      cout<<"Finish "<<100*(i-1)/nEvents <<" % ! \r"<<flush;
    }
    eve = i;    
    ne++;
    // initialization
    pdg[0]=-10;
    pdg[1]=-10;

    for(int pernum=0; pernum<2; pernum++){
      startx[pernum]=-1000;
      starty[pernum]=-1000;
      startz[pernum]=-1000;      
    }   
    lab_theta[0]=-1000;
    cms_theta[0]=-1000;
    //    lab_phi[1]=-1000;
    //    cms_phi[1]=-1000;
    v_gx=-1000;
    v_gy=-1000;
    v_gz=-1000;
    vxc[0]=-1000;
    vyc[0]=-1000;
    vzc[0]=-1000;
    vxc[1]=-1000;
    vyc[1]=-1000;
    vzc[1]=-1000;
    
    //get event
    tree->GetEvent(i);
    n=trackArray->GetEntries();
    np=pointArray->GetEntries();
    //    cout<<i <<" "<<n<<endl;

    //calculate energy, ang
    //set for beam event
    if (i%2==0 && n==1) {
      track=(AtMCTrack*)trackArray->At(0);
      pdg[0]=track->GetPdgCode();
      pdg[1]=0;
      if (pdg[0]==1000060120) {
	ne_b++;
	startx[0]=track->GetStartX();
	startx[1]=0;
	starty[0]=track->GetStartY();
	starty[1]=0;
	startz[0]=track->GetStartZ();
	startz[1]=0;
	px[0]=track->GetPx();
	py[0]=track->GetPy();
	pz[0]=track->GetPz();
	px[1]=0;
	py[1]=0;
	pz[1]=0;
	ener[0]=(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0])*GtoM2/c12_mass/2;
	ener[1]=0;
	ener_max=0;
	tener=0;
	vxc[0]=0;
	vyc[0]=0;
	vzc[0]=0;
	vxc[1]=0;
	vyc[1]=0;
	vzc[1]=0;
	lab_theta[0]=0;
	lab_theta[1]=0;
	cms_theta[0]=0;
	cms_theta[1]=0;
	point=(FairMCPoint*)pointArray->At(np-1);
	pbx=point->GetPx();
	pby=point->GetPy();
	pbz=point->GetPz();
	pba=sqrt(pbx*pbx+pby*pby+pbz*pbz);
	/*
	Double_t xb=point->GetX();
	Double_t yb=point->GetY();
	Double_t zb=point->GetZ();
	
	if(i==3284 || i==100){
	  cout<<i<<" pbx "<<pbx<<" pby "<<pby<<" pbz "<<pbz<<endl;
	  cout<<i<<" xb "<<xb<<" yb "<<yb<<" zb "<<zb<<endl;
	}
	*/
      }
      else {cout<<i<<" something is wrong with pdg at beam event !!!     "<<endl<<endl;}
    }

    //exclude beam and strange event(n>2)
    
    //    if (i%2==1 && n==2) {
    else if (i%2==1 && n==2 && n_d[i-1]==1) {
	for(int j=0;j<n ; j++){	
	track=(AtMCTrack*)trackArray->At(j);
	pdg[j]=track->GetPdgCode();
	startx[j]=track->GetStartX();
	starty[j]=track->GetStartY();
	startz[j]=track->GetStartZ();      
	px[j]=track->GetPx();
	py[j]=track->GetPy();
	pz[j]=track->GetPz();

	/*
	if(i==3285|| i==101){
	  if(j==0){cout<<i<<" px[0] "<<px[0]<<"  py[0]  "<<py[0]<<"  pz[0]  "<<pz[0]<<endl;}
	  else if(j==1){
	    cout<<i<<" x[0] "<<startx[0]<<"  y[0]  "<<starty[0]<<"  z[0]  "<<startz[0]<<endl;
	    cout<<i<<" x[1] "<<startx[1]<<"  y[1]  "<<starty[1]<<"  z[1]  "<<startz[1]<<endl;
	    cout<<i<<" px[1] "<<px[1]<<"  py[1]  "<<py[1]<<"  pz[1]  "<<pz[1]<<endl;
	    cout<<i<<" sum px "<<px[0]+px[1]<<"  sum py  "<<py[0]+py[1]<<" sum pz "<<pz[0]+pz[1]<<endl;
	    cout<<i<<" pbx-px "<<pbx-(px[0]+px[1])<<"  pby-py  "<<pby-(py[0]+py[1])<<" pbz-pz "<<pbz-(pz[0]+pz[1])<<endl<<endl;
	  }
	}
	*/
	if (pdg[j]==1000060120) {
	  mass[j]=11177.736;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/mass[j]/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	  //	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j] +pz[j]*pz[j])));
	  //	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));
	}
	else if (pdg[j]==2212) {
	  mass[j]=938.7669711;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/h1_mass/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	  //	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])));
	  //	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));	  
	}
	else {cout<<i<<" "<<j<<" something is wrong with pdg !!!     "<<endl<<endl;
	}
      }

      //total energy
      tener=(ener[0]+ener[1]);
      if(ener[0] >  ener[1]) ener_max = ener[0];
      if(ener[0] <= ener[1]) ener_max = ener[1];      

      v_gx=pbx/(mass[0]+mass[1]);
      v_gy=pby/(mass[0]+mass[1]);
      v_gz=pbz/(mass[0]+mass[1]);
      v_ga=sqrt(v_gx*v_gx+v_gy*v_gy+v_gz*v_gz); //abs vector in cms
      
      //about cms
      for(int k=0; k<n; k++){
	vxc[k]=px[k]/mass[k]-v_gx;
	vyc[k]=py[k]/mass[k]-v_gy;
	vzc[k]=pz[k]/mass[k]-v_gz;
	vac[k]=sqrt(vxc[k]*vxc[k]+vyc[k]*vyc[k]+vzc[k]*vzc[k]); //abs vector in cms
	vsp[k]=v_gx*vxc[k]+v_gy*vyc[k]+v_gz*vzc[k]; // scalar product
	
	cms_theta[k]=rad_to_deg*acos(vsp[k]/(v_ga*vac[k]));
	//	cms_theta[k]=rad_to_deg*acos(pz_cms[k]/(sqrt(px[k]*px[k]+py[k]*py[k]+pz_cms[k]*pz_cms[k])));
	//	cms_phi[k]=rad_to_deg*acos(px[k]/sqrt(px[k]*px[k]+py[k]*py[k]));	
      }
    }

    //display n>3 event    
    //exclude the beam event
    else if (i%2==1 && n>2) {
      eens=i;
      enps=n;
      ne_es++;
      expscat->Fill();
      continue;
    }
    else if (i%2==0) {
      eenb=i;
      enpb=n;
      ne_eb++;
      expbeam->Fill();
      continue;
    }
    else {
      eeno=i;
      enpo=n;
      ne_eo++;
      expothe->Fill();
      //      cout<<i <<" "<<n<<endl;
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

    if(i%2==1 && n<3) {
      ne_s++;
      anadata->Fill();
    }
    //    if(n<3) {anadata->Fill();// need to check; i%2==1 or no need    
    //    else if(i%2==0 && n_d<3) {anadata->Fill();}
    else if(i%2==0 && n<2) {anadata->Fill();}
    //  anadata->Fill();
    else{continue;}
  }
  anadata->AutoSave();
  expscat->AutoSave();
  expbeam->AutoSave();
  expothe->AutoSave();
  simmacro->Close();
  cout<<"Finish 100 % !"<<endl<<endl;
  cout<<"all event "<<ne<<"   beam event "<<ne_b<<"   scattaring event "<<ne_s<<endl;
  cout<<"except event "<<ne_es+ne_eb+ne_eo<<endl;
  cout<<"except scattering event "<<ne_es<<endl;
  cout<<"except beam event "<<ne_eb<<endl;
  cout<<"except other problem event "<<ne_eo<<endl;
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
