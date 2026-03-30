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

int sim_macro_0912_40matm(Int_t subnum=0){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set constants
  Double_t mass[2];
  Double_t c12_mass=11177.736;  
  Double_t h1_mass=938.7669711;
  Double_t alpha_mass=3728.33692;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t dz=1;
  int GtoM =1e3;
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
  Int_t massnum[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t atomnum[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t massstop=-1000;
  Int_t atomstop=-1000;
  Int_t massstar=-1000;
  Int_t atomstar=-1000;
  Int_t masstemp[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t atomtemp[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t massbefo[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t atombefo[7]={-100,-100,-100,-100,-100,-100,-100};
  Int_t ne=0;
  Int_t ne_b=0;
  Int_t ne_s=0;
  Int_t ne_es=0;
  Int_t ne_eb=0;
  Int_t ne_eo=0;
  Int_t trackid=-100;
  Double_t deltaz;
  Double_t rmax;
  Double_t cmax;
  Double_t cmin;
  Double_t r_tmp;
  Double_t ener[7]={0,0,0,0,0,0,0};
  Double_t ener_max;
  Double_t tener;
  Double_t lab_theta[7];
  Double_t cms_theta[7];
  Double_t bx;
  Double_t by;
  Double_t bz;
  Double_t pba;
  Double_t startx[7];
  Double_t starty[7];
  Double_t startz[7];
  Double_t stopx[7];
  Double_t stopy[7];
  Double_t stopz[7];
  Double_t p_gx;
  Double_t p_gy;
  Double_t p_gz;
  Double_t p_ga;
  Double_t v_gx;
  Double_t v_gy;
  Double_t v_gz;
  Double_t v_g;
  Double_t deltaE;
  Double_t been;
  Double_t momd;
  Double_t pxc[7];
  Double_t pyc[7];
  Double_t pzc[7];
  Double_t pac[7];
  Double_t psp[7];
  Double_t pbx;
  Double_t pby;
  Double_t pbz;
  Double_t px[7]={0,0,0,0,0,0,0};
  Double_t py[7]={0,0,0,0,0,0,0};
  Double_t pz[7]={0,0,0,0,0,0,0};
  Double_t stoppx[7]={0,0,0,0,0,0,0};
  Double_t stoppy[7]={0,0,0,0,0,0,0};
  Double_t stoppz[7]={0,0,0,0,0,0,0};
  Double_t pa[7]={0,0,0,0,0,0,0};
  Double_t sp[7]={0,0,0,0,0,0,0};
  Double_t r_nd[7]={0,0,0,0,0,0,0};
  Double_t rmax_nd=0;
  bool beam_start=true;
  bool scat_start[7]={true,true,true,true,true,true,true};
  bool check_alpha[7]={false,false,false,false,false,false,false};
  bool check_6alpha=false;

  //input
  //  TString infname = "data/attpcsim_12c12c_" + to_string(subnum) + ".root";
  TString infname = "data1/attpcsim_40matm_" + to_string(subnum) + ".root";
  TFile *file =new TFile(infname);  
  TTree *tree = (TTree*)file->Get("cbmsim");
  
  // digi data input
  TChain tpc_tree("cbmsim");
  //  TString digifname ="data/output_digi_rcnp_12c_" + to_string(subnum) +".root";
  TString digifname ="data1/output_digi_40matm_" + to_string(subnum) + ".root";
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
  
  // set reading data
  AtMCTrack* track=new AtMCTrack();
  TClonesArray *trackArray=0;
  tree= (TTree*) file->Get("cbmsim");
  tree->SetBranchAddress("MCTrack",&trackArray);

  //get AtMCPoint data
  AtMCPoint* point = new AtMCPoint();    
  TClonesArray *pointArray=0;
  tree->SetBranchAddress("AtTpcPoint",&pointArray);
  Int_t nEvents=tree->GetEntriesFast();
  Int_t n_d[nEvents];
  cout<<endl;

  //set output
  //  TString outfname = "data/sim_macro_12c_" + to_string(subnum) + ".root";
  TString outfname = "data1/sim_macro_40matm_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  TTree *anadata = new TTree("anadata","anadata");
  TTree *expbeam = new TTree("expbeam","expbeam");
  TTree *expscat = new TTree("expscat","expscat");
  TTree *expothe = new TTree("expothe","expothe");

  // new branch for simmacro.root
  anadata->Branch("event_number",&eve,"eve/I");
  anadata->Branch("startx",&startx,"startx[7]/D");
  anadata->Branch("starty",&starty,"starty[7]/D");
  anadata->Branch("startz",&startz,"startz[7]/D");
  anadata->Branch("stopx",&stopx,"stopx[7]/D");
  anadata->Branch("stopy",&stopy,"stopy[7]/D");
  anadata->Branch("stopz",&stopz,"stopz[7]/D");
  anadata->Branch("massnum",&massnum,"massnum[7]/I");
  anadata->Branch("atomnum",&atomnum,"atomnum[7]/I");
  anadata->Branch("startpx",&px,"px[7]/D");
  anadata->Branch("startpy",&py,"py[7]/D");
  anadata->Branch("startpz",&pz,"pz[7]/D");
  anadata->Branch("stoppx",&stoppx,"stoppx[7]/D");
  anadata->Branch("stoppy",&stoppy,"stoppy[7]/D");
  anadata->Branch("stoppz",&stoppz,"stoppz[7]/D");
  anadata->Branch("energy",&ener,"ener[7]/D");
  anadata->Branch("energy_max",&ener_max,"ener_max/D");
  anadata->Branch("beam_energy",&tener,"TotalE/D");
  //  anadata->Branch("theta_lab",&lab_theta,"theta_lab[7]/D");
  //  anadata->Branch("theta_cms",&cms_theta,"theta_cms[7]/D");
  anadata->Branch("n",&n,"n/I");
  anadata->Branch("rmax_nd", &rmax_nd, "rmax_nd/D");
  // the branch from run_digi.root
  anadata->Branch("rmax", &rmax, "rmax/D");
  anadata->Branch("cmax", &cmax, "cmax/D");
  anadata->Branch("cmin", &cmin, "cmin/D");
  // the branch in other tree
  expbeam->Branch("event_number_beam",&eenb,"eenb/I");
  expbeam->Branch("number_of_particle_beam",&enpb,"enpb/I");
  expscat->Branch("event_number_scat",&eens,"eens/I");
  expscat->Branch("number_of_particle_scat",&enps,"enps/I");
  expothe->Branch("event_number_other",&eeno,"eeno/I");
  expothe->Branch("number_of_particle_other",&enpo,"enpo/I");
  
  // analysis of attpcsim.root
  for(int i=0;i < nEvents; i++){
    tree->GetEvent(i);
    p=trackArray->GetEntries();
    n_d[i]=p;
    if(i%100==0) {
      //need
      cout<<"Get data "<<100*i/nEvents <<" % ! \r"<<flush;
    }
  }
  cout<<"Get data 100 % ! \r"<<flush;
  for(int i=0;i < nEvents; i++){
  //  for(int i=0;i < 100; i++){
    if(i%1000==0) {
      //need
      cout<<"Finish "<<100*i/nEvents <<" % !     \r"<<flush;
    }
    eve = i;    
    ne++;
    // initialization
    for(int pernum=0; pernum<2; pernum++){
      startx[pernum]=-1000;
      starty[pernum]=-1000;
      startz[pernum]=-1000;
      stopx[pernum]=-1000;
      stopy[pernum]=-1000;
      stopz[pernum]=-1000;
      r_nd[pernum]=-1000;
      massnum[pernum]=-100;
      atomnum[pernum]=-100;
    }
    //    lab_theta[0]=-1000;
    //    cms_theta[0]=-1000;
    p_gx=-1000;
    p_gy=-1000;
    p_gz=-1000;
    pxc[0]=-1000;
    pyc[0]=-1000;
    pzc[0]=-1000;
    pxc[1]=-1000;
    pyc[1]=-1000;
    pzc[1]=-1000;
    rmax_nd=-1000;
    
    //get event
    tree->GetEvent(i);
    n=trackArray->GetEntries();
    np=pointArray->GetEntries();
    //    cout<<i <<" "<<n<<endl;

    //calculate energy, ang
    //set for beam event
    if (i%2==0 && n==1) {
      ne_b++;
      point=(AtMCPoint*)pointArray->At(0); // get the point information
      trackid=(int)(point->GetTrackID()); //get track ID
      if(trackid==0){
	massstar=point->GetMassNum();
	atomstar=point->GetAtomicNum();
	//cout<<i<<"  mass number "<<massstar<<endl;
	if(massstar==12 && atomstar==6){
	  startx[0]=point->GetX();
	  starty[0]=point->GetY();
	  startz[0]=point->GetZ();
	  px[0]=point->GetPx();
	  py[0]=point->GetPy();
	  pz[0]=point->GetPz();
	  ener[0]=(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0])*GtoM2/c12_mass/2;
	  ener_max=0;
	  tener=0;
	  pxc[0]=0;
	  pyc[0]=0;
	  pzc[0]=0;
	  //	  lab_theta[0]=0;
	  //	  cms_theta[0]=0;
	  //	  for(int j=1; j<7; j++){
	  //7->6
	  for(int j=1; j<6; j++){
	    startx[j]=0;
	    starty[j]=0;
	    startz[j]=0;
	    px[j]=0;
	    py[j]=0;
	    pz[j]=0;
	    ener[j]=0;
	    pxc[j]=0;
	    pyc[j]=0;
	    pzc[j]=0;
	    //	    cms_theta[j]=0;
	    //	    lab_theta[j]=0;
	  }
	}
      }
      else {cout<<i<<" something is wrong with track ID at start point in beam event !!!     "<<endl<<endl;}

      for(int k=0; k<np; k++){
	point=(AtMCPoint*)pointArray->At(k); // get the point information
	trackid=(int)(point->GetTrackID()); //get track ID
	massstop=point->GetMassNum();
	atomstop=point->GetAtomicNum();
	if(trackid==0 && massstop==12 && atomstop==6 ){
	  /*
	    if(beam_start){
	    massbefo[j]=masstemp[j];
	    atombefo[j]=atomtemp[j];
	    }
	  */
	  //	cout<<i<<"  mass number "<<massstop<<endl;
	  //	  for(int j=1; j<7; j++){
	  //7->6
	  for(int j=1; j<6; j++){
	    stopx[j]=0;
	    stopy[j]=0;
	    stopz[j]=0;
	    stoppx[j]=0;
	    stoppy[j]=0;
	    stoppz[j]=0;
	    r_nd[j]=0;
	    if(massstar==massstop && atomstar==atomstop){
	      massnum[j]=0;
	      atomnum[j]=0;
	    }
	  }
	  stopx[0]=point->GetX();
	  stopy[0]=point->GetY();
	  stopz[0]=point->GetZ();
	  r_nd[0]=sqrt(stopx[0]*stopx[0]+stopy[0]*stopy[0]);
	  rmax_nd=r_nd[0];
	  bx=stopx[0];
	  by=stopy[0];
	  bz=stopz[0];
	  pbx=point->GetPx();
	  pby=point->GetPy();
	  pbz=point->GetPz();
	  stoppx[0]=pbx;
	  stoppy[0]=pby;
	  stoppz[0]=pbz;
	  pba=sqrt(pbx*pbx+pby*pby+pbz*pbz);
	  if(massstar==massstop && atomstar==atomstop){
	    massnum[0]=massstop;
	    atomnum[0]=atomstop;
	  }
	  else {cout<<i<<"  "<<k<<" something is wrong with mass and atomic number at start and stop in beam event !!!"<<endl<<endl;
	  }
	}
	else {cout<<i<<" something is wrong with track ID at end point in beam event !!!     "<<endl<<endl;}
      }
    }

    //exclude beam and strange event(n>2)
    //    else if (i%2==1 && n==7 && n_d[i-1]==1) {
    //7->6
    else if (i%2==1 && n==6 && n_d[i-1]==1) {
      for(int j=0; j<n ; j++){	
	scat_start[j]=true;
	check_alpha[j]=false;
	for(int k=0; k<np ; k++){
	  point=(AtMCPoint*)pointArray->At(k);
	  trackid=(int)(point->GetTrackID());
	  if(trackid==j){
	    masstemp[j]=point->GetMassNum();
	    atomtemp[j]=point->GetAtomicNum();
	    if(scat_start[j]){
	      massbefo[j]=masstemp[j];
	      atombefo[j]=atomtemp[j];
	      startx[j]=point->GetX();
	      starty[j]=point->GetY();
	      startz[j]=point->GetZ();
	      px[j]=point->GetPx();
	      py[j]=point->GetPy();
	      pz[j]=point->GetPz();
	      scat_start[j]=false;
	    }
	    stopx[j]=point->GetX();
	    stopy[j]=point->GetY();
	    stopz[j]=point->GetZ();
	    stoppx[j]=point->GetPx();
	    stoppy[j]=point->GetPy();
	    stoppz[j]=point->GetPz();
	    r_nd[j]=10*sqrt(stopx[j]*stopx[j]+stopy[j]*stopy[j]);
	    if(masstemp[j]==massbefo[j] && atomtemp[j]==atombefo[j]){
	      massnum[j]=masstemp[j];
	      atomnum[j]=atomtemp[j];
	    }
	    else{cout<<i<<"  "<<j<<" something wrong with mass and atomic number in scattering event !!!   "<<endl<<endl;}
	    massbefo[j]=masstemp[j];
	    atombefo[j]=atomtemp[j];
	  }
	}
	if (massnum[j]==12 && atomnum[j]==6) {
	  mass[j]=c12_mass;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/mass[j]/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	}
	else if (massnum[j]==1 && atomnum[j]==1) {
	  mass[j]=h1_mass;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/h1_mass/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	}
	else if (massnum[j]==4 && atomnum[j]==2) {
	  mass[j]=alpha_mass;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/alpha_mass/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	  check_alpha[j]=true;
	}
	else {cout<<i<<" "<<j<<" something is wrong with mass and atomic number !!!     "<<endl<<endl;}
      }
      rmax_nd=0;
      ener_max=0;
      for(Int_t j=1; j<n; j++){
	if(rmax_nd<r_nd[j]){rmax_nd=r_nd[j];}
	if(ener_max<ener[j]){ener_max=ener[j];}
      }
      //beam energy
      check_6alpha=true;
      for(int j=1; j<n; j++){
	if(check_alpha[j]==false){
	  check_6alpha=false;
	}
      }

      if(check_6alpha){
	p_gx=(px[1]+px[2]+px[3]+px[4]+px[5]+px[6])*GtoM/6;
	p_gy=(py[1]+py[2]+py[3]+py[4]+py[5]+py[6])*GtoM/6;
	p_gz=(pz[1]+pz[2]+pz[3]+pz[4]+pz[5]+pz[6])*GtoM/6;
	v_gx=p_gx/alpha_mass;
	v_gy=p_gy/alpha_mass;
	v_gz=p_gz/alpha_mass;
	v_g=sqrt(v_gx*v_gx+v_gy*v_gy+v_gz*v_gz);
	tener=(2*v_g)*(2*v_g)*c12_mass/2;
      }
      else{
	cout<<i<<" something is wrong with check 6alpha !!!     "<<endl;
	for(int j=1; j<n; j++){
	  if(check_alpha[j]==false){
	    cout<<i<<"  "<<j<<" check alpha false !!      "<<endl;
	  }
	}
      }

      /*      
      tener=(ener[1]+ener[2]+ener[3]
	     +ener[4]+ener[5]+ener[6]);
      p_gx=c12_mass*pbx/(mass[0]+mass[1]);
      p_gy=c12_mass*pby/(mass[0]+mass[1]);
      p_gz=c12_mass*pbz/(mass[0]+mass[1]);
      p_ga=sqrt(p_gx*p_gx+p_gy*p_gy+p_gz*p_gz); //abs vector in cms
      //about angle in cms
      for(int k=0; k<n; k++){
	pxc[k]=px[k]-p_gx;
	pyc[k]=py[k]-p_gy;
	pzc[k]=pz[k]-p_gz;
	pac[k]=sqrt(pxc[k]*pxc[k]+pyc[k]*pyc[k]+pzc[k]*pzc[k]); //abs vector in cms
	psp[k]=p_gx*pxc[k]+p_gy*pyc[k]+p_gz*pzc[k]; // scalar product
	cms_theta[k]=rad_to_deg*acos(psp[k]/(p_ga*pac[k]));
      }
      */
    }

    //display n>3 event    
    //exclude the beam event
    //    else if (i%2==1 && ((n>1 && n<7) || n>7)) {
    //7->6
    else if (i%2==1 && ((n>1 && n<6) || n>6)) {
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
    double cmax_tmp=-100;
    double cmin_tmp=1e7;

    for(int m = 0; m<numHits; m++) {
      // Get the padNum for the hit
      auto hit = eventPtr->GetHit(m);
      auto position = hit.GetPosition();      
      auto charge = hit.GetCharge();
      
      r_tmp=sqrt(position.X()*position.X() + position.Y()*position.Y());
      if(r_tmp>rmax_tmp) rmax_tmp = r_tmp;
      if(charge>cmax_tmp) cmax_tmp = charge;
      if(charge<cmin_tmp) cmin_tmp = charge;
    }
    rmax = rmax_tmp;
    cmax = cmax_tmp;
    cmin = cmin_tmp;

    //    if(i%2==1 && n<3) {
    //7->6
    //    if(i%2==1 && n==7) {
    if(i%2==1 && n==6) {
      ne_s++;
      anadata->Fill();
    }
    else if(i%2==0 && n<2) {anadata->Fill();}
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
