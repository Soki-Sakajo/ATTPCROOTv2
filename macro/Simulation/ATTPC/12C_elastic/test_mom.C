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

int test_mom(Int_t subnum=0){ 

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
  Double_t dz=1;
  
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
  //  Int_t pdg[2]={0,0};
  Int_t massnum[2]={-100,-100};
  Int_t atomnum[2]={-100,-100};
  Int_t massstop=-1000;
  Int_t atomstop=-1000;
  Int_t massstar=-1000;
  Int_t atomstar=-1000;
  Int_t masstemp[2]={-100,-100};
  Int_t atomtemp[2]={-100,-100};
  Int_t massbefo[2]={-100,-100};
  Int_t atombefo[2]={-100,-100};
  Int_t ne=0;
  Int_t ne_b=0;
  Int_t ne_s=0;
  Int_t ne_es=0;
  Int_t ne_eb=0;
  Int_t ne_eo=0;
  Int_t ne_p1=0;
  Int_t ne_p2=0;
  Int_t trackid=-100;
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
  Double_t bx;
  Double_t by;
  Double_t bz;
  Double_t pba;
  Double_t p_gx;
  Double_t p_gy;
  Double_t p_gz;
  Double_t p_ga;
  Double_t deltaE;
  Double_t been;
  Double_t momd;
  Double_t pxc[2];
  Double_t pyc[2];
  Double_t pzc[2];
  Double_t pac[2];
  Double_t psp[2];
  Double_t cms_theta[2];
  //  Double_t cms_phi[2];
  
  // get from data
  //  start position
  Double_t startx[2];
  Double_t starty[2];
  Double_t startz[2];
  //  stop position
  Double_t stopx[2];
  Double_t stopy[2];
  Double_t stopz[2];
  //  momentum
  Double_t px[2]={0,0};
  Double_t py[2]={0,0};
  Double_t pz[2]={0,0};
  Double_t pa[2]={0,0};
  Double_t sp[2]={0,0};
  Double_t r_nd[2]={0,0};
  Double_t rmax_nd=0;
  bool beam_start=true;
  bool scat_start[2]={true,true};
  Double_t deltaz;  
  
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

  // set reading data
  AtMCTrack* track=new AtMCTrack();
  TClonesArray *trackArray=0;
  tree= (TTree*) file->Get("cbmsim");
  tree->SetBranchAddress("MCTrack",&trackArray);
  //  tree->SetBranchAddress("AtMCTrack",&trackArray);

  //get AtMCPoint data
  AtMCPoint* point = new AtMCPoint();
  TClonesArray *pointArray=0;
  tree->SetBranchAddress("AtTpcPoint",&pointArray);

  Int_t nEvents=tree->GetEntriesFast();
  Int_t n_d[nEvents];

  Int_t n_bin=450;
  TH1F *h_z=new TH1F("h_z","delta-z :count",n_bin, -50, 400);
  TH2F *h_z2=new TH2F("h_z2","z:delta-z",400,0,400,n_bin, -50, 400);
  TH1F *h_mom=new TH1F("h_mom","sum(abs(pb-p))",100,0,2e-4);
  TH1F *h_ener=new TH1F("h_ener","delta-E",100,-5e-5,5e-5);
  TH2F *h_de=new TH2F("h_de","Ebeam:deltaE",50,0,0.05,100, -5e-5, 5e-5);

  //  cout << nEvents <<endl;
  cout<<endl;
  
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
    //    pdg[0]=-10;
    //    pdg[1]=-10;

    for(int pernum=0; pernum<2; pernum++){
      startx[pernum]=-1000;
      starty[pernum]=-1000;
      startz[pernum]=-1000;
      stopx[pernum]=-1000;
      stopy[pernum]=-1000;
      stopz[pernum]=-1000;
      r_nd[pernum]=-1000;
    }   
    lab_theta[0]=-1000;
    cms_theta[0]=-1000;
    //    lab_phi[1]=-1000;
    //    cms_phi[1]=-1000;
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
      
      //      track=(AtMCTrack*)trackArray->At(0);
      //      pdg[0]=track->GetPdgCode();
      //      pdg[1]=0;
      //      if (pdg[0]==1000060120) {
      //      ne_b++;
      point=(AtMCPoint*)pointArray->At(0); // get the point information
      trackid=(int)(point->GetTrackID()); //get track ID
      if(trackid==0){
	massstar=point->GetMassNum();
	atomstar=point->GetAtomicNum();
	//cout<<i<<"  mass number "<<massstar<<endl;
	if(massstar==12 && atomstar==6){
	  startx[0]=point->GetX();
	startx[1]=0;
	starty[0]=point->GetY();
	starty[1]=0;
	startz[0]=point->GetZ();
	startz[1]=0;
	px[0]=point->GetPx();
	py[0]=point->GetPy();
	pz[0]=point->GetPz();
	px[1]=0;
	py[1]=0;
	pz[1]=0;
	ener[0]=(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0])*GtoM2/c12_mass/2;
	ener[1]=0;
	ener_max=0;
	tener=0;
	pxc[0]=0;
	pyc[0]=0;
	pzc[0]=0;
	pxc[1]=0;
	pyc[1]=0;
	pzc[1]=0;
	lab_theta[0]=0;
	lab_theta[1]=0;
	cms_theta[0]=0;
	cms_theta[1]=0;
	}
      }
      else {cout<<i<<" something is wrong with track ID at start point in beam event !!!     "<<endl<<endl;}

      for(int k=0; k<np; k++){
	point=(AtMCPoint*)pointArray->At(k); // get the point information
      trackid=(int)(point->GetTrackID()); //get track ID
      massstop=point->GetMassNum();
      atomstop=point->GetAtomicNum();
      //      if(trackid==0){
     if(trackid==0 && massstop==12 && atomstop==6 ){
	/*
	if(beam_start){
	      massbefo[j]=masstemp[j];
	      atombefo[j]=atomtemp[j];
	}
	*/
	//	cout<<i<<"  mass number "<<massstop<<endl;
	stopx[0]=point->GetX();
	stopx[1]=0;
	stopy[0]=point->GetY();
	stopy[1]=0;
	stopz[0]=point->GetZ();
	stopz[1]=0;
	r_nd[0]=sqrt(stopx[0]*stopx[0]+stopy[0]*stopy[0]);
	r_nd[1]=0;
	rmax_nd=r_nd[0];
	bx=stopx[0];
	by=stopy[0];
	bz=stopz[0];
	pbx=point->GetPx();
	pby=point->GetPy();
	pbz=point->GetPz();
	pba=sqrt(pbx*pbx+pby*pby+pbz*pbz);
	if(massstar==massstop && atomstar==atomstop){
	  massnum[0]=massstop;
	  massnum[1]=0;
	  atomnum[0]=atomstop;
	  atomnum[1]=0;
	}
	else {cout<<i<<"  "<<k<<" something is wrong with mass and atomic number at start and stop in beam event !!!"<<endl<<endl;
	}

      }
      else {cout<<i<<" something is wrong with track ID at end point in beam event !!!     "<<endl<<endl;}
      /*
	if(i==3284 || i==100){
	cout<<i<<" pbx "<<pbx<<" pby "<<pby<<" pbz "<<pbz<<endl;
	cout<<i<<" xb "<<stopx[0]<<" yb "<<stopy[0]<<" zb "<<stopz[0]<<endl;
	}
      */
      //    }
      //      else {cout<<i<<" something is wrong with pdg at beam event !!!     "<<endl<<endl;}
      }
      /*
      if(i==12){
	  cout<<i<<" np "<<np<<"     "<<endl;
	    cout<<i<<" pbx "<<pbx<<" pby "<<pby<<" pbz "<<pbz<<endl;
	    cout<<i<<" xb "<<bx<<" yb "<<by<<" zb "<<bz<<endl;
	    cout<<endl;
	}
      */
    }

    //exclude beam and strange event(n>2)
    //    if (i%2==1 && n==2) {
    else if (i%2==1 && n==2 && n_d[i-1]==1) {
      for(int j=0; j<n ; j++){	
	scat_start[0]=true;
	scat_start[1]=true;
	//	track=(AtMCTrack*)trackArray->At(j);
	//	pdg[j]=track->GetPdgCode();
	
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
	if(j==1 && fabs(startz[0]-startz[1])<dz/10){
	  deltaz=(bz-startz[0])*10;
	  h_z->Fill(deltaz);
	  h_z2->Fill(startz[0]*10,deltaz);
	  }
	else if(j==0){}
	else{cout<<i<<" per partilce's startz is not same                                           "<<endl<<endl;}

	if((bz>startz[0] || bz>startz[1]) && j==1){
	//	if((bz>startz[0] || bz>startz[1]) && j==1 && bz>36.6){
	  ne_p1++;
	  //	if(j==1){
	  //	  cout<<i<<" bz "<<bz<<"   startz    "<<startz[0]<<"      "<<startz[1]<<endl;
	}

	if(j==1){
	  momd=fabs(pbx-(px[0]+px[1]))+fabs(pby-(py[0]+py[1]))+fabs(pbz-(pz[0]+pz[1]));
	  been=(pbx*pbx+pby*pby+pbz*pbz)/(2*c12_mass*1e-3);
	  deltaE=(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0]+px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1])/(2*c12_mass*1e-3)-(pbx*pbx+pby*pby+pbz*pbz)/(2*c12_mass*1e-3);
	  h_ener->Fill(deltaE,1);
	  h_de->Fill(been,deltaE);
	  h_mom->Fill(momd,1);
	  //	if(i==1 && j==1){
	  //	  if(fabs(pbx-(px[0]+px[1]))+fabs(pby-(py[0]+py[1]))+fabs(pbz-(pz[0]+pz[1]))>1.6e-4){
	  if(fabs(pbx-(px[0]+px[1]))+fabs(pby-(py[0]+py[1]))+fabs(pbz-(pz[0]+pz[1]))>1.6e-4){
	    ne_p2++;

	    cout<<i-1<<" xb "<<bx<<" yb "<<by<<" zb "<<bz<<endl;
	    cout<<i-1<<" pbx "<<pbx<<" pby "<<pby<<" pbz "<<pbz<<endl;
	    cout<<i-1<<" Eb "<<(pbx*pbx+pby*pby+pbz*pbz)/(2*c12_mass*1e-3)<<" GeV"<<endl;
	    cout<<i<<" x[0] "<<startx[0]<<"  y[0]  "<<starty[0]<<"  z[0]  "<<startz[0]<<endl;
	    cout<<i<<" x[1] "<<startx[1]<<"  y[1]  "<<starty[1]<<"  z[1]  "<<startz[1]<<endl;
	    cout<<i<<" px[0] "<<px[0]<<"  py[0]  "<<py[0]<<"  pz[0]  "<<pz[0]<<endl;
	    cout<<i<<" px[1] "<<px[1]<<"  py[1]  "<<py[1]<<"  pz[1]  "<<pz[1]<<endl;
	    cout<<i<<" sum px "<<px[0]+px[1]<<"  sum py  "<<py[0]+py[1]<<" sum pz "<<pz[0]+pz[1]<<endl;
	    cout<<i<<" E "<<(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0]+px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1])/(2*c12_mass*1e-3)<<" GeV"<<endl;
	    cout<<i<<" delta E "<<(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0]+px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1])/(2*c12_mass*1e-3)-(pbx*pbx+pby*pby+pbz*pbz)/(2*c12_mass*1e-3)<<" GeV"<<endl;
	    cout<<i<<" pbx-px "<<pbx-(px[0]+px[1])<<"  pby-py  "<<pby-(py[0]+py[1])<<" pbz-pz "<<pbz-(pz[0]+pz[1])<<endl;
	    cout<<i<<"  sum(abs(pb-px))  "<<fabs(pbx-(px[0]+px[1]))+fabs(pby-(py[0]+py[1]))+fabs(pbz-(pz[0]+pz[1]))<<endl<<endl;
	    cout<<"event "<<i<<" nmom "<<ne_p2<<endl<<endl;

	  }
	  }

	/*	
	if (massnum[j]==12 && atomnum[j]==6) {
	  mass[j]=c12_mass;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/mass[j]/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	  //	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j] +pz[j]*pz[j])));
	  //	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));
	}
	else if (massnum[j]==1 && atomnum[j]==1) {
	  mass[j]=h1_mass;
	  ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])*GtoM2/h1_mass/2;
	  pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
	  sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
	  //	  lab_theta[j]=rad_to_deg*acos(pz[j]/(sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])));
	  //	  lab_phi[j]=rad_to_deg*acos(px[j]/sqrt(px[j]*px[j]+py[j]*py[j]));	  
	}
	else {cout<<i<<" "<<j<<" something is wrong with mass and atomic number !!!     "<<endl<<endl;}
      }
      if(r_nd[0]>=r_nd[1]){rmax_nd=r_nd[0];}
      else if(r_nd[0]<r_nd[1]){rmax_nd=r_nd[1];}
      
      //total energy
      tener=(ener[0]+ener[1]);
      if(ener[0] >  ener[1]) ener_max = ener[0];
      if(ener[0] <= ener[1]) ener_max = ener[1];      
      
      p_gx=c12_mass*pbx/(mass[0]+mass[1]);
      p_gy=c12_mass*pby/(mass[0]+mass[1]);
      p_gz=c12_mass*pbz/(mass[0]+mass[1]);
      p_ga=sqrt(p_gx*p_gx+p_gy*p_gy+p_gz*p_gz); //abs vector in cms
      
      //about cms
      for(int k=0; k<n; k++){
	pxc[k]=px[k]-p_gx;
	pyc[k]=py[k]-p_gy;
	pzc[k]=pz[k]-p_gz;
	pac[k]=sqrt(pxc[k]*pxc[k]+pyc[k]*pyc[k]+pzc[k]*pzc[k]); //abs vector in cms
	psp[k]=p_gx*pxc[k]+p_gy*pyc[k]+p_gz*pzc[k]; // scalar product
	
	cms_theta[k]=rad_to_deg*acos(psp[k]/(p_ga*pac[k]));
	//	cms_theta[k]=rad_to_deg*acos(pz_cms[k]/(sqrt(px[k]*px[k]+py[k]*py[k]+pz_cms[k]*pz_cms[k])));
	//	cms_phi[k]=rad_to_deg*acos(px[k]/sqrt(px[k]*px[k]+py[k]*py[k]));	
	*/
      }
    }
  }

  for(int m=0; m<=n_bin; m++){
    double err_z=sqrt(h_z->GetBinContent(m));
    h_z->SetBinError(m,err_z);
  }
  TCanvas *can=new TCanvas("hist","hist");
  can->Divide(2,3);
  can->cd(1);
  gPad->SetLogy(1);
  h_z->Draw("hist");

  can->cd(2);
  //  gPad->SetLogy(1);
  h_z2->Draw("colz");

  can->cd(3);
  //  gPad->SetLogy(1);
  h_mom->Draw("hist");

  can->cd(4);
  //  gPad->SetLogy(1);
  h_ener->Draw("hist");

  can->cd(5);
  //  gPad->SetLogy(1);
  h_de->Draw("colz");

  cout<<"nz "<<ne_p1<<"  nmom "<<ne_p2<<"    "<<endl<<endl;
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;

}
