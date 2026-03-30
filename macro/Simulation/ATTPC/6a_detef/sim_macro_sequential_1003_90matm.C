#define SAVE_FILE
#define CANVAS

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

int sim_macro_sequential_1003_90matm(Int_t subnum=0){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  Double_t mass[2];
  //set situation
  const Int_t dz=1;           // mm; bin of depth
  const Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  const Int_t r_attpc=275;    // mm; radius of attpc
  const Int_t l_attpc=1000;   // mm; depth of attpc
  const Int_t div_l=(int)(l_attpc/dz+0.9);

  //set constants
  const Int_t ctom=10;
  const Int_t GtoM =1e3;
  const Int_t GtoM2=1e6;
  const Double_t amu=931.49432; //atomic mass unit in MeV
  const Double_t c12_mass=12*amu;     //MeV
  const Double_t mg24_mass=24*amu-13.933578;    //MeV
  const Double_t h1_mass=amu+7.28897;    //MeV
  const Double_t alpha_mass=4*amu+2.4249;  //MeV
  const Double_t pi=3.1415926535;
  const Double_t rad_to_deg=180/pi;
  
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
  Int_t massnum[6]={-100,-100,-100,-100,-100,-100};
  Int_t atomnum[6]={-100,-100,-100,-100,-100,-100};
  Int_t massstop=-1000;
  Int_t atomstop=-1000;
  Int_t massstar=-1000;
  Int_t atomstar=-1000;
  Int_t masstemp[6]={-100,-100,-100,-100,-100,-100};
  Int_t atomtemp[6]={-100,-100,-100,-100,-100,-100};
  Int_t massbefo[6]={-100,-100,-100,-100,-100,-100};
  Int_t atombefo[6]={-100,-100,-100,-100,-100,-100};
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
  Double_t ener[6]={0,0,0,0,0,0};
  Double_t ener_max;
  Double_t tener;
  Double_t lab_theta[6];
  Double_t cms_theta[6];
  Double_t bx;
  Double_t by;
  Double_t bz;
  Double_t pba;
  Double_t startx[6];
  Double_t starty[6];
  Double_t startz[6];
  Double_t stopx[6];
  Double_t stopy[6];
  Double_t stopz[6];
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
  Double_t pxc[6];
  Double_t pyc[6];
  Double_t pzc[6];
  Double_t pac[6];
  Double_t psp[6];
  Double_t pbx;
  Double_t pby;
  Double_t pbz;
  Double_t px[6]={0,0,0,0,0,0};
  Double_t py[6]={0,0,0,0,0,0};
  Double_t pz[6]={0,0,0,0,0,0};
  Double_t stoppx[6]={0,0,0,0,0,0};
  Double_t stoppy[6]={0,0,0,0,0,0};
  Double_t stoppz[6]={0,0,0,0,0,0};
  Double_t pa[6]={0,0,0,0,0,0};
  Double_t sp[6]={0,0,0,0,0,0};
  Double_t r_nd[6]={0,0,0,0,0,0};
  Double_t rmax_nd=0;
  Double_t Esum6=0;
  Double_t Esum12[2]={0,0};
  Double_t Ex12c[2]={0,0};
  Double_t k12c[2]={0,0};
  Double_t x[6];
  Double_t y[6];
  Double_t z[6];
  Double_t pxg[6];
  Double_t pyg[6];
  Double_t pzg[6];
  Double_t psum12x[6];
  Double_t psum12y[6];
  Double_t psum12z[6];
  Double_t psum12[6];
  Double_t p_alpha[6];
  Double_t palpha[6];
  Double_t px12[2]={0,0};
  Double_t py12[2]={0,0};
  Double_t pz12[2]={0,0};
  Double_t E[6];
  Double_t Ealpha[6];
  bool beam_start=true;
  bool scat_start[6]={true,true,true,true,true,true};
  bool check_alpha[6]={false,false,false,false,false,false};
  bool check_6alpha=false;
  bool check_6alpha_event=false;
  bool flag_dx6=false;
  bool flag_dy6=false;
  bool flag_dz6=false;

  //in and out dire set
  //  TString outdir= "./data/";
  TString outdir= "./data1/";

  //input
  //  TString infname = "data/attpcsim_12c12c_" + to_string(subnum) + ".root";
  //  TString infname = "data1/attpcsim_hoylehoyle_" + to_string(subnum) + ".root";
  //  TString infname = "data/attpcsim_sequential_hoylehoyle_" + to_string(subnum) + ".root";
  TString infname =
    outdir + "attpcsim_sequential_hoylehoyle_90matm_" + to_string(subnum) + ".root";
  TFile *file =new TFile(infname);  
  TTree *tree = (TTree*)file->Get("cbmsim");

  // digi data input
  TChain tpc_tree("cbmsim");
  //  TString digifname ="data/output_digi_rcnp_12c_nd_"+ to_string(subnum) +".root";
  TString digifname =
    //    "data1/output_digi_rcnp_hoylehoyle_nd_" + to_string(subnum) + ".root";
    //    "data/output_digi_rcnp_sequential_hoylehoyle_nd_" + to_string(subnum) + ".root";
    outdir + "output_digi_rcnp_sequential_hoylehoyle_90matm_" + to_string(subnum) + ".root";
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
  fAtMapPtr->GeneratePadPlane();
  auto fPadPlane = fAtMapPtr->GetPadPlane();
  
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
  #ifdef SAVE_FILE
  //  TString outfname = "data1/check_hoyle_dire_" + to_string(subnum) + ".root";
  //  TString outfname = "data1/sim_macro_hoylehoyle_nd_" + to_string(subnum) + ".root";
  //  TString outfname = "data/sim_macro_sequential_hoylehoyle_nd_" + to_string(subnum) + ".root";
  TString outfname =
    outdir + "sim_macro_sequential_hoylehoyle_90matm_" + to_string(subnum) + ".root";
  TFile *simmacro=new TFile(outfname, "recreate");
  TTree *anadata = new TTree("anadata","anadata");
  TTree *expbeam = new TTree("expbeam","expbeam");
  TTree *expscat = new TTree("expscat","expscat");
  TTree *expothe = new TTree("expothe","expothe");

  // new branch for simmacro.root
  //  the branch from attpcsim
  anadata->Branch("event_number",&eve,"eve/I");
  anadata->Branch("n",&n,"n/I");
  anadata->Branch("massnum",&massnum,"massnum[6]/I");
  anadata->Branch("atomnum",&atomnum,"atomnum[6]/I");
  anadata->Branch("startx",&startx,"startx[6]/D");
  anadata->Branch("starty",&starty,"starty[6]/D");
  anadata->Branch("startz",&startz,"startz[6]/D");
  anadata->Branch("startpx",&px,"px[6]/D");
  anadata->Branch("startpy",&py,"py[6]/D");
  anadata->Branch("startpz",&pz,"pz[6]/D");
  anadata->Branch("stopx",&stopx,"stopx[6]/D");
  anadata->Branch("stopy",&stopy,"stopy[6]/D");
  anadata->Branch("stopz",&stopz,"stopz[6]/D");
  anadata->Branch("stoppx",&stoppx,"stoppx[6]/D");
  anadata->Branch("stoppy",&stoppy,"stoppy[6]/D");
  anadata->Branch("stoppz",&stoppz,"stoppz[6]/D");
  anadata->Branch("rmax_nd", &rmax_nd, "rmax_nd/D");
  anadata->Branch("energy",&ener,"ener[6]/D");
  anadata->Branch("energy_max",&ener_max,"ener_max/D");
  anadata->Branch("beam_energy",&tener,"TotalE/D");
  //  the branch about reconstracted hoyle
  anadata->Branch("px_of_CM",&p_gx,"pgx/D");
  anadata->Branch("py_of_CM",&p_gy,"pgy/D");
  anadata->Branch("pz_of_CM",&p_gz,"pgz/D");
  anadata->Branch("sum_energy",&Esum6,"Esum6/D");
  anadata->Branch("sum_energy_of_hoyle",&Esum12,"Esum12c/D");
  anadata->Branch("Ex_energy",&Ex12c,"Ex12c/D");
  anadata->Branch("ki_energy",&k12c,"k12c/D");
  anadata->Branch("px_of_hoyle",&px12,"px12c/D");
  anadata->Branch("py_of_hoyle",&py12,"py12c/D");
  anadata->Branch("pz_of_hoyle",&pz12,"pz12c/D");
  anadata->Branch("theta_of_hoyle",&cms_theta,"cms_theta/D");
  //  the branch from run_digi.root
  anadata->Branch("rmax", &rmax, "rmax/D");
  anadata->Branch("cmax", &cmax, "cmax/D");
  anadata->Branch("cmin", &cmin, "cmin/D");
  //  the branch in other tree
  expbeam->Branch("event_number_beam",&eenb,"eenb/I");
  expbeam->Branch("number_of_particle_beam",&enpb,"enpb/I");
  expscat->Branch("event_number_scat",&eens,"eens/I");
  expscat->Branch("number_of_particle_scat",&enps,"enps/I");
  expothe->Branch("event_number_other",&eeno,"eeno/I");
  expothe->Branch("number_of_particle_other",&enpo,"enpo/I");
  #endif
  
  //set histgram
  TH2F *h_posi=new TH2F("h_posi","reaction position:energy of alpha",nbin_zn,0,400,150,0,15);
  TH2F *h_dep=new TH2F("h_dep","reaction position: total energy of 6 alpha",200,0,200,40,0,40);
  TH2F *h_beam=new TH2F("h_beam","reaction position:beam energy",200,0,200,520,0,52);
  //  TH2F *h_hoy=new TH2F("h_hoy","reaction position: exenergy of hoyle",200,0,200,50,7,8);
  TH2F *h_hoy=new TH2F("h_hoy","reaction position: exenergy of hoyle",200,0,200,50,7.64,7.71);
  //  TH2F *h_hang=new TH2F("h_hang","reaction angle: energy of hoyle",180,0,180,30,5,10);
  TH2F *h_hang=new TH2F("h_hang","reaction angle: energy of hoyle",180,0,180,50,7.64,7.71);
  TH2F *h_h1=new TH2F("h_h1","reaction position: reaction angle of hoyle",400,0,400,180,0,180);
  //  TH1F *h_ex=new TH1F("h_ex","excitation energy",100,7.6,7.8);
  TH1F *h_ex=new TH1F("h_ex","diff of excitation energy",100,-0.05,0.051);
  TH2F *h_ang=new TH2F("h_ang","reaction angle: kienergy of hoyle",180,0,180,350,0,35);
  TH2F *h_ang1=new TH2F("h_ang1","reaction angle: energy of hoyle",180,0,180,350,0,35);
  //  TH2F *h_ex_dep=new TH2F("h_ex_dep","reaction position: excitation energy in dep range",20,100,110,50,5,10);
  TH1F *h_ex_dep=new TH1F("h_ex_dep","excitation energy of hoyle1",100,7.64,7.74);
  //  TH2F *h_ex_ang=new TH2F("h_ex_ang","reaction angle: excitation energy in dep range",180,0,180,20,5,10);
  TH1F *h_ex_ang=new TH1F("h_ex_ang","excitation energy of hoyle2",100,7.64,7.74);

  //  TH2F *h_rmax=new TH2F("h_rmax","reaction position: rmax",nbin_zn,0,400,160,0,160);
  TH2F *h_rmax=new TH2F("h_rmax","reaction position: (Ebeam-Esum6alpha)/2",180,0,180,100,7.27,7.285);
  //  TH2F *h_rmaxnd=new TH2F("h_rmaxnd","reaction position: rmax_nd",nbin_zn,0,400,160,0,160);
  TH2F *h_rmaxnd=new TH2F("h_rmaxnd","reaction position: Ebeam-Esum6alpha",180,0,180,150,14.54,14.57);

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
    flag_dx6=true;
    flag_dy6=true;
    flag_dz6=true;

    if(i%1000==0) {
      //need
      cout<<"Finish "<<100*i/nEvents <<" % !     \r"<<flush;
    }
    eve = i;    
    ne++;
    // initialization
    for(int pernum=0; pernum<6; pernum++){
      startx[pernum]=-1000;
      starty[pernum]=-1000;
      startz[pernum]=-1000;
      stopx[pernum]=-1000;
      stopy[pernum]=-1000;
      stopz[pernum]=-1000;
      r_nd[pernum]=-1000;
      massnum[pernum]=-100;
      atomnum[pernum]=-100;
      pxc[pernum]=-1000;
      pyc[pernum]=-1000;
      pzc[pernum]=-1000;
    }
    p_gx=-1000;
    p_gy=-1000;
    p_gz=-1000;
    rmax_nd=-1000;
    
    //get event
    tree->GetEvent(i);
    n=trackArray->GetEntries();
    np=pointArray->GetEntries();

    //calculate energy, ang
    // getting data of beam event
    if (i%2==0 && n==1) {
      ne_b++;
      point=(AtMCPoint*)pointArray->At(0); // get the point information
      trackid=(int)(point->GetTrackID()); //get track ID
      if(trackid==0){
        massstar=point->GetMassNum();
        atomstar=point->GetAtomicNum();
        if(massstar==12 && atomstar==6){
          startx[0]=point->GetX()*ctom;
          starty[0]=point->GetY()*ctom;
          startz[0]=point->GetZ()*ctom;
          px[0]=point->GetPx()*GtoM;
          py[0]=point->GetPy()*GtoM;
          pz[0]=point->GetPz()*GtoM;
          ener[0]=(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0])/c12_mass/2;
          ener_max=0;
          tener=0;
          pxc[0]=0;
          pyc[0]=0;
          pzc[0]=0;
	  //          for(int j=1; j<7; j++){
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
	  //          for(int j=1; j<7; j++){
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
          stopx[0]=point->GetX()*ctom;
          stopy[0]=point->GetY()*ctom;
          stopz[0]=point->GetZ()*ctom;
          r_nd[0]=sqrt(stopx[0]*stopx[0]+stopy[0]*stopy[0]);
          rmax_nd=r_nd[0];
          bx=stopx[0];
          by=stopy[0];
          bz=stopz[0];
          pbx=point->GetPx()*GtoM;
          pby=point->GetPy()*GtoM;
          pbz=point->GetPz()*GtoM;
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

    // getting data of 6alpha event
    //    else if (i%2==1 && n==7 && n_d[i-1]==1) {
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
              startx[j]=point->GetX()*ctom;
              starty[j]=point->GetY()*ctom;
              startz[j]=point->GetZ()*ctom;
              px[j]=point->GetPx()*GtoM;
              py[j]=point->GetPy()*GtoM;
              pz[j]=point->GetPz()*GtoM;
              scat_start[j]=false;
            }
            stopx[j]=point->GetX()*ctom;
            stopy[j]=point->GetY()*ctom;
            stopz[j]=point->GetZ()*ctom;
            stoppx[j]=point->GetPx()*GtoM;
            stoppy[j]=point->GetPy()*GtoM;
            stoppz[j]=point->GetPz()*GtoM;
            r_nd[j]=sqrt(stopx[j]*stopx[j]+stopy[j]*stopy[j]);
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
          ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])/mass[j]/2;
          pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
          sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
          //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
        }
        else if (massnum[j]==1 && atomnum[j]==1) {
          mass[j]=h1_mass;
          ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])/h1_mass/2;
          pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
          sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
          //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
        }
        else if (massnum[j]==4 && atomnum[j]==2) {
          mass[j]=alpha_mass;
          ener[j]=(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j])/alpha_mass/2;
          pa[j]=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); //abs vector
          sp[j]=pbx*px[j]+pby*py[j]+pbz*pz[j];  //scalar product
          //	  lab_theta[j]=rad_to_deg*acos(sp[j]/(pba*pa[j]));
          check_alpha[j]=true;
        }
        else {cout<<i<<" "<<j<<" something is wrong with mass and atomic number !!!     "<<endl<<endl;}
      }
      rmax_nd=0;
      ener_max=0;
      for(Int_t j=0; j<n; j++){
        if(rmax_nd<r_nd[j]){rmax_nd=r_nd[j];}
        if(ener_max<ener[j]){ener_max=ener[j];}
      }
      //check 6alpha
      check_6alpha=true;
      for(int j=0; j<n; j++){
        if(check_alpha[j]==false){
          check_6alpha=false;
        }
      }

      //reconsutract hoyle state
      for(int l=0; l<6; l++){
        if(l<5){
          if(fabs(startx[l]-startx[l+1])>1){
            flag_dx6=false;
            cout<<i<<"  "<<l<<"  "<<startx[l]<<"  "<<fabs(startx[l]-startx[l+1])<<endl;
          }
          if(fabs(starty[l]-starty[l+1])>1){
            flag_dy6=false;
            cout<<i<<"  "<<l<<"  "<<starty[l]<<"  "<<fabs(starty[l]-starty[l+1])<<endl;
          }
          if(fabs(startz[l]-startz[l+1])>dz){
            flag_dz6=false;
            cout<<i<<"  "<<l<<"  "<<startz[l]<<"  "<<fabs(startz[l]-startz[l+1])<<endl;
          }
        }
        else if(l==5){
          if(fabs(startx[l]-startx[0])>1){
            flag_dx6=false;
          }
          if(fabs(starty[l]-starty[0])>1){
            flag_dy6=false;
          }
          if(fabs(startz[l]-startz[0])>dz){
            flag_dz6=false;
          }
        }
      }
      if(check_6alpha && flag_dx6 && flag_dy6 && flag_dz6){
        p_gx=(px[0]+px[1]+px[2]+px[3]+px[4]+px[5])/6;
        p_gy=(py[0]+py[1]+py[2]+py[3]+py[4]+py[5])/6;
        p_gz=(pz[0]+pz[1]+pz[2]+pz[3]+pz[4]+pz[5])/6;
        p_ga=sqrt(p_gx*p_gx+p_gy*p_gy+p_gz*p_gz);     //abs vector in cms
        v_gx=p_gx/alpha_mass;
        v_gy=p_gy/alpha_mass;
        v_gz=p_gz/alpha_mass;
        v_g=sqrt(v_gx*v_gx+v_gy*v_gy+v_gz*v_gz);
        tener=(2*v_g)*(2*v_g)*c12_mass/2;
        h_beam->Fill(startz[0],tener);

        for(int l=0; l<6; l++){
          E[l]=ener[l];
          h_posi->Fill(startz[l],E[l]);
          pxg[l]=px[l]-p_gx;
          pyg[l]=py[l]-p_gy;
          pzg[l]=pz[l]-p_gz;
          p_alpha[l]=sqrt(pxg[l]*pxg[l]+pyg[l]*pyg[l]+pzg[l]*pzg[l]);
          Ealpha[l]=E[l]+alpha_mass;
          palpha[l]=sqrt(Ealpha[l]*Ealpha[l]-alpha_mass*alpha_mass);
        }

        Esum6=E[0]+E[1]+E[2]+E[3]+E[4]+E[5];
	h_rmax->Fill(startz[1],(tener-Esum6)/2);
	h_rmaxnd->Fill(startz[1],tener-Esum6);
        //need to check the delta E of hoyle and 3alpha
        //hoyle 1
        Esum12[0]=Ealpha[0]+Ealpha[1]+Ealpha[2];
        psum12x[0]=(px[0]+px[1]+px[2]);
        psum12y[0]=(py[0]+py[1]+py[2]);
        psum12z[0]=(pz[0]+pz[1]+pz[2]);
        psum12[0]=sqrt(psum12x[0]*psum12x[0]+psum12y[0]*psum12y[0]+psum12z[0]*psum12z[0]);
        Ex12c[0]=sqrt(Esum12[0]*Esum12[0]-psum12[0]*psum12[0])-c12_mass;
        k12c[0]=Esum12[0]-sqrt(Esum12[0]*Esum12[0]-psum12[0]*psum12[0]);
        //hoyle 2
        Esum12[1]=Ealpha[3]+Ealpha[4]+Ealpha[5];
        psum12x[1]=(px[3]+px[4]+px[5]);
        psum12y[1]=(py[3]+py[4]+py[5]);
        psum12z[1]=(pz[3]+pz[4]+pz[5]);
        psum12[1]=sqrt(psum12x[1]*psum12x[1]+psum12y[1]*psum12y[1]+psum12z[1]*psum12z[1]);
        Ex12c[1]=sqrt(Esum12[1]*Esum12[1]-psum12[1]*psum12[1])-c12_mass;
        k12c[1]=Esum12[1]-sqrt(Esum12[1]*Esum12[1]-psum12[1]*psum12[1]);

        px12[0]=(px[0]+px[1]+px[2])/3/alpha_mass*c12_mass;
        py12[0]=(py[0]+py[1]+py[2])/3/alpha_mass*c12_mass;
        pz12[0]=(pz[0]+pz[1]+pz[2])/3/alpha_mass*c12_mass;
        px12[1]=(px[3]+px[4]+px[5])/3/alpha_mass*c12_mass;
        py12[1]=(py[3]+py[4]+py[5])/3/alpha_mass*c12_mass;
        pz12[1]=(pz[3]+pz[4]+pz[5])/3/alpha_mass*c12_mass;
        /*
        if((Ex12c[0]<7.65 || Ex12c[0]>7.66)||(Ex12c[1]<7.65 || Ex12c[1]>7.66)){
        //if((Ex12c[0]>7.65 && Ex12c[0]<7.66)&&(Ex12c[1]>7.65 && Ex12c[1]<7.66)){
          cout<<j<<"  excitation 0 "<<Ex12c[0]<<"  excitation 1 "<<Ex12c[1]<<endl;
        }
        */
        //	cout<<i<<"  "<<j<<"  "<<k12c[0]<<"     "<<Esum12[0]<<"     "<<psum12[0]<<endl;
        //about angle in cms
        for(int k=0; k<2; k++){
          pxc[k]=px12[k]-p_gx/alpha_mass*c12_mass;
          pyc[k]=py12[k]-p_gy/alpha_mass*c12_mass;
          pzc[k]=pz12[k]-p_gz/alpha_mass*c12_mass;
          pac[k]=sqrt(pxc[k]*pxc[k]+pyc[k]*pyc[k]+pzc[k]*pzc[k]);//abs vector in cms
          psp[k]=p_gx*pxc[k]+p_gy*pyc[k]+p_gz*pzc[k]; // scalar product
          cms_theta[k]=rad_to_deg*acos(psp[k]/(p_ga*pac[k]));
          h_ang->Fill(cms_theta[k],k12c[k]);
	  //          h_ex->Fill(Ex12c[k]);
	  //hist condition
	  /*
          if(tener>=45 && tener<50){
            h_ang1->Fill(cms_theta[k],k12c[k]);
          }
	  if(fabs(Ex12c[k]-7.654)<0.001){
	    h_hoy->Fill(startz[1],Ex12c[k]);
	    h_hang->Fill(cms_theta[k],Ex12c[k]);
	    h_ex_dep->Fill(Ex12c[k]);
	  }
	  if(fabs(Ex12c[k]-7.654)<0.001 && k==0){
	    h_ex->Fill(Ex12c[0]-Ex12c[1]);
	  }
	  */

	  h_ang1->Fill(cms_theta[k],k12c[k]);
	  h_hoy->Fill(startz[1],Ex12c[k]);
	  h_hang->Fill(cms_theta[k],Ex12c[k]);
	  h_ex_dep->Fill(Ex12c[k]);
	  h_ex->Fill(Ex12c[0]-Ex12c[1]);
	  
        }
        h_dep->Fill(startz[1],Esum6);
	h_ex_ang->Fill(Ex12c[0]);
        h_h1->Fill(startz[1],cms_theta[0]);
        h_h1->Fill(startz[1],cms_theta[1]);
	/*
        h_hang->Fill(cms_theta[0],Ex12c[0]);
        h_hang->Fill(cms_theta[1],Ex12c[1]);
	h_ex_dep->Fill(Ex12c[0]);
	h_ex_ang->Fill(Ex12c[1]);
	h_ex->Fill(Ex12c[0]-Ex12c[1]);
	*/
      }
      else{
        check_6alpha_event=true;
        cout<<i<<" something is wrong with check 6alpha !!!     "<<endl;
        for(int j=0; j<n; j++){
          if(check_alpha[j]==false){
            check_6alpha_event=false;
            cout<<i<<"  "<<j<<" check alpha false !!      "<<endl;
          }
        }
        if(check_6alpha_event) {
          if(flag_dx6==false || flag_dy6==false){
            cout<<i<<" this event is starting x or y position is not same "<<endl;
          }
          if(flag_dz6==false){
            cout<<i<<" this event is starting z position is not same  "<<endl;
            cout<<" z[0] "<<z[0]<<" z[1] "<<z[1]<<" z[2] "<<z[2]<<" z[3] "<<z[3]
            <<" z[4] "<<z[4]<<" z[5] "<<z[5]<<endl;
          }
	  else{
	    cout<<i<<" something is wrong  check_6alpha_event=true, flag_d(x,y,z)6=ture"<<endl;
	  }
        }
        continue;
      }	
    }
  
    //exclude 1<n<7 or n>7 event
    //    else if (i%2==1 && ((n>1 && n<7) || n>7)) {
    else if (i%2==1 && ((n>1 && n<6) || n>6)) {
      eens=i;
      enps=n;
      ne_es++;
      #ifdef SAVE_FILE
      expscat->Fill();
      #endif
      continue;
    }
    //exclude the beam event
    else if (i%2==0) {
      eenb=i;
      enpb=n;
      ne_eb++;
      #ifdef SAVE_FILE
      expbeam->Fill();
      #endif
      continue;
    }
    else {
      eeno=i;
      enpo=n;
      ne_eo++;
      #ifdef SAVE_FILE
      expothe->Fill();
      #endif
      continue;      
    }
    // analysis of attpcsim.root to here

    // analysis of digi out from here
    //   get hit pattern
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
    
    //    if(i%2==1 && n==7) {
    if(i%2==1 && n==6) {
      ne_s++;
      #ifdef SAVE_FILE
      anadata->Fill();
      #endif
    }
    else if(i%2==0 && n<2) {
      #ifdef SAVE_FILE
      anadata->Fill();
      #endif
    }
    else{continue;}
  }

  cout<<"Finish 100 % !"<<endl<<endl;
  cout<<"all event "<<ne<<"   beam event "<<ne_b<<"   scattaring event "<<ne_s<<endl;
  cout<<"except event "<<ne_es+ne_eb+ne_eo<<endl;
  cout<<"except scattering event "<<ne_es<<endl;
  cout<<"except beam event "<<ne_eb<<endl;
  cout<<"except other problem event "<<ne_eo<<endl;

  //need
  #ifdef CANVAS
  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(4,3);
    
  //  h_z->SetStats(0);
  h_posi->SetStats(1);
  can->cd(1);
  //  gPad->SetLogy(1);
  h_posi->Draw("colz");

  h_dep->SetStats(1);
  can->cd(2);
  //  gPad->SetLogy(1);
  h_dep->Draw("colz");

  h_beam->SetStats(1);
  can->cd(3);
  //  gPad->SetLogy(1);
  h_beam->Draw("colz");

  h_h1->SetStats(1);
  can->cd(4);
  //  gPad->SetLogy(1);
  h_h1->Draw("colz");

  h_hoy->SetStats(1);
  can->cd(5);
  //  gPad->SetLogy(1);
  h_hoy->Draw("colz");

  h_ex->SetStats(1);
  can->cd(6);
  //  gPad->SetLogy(1);
  h_ex->Draw("colz");

  h_hang->SetStats(1);
  can->cd(7);
  //  gPad->SetLogy(1);
  h_hang->Draw("colz");

  h_ang->SetStats(1);
  can->cd(8);
  //  gPad->SetLogy(1);
  h_ang->Draw("colz");

  h_ex_dep->SetStats(1);
  can->cd(9);
  //  gPad->SetLogy(1);
  h_ex_dep->Draw("colz");

  h_ex_ang->SetStats(1);
  can->cd(10);
  //  gPad->SetLogy(1);
  h_ex_ang->Draw("colz");

  h_rmax->SetStats(1);
  can->cd(11);
  //  gPad->SetLogy(1);
  h_rmax->Draw("colz");

  h_rmaxnd->SetStats(1);
  can->cd(12);
  //  gPad->SetLogy(1);
  h_rmaxnd->Draw("colz");

  //output for canvas
  //  TString outcan = "data/check_hoyle_state_"+ to_string(subnum) +".C";
  //  TString outcan = "data1/can_check_hoyle_dire_"+ to_string(subnum) +".C";
  //  TString outcan = "data1/can_sim_macro_hoyle_"+ to_string(subnum) +".C";
  //  TString outcan = "data1/can_sim_macro_hoyle_nd_"+ to_string(subnum) +".C";
  //  TString outcan = "data/can_sim_macro_sequential_hoyle_nd_"+ to_string(subnum) +".C";
  TString outcan =
    outdir + "can_sim_macro_sequential_hoyle_90matm_"+ to_string(subnum) +".C";
  can->SaveAs(outcan);

  //need to here
  /*
  //  h_rmax->SetStats(0);
  h_dep->SetStats(1);
  //  gPad->SetLogy(1);
  h_dep->Draw("colz");
  */
  #endif
  #ifdef SAVE_FILE
  anadata->AutoSave();
  expscat->AutoSave();
  expbeam->AutoSave();
  expothe->AutoSave();
  #endif

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;

  #if defined(SAVE_FILE) && defined(CANVAS)
  cout<<"output: "<< outfname <<" and "<< outcan << endl;
  #elif ! defined(SAVE_FILE) && defined(CANVAS)
  cout<<"output: "<< outcan << endl;
  #elif defined(SAVE_FILE) && ! defined(CANVAS)
  cout<<"output: "<< outfname << endl;
  #endif
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}
