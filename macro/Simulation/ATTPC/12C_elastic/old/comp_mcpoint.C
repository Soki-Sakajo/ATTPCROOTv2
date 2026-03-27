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
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TMath.h"
#include <fstream>

int comp_mcpoint(Int_t subnum = 0)
{

  // set timer
  TStopwatch timer;
  timer.Start();

  // set constants
  Double_t c12_mass = 11177.736;
  Double_t h1_mass = 938.7669711;
  Double_t pi = 3.1415926535;
  Double_t rad_to_deg = 180 / pi;
  int GtoM2 = 1e6;

  //for MCTrack
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
  Int_t ne = 0;
  Int_t ne_b = 0;
  Int_t ne_s = 0;
  Int_t ne_es = 0;
  Int_t ne_eb = 0;
  Int_t ne_eo = 0;
  Int_t trackid = -100;
  Double_t mass[2];
  Double_t rmax;
  Double_t r_tmp;
  Double_t ener[2] = {0, 0};
  Double_t ener_max;
  Double_t tener;
  Double_t lab_theta[2];
  Double_t cms_theta[2];
  Double_t pba;
  Double_t startx[2];
  Double_t starty[2];
  Double_t startz[2];
  Double_t stopx[2];
  Double_t stopy[2];
  Double_t stopz[2];
  Double_t p_gx;
  Double_t p_gy;
  Double_t p_gz;
  Double_t p_ga;
  Double_t pxc[2];
  Double_t pyc[2];
  Double_t pzc[2];
  Double_t pac[2];
  Double_t psp[2];
  Double_t pbx;
  Double_t pby;
  Double_t pbz;
  Double_t px[2] = {0, 0};
  Double_t py[2] = {0, 0};
  Double_t pz[2] = {0, 0};
  Double_t pa[2] = {0, 0};
  Double_t sp[2] = {0, 0};
  Double_t r_nd[2] = {0, 0};
  Double_t rmax_nd = 0;
  Int_t pdg[2] = {0, 0};

  // for MCPoint
  Int_t eve_p;
  Int_t eenb_p;
  Int_t eens_p;
  Int_t eeno_p;
  Int_t n_p;
  Int_t np_p;
  Int_t enpb_p;
  Int_t enps_p;
  Int_t enpo_p;
  Int_t s_p;
  Int_t p_p;
  Int_t massnum_p[2] = {-100, -100};
  Int_t atomnum_p[2] = {-100, -100};
  Int_t massstop_p = -1000;
  Int_t atomstop_p = -1000;
  Int_t massstar_p = -1000;
  Int_t atomstar_p = -1000;
  Int_t masstemp_p[2] = {-100, -100};
  Int_t atomtemp_p[2] = {-100, -100};
  Int_t massbefo_p[2] = {-100, -100};
  Int_t atombefo_p[2] = {-100, -100};
  Int_t ne_p = 0;
  Int_t ne_b_p = 0;
  Int_t ne_s_p = 0;
  Int_t ne_es_p = 0;
  Int_t ne_eb_p = 0;
  Int_t ne_eo_p = 0;
  Int_t trackid_p = -100;
  Double_t mass_p[2];
  Double_t deltaz_p;
  Double_t rmax_p;
  Double_t r_tmp_p;
  Double_t ener_p[2] = {0, 0};
  Double_t ener_max_p;
  Double_t tener_p;
  Double_t lab_theta_p[2];
  Double_t cms_theta_p[2];
  Double_t bx_p;
  Double_t by_p;
  Double_t bz_p;
  Double_t pba_p;
  Double_t startx_p[2];
  Double_t starty_p[2];
  Double_t startz_p[2];
  Double_t stopx_p[2];
  Double_t stopy_p[2];
  Double_t stopz_p[2];
  Double_t p_gx_p;
  Double_t p_gy_p;
  Double_t p_gz_p;
  Double_t p_ga_p;
  Double_t deltaE_p;
  Double_t been_p;
  Double_t momd_p;
  Double_t pxc_p[2];
  Double_t pyc_p[2];
  Double_t pzc_p[2];
  Double_t pac_p[2];
  Double_t psp_p[2];
  Double_t pbx_p;
  Double_t pby_p;
  Double_t pbz_p;
  Double_t px_p[2] = {0, 0};
  Double_t py_p[2] = {0, 0};
  Double_t pz_p[2] = {0, 0};
  Double_t pa_p[2] = {0, 0};
  Double_t sp_p[2] = {0, 0};
  Double_t r_nd_p[2] = {0, 0};
  Double_t rmax_nd_p = 0;
  bool beam_start_p = true;
  bool scat_start_p[2] = {true, true};

  // input
  TString infname = "data1/attpcsim_pp_" + to_string(subnum) + ".root";
  TFile *file = new TFile(infname);
  TTree *tree = (TTree *)file->Get("cbmsim");

  // digi data input
  TChain tpc_tree("cbmsim");
  TString digifname = "data1/output_digi_rcnp_p_nd_" + to_string(subnum) + ".root";
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
  fAtMapPtr->GeneratePadPlane();             // Adam suggest
  auto fPadPlane = fAtMapPtr->GetPadPlane(); // Adam suggest

  // set reading data
  AtMCTrack *track = new AtMCTrack();
  TClonesArray *trackArray = 0;
  tree = (TTree *)file->Get("cbmsim");
  tree->SetBranchAddress("MCTrack", &trackArray);

  // get AtMCPoint data for MCPoint
  AtMCPoint *point_p = new AtMCPoint();
  TClonesArray *pointArray_p = 0;
  tree->SetBranchAddress("AtTpcPoint", &pointArray_p);
  Int_t nEvents_p = tree->GetEntriesFast();
  Int_t n_d_p[nEvents_p];
  cout <<"set AtMCPoint"<< endl;

  // get FairMCPoint data for MCTrack
  FairMCPoint *point = new FairMCPoint();
  TClonesArray *pointArray = 0;
  tree->SetBranchAddress("AtTpcPoint", &pointArray);

  Int_t nEvents = tree->GetEntriesFast();
  Int_t n_d[nEvents];

  //  cout << nEvents <<endl;
  cout <<"set FairMCPoint"<< endl;

  // set output
  TString outfname = "data1/comp_mcpoint_" + to_string(subnum) + ".root";
  TFile *simmacro = new TFile(outfname, "recreate");

  //  if(subnum==0){TFile *simmacro = new TFile(outfname, "recreate");}
  //  else {
  if(subnum!=0){
    cout<<"need to write output file when subnum!=0"<<endl;
    //TFile *simmacro = new TFile(outfname, "update");
    exit(0);
  }
  //for MCTrack
  TTree *trackdata = new TTree("trackdata", "trackdata");
  TTree *trackexbeam = new TTree("trackexbeam", "trackexbeam");
  TTree *trackexscat = new TTree("trackexscat", "trackexscat");
  TTree *trackexoth = new TTree("trackexoth", "trackexoth");
  trackdata->Branch("event_number", &eve, "eve/I");
  trackdata->Branch("PdgCode", &pdg, "pdg[2]/I");
  trackdata->Branch("startx", &startx, "startx[2]/D");
  trackdata->Branch("starty", &starty, "starty[2]/D");
  trackdata->Branch("startz", &startz, "startz[2]/D");
  trackdata->Branch("stopx", &stopx, "stopx[2]/D");
  trackdata->Branch("stopy", &stopy, "stopy[2]/D");
  trackdata->Branch("stopz", &stopz, "stopz[2]/D");
  trackdata->Branch("Energy", &ener, "ener[2]/D");
  trackdata->Branch("ener_max", &ener_max, "ener_max/D");
  trackdata->Branch("TotalE", &tener, "TotalE/D");
  trackdata->Branch("theta_lab", &lab_theta, "theta_lab[2]/D");
  trackdata->Branch("theta_cms", &cms_theta, "theta_cms[2]/D");
  trackdata->Branch("n", &n, "n/I");
  trackdata->Branch("rmax_nd", &rmax_nd, "rmax_nd/D");
  trackdata->Branch("rmax", &rmax, "rmax/D");
  trackexbeam->Branch("event_number_beam", &eenb, "eenb/I");
  trackexbeam->Branch("number_of_particle_beam", &enpb, "enpb/I");
  trackexscat->Branch("event_number_scat", &eens, "eens/I");
  trackexscat->Branch("number_of_particle_scat", &enps, "enps/I");
  trackexoth->Branch("event_number_other", &eeno, "eeno/I");
  trackexoth->Branch("number_of_particle_other", &enpo, "enpo/I");

  //for MCPoint
  TTree *pointdata = new TTree("pointdata", "pointdata");
  TTree *pointexbeam = new TTree("pointexbeam", "pointexbeam");
  TTree *pointexscat = new TTree("pointexscat", "pointexscat");
  TTree *pointexoth = new TTree("pointexoth", "pointexoth");
  pointdata->Branch("event_number", &eve_p, "eve/I");
  pointdata->Branch("startx", &startx_p, "startx[2]/D");
  pointdata->Branch("starty", &starty_p, "starty[2]/D");
  pointdata->Branch("startz", &startz_p, "startz[2]/D");
  pointdata->Branch("stopx", &stopx_p, "stopx[2]/D");
  pointdata->Branch("stopy", &stopy_p, "stopy[2]/D");
  pointdata->Branch("stopz", &stopz_p, "stopz[2]/D");
  pointdata->Branch("massnum", &massnum_p, "massnum[2]/I");
  pointdata->Branch("atomnum", &atomnum_p, "atomnum[2]/I");
  pointdata->Branch("Energy", &ener_p, "ener[2]/D");
  pointdata->Branch("ener_max", &ener_max_p, "ener_max/D");
  pointdata->Branch("TotalE", &tener_p, "TotalE/D");
  pointdata->Branch("theta_lab", &lab_theta_p, "theta_lab[2]/D");
  pointdata->Branch("theta_cms", &cms_theta_p, "theta_cms[2]/D");
  pointdata->Branch("n", &n_p, "n/I");
  pointdata->Branch("rmax_nd", &rmax_nd_p, "rmax_nd/D");
  pointdata->Branch("rmax", &rmax_p, "rmax/D");
  pointexbeam->Branch("event_number_beam", &eenb_p, "eenb/I");
  pointexbeam->Branch("number_of_particle_beam", &enpb_p, "enpb/I");
  pointexscat->Branch("event_number_scat", &eens_p, "eens/I");
  pointexscat->Branch("number_of_particle_scat", &enps_p, "enps/I");
  pointexoth->Branch("event_number_other", &eeno_p, "eeno/I");
  pointexoth->Branch("number_of_particle_other", &enpo_p, "enpo/I");

  /*
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
  cout<<endl;
  */

  // analysis of attpcsim.root for MCTrack from here
  for (int i = 0; i < nEvents; i++)
  {
    tree->GetEvent(i);
    p = trackArray->GetEntries();
    n_d[i] = p;
    if (i % 100 == 0)
    {
      // need
      cout << "Get data for MCTrack " << 100 * i / nEvents << " % ! \r" << flush;
    }
  }
  cout << "Get data for MCTrack 100 % ! \r" << flush;
  for (int i = 0; i < nEvents; i++)
  {
    //  for(int i=0;i < 100; i++){
    if (i % 1000 == 0)
    {
      // need
      cout << "Finish for MCTrack " << 100 * i / nEvents << " % !     \r" << flush;
    }
    eve = i;
    ne++;
    // initialization
    pdg[0] = -10;
    pdg[1] = -10;
    for (int pernum = 0; pernum < 2; pernum++)
    {
      startx[pernum] = -1000;
      starty[pernum] = -1000;
      startz[pernum] = -1000;
      stopx[pernum] = -1000;
      stopy[pernum] = -1000;
      stopz[pernum] = -1000;
      r_nd[pernum] = -1000;
    }
    lab_theta[0] = -1000;
    cms_theta[0] = -1000;
    p_gx = -1000;
    p_gy = -1000;
    p_gz = -1000;
    pxc[0] = -1000;
    pyc[0] = -1000;
    pzc[0] = -1000;
    pxc[1] = -1000;
    pyc[1] = -1000;
    pzc[1] = -1000;
    rmax_nd = -1000;

    // get event
    tree->GetEvent(i);
    n = trackArray->GetEntries();
    np = pointArray->GetEntries();
    //    cout<<i <<" "<<n<<endl;

    // calculate energy, ang
    // set for beam event
    if (i % 2 == 0 && n == 1)
    {
      track = (AtMCTrack *)trackArray->At(0);
      pdg[0] = track->GetPdgCode();
      pdg[1] = 0;
      if (pdg[0] == 1000060120)
      {
        ne_b++;
        startx[0] = track->GetStartX();
        starty[0] = track->GetStartY();
        startz[0] = track->GetStartZ();
        startx[1] = 0;
        starty[1] = 0;
        startz[1] = 0;
        px[0] = track->GetPx();
        py[0] = track->GetPy();
        pz[0] = track->GetPz();
        px[1] = 0;
        py[1] = 0;
        pz[1] = 0;
        ener[0] = (px[0] * px[0] + py[0] * py[0] + pz[0] * pz[0]) * GtoM2 / c12_mass / 2;
        ener[1] = 0;
        ener_max = 0;
        tener = 0;
        pxc[0] = 0;
        pyc[0] = 0;
        pzc[0] = 0;
        pxc[1] = 0;
        pyc[1] = 0;
        pzc[1] = 0;
        lab_theta[0] = 0;
        lab_theta[1] = 0;
        cms_theta[0] = 0;
        cms_theta[1] = 0;
        point = (FairMCPoint *)pointArray->At(np - 1);
        stopx[0] = point->GetX();
        stopy[0] = point->GetY();
        stopz[0] = point->GetZ();
        stopx[1] = 0;
        stopy[1] = 0;
        stopz[1] = 0;
        r_nd[0] = sqrt(stopx[0] * stopx[0] + stopy[0] * stopy[0]);
        r_nd[1] = 0;
        rmax_nd = r_nd[0];
        pbx = point->GetPx();
        pby = point->GetPy();
        pbz = point->GetPz();
        pba = sqrt(pbx * pbx + pby * pby + pbz * pbz);
      }
      else
      {
        cout << i << " something is wrong with pdg at beam event !!!     " << endl
             << endl;
      }
    }

    // exclude beam and strange event(n>2)
    else if (i % 2 == 1 && n == 2 && n_d[i - 1] == 1)
    {
      for (int j = 0; j < n; j++)
      {
        track = (AtMCTrack *)trackArray->At(j);
        pdg[j] = track->GetPdgCode();
        startx[j] = track->GetStartX();
        starty[j] = track->GetStartY();
        startz[j] = track->GetStartZ();
        px[j] = track->GetPx();
        py[j] = track->GetPy();
        pz[j] = track->GetPz();
        for (int k = 0; k < np; k++)
        {
          point = (FairMCPoint *)pointArray->At(k);
          trackid = (int)(point->GetTrackID());
          stopx[j] = point->GetX();
          stopy[j] = point->GetY();
          stopz[j] = point->GetZ();
          r_nd[j] = 10 * sqrt(stopx[j] * stopx[j] + stopy[j] * stopy[j]);
        }
        if (pdg[j] == 1000060120)
        {
          mass[j] = c12_mass;
          ener[j] = (px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]) * GtoM2 / mass[j] / 2;
          pa[j] = sqrt(px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]); // abs vector
          sp[j] = pbx * px[j] + pby * py[j] + pbz * pz[j];             // scalar product
          lab_theta[j] = rad_to_deg * acos(sp[j] / (pba * pa[j]));
        }
        else if (pdg[j] == 2212)
        {
          mass[j] = h1_mass;
          ener[j] = (px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]) * GtoM2 / h1_mass / 2;
          pa[j] = sqrt(px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]) / mass[j];  // abs vector
          sp[j] = (pbx * px[j] + pby * py[j] + pbz * pz[j]) / c12_mass / mass[j]; // scalar product
          lab_theta[j] = rad_to_deg * acos(sp[j] / (pba * pa[j]));
        }
        else
        {
          cout << i << " " << j << " something is wrong with pdg !!!     " << endl
               << endl;
        }
      }
      if (r_nd[0] >= r_nd[1])
      {
        rmax_nd = r_nd[0];
      }
      else if (r_nd[0] < r_nd[1])
      {
        rmax_nd = r_nd[1];
      }
      // total energy
      tener = (ener[0] + ener[1]);
      if (ener[0] > ener[1])
        ener_max = ener[0];
      if (ener[0] <= ener[1])
        ener_max = ener[1];
      p_gx = c12_mass * pbx / (mass[0] + mass[1]);
      p_gy = c12_mass * pby / (mass[0] + mass[1]);
      p_gz = c12_mass * pbz / (mass[0] + mass[1]);
      p_ga = sqrt(p_gx * p_gx + p_gy * p_gy + p_gz * p_gz); // abs vector in cms

      // about angle in cms
      for (int k = 0; k < n; k++)
      {
        pxc[k] = px[k] - p_gx;
        pyc[k] = py[k] - p_gy;
        pzc[k] = pz[k] - p_gz;
        pac[k] = sqrt(pxc[k] * pxc[k] + pyc[k] * pyc[k] + pzc[k] * pzc[k]); // abs vector in cms
        psp[k] = p_gx * pxc[k] + p_gy * pyc[k] + p_gz * pzc[k];             // scalar product
        cms_theta[k] = rad_to_deg * acos(psp[k] / (p_ga * pac[k]));
      }
    }

    // display n>3 event
    // exclude the beam event
    else if (i % 2 == 1 && n > 2)
    {
      eens = i;
      enps = n;
      ne_es++;
      trackexscat->Fill();
      continue;
    }
    else if (i % 2 == 0)
    {
      eenb = i;
      enpb = n;
      ne_eb++;
      trackexbeam->Fill();
      continue;
    }
    else
    {
      eeno = i;
      enpo = n;
      ne_eo++;
      trackexoth->Fill();
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

    double rmax_tmp = -100;

    for (int m = 0; m < numHits; m++)
    {
      // Get the padNum for the hit
      auto hit = eventPtr->GetHit(m);
      auto position = hit.GetPosition();
      auto charge = hit.GetCharge();

      r_tmp = sqrt(position.X() * position.X() + position.Y() * position.Y());
      if (r_tmp > rmax_tmp)
        rmax_tmp = r_tmp;
    }
    rmax = rmax_tmp;

    if (i % 2 == 1 && n < 3)
    {
      ne_s++;
      trackdata->Fill();
    }
    else if (i % 2 == 0 && n < 2)
    {
      trackdata->Fill();
    }
    else
    {
      continue;
    }
  }
  cout << "Finish for MCTrack 100 % !" << endl;

  //for MCTrack to here

  // analysis of attpcsim.root for MCPoint from here
  for (int i = 0; i < nEvents_p; i++)
  {
    tree->GetEvent(i);
    p_p = trackArray->GetEntries();
    n_d_p[i] = p_p;
    if (i % 100 == 0)
    {
      // need
      cout << "Get data for MCPoint " << 100 * i / nEvents_p << " % ! \r" << flush;
    }
  }
  cout << "Get data for MCPoint 100 % ! \r" << flush;
  for (int i = 0; i < nEvents_p; i++)
  {
    //  for(int i=0;i < 100; i++){
    if (i % 1000 == 0)
    {
      // need
      cout << "Finish for MCPoint " << 100 * i / nEvents_p << " % !     \r" << flush;
    }
    eve_p = i;
    ne_p++;
    // initialization
    for (int pernum = 0; pernum < 2; pernum++)
    {
      startx_p[pernum] = -1000;
      starty_p[pernum] = -1000;
      startz_p[pernum] = -1000;
      stopx_p[pernum] = -1000;
      stopy_p[pernum] = -1000;
      stopz_p[pernum] = -1000;
      r_nd_p[pernum] = -1000;
    }
    lab_theta_p[0] = -1000;
    cms_theta_p[0] = -1000;
    p_gx_p = -1000;
    p_gy_p = -1000;
    p_gz_p = -1000;
    pxc_p[0] = -1000;
    pyc_p[0] = -1000;
    pzc_p[0] = -1000;
    pxc_p[1] = -1000;
    pyc_p[1] = -1000;
    pzc_p[1] = -1000;
    rmax_nd_p = -1000;

    // get event
    tree->GetEvent(i);
    n_p = trackArray->GetEntries();
    np_p = pointArray_p->GetEntries();
    //    cout<<i <<" "<<n<<endl;

    // calculate energy, ang
    // set for beam event
    if (i % 2 == 0 && n_p == 1)
    {
      ne_b_p++;
      point_p = (AtMCPoint *)pointArray_p->At(0); // get the point information
      trackid_p = (int)(point_p->GetTrackID());   // get track ID
      if (trackid_p == 0)
      {
        massstar_p = point_p->GetMassNum();
        atomstar_p = point_p->GetAtomicNum();
        // cout<<i<<"  mass number "<<massstar<<endl;
        if (massstar_p == 12 && atomstar_p == 6)
        {
          startx_p[0] = point_p->GetX();
          starty_p[0] = point_p->GetY();
          startz_p[0] = point_p->GetZ();
          startx_p[1] = 0;
          starty_p[1] = 0;
          startz_p[1] = 0;
          px_p[0] = point_p->GetPx();
          py_p[0] = point_p->GetPy();
          pz_p[0] = point_p->GetPz();
          px_p[1] = 0;
          py_p[1] = 0;
          pz_p[1] = 0;
          ener_p[0] = (px_p[0] * px_p[0] + py_p[0] * py_p[0] + pz_p[0] * pz_p[0]) * GtoM2 / c12_mass / 2;
          ener_p[1] = 0;
          ener_max_p = 0;
          tener_p = 0;
          pxc_p[0] = 0;
          pyc_p[0] = 0;
          pzc_p[0] = 0;
          pxc_p[1] = 0;
          pyc_p[1] = 0;
          pzc_p[1] = 0;
          lab_theta_p[0] = 0;
          lab_theta_p[1] = 0;
          cms_theta_p[0] = 0;
          cms_theta_p[1] = 0;
        }
      }
      else
      {
        cout << i << " something is wrong with track ID at start point in beam event !!!     " << endl
             << endl;
      }

      for (int k = 0; k < np_p; k++)
      {
        point_p = (AtMCPoint *)pointArray_p->At(k); // get the point information
        trackid_p = (int)(point_p->GetTrackID());   // get track ID
        massstop_p = point_p->GetMassNum();
        atomstop_p = point_p->GetAtomicNum();
        if (trackid_p == 0 && massstop_p == 12 && atomstop_p == 6)
        {
          /*
            if(beam_start){
            massbefo[j]=masstemp[j];
            atombefo[j]=atomtemp[j];
            }
          */
          //	cout<<i<<"  mass number "<<massstop<<endl;
          stopx_p[0] = point_p->GetX();
          stopy_p[0] = point_p->GetY();
          stopz_p[0] = point_p->GetZ();
          stopx_p[1] = 0;
          stopy_p[1] = 0;
          stopz_p[1] = 0;
          r_nd_p[0] = sqrt(stopx_p[0] * stopx_p[0] + stopy_p[0] * stopy_p[0]);
          r_nd_p[1] = 0;
          rmax_nd_p = r_nd_p[0];
          bx_p = stopx_p[0];
          by_p = stopy_p[0];
          bz_p = stopz_p[0];
          pbx_p = point_p->GetPx();
          pby_p = point_p->GetPy();
          pbz_p = point_p->GetPz();
          pba_p = sqrt(pbx_p * pbx_p + pby_p * pby_p + pbz_p * pbz_p);
          if (massstar_p == massstop_p && atomstar_p == atomstop_p)
          {
            massnum_p[0] = massstop_p;
            massnum_p[1] = 0;
            atomnum_p[0] = atomstop_p;
            atomnum_p[1] = 0;
          }
          else
          {
            cout << i << "  " << k << " something is wrong with mass and atomic number at start and stop in beam event !!!" << endl
                 << endl;
          }
        }
        else
        {
          cout << i << " something is wrong with track ID at end point in beam event !!!     " << endl
               << endl;
        }
      }
    }

    // exclude beam and strange event(n>2)
    else if (i % 2 == 1 && n_p == 2 && n_d_p[i - 1] == 1)
    {
      for (int j = 0; j < n_p; j++)
      {
        scat_start_p[0] = true;
        scat_start_p[1] = true;
        for (int k = 0; k < np_p; k++)
        {
          point_p = (AtMCPoint *)pointArray_p->At(k);
          trackid_p = (int)(point_p->GetTrackID());
          if (trackid_p == j)
          {
            masstemp_p[j] = point_p->GetMassNum();
            atomtemp_p[j] = point_p->GetAtomicNum();
            if (scat_start_p[j])
            {
              massbefo_p[j] = masstemp_p[j];
              atombefo_p[j] = atomtemp_p[j];
              startx_p[j] = point_p->GetX();
              starty_p[j] = point_p->GetY();
              startz_p[j] = point_p->GetZ();
              px_p[j] = point_p->GetPx();
              py_p[j] = point_p->GetPy();
              pz_p[j] = point_p->GetPz();
              scat_start_p[j] = false;
            }
            stopx_p[j] = point_p->GetX();
            stopy_p[j] = point_p->GetY();
            stopz_p[j] = point_p->GetZ();
            r_nd_p[j] = 10 * sqrt(stopx_p[j] * stopx_p[j] + stopy_p[j] * stopy_p[j]);
            if (masstemp_p[j] == massbefo_p[j] && atomtemp_p[j] == atombefo_p[j])
            {
              massnum_p[j] = masstemp_p[j];
              atomnum_p[j] = atomtemp_p[j];
            }
            else
            {
              cout << i << "  " << j << " something wrong with mass and atomic number in scattering event !!!   " << endl
                   << endl;
            }
            massbefo_p[j] = masstemp_p[j];
            atombefo_p[j] = atomtemp_p[j];
          }
        }
        if (massnum_p[j] == 12 && atomnum_p[j] == 6)
        {
          mass_p[j] = c12_mass;
          ener_p[j] = (px_p[j] * px_p[j] + py_p[j] * py_p[j] + pz_p[j] * pz_p[j]) * GtoM2 / mass_p[j] / 2;
          pa_p[j] = sqrt(px_p[j] * px_p[j] + py_p[j] * py_p[j] + pz_p[j] * pz_p[j]); // abs vector
          sp_p[j] = pbx_p * px_p[j] + pby_p * py_p[j] + pbz_p * pz_p[j];             // scalar product
          lab_theta_p[j] = rad_to_deg * acos(sp_p[j] / (pba_p * pa_p[j]));
        }
        else if (massnum_p[j] == 1 && atomnum_p[j] == 1)
        {
          mass_p[j] = h1_mass;
          ener_p[j] = (px_p[j] * px_p[j] + py_p[j] * py_p[j] + pz_p[j] * pz_p[j]) * GtoM2 / h1_mass / 2;
          pa_p[j] = sqrt(px_p[j] * px_p[j] + py_p[j] * py_p[j] + pz_p[j] * pz_p[j]); // abs vector
          sp_p[j] = pbx_p * px_p[j] + pby_p * py_p[j] + pbz_p * pz_p[j];             // scalar product
          lab_theta_p[j] = rad_to_deg * acos(sp_p[j] / (pba_p * pa_p[j]));
        }
        else
        {
          cout << i << " " << j << " something is wrong with mass and atomic number !!!     " << endl
               << endl;
        }
      }
      if (r_nd_p[0] >= r_nd_p[1])
      {
        rmax_nd_p = r_nd_p[0];
      }
      else if (r_nd_p[0] < r_nd_p[1])
      {
        rmax_nd_p = r_nd_p[1];
      }
      // total energy
      tener_p = (ener_p[0] + ener_p[1]);
      if (ener_p[0] > ener_p[1])
        ener_max_p = ener_p[0];
      if (ener_p[0] <= ener_p[1])
        ener_max_p = ener_p[1];
      p_gx_p = c12_mass * pbx_p / (mass_p[0] + mass_p[1]);
      p_gy_p = c12_mass * pby_p / (mass_p[0] + mass_p[1]);
      p_gz_p = c12_mass * pbz_p / (mass_p[0] + mass_p[1]);
      p_ga_p = sqrt(p_gx_p * p_gx_p + p_gy_p * p_gy_p + p_gz_p * p_gz_p); // abs vector in cms

      // about angle in cms
      for (int k = 0; k < n_p; k++)
      {
        pxc_p[k] = px_p[k] - p_gx_p;
        pyc_p[k] = py_p[k] - p_gy_p;
        pzc_p[k] = pz_p[k] - p_gz_p;
        pac_p[k] = sqrt(pxc_p[k] * pxc_p[k] + pyc_p[k] * pyc_p[k] + pzc_p[k] * pzc_p[k]); // abs vector in cms
        psp_p[k] = p_gx_p * pxc_p[k] + p_gy_p * pyc_p[k] + p_gz_p * pzc_p[k];             // scalar product
        cms_theta_p[k] = rad_to_deg * acos(psp_p[k] / (p_ga_p * pac_p[k]));
      }
    }

    // display n>3 event
    // exclude the beam event
    else if (i % 2 == 1 && n_p > 2)
    {
      eens_p = i;
      enps_p = n_p;
      ne_es_p++;
      pointexscat->Fill();
      continue;
    }
    else if (i % 2 == 0)
    {
      eenb_p = i;
      enpb_p = n_p;
      ne_eb_p++;
      pointexbeam->Fill();
      continue;
    }
    else
    {
      eeno_p = i;
      enpo_p = n_p;
      ne_eo_p++;
      pointexoth->Fill();
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

    double rmax_tmp_p = -100;

    for (int m = 0; m < numHits; m++)
    {
      // Get the padNum for the hit
      auto hit = eventPtr->GetHit(m);
      auto position = hit.GetPosition();
      auto charge = hit.GetCharge();

      r_tmp_p = sqrt(position.X() * position.X() + position.Y() * position.Y());
      if (r_tmp_p > rmax_tmp_p)
        rmax_tmp_p = r_tmp_p;
    }
    rmax_p = rmax_tmp_p;

    if (i % 2 == 1 && n_p < 3)
    {
      ne_s_p++;
      pointdata->Fill();
    }
    else if (i % 2 == 0 && n_p < 2)
    {
      pointdata->Fill();
    }
    else
    {
      continue;
    }
  }
  cout << "Finish for MCPoint 100 % !" << endl << endl;

  //for MCPoint to here

  //save data
  trackdata->AutoSave();
  trackexscat->AutoSave();
  trackexbeam->AutoSave();
  trackexoth->AutoSave();
  pointdata->AutoSave();
  pointexscat->AutoSave();
  pointexbeam->AutoSave();
  pointexoth->AutoSave();
  simmacro->Close();

  //output about exepted data
  cout << "For MCTrack data "<< endl;
  cout << "all event " << ne << "   beam event " << ne_b << "   scattaring event " << ne_s << endl;
  cout << "except event " << ne_es + ne_eb + ne_eo << endl;
  cout << "except scattering event " << ne_es << endl;
  cout << "except beam event " << ne_eb << endl;
  cout << "except other problem event " << ne_eo << endl<< endl;

  cout << "For MCPoint data "<< endl;
  cout << "all event " << ne_p << "   beam event " << ne_b_p << "   scattaring event " << ne_s_p << endl;
  cout << "except event " << ne_es_p + ne_eb_p + ne_eo_p << endl;
  cout << "except scattering event " << ne_es_p << endl;
  cout << "except beam event " << ne_eb_p << endl;
  cout << "except other problem event " << ne_eo_p << endl;

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl
       << endl;

  return 0;
}
