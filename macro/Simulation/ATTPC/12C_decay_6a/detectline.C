//#define bins_set_ver

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
#include <fstream>
#include <vector>
#include"TString.h"
#include"TFile.h"
#include"TTree.h"
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TClonesArray.h"
#include"TCanvas.h"
#include"TMath.h"
#include"TChain.h"
#include"THnSparse.h"
#include<algorithm>

int SetParameter(Int_t& dim, Double_t* start_val, Double_t* stop_val, Double_t* step_val, Int_t* bins_set,Int_t& n_get){
  dim=4;
  n_get=15; // getting combinations of parameter of line
  //psi
  start_val[0]=0;
  stop_val[0] =360;
  step_val[0] =360;
  bins_set[0]=-10;
  //theta
  start_val[1]=0;
  stop_val[1] =90;
  step_val[1] =90;
  bins_set[1]=-10;
  //phi
  start_val[2]=0;
  stop_val[2] =360;
  step_val[2] =360;
  bins_set[2]=-10;
  //rho
  start_val[3]=-1076;
  stop_val[3] =1076;
  step_val[3] =1076;
  bins_set[3]=-10;

  return 0;
};

int GetHit(Int_t& eventNumber,Int_t& hitnum, vector<vector<Double_t>>& hitpoint){

  double fThreshold = 0; // Only save hits if they're above this threshold
  
  TChain tpc_tree("cbmsim");
  tpc_tree.Add("./data1/output_digi_rcnp_sequential_hoylehoyle_nd_99.root");
  
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
  
  std::vector<std::vector<Short_t>> traces; // traces[event][time]
  
  // Get the event to read
  reader.SetEntry(eventNumber);
  // Get the event
  AtEvent *eventPtr = (AtEvent *)(event->At(0));
  auto numHits = eventPtr->GetNumHits();
  
  hitnum=numHits;
  hitpoint.resize(hitnum, vector<Double_t>(5,0));

  //if multi event, need to modify about vector "hitpoint"
  for (int i = 0; i < numHits; ++i) {
    // Get the padNum for the hit
    auto hit = eventPtr->GetHit(i);
    auto position = hit.GetPosition();
    auto charge = hit.GetCharge();

    hitpoint.at(i).at(0)=eventNumber;
    hitpoint.at(i).at(1)=position.X();
    hitpoint.at(i).at(2)=position.Y();
    hitpoint.at(i).at(3)=position.Z();
    hitpoint.at(i).at(4)=charge;
    
    //    traceFile << position.X() << " " << position.Y() << " " << position.Z() << " " << charge << std::endl;
  }
  return 0;
};

int ShiftPoint(vector<vector<Double_t>>& hitpoint, Double_t& max_z){

  Double_t temp_z, shift_z;
  
  for(int i=0; i<hitpoint.size(); i++){
    if(max_z<hitpoint.at(i).at(3)){
      max_z=hitpoint.at(i).at(3);
    }
  }
  for(int i=0; i<hitpoint.size(); i++){
    temp_z=hitpoint.at(i).at(3);
    shift_z=-temp_z+max_z;
    hitpoint.at(i).at(3)=shift_z;
  }

  return 0;  
};

Double_t rho_fun(const Double_t* coord){

  //set rho function
  Double_t a=coord[0];
  Double_t b=coord[1];
  Double_t c=coord[2];
  Double_t x=coord[3];
  Double_t y=coord[4];
  Double_t z=coord[5];

  Double_t A=cos(a)*cos(c)-sin(a)*cos(b)*sin(c);
  Double_t B=sin(a)*cos(c)+cos(a)*cos(b)*sin(c);
  Double_t C=sin(b)*sin(c);

  return x*A+y*B+z*C;
};

int CalculateHoughSpace(Int_t* bins, Int_t bins_dim, Double_t* start_val, Double_t* stop_val, Double_t* step_val, Int_t* bins_set, Bool_t check_bins_set_ver, Double_t* end_val){

  if(check_bins_set_ver){
    // need to modify
    Double_t deff_val[4];
    for(Int_t i=0; i<bins_dim; i++){
      deff_val[i]=stop_val[i]-start_val[i];
      bins[i]=bins_set[i];
      end_val[i]=start_val[i]+bins[i]*step_val[i];
    }
    cout<<endl;
    cout<<"information of hough space"<<endl;
    cout<<"\"start\" is the start point of calculation."<<endl;
    cout<<"\"stop\" is the goal point of calculation."<<endl;
    cout<<"\"bin\" is the number of calculation points."<<endl<<endl;
    cout<<"start (psi,theta,phi,x,y,z)= ("<<start_val[0]<<", "<<start_val[1]<<", "
	<<start_val[2]<<", "<<start_val[3]<<")"<<endl;
    cout<<"stop  (psi,theta,phi,x,y,z)= ("<<stop_val[0]<<", "<<stop_val[1]<<", "
	<<stop_val[2]<<", "<<stop_val[3]<<")"<<endl;
    cout<<"step  (psi,theta,phi,x,y,z)= ("<<step_val[0]<<", "<<step_val[1]<<", "
	<<step_val[2]<<", "<<step_val[3]<<")"<<endl;
    cout<<"bin   (psi,theta,phi,rho)=("<<bins[0]<<", "<<bins[1]<<", "<<bins[2]<<", "<<bins[3]<<") "<<endl;    
  }
  
  else{
    Double_t deff_val[4];
    Double_t end_val[4];
    for(Int_t i=0; i<bins_dim; i++){
      deff_val[i]=stop_val[i]-start_val[i];
      bins[i]=(int)(deff_val[i]/step_val[i]+0.9);
      end_val[i]=start_val[i]+bins[i]*step_val[i];
    }
    cout<<endl;
    cout<<"information of hough space"<<endl;
    cout<<"\"start\" is the start point of calculation."<<endl;
    cout<<"\"stop\" is the goal point of calculation."<<endl;
    cout<<"\"bin\" is the number of calculation points."<<endl<<endl;
    cout<<"start (psi,theta,phi,x,y,z)= ("<<start_val[0]<<", "<<start_val[1]<<", "
	<<start_val[2]<<", "<<start_val[3]<<")"<<endl;
    cout<<"stop  (psi,theta,phi,x,y,z)= ("<<stop_val[0]<<", "<<stop_val[1]<<", "
	<<stop_val[2]<<", "<<stop_val[3]<<")"<<endl;
    cout<<"step  (psi,theta,phi,x,y,z)= ("<<step_val[0]<<", "<<step_val[1]<<", "
	<<step_val[2]<<", "<<step_val[3]<<")"<<endl;
    cout<<"bin   (psi,theta,phi,rho)=("<<bins[0]<<", "<<bins[1]<<", "<<bins[2]<<", "<<bins[3]<<") "<<endl;
  }
  return 0;
};

int GenerateHoughspace(Int_t dim, Int_t* bins, Double_t* start_val, Double_t* end_val, Double_t* step_val,THnSparseD* houghspace, vector<vector<Double_t>> hitpoint, Int_t hitnum){

  const Double_t pi=3.141592653589;
  const Double_t deg_to_rad=pi/180;
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t hough_val[4];
  Double_t psi;
  Double_t theta;
  Double_t phi;
  Double_t psi_r;
  Double_t theta_r;
  Double_t phi_r;
  Double_t rho;
  Double_t rho_coord[6];
  
  for(Int_t i=0; i<hitnum; i++){
    x=hitpoint.at(i).at(1);
    y=hitpoint.at(i).at(2);
    z=hitpoint.at(i).at(3);
    rho_coord[3]=x;
    rho_coord[4]=y;
    rho_coord[5]=z;
    for(Int_t j=0; j<bins[0]; j++){
      psi=start_val[0]+j*step_val[0];
      psi_r=psi*deg_to_rad;
      hough_val[0]=psi;
      rho_coord[0]=psi_r;
      for(Int_t k=0; k<bins[1]; k++){
	theta=start_val[1]+k*step_val[1];
	theta_r=theta*deg_to_rad;
	hough_val[1]=theta;
	rho_coord[1]=theta_r;
	for(Int_t l=0; l<bins[2]; l++){
	  phi=start_val[2]+l*step_val[2];
	  phi_r=phi*deg_to_rad;
	  hough_val[2]=phi;
	  rho_coord[2]=phi_r;

	  rho=rho_fun(rho_coord);
	  hough_val[3]=rho;
	  houghspace->Fill(hough_val);
	}
      }
    }
  }

  //  TString outcan="data1/hough_space.C";
  //  houghspace->SaveAs(outcan);

  return 0;
};

int GetBinInfo(THnSparseD* houghspace,Int_t* bins, Long64_t* n_bin, vector<vector<Double_t>>& binlist, vector<vector<Double_t>>& binlist_raw){

  Int_t frequency;
  Int_t check_bins[4];
  Double_t coord[4];
  Bool_t check_number_bins=true;
  binlist.resize(6,vector<Double_t>(0));
  
  n_bin=houghspace->GetNbins();
  houghspace->GetBins(check_bins);
  //need to check number of bin

  for(Int_t i=0; i<4; i++){
    if(bins[i]!=check_bins[i]){
      check_number_bins=false;
    }
  }
  if(check_number_bins==false){
    cout<<endl;
    cout<<"Something worng about bins !!!! need to check bins in this macro!!!"<<endl;
    cout<<"end this macro!!"<<endl;
    hough->AutoSave();
    outfile->Close();
    gSystem->Exit(0);
  }

  for(Int_t bin_number=0; bin_number<n_bin; bin_number++){
    frequency=houghspace->GetBinContent(bin_number);
    //    if(frequency>0){
    if(frequency>5){
      houghspace->GetBinContent(bin_number,coord);

      binlist[0].push_back(bin_number);
      binlist[1].push_back(frequency);
      binlist[2].push_back(coord[0]);
      binlist[3].push_back(coord[1]);
      binlist[4].push_back(coord[2]);
      binlist[5].push_back(coord[3]);
      
    }
  }

  binlist_raw=binlist;
  sort(binlist.begin(),binlist.end(),[](const vector<Int_t> &a, const vector<Int_t> &b){return a[1]<b[1]});

  return 0;
};

int GetLineParameter(vector<vector<Double_t>> binlist, vector<vector<Double_t>>& line_par, Int_t n_bin, Int_t n_get){

  const Double_t pi=3.141592653589;
  const Double_t deg_to_rad=pi/180;
  Double_t ax;
  Double_t ay;
  Double_t az;
  Double_t bx;
  Double_t by;
  Double_t bz;
  Double_t psi_tem;
  Double_t theta_tem;
  Double_t phi_tem;
  Double_t rho_tem;
  Double_t psi_tem_r;
  Double_t theta_tem_r;
  Double_t phi_tem_r;

  n_bin_check=binlist.size();
  //need to modify
  if(n_bin_check!=n_bin){
    cout<<endl;
    cout<<"Something worng about n_bin !!!! need to check about n_bin in this macro!!!"<<endl;
    cout<<"end this macro!!"<<endl;
    hough->AutoSave();
    outfile->Close();
    gSystem->Exit(0);
  }

  if(n_get>n_bin){
    n_get=n_bin;
  }

  line_par.resize(12,vector<Double_t>(n_get,0));

  //get parameter of line 
  for(Int_t i=0; i<n_get; i++){
    alph=binlist.at(2).at(i)*deg_to_rad;
    beta=binlist.at(3).at(i)*deg_to_rad;
    gamm=binlist.at(4).at(i)*deg_to_rad;
    delt=binlist.at(5).at(i);
    
    ax=delt*(cos(alph)*cos(gamm)-sin(alph)*cos(beta)*sin(gamm));
    ay=delt*(sin(alph)*cos(gamm)+cos(alph)*cos(beta)*sin(gamm));
    az=delt*sin(beta)*sin(gamm);
    bx=sin(alph)*sin(beta);
    by=-cos(alph)*sin(beta);
    bz=cos(beta);

    line_par.at(0).at(i)=  binlist.at(0).at(i);
    line_par.at(1).at(i)=  binlist.at(1).at(i);
    line_par.at(2).at(i)=  binlist.at(2).at(i);
    line_par.at(3).at(i)=  binlist.at(3).at(i);
    line_par.at(4).at(i)=  binlist.at(4).at(i);
    line_par.at(5).at(i)=  binlist.at(5).at(i);
    line_par.at(6).at(i)=  ax;
    line_par.at(7).at(i)=  ay;
    line_par.at(8).at(i)=  az;
    line_par.at(9).at(i)=  bx;
    line_par.at(10).at(i)= by;
    line_par.at(11).at(i)= bz;

  }

};

int detectline(int eventNumber=1){

  //set timer
  TStopwatch timer;
  timer.Start();
  
  Int_t dim=0;
  Int_t bins_set[4];
  Double_t start_val[4];
  Double_t stop_val[4];
  Double_t step_val[4];
  Double_t end_val[4];
  SetParameter(dim,start_val,stop_val_step_val,bins_set,n_get);

  Int_t bins_dim=4;
  Int_t bins[4];
  Int_t n_get;
  Int_t hitnum;
  Double_t max_z=0;
  //need to modify
  Long64_t n_bin;

  //if loop multi event, need to modify about vector "hitpoint"
  vector<vector<Double_t>> hitpoint;
  vector<vector<Double_t>> hitpoint_raw;
  vector<vector<Double_t>> binlist;
  vector<vector<Double_t>> binlist_raw;
  vector<vector<Double_t>> line_par;

  TString check_input;
  Bool_t check_bins_set_ver=false;
  #ifdef bins_set_ver
  check_bins_set_ver=true;
  #endif
  Bool_t check_bins=true;

  //generate hough space
  THnSparseD* houghspace =new THnSpareseD("houghspace","hough space",dim,bins,start_val,end_val);


  //set output
  TString outfname="data1/out_detectline"+ to_string(eventNumber) +".root";
  TFile *outfile= new TFile(outfname,"receate");
  TTree *hough= new TTree("hough","hough");
  // new branch for new file
  hough->Branch("hitpoint_raw",&hitpoint_raw);
  hough->Branch("hitpoint",&hitpoint);
  hough->Branch("line_par",&line_par);

  //for malti event, use after this macro updated.
  //  TTree *padhit= new TTree("padhit","padhit");
  //  padhit->Branch("event_number",&eve,"eve/I");
  //  padhit->Branch("hit_number",&hitnum,"hitnum/I");

  //get hits
  GetHit(eventNumber,hitnum,hitpoint);
  cout<<"finish get hits"<<endl;
  hitpoint_raw=hitpoint;
  
  //shift point
  ShiftPoint(hitpoint,max_z);
  cout<<"finish shift point. shift value:"<<max_z<<endl;

  //calculate for hough space
  CalculateHoughSpace(bins,bins_dim,startval,stop_val,bins_set,check_bins_set_ver);
  for(i=0; i<bins_dim; i++){
    if(bins[i]<=0){check_bins=false}
  }
  cout<<endl;
  cout<<"If you want to stop hough tranceform, enter \"NO\""<<endl;
  cin>>check_input;
  
  if(check_bins==false || check_input="NO"){
    if(check_bins==false){
      cout<<"Something wrong about check_bins!! bins<=0 !! need to check bin!!!"<<endl;
    }
    else{
      cout<<"input is \"NO\"!!!"<<endl;
    }
    cout<<"end this macro!!"<<endl;
    hough->AutoSave();
    outfile->Close();
    gSystem->Exit(0);
  }

  cout<<"making hough space..."<<endl;
  GenerateHoughspace(dim,bins,start_cal,end_val,step_val,houghspace,hitpoint);

  GetBinInfo(houghspace,n_bin,binlist,binlist_raw);
  delete houghspace;
  cout<<"erased hough space..."<<endl;


  GetLineParameter(binlist, line_par, n_bin, n_get);

  cout<<"number of hit:"<<hitnum<<endl;
  hough->Fill();
  hough->AutoSave();
  outfile->Close();

  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  //  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "output: " << outdata << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  
  return 0;
}

