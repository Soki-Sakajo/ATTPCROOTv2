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

int sim_macro_furuno(Int_t subnum=0){ 

  //set timer
  TStopwatch timer;
  timer.Start();
  
  int ne=0;
  
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
  
  tree= (TTree*) file->Get("cbmsim");
  
  // set reading data
  AtMCTrack* track=new AtMCTrack();
  TClonesArray *trackArray=0;
  tree->SetBranchAddress("MCTrack",&trackArray);

  // furuno
  //  AtTpc* point = new AtTpc();
  //  AtMCPoint* point = new AtMCPoint();
  FairMCPoint* point = new FairMCPoint();    
  TClonesArray *pointArray=0;
  tree->SetBranchAddress("AtTpcPoint",&pointArray);
  
  Int_t nEvents=tree->GetEntriesFast();
  
  //  cout << nEvents <<endl;

  // analysis of attpcsim.root
  for(int i=0;i < nEvents; i++){
    if((i-1)%1000==0) {
      cout<<"Finish "<<100*(i-1)/nEvents <<" % ! \r"<<flush;
    }
    int eve = i;    
    ne++;

    
    //get event
    tree->GetEvent(i);
    int n=pointArray->GetEntries();

    if(i<10) printf("npoint = %d\n", n);
    
    Double_t px[2];
    Double_t py[2];
    Double_t pz[2]; 

    Double_t x[2];
    Double_t y[2];
    Double_t z[2];

    int trackid;
    
    if(i%2==1){  // analyze only the scattering event
      for(int j=0; j<n; j++){
	point = (FairMCPoint*)pointArray->At(j); // get the point information
	
	trackid = (int)(point->GetTrackID());  // get track ID
	
	if(trackid==0 || trackid==1){
	  // positions
	  x[trackid]=point->GetX();
	  y[trackid]=point->GetY();
	  z[trackid]=point->GetZ();      
	  
	  // momenta
	  px[trackid]=point->GetPx();
	  py[trackid]=point->GetPy();
	  pz[trackid]=point->GetPz();      
	}
      }
      
      if(i==1){
	printf("ID0: X:%f, Y:%f, Z:%f\n", x[0], y[0], z[0]);
	printf("ID1: X:%f, Y:%f, Z:%f\n", x[1], y[1], z[1]);      
      }
    }
  }
  
  return 0;
}
