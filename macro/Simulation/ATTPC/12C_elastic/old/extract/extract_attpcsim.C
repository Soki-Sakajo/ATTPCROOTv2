#include <TTree.h>
#include <TFile.h>
#include"TString.h"
#include"TFile.h"
#include"TTree.h"

int extract(Int_t subnum=0, Int_t nEvent=0){ 

  //extract simfile
  //  input
  TString insim = "data1/attpcsim_12c12c_" + to_string(subnum) + ".root";
  TFile *simfile =new TFile(insim);  
  TTree *simtree = (TTree*)simfile->Get("cbmsim");
  Int_t evesim=simtree->GetEntries();

  //  set output
  TString outsim = "data1/extract_c12_" + to_string(subnum)+"_"+ to_string(nEvent)+ ".root";
  TFile *exsim=new TFile(outsim, "recreate");
  TTree *exsimtree=0;  
  
  exsimtree =simtree->CloneTree(0);
  simtree->GetEntry(nEvent-1);
  exsimtree->Fill();
  simtree->GetEntry(nEvent);
  exsimtree->Fill();
  simtree->GetEntry(nEvent+1);
  exsimtree->Fill();  

  exsim->Write();
  simfile->Close();
  exsim->Close();

  /*
  //extract parfile
  //  input
  TString inpar = "data1/attpcsim_12c12c_" + to_string(subnum) + ".root";
  TFile *parfile =new TFile(infname);  
  TTree *partree = (TTree*)file->Get("cbmsim");
  Int_t evepar=partree->GetEntries();

  //  set output
  TString outpar = "data1/extract_c12_" + to_string(subnum)+"_"+ to_string(nEvent)+ ".root";
  TFile *exsim=new TFile(outfname, "recreate");
  TTree *exsitree=0;  

  extree = tree->CloneTree(0);
  tree->GetEntry(nEvent-1);
  extree->Fill();
  tree->GetEntry(nEvent);
  extree->Fill();
  tree->GetEntry(nEvent+1);
  extree->Fill();

  expar->Write();
  parfile->Close();
  expar->Close();
  */

  return 0;
}
