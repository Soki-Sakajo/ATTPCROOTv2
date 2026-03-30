#include <TCanvas.h>
#include "TCanvas.h"

void pad_plot_time(Int_t eve_num=17){

  //  double fThreshold =0; //Only save hits if they're above this threshold
  Double_t q_t=1; //only fill hists if they're above this thereshold
  TChain tpc_tree("cbmsim");

  //  tpc_tree.Add("./data1/output_digi_50matm_0.root");
  tpc_tree.Add("./data1/output_digi_rcnp_sequential_hoylehoyle_0.root");
  //  tpc_tree.Add("./data1/output_digi_rcnp_sequential_hoylehoyle_70matm_0.root");

  //  TString outcanp="test.pdf";
  //  TString outcanp="can_padplane_12c_130matm_0-"+ to_string(eve_num) +".pdf";
  //  TString outcanp="./can_padpane/can_padplane_part_sequential_hoylehoyle_130matm_0-"+ to_string(eve_num) +".pdf";
  //  TString outcanpt="test.pdf";
  //  TString outcanpt="test_time.pdf";
  //  TString outcanpt="can_padplane_12c_130matm_0-"+ to_string(eve_num) +".pdf";
  //  TString outcanpt="./pdf_output/can_padplane_time_sequential_hoylehoyle_130matm_0-"+ to_string(eve_num) +".pdf";
  //  TString outcan3="test.pdf";
  //  TString outcan3="./pdf_output/can_padplane_part_sequential_hoylehoyle_130matm_0-"+ to_string(eve_num) +".pdf";
  TString outcan3="./pdf_output/can_padplane_part_sequential_hoylehoyle_70matm_0-"+ to_string(eve_num) +".pdf";

  TTreeReader reader(&tpc_tree);
  TTreeReaderValue<TClonesArray> event(reader,"AtEventH");

  //Create pad plane and load map
  TString mapFile = "Lookup20150611.xml";
  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString mapDir = dir + "/scripts/"+mapFile;
  auto fAtMapPtr = new AtTpcMap();
  fAtMapPtr->ParseXMLMap(mapDir.Data());
  fAtMapPtr->GeneratePadPlane();
  auto fPadPlane = fAtMapPtr->GetPadPlane();
  
  // generate the pad plane histogram
  gSystem->Load("/home/attpc/fair_install/ATTPCROOTv2/build/lib/libAtMap.so");
  AtMap *maptest = new AtTpcMap();
  maptest->GeneratePadPlane();
  cout<<endl;

  //  Int_t n_hist=6;
  Int_t n_hist=3;
  //  TH2Poly *padplane = maptest->GetPadPlane();
  TH2Poly *padplane1 = maptest->GetPadPlane();
  TH2Poly *padplane2 = maptest->GetPadPlane();
  TH2Poly *padplane3 = maptest->GetPadPlane();
  TH2Poly *padplane4 = maptest->GetPadPlane();
  TH2Poly *padplane5 = maptest->GetPadPlane();
  TH2Poly *padplane6 = maptest->GetPadPlane();

  Double_t x,y,z,q;
  Double_t z_max=-100;
  Double_t z_max_tem=-100;
  Double_t z_min=1e5;
  Double_t z_min_tem=1e5;
  
  //get hit pattern data from digi
  std::vector<std::vector<Short_t>> traces; // traces[event][time]

  //get the event to read
  reader.SetEntry(eve_num);
  //Get the event
  AtEvent *eventPtr =(AtEvent *)(event->At(0));
  auto numHits = eventPtr->GetNumHits();
  Int_t n_skip=0;
  Int_t n_sz=0;
  vector<vector<Double_t>> hitpoint(5,vector<Double_t>(numHits,0));
  vector<Int_t> z_val(0);
  Int_t z_point=-10000;
  Int_t z_tem=-10000;
  Int_t z_before=0;
  bool z_not_first=false;

  for(Int_t i=0; i<numHits; i++){
    auto hit=eventPtr->GetHit(i);
    auto position=hit.GetPosition();
    auto charge=hit.GetCharge();
    if(charge<q_t){
      n_skip++;
      cout<<"hitnum: "<<i<<" charge is under "<<q_t<<". skip this hit point."<<endl;
      continue;
    }
    hitpoint.at(0).at(i)=eve_num;
    hitpoint.at(1).at(i)=position.X();
    hitpoint.at(2).at(i)=position.Y();
    hitpoint.at(3).at(i)=position.Z();
    hitpoint.at(4).at(i)=charge;

    z_point=hitpoint.at(3).at(i);
    if(z_max_tem<z_point){
      z_max_tem=z_point;
    }
    if(z_min_tem>z_point){
      z_min_tem=z_point;
    }
    if(z_tem!=z_point){
      z_val.push_back(z_point);
      z_tem=z_point;
      z_not_first=true;
    }
  }
  sort(z_val.rbegin(),z_val.rend());
  z_max=z_max_tem;
  z_min=z_min_tem;

  /*
  for(Int_t i=0; i<z_val.size(); i++){
    cout<<"z value: "<<z_val.at(i)<<endl;
  }
  */
  Int_t n_z=z_val.size();
  Int_t n_p=hitpoint.at(0).size();
  //  vector<vector<Int_t>> z_hit(n_z,vector<Int_t> (3,0));
  vector<vector<Int_t>> z_hit(n_z,vector<Int_t> (4,0));
  cout<<endl;
  cout<<"skip point: "<<n_skip<<", kinds of z: "<<n_z<<", recoding point: "<<n_p<<endl;
  cout<<"range of z: "<<z_min<<" -- "<<z_max<<endl;
  cout<<"z: ";
  for(Int_t i=0; i<z_val.size(); i++){
    cout<<z_val.at(i)<<"  ";
  }
  cout<<endl;

  for(Int_t i=0; i<n_z; i++){
    z_hit.at(i).at(0)=z_val.at(i);
    for(Int_t j=0; j<n_p; j++){
      if(hitpoint.at(3).at(j)==z_val.at(i)){
	z_hit.at(i).at(1)++;
      }
    }
    if(i==0){
      z_hit.at(i).at(2)=z_hit.at(i).at(1);
      z_hit.at(i).at(3)=int(100*z_hit.at(i).at(2)/n_p+0.9);
    }
    else{
      z_hit.at(i).at(2)=z_hit.at(i-1).at(2)+z_hit.at(i).at(1);
      z_hit.at(i).at(3)=int(100*z_hit.at(i).at(2)/n_p+0.9);
    }
  }
  /*
  for(Int_t i=0; i<z_hit.size(); i++){
    cout<<"z value: "<<z_hit.at(i).at(0)<<", number of points: "<<z_hit.at(i).at(1)
	<<", sum of points from this z: "<<z_hit.at(i).at(2)
	<<", rate: "<<z_hit.at(i).at(3)<<endl;
  }
  */

  Int_t drate=int(100/n_hist+1);
  // fill in the histogram
  for(Int_t i=0; i<n_p; i++){
    x=hitpoint.at(1).at(i);
    y=hitpoint.at(2).at(i);
    z=hitpoint.at(3).at(i);
    q=hitpoint.at(4).at(i);
    padplane4->Fill(x,y,q);
    for(Int_t j=0; j<n_z; j++){
      if(z==z_hit.at(j).at(0)){
	if(z_hit.at(j).at(3)<drate){
	  padplane1->Fill(x,y,q);
	  //	  cout<<"fill padplane1 "<<z<<endl;
	}
	else if(z_hit.at(j).at(3)<drate*2){
	  padplane2->Fill(x,y,q);
	  //	  cout<<"fill padplane2 "<<z<<endl;
	}
	else if(z_hit.at(j).at(3)<drate*3){
	  padplane3->Fill(x,y,q);
	  //	  cout<<"fill padplane3 "<<z<<endl;
	}
	/*
	else if(z_hit.at(j).at(3)<drate*4){
	  padplane4->Fill(x,y,q);	  
	}
	else if(z_hit.at(j).at(3)<drate*5){
	  padplane5->Fill(x,y,q);	  
	}
	else {
	  padplane6->Fill(x,y,q);
	}
	*/
      }
    }
  }
  cout<<endl;
  /*
  //  TCanvas *can1=new TCanvas("hist1","hist1",1300,1000);
  TCanvas *can1=new TCanvas("hist1","hist1",1250,1000);
  can1->cd(1);
  can1->SetLeftMargin(0.14);
  can1->SetRightMargin(0.23);
  can1->SetTopMargin(0.08);
  can1->SetBottomMargin(0.13);
  padplane4->SetStats(0);
  padplane4->GetXaxis()->SetRangeUser(-100,100);
  padplane4->GetYaxis()->SetRangeUser(-100,100);
  padplane4->GetXaxis()->SetTitle("x [mm]");
  padplane4->GetYaxis()->SetTitle("y [mm]");
  padplane4->GetZaxis()->SetTitle("charge");
  padplane4->Draw("colz");

  padplane4->GetXaxis()->SetLabelSize(0.06);
  padplane4->GetYaxis()->SetLabelSize(0.06);
  padplane4->GetZaxis()->SetLabelSize(0.06);
  padplane4->GetXaxis()->SetTitleSize(0.06);
  padplane4->GetYaxis()->SetTitleSize(0.06);
  padplane4->GetYaxis()->SetTitleOffset(1.1);
  padplane4->GetZaxis()->SetTitleSize(0.06);
  padplane4->GetZaxis()->SetTitleOffset(1.6);
  gStyle->SetOptTitle(0);
  gPad->Update();

  can1->SaveAs(outcanp);

  //  TCanvas *can2=new TCanvas("hist","hist",3750,1000);
  TCanvas *can2=new TCanvas("hist","hist",3800,1000);
  can2->Divide(3,1);
  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetTitleX(0.45);

  can2->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.23);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.13);
  padplane1->SetStats(0);
  padplane1->SetTitle("(a) 0 < z #leq 60");
  padplane1->GetXaxis()->SetTitle("x [mm]");
  padplane1->GetYaxis()->SetTitle("y [mm]");
  padplane1->GetZaxis()->SetTitle("charge");
  padplane1->Draw("colz");
  padplane1->GetXaxis()->SetLabelSize(0.06);
  padplane1->GetYaxis()->SetLabelSize(0.06);
  padplane1->GetZaxis()->SetLabelSize(0.06);
  padplane1->GetXaxis()->SetTitleSize(0.06);
  padplane1->GetYaxis()->SetTitleSize(0.06);
  padplane1->GetYaxis()->SetTitleOffset(1.1);
  padplane1->GetZaxis()->SetTitleSize(0.06);
  padplane1->GetZaxis()->SetTitleOffset(1.6);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can2->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.23);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.13);
  padplane2->SetStats(0);
  padplane2->SetTitle("(b) 60 < z #leq 80");
  padplane2->GetXaxis()->SetTitle("x [mm]");
  padplane2->GetYaxis()->SetTitle("y [mm]");
  padplane2->GetZaxis()->SetTitle("charge");
  padplane2->Draw("colz");
  padplane2->GetXaxis()->SetLabelSize(0.06);
  padplane2->GetYaxis()->SetLabelSize(0.06);
  padplane2->GetZaxis()->SetLabelSize(0.06);
  padplane2->GetXaxis()->SetTitleSize(0.06);
  padplane2->GetYaxis()->SetTitleSize(0.06);
  padplane2->GetYaxis()->SetTitleOffset(1.1);
  padplane2->GetZaxis()->SetTitleSize(0.06);
  padplane2->GetZaxis()->SetTitleOffset(1.6);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can2->cd(3);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.23);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.13);
  padplane3->SetStats(0);
  padplane3->SetTitle("(c) 80 < z #leq 150");
  padplane3->GetXaxis()->SetTitle("x [mm]");
  padplane3->GetYaxis()->SetTitle("y [mm]");
  padplane3->GetZaxis()->SetTitle("charge");
  padplane3->Draw("colz");
  padplane3->GetXaxis()->SetLabelSize(0.06);
  padplane3->GetYaxis()->SetLabelSize(0.06);
  padplane3->GetZaxis()->SetLabelSize(0.06);
  padplane3->GetXaxis()->SetTitleSize(0.06);
  padplane3->GetYaxis()->SetTitleSize(0.06);
  padplane3->GetYaxis()->SetTitleOffset(1.1);
  padplane3->GetZaxis()->SetTitleSize(0.06);
  padplane3->GetZaxis()->SetTitleOffset(1.6);
  //  gStyle->SetOptTitle(0);
  gPad->Update();

  can2->SaveAs(outcanpt);
  */
  TCanvas *can1=new TCanvas("hist1","hist1",1000,1000);
  can1->cd(1);
  /*
  can1->SetLeftMargin(0.14);
  can1->SetRightMargin(0.23);
  can1->SetTopMargin(0.08);
  can1->SetBottomMargin(0.13);
  */
  padplane4->SetStats(0);
  padplane4->GetXaxis()->SetRangeUser(-150,150);
  padplane4->GetYaxis()->SetRangeUser(-150,150);
  //  padplane4->GetXaxis()->SetRangeUser(-100,100);
  //  padplane4->GetYaxis()->SetRangeUser(-100,100);
  padplane4->GetXaxis()->SetTitle("");
  padplane4->GetYaxis()->SetTitle("");
  padplane4->GetZaxis()->SetTitle("");
  padplane4->GetXaxis()->SetLabelSize(0);
  padplane4->GetYaxis()->SetLabelSize(0);
  padplane4->GetZaxis()->SetLabelSize(0);
  padplane4->GetXaxis()->SetTickLength(0);
  padplane4->GetYaxis()->SetTickLength(0);
  padplane4->GetZaxis()->SetTickLength(0);
  padplane4->GetZaxis()->SetAxisColor(0);
  gStyle->SetOptTitle(0);
  padplane4->Draw("col");
  gPad->Update();

  can1->SaveAs(outcan3);

}
