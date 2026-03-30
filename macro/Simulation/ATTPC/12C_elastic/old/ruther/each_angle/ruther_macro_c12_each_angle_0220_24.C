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
#include<vector>
#include"TChain.h"

int ruther_macro_c12_each_angle_0220_24(Int_t subnum=1){ 
  
  //This file is only 12C12C elastic scattering
  
  //set timer
  TStopwatch timer;
  timer.Start();
  
  //set situation
  Int_t z_s=6;          // Z of scatterd particl in simfile
  Int_t z_r1=6;         // Z of an recoil particl in simfile
  Int_t z_r2=1;         // Z of another recoil particl in simfile
  Int_t dz=1;           // mm; bin of depth
  Int_t n_b=1e5;        // /s; beam intensity
  Double_t T=300;       // K; thermodynamic temparature
  Double_t P=0.13;      // atm; gas pressure
  Double_t da=1;        // deg; bin of angle
  Double_t dr=1;        // mm; bin of reaction rate
  Int_t nbin=int(300/dr);    // bin of rmax-reaction_rate hist
  Int_t nbin_zn=int(400/dz); // bin of z-reaction_rate hist
  Int_t N_rate1=4;      // (number of scatterd particle)/molecular
  Int_t N_rate2=10;     // (number of recoli partcle)/molecular  
  Int_t theta_min=3;    // deg; rutherford scattar angle
  Int_t r_attpc=275;    // mm; radius of attpc
  Int_t l_attpc=1000;   // mm; depth of attpc
  
  //set constants
  int ctom=10;
  Double_t pi=3.1415926535;
  Double_t rad_to_deg=180/pi;
  Double_t deg_to_rad=pi/180;
  Double_t atmtoPa=101325; // Pa
  Double_t alfa=7.297352569e-3; // fine structure constant
  Double_t hc=197.463587*1e-15; // MeV•m
  Double_t k_b=1.38065e-23; // J/K ; Bolzman constant
  Double_t rho=P*atmtoPa/(k_b*T);
  Double_t m12=11177.73600; // MeV/c^2 ; 12C mass
  
  //constants in ruthford scatter
  Double_t a=alfa*hc*z_s*z_r1/4;  
  
  //set variable
  Int_t div_l=(int)(l_attpc/dz+0.9);
  Int_t div_r=(int)(r_attpc/dr+0.9);
  Int_t div_a=(int)((180-theta_min)/da+0.9);
  Int_t eve;
  Int_t n_data;
  Int_t n_p1;
  Int_t n_fi1;
  Int_t n_ev1;
  Int_t n_ef1p=0;
  Int_t n_ee1p=0;
  Int_t n_ef1l=0;
  Int_t n_ee1l=0;
  Int_t n_ef1a=0;
  Int_t n_ee1a=0;
  Int_t n_ef1b=0;
  Int_t n_ee1b=0;
  Int_t n_ef1o=0;
  Int_t n_ee1o=0;
  Int_t n_fi_temp=0;
  Int_t eve_temp=0;
  //  Int n_k1;

  Double_t rmax1;
  Double_t E1;
  Double_t startz1[2];  
  Double_t z1[2];
  Double_t theta1_cms[2];
  Double_t theta1_cms_r;
  Double_t Y1_rmax=0;
  Double_t da_r=da*deg_to_rad;
  Double_t l1[div_l];
  Double_t ang1[div_a];
  Double_t ang_r1[div_a];
  Double_t n1[div_l];
  Double_t Esum1[div_l];
  Double_t Eave1[div_l];
  
  vector<vector<vector<Int_t>>> fi_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Int_t>>> ev_ve1(div_l,vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  vector<vector<vector<Double_t>>> E_ve(div_l,vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Double_t>>> r_ve(div_l, vector<vector<Double_t>>(div_a,vector<Double_t>(0)));
  vector<vector<vector<Int_t>>> np_ve(div_l, vector<vector<Int_t>>(div_a,vector<Int_t>(0)));
  //  vector<vector<vector<Double_t>>> Y1_ve(div_l, vector<vector<Int_t>>(div_a,vector<Int_t>(0)));  
  
  Int_t z_n=0;
  Int_t a_n=0;
  //  Int_t nEntry=0;
  
  vector<vector<Double_t>> n_a1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> I(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> n_beam1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> sigma1(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_da(div_l, vector<Double_t>(div_a,0));
  vector<vector<Double_t>> Y1_nom(div_l, vector<Double_t>(div_a,0));
  
  Double_t n_react1[div_l];
  Int_t n_rmax1=0;
  Double_t Y1[div_l];
  Double_t Y1_sum=0;
  Double_t z[div_r];
  Double_t z1_max;
  Double_t i1_max;
  Double_t j1_max;
  Double_t dep1;
  Double_t dir1;
  Double_t Y1_eve;
  Double_t theta1;

  bool setbeam=true;
  bool fin_beam1=false;
  bool fin_rea1=false;
  
  //for 12Cp
  Double_t b=alfa*hc*z_s*z_r2/4;  
  Double_t rmax2;
  Double_t E2;
  Double_t startz2[2];  
  Double_t z2[2];
  Double_t theta2_cms[2];
  Int_t n_p2;
  Int_t n_rmax2[div_l];   
  Double_t n2[div_l];
  Double_t n_a2[div_a];
  Double_t Esum2[div_l];
  Double_t Eave2[div_l];
  Double_t Y2[div_l];

  
  //input multi-file by using TChain
  // Using this or making root files
  // Tchain of c12

  TChain *tree1=new TChain("anadata","anadata");
    for(int i=0; i<subnum; i++){
      //      TString refi="data/sim_macro_c12_" + to_string(i) + ".root";
      TString dir="data1";
      TString refi="sim_macro_c12_" + to_string(i) + ".root";
      TString direfi=dir+"/"+refi;
      tree1->Add(direfi);
      if(i==0){cout<<"read from  "<<refi<<flush;}
      if(i==subnum-1){cout<<"  to  "<<refi<<endl;}
    }
  cout<<"end reading files"<<endl<<endl;
  if(subnum==1){
    cout<<endl;
    cout<<"ATTENTION!! read ONLY one file!!"<<endl<<endl<<endl;
  }
  
  // Tchain of p for 12Cp
  /*
    TChain *trp=new TChain("anadata","anadata");
    for(i=0, i<subnum){
    trp->Add("data1/sim_macro_p_" + to_string(subnum) + ".root");
    }
    //trp->Add("data1/sim_macro_p_0.root");
   */  

  // set reading data
  Int_t Entry1=tree1->GetEntries();
  tree1->SetBranchAddress("event_number",&eve);
  tree1->SetBranchAddress("startz",&startz1);
  tree1->SetBranchAddress("rmax",&rmax1);
  tree1->SetBranchAddress("TotalE",&E1);
  tree1->SetBranchAddress("theta_cms",&theta1_cms);
  tree1->SetBranchAddress("n",&n_p1);

  //  cout<<"set branch adress"<<endl;
 
  Int_t n_f=Entry1/1e5;
  
  //  cout<<"before set histgram"<<endl;
  
  //set histgram
  TH1F *h_rmax=new TH1F("h_rmax","rmax:reaction_rate",nbin,0,300);
  TH1F *h_z=new TH1F("h_z","reaction position:reaction_rate",nbin_zn,0,400);
  TH2F *h_ztheta=new TH2F("h_ztheta","reaction position:angle",400,0,400,100,0,100);
  TH1F *h_thetaY1=new TH1F("h_thetaY1","angle :reaction_rate",100,0,100);
  TH1F *h_thetaY2=new TH1F("h_thetaY2","angle :reaction_rate",100,0,100);
  TH1F *h_thetaY3=new TH1F("h_thetaY3","angle :reaction_rate",100,0,100);
  
  //set output for making ratdata by each subnum
  //  TString outdata = "data/rutdata_c12_vector.root";
  //  TString outdata = "data/rutdata_c12_" + to_string(subnum) + ".root";
  TString outdata = "data1/rutdata_c12_" + to_string(subnum) + ".root";
  TFile *rutfile=new TFile(outdata, "recreate");
  TTree *rutdata = new TTree("rutdata","rutdata");
  /*
  TTree *exp12cp = new TTree("exp12cp","exp12cp");
  TTree *exp12cl = new TTree("exp12cl","exp12cl");
  TTree *exp12ca = new TTree("exp12ca","exp12ca");
  TTree *exp12cb = new TTree("exp12cb","exp12cb");
  TTree *exp12co = new TTree("exp12co","exp12co");
  */
  
  rutdata->Branch("file_number",&n_fi1,"eve/I");
  rutdata->Branch("event_number",&n_ev1,"eve/I");
  rutdata->Branch("reaction_position",&dep1,"dep1/D");
  rutdata->Branch("angle",&dir1,"dir1/D");
  rutdata->Branch("rmax_c12",&rmax1,"rmax1/D");
  rutdata->Branch("reaction_rate",&Y1_eve,"Y1_eve/D");
  //  rutdata->Branch("depth",&l1,"l1/D");
  rutdata->Branch("Eave",&Eave1,"Eave1/D");

  /*
  exp12cp->Branch("file_particle",&n_ef1p,"n_ef1p/I");
  exp12cp->Branch("event_number_particle",&n_ee1p,"n_ee1p/I");
  exp12cl->Branch("file_startz",&n_ef1l,"n_ef1l/I");
  exp12cl->Branch("event_number_startz",&n_ee1l,"n_ee1l/I");
  exp12ca->Branch("file_angle",&n_ef1a,"n_fi1a/I");
  exp12ca->Branch("event_number_angle",&n_ee1a,"n_ee1a/I");
  exp12cb->Branch("file_angle",&n_ef1b,"n_fi1b/I");
  exp12cb->Branch("event_number_angle",&n_ee1b,"n_ee1b/I");
  exp12co->Branch("file_other",&n_ef1o,"n_fi1o/I");
  exp12co->Branch("event_number_other",&n_ee1o,"n_ee1o/I");
  */
  
  //make event lists  
  for(int j=0; j<Entry1; j++){
    tree1->GetEntry(j);
    z1[0]=startz1[0]*ctom;
    z1[1]=startz1[1]*ctom;
    if(eve_temp>eve){n_fi_temp++;}

    if(z1_max<z1[0]){z1_max=z1[0];}
    if(z1_max<z1[1]){z1_max=z1[1];}
    if(theta1_cms[0]<=theta1_cms[1]){theta1=theta1_cms[0];}
    else if(theta1_cms[0]>theta1_cms[1]){theta1=theta1_cms[1];}

    if(n_p1==2 && fabs(z1[0]-z1[1])<dz && theta1-theta_min>0){
      z_n=int(z1[0]/dz);
      a_n=int((theta1-theta_min)/da);
      ev_ve1[z_n][a_n].push_back(eve);
      E_ve[z_n][a_n].push_back(E1);
      r_ve[z_n][a_n].push_back(rmax1);
      np_ve[z_n][a_n].push_back(n_p1);
      //      cout<<j<<"  "<<n_ve.at(z_n).at(a_n).at(0)<<"   "<<E_ve.at(z_n).at(a_n).at(0)<<endl;
      //      nEntry++;

      fi_ve1[z_n][a_n].push_back(n_fi_temp);

    }
    /*
    else if(fabs(z1[0]-z1[1])<dz && theta1-theta_min>0){
      n_ef1p=n_fi_temp;
      n_ee1p=eve;
      exp12cp->Fill();
      //      cout<<"particle  "<<eve<<endl;
    }
    else if(n_p1==2 && theta1-theta_min>0){
      n_ef1l=n_fi_temp;
      n_ee1l=eve;
      exp12cl->Fill();
      //      cout<<"length  "<<eve<<endl;
    }
    else if(n_p1==2 && fabs(z1[0]-z1[1])<dz){
      n_ef1a=n_fi_temp;
      n_ee1a=eve;
      exp12ca->Fill();
      //      cout<<"angle  "<<eve<<endl;
    }
    else if(eve%2==0){
      n_ef1b=n_fi_temp;
      n_ee1b=eve;
      exp12cb->Fill();
      //      cout<<"beam  "<<eve<<endl;
    }
    else{
      n_ef1o=n_fi_temp;
      n_ee1o=eve;
      exp12co->Fill();
      //      cout<<"else  "<<eve<<endl;
    }
    */
    
    eve_temp=eve;
  }
  
  i1_max=int(z1_max/dz);
  
  //12c12c calculate
  // devide depth
  
  //  for(int i=0;i < div_l ; i++){
  for(int i=0;i < ev_ve1.size(); i++){
    l1[i]=i*dz;
    n1[i]=0;
    Esum1[i]=0;    
    Eave1[i]=0;
    
    //calculate Eave
    for(int j=0; j< ev_ve1.at(i).size(); j++){
      for(int k=0; k< ev_ve1.at(i).at(j).size(); k++){
	//	if(np_ve.at(i).at(j).at(k)==1){continue;}
	n1[i]++;
	Esum1[i]+=E_ve.at(i).at(j).at(k);
	if(i==i1_max){
	  if(ev_ve1.at(i).at(j).size()>0){j1_max=j;}
	}
      }
    }
    if(n1[i]<1){continue;}
    else{Eave1[i]=Esum1[i]/(2*n1[i]);}
    
    //calculate reaction rate
    Y1[i]=0;
    n_react1[i]=0;
    
    //divide angle
    for(int j=0; j< ev_ve1.at(i).size(); j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;

      //need
      cout<<"calculating reaction rate at "<<l1[i]<<" mm "<<ang1[j]<<" deg  \r"<<flush;
      
      n_a1.at(i).at(j)++;
      //	    cout<<j<<"  "<<i<<"   "<<k<<"  "<<n_a1.at(i).at(k)<<endl;

      if(ev_ve1.at(i).at(j).size()<1){continue;}
      
      I.at(i).at(j)=2*pi*a*a*da_r*sin(ang_r1[j])
	/(Eave1[i]*Eave1[i]
	  *sin(ang_r1[j]/2)*sin(ang_r1[j]/2)
	  *sin(ang_r1[j]/2)*sin(ang_r1[j]/2))

	+2*pi*a*a*da_r*sin(ang_r1[j])
	/(Eave1[i]*Eave1[i]
	  *cos(ang_r1[j]/2)*cos(ang_r1[j]/2)
	  *cos(ang_r1[j]/2)*cos(ang_r1[j]/2));

	//	+4*pi*a*a*cos(sqrt(2*m12*Eave1[i])/hc*log(tan(ang_r1[j]/2)*tan(ang_r1[j]/2)))/sin(ang_r1[j]);
	//	+4*pi*a*a*cos(sqrt(2*m12/Eave1[i])*a/hc*log(tan(ang_r1[j]/2)*tan(ang_r1[j]/2)))/sin(ang_r1[j]);
       
      if(setbeam){
	n_beam1.at(i).at(j)=n_b;
	setbeam=false;
      }
      else{n_beam1.at(i).at(j)=n_b-Y1_sum;}
      //      cout<<i<<"       "<<j<<"      "<<n_beam1.at(i).at(j)<<endl;
      if(n_beam1.at(i).at(j)<0){
      cout<<"                                                                 "<<endl;
      cout<<"end the loop at"<<l1[i]<<" mm~"<<l1[i]+dz<<" mm, "<<ang1[j]<<"deg~"<<ang1[j]+da<<"deg, because beam particle is 0."<<endl<<endl;
	fin_beam1=true;
	break;
      }
      else{
	
	//	sigma1.at(i).at(j)=2*pi*a*a*I.at(i).at(j)/(Eave1[i]*Eave1[i]);
	//	sigma1.at(i).at(j)=a*a*I.at(i).at(j)/(Eave1[i]*Eave1[i]);
	//	Y1_da.at(i).at(j)=sigma1.at(i).at(j)*n_beam1.at(i).at(j)*rho*dz*1e-3*N_rate1;
	Y1_da.at(i).at(j)=I.at(i).at(j)*n_beam1.at(i).at(j)*rho*dz*1e-3*N_rate1;
	//	Y1_nom.at(i).at(j)+=Y1_da.at(i).at(j)/n_a1.at(i).at(j);
	Y1[i]+=Y1_da.at(i).at(j);
	Y1_sum+=Y1_da.at(i).at(j);
	//	cout<<i<<"   "<<k<<"          "<<Y1_da.at(i).at(k)<<"          "<<Y1[i]<<"          "<<Y1_sum<<endl;
	
      }
      //      cout<<i<<"       "<<j<<"      "<<ev_ve.at(i).at(j).size()<<endl;
      if(i==i1_max && j==j1_max){
      cout<<"                                                                 "<<endl;
      cout<<"end the loop at "<<l1[i]<<" mm~"<<l1[i]+dz<<" mm "
	  <<ang1[j]<<" deg~"<<ang1[j]+da<<" deg, because of no reaction in over "
	  <<l1[i]+dz<<" mm "<<ang1[j]+da<<" deg"<<endl<<endl;
      //      cout<<i<<"    Y1sum"<<Y1_sum<<endl;
      fin_rea1=true;
      break;
    }

    }
    
    //z-reaction rate hist filling data
    if(n1[i]>0){
      //      h_z->Fill(l11[i],Y1_nom[i]);
      h_z->Fill(l1[i],Y1[i]);
    }
    if (i%100==0){
      //need
      cout<<"12C12C finish "<<l1[i]<<" mm !!                     "<<endl;
    }
    if(fin_beam1){break;}
    if(fin_rea1){break;}
    //    if(z1_max<l1[i]+dz){


  }

  //make conections between rmax and reaction rate
  for(int i=0; i< ev_ve1.size(); i++){
    l1[i]=i*dz;
    
    for(int j=0;j < ev_ve1.at(i).size() ; j++){
      ang1[j]=j*da+theta_min;
      ang_r1[j]=ang1[j]*deg_to_rad;
      //      n_k1=0;

      for(int k=0; k< ev_ve1.at(i).at(j).size() ; k++){
	//	n_k1++;
	n_fi1=fi_ve1.at(i).at(j).at(k);
	n_ev1=ev_ve1.at(i).at(j).at(k);

	dep1=l1[i];
	dir1=ang1[j];
	Y1_eve=Y1_da.at(i).at(j)/ev_ve1.at(i).at(j).size();
	rmax1=r_ve.at(i).at(j).at(k);
	h_rmax->Fill(rmax1,Y1_eve);
	n_rmax1++;
	Y1_rmax+=Y1_eve;
	h_ztheta->Fill(l1[i],ang1[j]);
	rutdata->Fill();

	if(l1[i]==5){
	  h_thetaY1->Fill(ang1[j],Y1_eve);
	    }
	if(l1[i]==20){
	  h_thetaY2->Fill(ang1[j],Y1_eve);
	    }
	if(l1[i]==30){
	  h_thetaY3->Fill(ang1[j],Y1_eve);
	    }
      }
    }
    //need
    cout<<"Filling data "<<100*i/ev_ve1.size() <<" % !  \r"<<flush;
  }
  
  //  cout<<"Entry "<<nEntry<<"   n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax<<"                     "<<endl;
    cout<<"n_rmax1 "<<n_rmax1<<"    Y1sum "<<Y1_rmax<<"                       "<<endl;
  
  //search rmax when over 100Hz
  double hist_cnt = 0;
  int flag_100hz=0;
  for(int ihist=nbin; ihist>0; ihist--){
    hist_cnt += h_rmax->GetBinContent(ihist);
    //    printf("cnt=%f\n", h_rmax->GetBinContent(nbin));
    if(hist_cnt > 100 && flag_100hz==0 ){
      printf("over 100Hz %d mm \n\n", ihist);
      flag_100hz=1;
    }
  }
  
  //make error bar
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_rmax->GetBinContent(m));
    h_rmax->SetBinError(m,err);
  }
  for(int m=0; m<=nbin_zn;m++){
    double err_z=sqrt(h_z->GetBinContent(m));
    h_z->SetBinError(m,err_z);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_ztheta->GetBinContent(m));
    h_ztheta->SetBinError(m,err);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_thetaY1->GetBinContent(m));
    h_thetaY1->SetBinError(m,err);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_thetaY2->GetBinContent(m));
    h_thetaY2->SetBinError(m,err);
  }
  for(int m=0; m<=nbin;m++){
    double err=sqrt(h_thetaY3->GetBinContent(m));
    h_thetaY3->SetBinError(m,err);
  }

  //make and save the hist 
  TCanvas *can=new TCanvas("hists","hists");
  can->Divide(2,3);
  
  // z-angle
  //  h_ztheta->SetStats(0);
  h_ztheta->SetStats(1);
  can->cd(1);
  h_ztheta->Draw("colz");

  // z-reaction rate
  //  h_z->SetStats(0);
  h_z->SetStats(1);
  can->cd(2);
  gPad->SetLogy(1);
  h_z->Draw("hist");

  // rmax-reaction rate
  //  h_rmax->SetStats(0);
  h_rmax->SetStats(1);
  can->cd(3);
  gPad->SetLogy(1);
  h_rmax->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY1->SetStats(1);
  can->cd(4);
  gPad->SetLogy(1);
  h_thetaY1->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY2->SetStats(1);
  can->cd(5);
  gPad->SetLogy(1);
  h_thetaY2->Draw("hist");

  // theta-reaction rate
  //  h_thetaY->SetStats(0);
  h_thetaY3->SetStats(1);
  can->cd(6);
  gPad->SetLogy(1);
  h_thetaY3->Draw("hist");

  //output for canvas
  //  TString outcan = "data/ruther_12c_can_vector.C";
  //  TString outcan = "data/ruther_12c_can_" + to_string(subnum) + ".C";
  TString outcan = "data1/ruther_12c_can_" + to_string(subnum) + ".C";
  can->SaveAs(outcan);
  h_ztheta->Write();
  h_z->Write();
  h_rmax->Write();
  rutdata->AutoSave();
  /*
  exp12cp->AutoSave();
  exp12cl->AutoSave();
  exp12ca->AutoSave();
  exp12cb->AutoSave();
  exp12co->AutoSave();
  */
  rutfile->Close();
  
  // stop timer
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl;
  cout << "output: " << outdata <<" and "<< outcan << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
  return 0;
}
