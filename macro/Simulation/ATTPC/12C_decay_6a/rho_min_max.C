//#define check_psi
#define check_theta
//#define check_phi

#include <iostream>
#include "TMath.h"
#include "TCanvas.h"

Double_t rho(const Double_t *vars){
  Double_t a = vars[0]; //psi
  Double_t b = vars[1]; //theta
  Double_t c = vars[2]; //phi
  Double_t x = vars[3];
  Double_t y = vars[4];
  Double_t z = vars[5];
  
  Double_t A = cos(a)*cos(c)-sin(a)*cos(b)*sin(c);
  Double_t B = sin(a)*cos(c)+cos(a)*cos(b)*sin(c);
  Double_t C = sin(b)*sin(c);
  
  return x*A+y*B+z*C;
 };

int rho_min_max(Int_t subnum=0){

  //set timer
  TStopwatch timer;
  timer.Start();

  const Double_t pi=3.141592653589;
  const Double_t rad_to_deg=180/pi;
  const Double_t deg_to_rad=pi/180;

  //set condition
  const Double_t start_val[6]={0,0,0,-280,-280,0};
  const Double_t stop_val[6]={360,90,360,280,280,1000};
  const Double_t check_val_max=1075.54;
  const Double_t check_val_min=-1075.54;
  //  const Double_t step_val[6]={360,90,360,560,560,1000};
  //  const Double_t step_val[6]={180,45,180,280,280,500};
  //  const Double_t step_val[6]={36,9,36,56,56,100};
  //  const Double_t step_val[6]={18,4.5,18,28,28,50};
  //  const Double_t step_val[6]={3.6,0.9,3.6,5.6,5.6,10};
  //  const Double_t step_val[6]={1,0.25,1,560,560,1000};
  const Double_t step_val[6]={0.36,0.09,0.36,560,560,1000};

  
  Double_t deff_val[6];
  for(int i=0; i<6; i++){
    deff_val[i]=stop_val[i]-start_val[i];
  }
  const Int_t div_psi=  (int)(deff_val[0]/step_val[0]+0.9);
  const Int_t div_theta=(int)(deff_val[1]/step_val[1]+0.9);
  const Int_t div_phi=  (int)(deff_val[2]/step_val[2]+0.9);
  const Int_t div_x=    (int)(deff_val[3]/step_val[3]+0.9);
  const Int_t div_y=    (int)(deff_val[4]/step_val[4]+0.9);
  const Int_t div_z=    (int)(deff_val[5]/step_val[5]+0.9);
  Double_t rho_tem=0;
  Double_t rho_max=0;
  Double_t rho_min=0;
  Double_t tem_coord[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t tem_coo_r[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t max_coord[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t max_coo_r[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t min_coord[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t min_coo_r[6]={-1,-1,-1,-1000,-1000,-1000};
  Double_t psi[div_psi+1];
  Double_t psi_r[div_psi+1];
  Double_t theta[div_theta+1];
  Double_t theta_r[div_theta+1];
  Double_t phi[div_phi+1];
  Double_t phi_r[div_phi+1];
  Double_t x[div_x+1];
  Double_t y[div_y+1];
  Double_t z[div_z+1];
  Double_t v_rho_max=0;
  Double_t v_rho_min=0;
  Double_t v_psi_max,v_theta_max,v_phi_max,v_x_max,v_y_max,v_z_max;
  Double_t v_psi_min,v_theta_min,v_phi_min,v_x_min,v_y_min,v_z_min;

  //set output
  TString outdata="data1/cal_rho_"+to_string(subnum)+".root";
  TFile *calfile= new TFile(outdata,"recreate");
  TTree *calmax= new TTree("calmax","calmax");
  calmax->Branch("max_rho",&v_rho_max,"v_rho_max/D");
  calmax->Branch("max_psi",&v_psi_max,"v_psi_max/D");
  calmax->Branch("max_theta",&v_theta_max,"v_theta_max/D");
  calmax->Branch("max_phi",&v_phi_max,"v_phi_max/D");
  calmax->Branch("max_x",&v_x_max,"v_x_max/D");
  calmax->Branch("max_y",&v_y_max,"v_y_max/D");
  calmax->Branch("max_z",&v_z_max,"v_z_max/D");
  TTree *calmin= new TTree("calmin","calmin");
  calmin->Branch("min_rho",&v_rho_min,"v_rho_min/D");
  calmin->Branch("min_psi",&v_psi_min,"v_psi_min/D");
  calmin->Branch("min_theta",&v_theta_min,"v_theta_min/D");
  calmin->Branch("min_phi",&v_phi_min,"v_phi_min/D");
  calmin->Branch("min_x",&v_x_min,"v_x_min/D");
  calmin->Branch("min_y",&v_y_min,"v_y_min/D");
  calmin->Branch("min_z",&v_z_min,"v_z_min/D");

  cout<<endl;
  cout<<"\"start\" is the start point of calculation."<<endl;
  cout<<"\"stop\" is the goal point of calculation."<<endl;
  cout<<"\"end\" is actually the end point of calculation."<<endl;
  cout<<"\"bin\" is the number of calculation points."<<endl<<endl;

  cout<<"start (psi,theta,phi,x,y,z)= ("<<start_val[0]<<", "<<start_val[1]<<", "<<start_val[2]<<", "
      <<start_val[3]<<", "<<start_val[4]<<", "<<start_val[5]<<")"<<endl;
  cout<<"stop  (psi,theta,phi,x,y,z)= ("<<stop_val[0]<<", "<<stop_val[1]<<", "<<stop_val[2]<<", "
      <<stop_val[3]<<", "<<stop_val[4]<<", "<<stop_val[5]<<")"<<endl;
  cout<<"end   (psi,theta,phi,x,y,z)= ("<<div_psi*step_val[0]+start_val[0]<<", "<<div_theta*step_val[1]+start_val[1]
      <<", "<<div_phi*step_val[2]+start_val[2]<<", "<<div_x*step_val[3]+start_val[3]
      <<", "<<div_y*step_val[4]+start_val[4]<<", "<<div_z*step_val[5]+start_val[5]<<")"<<endl;
  cout<<"step  (psi,theta,phi,x,y,z)= ("<<step_val[0]<<", "<<step_val[1]<<", "<<step_val[2]<<", "
      <<step_val[3]<<", "<<step_val[4]<<", "<<step_val[5]<<")"<<endl;
  cout<<"bin   (psi,theta,phi,x,y,z)= ("<<div_psi+1<<", "<<div_theta+1<<", "<<div_phi+1<<", "
      <<div_x+1<<", "<<div_y+1<<", "<<div_z+1<<")"<<endl<<endl;
  
  cout<<"checking value under "<<check_val_min<<" and over "<<check_val_max<<" for output: "<<outdata<<endl<<endl;

  #ifdef check_psi
  cout<<"calculating ... psi 0% compleated now      \r"<<flush;
  #endif
  #ifdef check_theta
  cout<<"calculating ... psi 0%, theta 0% compleated now      \r"<<flush;
  #endif
  #ifdef check_phi
  cout<<"calculating ... psi 0%, theta 0%, phi 0% compleated now      \r"<<flush;
  #endif

  for(int i=0; i<div_psi+1; i++){
    psi[i]=i*step_val[0]+start_val[0];
    psi_r[i]=psi[i]*deg_to_rad;
    tem_coord[0]=psi[i];
    tem_coo_r[0]=psi_r[i];
    for(int j=0; j<div_theta+1; j++){
      theta[j]=j*step_val[1]+start_val[1];
      theta_r[j]=theta[j]*deg_to_rad;
      tem_coord[1]=theta[j];
      tem_coo_r[1]=theta_r[j];
      for(int k=0; k<div_phi+1; k++){
	phi[k]=k*step_val[2]+start_val[2];
	phi_r[k]=phi[k]*deg_to_rad;
	tem_coord[2]=phi[k];
	tem_coo_r[2]=phi_r[k];
	for(int l=0; l<div_x+1; l++){
	  x[l]=l*step_val[3]+start_val[3];
	  tem_coord[3]=x[l];
	  tem_coo_r[3]=x[l];
	  for(int m=0; m<div_y+1; m++){
	    y[m]=m*step_val[4]+start_val[4];
	    tem_coord[4]=y[m];
	    tem_coo_r[4]=y[m];
	    for(int n=0; n<div_z+1; n++){
	      z[n]=n*step_val[5]+start_val[5];
	      tem_coord[5]=z[n];
	      tem_coo_r[5]=z[n];
	      rho_tem=rho(tem_coo_r);
	      if(rho_max<rho_tem){
		rho_max=rho_tem;
		for(int p=0; p<6; p++){
		  max_coo_r[p]=tem_coo_r[p];
		  max_coord[p]=tem_coord[p];
		}
	      }
	      if(check_val_max<rho_tem){
		v_rho_max=rho_tem;
		v_psi_max=tem_coord[0];
		v_theta_max=tem_coord[1];
		v_phi_max=tem_coord[2];
		v_x_max=tem_coord[3];
		v_y_max=tem_coord[4];
		v_z_max=tem_coord[5];
		calmax->Fill();
	      }
	      if(rho_min>rho_tem){
		rho_min=rho_tem;
		for(int q=0; q<6; q++){
		  min_coo_r[q]=tem_coo_r[q];
		  min_coord[q]=tem_coord[q];
		}
	      }
	      if(check_val_min>rho_tem){
		v_rho_min=rho_tem;
		v_psi_min=tem_coord[0];
		v_theta_min=tem_coord[1];
		v_phi_min=tem_coord[2];
		v_x_min=tem_coord[3];
		v_y_min=tem_coord[4];
		v_z_min=tem_coord[5];
		calmin->Fill();
	      }
	    }
	  }
	}
        #ifdef check_phi
	if(k%100==0){
	  cout<<"calculating ... psi "<<100*i/div_psi<<"%, theta "<<100*j/div_theta<<"%, phi "
	      <<100*j/div_theta<<"% compleated now     \r"<<flush;
	}
	if(i==div_psi && j==div_theta && k==div_phi){
	  cout<<"calculating ... psi 100%, theta 100%, phi 100% compleated now      "<<endl;
	}
        #endif
      }
      #ifdef check_theta
      if(j%100==0){
	cout<<"calculating ... psi "<<100*i/div_psi<<"%, theta "<<100*j/div_theta<<"% compleated now    \r"<<flush;
      }
      if(i==div_psi && j==div_theta){
	cout<<"calculating ... psi 100%, theta 100% compleated now      "<<endl;
      }
      #endif
    }
    #ifdef check_psi
    if(i%10==0){
      cout<<"calculating ... psi "<<100*i/div_psi<<"% compleated now    \r"<<flush;
    }
    if(i==div_psi){
      cout<<"calculating ... psi 100% compleated now      "<<endl;
    }
    #endif
  }
  /*
  //val check
  cout<<"start val check "<<psi[0]<<"  "<<theta[0]<<"  "<<phi[0]<<"  "
      <<x[0]<<"  "<<y[0]<<"  "<<z[0]<<endl;
  cout<<"     val check "<<psi[1]<<"  "<<theta[1]<<"  "<<phi[1]<<"  "
      <<x[1]<<"  "<<y[1]<<"  "<<z[1]<<endl;
  cout<<"end val check "<<psi[div_psi]<<"  "<<theta[div_theta]<<"  "<<phi[div_phi]<<"  "
      <<x[div_x]<<"  "<<y[div_y]<<"  "<<z[div_z]<<endl;
  */
  //out put
  cout<<endl;
  cout<<"min point"<<endl;
  cout<<" rho: "<<rho_min<<endl;
  cout<<"   (psi, theta, phi, x, y, z)=("
      <<min_coord[0]<<", "<<min_coord[1]<<", "<<min_coord[2]<<", "
      <<min_coord[3]<<", "<<min_coord[4]<<", "<<min_coord[5]<<")"<<endl;

  cout<<"max point"<<endl;
  cout<<" rho: "<<rho_max<<endl;
  cout<<"   (psi, theta, phi, x, y, z)=("
      <<max_coord[0]<<", "<<max_coord[1]<<", "<<max_coord[2]<<", "
      <<max_coord[3]<<", "<<max_coord[4]<<", "<<max_coord[5]<<")"<<endl<<endl;

  /*
  //check value
  Double_t test1[6]={max_coo_r[0],max_coo_r[1],max_coo_r[2],280,280,1000};
  Double_t rho_test1=rho(test1);
  Double_t test2[6]={max_coo_r[0],max_coo_r[1],max_coo_r[2],-280,-280,0};
  Double_t rho_test2=rho(test2);
  Double_t test3[6]={min_coord[0],min_coord[1],min_coord[2],-280,-280,0};
  Double_t rho_test3=rho(test3);
  Double_t test4[6]={min_coord[0],min_coord[1],min_coord[2],280,280,1000};
  Double_t rho_test4=rho(test4);

  cout<<"test1 point (psi, theta, phi, x, y, z)=("
      <<test1[0]*180/pi<<", "<<test1[1]*180/pi<<", "<<test1[2]*180/pi<<", "
      <<test1[3]<<", "<<test1[4]<<", "<<test1[5]<<")"<<endl;
  cout<<"rho_test1 "<<rho_test1<<endl<<endl;

  cout<<"test2 point (psi, theta, phi, x, y, z)=("
      <<test2[0]*180/pi<<", "<<test2[1]*180/pi<<", "<<test2[2]*180/pi<<", "
      <<test2[3]<<", "<<test2[4]<<", "<<test2[5]<<")"<<endl;
  cout<<"rho_test2 "<<rho_test2<<endl<<endl;

  cout<<"test3 point (psi, theta, phi, x, y, z)=("
      <<test3[0]*180/pi<<", "<<test3[1]*180/pi<<", "<<test3[2]*180/pi<<", "
      <<test3[3]<<", "<<test3[4]<<", "<<test3[5]<<")"<<endl;
  cout<<"rho_test3 "<<rho_test3<<endl<<endl;

  cout<<"test4 point (psi, theta, phi, x, y, z)=("
      <<test4[0]*180/pi<<", "<<test4[1]*180/pi<<", "<<test4[2]*180/pi<<", "
      <<test4[3]<<", "<<test4[4]<<", "<<test4[5]<<")"<<endl;
  cout<<"rho_test4 "<<rho_test4<<endl<<endl;
  */

  //save branch
  calmax->AutoSave();
  calmin->AutoSave();
  calfile->Close();
  cout<<"output: "<<outdata<<endl;

  // stop timer
  timer.Stop();
  Double_t rtime=timer.RealTime();
  Double_t ctime=timer.CpuTime();
  cout<<endl;
  cout<<"Real time "<<rtime<<" s, CPU time "<<ctime<<" s"<<endl<<endl;
  
  return 0;
}
