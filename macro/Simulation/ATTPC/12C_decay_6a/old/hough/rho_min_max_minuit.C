#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

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

void rho_val(Int_t &npar, Double_t *gin, Double_t &fval, Double_t *vars, Int_t iflag){
  fval=rho(vars);
}
void rho_max_val(Int_t &npar, Double_t *gin, Double_t &fval, Double_t *vars, Int_t iflag){
  fval=-rho(vars);
}

void rho_min_max_minuit(){
  const Double_t pi=3.141592653589;
  Double_t ang_step=0.1;
  Double_t xy_step=1;
  Double_t z_step=1;
  Double_t rho_max;
  Double_t rho_min;
  Double_t max_coord[6]={0,0,0,0,0,0};
  Double_t max_coord_err[6]={0,0,0,0,0,0};
  Double_t min_coord[6]={0,0,0,0,0,0};
  Double_t min_coord_err[6]={0,0,0,0,0,0};
  Double_t a_min;
  Double_t b_min;
  Double_t c_min;
  Double_t a_min_err;
  Double_t b_min_err;
  Double_t c_min_err;
  Double_t a_max;
  Double_t b_max;
  Double_t c_max;
  Double_t a_max_err;
  Double_t b_max_err;
  Double_t c_max_err;
  /*
  Double_t test1[6]={0,0,0,0,0,0};
  Double_t rho_test1;
  Double_t test2[6]={0,0,0,0,0,0};
  Double_t rho_test2;
  */
  
  TMinuit minuit(6);
  minuit.SetFCN(rho_val);

  //set condition
  string name[6]={"var1","var2","var3","var4","var5","var6"};
  Double_t start_val[6]={pi,pi/4,pi,0,0,500};
  //  Double_t start_val[6]={0,0,0,0,0,0};
  Double_t step[6]={ang_step,ang_step,ang_step,xy_step,xy_step,z_step};
  Double_t min[6]={0,0,0,-280,-280,0};
  Double_t max[6]={2*pi,pi/2,2*pi,280,280,1000};

  for(int i=0; i<6; i++){
    minuit.DefineParameter(i,name[i].c_str(),start_val[i],step[i],min[i],max[i]);
  }

  //min calculation
  minuit.Migrad();
  for(int i=0; i<6; i++){
    minuit.GetParameter(i, min_coord[i], min_coord_err[i]);
  }
  rho_min=rho(min_coord);
  a_min=min_coord[0]*180/pi;
  b_min=min_coord[1]*180/pi;
  c_min=min_coord[2]*180/pi;
  a_min_err=min_coord_err[0]*180/pi;
  b_min_err=min_coord_err[1]*180/pi;
  c_min_err=min_coord_err[2]*180/pi;

  /*
  //max calculation
  minuit.SetFCN(rho_max_val);
  minuit.Migrad();
  for(int i=0; i<6; i++){
    minuit.GetParameter(i, max_coord[i], max_coord_err[i]);
  }
  rho_max=rho(max_coord);
  a_max=max_coord[0]*180/pi;
  b_max=max_coord[1]*180/pi;
  c_max=max_coord[2]*180/pi;
  a_max_err=max_coord_err[0]*180/pi;
  b_max_err=max_coord_err[1]*180/pi;
  c_max_err=max_coord_err[2]*180/pi;
  */
  
  //out put
  cout<<endl<<endl;
  cout<<"min_point (psi, theta, phi, x, y, z)=("
      <<a_min<<"± "<<a_min_err<<", "<<b_min<<"± "<<b_min_err<<", "<<c_min<<"± "<<c_min_err<<", "
      <<min_coord[3]<<", "<<min_coord[4]<<", "<<min_coord[5]<<")"<<endl;
  cout<<"rho_min "<<rho_min<<endl<<endl;
  /*
  cout<<"max_point (psi, theta, phi, x, y, z)=("
      <<a_max<<"± "<<a_max_err<<", "<<b_max<<"± "<<b_max_err<<", "<<c_max<<"± "<<c_max_err<<", "
      <<max_coord[3]<<", "<<max_coord[4]<<", "<<max_coord[5]<<")"<<endl;
  cout<<"rho_max "<<rho_max<<endl<<endl<<endl;

  //check value
  Double_t test1[6]={max_coord[0],max_coord[1],max_coord[2],280,280,1000};
  Double_t rho_test1=rho(test1);
  Double_t test2[6]={max_coord[0],max_coord[1],max_coord[2],-280,-280,0};
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

  return 0;
}
