#include "TMath.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

static double sq(double x){ return x*x; }

void kinema_12c12c(const char *states = "gsgs"){

   TString outdi="./out_kinema_69.5/";
   TString oufi="";
   Int_t n_data =0;
   Int_t n_file =0;
   Int_t verz =0;
   Double_t n=1800;
   Double_t m = 11174.863254; // MeV/c^2; 12C mass
   Double_t Ex3 = 0.0;
   Double_t Ex4 = 0.0;
   Double_t K1_lab = 0.0;
   Double_t dtheta= TMath::Pi()/n;
   //  cout<<"dtheta:"<<dtheta<<endl;

   Double_t m1 = m;
   Double_t m2 = m;
   Double_t m3 = m;
   Double_t m4 = m;

   if (strcmp(states, "gsgs") == 0){
      Ex3 = 0.0;
      Ex4 = 0.0;
      std::cout << "You chose reaction: gsgs; Ex3 = " << Ex3 << " MeV, Ex4 = " << Ex4 << " MeV" << std::endl;
   }
   else if (strcmp(states, "gsex") == 0){
      Ex3 = 0.0;
      Ex4 = 4.44;
      std::cout << "You chose reaction: gsex; Ex3 = " << Ex3 << " MeV, Ex4 = " << Ex4 << " MeV" << std::endl;
   }
   else if (strcmp(states, "exex") == 0){
      Ex3 = 4.44;
      Ex4 = 4.44;
      std::cout << "You chose reaction: exex; Ex3 = " << Ex3 << " MeV, Ex4 = " << Ex4 << " MeV" << std::endl;
   }
   else{
      std::cerr << "Please input the correct states: gsgs, gsex, exex" << std::endl;
      gSystem->Exit(1);
   }

   //calculate mass include Ex
   Double_t m3p = m3 + Ex3;
   Double_t m4p = m4 + Ex4;

   /*
   //  69.0 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 69.027}, {50, 66.411}, {100, 63.727}, {150, 60.949}, {200, 58.083}, {250, 55.105}, {300, 52.061},
         {350, 48.784},{400, 45.411},{450, 41.844}, {500, 38.084}, {550, 34.062},{600, 29.719}, {650, 24.984},
         {700, 19.708}, {750, 13.621}, {800, 6.194}, {830, 0.832}, {835, 0.271}
   };
   */

   //  69.5 MeV injection energy to ATTPC
   vector<pair<Double_t, Double_t>> lise_data={
         {0, 69.500}, {50, 66.899}, {100, 64.228}, {150, 61.468}, {200, 58.618}, {250, 55.662}, {300, 52.594},
         {350, 49.391},{400, 46.046},{450, 42.514}, {500, 38.797}, {550, 34.822},{600, 30.545}, {650, 25.892},
         {700, 20.728}, {750, 14.821}, {800, 7.705}, {840, 0.723}, {848, 0.028}
   };

   Int_t l = lise_data.size();
   std::vector<std::vector<Double_t>> Ebeam(2, std::vector<Double_t>(l, 0));
   for (Int_t i = 0; i < l; i++){
      Ebeam.at(0).at(i) = lise_data[i].first;
      Ebeam.at(1).at(i) = lise_data[i].second;
      if (Ebeam.at(1).at(i) > 1.0){
         n_data ++;
      }
   }

   for (Int_t i = 0; i < l; i++){
      if (Ebeam.at(1).at(i) < 1.0 ){
         continue;
      }
      else {
         K1_lab = Ebeam.at(1).at(i);
         verz = int(Ebeam.at(0).at(i));
      }
      //starting condition
      Double_t E1_lab = K1_lab + m1;
      Double_t p1_lab = sqrt(max(0.0, E1_lab*E1_lab - m1*m1));
      Double_t s = m1*m1 + m2*m2 + 2.0*m2*E1_lab; //Lorentz invariant
      Double_t sqrt_s = sqrt(s);

      if(s < sq(m3p+m4p)){
         std::cout << std::endl;
         std::cout << "Depth:" << verz << " mm, reaction forbidden: sqrt(s) < m3'+m4'" << std::endl;
         //         return;
         break;
      }
      oufi = outdi + Form("kine_12c12c_%.4s_%.1f_verz_%d.txt", states, K1_lab, verz);
      n_file ++;

      Double_t term1 = s - sq(m3p + m4p);
      Double_t term2 = s - sq(m3p - m4p);
      Double_t p_cm = 0.5 * sqrt(max(0.0, term1*term2)) / sqrt_s;

      Double_t E3_cm = sqrt(sq(p_cm) + sq(m3p));
      Double_t E4_cm = sqrt(sq(p_cm) + sq(m4p));

      //CM to Lab
      Double_t P_lab = p1_lab;
      Double_t E_lab_total = E1_lab + m2;
      Double_t beta = P_lab / E_lab_total;
      Double_t gamma = 1.0 / sqrt(1.0 - beta*beta);

      //set output
      std::ofstream fout(oufi);
      if(!fout){
         std::cout << "Cannot open output file." << std::endl;
         return;
      }

      //calculation
      fout << "#kinematics calculation of " << states << "\n";
      fout << "#vertex z = " << verz << " mm, E_reaction = " << K1_lab << " MeV" << "\n";
      fout << "#theta_cm, theta3_lab, Ek3_lab, theta4_lab, Ek4_lab" << "\n";
      for(int k = 0; k < n+1;k++){
         if(k%100 == 0){
            std::cout << "Finish " << verz << "mm: " << 100*k/n << " % !      \r" << std::flush;
         }

         Double_t theta_cm_rad = k * dtheta;
         Double_t theta_cm_deg = theta_cm_rad * 180.0 / TMath::Pi();
         Double_t ct = cos(theta_cm_rad);
         Double_t st = sin(theta_cm_rad);
         //    std::cout << "i:" << i << ", theta:"<<theta_cm/TMath::Pi() << std::endl;

         //particle 3
         Double_t E3_lab = gamma*(E3_cm + beta*p_cm*ct);
         Double_t p3_para = gamma*(p_cm*ct + beta*E3_cm);
         Double_t p3_perp = p_cm*st;
         //    Double_t theta3_lab = atan2(p3_perp, p3_para);
         Double_t theta3_lab = atan2(p3_perp, p3_para);
         Double_t theta3_deg = theta3_lab * 180.0 / TMath::Pi();
         //    if(theta3_deg < 0) theta3_deg += 180.0;

         Double_t Ek3_lab = E3_lab - m3p;
         if(Ek3_lab < 0) Ek3_lab = 0;

         //particle 4  (cm angle to pi-theta_cm)
         Double_t ct4 = -ct;
         Double_t st4 =  st;
         Double_t E4_lab = gamma*(E4_cm + beta*p_cm*ct4);
         Double_t p4_para = gamma*(p_cm*ct4 + beta*E4_cm);
         Double_t p4_perp = p_cm*st4;
         //    Double_t theta4_lab = atan2(p4_perp,p4_para);
         Double_t theta4_lab = atan2(p4_perp, p4_para);
         Double_t theta4_deg = theta4_lab * 180.0 / TMath::Pi();
         //    if(theta4_deg < 0) theta4_deg += 180.0;

         Double_t Ek4_lab = E4_lab - m4p;
         if(Ek4_lab < 0) Ek4_lab = 0;

         //output (cm angle, theta3, Ek3, theta4, Ek4)
         fout << theta_cm_deg << " " << theta3_deg << " " << Ek3_lab << " " << theta4_deg << " " << Ek4_lab << "\n";
      }
      fout.close();
   }
   if (n_data != n_file){
      std::cout << std::endl;
      std::cout << "Warning!! Data and file mismatch for some reason." << std::endl;
      std::cout << "  n_data: " << n_data << ", n_file: " << n_file << std::endl;
   }

   std::cout << std::setprecision(3) << std::endl;
   std::cout << "Finish calculation, E: " << Ebeam.at(1).at(0) << " MeV, Depth: " << Ebeam.at(0).at(n_file - 1) << "mm " << std::endl;
   std::cout << "Output saved. n_file: " << n_file << std::endl;
   gSystem->Exit(0);

}
