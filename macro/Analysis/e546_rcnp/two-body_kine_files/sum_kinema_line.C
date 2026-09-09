//#define debug_mode
#include <string>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TString.h"
#include "TObject.h"
#include "TCanvas.h"
#include "Math/Point3D.h"

TGraph* ReadKinematics(TString kineFile, TString ReturnType);

void sum_kinema_line(){

   TFile * Results = new TFile("../data2/c12_12c12c/sum_kine_lines_E69.5.root","recreate");

   // Kinematic lines.
   TGraph *kine_gsgs_0   = ReadKinematics("gsgs_69.5_verz_0",   "kine");
   TGraph *kine_gsgs_50  = ReadKinematics("gsgs_66.9_verz_50",  "kine");
   TGraph *kine_gsgs_100 = ReadKinematics("gsgs_64.2_verz_100", "kine");
   TGraph *kine_gsgs_150 = ReadKinematics("gsgs_61.5_verz_150", "kine");
   TGraph *kine_gsgs_200 = ReadKinematics("gsgs_58.6_verz_200", "kine");
   TGraph *kine_gsgs_250 = ReadKinematics("gsgs_55.7_verz_250", "kine");
   TGraph *kine_gsgs_300 = ReadKinematics("gsgs_52.6_verz_300", "kine");
   TGraph *kine_gsgs_350 = ReadKinematics("gsgs_49.4_verz_350", "kine");
   TGraph *kine_gsgs_400 = ReadKinematics("gsgs_46.0_verz_400", "kine");
   TGraph *kine_gsgs_450 = ReadKinematics("gsgs_42.5_verz_450", "kine");
   TGraph *kine_gsgs_500 = ReadKinematics("gsgs_38.8_verz_500", "kine");
   TGraph *kine_gsgs_550 = ReadKinematics("gsgs_34.8_verz_550", "kine");
   TGraph *kine_gsgs_600 = ReadKinematics("gsgs_30.5_verz_600", "kine");
   TGraph *kine_gsgs_650 = ReadKinematics("gsgs_25.9_verz_650", "kine");
   TGraph *kine_gsgs_700 = ReadKinematics("gsgs_20.7_verz_700", "kine");
   TGraph *kine_gsgs_750 = ReadKinematics("gsgs_14.8_verz_750", "kine");
   TGraph *kine_gsgs_800 = ReadKinematics( "gsgs_7.7_verz_800", "kine");

   TGraph *kine_gsex_0   = ReadKinematics("gsex_69.5_verz_0",   "kine");
   TGraph *kine_gsex_50  = ReadKinematics("gsex_66.9_verz_50",  "kine");
   TGraph *kine_gsex_100 = ReadKinematics("gsex_64.2_verz_100", "kine");
   TGraph *kine_gsex_150 = ReadKinematics("gsex_61.5_verz_150", "kine");
   TGraph *kine_gsex_200 = ReadKinematics("gsex_58.6_verz_200", "kine");
   TGraph *kine_gsex_250 = ReadKinematics("gsex_55.7_verz_250", "kine");
   TGraph *kine_gsex_300 = ReadKinematics("gsex_52.6_verz_300", "kine");
   TGraph *kine_gsex_350 = ReadKinematics("gsex_49.4_verz_350", "kine");
   TGraph *kine_gsex_400 = ReadKinematics("gsex_46.0_verz_400", "kine");
   TGraph *kine_gsex_450 = ReadKinematics("gsex_42.5_verz_450", "kine");
   TGraph *kine_gsex_500 = ReadKinematics("gsex_38.8_verz_500", "kine");
   TGraph *kine_gsex_550 = ReadKinematics("gsex_34.8_verz_550", "kine");
   TGraph *kine_gsex_600 = ReadKinematics("gsex_30.5_verz_600", "kine");
   TGraph *kine_gsex_650 = ReadKinematics("gsex_25.9_verz_650", "kine");
   TGraph *kine_gsex_700 = ReadKinematics("gsex_20.7_verz_700", "kine");
   TGraph *kine_gsex_750 = ReadKinematics("gsex_14.8_verz_750", "kine");

   TGraph *kine_exex_0   = ReadKinematics("exex_69.5_verz_0",   "kine");
   TGraph *kine_exex_50  = ReadKinematics("exex_66.9_verz_50",  "kine");
   TGraph *kine_exex_100 = ReadKinematics("exex_64.2_verz_100", "kine");
   TGraph *kine_exex_150 = ReadKinematics("exex_61.5_verz_150", "kine");
   TGraph *kine_exex_200 = ReadKinematics("exex_58.6_verz_200", "kine");
   TGraph *kine_exex_250 = ReadKinematics("exex_55.7_verz_250", "kine");
   TGraph *kine_exex_300 = ReadKinematics("exex_52.6_verz_300", "kine");
   TGraph *kine_exex_350 = ReadKinematics("exex_49.4_verz_350", "kine");
   TGraph *kine_exex_400 = ReadKinematics("exex_46.0_verz_400", "kine");
   TGraph *kine_exex_450 = ReadKinematics("exex_42.5_verz_450", "kine");
   TGraph *kine_exex_500 = ReadKinematics("exex_38.8_verz_500", "kine");
   TGraph *kine_exex_550 = ReadKinematics("exex_34.8_verz_550", "kine");
   TGraph *kine_exex_600 = ReadKinematics("exex_30.5_verz_600", "kine");
   TGraph *kine_exex_650 = ReadKinematics("exex_25.9_verz_650", "kine");
   TGraph *kine_exex_700 = ReadKinematics("exex_20.7_verz_700", "kine");

   TGraph *angle_gsgs_0 = ReadKinematics("gsgs_69.5_verz_0", "angle");
   TGraph *angle_gsex_0 = ReadKinematics("gsex_69.5_verz_0", "angle");
   TGraph *angle_exex_0 = ReadKinematics("exex_69.5_verz_0", "angle");

   // save lines
   Results->cd();
   kine_gsgs_0  ->Write("kine_12c_gsgs_z0"  );
   kine_gsgs_50 ->Write("kine_12c_gsgs_z50" );
   kine_gsgs_100->Write("kine_12c_gsgs_z100");
   kine_gsgs_150->Write("kine_12c_gsgs_z150");
   kine_gsgs_200->Write("kine_12c_gsgs_z200");
   kine_gsgs_250->Write("kine_12c_gsgs_z250");
   kine_gsgs_300->Write("kine_12c_gsgs_z300");
   kine_gsgs_350->Write("kine_12c_gsgs_z350");
   kine_gsgs_400->Write("kine_12c_gsgs_z400");
   kine_gsgs_450->Write("kine_12c_gsgs_z450");
   kine_gsgs_500->Write("kine_12c_gsgs_z500");
   kine_gsgs_550->Write("kine_12c_gsgs_z550");
   kine_gsgs_600->Write("kine_12c_gsgs_z600");
   kine_gsgs_650->Write("kine_12c_gsgs_z650");
   kine_gsgs_700->Write("kine_12c_gsgs_z700");
   kine_gsgs_750->Write("kine_12c_gsgs_z750");
   kine_gsgs_800->Write("kine_12c_gsgs_z800");

   kine_gsex_0  ->Write( "kine_12c_gsex_z0" );
   kine_gsex_50 ->Write("kine_12c_gsex_z50" );
   kine_gsex_100->Write("kine_12c_gsex_z100");
   kine_gsex_150->Write("kine_12c_gsex_z150");
   kine_gsex_200->Write("kine_12c_gsex_z200");
   kine_gsex_250->Write("kine_12c_gsex_z250");
   kine_gsex_300->Write("kine_12c_gsex_z300");
   kine_gsex_350->Write("kine_12c_gsex_z350");
   kine_gsex_400->Write("kine_12c_gsex_z400");
   kine_gsex_450->Write("kine_12c_gsex_z450");
   kine_gsex_500->Write("kine_12c_gsex_z500");
   kine_gsex_550->Write("kine_12c_gsex_z550");
   kine_gsex_600->Write("kine_12c_gsex_z600");
   kine_gsex_650->Write("kine_12c_gsex_z650");
   kine_gsex_700->Write("kine_12c_gsex_z700");
   kine_gsex_750->Write("kine_12c_gsex_z750");
   //   kine_gsex_800->Write("kine_12c_gsex_z800");

   kine_exex_0  ->Write( "kine_12c_exex_z0" );
   kine_exex_50 ->Write("kine_12c_exex_z50" );
   kine_exex_100->Write("kine_12c_exex_z100");
   kine_exex_150->Write("kine_12c_exex_z150");
   kine_exex_200->Write("kine_12c_exex_z200");
   kine_exex_250->Write("kine_12c_exex_z250");
   kine_exex_300->Write("kine_12c_exex_z300");
   kine_exex_350->Write("kine_12c_exex_z350");
   kine_exex_400->Write("kine_12c_exex_z400");
   kine_exex_450->Write("kine_12c_exex_z450");
   kine_exex_500->Write("kine_12c_exex_z500");
   kine_exex_550->Write("kine_12c_exex_z550");
   kine_exex_600->Write("kine_12c_exex_z600");
   kine_exex_650->Write("kine_12c_exex_z650");
   kine_exex_700->Write("kine_12c_exex_z700");
   //   kine_exex_750->Write("kine_12c_exex_z750");
   //   kine_exex_800->Write("kine_12c_exex_z800");

   angle_gsgs_0->Write("angle_12c_gsgs_z0");
   angle_gsex_0->Write("angle_12c_gsex_z0");
   angle_exex_0->Write("angle_12c_exex_z0");

   Results->Close();
   gSystem->Exit(0);

}

TGraph* ReadKinematics(TString kineFileName, TString ReturnType){

   TString kineder = "./out_kinema_69.5/";
   TString kineFile = kineder + "kine_12c12c_" + kineFileName + ".txt";
   //   std::vector<Double_t> ThetacmS, ThetalabR, ElabR, ThetalabS, ElabS, MomlabR;
   std::vector<Double_t> Thetacm, Theta3lab, E3lab, Theta4lab, E4lab;
#ifdef debug_mode
   // Debug: Check if file exists and current working directory
   std::cout << "DEBUG ReadKinematics: pwd = " << gSystem->pwd() << std::endl;
   std::cout << "DEBUG ReadKinematics: Attempting to read: " << kineFile << std::endl;
   if (gSystem->AccessPathName(kineFile.Data(), kFileExists)) {
      std::cout << "DEBUG ReadKinematics: FILE NOT FOUND!" << std::endl;
   } else {
      std::cout << "DEBUG ReadKinematics: File found" << std::endl;
   }
#endif

   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr -> is_open() || !kineStr -> good()){
      std::cout << "Warning : No kinematics file found or cannot read: " << kineFile << std::endl;
      return new TGraph();
   }
   /*
   while (!kineStr->eof()){
   //	      *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >> EnerLabSca[numKin];
	   *kineStr >> ThetaLabRec[numKin] >> EnerLabRec[numKin];
      numKin++;
      }
   */
   double tcm, t3lab, e3lab, t4lab, e4lab;
   string line;
   while (getline(*kineStr, line)){
      if(line.empty() || line[0] == '#'){
         continue;
      }
      istringstream iss(line);
      if (iss >> tcm >> t3lab >> e3lab >> t4lab >> e4lab){
         Thetacm.push_back(tcm);
         Theta3lab.push_back(t3lab);
         E3lab.push_back(e3lab);
         Theta4lab.push_back(t4lab);
         E4lab.push_back(e4lab);
      }
   }

   /*
   while ( *kineStr >> tlabr >> elabr){
      ThetalabR.push_back(tlabr);
      ElabR.push_back(elabr);
   }
   */
#ifdef debug_mode
      std::cout << "DEBUG ReadKinematics: Successfully read " << Theta3lab.size() << " points" << std::endl;
#endif
   if (Theta3lab.size() == 0){
      std::cout << " Warning : No data read from kinematics file: " << kineFile << std::endl;
   }
   //   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   TGraph *kine  = new TGraph(Theta3lab.size(), Theta3lab.data(), E3lab.data());
   TGraph *angle = new TGraph(Theta3lab.size(), Theta3lab.data(), Theta4lab.data());
   //   kine->Sort();

   if (ReturnType == "kine"){
      return kine;
   }
   else if (ReturnType == "angle"){
      return angle;
   }
   else {
      std::cout << "Warning: Invalid ReturnType specified. Returning empty TGraph." << std::endl;
      return new TGraph();
   }

}
