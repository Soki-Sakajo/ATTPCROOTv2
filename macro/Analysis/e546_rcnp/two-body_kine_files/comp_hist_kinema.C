void hist_kinema(const char* filename = "output.txt"){
   //read file (theta_cm, theta3_lab, E3. theta4_lab, E4)

   TString dire = "./out_file/";
   //  TString infi=dire+filename;
   TString infi = "./out_kinema_69/kine_12c12c_exex_69.0_verz_0.txt";
   TString infi2 = "./Einjection69/angle_12c12c_exex_4.44_69.0.txt";
   TString infi3 = "./Einjection69/kine_12c12c_exex_4.44_69.0.txt";
   ifstream fin(infi);
   ifstream fin2(infi2);
   ifstream fin3(infi3);
   if(!fin || !fin2 || !fin3){
   //      std::cerr << "Input file not found: " << filename << std::endl;
      std::cerr << "Some input files not found" << std::endl;
      if(!fin){
         std::cerr << "file1 not found" << std::endl;
      }
      if(!fin2){
         std::cerr << "file2 not found" << std::endl;
      }
      if(!fin3){
         std::cerr << "file3 not found" << std::endl;
      }
      return;
   }

   //hist difinition
   TH2D* hthetaE_m3 =
      new TH2D("hThetaE_m3","m3 Scattering;#theta_{lab} (deg);Energy (MeV)",180,0,180,500,0,50);

   TH2D* hthetaE_m4 =
      new TH2D("hThetaE_m4","m4 Scattering;#theta_{lab} (deg);Energy (MeV)",180,0,180,500,0,50);

   TH2D* hthetacmtheta3 =
      new TH2D("hthetacmtheta","cm-m4;#theta_{lab}(deg);m4 #theta_{lab}(deg)",180,0,180,360,0,180);

   TH2D* htheta3theta4 =
      new TH2D("htheta3theta4","m3-m4;#theta_{lab}(deg);m4 #theta_{lab}(deg)",180,0,180,360,0,180);

   TH2D* hthetaE_m3_2 =
      new TH2D("hThetaE_m3_2","m3 Scattering_from_site;#theta_{lab} (deg);Energy (MeV)",180,0,180,500,0,50);

   TH2D* htheta3theta4_2 =
      new TH2D("htheta3theta4_2","m3-m4_from_site;#theta_{lab}(deg);m4 #theta_{lab}(deg)",180,0,180,360,0,180);

   std::string line;
   Int_t i=0;
   Int_t n_line=0;
   Double_t th_cm, th3, E3, th4, E4;
   vector<vector<Double_t>> info(5,vector<Double_t>(0));
   vector<vector<Double_t>> info2(2,vector<Double_t>(0));
   vector<vector<Double_t>> info3(2,vector<Double_t>(0));

   //read file
   // line count
   while(getline(fin,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      n_line++;
   }
   cout<<"Number of lines: "<<n_line<<endl;
   fin.clear();
   fin.seekg(0,std::ios::beg);
   // read value
   while(getline(fin,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      cout<<"Finish "<<100*i/n_line<<" %!       \r"<<flush;
      //    cout<<i<<" %!       \r"<<flush;
      i++;
      istringstream iss(line);
      iss >> th_cm >> th3 >> E3 >> th4 >> E4;
      //    cout<<th_cm<<" "<<th3<<" "<<E3<<" "<<th4<<" "<<E4<<endl;
      info[0].push_back(th_cm);
      info[1].push_back(th3);
      info[2].push_back(E3);
      info[3].push_back(th4);
      info[4].push_back(E4);
      hthetaE_m3->Fill(th3,E3);
      hthetaE_m4->Fill(th4,E4);
      htheta3theta4->Fill(th3,th4);
      hthetacmtheta3->Fill(th_cm,th4);
   }

   while(getline(fin2,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      n_line++;
   }
   cout<<"Number of lines: "<<n_line<<endl;
   fin2.clear();
   fin2.seekg(0,std::ios::beg);
   // read value
   while(getline(fin2,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      cout<<"Finish "<<100*i/n_line<<" %!       \r"<<flush;
      //    cout<<i<<" %!       \r"<<flush;
      i++;
      istringstream iss(line);
      iss >> th3 >> th4;
      //    cout<<th_cm<<" "<<th3<<" "<<E3<<" "<<th4<<" "<<E4<<endl;
      info2[0].push_back(th3);
      info2[1].push_back(th4);
      htheta3theta4_2->Fill(th3,th4);
   }

   while(getline(fin3,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      n_line++;
   }
   cout<<"Number of lines: "<<n_line<<endl;
   fin3.clear();
   fin3.seekg(0,std::ios::beg);
   // read value
   while(getline(fin3,line)){
      if(line.empty()||line[0]=='#'){
         continue;
      }
      cout<<"Finish "<<100*i/n_line<<" %!       \r"<<flush;
      //    cout<<i<<" %!       \r"<<flush;
      i++;
      istringstream iss(line);
      iss >> th3 >> E3;
      //    cout<<th_cm<<" "<<th3<<" "<<E3<<" "<<th4<<" "<<E4<<endl;
      info3[0].push_back(th3);
      info3[1].push_back(E3);
      hthetaE_m3_2->Fill(th3,E3);
   }

   TGraph *theta_cm=new TGraph(info.at(0).size(),info.at(0).data(),info.at(1).data());
   TGraph *kine3=new TGraph(info.at(1).size(),info.at(1).data(),info.at(2).data());
   TGraph *theta=new TGraph(info.at(1).size(),info.at(1).data(),info.at(3).data());
   TGraph *kine4=new TGraph(info.at(3).size(),info.at(3).data(),info.at(4).data());
   TGraph *theta_2=new TGraph(info2.at(0).size(),info2.at(0).data(),info2.at(1).data());
   TGraph *kine3_2=new TGraph(info3.at(0).size(),info3.at(0).data(),info3.at(1).data());

   fin.close();
   cout<<"Finish 100 %!       "<<endl;
   //write
   TCanvas* can = new TCanvas("can","can");
   can->Divide(2,2);
   //  can->Divide(3,3);
   can->cd(1);
   theta_cm->SetTitle("cm vs lab;#theta_{3,cm} (deg);#theta_{3,lab} (deg)");
   theta_cm->Draw("AP");
   can->cd(2);
   kine3->SetTitle("m3 Kinetic Energy;#theta_{3,lab} (deg);Energy (MeV)");
   kine3->Draw("AP");
   can->cd(3);
   kine4->SetTitle("m4 Kinetic Energy;#theta_{4,lab} (deg);Energy (MeV)");
   kine4->Draw("AP");
   can->cd(4);
   theta->SetTitle("m3 vs m4;#theta_{3,lab} (deg);#theta_{4,lab} (deg)");
   theta->Draw("AP");
   can->Update();

   TCanvas* c1 = new TCanvas("c1","c1");
   c1->Divide(2,2);
   //  can->Divide(3,3);
   c1->cd(1);
   kine3->SetTitle("m3 Kinetic Energy;#theta_{3,lab} (deg);Energy (MeV)");
   kine3->Draw("AP");
   c1->cd(2);
   kine3_2->SetTitle("m3 Kinetic Energy from site;#theta_{3,lab} (deg);Energy (MeV)");
   kine3_2->Draw("AP");
   c1->cd(3);
   theta->SetTitle("m3 vs m4;#theta_{3,lab} (deg);#theta_{4,lab} (deg)");
   theta->Draw("AP");
   c1->cd(4);
   theta_2->SetTitle("m3 vs m4 from site;#theta_{3,lab} (deg);#theta_{4,lab} (deg)");
   theta_2->Draw("AP");
   c1->Update();
}
