{
gStyle->SetPadGridX(1);
gStyle->SetPadGridY(1);
}

#include "userdef.h++"

#define N_LOOP 10000000
#define OFFSET_LOOP 1
#define N_EVE 1000000
//#define N_EVE 100000

#define N_FADC_BOARD 8
#define N_CH 16
#define N_MADC_CH 64

struct detector {
  double angle;
  double pos[3][3]; // position of corner (left,right,top)
  double vec[3][3]; // vec reac_pos -> det_corner
  double n_vec[4][3]; // n_vec for face of detector_pyramid
  double plane[4]; // n_vec [a,b,c,d] for detector plane ax+by+cz+d=0
};

struct detector_ROOT{
  TVector3 pos[3];
  TVector3 normal;
  Double_t d;
};

struct func_return{
  Int_t all_event;
  Int_t count[7];
  Int_t count_6a[7];
};

struct hit_part{
  Int_t parent,chF,chR,PID;
  Double_t E;
  int hit;
};

struct hit_seg{
  vector<hit_part> part;
  map<int,int> hit_F,hit_R;
};

int ch2seg(int ch){
  if(ch<0 || 40<=ch){
    return -1;
  }
  int seg;
  if(ch<8) seg=3; //Si4
  else if(ch<16) seg=4; //Si5
  else if(ch<24) seg=1; //Si2
  else if(ch<32) seg=0; //Si1
  else seg=2; //Si3
  return seg;
}

#include "veclib.c"

const Double_t AMU = 0.93149432;      // atomic mass unit [GeV]
const Double_t C = 299792458;       // light speed [m/s]
const Double_t M_alpha = 4*AMU+2.4249/1000;     // alpha (= 4He) mass [GeV]
const Double_t M_4He = M_alpha;     // alpha (= 4He) mass [GeV]
const Double_t M_12C = 12*AMU;      // 12C mass [GeV]
const Double_t M_8Be = 8*AMU+4.9416/1000; // 8Be mass [GeV]
const Double_t M_16O = 16*AMU-4.7370/1000; // 16O mass [GeV]
const Double_t weight_max[3] = {3, 3, 3};
const Double_t Ex_3a = 0.00765; // 12C(0+2) excitation energy [GeV]
const Double_t Ex_12C[]={0.00/1000, 4.44/1000, 7.65/1000, 9.64/1000, 10000000};
const Double_t Ex_4a = 0.0151; // 16O(0+6) excitation energy [GeV]

const Double_t Density_C=2.253; // [mg/cm^3]

#include "target_effect.cpp"

TCanvas *c_ene_ang;
TCanvas *c_det_xy;
TCanvas *c_det_xz;
TCanvas *c_ene;
TCanvas *c_dist;
TCanvas *c_count;
TCanvas *c_ene_3a;
TCanvas *c_pol;
TCanvas *c_track;

void h_Write(double _E_beam=-1, int i_mon=-1, int i_ene=-1, int mode=-1){
  string root_dir="root_file";
  string of_name;
  if(_E_beam==-1) of_name=Form("%s/simu_%s_%02d_%d_hist.root",root_dir.c_str(),s_Mon[i_mon].c_str(),i_ene,mode);
  else of_name=Form("%s/simu_%2.1f_hist.root",root_dir.c_str(),_E_beam);

  auto objs=gDirectory->GetList();
  unsigned int nObjs=objs->GetSize();

  TFile *fout=new TFile(of_name.c_str(),"recreate");
  for(int i=0;i<nObjs;++i) objs->At(i)->Write();
  fout->Close();
  
  std::cout << Form("histgram written in %s",of_name.c_str()) << std::endl;

  return;
}

//func_return phase_space_2hoyle(Double_t _E_beam, Double_t offset=0){
func_return phase_space_2hoyle(int i_mon, int i_ene, int mode){
  double _E_beam=ene_list[i_mon][i_ene];
  double offset=0;
  //string fname=Form("Multi6_%d_simu.gif",37);
  //if(gROOT->IsBatch()) gSystem->Unlink(fname.c_str());
  gErrorIgnoreLevel=kWarning;
  /* Canvas */
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  //gStyle->SetOptLogz(1);
  if(c_ene_ang!=NULL){
    c_ene_ang->Close();
    c_det_xy->Close();
    c_det_xz->Close();
    c_ene->Close();
    c_dist->Close();
    c_count->Close();
    c_ene_3a->Close();
    c_pol->Close();
    c_track->Close();
  }
  c_ene_ang = new TCanvas("c_ene_ang","Alpha ene and ang",1400,500);
  c_ene_ang->Divide(2,1);
  c_det_xy = new TCanvas("c_det_xy","x-y @ det",700,500);
  c_det_xz = new TCanvas("c_det_xz","x-z @ det",700,500);
  c_ene = new TCanvas("c_ene","Kinetic energy 12C+12C",700,500);
  c_ene->Divide(3,2);
  c_dist = new TCanvas("c_dist","Angler Diatribution",700,500);
  c_dist->Divide(1,2);
  c_count = new TCanvas("c_count","Hit pattern",700,500);
  c_count->Divide(sqrt(OFFSET_LOOP),sqrt(OFFSET_LOOP));
  c_ene_3a=new TCanvas("c_ene_3a","c_ene_3a",700,500);
  c_pol=new TCanvas("c_pol","c_pol",600,600);
  gPad->DrawFrame(-150., -150., 150., 150.);
  c_track=new TCanvas("c_track","c_track",600,600);

  /* Histograms */
  auto objs=gDirectory->GetList();
  unsigned int nObjs=objs->GetSize();
  while(objs->GetSize()) objs->At(0)->Delete();

  TH1D *h_ene[6];
  for(int i=0; i<6; i++)
    h_ene[i]= new TH1D(Form("h_ene_%d",i), Form("particle%d",i), 100, 0.0, 30.0);
  TH1D *h_ene_6a=new TH1D("h_ene_6a","h_ene_6a" , 100, 0.0, 30.0);
  TH1D *h_ene_3a[2];
  //for(int i=0;i<2;++i) h_ene_3a[i]=new TH1D(Form("h_ene_3a_%d",i),Form("h_ene_3a %d",i), 100, 0.0, 20.0);
  for(int i=0;i<2;++i) h_ene_3a[i]=new TH1D(Form("h_ene_3a_%d",i),Form("h_ene_3a %d",i), 360, 0,360.0);
  h_ene_3a[1]->SetLineColor(kRed);
  TH1D *h_inv2C=new TH1D("h_inv2C","inv2C",500,0,30);
  TH2D *h_th_phi = new TH2D("h_th_phi", "h2 title", 90, 0., 180., 90, -180, 180);
  TH1D *h_th = new TH1D("h_th", "h2 title", 90, 0., 180.);
  TH2D *h_xy = new TH2D("h_xy", "x-y", 100, -1.0, 1.0, 100, -1.0, 1.0);
  TH2D *h_xz = new TH2D("h_xz", "x-z", 100, -1.0, 1.0, 100, -1.0, 1.0);
  TH1D *h_ang[5];
  for(int i=0;i<5;i++)
    h_ang[i]=new TH1D(Form("h_ang_%d",i), Form("particle%d_1",i), 360, 0.0, 360.);
  TH2D *h_K_theta[2];
  for(int i=0;i<2;++i) h_K_theta[i]=new TH2D(Form("h_K_theta_%d",i),Form("h_K_theta_%d",i),200,0,50,200,0,50);
  h_K_theta[0]->SetTitle("12C K vs theta, 3a on same Si;K[MeV];Theat[deg]");
  h_K_theta[1]->SetTitle("12C K vs theta, 3a;K[MeV];Theat[deg]");
  TH2D *h_det_xy=new TH2D("h_det_xy","x-y @ det ",200, -250.0, 250.0, 200, -250.0, 250.0);
  TH2D *h_det_xz=new TH2D("h_det_xz","x-z @ det ",200, -250.0, 250.0, 200, -50.0-offset, 200.0-offset);
  TH1D *h_r = new TH1D("h_r","hit R", 1000, 0.0, 250.);
  TH1D *h_theta = new TH1D("h_theta","hit theta", 900, 0.0, 90.);
  TH1D *h_count = new TH1D("h_count","hit count",7,0,7);
  TH1D *h_hit_3a=new TH1D("h_hit_3a","h_hit_3a",5,0,5);
  TH1D *h_TDC=new TH1D("h_TDC","h_TDC",500,0,2500);
  TH2D *h_Esum_FR[5];
  for(int seg=0;seg<5;++seg) h_Esum_FR[seg]=new TH2D(Form("h_Esum_FR_%d",seg),Form("h_Esum_FR_%d",seg),200,0,30,200,0,30);
  TH1D *h_angle[2];
  for(int i=0;i<2;++i) h_angle[i]=new TH1D(Form("h_angle_%d",i),Form("h_angle_%d",i),360,0,360);
  h_angle[1]->SetLineColor(kRed);
  TH1D *h_angle_min=new TH1D("h_angle_min","h_angle_min",3,-1,2);
  TH1D *h_NnotAssigned=new TH1D("h_NnotAssigned","h_NnotAssigned",5,0,5);
  TH1D *h1_Ex=new TH1D("h1_Ex","h1_Ex",200,5,20);
  TH2D *h2_Ex[2],*h2_Ex_id[2];
  TH2D *h2_theta1_vs_theta2[2];
  TH2D *h2_12C_theta_theta[2];
  TH2D *h2_E1_vs_E2[2];
  TH1D *h1_delta_phi[2],*h1_delta_theta[2],*h1_delta_E[2],*h1_angle_cm[2];
  for(int i=0;i<2;++i){
    h2_Ex[i]=new TH2D(Form("h2_Ex_%d",i),Form("h2_Ex_%d",i),110,-1,10,110,-1,10);
    h2_Ex_id[i]=new TH2D(Form("h2_Ex_id_%d",i),Form("h2_Ex_id_%d",i),5,0,5,5,0,5);
    h2_theta1_vs_theta2[i]=new TH2D(Form("h2_theta1_vs_theta2_%d",i),Form("h2_theta1_vs_theta2_%d",i),360,0,180,360,0,180);
    h2_E1_vs_E2[i]=new TH2D(Form("h2_E1_vs_E2_%d",i),Form("h2_E1_vs_E2_%d",i),360,0,40,360,0,40);
    h1_delta_phi[i]=new TH1D(Form("h1_delta_phi_%d",i),Form("h1_delta_phi_%d",i),3600,-360,360);
    h1_delta_phi[i]->SetLineColor(colors[i]);
    h1_delta_theta[i]=new TH1D(Form("h1_delta_theta_%d",i),Form("h1_delta_theta_%d",i),3600,-180,180);
    h1_delta_theta[i]->SetLineColor(colors[i]);
    h1_delta_E[i]=new TH1D(Form("h1_delta_E_%d",i),Form("h1_delta_E_%d",i),5000,-50,50);
    h1_delta_E[i]->SetLineColor(colors[i]);
    h1_angle_cm[i]=new TH1D(Form("h1_angle_cm_%d",i),Form("h1_angle_cm_%d",i),360,0,200);
    h1_angle_cm[i]->SetLineColor(colors[i]);
    h2_12C_theta_theta[i]=new TH2D(Form("h2_12C_theta_theta_%d",i),Form("h2_12C_theta_theta_%d;#theta_{cm};#theta_{cm}",i),180,0,180,180,0,180);
  }
  TH1D *h1_theta_12C[2][2][4]; //[all, in one seg][7.65, 9.64][gs,4.44,7.65,9.64]
  TH1D *h1_theta_12C_ratio[2][4];
  for(int i=0;i<2;++i){
    for(int j=0;j<2;++j){
      for(int k=0;k<4;++k){
	h1_theta_12C[i][j][k]=new TH1D(Form("h1_theta_12C_%d_%d_%d",i,j,k),Form("h1_theta_12C_%d_%d_%d",i,j,k),90,0,180);
	h1_theta_12C[i][j][k]->SetLineColor(colors[i]);
	if(i==0) h1_theta_12C_ratio[j][k]=new TH1D(Form("h1_theta_12C_ratio_%d_%d",j,k),Form("h1_theta_12C_ratio_%d_%d",j,k),90,0,180);
      }
    }
  }
  TH1D *h1_Esum[6];
  for(int i=0;i<6;++i) h1_Esum[i]=new TH1D(Form("h1_Esum_%d",i),Form("h1_Esum_%d",i),1000,-50,50);
  
  TH3D *h3_p=new TH3D("h3_p","h3_p",50,-1,1,50,-1,1,50,-1,1);
  TH1D *h1_p=new TH1D("h1_p","h1_p",500,-0.5,0.5);
  
  TPolyLine *det_pline_xz[5];
  TPolyLine *det_pline_xy[5];
  TH3D *h3=new TH3D("h3",";x;z;y",200,-150,150,200,0,150,200,-150,150);
  h3->SetStats(0);
  h3->GetXaxis()->SetTitleOffset(2);
  h3->GetYaxis()->SetTitleOffset(2);
  h3->GetZaxis()->SetTitleOffset(1.5);

  TH1D *h1_E_beam=new TH1D("h1_E_beam","h1_E_beam",500,30,60);
  TH2D *h2_reac_pos=new TH2D("h2_reac_pos","h2_reac_pos",100,-0.01,0.01,100,-0.01,0.01);
  TH1D *h1_3aInSeg_resi=new TH1D("h1_3aInSeg_resi","h1_3aInSeg_resi",500,-15,25);
  TH2D *h2_theta_E[5][5];
  for(int i=0;i<5;++i){
    for(int j=0;j<5;++j){
      h2_theta_E[i][j]=new TH2D(Form("h2_theta_E_%d_%d",i,j),Form("h2_theta_E_%d_%d;Ecm;theta_cm",i,j),200,0,50,180,0,180);
    }
  }
  TH2D *h2_theta_phi_incl[2];
  for(int i=0;i<2;++i) h2_theta_phi_incl[i]=new TH2D(Form("h2_theta_phi_incl_%d",i),Form("h2_theta_phi_incl_%d",i),180,-180,180,180,0,180);
  
  TH2D *polar_hist=new TH2D("polar_hist","polar_hist",
			    3*360/60*8, -2*TMath::Pi(), 4*TMath::Pi(),
			    N_RING, R_INNER, R_OUTER);
  // To draw the frame of the histogram,
  // Fill one event to each bin
  for(int board=0;board<5;++board){
    for(int i=0;i<16;++i){
      for(int j=0;j<8;++j){
	double r = R_INNER + (R_PITCH/2.0) + (double)i*R_PITCH;
	double theta = (60/2.0-60.0/8/2.0)/180.0*TMath::Pi()
	  - (double)j*60.0/8/180.0*TMath::Pi()
	  + (double)(board*Si_phi+Si_phi_off)/180*TMath::Pi();
	
	polar_hist->Fill(theta, r);
      }
    }
  }

  c_pol->cd();
  polar_hist->Draw("same colz pol");
  TText *text_pol=new TText(-150,155,Form("E_beam=%.1f MeV",_E_beam));
  text_pol->Draw();
  TText *t_si[5];
  t_si[0]=new TText(-110,110,"Si1");
  t_si[1]=new TText(-140,-100,"Si2");
  t_si[2]=new TText(70,-140,"Si3");
  t_si[3]=new TText(110,-100,"Si4");
  t_si[4]=new TText(90,110,"Si5");
  for(int i=0;i<5;++i) t_si[i]->Draw();

  // vvvvv Tree vvvvv //
  //string tree_file=Form("root_file/simu_%.1f.root",_E_beam);
  string tree_file=Form("root_file/simu_%s_%02d_%d.root",s_Mon[i_mon].c_str(),i_ene,mode);
  TFile *tfile=new TFile(tree_file.c_str(),"recreate");
  TTree *tree=new TTree("tree","tree");
  int entry;
  int Ex_count[4][4]={}; // gs, 4.44, 7.65, 9.64  
  int ADC[N_FADC_BOARD][N_CH];
  int PeakClk[N_FADC_BOARD][N_CH];
  int sum[N_FADC_BOARD][N_CH];
  int baseline[N_FADC_BOARD][N_CH];
  Double_t Energy[N_FADC_BOARD][N_CH];
  Double_t Energy_MADC[N_MADC_CH];
  vector<vector<int>> *TDC=new vector<vector<int>>(128);
  int Ans_Nhit,Ans_Npix,Ans_Nfront,Ans_Nrear;
  vector<int> *Ans_parent=new vector<int>();
  vector<int> *Ans_chF=new vector<int>();
  vector<int> *Ans_chR=new vector<int>();
  vector<int> *Ans_hit=new vector<int>();
  vector<double> *Ans_Energy=new vector<double>();
  vector<int> *Ans_PID=new vector<int>();
  tree->Branch("entry",&entry,"entry/I");
  tree->Branch("Ex_count",Ex_count,"Ex_count[4][4]/I");
  tree->Branch("ADC",ADC,Form("ADC[%d][%d]/I",N_FADC_BOARD,N_CH));
  tree->Branch("PeakClk",PeakClk,Form("PeakClk[%d][%d]/I",N_FADC_BOARD,N_CH));
  tree->Branch("sum",sum,Form("sum[%d][%d]/I",N_FADC_BOARD,N_CH));
  tree->Branch("baseline",baseline,Form("baseline[%d][%d]/I",N_FADC_BOARD,N_CH));
  tree->Branch("Energy",Energy,Form("Energy[%d][%d]/D",N_FADC_BOARD,N_CH));
  tree->Branch("Energy_MADC",Energy_MADC,Form("Energy_MADC[%d]/D",N_MADC_CH));
  tree->Branch("TDC","std::vector<std::vector<int>>",TDC);
  tree->Branch("Ans_Nhit",&Ans_Nhit,"Ans_Nhit/I");
  tree->Branch("Ans_Npix",&Ans_Npix,"Ans_Npix/I");
  tree->Branch("Ans_Nfront",&Ans_Nfront,"Ans_Nfront/I");
  tree->Branch("Ans_Nrear",&Ans_Nrear,"Ans_Nrear/I");
  tree->Branch("Ans_parent","std::vector<int>",Ans_parent);
  tree->Branch("Ans_chF","std::vector<int>",Ans_chF);
  tree->Branch("Ans_chR","std::vector<int>",Ans_chR);
  tree->Branch("Ans_hit","std::vector<int>",Ans_hit);
  tree->Branch("Ans_Energy","std::vector<double>",Ans_Energy);
  tree->Branch("Ans_PID","std::vector<int>",Ans_PID);
  // ^^^^^ Tree ^^^^^ //

  TH2D *h2_EMADC_eff;
  TFile *ifile_sum=new TFile(Form("root_file%d/inv/inv_vs_run_%s.root",i_mon,s_Mon[i_mon].c_str()));
  gDirectory->cd("Rint:/");
  h2_EMADC_eff=(TH2D*)((TH2D*)ifile_sum->Get(Form("h2_EMADC_ch_ene_%d_%d_1",i_mon,i_ene)))->Clone();
  h2_EMADC_eff->SetName("h2_EMADC_eff");
  TH2D h2_all=*(TH2D*)ifile_sum->Get(Form("h2_EMADC_ch_ene_%d_%d_0",i_mon,i_ene));
  h2_EMADC_eff->Divide(&h2_all);
  ifile_sum->Close();

  double total_cs[2]={};
  double integ_cs[2][4]={};
  int pointX=-1;
  double angle_cs_max[2][4]={};
  TGraphErrors *g_angle_cs_ans_excl[N_MON][N_ENERGY][N_Ex_INDEX][N_Ex_INDEX];
  vector<bitset<4>> use_state(4);
  double anlge_low[N_Ex_INDEX][N_Ex_INDEX]={},angle_high[N_Ex_INDEX][N_Ex_INDEX]={};
  TGraphMultiErrors *g_cs_ans_integ_excl[N_MON][N_Ex_INDEX][N_Ex_INDEX];
  TH1D *h1_3aInSeg_angle_cs_ans_excl[N_MON][N_ENERGY][N_Ex_INDEX][N_Ex_INDEX];
  TFile *ifile_excl=new TFile("root_file_res/angle_dist_excl.root");
  TFile *ifile_excl_total=new TFile("root_file_res/excl_total.root");
  TGraphMultiErrors *g_cs_ans_excl[N_MON][N_Ex_INDEX][N_Ex_INDEX];
  gDirectory->cd("Rint:/");
  for(int i_inv=I_Ex_765;i_inv<=I_Ex_964;++i_inv){
    for(int i_mis=I_Ex_000;i_mis<=I_Ex_964;++i_mis){
      h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]=(TH1D*)ifile_excl->Get(Form("h1_3aInSeg_angle_cs_ans_excl_%d_%d_%d_%d",i_mon,i_ene,i_inv,i_mis));
      int tmp_bin=h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]->GetMaximumBin();
      double tmp_theta=h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]->GetBinCenter(tmp_bin);
      double bin_sin=sin(tmp_theta*D_TO_R);
      double dtheta=180./h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]->GetNbinsX()*D_TO_R;
      angle_cs_max[i_inv-I_Ex_765][i_mis]=h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]->GetMaximum()*2*M_PI*bin_sin*dtheta;
      if(angle_cs_max[i_inv-I_Ex_765][i_mis]!=0) use_state[i_inv].set(i_mis);
      else use_state[i_inv].reset(i_mis);      
      //g_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]=(TGraphErrors*)ifile_excl->Get(Form("g_angle_cs_ans_excl_%d_%d_%d_%d",i_mon,i_ene,i_inv,i_mis));
      g_angle_cs_ans_excl[i_mon][i_ene][i_inv][i_mis]=(TGraphErrors*)ifile_excl->Get(Form("g_angle_cs_ans_Ngate_excl_%d_%d_%d_%d",i_mon,i_ene,i_inv,i_mis));
      g_cs_ans_integ_excl[i_mon][i_inv][i_mis]=(TGraphMultiErrors*)ifile_excl->Get(Form("g_cs_ans_integ_excl_%d_%d_%d",i_mon,i_inv,i_mis));
      g_cs_ans_excl[i_mon][i_inv][i_mis]=(TGraphMultiErrors*)ifile_excl_total->Get(Form("g_cs_ans_excl_%d_%d_%d",i_mon,i_inv,i_mis));
      if(pointX==-1){
	for(int i=0;i<16;++i){
	  if(_E_beam==g_cs_ans_excl[i_mon][i_inv][i_mis]->GetPointX(i)*2.){
	    pointX=i;
	    break;
	  }
	}
      }
      total_cs[i_inv-I_Ex_765]+=g_cs_ans_excl[i_mon][i_inv][i_mis]->GetPointY(pointX);
      integ_cs[i_inv-I_Ex_765][i_mis]=total_cs[i_inv-I_Ex_765];
    }
  }
  for(int i=0;i<2;++i){
    for(int j=0;j<4;++j){
      integ_cs[i][j]/=total_cs[i];
    }
  }

  // vvvvv detector setting vvvvv //
  TPolyLine3D *det_line[5];
  struct detector_ROOT det[5];
  for(int seg=0;seg<5;++seg){
    det_line[seg]=new TPolyLine3D();
    TVector3 det_pos[3];
    for(int i=0;i<3;++i){
      det[seg].pos[i].SetXYZ(0,R_OUTER,0);      
      det[seg].pos[i].RotateZ(pow(-1,i)*TH_PITCH*4/180.*TMath::Pi());
      det[seg].pos[2].SetXYZ(0,0,0);
      det[seg].pos[i].SetY(det[seg].pos[i].Y()+Si_offset);
      det[seg].pos[i].RotateX(-(90-Si_theta)/180.*TMath::Pi());
      det[seg].pos[i].RotateZ(Si_phi*(seg+0.5)/180.*TMath::Pi());
      det[seg].pos[i].SetZ(det[seg].pos[i].Z()+SAKRA_offset+offset);
    }    
    form_3d_plane(&det[seg]);
    det_pline_xy[seg]=new TPolyLine();
    det_pline_xz[seg]=new TPolyLine();
    for(int i=0;i<4;++i){
      det_pline_xy[seg]->SetNextPoint(det[seg].pos[i%3].X(),det[seg].pos[i%3].Y());
      det_pline_xz[seg]->SetNextPoint(det[seg].pos[i%3].X(),det[seg].pos[i%3].Z());
      det_line[seg]->SetNextPoint(det[seg].pos[i%3].X(),det[seg].pos[i%3].Z(),det[seg].pos[i%3].Y());
    }    
  }
  // ^^^^^ detector setting ^^^^^ //
  
  /* Physics setting */
  Double_t T_beam = _E_beam/1000;        // beam (12C) kinetic energy [GeV]
  Double_t E_beam = M_12C + T_beam;  // beam (12C) total energy [GeV]
  Double_t Pz_beam = TMath::Sqrt(E_beam*E_beam - M_12C*M_12C);  // beam (12C) Pz [GeV]
  TLorentzVector beam(0.0, 0.0, Pz_beam, E_beam);  // 4-mom. vector for the beam
  TLorentzVector target(0.0, 0.0, 0.0, M_12C);     // 4-mom. vector for the target
  TLorentzVector W = beam + target;
  TVector3 reac_pos(0,0,0);
  double E_reac=_E_beam-LISE_Eloss_12C(_E_beam,50);
  
  /* event generator setting */
  TGenPhaseSpace event[6];
  TGenPhaseSpace event_8BeAlpha[2],event_2alpha[2];
  // 0:12C(hoyle) + 12C(hoyle) 1,2:12C(hoyle)->3-alpha
  Double_t mass_12c12c[2] = {M_12C + Ex_3a, M_12C + Ex_3a};
  Double_t mass_8BeAlpha[2] ={M_8Be, M_alpha};
  Double_t mass_2alpha[2] = {M_alpha, M_alpha};
  Double_t mass_3alpha[3] = {M_alpha, M_alpha, M_alpha};
  event[0].SetDecay(W, 2, mass_12c12c);
  TLorentzVector *pAlpha[6];

  Double_t E_thre=1.0;
  Int_t hit_det_all=0;
  Int_t count_6a[7]={};

  string fnam_track=Form("fig/SAKRA_beam%.1f_off%.1f_track.pdf",T_beam*1000,offset);
  //c_track->Print((fnam_track+"[").c_str());
  // vvvvv event generation vvvvv //
  Double_t weight, uniform_ran;
  Int_t n=0;
  cout << Form("start E=%.1f MeV",_E_beam) << endl;
  cout << "total    event : " << setw(10) << N_EVE << endl;
  for (n=0; n<N_LOOP; n++) {
    if(tree->GetEntries()==N_EVE) break;
    if(n%(N_EVE/100)==0){
      cout << "\rGenerate event : " << setw(10) << n;
      cout << flush;
    }
    /* create 12C+12C->24Mg */
    
    double coll_depth=gRandom->Uniform(0,1)*100; // [ug/cm^2]
    TLorentzVector beam_straggled;
    create24Mg(_E_beam,coll_depth,&beam_straggled,&reac_pos);
    h1_E_beam->Fill((beam_straggled.E()-M_12C)*1000);
    h2_reac_pos->Fill(reac_pos[0],reac_pos[1]);
    W=beam_straggled+target;
    
    /* event generation 12C(hoyle) + 12C(hoyle) */
    int Ex_index[2];
    do {      
      //Ex_index[0]=2; // 12C-1 : 7.65
      //Ex_index[0]=(int)gRandom->Uniform(2,4); // 12C-1 : 7.65 or 9.64
      double tmp=gRandom->Uniform(0,1);
      if(tmp<0.3) Ex_index[0]=2;
      else Ex_index[0]=3;
      if(mode==2 || mode==3){
	if(!use_state[3].any()) Ex_index[0]=2;
      }
      else if(mode==5) Ex_index[0]=0;
      mass_12c12c[0]=M_12C+Ex_12C[Ex_index[0]];
      if(mode==1){
	//Ex_index[1]=(int)gRandom->Uniform(0,5); // 12C-2 : gs or 4.44 or 7.65 or 9.64
	double tmp=gRandom->Uniform(0,1);
	double tmp_ratio[4]={0.1,0.2,0.4,1};
	for(int i=0;i<4;++i){
	  if(tmp<tmp_ratio[i]){
	    Ex_index[1]=i;
	    break;
	  }	  
	}
      }
      else if(mode==2){
	while(1){
	  int tmp=(int)gRandom->Uniform(0,4);
	  if(use_state[Ex_index[0]][tmp]){
	    Ex_index[1]=tmp;
	    break;
	  }
	}
      }
      else if(mode==3){
	double tmp=gRandom->Uniform(0,1);
	for(int i=0;i<4;++i){
	  if(tmp<integ_cs[Ex_index[0]-2][i]){
	    Ex_index[1]=i;
	    break;
	  }
	}
      }
      else if(mode==4){
	Ex_index[1]=0;
      }
      else if(mode==5){
	Ex_index[1]=0;
      }
      mass_12c12c[1]=M_12C+Ex_12C[Ex_index[1]];
    }
    while(!event[0].SetDecay(W, 2, mass_12c12c));
    //while(!event[0].SetDecay(W, 2, mass_12c12c) || Ex_index[0]<Ex_index[1]);
    Ex_count[Ex_index[0]][Ex_index[1]]++;
  reGene:
    do {
      weight = event[0].Generate();
      if (weight > weight_max[0]) {
	std::cerr << "weight_max[0] is too small!" << std::endl;
	exit(0);
      }
      uniform_ran = gRandom->Uniform(0.,weight_max[0]);
    } while (uniform_ran > weight);
    if(mode==2){
      TLorentzVector c1=*event[0].GetDecay(0);
      TVector3 tmp_boost=W.BoostVector();
      c1.Boost(-tmp_boost);
      double tmp_theta=c1.Theta()*R_TO_D;
      double bin_sin=sin(tmp_theta*D_TO_R);
      double dtheta=180./h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][Ex_index[0]][Ex_index[1]]->GetNbinsX()*D_TO_R;
      int tmp_bin=h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][Ex_index[0]][Ex_index[1]]->FindBin(tmp_theta);
      
      if(h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][Ex_index[0]][Ex_index[1]]->GetBinContent(tmp_bin)==0)
	goto reGene;
            
      //double tmp_cs=h1_3aInSeg_angle_cs_ans_excl[i_mon][i_ene][Ex_index[0]][Ex_index[1]]->GetBinContent(tmp_bin)*2*M_PI*bin_sin*dtheta;
      double tmp_cs=g_angle_cs_ans_excl[i_mon][i_ene][Ex_index[0]][Ex_index[1]]->Eval(tmp_theta,0,"S")*2*M_PI*bin_sin*dtheta;
      if(tmp_cs<gRandom->Uniform(0.,angle_cs_max[Ex_index[0]-I_Ex_765][Ex_index[1]])){
	goto reGene;
      }
    }
    
    if(mode==5){
      TLorentzVector c1=*event[0].GetDecay(0);
      TVector3 tmp_boost=W.BoostVector();
      c1.Boost(-tmp_boost);
      double tmp_theta=c1.Theta()*R_TO_D;
      if(tmp_theta<85 || 95<=tmp_theta) goto reGene;
    }
    
    TLorentzVector C12_1 = *event[0].GetDecay(0);
    TLorentzVector C12_2 = *event[0].GetDecay(1);
    h2_theta1_vs_theta2[0]->Fill(C12_1.Theta()*R_TO_D,C12_2.Theta()*R_TO_D);
    h2_E1_vs_E2[0]->Fill((C12_1.E()-M_12C)*1000,(C12_2.E()-M_12C)*1000);
    h1_Ex->Fill((C12_1.M()-M_12C)*1000);
    h1_Ex->Fill((C12_2.M()-M_12C)*1000);
    h2_Ex[0]->Fill((C12_1.M()-M_12C)*1000,(C12_2.M()-M_12C)*1000);
    h2_Ex_id[0]->Fill(Ex_index[0],Ex_index[1]);
    TVector3 tmp_boost=W.BoostVector();
    TLorentzVector c1_cm=C12_1,c2_cm=C12_2;
    c1_cm.Boost(-tmp_boost); c2_cm.Boost(-tmp_boost);
    h1_theta_12C[0][Ex_index[0]-2][Ex_index[1]]->Fill(c1_cm.Theta()*R_TO_D);
    h2_theta_E[Ex_index[0]][Ex_index[1]]->Fill((c1_cm.E()-M_12C)*1000,c1_cm.Theta()*R_TO_D);
    h2_theta_E[Ex_index[0]][4]->Fill((c1_cm.E()-M_12C)*1000,c1_cm.Theta()*R_TO_D);
    h2_theta_E[4][Ex_index[1]]->Fill((c1_cm.E()-M_12C)*1000,c1_cm.Theta()*R_TO_D);
    if(Ex_index[1]>1){
      h1_theta_12C[0][Ex_index[0]-2][Ex_index[1]]->Fill(c2_cm.Theta()*R_TO_D);      
      //h2_theta_E[Ex_index[0]][4]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);
      //h2_theta_E[4][Ex_index[1]]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);      
    }
    if(Ex_index[0]==Ex_index[1]){
      h2_theta_E[Ex_index[0]][Ex_index[1]]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);
      h2_theta_E[Ex_index[0]][4]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);
    }
    if(Ex_index[0]==3 && Ex_index[1]==2){
      //h2_theta_E[Ex_index[1]][Ex_index[0]]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);
      //h2_theta_E[Ex_index[1]][4]->Fill((c2_cm.E()-M_12C)*1000,c2_cm.Theta()*R_TO_D);
    }
    if(Ex_index[0]==2 && Ex_index[1]==2){
      h1_p->Fill(c1_cm.P()); h1_p->Fill(c2_cm.P());
      h3_p->Fill(c1_cm.Px(),c1_cm.Py(),c1_cm.Pz());
    }
    for(int i=0;i<2;++i){
      if(Ex_index[0]==i+2) h2_theta_phi_incl[i]->Fill(c1_cm.Phi()*R_TO_D,c1_cm.Theta()*R_TO_D);
      if(Ex_index[1]==i+2) h2_theta_phi_incl[i]->Fill(c2_cm.Phi()*R_TO_D,c2_cm.Theta()*R_TO_D);
    }
    if(mode==5){
      h2_12C_theta_theta[0]->Fill(c1_cm.Theta()*R_TO_D,c2_cm.Theta()*R_TO_D);
    }
    /* event generation 12C(hoyle) -> 3-alpha */
    int Ndecay=0;
    if(event_8BeAlpha[0].SetDecay(C12_1, 2, mass_8BeAlpha)){
      do {
	weight = event_8BeAlpha[0].Generate();
	if (weight > weight_max[1]) {
	  std::cerr << "weight_max[1] is too small!" << std::endl;
	  exit(0);
	}
	uniform_ran = gRandom->Uniform(0.,weight_max[1]);
      } while (uniform_ran > weight);      
      TLorentzVector tmp_8Be=*event_8BeAlpha[0].GetDecay(0);
      pAlpha[0]=event_8BeAlpha[0].GetDecay(1);
      Ndecay+=1;

      if(event_2alpha[0].SetDecay(tmp_8Be,2,mass_2alpha)){
	do {
	  weight = event_2alpha[0].Generate();
	  if (weight > weight_max[1]) {
	    std::cerr << "weight_max[1] is too small!" << std::endl;
	    exit(0);
	  }
	  uniform_ran = gRandom->Uniform(0.,weight_max[1]);
	} while (uniform_ran > weight);
	pAlpha[1]=event_2alpha[0].GetDecay(0);
	pAlpha[2]=event_2alpha[0].GetDecay(1);
	Ndecay+=2;
      }
    }
    
    if(Ex_index[0]==Ex_index[1]){
      if(event_8BeAlpha[1].SetDecay(C12_2, 2, mass_8BeAlpha)){
	do {
	  weight = event_8BeAlpha[1].Generate();
	  if (weight > weight_max[2]) {
	    std::cerr << "weight_max[2] is too small!" << std::endl;
	    exit(0);
	  }
	  uniform_ran = gRandom->Uniform(0.,weight_max[2]);
	} while (uniform_ran > weight);      
	TLorentzVector tmp_8Be=*event_8BeAlpha[1].GetDecay(0);
	pAlpha[3]=event_8BeAlpha[1].GetDecay(1);
	Ndecay+=1;

	if(event_2alpha[1].SetDecay(tmp_8Be,2,mass_2alpha)){
	  do {
	    weight = event_2alpha[1].Generate();
	    if (weight > weight_max[2]) {
	      std::cerr << "weight_max[2] is too small!" << std::endl;
	      exit(0);
	    }
	    uniform_ran = gRandom->Uniform(0.,weight_max[2]);
	  } while (uniform_ran > weight);
	  pAlpha[4]=event_2alpha[1].GetDecay(0);
	  pAlpha[5]=event_2alpha[1].GetDecay(1);
	  Ndecay+=2;
	}
      }
    }
    // ^^^^^ event generation ^^^^^ //

    // vvvvv initialization vvvvv //
    entry=-1;
    for(int seg=0;seg<N_FADC_BOARD;++seg){
      for(int ch=0;ch<N_CH;++ch){
	ADC[seg][ch]=-1;
	PeakClk[seg][ch]=-1;
	sum[seg][ch]=-1;
	baseline[seg][ch]=-1;
	Energy[seg][ch]=-1;
      }
    }
    for(int ch=0;ch<N_MADC_CH;++ch) Energy_MADC[ch]=-1;
    for(int i=0;i<TDC->size();++i) TDC->at(i).clear();

    Ans_Nhit=0; Ans_Npix=0; Ans_Nfront=0; Ans_Nrear=0;
    Ans_parent->clear();
    Ans_chF->clear();
    Ans_chR->clear();
    Ans_hit->clear();
    Ans_Energy->clear();
    Ans_PID->clear();
    // ^^^^^ initialization ^^^^^ //
    for (int i=0; i<Ndecay; i++) {
      deflect(100-coll_depth,pAlpha[i],reac_pos);
      h_ene[i]->Fill((pAlpha[i]->E() - M_alpha)*1000); // kinetic energy1 [MeV]
      /* particle distribution at r=1 */
      h_xy->Fill(pAlpha[i]->Px()/pAlpha[i]->P(),pAlpha[i]->Py()/pAlpha[i]->P());
      h_xz->Fill(pAlpha[i]->Px()/pAlpha[i]->P(),pAlpha[i]->Pz()/pAlpha[i]->P());
      // angular distribution
      h_th_phi->Fill(pAlpha[i]->Theta()*R_TO_D, pAlpha[i]->Phi()*R_TO_D); 
      h_th->Fill(pAlpha[i]->Theta()*R_TO_D);    
    }
    if(Ndecay==3){
      pAlpha[3]=&C12_2;
      deflect_12C(100-coll_depth,pAlpha[3],reac_pos);
    }
    if(mode==5){
      for(int i=0;i<2;++i){
	pAlpha[i] = event[0].GetDecay(i);
	deflect_12C(100-coll_depth,pAlpha[i],reac_pos);
      }
    }
    Int_t hit_count[5][16][8]={};
    Int_t hit_count_add=0;
    Int_t N_alpha=0; // number of alpha particle which over threshold      
    hit_det_all=0;
    N_alpha=0;
    Int_t flag[6]={-1,-1,-1,-1,-1,-1};
    bitset<6> hit_bit;
    bitset<5*16> hitF_bit;
    bitset<5*8>  hitR_bit;
    vector<int> Nhit_chF(5*16,0),Nhit_chR(5*8,0);
    vector<hit_seg> segment(5);
    vector<TVector3> pos;
    double Esum=0,Esum_f[5]={},Esum_r[5]={};
    int Npart=Ndecay==6?6:4;
    if(mode==5) Npart=2;    
    vector<double> mass_part(Npart,M_alpha);
    if(Npart==4) mass_part[3]=M_12C+Ex_12C[Ex_index[1]];
    if(mode==5){
      mass_part[0]=M_12C; mass_part[1]=M_12C; 
    }
    /* detector hit detection */
    for(int i=0;i<Npart;i++){
      //if((pAlpha[i]->E()-M_alpha)*1000<E_thre || 15<(pAlpha[i]->E()-M_alpha)*1000) continue;
      if((pAlpha[i]->E()-mass_part[i])*1000>=E_thre) N_alpha++;
      if(pAlpha[i]->Theta()>90*D_TO_R) continue;
      Double_t tmp_phi=pAlpha[i]->Phi()*R_TO_D-90;
      if(tmp_phi<0) tmp_phi+=360;
      for(int seg=0;seg<5;seg++){
	if(72*seg<tmp_phi && tmp_phi<72*(seg+1)){
	  TVector3 hit_pos=det_hit_pos(&det[seg],&reac_pos,pAlpha[i]);
	  pos.push_back(hit_pos);
	  Double_t l=(hit_pos-reac_pos).Mag()*1e-3;
	  Double_t TOF=1e9*l/(pAlpha[i]->Beta()*C)*10; // 100ps unit
	  h_TDC->Fill(TOF);
	  TVector3 hit_pos_local=hit_pos-det[seg].pos[2];
	  Double_t hit_r=hit_pos_local.Mag();
	  if(hit_r<R_INNER || R_OUTER<hit_r) continue;
	  Double_t hit_theta=hit_pos_local.Angle(det[seg].pos[0]-det[seg].pos[2])*R_TO_D;
	  Double_t hit_theta2=hit_pos_local.Angle(det[seg].pos[1]-det[seg].pos[2])*R_TO_D;
	  if(hit_theta>TH_COVERAGE || hit_theta2>TH_COVERAGE) continue;
	  h_r->Fill(hit_r);
	  h_theta->Fill(hit_theta);	  
	  for(int r=0;r<16;r++){
	    if(R_OUTER-R_PITCH*(r+1)<hit_r && hit_r<R_OUTER-R_PITCH*r){
	      for(int th=0;th<8;th++){
		if(TH_PITCH*th<hit_theta && hit_theta<TH_PITCH*(th+1)){
		  hit_bit.set(i);
		  int front_ch=16*seg+r;
		  hitF_bit.set(front_ch);
		  int rear_ch;
		  if(seg==0) rear_ch=24+th;
		  else if(seg==1) rear_ch=16+th;
		  else if(seg==2) rear_ch=32+th;
		  else if(seg==3) rear_ch=th;
		  else rear_ch=8+th;
		  hitR_bit.set(rear_ch);
		  // vvvvv data fill for Tree vvvvv //
		  entry=n;
		  PeakClk[seg][r]=500;
		  sum[seg][r]=300e3;
		  ADC[seg][r]=300;
		  baseline[seg][r]=600;
		  double dE=gRandom->Gaus(0,3e-2);
		  //dE=0;
		  if(Energy[seg][r]==-1) Energy[seg][r]=(pAlpha[i]->E()-mass_part[i])*1000+dE;
		  else Energy[seg][r]+=(pAlpha[i]->E()-mass_part[i])*1000+dE;
		  Esum+=(pAlpha[i]->E()-mass_part[i])*1000+dE;
		  dE=gRandom->Gaus(0,4e-2);
		  if(Energy_MADC[rear_ch]==-1) Energy_MADC[rear_ch]=(pAlpha[i]->E()-mass_part[i])*1000+dE;
		  else Energy_MADC[rear_ch]+=(pAlpha[i]->E()-mass_part[i])*1000+dE;
		  double random_eff=gRandom->Uniform(0,1);
		  if(mode==5) random_eff=-1;
		  if(i_mon==0){
		    if(14.8<(pAlpha[i]->E()-mass_part[i])*1000+dE) TDC->at(rear_ch).push_back((int)TOF);
		    else if(random_eff<h2_EMADC_eff->GetBinContent(h2_EMADC_eff->FindBin(rear_ch,(pAlpha[i]->E()-mass_part[i])*1000+dE))){
		      TDC->at(rear_ch).push_back((int)TOF);
		    }
		  }
		  else{
		    if(random_eff<h2_EMADC_eff->GetBinContent(h2_EMADC_eff->FindBin(rear_ch,(pAlpha[i]->E()-mass_part[i])*1000+dE))){
		      TDC->at(rear_ch).push_back((int)TOF);
		    }
		  }
		  if(TDC->at(46).size()==0) TDC->at(46).push_back(120);
		  Ans_Nhit++;
		  hit_part tmp_part={3<=i,front_ch,rear_ch,mass_part[i]==M_alpha?1:2,(pAlpha[i]->E()-mass_part[i])*1000,F3R3};
		  segment[seg].part.push_back(tmp_part);
		  // ^^^^^ data fill for Tree ^^^^^ //
		  if(hit_count[seg][r][th]==0){
		    flag[i]=seg;
		    segment[seg].hit_F[front_ch%16]++;
		    segment[seg].hit_R[rear_ch%8]++;
		    Nhit_chF[front_ch]++;
		    Nhit_chR[rear_ch]++;		  
		  }
		  hit_count[seg][r][th]=1;
		  double tmp_r=R_OUTER-R_PITCH*r-R_PITCH/2.0;
		  double theta=(60.0/8*(4-th)-60.0/8/2+Si_phi*seg+Si_phi_off)*D_TO_R;
		  polar_hist->Fill(theta, tmp_r);		    
		  break;
		}
	      }
	      break;
	    }
	  }
	  h_det_xy->Fill(hit_pos[0],hit_pos[1]);
	  h_det_xz->Fill(hit_pos[0],hit_pos[2]);
	  hit_det_all++;
	  h_ang[seg]->Fill(tmp_phi);
	  break;
	}
      }
    }

    Ans_Nfront=hitF_bit.count();
    Ans_Nrear=hitR_bit.count();

    for(int seg=0;seg<5;++seg){
      int Nhit=segment[seg].part.size();
      int NnotAssigned=Nhit;
      for(int i=0;i<Nhit;++i){
	for(int j=i+1;j<Nhit;++j){
	  if(segment[seg].part[i].chF==segment[seg].part[j].chF && segment[seg].part[i].chR==segment[seg].part[j].chR){
	    if(segment[seg].part[j].hit==F3R3){
	      segment[seg].part[j].hit=SAME;
	      NnotAssigned--;
	    }
	  }
	}
      }
      int Nfront=segment[seg].hit_F.size(),Nrear=segment[seg].hit_R.size();
      int hit=0;
      vector<int> v_hit(Nhit,F3R3);
      if(Nhit==Nfront && Nhit==Nrear){
	for(int i=0;i<Nhit;++i) segment[seg].part[i].hit=FR;
	NnotAssigned=0;
      }
      else{
	for(int i=0;i<Nhit;++i){
	  int chF=segment[seg].part[i].chF,chR=segment[seg].part[i].chR;
	  if(segment[seg].hit_F[chF%16]==1 && segment[seg].hit_R[chR%8]==1){
	    segment[seg].part[i].hit=FR;
	    segment[seg].hit_F.erase(chF%16); segment[seg].hit_R.erase(chR%8);
	    NnotAssigned--;
	  }
	}
      }
      if(NnotAssigned==2){
	if(segment[seg].hit_F.size()==2 && segment[seg].hit_R.size()==1){	  
	  for(int i=0;i<Nhit;++i){
	    if(segment[seg].part[i].hit==F3R3){
	      segment[seg].part[i].hit=FFR;
	      NnotAssigned--;
	    }
	  }
	}
	else if(segment[seg].hit_F.size()==1 && segment[seg].hit_R.size()==2){	  
	  for(int i=0;i<Nhit;++i){
	    if(segment[seg].part[i].hit==F3R3){
	      segment[seg].part[i].hit=FRR;
	      NnotAssigned--;
	    }
	  }
	}
      }
      if(NnotAssigned==3){
	if(segment[seg].hit_F.size()==3 && segment[seg].hit_R.size()==1){	  
	  for(int i=0;i<Nhit;++i){
	    if(segment[seg].part[i].hit==F3R3){
	      segment[seg].part[i].hit=FFFR;
	      NnotAssigned--;
	    }
	  }
	}
	else if(segment[seg].hit_F.size()==1 && segment[seg].hit_R.size()==3){	  
	  for(int i=0;i<Nhit;++i){
	    if(segment[seg].part[i].hit==F3R3){
	      segment[seg].part[i].hit=FRRR;
	      NnotAssigned--;
	    }
	  }
	}
	else if(segment[seg].hit_F.size()==2 && segment[seg].hit_R.size()==2){
	  for(int i=0;i<Nhit;++i){
	    if(segment[seg].part[i].hit==F3R3){
	      segment[seg].part[i].hit=FFRR;
	      NnotAssigned--;
	    }
	  }
	}
      }
      h_NnotAssigned->Fill(NnotAssigned);
      for(int i=0;i<segment[seg].part.size();++i){
	Ans_parent->push_back(segment[seg].part[i].parent);
	Ans_chF->push_back(segment[seg].part[i].chF);
	Ans_chR->push_back(segment[seg].part[i].chR);
	Ans_Energy->push_back(segment[seg].part[i].E);
	Ans_hit->push_back(segment[seg].part[i].hit);
	Ans_PID->push_back(segment[seg].part[i].PID);
      }
    }
    
    if(hit_bit.count()==5){
      hit_bit.flip();
      double ang_sum[2]={0,0};
      int i_off=0;
      int no_hit=hit_bit._Find_first();
      double angle_min=360*4;
      int min_is_same=-1;
      for(int i=0;i<6;++i){
	if(i==no_hit) continue;
	for(int j=i+1;j<6;++j){
	  if(j==no_hit) continue;
	  for(int k=j+1;k<6;++k){
	    if(k==no_hit) continue;
	    int isSame;
	    if(i+j+k==0+1+2 || i+j+k==3+4+5) isSame=1;
	    else isSame=0;
	    ang_sum[isSame]=0;
	    ang_sum[isSame]+=pAlpha[i]->Angle(pAlpha[j]->Vect());
	    ang_sum[isSame]+=pAlpha[j]->Angle(pAlpha[k]->Vect());
	    ang_sum[isSame]+=pAlpha[k]->Angle(pAlpha[i]->Vect());
	    h_angle[isSame]->Fill(ang_sum[isSame]*R_TO_D);	    
	    if(ang_sum[isSame]<angle_min){
	      angle_min=ang_sum[isSame];
	      min_is_same=isSame;
	    }
	  }
	}
      }
      h_angle_min->Fill(min_is_same);
      hit_bit.flip();
    }

    if(hit_bit.count()>0) h1_Esum[hit_bit.count()-1]->Fill(Esum-E_reac);
    if(hit_bit.count()==6){
      //h1_Esum[5]->Fill(Esum-E_reac);
      h2_Ex[1]->Fill((C12_1.M()-M_12C)*1000,(C12_2.M()-M_12C)*1000);
      h2_theta1_vs_theta2[1]->Fill(C12_1.Theta()*R_TO_D,C12_2.Theta()*R_TO_D);
      h2_E1_vs_E2[1]->Fill((C12_1.E()-M_12C)*1000,(C12_2.E()-M_12C)*1000);
      for(int i=1;i<6;++i){
	for(int j=i+1;j<6;++j){
	  TLorentzVector c1=(*pAlpha[0])+(*pAlpha[i])+(*pAlpha[j]),c2;
	  set<int> index_c1={0,i,j};
	  for(int i_c2=1;i_c2<6;++i_c2){
	    if(index_c1.count(i_c2)) continue;
	    c2+=(*pAlpha[i_c2]);
	  }
	  TLorentzVector redi_c=W-c1;
	  bool l=(i==1 && j==2);
	  h1_delta_phi[l]->Fill((redi_c.Phi()-c2.Phi())*180/TMath::Pi());
	  h1_delta_theta[l]->Fill((redi_c.Theta()-c2.Theta())*180/TMath::Pi());
	  h1_delta_E[l]->Fill(redi_c.E()-c2.E());
	  TVector3 b=W.BoostVector();
	  TLorentzVector c1_cm=c1,c2_cm=c2;
	  c1_cm.Boost(-b); c2_cm.Boost(-b);
	  double phi=fabs(c1_cm.Phi()-c2_cm.Phi())*R_TO_D;
	  phi=min(phi,360-phi);
	  h1_angle_cm[l]->Fill(c1_cm.Theta()*R_TO_D);
	  if(l){	    
	    TPolyLine3D *track[6];
	    TPolyMarker3D *hit[6];
	    c_track->Clear();
	    c_track->SetTheta(10);
	    c_track->SetPhi(20);
	    h3->SetTitle(Form("Red:%1.2f MeV Blue:%1.2f MeV;x;z;y",(C12_1.M()-M_12C)*1000,(C12_2.M()-M_12C)*1000));
	    h3->Draw();
	    for(int seg=0;seg<5;++seg) det_line[seg]->Draw("same");
	    //TAxis3D rulers;
	    //rulers.Draw();
	    for(int i=0;i<6;++i){
	      track[i]= new TPolyLine3D();
	      hit[i]=new TPolyMarker3D(1,8);
	      if(i==0 || i==1 || i==2){
		track[i]->SetLineColor(kRed);
		hit[i]->SetMarkerColor(kRed);
	      }
	      else{
		track[i]->SetLineColor(kBlue);
		hit[i]->SetMarkerColor(kBlue);
	      }
	      track[i]->SetNextPoint(0,0,0);
	      double P=(*pAlpha[i]).P();
	      //track[i]->SetNextPoint((*pAlpha[i]).Px()*10000,(*pAlpha[i]).Pz()*10000,(*pAlpha[i]).Py()*10000);
	      track[i]->SetNextPoint(pos[i].X(),pos[i].Z(),pos[i].Y());
	      hit[i]->SetPoint(0,pos[i].X(),pos[i].Z(),pos[i].Y());
	      if(1) track[i]->Draw("same");
	      else track[i]->Draw("A");
	      hit[i]->Draw("same");
	    }
	    /*
	    c_track->Modified();
	    c_track->Update();
	    c_track->Print(fnam_track.c_str());
	    */
	    //sleep(1);
	  }
	}
      }
    }
    
    for(int seg=0;seg<5;++seg){
      for(int r=0;r<16;++r){
	if(Energy[seg][r]!=-1) Esum_f[seg]+=Energy[seg][r];
      }
    }
    for(int ch=0;ch<40;++ch){
      int seg=ch2seg(ch);
      if(Energy_MADC[ch]!=-1) Esum_r[seg]+=Energy_MADC[ch];
    }
    for(int seg=0;seg<5;++seg) h_Esum_FR[seg]->Fill(Esum_f[seg],Esum_r[seg]);
    
    if(flag[0]!=-1 && flag[0]==flag[1] && flag[1]==flag[2]){
      h_K_theta[0]->Fill((C12_1.E()-M_12C)*1000,C12_1.Theta()*R_TO_D);
      h_hit_3a->Fill(flag[0]);
      h1_theta_12C[1][Ex_index[0]-2][Ex_index[1]]->Fill(c1_cm.Theta()*R_TO_D);
      h1_theta_12C_ratio[Ex_index[0]-2][Ex_index[1]]->Fill(c1_cm.Theta()*R_TO_D);
      h1_3aInSeg_resi->Fill(((W-(*pAlpha[0]+*pAlpha[1]+*pAlpha[2])).M()-M_12C)*1000);
    }
    if(flag[3]!=-1 && flag[3]==flag[4] && flag[4]==flag[5]){
      h_K_theta[0]->Fill((C12_2.E()-M_12C)*1000,C12_2.Theta()*R_TO_D);
      h_hit_3a->Fill(flag[3]);
      h1_theta_12C[1][Ex_index[0]-2][Ex_index[1]]->Fill(c2_cm.Theta()*R_TO_D);
      h1_theta_12C_ratio[Ex_index[0]-2][Ex_index[1]]->Fill(c2_cm.Theta()*R_TO_D);
      h1_3aInSeg_resi->Fill(((W-(*pAlpha[3]+*pAlpha[4]+*pAlpha[5])).M()-M_12C)*1000);
    }

    if(flag[0]!=-1 && flag[1]!=-1 && flag[2]!=-1) h_K_theta[1]->Fill((C12_1.E()-M_12C)*1000,C12_1.Theta()*R_TO_D);
    if(flag[3]!=-1 && flag[4]!=-1 && flag[5]!=-1) h_K_theta[1]->Fill((C12_2.E()-M_12C)*1000,C12_2.Theta()*R_TO_D);
      
    for(int seg=0;seg<5;seg++){
      for(int r=0;r<16;r++){
	for(int th=0;th<8;th++){
	  hit_count_add+=hit_count[seg][r][th];
	  Ans_Npix+=hit_count[seg][r][th];
	}
      }
    }
    h_count->Fill(hit_count_add);
    if(hit_count_add==6){
      count_6a[N_alpha]++;
      for(int k=0;k<6;++k) h_ene_6a->Fill((pAlpha[k]->E() - M_alpha)*1000); // kinetic energy1 [MeV]
      TLorentzVector temp_c12[2];
      TVector3 vect(0,0,1);
      for(int i=0;i<2;++i){
	/*
	  vector<double> ene_a;
	  for(int j=0;j<3;++j) ene_a.push_back((pAlpha[3*i+j]->E() - M_alpha)*1000);
	  sort(ene_a.begin(),ene_a.end());
	  for(int j=1;j<3;++j) h_ene_3a[i]->Fill(ene_a[j]-ene_a[0]);
	*/
	  
	for(int j=0;j<3;++j) temp_c12[i]+=*pAlpha[3*i+j];
	//h_ene_3a[i]->Fill(temp_c12.Phi()*R_TO_D);
	h_inv2C->Fill((temp_c12[i].M()-M_12C)*1000);
      }
      h_ene_3a[0]->Fill(fabs(temp_c12[0].Phi()-temp_c12[1].Phi())*R_TO_D);      
    }

    int N_TDC=0;
    for(int i_ch=0;i_ch<40;++i_ch){
      if(TDC->at(i_ch).size()>0){
	/*
	double random_eff=gRandom->Uniform(0,1);
	if(h2_EMADC_eff->GetBinContent(h2_EMADC_eff->FindBin(i_ch,Energy_MADC[i_ch]))<random_eff){
	  if(i_mon==1 || (0<Energy_MADC[i_ch] && Energy_MADC[i_ch]<14.8)){
	    TDC->at(i_ch).clear();
	  }
	}
	*/
	if(TDC->at(i_ch).size()>0) N_TDC++;		
      }
    }
    //if(Ans_Npix>0) tree->Fill();
    if(mode==5){
      if(N_TDC>0) tree->Fill();
    }
    else if(N_TDC>=3) tree->Fill();
  } // end N_LOOP
  cout << "\rGenerate event : " << setw(10) << n << endl;

  tree->AutoSave();
  tfile->Close();

  func_return ret={};
  ret.all_event=n;
  for(int i=0;i<7;++i){
    ret.count[i]=h_count->GetBinContent(i+1);
    ret.count_6a[i]=count_6a[i];
  }
  
  /* draw multiplicity graph */
  TGraph *g_count_sum_6a[7];


  TText *text_sum_6a[7];
  TMarker *marker[7];
  TMarker *marker_sum[7];
  
  TMultiGraph *mg_sum_6a = new TMultiGraph("multi_g_sum_6a",
					   ";offset[mm];efficiency");
  /*
  for(int i=0;i<7;i++){
    text_sum_6a[i] = new TText(offset_z[OFFSET_LOOP-1]+15,0.3-i*0.05,Form("Multi. >%d",i));
    text_sum_6a[i]->SetTextColor(i+1);
    marker[i] = new TMarker(offset_z[OFFSET_LOOP-1]+30,0.4-i*0.05,i+20);
    marker_sum[i] = new TMarker(offset_z[OFFSET_LOOP-1]+30,1-i*0.15,i+20);
    marker_sum[i]->SetMarkerColor(i+1);
    marker[i]->SetMarkerColor(i+1);
    
    g_count_sum_6a[i] = new TGraph(OFFSET_LOOP,offset_z,eff_sum_6a[i]);
    g_count_sum_6a[i]->SetLineColor(i+1);
    g_count_sum_6a[i]->SetMarkerStyle(i+20);
    g_count_sum_6a[i]->SetMarkerColor(i+1);
    if(i==4){
      text_sum[i]->SetTextColor(kOrange+10);
      text_sum_6a[i]->SetTextColor(kOrange+10);
      g_count_sum[i]->SetLineColor(kOrange+10);
      g_count_sum[i]->SetMarkerColor(kOrange+10);
      g_count_sum_6a[i]->SetLineColor(kOrange+10);
      g_count_sum_6a[i]->SetMarkerColor(kOrange+10);
      text[i]->SetTextColor(kOrange+10);
      marker_sum[i]->SetMarkerColor(kOrange+10);
      marker[i]->SetMarkerColor(kOrange+10);
    }
    
    g_count_sum_6a[i]->SetTitle(Form("Multi.> %d;offset[mm];efficiency",i));
    g_count[i]->SetTitle(Form("Multi. %d;offset[mm];efficiency",i));
    mg_sum->Add(g_count_sum[i]);
    mg_sum_6a->Add(g_count_sum_6a[i]);
    
  }
  
  c_graph_sum->cd()->SetMargin(0.12,0.15,0.1,0.1);
  mg_sum->Draw("APC");
  mg_sum->GetXaxis()->SetTitleSize(0.05);
  mg_sum->GetYaxis()->SetTitleSize(0.05);
  mg_sum->GetXaxis()->SetLabelSize(0.05);
  mg_sum->GetYaxis()->SetLabelSize(0.05);
  for(int i=0;i<7;i++){
    text_sum[i]->Draw();
    //marker_sum[i]->Draw();
  }
  
  c_graph_sum_6a->cd()->SetMargin(0.15,0.15,0.1,0.1);
  mg_sum_6a->Draw("APC");
  mg_sum_6a->GetXaxis()->SetTitleSize(0.05);
  mg_sum_6a->GetYaxis()->SetTitleSize(0.05);
  mg_sum_6a->GetXaxis()->SetLabelSize(0.05);
  mg_sum_6a->GetYaxis()->SetLabelSize(0.05);
  for(int i=0;i<7;i++){
    text_sum_6a[i]->Draw();
    //marker_sum[i]->Draw();
  }

  */
  
  /* draw the histograms */
  h_ene_6a->SetTitle("Energy of alpha;Energy[MeV];count");
  h_ene_6a->SetTitleSize(0.05,"xy");
  h_ene_6a->SetLabelSize(0.05,"xy");
  //h_ene[0]->SetStats(0);
  h_th_phi->SetTitle(";Theta[deg];Phi[deg]");
  h_th_phi->SetTitleSize(0.05,"xy");
  h_th_phi->SetLabelSize(0.05,"xy");
  h_th_phi->SetStats(0);
  h_th->SetTitle(";Theta[deg];count");
  h_th->SetTitleSize(0.05,"xy");
  h_th->SetLabelSize(0.05,"xy");
  h_th->SetStats(0);
  c_ene_ang->cd(1)->SetMargin(0.15,0.05,0.12,0.1);
  //h_ene[0]->Draw();
  h_ene_6a->Draw();
  c_ene_ang->cd(2)->SetMargin(0.15,0.05,0.12,0.1);
  //  h_th_phi->Draw("col");
  h_th->Draw();
  c_det_xy->cd(0); h_det_xy->Draw("zcol");
  for(int i=0;i<5;i++){
    det_pline_xy[i]->Draw("same");
  }
  //h_hit->Draw();
  c_det_xz->cd(0); h_det_xz->Draw("zcol");
  for(int i=0;i<5;++i) det_pline_xz[i]->Draw("same");
  //h_hit->Draw();  
  c_dist->cd(1); h_xy->Draw("zcol");
  c_dist->cd(2); h_xz->Draw("zcol");
  for(int i=0; i<6; i++){
    c_ene->cd(i+1);
    h_ene[i]->Draw();
  }

  for(int i=0;i<OFFSET_LOOP;i++){
    c_count->cd(i+1);
    h_count->Draw();
  }
  c_ene_3a->cd(0);
  h_ene_3a[0]->Draw();
  h_ene_3a[1]->Draw("same");
  h_hit_3a->Draw();
  h_hit_3a->GetYaxis()->SetRangeUser(0,h_hit_3a->GetMaximum()*1.1);
  
  /*  
  c2->cd(0);
  for(int i=0;i<5;i++){
    h_ang[i]->SetLineColor(i+1);
    h_ang[i]->Draw("same");
  }
  h_r->Draw();
  h_theta->Draw();
  h_count[0]->Draw();
  */
  
  c_dist->Clear();
  c_dist->Divide(1,2);
  c_dist->cd(1)->SetGridx(1);
  c_dist->cd(1)->SetGridy(1);
  h_K_theta[0]->Draw("colz");
  c_dist->cd(2)->SetGridx(1);
  c_dist->cd(2)->SetGridy(1);
  h_K_theta[1]->Draw("colz");
  
  string fnam=Form("fig/SAKRA_beam%.1f_off%.1f.pdf",T_beam*1000,offset);
  c_pol->Print((fnam+"[").c_str());
  c_pol->Print(fnam.c_str());
  c_det_xy->Print(fnam.c_str());
  c_det_xz->Print(fnam.c_str());
  c_ene_ang->Print(fnam.c_str());
  c_ene->Print(fnam.c_str());
  c_ene_3a->Print(fnam.c_str());
  c_count->Print(fnam.c_str());
  c_dist->Print(fnam.c_str());
  //c_graph_sum->Print(fnam.c_str());
  //c_graph->Print(fnam.c_str());
  c_pol->Print((fnam+"]").c_str());
  //c_track->Print((fnam_track+"]").c_str());
  
  c_ene_3a->Clear();
  c_ene_3a->Divide(1,2);
  for(int i=0;i<2;++i){
    c_ene_3a->cd(i+1);
    h2_E1_vs_E2[i]->Draw("colz");
  }

  for(int j=0;j<2;++j){
    for(int k=0;k<4;++k)
      h1_theta_12C_ratio[j][k]->Divide(h1_theta_12C[0][j][k]);
  }
  h_Write(-1,i_mon,i_ene,mode);
  
  return ret;
}

void all_ene(){
  double ene[N_ENERGY];

  double *eff[7],*eff_sum[7];
  for(int i=0;i<7;++i){
    eff[i]=new double[N_ENERGY];
    eff_sum[i]=new double[N_ENERGY];  
  }  
  func_return ret[50]={};
  for(int i=0;i<N_ENERGY;++i){
    ene[i]=ene_list[0][i];
    //ret[i]=phase_space_2hoyle(0,i);
    cout << endl;
    for(int j=0;j<7;++j){
      eff[j][i]=(double)ret[i].count[j]/ret[i].all_event;
      long long tmp=0;
      for(int k=j;k<7;++k) tmp+=ret[i].count[k];
      eff_sum[j][i]=(double)tmp/ret[i].all_event;
    }
  }

  int colors[]={1,2,3,4,kOrange+10,6,7,8,9,10};
  TGraph *g_count[7];
  TMultiGraph *mg_count=new TMultiGraph("multi_g",Form(";E_beam[MeV];efficiency"));
  TText *text[7];
  TGraph *g_count_sum[7];
  TMultiGraph *mg_sum = new TMultiGraph("multi_g_sum",Form(";E_beam[MeV];efficiency"));
  TText *text_sum[7];
  for(int i=0;i<7;++i){
    g_count[i]=new TGraph(N_ENERGY,ene,eff[i]);
    g_count[i]->SetLineColor(colors[i]);
    g_count[i]->SetMarkerStyle(i+20);
    g_count[i]->SetMarkerColor(colors[i]);
    g_count[i]->SetTitle(Form("Multi. %d;E_beam[MeV];efficiency",i));
    mg_count->Add(g_count[i]);
    text[i] = new TText(ene[N_ENERGY-1],0.3-i*0.05,Form("Multi. %d",i));
    text[i]->SetTextColor(colors[i]);
    
    g_count_sum[i]=new TGraph(N_ENERGY,ene,eff_sum[i]);
    g_count_sum[i]->SetLineColor(colors[i]);
    g_count_sum[i]->SetMarkerStyle(i+20);
    g_count_sum[i]->SetMarkerColor(colors[i]);
    g_count_sum[i]->SetTitle(Form("Multi.> %d;E_beam[MeV];efficiency",i));
    mg_sum->Add(g_count_sum[i]);
    text_sum[i] = new TText(ene[N_ENERGY-1],1-i*0.1,Form("Multi. >%d",i));
    text_sum[i]->SetTextColor(colors[i]);
  }
  
  TCanvas *c_graph = new TCanvas("c_graph","eff vs offset",700,500);
  c_graph->SetRightMargin(0.15);
  mg_count->Draw("APC");
  c_graph->Modified();
  c_graph->Update();
  double x_max=c_graph->GetUxmax();
  for(int i=0;i<7;i++){
    text[i]->SetX(x_max);
    text[i]->Draw();
  }
  
  TCanvas *c_graph_sum = new TCanvas("c_graph_sum","eff(sum) vs offset",700,500);  
  c_graph_sum->SetRightMargin(0.15);
  mg_sum->Draw("APC");
  c_graph_sum->Modified();
  c_graph_sum->Update();
  x_max=c_graph_sum->GetUxmax();
  for(int i=0;i<7;++i){
    text_sum[i]->SetX(x_max);
    text_sum[i]->Draw();
  }
  
  //TCanvas *c_graph_sum_6a = new TCanvas("c_graph_sum_6a","eff(sum_6a) vs offset",700,500);
  
  string fnam=Form("fig/efficiency_vs_Ebeam.pdf");
  c_graph->Print((fnam+"[").c_str());
  c_graph->Print(fnam.c_str());
  c_graph_sum->Print(fnam.c_str());
  c_graph->Print((fnam+"]").c_str());
  
}

void all_offset(Double_t ene=39){
  int N_offset=20;
  double *offset;
  offset=new double[N_offset];
  double *eff[7],*eff_sum[7];
  for(int i=0;i<7;++i){
    eff[i]=new double[N_offset];
    eff_sum[i]=new double[N_offset];  
  }
  func_return ret[50]={};
  for(int i=0;i<N_offset;++i){
    offset[i]=-20+10*i;
    //ret[i]=phase_space_2hoyle(ene,offset[i]);
    for(int j=0;j<7;++j){
      eff[j][i]=(double)ret[i].count[j]/ret[i].all_event;
      long long tmp=0;
      for(int k=j;k<7;++k) tmp+=ret[i].count[k];
      eff_sum[j][i]=(double)tmp/ret[i].all_event;
    }
  }
  
  int colors[]={1,2,3,4,kOrange+10,6,7,8,9,10};
  TGraph *g_count[7];
  TMultiGraph *mg_count=new TMultiGraph("multi_g",Form("E_beam=%.1f MeV;offset[mm];efficiency",ene));
  TText *text[7];
  TGraph *g_count_sum[7];
  TMultiGraph *mg_sum = new TMultiGraph("multi_g_sum",Form("E_beam=%.1f MeV;offset[mm];efficiency",ene));
  TText *text_sum[7];
  for(int i=0;i<7;++i){
    g_count[i]=new TGraph(N_offset,offset,eff[i]);
    g_count[i]->SetLineColor(colors[i]);
    g_count[i]->SetMarkerStyle(i+20);
    g_count[i]->SetMarkerColor(colors[i]);
    g_count[i]->SetTitle(Form("Multi. %d;offset[mm];efficiency",i));
    mg_count->Add(g_count[i]);
    text[i] = new TText(offset[N_offset-1],0.4-i*0.05,Form("Multi. %d",i));
    text[i]->SetTextColor(colors[i]);

    g_count_sum[i]=new TGraph(N_offset,offset,eff_sum[i]);
    g_count_sum[i]->SetLineColor(colors[i]);
    g_count_sum[i]->SetMarkerStyle(i+20);
    g_count_sum[i]->SetMarkerColor(colors[i]);
    g_count_sum[i]->SetTitle(Form("Multi.> %d;offset[mm];efficiency",i));
    mg_sum->Add(g_count_sum[i]);
    text_sum[i] = new TText(offset[N_offset-1],1-i*0.1,Form("Multi. >%d",i));
    text_sum[i]->SetTextColor(colors[i]);
  }
  
  TCanvas *c_graph = new TCanvas("c_graph","eff vs offset",700,500);
  c_graph->SetRightMargin(0.15);
  mg_count->Draw("APC");
  c_graph->Modified();
  c_graph->Update();
  double x_max=c_graph->GetUxmax();
  for(int i=0;i<7;i++){
    text[i]->SetX(x_max);
    text[i]->Draw();
  }
  
  TCanvas *c_graph_sum = new TCanvas("c_graph_sum","eff(sum) vs offset",700,500);  
  c_graph_sum->SetRightMargin(0.15);
  mg_sum->Draw("APC");
  c_graph_sum->Modified();
  c_graph_sum->Update();
  x_max=c_graph_sum->GetUxmax();
  for(int i=0;i<7;++i){
    text_sum[i]->SetX(x_max);
    text_sum[i]->Draw();
  }
  
  //TCanvas *c_graph_sum_6a = new TCanvas("c_graph_sum_6a","eff(sum_6a) vs offset",700,500);

  string fnam=Form("fig/efficiency_E%.1fMeV.pdf",ene);
  c_graph->Print((fnam+"[").c_str());
  c_graph->Print(fnam.c_str());
  c_graph_sum->Print(fnam.c_str());
  c_graph->Print((fnam+"]").c_str());
}
