#ifndef DataCards_h
#define DataCards_h

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;


template<class T1, class T2>
vector<T1> extract_first(const map<T1, T2>& v) {
    vector<T1> vFirst;
    for (auto element : v) {
        vFirst.push_back(element.first);
    }
    return vFirst;
};

template<class T1, class T2>
vector<T2> extract_second(const vector<pair<T1, T2> >& v) {
    vector<T2> vSecond;
    for (auto element : v) {
        vSecond.push_back(element.second);
    }
    return vSecond;
};

struct params {
  TString weights;
  TString cuts;
  bool hist2D;
  int nbins;
  double xmin;
  double xmax;
  vector<double> xDNN;
  TString varToPlot;
};

class DataCards {

 public:

  DataCards(TString ditauchannel,
	    TString era,
	    bool embedded, 
	    bool FFmethod, 
	    TString variableCP,
	    map<TString,int> binsperchannel,
	    double xmin,
	    double xmax,
	    vector<double> xDNNSig,
	    vector<double> xDNNZtt,
	    vector<double> xDNNFakes,
	    bool splitBkg,
	    bool useTH1forHiggs,
	    bool useTH2forZtt,
	    bool useTH2forFakes,
	    bool mvaDM,
	    bool applyIPcut,
	    bool applyIPcutOnBkg,
	    bool runSystematic,
	    bool checkPhiModulation); 

  void SetInputDirectory(TString input_dir);
  void SetOutputDirectory(TString output_dir);
  void SetOutputFileName(TString output_filename);

  void SetCutIPmuon(TString CutIP_muon) {
    CutIP_muon_ = CutIP_muon;
  };

  void SetCutIPpion(TString CutIP_pion) {
    CutIP_pion_ = CutIP_pion;
  };



  bool Run(int classIndex, TString channel);

  ~DataCards();

  TH1D * Unfold(TH2D * histInput) {
    
    int nBinsX = histInput->GetNbinsX();
    int nBinsY = histInput->GetNbinsY();
    
    int nBins = nBinsX * nBinsY;
    
    TString NameNew = TString(histInput->GetName())+TString("_unfolded");
    
    TH1D * histOutput = new TH1D(NameNew,"",nBins,0,float(nBins));
    int iBin = 1;
    for (int j=1; j<=nBinsY; ++j) {
      for (int i=1; i<=nBinsX; ++i) {
	histOutput->SetBinContent(iBin,histInput->GetBinContent(i,j));
	histOutput->SetBinError(iBin,histInput->GetBinError(i,j));
	iBin++;
      }
    }
    
    return histOutput;
    
  };


 private:
  bool embedded_;
  bool fakeFactor_;
  bool useTH2forZtt_;
  bool useTH1forHiggs_;
  bool useTH2forFakes_;
  bool splitBkg_;
  bool mvaDM_;
  bool applyIPcut_;
  bool applyIPcutOnBkg_;
  bool runSystematics_;
  bool checkPhiModulation_;
  TString era_;
  TString ditauchannel_;
  TString variableCP_;
  int nbins_;
  double xmin_;
  double xmax_;
  vector<double> xDNNSig_;  
  vector<double> xDNNZtt_;  
  vector<double> xDNNFakes_;  
  map<TString,int> binsperchannel_;
  vector<TString> FFSystematics_;

  TString CutIP_muon_;
  TString CutIP_pion_;

  const TString prefix_ = "-NOMINAL_ntuple_";

  TString EmbedCut_ = "gen_match_2==5";
  TString input_dir_;
  TString output_dir_;
  TString output_filename_;

  bool loadFiles();
  void createOutputFile(int classIndex, TString channel);
  
  void closeOutputFile();
  
  TFile * outputFile_;  
  
  vector<TH1D*> CreateCardsSample(TString sampleName, 
				  params param,
				  bool runSystematics);
  vector<TH1D*> CreateCardsSample(TString sampleName, 
				  params param,
				  TString systematicName);

  TH1D* CreateCardsFakesOrQCD(TString sampleName, params param, TString weightFForQCD);
  vector<TH1D*> CreateCardsFakes(TString sampleName, params param, TString weightFForQCD,bool runSystematics);
  void RunOnCategory(TString category);
  vector<TH1D*> CreateCardsEmbedSyst(params param);

  vector<TString> SystematicsNames = {
    "",
    "CMS_scale_t_1prong_13TeVUp",
    "CMS_scale_t_1prong_13TeVDown",
    "CMS_scale_t_1prong1pizero_13TeVUp",
    "CMS_scale_t_1prong1pizero_13TeVDown",
    "CMS_scale_t_3prong_13TeVUp",
    "CMS_scale_t_3prong_13TeVDown",
    "CMS_scale_t_3prong1pizero_13TeVUp",
    "CMS_scale_t_3prong1pizero_13TeVDown",
    "CMS_shape_dyShape_13TeVUp",
    "CMS_shape_dyShape_13TeVDown",
    "CMS_scale_met_unclustered_13TeVUp",
    "CMS_scale_met_unclustered_13TeVDown",
    "CMS_htt_boson_reso_met_13TeVUp",
    "CMS_htt_boson_reso_met_13TeVDown",
    "CMS_htt_boson_scale_met_13TeVUp",
    "CMS_htt_boson_scale_met_13TeVDown",
    "CMS_htt_ZLShape_mt_1prong_13TeVUp",
    "CMS_htt_ZLShape_mt_1prong_13TeVDown",
    "CMS_htt_ZLShape_mt_1prong1pi_13TeVUp",
    "CMS_htt_ZLShape_mt_1prong1pi_13TeVDown",
    "CMS_htt_ZLShape_et_1prong_13TeVUp",
    "CMS_htt_ZLShape_et_1prong_13TeVDown",
    "CMS_htt_ZLShape_et_1prong1pi_13TeVUp",
    "CMS_htt_ZLShape_et_1prong1pi_13TeVDown",
    "CMS_htt_ZLShape_1prong_13TeVUp",
    "CMS_htt_ZLShape_1prong_13TeVDown",
    "CMS_htt_ZLShape_1prong1pi_13TeVUp",
    "CMS_htt_ZLShape_1prong1pi_13TeVDown",
    "CMS_scale_j_FlavorQCD_13TeVUp",
    "CMS_scale_j_FlavorQCD_13TeVDown",
    "CMS_scale_j_RelativeBal_13TeVUp",
    "CMS_scale_j_RelativeBal_13TeVDown",
    "CMS_scale_j_HF_13TeVUp",
    "CMS_scale_j_HF_13TeVDown",
    "CMS_scale_j_BBEC1_13TeVUp",
    "CMS_scale_j_BBEC1_13TeVDown",
    "CMS_scale_j_EC2_13TeVUp",
    "CMS_scale_j_EC2_13TeVDown",
    "CMS_scale_j_Absolute_13TeVUp",
    "CMS_scale_j_Absolute_13TeVDown",
    "CMS_scale_j_Absolute_2016_13TeVUp",
    "CMS_scale_j_Absolute_2016_13TeVDown",
    "CMS_scale_j_HF_2016_13TeVUp",
    "CMS_scale_j_HF_2016_13TeVDown",
    "CMS_scale_j_EC2_2016_13TeVUp",
    "CMS_scale_j_EC2_2016_13TeVDown",
    "CMS_scale_j_RelativeSample_2016_13TeVUp",
    "CMS_scale_j_RelativeSample_2016_13TeVDown",
    "CMS_scale_j_BBEC1_2016_13TeVUp",
    "CMS_scale_j_BBEC1_2016_13TeVDown",
    "CMS_scale_j_Absolute_2017_13TeVUp",
    "CMS_scale_j_Absolute_2017_13TeVDown",
    "CMS_scale_j_HF_2017_13TeVUp",
    "CMS_scale_j_HF_2017_13TeVDown",
    "CMS_scale_j_EC2_2017_13TeVUp",
    "CMS_scale_j_EC2_2017_13TeVDown",
    "CMS_scale_j_RelativeSample_2017_13TeVUp",
    "CMS_scale_j_RelativeSample_2017_13TeVDown",
    "CMS_scale_j_BBEC1_2017_13TeVUp",
    "CMS_scale_j_BBEC1_2017_13TeVDown",
    "CMS_scale_j_Absolute_2018_13TeVUp",
    "CMS_scale_j_Absolute_2018_13TeVDown",
    "CMS_scale_j_HF_2018_13TeVUp",
    "CMS_scale_j_HF_2018_13TeVDown",
    "CMS_scale_j_EC2_2018_13TeVUp",
    "CMS_scale_j_EC2_2018_13TeVDown",
    "CMS_scale_j_RelativeSample_2018_13TeVUp",
    "CMS_scale_j_RelativeSample_2018_13TeVDown",
    "CMS_scale_j_BBEC1_2018_13TeVUp",
    "CMS_scale_j_BBEC1_2018_13TeVDown",
    "CMS_scale_mu_13TeVUp",
    "CMS_scale_mu_13TeVDown",
    "CMS_scale_e_13TeVUp",
    "CMS_scale_e_13TeVDown",
    "CMS_res_j_13TeVUp",
    "CMS_res_j_13TeVDown",
    "CMS_eff_b_13TeVUp",
    "CMS_eff_b_13TeVDown",
  };


  vector<TString> FFSystematics_mt = {
        "",
	"ff_mt_sub_systUp",
	"ff_mt_sub_systDown",
	//stat unc1
	"ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_ltUp",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_ltUp",
		
	"ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_ltUp",
	

	"ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc1_njets0_mvadm1Up",
	"ff_mt_qcd_stat_unc1_njets0_mvadm1Up",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm1Up",
	"ff_mt_qcd_stat_unc1_njets1_mvadm1Up",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm1Up",
	"ff_mt_qcd_stat_unc1_njets2_mvadm1Up",
	
	
	"ff_mt_wjets_stat_unc1_njets0_mvadm2Up",
	"ff_mt_qcd_stat_unc1_njets0_mvadm2Up",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm2Up",
	"ff_mt_qcd_stat_unc1_njets1_mvadm2Up",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm2Up",
	"ff_mt_qcd_stat_unc1_njets2_mvadm2Up",


	"ff_mt_wjets_stat_unc1_njets0_mvadm10Up",
	"ff_mt_qcd_stat_unc1_njets0_mvadm10Up",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm10Up",
	"ff_mt_qcd_stat_unc1_njets1_mvadm10Up",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm10Up",
	"ff_mt_qcd_stat_unc1_njets2_mvadm10Up",


	//"ff_mt_wjets_stat_unc1_njets0_mvadm11Up",
	//"ff_mt_qcd_stat_unc1_njets0_mvadm11Up",
	
	//"ff_mt_wjets_stat_unc1_njets1_mvadm11Up",
	//"ff_mt_qcd_stat_unc1_njets1_mvadm11Up",
	
	//"ff_mt_wjets_stat_unc1_njets2_mvadm11Up",
	//"ff_mt_qcd_stat_unc1_njets2_mvadm11Up",

	"ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_ltDown",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_ltDown",
		
	"ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_ltDown",
	

	"ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc1_njets0_mvadm1Down",
	"ff_mt_qcd_stat_unc1_njets0_mvadm1Down",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm1Down",
	"ff_mt_qcd_stat_unc1_njets1_mvadm1Down",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm1Down",
	"ff_mt_qcd_stat_unc1_njets2_mvadm1Down",
	
	
	"ff_mt_wjets_stat_unc1_njets0_mvadm2Down",
	"ff_mt_qcd_stat_unc1_njets0_mvadm2Down",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm2Down",
	"ff_mt_qcd_stat_unc1_njets1_mvadm2Down",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm2Down",
	"ff_mt_qcd_stat_unc1_njets2_mvadm2Down",


	"ff_mt_wjets_stat_unc1_njets0_mvadm10Down",
	"ff_mt_qcd_stat_unc1_njets0_mvadm10Down",
	
	"ff_mt_wjets_stat_unc1_njets1_mvadm10Down",
	"ff_mt_qcd_stat_unc1_njets1_mvadm10Down",
	
	"ff_mt_wjets_stat_unc1_njets2_mvadm10Down",
	"ff_mt_qcd_stat_unc1_njets2_mvadm10Down",


	//"ff_mt_wjets_stat_unc1_njets0_mvadm11Down",
	//"ff_mt_qcd_stat_unc1_njets0_mvadm11Down",
	
	//"ff_mt_wjets_stat_unc1_njets1_mvadm11Down",
	//"ff_mt_qcd_stat_unc1_njets1_mvadm11Down",
	
	//"ff_mt_wjets_stat_unc1_njets2_mvadm11Down",
	//"ff_mt_qcd_stat_unc1_njets2_mvadm11Down",

	//stat unc2
	"ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_ltUp",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_ltUp",
		
	"ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_ltUp",
	"ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_ltUp",
	

	"ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_gtUp",
	"ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_gtUp",
	
	"ff_mt_wjets_stat_unc2_njets0_mvadm1Up",
	"ff_mt_qcd_stat_unc2_njets0_mvadm1Up",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm1Up",
	"ff_mt_qcd_stat_unc2_njets1_mvadm1Up",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm1Up",
	"ff_mt_qcd_stat_unc2_njets2_mvadm1Up",
	
	
	"ff_mt_wjets_stat_unc2_njets0_mvadm2Up",
	"ff_mt_qcd_stat_unc2_njets0_mvadm2Up",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm2Up",
	"ff_mt_qcd_stat_unc2_njets1_mvadm2Up",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm2Up",
	"ff_mt_qcd_stat_unc2_njets2_mvadm2Up",


	"ff_mt_wjets_stat_unc2_njets0_mvadm10Up",
	"ff_mt_qcd_stat_unc2_njets0_mvadm10Up",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm10Up",
	"ff_mt_qcd_stat_unc2_njets1_mvadm10Up",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm10Up",
	"ff_mt_qcd_stat_unc2_njets2_mvadm10Up",


	//"ff_mt_wjets_stat_unc2_njets0_mvadm11Up",
	//"ff_mt_qcd_stat_unc2_njets0_mvadm11Up",
	
	//"ff_mt_wjets_stat_unc2_njets1_mvadm11Up",
	//"ff_mt_qcd_stat_unc2_njets1_mvadm11Up",
	
	//"ff_mt_wjets_stat_unc2_njets2_mvadm11Up",
	//"ff_mt_qcd_stat_unc2_njets2_mvadm11Up",

	"ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_ltDown",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_ltDown",
		
	"ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_ltDown",
	"ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_ltDown",
	

	"ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_gtDown",
	"ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_gtDown",
	
	"ff_mt_wjets_stat_unc2_njets0_mvadm1Down",
	"ff_mt_qcd_stat_unc2_njets0_mvadm1Down",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm1Down",
	"ff_mt_qcd_stat_unc2_njets1_mvadm1Down",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm1Down",
	"ff_mt_qcd_stat_unc2_njets2_mvadm1Down",
	
	
	"ff_mt_wjets_stat_unc2_njets0_mvadm2Down",
	"ff_mt_qcd_stat_unc2_njets0_mvadm2Down",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm2Down",
	"ff_mt_qcd_stat_unc2_njets1_mvadm2Down",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm2Down",
	"ff_mt_qcd_stat_unc2_njets2_mvadm2Down",


	"ff_mt_wjets_stat_unc2_njets0_mvadm10Down",
	"ff_mt_qcd_stat_unc2_njets0_mvadm10Down",
	
	"ff_mt_wjets_stat_unc2_njets1_mvadm10Down",
	"ff_mt_qcd_stat_unc2_njets1_mvadm10Down",
	
	"ff_mt_wjets_stat_unc2_njets2_mvadm10Down",
	"ff_mt_qcd_stat_unc2_njets2_mvadm10Down",


	//"ff_mt_wjets_stat_unc2_njets0_mvadm11Down",
	//"ff_mt_qcd_stat_unc2_njets0_mvadm11Down",
	
	//"ff_mt_wjets_stat_unc2_njets1_mvadm11Down",
	//"ff_mt_qcd_stat_unc2_njets1_mvadm11Down",
	
	//"ff_mt_wjets_stat_unc2_njets2_mvadm11Down",
	//"ff_mt_qcd_stat_unc2_njets2_mvadm11Down",

	//met_var_qcd and met_var_w non-closure corrections

	"ff_mt_qcd_met_closure_systUp",
	"ff_mt_wjets_met_closure_systUp",
	"ff_mt_qcd_met_closure_systDown",
	"ff_mt_wjets_met_closure_systDown",
	"ff_mt_qcd_met_closure_syst_njets0Up",
	"ff_mt_wjets_met_closure_syst_njets0Up",
	"ff_mt_qcd_met_closure_syst_njets0Down",
	"ff_mt_wjets_met_closure_syst_njets0Down",
	"ff_mt_qcd_met_closure_syst_njets1Up",
	"ff_mt_wjets_met_closure_syst_njets1Up",
	"ff_mt_qcd_met_closure_syst_njets1Down",
	"ff_mt_wjets_met_closure_syst_njets1Down",
	"ff_mt_wjets_met_closure_syst_njets2Up",
	"ff_mt_wjets_met_closure_syst_njets2Down",
	"ff_mt_ttbar_met_closure_systUp",
	"ff_mt_ttbar_met_closure_systDown",

	//m_pt non-closure corrections

	"ff_mt_qcd_l_pt_closure_systUp",
	"ff_mt_qcd_l_pt_closure_systDown",
	"ff_mt_wjets_l_pt_closure_systUp",
	"ff_mt_wjets_l_pt_closure_systDown",
	"ff_mt_qcd_l_pt_closure_syst_njets0Up",
	"ff_mt_qcd_l_pt_closure_syst_njets0Down",
	"ff_mt_wjets_l_pt_closure_syst_njets0Up",
	"ff_mt_wjets_l_pt_closure_syst_njets0Down",
	"ff_mt_qcd_l_pt_closure_syst_njets1Up",
	"ff_mt_qcd_l_pt_closure_syst_njets1Down",
	"ff_mt_wjets_l_pt_closure_syst_njets1Up",
	"ff_mt_wjets_l_pt_closure_syst_njets1Down",
	"ff_mt_wjets_l_pt_closure_syst_njets2Up",
	"ff_mt_wjets_l_pt_closure_syst_njets2Down",

	//extrapolations from DR to SR
	"ff_mt_qcd_systUp",
	"ff_mt_qcd_systDown",
	"ff_mt_qcd_syst_njets0Up",
	"ff_mt_qcd_syst_njets0Down",
	"ff_mt_qcd_syst_njets1Up",
	"ff_mt_qcd_syst_njets1Down",
	"ff_mt_wjets_systUp",
	"ff_mt_wjets_systDown",
	"ff_mt_wjets_syst_njets0Up",
	"ff_mt_wjets_syst_njets0Down",
	"ff_mt_wjets_syst_njets1Up",
	"ff_mt_wjets_syst_njets1Down",
	"ff_mt_wjets_syst_njets2Up",
	"ff_mt_wjets_syst_njets2Down",
	"ff_mt_ttbar_systUp",
	"ff_mt_ttbar_systDown"

  };

  vector<TString> WeightSystematics = {
    "CMS_eff_Xtrigger_mt_MVADM0_13TeVUp",
    "CMS_eff_Xtrigger_mt_MVADM1_13TeVUp",
    "CMS_eff_Xtrigger_mt_MVADM2_13TeVUp",
    "CMS_eff_Xtrigger_mt_MVADM10_13TeVUp",
    //"CMS_eff_Xtrigger_mt_MVADM11_13TeVUp",
    "CMS_eff_Xtrigger_mt_MVADM0_13TeVDown",
    "CMS_eff_Xtrigger_mt_MVADM1_13TeVDown",
    "CMS_eff_Xtrigger_mt_MVADM2_13TeVDown",
    "CMS_eff_Xtrigger_mt_MVADM10_13TeVDown",
    //"CMS_eff_Xtrigger_mt_MVADM11_13TeVDown",
    "CMS_eff_Xtrigger_et_MVADM0_13TeVUp",
    "CMS_eff_Xtrigger_et_MVADM1_13TeVUp",
    "CMS_eff_Xtrigger_et_MVADM2_13TeVUp",
    "CMS_eff_Xtrigger_et_MVADM10_13TeVUp",
    //"CMS_eff_Xtrigger_et_MVADM11_13TeVUp",
    "CMS_eff_Xtrigger_et_MVADM0_13TeVDown",
    "CMS_eff_Xtrigger_et_MVADM1_13TeVDown",
    "CMS_eff_Xtrigger_et_MVADM2_13TeVDown",
    "CMS_eff_Xtrigger_et_MVADM10_13TeVDown",
    //"CMS_eff_Xtrigger_et_MVADM11_13TeVDown",
    "CMS_efake_et_MVADM0_13TeVUp",
    "CMS_efake_et_MVADM1_13TeVUp",
    "CMS_efake_et_MVADM2_13TeVUp",
    "CMS_efake_et_MVADM10_13TeVUp",
    //"CMS_efake_et_MVADM11_13TeVUp",
    "CMS_efake_et_MVADM0_13TeVDown",
    "CMS_efake_et_MVADM1_13TeVDown",
    "CMS_efake_et_MVADM2_13TeVDown",
    "CMS_efake_et_MVADM10_13TeVDown",
    //"CMS_efake_mt_MVADM11_13TeVDown",
    "CMS_mufake_mt_MVADM0_13TeVUp",
    "CMS_mufake_mt_MVADM1_13TeVUp",
    "CMS_mufake_mt_MVADM2_13TeVUp",
    "CMS_mufake_mt_MVADM10_13TeVUp",
    //"CMS_mufake_mt_MVADM11_13TeVUp",
    "CMS_mufake_mt_MVADM0_13TeVDown",
    "CMS_mufake_mt_MVADM1_13TeVDown",
    "CMS_mufake_mt_MVADM2_13TeVDown",
    "CMS_mufake_mt_MVADM10_13TeVDown",
    //"CMS_mufake_mt_MVADM11_13TeVDown",
    "CMS_eff_t_pTlow_MVADM0_13TeVUp", 
    "CMS_eff_t_pTlow_MVADM1_13TeVUp", 
    "CMS_eff_t_pTlow_MVADM2_13TeVUp", 
    "CMS_eff_t_pTlow_MVADM10_13TeVUp",
    //"CMS_eff_t_pTlow_MVADM11_13TeVUp",
    "CMS_eff_t_pThigh_MVADM0_13TeVUp",
    "CMS_eff_t_pThigh_MVADM1_13TeVUp",
    "CMS_eff_t_pThigh_MVADM2_13TeVUp",
    "CMS_eff_t_pThigh_MVADM10_13TeVUp", 
    //"CMS_eff_t_pThigh_MVADM11_13TeVUp", 
    "CMS_eff_t_pTlow_MVADM0_13TeVDown", 
    "CMS_eff_t_pTlow_MVADM1_13TeVDown", 
    "CMS_eff_t_pTlow_MVADM2_13TeVDown", 
    "CMS_eff_t_pTlow_MVADM10_13TeVDown", 
    //"CMS_eff_t_pTlow_MVADM11_13TeVDown", 
    "CMS_eff_t_pThigh_MVADM0_13TeVDown", 
    "CMS_eff_t_pThigh_MVADM1_13TeVDown", 
    "CMS_eff_t_pThigh_MVADM2_13TeVDown", 
    "CMS_eff_t_pThigh_MVADM10_13TeVDown",
    //"CMS_eff_t_pThigh_MVADM11_13TeVDown",
    "CMS_PreFire_13TeVUp",
    "CMS_PreFire_13TeVDown",
  };


  vector<TString> sampleNames = {
    "data_obs",
    "EmbedZTT",
    "ZTT",
    "ZL",
    //"ST",
    "VVT",
    "TTT",
    "QCD",
    "jetFakes",
    "W",
    "HToWW",
    "ZH_sm_htt125",
    "ZH_ps_htt125",
    "ZH_mm_htt125",
    "ZH_flat_htt125",
    "WH_sm_htt125",
    "WH_ps_htt125",
    "WH_mm_htt125",
    "WH_flat_htt125",
    "qqH_sm_htt125",
    "qqH_ps_htt125",
    "qqH_mm_htt125",
    "qqH_flat_htt125",
    "ggH_sm_htt125",
    "ggH_ps_htt125",
    "ggH_mm_htt125",
    "ggH_flat_htt125"
  };

  vector<TString> samplesToSubtract = {
    "EmbedZTT",
    "ZTT",
    "ZL",
    //"ST",
    "VVT",
    "TTT",
    "W",
//    "HToWW"
  };

  vector<TString> fileNames;
  vector<TString> fileNamesOutputDNN = {
    "data",
    "EMB",
    "DY",
    "HToWW",
    "VV",
    "W",
    //"ST",
    "TT",
    "ggH125",
    "qqH125",
    "ZH125",
    "WH125"
  }; 
  vector<TString> fileNamesInputDNN = {
    "SingleMuon",
    "EmbeddedMuTau",
    "DYJets",
    "HToWW",
    "Diboson",
    "WJets",
    "TTbar",
    "GluGluHToUncorrTauTau",
    "VBFHToUncorrTauTau",
    "ZHToUncorrTauTau",
    "WHToUncorrTauTau"
  }; 

  vector<TFile*> filePointer;

  map<TString, TString> mapSampleFileName;

  map<TString, TString> mapSampleFileNameOutputDNN ={
    {"data_obs","data"},
    {"EmbedZTT","EMB"},
    {"ZL","DY"},
    {"ZTT","DY"},
    {"EWKZ","EWKZ"},
    //{"ST","ST"},
    {"TTT","TT"},
    {"VVT","VV"},
    {"W","W"},
    {"QCD","data"},
    {"jetFakes","data"},
    {"HToWW","HToWW"},
    {"ggH_sm_htt125","ggH125"},
    {"ggH_ps_htt125","ggH125"},
    {"ggH_mm_htt125","ggH125"},
    {"ggH_flat_htt125","ggH125"},
    {"qqH_sm_htt125","qqH125"},
    {"qqH_ps_htt125","qqH125"},
    {"qqH_mm_htt125","qqH125"},
    {"qqH_flat_htt125","qqH125"},
    {"ZH_sm_htt125","ZH125"},
    {"ZH_ps_htt125","ZH125"},
    {"ZH_mm_htt125","ZH125"},
    {"ZH_flat_htt125","ZH125"},
    {"WH_sm_htt125","WH125"},
    {"WH_ps_htt125","WH125"},
    {"WH_mm_htt125","WH125"},
    {"WH_flat_htt125","WH125"},
  };
  map<TString, TString> mapSampleFileNameInputDNN ={
    {"data_obs","SingleMuon"},
    {"EmbedZTT","EmbeddedMuTau"},
    {"ZL","DYJets"},
    {"ZTT","DYJets"},
    {"EWKZ","EWKZ"},
    {"TTT","TTbar"},
    {"VVT","Diboson"},
    {"W","WJets"},
    {"QCD","SingleMuon"},
    {"jetFakes","SingleMuon"},
    {"HToWW","HToWW"},
    {"ggH_sm_htt125","GluGluHToUncorrTauTau"},
    {"ggH_ps_htt125","GluGluHToUncorrTauTau"},
    {"ggH_mm_htt125","GluGluHToUncorrTauTau"},
    {"ggH_flat_htt125","GluGluHToUncorrTauTau"},
    {"qqH_sm_htt125","VBFHToUncorrTauTau"},
    {"qqH_ps_htt125","VBFHToUncorrTauTau"},
    {"qqH_mm_htt125","VBFHToUncorrTauTau"},
    {"qqH_flat_htt125","VBFHToUncorrTauTau"},
    {"ZH_sm_htt125","ZHToUncorrTauTau"},
    {"ZH_ps_htt125","ZHToUncorrTauTau"},
    {"ZH_mm_htt125","ZHToUncorrTauTau"},
    {"ZH_flat_htt125","ZHToUncorrTauTau"},
    {"WH_sm_htt125","WHToUncorrTauTau"},
    {"WH_ps_htt125","WHToUncorrTauTau"},
    {"WH_mm_htt125","WHToUncorrTauTau"},
    {"WH_flat_htt125","WHToUncorrTauTau"},
  };

  map<TString, TFile*> mapSampleFile;

  vector<TString> categories;/* = {
    "mt_mupi_sig",
    "mt_mupi_ztt",
    "mt_mupi_fakes",
    //    "mt_mupi_misc",
    "mt_murho_sig",
    "mt_murho_ztt",
    "mt_murho_fakes",
    //    "mt_murho_misc",
    "mt_mua1_sig",
    "mt_mua1_ztt",
    "mt_mua1_fakes",
    //    "mt_ztt",
    //    "mt_fakes"
    };*/


  map<int,TString> classNames;
  map<int,TString> standardClassNames = {
    {0,"sig"},
    {1,"ztt"},
    {2,"fakes"},
  };
  map<int,TString> extraClassNames = {
    {-2,"alphaLtPiOver4"},
    {-3,"alphaGtPiOver4"},
  };
  vector<TString> channelNames = {"mupi","murho","mua1","mu0a1"};

  map<TString,TString> mapCategoryCut;

  void createCategoryList(int classIndex, TString channel);
  void setCategoryCuts();


};



#endif
