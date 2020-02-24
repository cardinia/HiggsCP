#ifndef DataCards_h
#define DataCards_h

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <vector>
#include <iostream>

using namespace std;


template<class T1, class T2>
vector<T1> extract_first(const vector<pair<T1, T2> >& v) {
    vector<T1> vFirst;
    for (size_t i = 0; i < v.size(); i++) {
        vFirst.push_back(v[i].first);
    }
    return vFirst;
};

template<class T1, class T2>
vector<T2> extract_second(const vector<pair<T1, T2> >& v) {
    vector<T2> vSecond;
    for (size_t i = 0; i < v.size(); i++) {
        vSecond.push_back(v[i].second);
    }
    return vSecond;
};

struct params {
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

  DataCards(TString era,
	    bool embedded, 
	    TString variableCP,
	    int nbins,
	    double xmin,
	    double xmax,
	    vector<double> xDNN,
	    bool useTH2forZtt,
	    bool mvaDM,
	    bool applyIPcut,
	    bool runSystematic); 

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
  bool mvaDM_;
  bool applyIPcut_;
  bool runSystematics_;
  TString era_;
  TString variableCP_;
  int nbins_;
  double xmin_;
  double xmax_;
  vector<double> xDNN_;  

  TString CutIP_muon_;
  TString CutIP_pion_;

  const TString prefix_ = "mt-NOMINAL_ntuple_";

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

  TH1D* CreateCardsFakesOrQCD(TString sampleName, params param, TString weightFForQCD);
  void RunOnCategory(TString category);
  
  vector<TString> SystematicsNames = {
    "",
    "CMS_scale_t_1prong_13TeVUp",
    "CMS_scale_t_1prong_13TeVDown",
    "CMS_scale_t_1prong1pizero_13TeVUp",
    "CMS_scale_t_1prong1pizero_13TeVDown",
    "CMS_scale_t_3prong_13TeVUp",
    "CMS_scale_t_3prong_13TeVDown",
    "CMS_shape_dyShape_13TeVUp",
    "CMS_shape_dyShape_13TeVDown",
    "topPtWeightUp",
    "topPtWeightDown",
    "CMS_scale_met_unclustered_13TeVUp",
    "CMS_scale_met_unclustered_13TeVDown",
    "CMS_scale_met_boson_resolution_13TeVUp",
    "CMS_scale_met_boson_resolution_13TeVDown",
    "CMS_scale_met_boson_response_13TeVUp",
    "CMS_scale_met_boson_response_13TeVDown",
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
    "CMS_res_j_13TeVUp",
    "CMS_res_j_13TeVDown"
  };

  vector<TString> sampleNames = {
    "data_obs",
    "EMB",
    "ZTT",
    "ZLL",
    "ST",
    "VV",
    "TT",
    "QCD",
    "fakes",
    "W",
    "qqH_sm_htt125",
    "qqH_ps_htt125",
    "qqH_mm_htt125",
    "ggH_sm_htt125",
    "ggH_ps_htt125",
    "ggH_mm_htt125"
  };

  vector<TString> samplesToSubtract = {
    "EMB",
    "ZTT",
    "ZLL",
    "ST",
    "VV",
    "TT",
    "W"
  };

  vector<TString> fileNames = {
    "data",
    "EMB",
    "DY",
    "VV",
    "W",
    "ST",
    "TT",
    "ggH125",
    "qqH125"
  }; 

  vector<TFile*> filePointer;

  map<TString, TString> mapSampleFileName ={
    {"data_obs","data"},
    {"EMB","EMB"},
    {"ZLL","DY"},
    {"ZTT","DY"},
    {"ST","ST"},
    {"TT","TT"},
    {"VV","VV"},
    {"W","W"},
    {"QCD","data"},
    {"fakes","data"},
    {"ggH_sm_htt125","ggH125"},
    {"ggH_ps_htt125","ggH125"},
    {"ggH_mm_htt125","ggH125"},
    {"qqH_sm_htt125","qqH125"},
    {"qqH_ps_htt125","qqH125"},
    {"qqH_mm_htt125","qqH125"}
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

  vector<pair<int,TString>> classNames = {
    {0,"sig"},
    {1,"ztt"},
    {2,"fakes"},
  };
  vector<TString> channelNames = {"mupi","murho","mua1"};

  map<TString,TString> mapCategoryCut;

  void createCategoryList(int classIndex, TString channel);
  void setCategoryCuts();


};



#endif
