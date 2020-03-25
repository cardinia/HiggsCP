#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "HiggsCP/Inputs/macros/Plot_lept_mutau_NNNTuples.C"


void   Plot_lept_mutau_NNScore(TString directory = "/nfs/dust/cms/user/rasp/storage/cardinia/2018/OutputDNN/March18/predictions_2018/",
			       TString outputDir = "./figures_March18_DNN/2018/categories/",
			       int year=2018,
			       bool FFmethod = true,  
			       bool useEmbedded = true,
			       bool LargeScale = false,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = true,
			       int scaleSignal = 100.,
			       bool blindData = true,
			       bool FORCE = false
			       )
{
  
  const int nCategories = 3;
  const int nSigCategories = 1;
  
  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_refitbs_00",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "weight*",
			      "(pt_1>21&&pt_2>20&&m_vis>40&&TMath::Abs(eta_1)<2.1&&puppimt_1<50)*(IP_signif_RefitV_with_BS_1>1.5)*(dmMVA_2==0)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir + "mupi/",
			      year,
			      FFmethod,
			      useEmbedded,
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }

  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_refitbs_01",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "weight*",
			      "(pt_1>21&&pt_2>20&&m_vis>40&&TMath::Abs(eta_1)<2.1&&puppimt_1<50)*(IP_signif_RefitV_with_BS_1>1.5)*(dmMVA_2==1)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir + "murho/",
			      year,
			      FFmethod,
			      useEmbedded,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
	
  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_refitbs_01",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "weight*",
			      "(pt_1>21&&pt_2>20&&m_vis>40&&TMath::Abs(eta_1)<2.1&&puppimt_1<50)*(IP_signif_RefitV_with_BS_1>1.5)*(dmMVA_2==10)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir + "mua1/",
			      year,
			      FFmethod,
			      useEmbedded,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
  
  bool _logY = false;
  bool _largeScale = false;
  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    if(categoryIndex<nSigCategories){
      blindData=true;
      _logY=true;
      _largeScale=true;
    }
    else {
      _logY=logY;
      _largeScale=LargeScale;
    }
    Plot_lept_mutau_NNNTuples("predicted_prob",
			      "NN Score",
			      50,0.,1.,
			      "weight*",
			      "(pt_1>21&&pt_2>20&&m_vis>40&&TMath::Abs(eta_1)<2.1&&puppimt_1<50)*(IP_signif_RefitV_with_BS_1>1.5)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      FFmethod,
			      useEmbedded,  
			      _largeScale,  
			      _logY,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      FORCE
			      );
  }
}
