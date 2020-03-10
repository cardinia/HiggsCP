#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "HiggsCP/Inputs/macros/Plot_lept_mutau_NNNTuples.C"


void   Plot_lept_mutau_NNScore(TString directory = "/nfs/dust/cms/user/rasp/storage/cardinia/2018/OutputDNN/March7/predictions_2018/",
			       TString outputDir = "./figures_March10/",
			       int year=2018,
			       bool FFmethod = true,  
			       bool useEmbedded = true,
			       bool LargeScale = false,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = true,
			       int scaleSignal = 1.,
			       bool blindData = false,
			       bool FORCE = false
			       )
{
  
  const int nCategories = 3;
  const int nSigCategories = 1;
  
  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_00",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "weight*",
			      "(puppimt_1<50&&pt_1>20)*(tau_decay_mode_2==0)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
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

  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_01",
			      "#phi_{CP} vs NN score",
			      5,0.,2*TMath::Pi(),
			      "weight*",
			      "(puppimt_1<50&&pt_1>20)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
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
			      10,0.,1.,
			      "weight*",
			      "(puppimt_1<50&&pt_1>20)*",
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
