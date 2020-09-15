#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "HiggsCP/Inputs/macros/Plot_lept_mutau_NNNTuples.C"



void   Plot_lept_mutau_NNScore(TString directory = "/nfs/dust/cms/user/rasp/storage/cardinia/2017/OutputDNN/May16/predictions_2017/",
			       TString outputDir = "./figures_June16/old/",
			       int year=2017,
			       bool FFmethod = true,  
			       bool useEmbedded = true,
			       bool LargeScale = false,  
			       bool logY = false,
			       bool showSignal = true,
			       bool compareCP = true,
			       int scaleSignal = 10.,
			       bool blindData = false,
			       bool FORCE = false
			       )
{
  
	TString cuts = "(iso_1<0.15&&pt_1>25&&pt_2>30&&abs(eta_1)<2.1&&abs(eta_2)<2.3)*(is_SingleLepTrigger)*"; // m_vis>40
	TString weights = "weight*trigweight_1/trigweight*";

  const int nCategories = 3;
  const int nSigCategories = 1;
  
  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_bs_00",
			      "#phi_{CP} vs DNN score",
			      8,0.,2*TMath::Pi(),
			      weights,
			      cuts + "(IP_signif_PV_with_BS_1>1.5&&IP_signif_PV_with_BS_2>1.5)*(dmMVA_2==0)*", //&&IP_signif_RefitV_with_BS_1>1.5&&IP_signif_RefitV_with_BS_2>1.5
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
						channel,
			      FFmethod,
			      useEmbedded,
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      "_mupi",
			      FORCE
			      );
  }

  for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_bs_01",
			      "#phi_{CP} vs DNN score",
			      8,0.,2*TMath::Pi(),
			      "weight*",
			      cuts + "(IP_signif_PV_with_BS_1>1.5)*(dmMVA_2==1)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
						channel,
			      FFmethod,
			      useEmbedded,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
						"_murho",
			      FORCE
			      );
  }

	for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
    Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_bs_01",
			      "#phi_{CP} vs DNN score",
			      8,0.,2*TMath::Pi(),
			      "weight*",
			      cuts + "(IP_signif_PV_with_BS_1>1.5)*(dmMVA_2==2)*",
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
						channel,
			      FFmethod,
			      useEmbedded,  
			      true,  
			      true,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
						"_mu0a1",
			      FORCE
			      );
  }
	
	for(int categoryIndex=0;categoryIndex<nSigCategories;categoryIndex++){
		Plot_lept_mutau_NNNTuples("predicted_prob:acotautau_bs_01",
						"#phi_{CP} vs DNN score",
						8,0.,2*TMath::Pi(),
						"weight*",
						cuts + "(IP_signif_PV_with_BS_1>1.5)*(dmMVA_2==10)*",
						"Events",
						categoryIndex,
						directory,
						outputDir,
						year,
						channel,
						FFmethod,
						useEmbedded,  
						true,  
						true,
						showSignal,
						compareCP,
						scaleSignal,
						blindData,
						"_mua1",
						FORCE
						);
	}
	  
  bool _logY = false;
  bool _largeScale = false;
  for(int categoryIndex=0;categoryIndex<nCategories;categoryIndex++){
    if(categoryIndex<nSigCategories){
      blindData=true;
      _logY=true;
    }
    else {
      _logY=logY;
    }
    int nbins = 7;
    float xmin = 0.3;
    float xmax = 1.;
    if (categoryIndex==1){
      nbins=5;
      xmax=0.8;
    }

    TString cuts = "(pt_1>21&&puppimt_1<50&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5&&abs(eta_1)<2.1&&m_vis>40)*";
    if (categoryIndex==0) cuts+= "(IP_signif_RefitV_with_BS_1>1.5)*";
    Plot_lept_mutau_NNNTuples("predicted_prob",
			      "NN Score",
			      nbins,xmin,xmax,
			      "weight*",
			      cuts,
			      "Events",
			      categoryIndex,
			      directory,
			      outputDir,
			      year,
			      channel,
			      //false,
			      FFmethod,
			      //false,  
			      useEmbedded,  
			      _largeScale,  
			      _logY,
			      showSignal,
			      compareCP,
			      scaleSignal,
			      blindData,
			      "",
			      FORCE
			      );
    
}
 
