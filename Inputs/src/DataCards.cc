#define DataCards_ccx

#include "HiggsCP/Inputs/interface/DataCards.h"

using namespace std;

DataCards::DataCards(TString era,
		     bool embedded,
		     int nbins,
		     double xmin,
		     double xmax,
		     std::vector<double> xDNN,
		     bool useTH2forZtt,
		     bool mvaDM,
		     bool applyIPcut,
		     bool runSystematic) {

  era_ = era;
  embedded_ = embedded;
  nbins_ = nbins;
  xmin_ = xmin;
  xmax_ = xmax;
  xDNN_ = xDNN;
  useTH2forZtt_ = useTH2forZtt;
  mvaDM_ = mvaDM;
  applyIPcut_ = applyIPcut;
  runSystematics_ = runSystematic;

}

DataCards::~DataCards() {
  for (auto file : filePointer) {
    file->Close();
    delete file;
  }
    
}

void DataCards::SetInputDirectory(TString input_dir) {
  input_dir_ = input_dir;
}

void DataCards::SetOutputDirectory(TString output_dir) {
  output_dir_ = output_dir;
}

void DataCards::SetOutputFileName(TString output_filename) {
  output_filename_ = output_filename;
}

bool DataCards::loadFiles() {

  bool allFilesPresent = true;
  for (auto fileName : fileNames) {
    TString fullName = input_dir_+"/"+prefix_+fileName+".root";
    TFile * file = new TFile(fullName);
    filePointer.push_back(file);
    if (file->IsZombie()) {
      cout << " file " << fullName << " does not exist " << endl;
      allFilesPresent = false;
      continue;
    }
    for (auto sampleName : sampleNames) {
      if (mapSampleFileName[sampleName]==fileName)
	mapSampleFile[sampleName] = file;
    }
  }

  return allFilesPresent;

}

void DataCards::createOutputFile() {

  outputFile_ = new TFile(output_dir_+"/"+output_filename_,"recreate");
  for (auto category : categories) {
    outputFile_->mkdir(category);
  } 

}

void DataCards::closeOutputFile() {
  
  outputFile_->Close();
  delete outputFile_;

}


void DataCards::setCategoryCuts() {
    
  for (auto category : categories) {
    TString catString = "";
    TString dmString = ""; 
    if (category.Contains("_sig"))
      catString = "predicted_class==0";
    else if (category.Contains("_ztt"))
      catString = "predicted_class==1";
    else 
      catString = "predicted_class==2";
    //    else 
    //      catString = "DNN==3";

    TString DM("tau_decay_mode_2");
    if (mvaDM_) 
      DM = "dmMVA_2";
    if (category.Contains("_mupi"))
      dmString = DM+"==0";
    else if (category.Contains("_murho"))
      dmString = DM+"==1";
    else 
      dmString = DM+"==10";

    TString cut = catString + "&&" + dmString;
    mapCategoryCut[category] = cut;
  }

}

vector<TH1D*> DataCards::CreateCardsSample(TString sampleName, params param, bool runSystematics) {

  TFile * fileSample = mapSampleFile[sampleName];

  vector<TH1D*> histos;

  cout << "   Running on sample " << sampleName << " : " << fileSample << endl;

  for (auto systematicName : SystematicsNames ) {
    cout << "     Running on systematics " << systematicName << std::endl;

    if (!runSystematics && systematicName!="") continue;
    TString TreeName = "TauCheck";
    TString HistName = sampleName;
    if (systematicName!="") {
      TreeName = "TauCheck_" + systematicName;
      HistName = sampleName + "_" + systematicName;
    }
    TTree * tree = (TTree*)fileSample->Get(TreeName);
    if (tree==NULL) continue;
    TH2D * hist2D;
    TH1D * hist1D;
    double xbins[10];
    double ybins[10];
    int nbinsx = param.nbins;
    double widthx = (param.xmax-param.xmin)/double(nbinsx);
    for (int i=0; i<=nbinsx; ++i)
      xbins[i] = param.xmin + double(i)*widthx;
    int nbinsy = param.xDNN.size() - 1;
    for (int i=0; i<=nbinsy; ++i)
      ybins[i] = param.xDNN.at(i);
    if (param.hist2D) {      
      hist2D = new TH2D("hist","",
			nbinsx,xbins,nbinsy,ybins);
    }
    else {
      hist1D = new TH1D("hist","",
                        nbinsy,ybins); 
    }
    //    cout << "      plot : " << param.varToPlot << endl;
    //    cout << "       cut : " << param.cuts << endl;
    tree->Draw(param.varToPlot+">>hist",param.cuts);
    if (param.hist2D) {
      TH1D * hist = Unfold(hist2D);
      histos.push_back((TH1D*)hist->Clone(HistName));
      delete hist2D;
      delete hist;
    }
    else {
      histos.push_back((TH1D*)hist1D->Clone(HistName));
      delete hist1D;
    }
  }

  return histos;

}

TH1D * DataCards::CreateCardsFakesOrQCD(TString sample, params parameters, TString weightFForQCD) {
  
  TString weight = "weight*"+weightFForQCD;
  TString cuts = parameters.cuts;
  parameters.cuts = weight + "(" + cuts + ")";
  bool runSystematics = false;

  vector<TH1D*> histFF = CreateCardsSample(sample,parameters,runSystematics);
  for (auto sample : samplesToSubtract) {
    if (sample=="ZTT"&&embedded_) continue;
    if (sample=="EMB"&&!embedded_) continue;

    TString cutsSample = cuts;
    params parametersSample = parameters;

    if (sample=="ZTT")
      cutsSample += "&&gen_match_1==4&&gen_match_2==5";
    else if (sample=="ZLL")
      cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";
    else if (embedded_)
      cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";

    cutsSample += "&&gen_match_2!=6";

    parametersSample.cuts = weight + "(" + cutsSample + ")";

    vector<TH1D*> histos = CreateCardsSample(sample,parametersSample,runSystematics);
    histFF[0]->Add(histFF[0],histos[0],1,-1);

  }

  return histFF[0];

}

void DataCards::RunOnCategory(TString category) {

  cout << "Running on category " << endl;

  vector<TH1D*> allHists; allHists.clear();

  for (auto sampleName : sampleNames) {
    params parameters;
    TString cuts = "";
    TString weight = "weight*";
    TString acotautau = "acotautau_refitbs_00";
    bool runSystematics = true;

    cuts = mapCategoryCut[category] + "&&pt_1>21&&os>0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2>0.5";    


    TString cutsFF = mapCategoryCut[category]  + "&&pt_1>21&&os>0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5";
    TString cutsQCD = mapCategoryCut[category] + "&&pt_1>21&&os<0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2>0.5";

    TString IPCut("");
    if (applyIPcut_) {
      if (category.Contains("_mupi_"))
	IPCut = "&&IP_signif_PV_with_BS_1>1.0&&IP_signif_PV_with_BS_2>1.0";
      else
	IPCut = "&&IP_signif_PV_with_BS_1>1.0";
    }
    cuts += IPCut;
    cutsFF += IPCut;
    cutsQCD += IPCut;
    
    if (sampleName=="ZTT") 
      cuts += "&&gen_match_1==4&&gen_match_2==5";
    else if (sampleName=="ZLL")
      cuts += "&&!(gen_match_1==4&&gen_match_2==5)";
    else if ((sampleName=="TT"||sampleName=="ST"||sampleName=="W"||sampleName=="VV") && embedded_) 
      cuts += "&&!(gen_match_1==4&&gen_match_2==5)";
    

    if (sampleName.Contains("_sm_htt125"))
      weight += "gen_sm_htt125*";
    else if (sampleName.Contains("_ps_htt125"))
      weight += "gen_ps_htt125*";
    else if (sampleName.Contains("_mm_htt125"))
      weight += "gen_mm_htt125*";

    parameters.cuts = weight + "(" + cuts + ")"; 
   
    
    if (category.Contains("_murho_")||category.Contains("_mua1_"))
      acotautau = "acotautau_refitbs_01";    

    if (category.Contains("_sig")) {
      parameters.hist2D = true;
      parameters.varToPlot = "predicted_prob:"+acotautau;
    }
    else if (category.Contains("_ztt")&&useTH2forZtt_) {
      parameters.hist2D = true;
      parameters.varToPlot = "predicted_prob:"+acotautau;      
    }
    else {
      parameters.hist2D = false;
      parameters.varToPlot = "predicted_prob";
    }

    parameters.nbins = nbins_;
    parameters.xmin = xmin_;
    parameters.xmax = xmax_;
    parameters.xDNN = xDNN_;


    if (sampleName=="fakes") {
      params parametersFF = parameters;
      parametersFF.cuts = cutsFF;
      TString weightFF("ff_nom*");
      TH1D * hist = CreateCardsFakesOrQCD(sampleName,parametersFF,weightFF);
      allHists.push_back(hist);
    }
    else if (sampleName=="QCD") {
      params parametersQCD = parameters;
      parametersQCD.cuts = cutsQCD;
      TString weightQCD("1.2*");
      TH1D * hist = CreateCardsFakesOrQCD(sampleName,parametersQCD,weightQCD);
      allHists.push_back(hist); 
    }
    else {
      if (sampleName=="data_obs") runSystematics = false;
      else runSystematics = runSystematics_;
      vector<TH1D*> hists = CreateCardsSample(sampleName, 
					      parameters, 
					      runSystematics);
      for (auto hist : hists)
	allHists.push_back(hist);
    }
  }  

  outputFile_->cd(category);
  for (auto hist : allHists) 
    hist->Write(hist->GetName());
  
} 

bool DataCards::Run() {

  bool isOK = loadFiles();
  if (!isOK) return isOK;

  createOutputFile();

  setCategoryCuts();

  for (auto category : categories) 
    RunOnCategory(category);

  closeOutputFile();

  isOK = true;
  return isOK;

}
