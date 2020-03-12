#define DataCards_ccx

#include "HiggsCP/Inputs/interface/DataCards.h"

using namespace std;

DataCards::DataCards(TString era,
		     bool embedded,
		     TString variableCP,
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
  variableCP_ = variableCP;

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

void DataCards::createOutputFile(int classIndex=-1, TString channel="") {
  if (classIndex==-1&&channel=="") outputFile_ = new TFile(output_dir_+"/"+output_filename_+".root","recreate");
  else outputFile_ = new TFile(output_dir_+"/"+output_filename_+"_"+channel+"_"+classNames[classIndex].second+".root","recreate");
  for (auto category : categories) {
    outputFile_->mkdir(category);
  } 
  cout << "Created file " << outputFile_->GetName() <<endl;
}

void DataCards::closeOutputFile() {
  
  outputFile_->Close();
  delete outputFile_;

}

void DataCards::createCategoryList(int classIndex=-1, TString channel="") {

  TString catNamePrefix = "mt";
  vector<int> classIndices;
  vector<TString> channels;
  if (channel=="") channels = channelNames;
  else channels.push_back(channel);
  if (classIndex==-1) classIndices = extract_first(classNames);
  else classIndices.push_back(classIndex);
  for (auto ch : channels){
    for (auto cl : classIndices){
      categories.push_back(catNamePrefix+"_"+ch+"_"+classNames[cl].second+"_"+era_);
    }
  }
  cout << "***************************************" << endl << "Running on the following categories: " <<endl;
  for (auto category : categories) cout << category << endl;
  cout << "***************************************" << endl << endl;

}



void DataCards::setCategoryCuts() {
    
  for (auto category : categories) {
    TString catString = "predicted_class==";
    TString dmString = ""; 
    for (auto className : classNames){
      if(category.Contains(className.second)) catString += TString::Itoa(className.first,10);
    }/*
    if (category.Contains("_sig"))
      catString = "predicted_class==0";
    else if (category.Contains("_ztt"))
      catString = "predicted_class==1";
    else 
      catString = "predicted_class==2";
    //    else 
    //      catString = "DNN==3";
    */
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


vector<TH1D*> DataCards::CreateCardsEmbedSyst(params param) {

  vector<TH1D*> histos;
  //create nominal embedded template
  vector<TH1D*> histsEmbed = CreateCardsSample("EmbedZTT", param, false);

  TH1D * histEmbedUp = (TH1D *) histsEmbed[0]->Clone("EmbedZTT_CMS_ttbar_embeded_13TeVUp");
  TH1D * histEmbedDown = (TH1D *) histsEmbed[0]->Clone("EmbedZTT_CMS_ttbar_embeded_13TeVDown");

  params paramShift = param;
  //get 10% of TT and VV true taus contributions
  paramShift.weights += "0.1*";
  paramShift.cuts += "&&(gen_match_1==4&&gen_match_2==5)";

  vector<TH1D*> histsTTT = CreateCardsSample("TTT", paramShift, false);
  vector<TH1D*> histsVVT = CreateCardsSample("VVT", paramShift, false);

  histsTTT[0]->Add(histsVVT[0]);

  histEmbedUp->Add(histEmbedUp,histsTTT[0],1,1);
  histEmbedDown->Add(histEmbedDown,histsTTT[0],1,-1);

  histos.push_back(histEmbedUp);
  histos.push_back(histEmbedDown);

  return histos;

}

vector<TH1D*> DataCards::CreateCardsSample(TString sampleName, params param, bool runSystematics) {

  TFile * fileSample = mapSampleFile[sampleName];

  vector<TH1D*> histos;

  cout << "   Running on sample " << sampleName << " : " << fileSample->GetName() << endl;
  vector<TString> SystematicsPerSample = SystematicsNames ;
  SystematicsPerSample.insert(SystematicsPerSample.end(),WeightSystematics.begin(),WeightSystematics.end());
  if(sampleName.Contains("htt")){
    cout << "Adding theory shape systematics" <<endl;
    SystematicsPerSample.push_back("CMS_scale_gg_13TeVUp");
    SystematicsPerSample.push_back("CMS_scale_gg_13TeVDown");
  }
  else if(sampleName=="ZTT"||sampleName=="ZL"){
    cout << "Adding DY shape systematics" <<endl;
    SystematicsPerSample.push_back("CMS_htt_dyShape_13TeVUp");
    SystematicsPerSample.push_back("CMS_htt_dyShape_13TeVDown");
  }

  for (auto systematicName : SystematicsPerSample ) {
    if (!runSystematics && systematicName!="") continue;
    if (systematicName!="") cout << "     Running on systematics " << systematicName << std::endl;
    
    TString TreeName = "TauCheck";
    TString HistName = sampleName;
    params parameterSystematic = param;
    if (systematicName!="") {
      if(find(WeightSystematics.begin(),WeightSystematics.end(),systematicName)==WeightSystematics.end()&&(!(systematicName.Contains("CMS_scale_gg_13TeV")))&&(!(systematicName.Contains("CMS_htt_dyShape_13TeV"))))TreeName = "TauCheck_" + systematicName; 
      else parameterSystematic.weights = param.weights + "weight_" +systematicName+"*";
      HistName = sampleName + "_" + systematicName;
    }

    TTree * tree = (TTree*)fileSample->Get(TreeName);
    //cout << TreeName <<endl;
    if (tree==NULL) continue;
    TH2D * hist2D;
    TH1D * hist1D;
    double xbins[10];
    double ybins[10];
    int nbinsx = parameterSystematic.nbins;
    double widthx = (parameterSystematic.xmax-parameterSystematic.xmin)/double(nbinsx);
    for (int i=0; i<=nbinsx; ++i)
      xbins[i] = parameterSystematic.xmin + double(i)*widthx;
    int nbinsy = parameterSystematic.xDNN.size() - 1;
    for (int i=0; i<=nbinsy; ++i)
      ybins[i] = parameterSystematic.xDNN.at(i);
    if (parameterSystematic.hist2D) {      
      hist2D = new TH2D("hist","",
			nbinsx,xbins,nbinsy,ybins);
    }
    else {
      hist1D = new TH1D("hist","",
                        nbinsy,ybins); 
    }
    //    cout << "      plot : " << parameterSystematic.varToPlot << endl;
    //    cout << "       cut : " << parameterSystematic.cuts << endl;
    tree->Draw(parameterSystematic.varToPlot+">>hist",parameterSystematic.weights+"("+parameterSystematic.cuts+")");
    //cout << parameterSystematic.weights << "(" << parameterSystematic.cuts << ")" << endl;
    if (parameterSystematic.hist2D) {
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

vector<TH1D*> DataCards::CreateCardsSample(TString sampleName, params param, TString systematicName) {

  TFile * fileSample = mapSampleFile[sampleName];

  vector<TH1D*> histos;

  cout << "   Running on sample " << sampleName << " : " << fileSample->GetName() << endl;
  
  if (systematicName!="") cout << "     Running on systematics " << systematicName << std::endl;
    
    TString TreeName = "TauCheck";
    TString HistName = sampleName;
    params parameterSystematic = param;
    //if systematic corresponds to a shift in event weight the histogram should be renamed, but not the tree
    if (systematicName!="") {
      HistName = sampleName + "_" + systematicName;
      if(!systematicName.Contains("jetFakes_ff_mt_sub_syst"))parameterSystematic.weights = param.weights + "weight_" +systematicName+"*";
      else if(sampleName=="fakes")parameterSystematic.weights = param.weights; // small exception: for jetFakes_ff_mt_sub_syst we want the fake bkg to not be scaled, but just the contribution coming from genuine taus
      else if(systematicName=="jetFakes_ff_mt_sub_systUp")parameterSystematic.weights = param.weights + "1.1*";
      else if(systematicName=="jetFakes_ff_mt_sub_systDown")parameterSystematic.weights = param.weights + "0.9*";
    }
    TTree * tree = (TTree*)fileSample->Get(TreeName);
    TH2D * hist2D;
    TH1D * hist1D;
    double xbins[10];
    double ybins[10];
    int nbinsx = parameterSystematic.nbins;
    double widthx = (parameterSystematic.xmax-parameterSystematic.xmin)/double(nbinsx);
    for (int i=0; i<=nbinsx; ++i)
      xbins[i] = parameterSystematic.xmin + double(i)*widthx;
    int nbinsy = parameterSystematic.xDNN.size() - 1;
    for (int i=0; i<=nbinsy; ++i)
      ybins[i] = parameterSystematic.xDNN.at(i);
    if (parameterSystematic.hist2D) {      
      hist2D = new TH2D("hist","",
			nbinsx,xbins,nbinsy,ybins);
    }
    else {
      hist1D = new TH1D("hist","",
                        nbinsy,ybins); 
    }
    //    cout << "      plot : " << parameterSystematic.varToPlot << endl;
    //    cout << "       cut : " << parameterSystematic.cuts << endl;
    tree->Draw(parameterSystematic.varToPlot+">>hist","("+parameterSystematic.weights+parameterSystematic.cuts+")");
    if (parameterSystematic.hist2D) {
      TH1D * hist = Unfold(hist2D);
      histos.push_back((TH1D*)hist->Clone(HistName));
      delete hist2D;
      delete hist;
    }
    else {
      histos.push_back((TH1D*)hist1D->Clone(HistName));
      delete hist1D;
    }
  
    
  return histos;

}

TH1D* DataCards::CreateCardsFakesOrQCD(TString FakeOrQCD, params parameters, TString weightFForQCD) {
  /*  
  TString weight = "weight*"+weightFForQCD;
  TString cuts = parameters.cuts;
  parameters.cuts = weight + "(" + cuts + ")";
  */
  TString cuts = parameters.cuts;
  parameters.weights = "weight*"+weightFForQCD;
  bool runSystematics = false;
  vector<TH1D*> histFF = CreateCardsSample(FakeOrQCD,parameters,runSystematics);
  for (auto sample : samplesToSubtract) {
    if (sample=="ZTT"&&embedded_) continue;
    if (sample=="EmbedZTT"&&!embedded_) continue;

    TString cutsSample = cuts;
    params parametersSample = parameters;

    if (sample=="ZTT")
      cutsSample += "&&gen_match_1==4&&gen_match_2==5";
    else if (sample=="ZL")
      cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";
    else if (embedded_)
      cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";

    if (sample=="EmbedZTT" && era_=="2016") cuts += "&&(weight<1000)";

    cutsSample += "&&gen_match_2!=6";

    //parametersSample.cuts = weight + "(" + cutsSample + ")";
    parametersSample.cuts = cutsSample;

    cout << "  Subtracting " << sample << " for data-driven estimation" << endl;
    vector<TH1D*> histos = CreateCardsSample(sample,parametersSample,runSystematics);
    histFF[0]->Add(histFF[0],histos[0],1,-1);

  }

  return histFF[0];

}

vector<TH1D*> DataCards::CreateCardsFakes(TString FakeOrQCD, params parameters, TString weightFForQCD, bool runSystematics) {
  /*  
  TString weight = "weight*"+weightFForQCD;
  TString cuts = parameters.cuts;
  parameters.cuts = weight + "(" + cuts + ")";
  */
  vector<TH1D*> all_histFF;
  TString cuts = parameters.cuts;
  for (auto systematicName : FFSystematics ) {
    parameters.weights = "weight*"+weightFForQCD;
    if (!runSystematics && systematicName!="") continue;
    if (systematicName!="") cout << "     Running on systematics " << systematicName << std::endl;

    vector<TH1D*> histFF = CreateCardsSample(FakeOrQCD,parameters,systematicName);
    for (auto sample : samplesToSubtract) {
      if (sample=="ZTT"&&embedded_) continue;
      if (sample=="EmbedZTT"&&!embedded_) continue;
      
      TString cutsSample = cuts;
      params parametersSample = parameters;
      
      if (sample=="ZTT")
	cutsSample += "&&gen_match_1==4&&gen_match_2==5";
      else if (sample=="ZL")
	cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";
      else if (embedded_)
	cutsSample += "&&!(gen_match_1==4&&gen_match_2==5)";
      
      if (sample=="EmbedZTT" && era_=="2016") cuts += "&&(weight<1000)";

      cutsSample += "&&gen_match_2!=6";
      
      //parametersSampleUp.cuts = weight + "(" + cutsSample + ")";
      cout << "  Subtracting " << sample << " for data-driven estimation" << endl;
      vector<TH1D*> histos = CreateCardsSample(sample,parametersSample,systematicName);
      histFF[0]->Add(histFF[0],histos[0],1,-1);
      
    }
    all_histFF.push_back(histFF[0]);

  }
  return all_histFF;
  
}
  
void DataCards::RunOnCategory(TString category) {

  cout << "Running on category " << category << endl;

  vector<TH1D*> allHists; allHists.clear();

  for (auto sampleName : sampleNames) {
    params parameters;
    TString cuts = "";
    TString weight = "weight*"; 
    TString acotautau = variableCP_+"_00";
    bool runSystematics = true;

    cuts = mapCategoryCut[category] + "&&pt_1>21&&pt_2>20&&TMath::Abs(eta_1)<2.1&&os>0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2>0.5";    


    TString cutsFF = mapCategoryCut[category]  + "&&pt_1>21&&pt_2>20&&TMath::Abs(eta_1)<2.1&&os>0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5";
    TString cutsQCD = mapCategoryCut[category] + "&&pt_1>21&&pt_2>20&&TMath::Abs(eta_1)<2.1&&os<0.5&&puppimt_1<50&&byMediumDeepTau2017v2p1VSjet_2>0.5";

    TString IPCut("");
    if (applyIPcut_) {
      if (category.Contains("_mupi_"))
	IPCut = "&&"+CutIP_muon_+"&&"+CutIP_pion_;
      else
	IPCut = "&&"+CutIP_muon_;
    }
    cuts += IPCut;
    cutsFF += IPCut;
    cutsQCD += IPCut;
    
    if (sampleName=="ZTT") 
      cuts += "&&gen_match_1==4&&gen_match_2==5";
    else if (sampleName=="ZL")
      cuts += "&&!(gen_match_1==4&&gen_match_2==5)";
    else if ((sampleName=="TTT"||sampleName=="ST"||sampleName=="W"||sampleName=="VVT") && embedded_) 
      cuts += "&&!(gen_match_1==4&&gen_match_2==5)";
    
    if (sampleName=="EmbedZTT" && era_=="2016") cuts += "&&(weight<1000)";


    if (sampleName.Contains("_sm_htt125"))
      weight += "gen_sm_htt125*";
    else if (sampleName.Contains("_ps_htt125"))
      weight += "gen_ps_htt125*";
    else if (sampleName.Contains("_mm_htt125"))
      weight += "gen_mm_htt125*";


    parameters.weights = weight ; 
    parameters.cuts = cuts ; 
   
    
    if (category.Contains("_murho_")||category.Contains("_mua1_"))
      acotautau = variableCP_+"_01";    

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

    if (sampleName=="EmbedZTT" && runSystematics  && embedded_){ 
      vector<TH1D*> hists = CreateCardsEmbedSyst(parameters);
      for (auto hist : hists)
	allHists.push_back(hist);
    }
    if (sampleName=="fakes") {
      params parametersFF = parameters;
      parametersFF.cuts = cutsFF;
      TString weightFF("ff_nom*");
      vector<TH1D*> hists = CreateCardsFakes(sampleName,parametersFF,weightFF, 
					      runSystematics);
      for (auto hist : hists)
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

bool DataCards::Run(int classIndex=-1, TString channel="") {

  bool isOK = loadFiles();
  if (!isOK) return isOK;

  createCategoryList(classIndex,channel);
  setCategoryCuts();
  createOutputFile(classIndex,channel);

  for (auto category : categories) 
    RunOnCategory(category);

  closeOutputFile();

  isOK = true;
  return isOK;

}
