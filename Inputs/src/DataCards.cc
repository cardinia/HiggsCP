#define DataCards_ccx

#include "HiggsCP/Inputs/interface/DataCards.h"

using namespace std;

DataCards::DataCards(TString ditauchannel,
		     TString era,
		     bool embedded,
		     bool FFmethod,
		     TString variableCP,
		     map<TString,int> binsperchannel,
		     double xmin,
		     double xmax,
		     std::vector<double> xDNNSig,
		     std::vector<double> xDNNZtt,
		     std::vector<double> xDNNFakes,
		     bool splitBkgCat,
		     bool useTH1forHiggs,
		     bool useTH2forZtt,
		     bool useTH2forFakes,
		     bool mvaDM,
		     bool applyIPcut,
		     bool applyIPcutOnBkg,
		     bool runSystematic,
		     bool checkPhiModulation) {

  ditauchannel_ = ditauchannel;
  era_ = era;
  embedded_ = embedded;
  fakeFactor_ = FFmethod;
  binsperchannel_ = binsperchannel;
  xmin_ = xmin;
  xmax_ = xmax;
  xDNNSig_ = xDNNSig;
  xDNNZtt_ = xDNNZtt;
  xDNNFakes_ = xDNNFakes;
  splitBkg_ = splitBkgCat;
  useTH1forHiggs_ = useTH1forHiggs;
  useTH2forZtt_ = useTH2forZtt;
  useTH2forFakes_ = useTH2forFakes;
  mvaDM_ = mvaDM;
  applyIPcut_ = applyIPcut;
  applyIPcutOnBkg_ = applyIPcutOnBkg;
  runSystematics_ = runSystematic;
  variableCP_ = variableCP;
  checkPhiModulation_ = checkPhiModulation;
  if(ditauchannel_=="mt"){
    EmbedCut_=EmbedCut_+"&&gen_match_1==4";
    FFSystematics_ = FFSystematics_mt;
  }else if(ditauchannel_=="et"){
    EmbedCut_=EmbedCut_+"&&gen_match_1==3";
    for (auto systName : FFSystematics_mt){
      systName.ReplaceAll("ff_mt_","ff_et_");
      FFSystematics_.push_back(systName);
    }
  }
  if(checkPhiModulation_){
    classNames = extraClassNames;
    fileNames = fileNamesInputDNN;
    mapSampleFileName = mapSampleFileNameInputDNN;
  }else{
    classNames = standardClassNames;
    fileNames = fileNamesOutputDNN;
    mapSampleFileName = mapSampleFileNameOutputDNN;
  }
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
    TString postfix="";
    if(checkPhiModulation_) postfix=TString("_")+era_;
    TString fullName = input_dir_+"/"+ditauchannel_+prefix_+fileName+postfix+".root";
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
  else if(channel=="") outputFile_ = new TFile(output_dir_+"/"+output_filename_+"_"+classNames[classIndex]+".root","recreate");
  else outputFile_ = new TFile(output_dir_+"/"+output_filename_+"_"+channel+"_"+classNames[classIndex]+".root","recreate");
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
  TString catNamePrefix = ditauchannel_;
  vector<int> classIndices;
  vector<TString> channels;
  if (channel==""&&((splitBkg_&&classIndex!=0)||classIndex==-1||(classIndex==0&&!useTH1forHiggs_))) channels = channelNames;
  else channels.push_back(channel);
  if (classIndex==-1) classIndices = extract_first(classNames);
  else classIndices.push_back(classIndex);
  for (auto cl : classIndices){
    for (auto ch : channels){
      if(ch!="") categories.push_back(catNamePrefix+"_"+ch+"_"+classNames[cl]+"_"+era_);
      else categories.push_back(catNamePrefix+"_"+classNames[cl]+"_"+era_);
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
    }
    TString minusOrplane = "";
    if(category.Contains("_mupi"))
      minusOrplane="minus";
    else minusOrplane="_plane_2";
    if(checkPhiModulation_){
      catString = "alpha";
      catString += minusOrplane;
      if(variableCP_.Contains("uncor")&&category.Contains("_mupi")) catString += "_uncorr";
      if(category.Contains("LtPiOver4")) catString += "<(TMath::Pi()/4.)";
      else if(category.Contains("GtPiOver4")) catString += ">=(TMath::Pi()/4.)";
      else exit(EXIT_FAILURE);
    }
    /*
    if (category.Contains("_sig"))
      catString = "predicted_class==0";
    else if (category.Contains("_ztt"))
      catString = "predicted_class==1";
    else 
      catString = "predicted_class==2";
    //    else 
    //      catString = "DNN==3";
    */

    TString cut = catString;

    TString DM("tau_decay_mode_2");
    if (mvaDM_) 
      DM = "dmMVA_2";
    if(category.Contains("sig")||(splitBkg_&&(category.Contains("ztt")||category.Contains("fakes")))||category.Contains("alpha")){
      if (category.Contains("_mupi"))
	dmString = DM+"==0";
      else if (category.Contains("_murho"))
      dmString = DM+"==1";
      else if (category.Contains("_mu0a1")&&mvaDM_)
	dmString = DM+"==2";
      else 
	dmString = DM+"==10";
    
      cut = cut + "&&" + dmString;
    }
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
    SystematicsPerSample.push_back("CMS_PS_ISR_ggH_13TeVUp");
    SystematicsPerSample.push_back("CMS_PS_ISR_ggH_13TeVDown");
    SystematicsPerSample.push_back("CMS_PS_FSR_ggH_13TeVUp");
    SystematicsPerSample.push_back("CMS_PS_FSR_ggH_13TeVDown");
  }
  else if(sampleName=="ZTT"||sampleName=="ZL"||sampleName=="EWKZ"){
    cout << "Adding DY shape systematics" <<endl;
    SystematicsPerSample.push_back("CMS_htt_dyShape_13TeVUp");
    SystematicsPerSample.push_back("CMS_htt_dyShape_13TeVDown");
    SystematicsPerSample.push_back("CMS_IPsignifCalib_13TeVDown");
  }
  else if(sampleName=="TTT"){
    cout << "Adding ttbar shape systematics" <<endl;
    SystematicsPerSample.push_back("CMS_htt_ttbarShape_13TeVUp");
    SystematicsPerSample.push_back("CMS_htt_ttbarShape_13TeVDown");
  }

  for (auto systematicName : SystematicsPerSample ) {
    if (!runSystematics && systematicName!="") continue;
    if (systematicName!="") cout << "     Running on systematics " << systematicName << std::endl;
    
    TString TreeName = "TauCheck";
    TString HistName = sampleName;
    params parameterSystematic = param;
    if (systematicName!="") {
      if(find(WeightSystematics.begin(),WeightSystematics.end(),systematicName)==WeightSystematics.end()&&(!(systematicName.Contains("CMS_PS_")))&&(!(systematicName.Contains("CMS_scale_gg_13TeV")))&&(!(systematicName.Contains("CMS_htt_dyShape_13TeV")))&&(!(systematicName.Contains("CMS_htt_ttbarShape_13TeV")))&&(!(systematicName.Contains("CMS_IPsignifCalib_13TeV"))))TreeName = "TauCheck_" + systematicName; 
      else if(systematicName.Contains("CMS_IPsignifCalib_13TeVDown")&&!(parameterSystematic.varToPlot.Contains("uncorr"))){
	parameterSystematic.varToPlot.ReplaceAll("refitbs","refitbs_uncorr");
	parameterSystematic.cuts.ReplaceAll("RefitV_with_BS","RefitV_with_BS_uncorr");
	if(parameterSystematic.cuts.Contains("alphaminus"))parameterSystematic.cuts.ReplaceAll("alphaminus","alphaminus_uncorr");
      }
      else parameterSystematic.weights = param.weights + "weight_" +systematicName+"*";
      HistName = sampleName + "_" + systematicName;
    }

    TTree * tree = (TTree*)fileSample->Get(TreeName);
    //cout << TreeName <<endl;
    if (tree==NULL){ //cout << TreeName << " not found in " << fileSample << endl;
      continue;}
    TH2D * hist2D;
    TH1D * hist1D;
    double xbins[20];
    double ybins[20];
    int nbinsx = parameterSystematic.nbins;
    double widthx = (parameterSystematic.xmax-parameterSystematic.xmin)/double(nbinsx);
    for (int i=0; i<=nbinsx; ++i)
      xbins[i] = parameterSystematic.xmin + double(i)*widthx;
    int nbinsy = parameterSystematic.xDNN.size() - 1;
    //cout << parameterSystematic.hist2D << endl;
    for (int i=0; i<=nbinsy; ++i){
      ybins[i] = parameterSystematic.xDNN.at(i);
      //cout << ybins[i] << " ";
    }//cout << endl;
    //for (int i=0; i<=nbinsx; ++i){
    //  cout << xbins[i] << " ";
    //}cout << endl;
    if (parameterSystematic.hist2D) {      
      hist2D = new TH2D("hist","",
			nbinsx,xbins,nbinsy,ybins);
    }
    else if(checkPhiModulation_){
      hist1D = new TH1D("hist","",
                        nbinsx,xbins); 
    }else{
      hist1D = new TH1D("hist","",
                        nbinsy,ybins); 
    }
    //cout << "      plot : " << parameterSystematic.varToPlot << endl;
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
      if(find(SystematicsNames.begin(),SystematicsNames.end(),systematicName)==SystematicsNames.end()){
	if(!systematicName.Contains("sub_syst")&&!systematicName.Contains("syst_njets"))
	  parameterSystematic.weights = param.weights + "weight_" +systematicName+"*";
	else if(systematicName.Contains("syst_njets")){
	  TString storedweight = systematicName;
	  if(storedweight.Contains("njets0")){
	    storedweight.ReplaceAll("_njets0","");
	    parameterSystematic.weights = param.weights + "((weight_" +storedweight+"*(njets==0))+(njets!=0)*ff_mva)*";
	  }else if(storedweight.Contains("njets1")){
	    storedweight.ReplaceAll("_njets1","");
	    if(storedweight.Contains("qcd"))
	      parameterSystematic.weights = param.weights + "((weight_" +storedweight+"*(njets>=1))+(njets<1)*ff_mva)*";
	    else
	      parameterSystematic.weights = param.weights + "((weight_" +storedweight+"*(njets==1))+(njets!=1)*ff_mva)*";
	  }else if(storedweight.Contains("njets2")){
	    storedweight.ReplaceAll("_njets2","");
	    parameterSystematic.weights = param.weights + "((weight_" +storedweight+"*(njets>=2))+(njets<2)*ff_mva)*";
	  }
	}
	else if(sampleName=="jetFakes")
	  parameterSystematic.weights = param.weights; // small exception: for jetFakes_ff_mt_sub_syst we want the fake bkg to not be scaled, but just the contribution coming from genuine taus
	else if(systematicName.Contains("sub_systUp"))
	  parameterSystematic.weights = param.weights + "1.1*";
	else if(systematicName.Contains("sub_systDown"))
	  parameterSystematic.weights = param.weights + "0.9*";
      }else TreeName = "TauCheck_" + systematicName; 

    }
    TTree * tree = (TTree*)fileSample->Get(TreeName);
    TH2D * hist2D;
    TH1D * hist1D;
    double xbins[20];
    double ybins[20];
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
    else if(checkPhiModulation_){
      hist1D = new TH1D("hist","",
                        nbinsx,xbins); 
    } else {
      hist1D = new TH1D("hist","",
                        nbinsy,ybins); 
      //for(int i=0; i<=nbinsy; ++i)cout <<ybins[i]<<" ";
      //cout<<endl;
    
    }
    //cout << "      plot : " << parameterSystematic.varToPlot << endl;
    //cout << "        2D : " <<  parameterSystematic.hist2D << endl;
    //cout << tree << endl;
	//    cout << "       cut : " << parameterSystematic.cuts << endl;
    //cout << sampleName << ":      " << TreeName << "   " << (parameterSystematic.weights+"("+parameterSystematic.cuts+")") << endl;

    tree->Draw(parameterSystematic.varToPlot+">>hist",parameterSystematic.weights+"("+parameterSystematic.cuts+")");
    //cout << hist1D->Integral() <<endl;
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
      cutsSample = cutsSample +"&&"+EmbedCut_;
    else if (sample=="ZL")
      cutsSample = cutsSample +"&&!("+EmbedCut_+")";
    else if (embedded_&&sample!="EmbedZTT")
      cutsSample = cutsSample +"&&!("+EmbedCut_+")";

    if (sample=="EmbedZTT" && era_=="2016") cutsSample += "&&(weight<1000)";

    if(fakeFactor_) cutsSample += "&&gen_match_2!=6";

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
  TString weights = parameters.weights;
  //runSystematics=false;
  for (auto systematicName : FFSystematics_ ) {
    params parametersSyst = parameters;
    TString weightsSyst = weights;
    if (systematicName==""||systematicName.Contains("sub_syst")) weightsSyst = weights+weightFForQCD;
    if (!runSystematics && systematicName!="") continue;
    if (systematicName!="") cout << "     Running on systematics " << systematicName << std::endl;
    //if(!systematicName.Contains("syst")&&systematicName!="") continue;

    parametersSyst.weights = weightsSyst;

    vector<TH1D*> histFF = CreateCardsSample(FakeOrQCD,parametersSyst,systematicName);

    for (auto sample : samplesToSubtract) {
      
      TString cutsSample = cuts;

      params parametersSample = parametersSyst;    

      if (sample=="ZTT"&&embedded_) continue;
      if (sample=="EmbedZTT"&&!embedded_) continue;
      
      if (sample=="ZTT" || sample=="EmbedZTT") 
        cutsSample = cutsSample+"&&"+EmbedCut_+"&&gen_match_2!=6";
      else if (sample=="ZL"){
        cutsSample = cutsSample+"&&!("+EmbedCut_+")";
	if(fakeFactor_) cutsSample += "&&gen_match_2!=6";
      }
      else if (sample=="TTT"||sample=="ST"||sample=="W"||sample=="VVT"){
	if (fakeFactor_) cutsSample += "&&gen_match_2!=6";
        if (embedded_) cutsSample = cutsSample+"&&!("+EmbedCut_+")";
      } 
      
      if (sample=="EmbedZTT" && era_=="2016")
          cutsSample += "&&(weight<1000)";
      
      parametersSample.cuts = cutsSample;
      cout << "      Subtracting " << sample << " for data-driven estimation" << endl << "         ";
      vector<TH1D*> histos = CreateCardsSample(sample,parametersSample,systematicName);
      //cout << histFF[0]->Integral() << " -----> ";
      histFF[0]->Add(histFF[0],histos[0],1,-1);
      //cout << histFF[0]->Integral() <<endl;
      
    }
    all_histFF.push_back(histFF[0]);

  }
  return all_histFF;
  
}
  
void DataCards::RunOnCategory(TString category) {

  cout << "Running on category " << category << endl;

  vector<TH1D*> allHists; allHists.clear();

  for (auto sampleName : sampleNames) {
    if(sampleName=="EWKZ"&&embedded_) continue;
    //if(sampleName!="jetFakes")continue;
    if(sampleName=="QCD"&&fakeFactor_==true)continue;
    else if(sampleName=="jetFakes"&&fakeFactor_==false)continue;

    if (sampleName=="ZTT"&&embedded_) continue;
    if (sampleName=="EmbedZTT"&&!embedded_) continue;
    params parameters;
    TString cuts = "";
    TString weight = "weight*"; 
    //if(sampleName.Contains("qqH")&&era_=="2018") weight= "weight*(4.110e-05/xsec_lumi_weight)*";
    //if(sampleName.Contains("qqH")&&era_=="2017") weight= "weight*(3.403e-05/xsec_lumi_weight)*";

    TString acotautau = variableCP_+"_01"; //We use 00 for mupi, and 01 for murho, mu0a1 and mua1. The s
    if (category.Contains("_mupi_"))
      acotautau = variableCP_+"_00";  
    bool runSystematics = runSystematics_;

    cuts = mapCategoryCut[category] + "&&pt_1>21&&pt_2>30&&m_vis>40&&TMath::Abs(eta_1)<2.1&&puppimt_1<50&&(dmMVA_2>-1&&dmMVA_2<11)";    
    TString elepTcut = "25";
    if(era_=="2017")elepTcut = "28";
    if(era_=="2018")elepTcut = "33";
    if(ditauchannel_=="et") cuts = mapCategoryCut[category] + "&&pt_1>25&&TMath::Abs(eta_1)<2.1&&pt_2>20&&TMath::Abs(eta_2)<2.3&&os>0.5&&puppimt_1<50&&((trg_singleelectron>0.5&&pt_1>"+elepTcut+")||(trg_etaucross>0.5&&pt_1>25&&pt_2>35&&TMath::Abs(eta_2)<2.1))&&(dmMVA_2>-1&&dmMVA_2<11)&&trg_singleelectron>0.5";
    if(checkPhiModulation_) cuts += "&&m_vis<85";
    //if(ditauchannel_=="et") cuts += "&&is_SingleLepTrigger>0.5";
    TString cut_FF_AR = "&&os>0.5&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5";
    TString cut_FF_SR = "&&os>0.5&&byMediumDeepTau2017v2p1VSjet_2>0.5";
    TString cut_QCD_SS = "&&os<0.5&&byMediumDeepTau2017v2p1VSjet_2>0.5";

    TString IPCut("");
    if (applyIPcut_ && ((applyIPcutOnBkg_ && !category.Contains("sig")) || category.Contains("sig"))) IPCut = "&&"+CutIP_muon_+"&&((dmMVA_2==0&&"+CutIP_pion_+")||dmMVA_2!=0)";
    cuts += IPCut;
    
    if (sampleName=="ZTT" || sampleName=="EmbedZTT") 
      cuts = cuts+"&&"+EmbedCut_+"&&gen_match_2!=6";
    else if (sampleName=="ZL")
      cuts = cuts+"&&!("+EmbedCut_+")&&gen_match_2!=6";
    else if (sampleName=="TTT"||sampleName=="ST"||sampleName=="W"||sampleName=="VVT"){
      if(fakeFactor_) cuts += "&&gen_match_2!=6";
      if (embedded_) cuts = cuts+ "&&!("+EmbedCut_+")";
    } 
    
    if (sampleName=="EmbedZTT" && era_=="2016")
        cuts += "&&(weight<1000)";


    if (sampleName.Contains("_sm_htt125"))
      weight += "gen_sm_htt125*";
    else if (sampleName.Contains("_ps_htt125"))
      weight += "gen_ps_htt125*";
    else if (sampleName.Contains("_mm_htt125"))
      weight += "gen_mm_htt125*";


    parameters.weights = weight ; 
    parameters.cuts = cuts + cut_FF_SR; // default approach, overwriten for QCD and jetFakes samples
     

    if (category.Contains("_mupi")){
	nbins_=binsperchannel_.at("mupi");
    }
    else if(category.Contains("_murho")){
	nbins_=binsperchannel_.at("murho");
    }
    else if(category.Contains("_mua1")){
	nbins_=binsperchannel_.at("mua1");
    }
    else if(category.Contains("_mu0a1")){
	nbins_=binsperchannel_.at("mu0a1");
	xDNNSig_={0.0, 0.45, 0.6, 0.8, 1.0};
    } else nbins_= 7;

    if (category.Contains("_sig")) {
      if(!useTH1forHiggs_){
	parameters.hist2D = true;
	parameters.varToPlot = "predicted_prob:"+acotautau;
      }else{
	parameters.varToPlot = "predicted_prob";
	parameters.hist2D = false;
	xDNNSig_={0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
      }
      parameters.xDNN = xDNNSig_;
    }
    else if (category.Contains("_ztt")) {
      if(useTH2forZtt_){
	parameters.hist2D = true;
	parameters.varToPlot = "predicted_prob:"+acotautau;
      }else{
	parameters.varToPlot = "predicted_prob";
	parameters.hist2D = false;
      }
      parameters.xDNN = xDNNZtt_; 
    }
    else if (category.Contains("_fakes")) {
      if(useTH2forFakes_){
	parameters.hist2D = true;
	parameters.varToPlot = "predicted_prob:"+acotautau;
      }else{
	parameters.varToPlot = "predicted_prob";
	parameters.hist2D = false;
      }
      parameters.xDNN = xDNNFakes_; 
    }
    else if (category.Contains("alpha")){
      parameters.hist2D = false;
      parameters.varToPlot = acotautau;
    }
      

    parameters.nbins = nbins_;
    parameters.xmin = xmin_;
    parameters.xmax = xmax_;


    if (sampleName=="EmbedZTT" && runSystematics  && embedded_){ 
      vector<TH1D*> hists = CreateCardsEmbedSyst(parameters);
      for (auto hist : hists)
	allHists.push_back(hist);
    }
    if (sampleName=="jetFakes") {
      params parametersFF = parameters;
      parametersFF.cuts = cuts + cut_FF_AR;
      TString weightFF("ff_mva*");
      vector<TH1D*> hists = CreateCardsFakes(sampleName,parametersFF,weightFF, 
					      runSystematics);
      for (auto hist : hists)
	allHists.push_back(hist);
    }
    else if (sampleName=="QCD") {
      params parametersQCD = parameters;
      parametersQCD.cuts = cuts + cut_QCD_SS;
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
