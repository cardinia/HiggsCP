#include <iostream>
#include <map>
#include <algorithm>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TSystem.h"
#include "settings_for_eras_newNTuples.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"

bool applyPreselection = true;

void create_dnn_ntuples_mt( TString era = "2017" , TString channel="mt"){
  if(channel!="mt"&&channel!="et"){
    cout << "ERROR: macro only for lep tau channels, please specify et or mt" << endl;
    exit(EXIT_FAILURE);
  }
  // Some definitions
  double luminosity = 0;
  float qcd_ss_os_iso_relaxed_ratio = 0;
  float trigger_filter_efficiency = 1;
  float embedded_trigger_weight = 1.0;
  float embedded_tracking_weight = 1.0;
  TString input_dir;

  // Mapping of subsamples to output root-file
  map< TString , vector<TString> > samples_map;
  const map<TString, double>  *xsec_map    = 0;
  const map<TString, TString> *process_map = 0;

  if (era == "2018"){
     xsec_map    = &xsec_map_2018;
     process_map = &process_map_2018;
     luminosity  = 59740;
     trigger_filter_efficiency = 1.0;
     qcd_ss_os_iso_relaxed_ratio = 1.89; //number from Janek's talk in TauPOG meeting (10.04.19)
     embedded_trigger_weight  = 1.00;
     embedded_tracking_weight = 1.00;
     if(channel=="mt"){
       samples_map[channel + "-NOMINAL_ntuple_Data"   ] = SingleMuon_Run2018;
       samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedMuTau_2018;
     }else{
       samples_map[channel + "-NOMINAL_ntuple_Data"   ] = SingleElectron_Run2018;
       samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedElTau_2018;
     }
     samples_map[channel + "-NOMINAL_ntuple_DYJets"   ] = DYJets_2018;
     samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2018;
     samples_map[channel + "-NOMINAL_ntuple_TTbar"    ] = TTbar_2018;
     samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2018;
     samples_map[channel + "-NOMINAL_ntuple_Diboson"  ] = Diboson_2018;
     input_dir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu_2018_all_eras/";  
  }
  else if(era == "2017"){
    xsec_map    = &xsec_map_2017; 
    process_map = &process_map_2017; 
    luminosity  = 41900;      
    trigger_filter_efficiency = 1.0; 
    qcd_ss_os_iso_relaxed_ratio = 1.; 
    embedded_trigger_weight  = 1.00;
    embedded_tracking_weight = 0.99;
    
    if(channel=="mt"){
      samples_map[channel + "-NOMINAL_ntuple_Data"     ] = SingleMuon_Run2017; 
      //samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedMuTau_2017;
    }else{
      samples_map[channel + "-NOMINAL_ntuple_Data"     ] = SingleElectron_Run2017; 
      samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedMuTau_2017;
    }
    samples_map[channel + "-NOMINAL_ntuple_DY"       ] = DYJets_2017;  		
    samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2017;  		
    samples_map[channel + "-NOMINAL_ntuple_TT"       ] = TTbar_2017;  		
    samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2017;  	
    samples_map[channel + "-NOMINAL_ntuple_VV"       ] = Diboson_2017;  		
    samples_map[channel + "-NOMINAL_ntuple_ggH"      ] = GluGluHToTauTau_2017; 	
    samples_map[channel + "-NOMINAL_ntuple_VBF"      ] = VBFHToTauTau_2017;
    //samples_map[channel + "-NOMINAL_ntuple_CPodd"    ] = SUSYGluGluToHToTauTau_2017;
    //samples_map[channel + "-NOMINAL_ntuple_ZH"       ] = ZHToTauTau_2017; 
    //samples_map[channel + "-NOMINAL_ntuple_WH"       ] = WHToTauTau_2017; 
    //samples_map[channel + "-NOMINAL_ntuple_ggHWW"    ] = ggHToWW_2017;   
    //samples_map[channel + "-NOMINAL_ntuple_VBFHWW"   ] = VBFHToWW_2017; 
    //samples_map[channel + "-NOMINAL_ntuple_ttH"      ] = ttH_2017; 
      
    input_dir="/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/mutau/Oleg";
  }  
  else if(era == "2016"){
    xsec_map    = &xsec_map_2016;
    process_map = &process_map_2016;
    luminosity  = 35866;                   
    trigger_filter_efficiency = 0.979;
    qcd_ss_os_iso_relaxed_ratio = 2.3;
    embedded_trigger_weight  = 1.03;
    embedded_tracking_weight = 0.98;
    if(channel=="mt"){
      //samples_map[channel + "-NOMINAL_ntuple_Data"   ] = SingleMuon_Run2016;
      //samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedMuTau_2016;
    }else{
	
      //samples_map[channel + "-NOMINAL_ntuple_Data"   ] = SingleElectron_Run2016;
      //samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = EmbeddedElTau_2016;
    }
    samples_map[channel + "-NOMINAL_ntuple_DYJets"   ] = DYJets_2016;
    samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2016;
    samples_map[channel + "-NOMINAL_ntuple_TTbar"    ] = TTbar_2016;
    samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2016;
    samples_map[channel + "-NOMINAL_ntuple_Diboson"  ] = Diboson_2016;
    samples_map[channel + "-NOMINAL_ntuple_ggH"      ] = GluGluHToTauTau_2016;
    samples_map[channel + "-NOMINAL_ntuple_VBFH"     ] = VBFHToTauTau_2016;
    samples_map[channel + "-NOMINAL_ntuple_ZH"       ] = ZHToTauTau_2016;
    samples_map[channel + "-NOMINAL_ntuple_WH"       ] = WHToTauTau_2016;
    samples_map[channel + "-NOMINAL_ntuple_ggHWW"    ] = ggHToWW_2016;
    samples_map[channel + "-NOMINAL_ntuple_VBFHWW"   ] = VBFHToWW_2016;
    samples_map[channel + "-NOMINAL_ntuple_ttH"      ] = ttH_2016;
    input_dir="/nfs/dust/cms/user/tlenz/13TeV/2017/CMSSW/2016_legacy/CMSSW_8_0_29/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu_2016/NTuples/ntuples_v5/";
  
  }

  // Needed for stitching
  double xsecWIncl      = xsec_map->at(process_map->at("WJets"));
  double xsecW1Jets     = xsec_map->at(process_map->at("W1Jets"));
  double xsecW2Jets     = xsec_map->at(process_map->at("W2Jets"));
  double xsecW3Jets     = xsec_map->at(process_map->at("W3Jets"));
  double xsecW4Jets     = xsec_map->at(process_map->at("W4Jets"));
  double xsecDYIncl     = xsec_map->at(process_map->at("DYJets"));
  double xsecDY1Jets    = xsec_map->at(process_map->at("DY1Jets"));
  double xsecDY2Jets    = xsec_map->at(process_map->at("DY2Jets"));
  double xsecDY3Jets    = xsec_map->at(process_map->at("DY3Jets"));
  double xsecDY4Jets    = xsec_map->at(process_map->at("DY4Jets"));
  double neventsWIncl   = getNEventsProcessed(input_dir,process_map->at("WJets"));
  double neventsW1Jets  = getNEventsProcessed(input_dir,process_map->at("W1Jets"));
  double neventsW2Jets  = getNEventsProcessed(input_dir,process_map->at("W2Jets"));
  double neventsW3Jets  = getNEventsProcessed(input_dir,process_map->at("W3Jets"));
  double neventsW4Jets  = getNEventsProcessed(input_dir,process_map->at("W4Jets"));
  double neventsDYIncl  = getNEventsProcessed(input_dir,process_map->at("DYJets"));
  double neventsDY1Jets = getNEventsProcessed(input_dir,process_map->at("DY1Jets"));
  double neventsDY2Jets = getNEventsProcessed(input_dir,process_map->at("DY2Jets"));
  double neventsDY3Jets = getNEventsProcessed(input_dir,process_map->at("DY3Jets"));
  double neventsDY4Jets = getNEventsProcessed(input_dir,process_map->at("DY4Jets"));

  TString output_dir = "test/NTuples_"+channel+"_" + era;
  gSystem -> Exec("mkdir " + output_dir);

  TFile* ff_file = TFile::Open("/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HTTutilities/Jet2TauFakes/data/Jet2TauFakesFiles-"+era+"-SM"+era+"/SM"+era+"/tight/vloose/mt/fakeFactors.root");
  FakeFactor* ff = (FakeFactor*)ff_file->Get("ff_comb");
  ff_file->Close();
  TFile* fake_frac_file = TFile::Open("/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/data/FakeFractions_mvis-njetsbinned.root");
  TH2D* h_ff_QCD = (TH2D*)fake_frac_file->Get("ff_QCD");
  TH2D* h_ff_W   = (TH2D*)fake_frac_file->Get("ff_W");
  TH2D* h_ff_tt  = (TH2D*)fake_frac_file->Get("ff_tt");
  h_ff_QCD->SetDirectory(0);
  h_ff_W->SetDirectory(0);
  h_ff_tt->SetDirectory(0);
  fake_frac_file->Close();

  ///////////////////////////
  ///  Loop over all samples  
  ///////////////////////////
  for (auto const& sample : samples_map){
    cout << endl << sample.first << "  :  " << endl ;
   
    TFile *outFile = new TFile(output_dir + "/" + sample.first + ".root","RECREATE");
    TTree *outTree = new TTree("TauCheck", "tree created as DNN input");
    bool firstTree = true;
    TList* treeList = new TList();

    for(TString const& subsample: sample.second) {

      cout << "  - " << subsample << " : ";

      TFile *inFile  = new TFile( input_dir + "/" + subsample + ".root" ,"READ");
      TTree *inTree  = (TTree*) inFile -> Get("TauCheck");
      double nevents = getNEventsProcessed( input_dir, subsample); 
      
      // SetBranchAddress for variables that need are needed for preselection or stitching
      //variables below are for preselection
      float iso_1;   
      int extraelec_veto; 
      int extramuon_veto; 
      int dilepton_veto;  
      float pt_1;  
      float pt_2;
      
      float mva17_2;
      float mt_1;
      float againstMuonTight3_2;
      float againstElectronVLooseMVA6_2;
      float againstMuonLoose3_2;
      float againstElectronTightMVA6_2;
      bool singleLepTrigger;
      bool xTrigger;
      //end declaration preselection variables

      int gen_noutgoing;
      uint run; //only used for embedded 2016

      //Merijn: vars below used for stxs. prefiring_weight is set for different stxs bins.. this will need special attention
      int njets;
      float mjj;
      float jdeta;
      float pt_tt;
      float dijetpt;
      float jpt_1;
      float jpt_2;
      float jeta_1;
      float jeta_2;

      //Variables used for FF method
      Int_t tau_decay_mode_2;
      Float_t m_vis;
      Float_t ff_nom;
      Float_t ff_sys;
      Float_t byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;

      //branches for preselection
      inTree->SetBranchAddress("iso_1",&iso_1); 
      inTree->SetBranchAddress("extraelec_veto",&extraelec_veto);
      inTree->SetBranchAddress("extramuon_veto",&extramuon_veto);
      inTree->SetBranchAddress("dilepton_veto",&dilepton_veto);
      inTree->SetBranchAddress("pt_1",&pt_1);
      inTree->SetBranchAddress("pt_2",&pt_2);

      inTree->SetBranchAddress("mva17_2",&mva17_2);
      inTree->SetBranchAddress("mt_1",&mt_1);
      inTree->SetBranchAddress("againstMuonTight3_2",&againstMuonTight3_2);
      inTree->SetBranchAddress("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2);
      inTree->SetBranchAddress("againstMuonLoose3_2",&againstMuonLoose3_2);
      inTree->SetBranchAddress("againstElectronTightMVA6_2",&againstElectronTightMVA6_2);
      inTree->SetBranchAddress("singleLepTrigger",&singleLepTrigger);
      inTree->SetBranchAddress("xTrigger",&xTrigger);

      //branches for stichting
      inTree->SetBranchAddress("gen_noutgoing",&gen_noutgoing);
      inTree->SetBranchAddress("run",&run);      
      inTree->SetBranchAddress("njets",&njets);
      inTree->SetBranchAddress("mjj",&mjj);
      inTree->SetBranchAddress("jdeta",&jdeta);
      inTree->SetBranchAddress("pt_tt",&pt_tt);
      inTree->SetBranchAddress("dijetpt",&dijetpt);
      inTree->SetBranchAddress("jpt_1",&jpt_1);
      inTree->SetBranchAddress("jpt_2",&jpt_2);
      inTree->SetBranchAddress("jeta_1",&jeta_1);
      inTree->SetBranchAddress("jeta_2",&jeta_2);

      //branches for FF
      inTree->SetBranchAddress("m_vis",&m_vis);      
      inTree->SetBranchAddress("tau_decay_mode_2",&tau_decay_mode_2); 
      inTree->SetBranchAddress("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2",&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);      

      outFile->cd();
      TTree *currentTree = new TTree(subsample,"temporary tree");

      // New branches
      float xsec_lumi_weight;      
      float qcd_correction;
      float trigger_filter_weight;
      float embedded_stitching_weight;
      float embedded_rate_weight;
      float prefiring_weight;
      int htxs_reco_flag_ggh;
      int htxs_reco_flag_qqh;
      if(firstTree){
	outTree    = inTree->CloneTree(0);
	outTree->Branch("xsec_lumi_weight", &xsec_lumi_weight, "xsec_lumi_weight/F");
	outTree->Branch("qcd_correction", &qcd_correction, "qcd_correction/F");
	outTree->Branch("trigger_filter_weight", &trigger_filter_weight, "trigger_filter_weight/F");
	outTree->Branch("embedded_stitching_weight", &embedded_stitching_weight, "embedded_stitching_weight/F");
	outTree->Branch("embedded_rate_weight", &embedded_rate_weight, "embedded_rate_weight/F");
	outTree->Branch("prefiring_weight", &prefiring_weight, "prefiring_weight/F");
	outTree->Branch("htxs_reco_flag_ggh", &htxs_reco_flag_ggh, "htxs_reco_flag_ggh/I");
	outTree->Branch("htxs_reco_flag_qqh", &htxs_reco_flag_qqh, "htxs_reco_flag_qqh/I");
	outTree->Branch("ff_nom", &ff_nom, "ff_nom/F");
	outTree->Branch("ff_sys", &ff_sys, "ff_sys/F");
	firstTree  = false;
      }
      currentTree = inTree->CloneTree(0);
      currentTree->Branch("xsec_lumi_weight", &xsec_lumi_weight, "xsec_lumi_weight/F");
      currentTree->Branch("qcd_correction", &qcd_correction, "qcd_correction/F");
      currentTree->Branch("trigger_filter_weight", &trigger_filter_weight, "trigger_filter_weight/F");
      currentTree->Branch("embedded_stitching_weight", &embedded_stitching_weight, "embedded_stitching_weight/F");
      currentTree->Branch("embedded_rate_weight", &embedded_rate_weight, "embedded_rate_weight/F");
      currentTree->Branch("prefiring_weight", &prefiring_weight, "prefiring_weight/F");
      currentTree->Branch("htxs_reco_flag_ggh", &htxs_reco_flag_ggh, "htxs_reco_flag_ggh/I");
      currentTree->Branch("htxs_reco_flag_qqh", &htxs_reco_flag_qqh, "htxs_reco_flag_qqh/I");
      currentTree->Branch("ff_nom", &ff_nom, "ff_nom/F");
      currentTree->Branch("ff_sys", &ff_sys, "ff_sys/F");

      // lumi-xsec-weight added
      if( xsec_map->find(subsample) == xsec_map->end() && !sample.first.Contains("Data")  && !sample.first.Contains("Embedded")){
	cout << endl << endl << "Sample " << subsample << " is missing in xsec_map. Exit code." << endl << endl ;

	exit(-1);
      }
      float xsec = 1;
      if(!sample.first.Contains("Data") && !sample.first.Contains("Embedded")) xsec = xsec_map->at(subsample);
      int counter=0;
      for (int i=0; i<inTree->GetEntries(); i++) {
	//for (int i=0; i<1000; i++) {
	inTree->GetEntry(i);
	
	if(applyPreselection){
	  if(channel=="mt"){
	    if( iso_1 > 0.15 )               continue;
	    else if( pt_1 < 20 )                  continue; 
	    else if( pt_2 < 30 )                  continue; 
	    //else if( mt_1>50 )                    continue;// kept for W+jets reweighting
	    else if( againstMuonTight3_2 < 0.5 )  continue;
	    else if( againstElectronVLooseMVA6_2<0.5) continue;
	  }else{
	    if( iso_1 > 0.15 )               continue;
	    else if( pt_1 < 20 )                  continue; 
	    else if( pt_2 < 30 )                  continue; 
	    //else if( mt_1>50 )                    continue;// kept for W+jets reweighting
	    else if( againstMuonLoose3_2 < 0.5 )  continue;
	    else if( againstElectronTightMVA6_2<0.5) continue;
	  }
	  if( extraelec_veto > 0.5 )       continue;
	  else if( extramuon_veto > 0.5 )       continue;
	  else if( dilepton_veto  > 0.5 )       continue;
	  else if( byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 < 0.5) continue;
	  else if( singleLepTrigger < 0.5&&xTrigger<0.5 ) continue;
	  
	}
	ff_nom=1.;

	if(mva17_2<0.5){
	  //FF method: the input variables are stored in an array, the order MUST be (tau pt, tau decay mode, njets, visible mass, muon isolation, ff_QCD, ff_W, ff_tt)
	  int bin=h_ff_QCD->FindBin(m_vis,njets);
	  std::vector<string> inputNames( ff->inputs());
	  std::vector<double> inputs={pt_2,                                  
				      static_cast<double>(tau_decay_mode_2), 
				      static_cast<double>(njets),            
				      m_vis,                                 
				      mt_1,                                  
				      iso_1,                                 
				      h_ff_QCD -> GetBinContent(bin),
				      h_ff_W   -> GetBinContent(bin),
				      h_ff_tt  -> GetBinContent(bin)};

	    if(abs(1-(inputs[6]+inputs[7]+inputs[8]))>1e-5){
	    cout << "ERROR: " << inputs[6]+inputs[7]+inputs[8] << "!=1 : ";
	    for(auto input: inputs) cout<< input << " ";
	    cout << endl;
	  }
	  ff_nom = ff->value(inputs);
	}
	ff_sys = ff_nom; // TO DO: fix systematics

	xsec_lumi_weight = xsec*luminosity/nevents;
	qcd_correction = qcd_ss_os_iso_relaxed_ratio;
	trigger_filter_weight = trigger_filter_efficiency;

	// Initialize variables for jet observables to -10 for NN
	if(njets < 2){
	  jdeta   = -10;
	  mjj     = -10;
	  dijetpt = -10;
	  //pt_ttjj = -10; 
	  jpt_2   = -10;
	  jeta_2  = -10;
	  if(njets < 1){
	    jpt_1 = -10;
	    jeta_1= -10;
	  }
	}
	
	// Stitching only for wjets MC in n-jet binned samples in gen_noutgoing
	if( subsample.Contains("W") && subsample.Contains("JetsToLNu") ){
	  if(gen_noutgoing == 1)      xsec_lumi_weight = luminosity / ( neventsW1Jets/xsecW1Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing == 2) xsec_lumi_weight = luminosity / ( neventsW2Jets/xsecW2Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing == 3) xsec_lumi_weight = luminosity / ( neventsW3Jets/xsecW3Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing == 4) xsec_lumi_weight = luminosity / ( neventsW4Jets/xsecW4Jets + neventsWIncl/xsecWIncl );
	  else                   xsec_lumi_weight = luminosity / ( neventsWIncl/xsecWIncl );
	}
	else if( subsample.Contains("DY") && subsample.Contains("JetsToLL_M-50") ){
	  //	  cout<<"Number partons: "<<gen_noutgoing<<endl;
	  if(gen_noutgoing == 1)      xsec_lumi_weight = luminosity / ( neventsDY1Jets/xsecDY1Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing == 2) xsec_lumi_weight = luminosity / ( neventsDY2Jets/xsecDY2Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing == 3) xsec_lumi_weight = luminosity / ( neventsDY3Jets/xsecDY3Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing == 4) xsec_lumi_weight = luminosity / ( neventsDY4Jets/xsecDY4Jets + neventsDYIncl/xsecDYIncl );
	  else                   xsec_lumi_weight = luminosity / ( neventsDYIncl/xsecDYIncl );
	}
	
	if( sample.first.Contains("Data") || sample.first.Contains("Embedded")){
	  xsec_lumi_weight = 1.;
	  trigger_filter_weight = 1.;
	}
	if( sample.first.Contains("Embedded") && era == "2016"){
	  embedded_stitching_weight = 
	    ((run >= 272007) && (run < 275657))*(1.0/0.891)
	    +((run >= 275657) && (run < 276315))*(1.0/0.910)
	    +((run >= 276315) && (run < 276831))*(1.0/0.953)
	    +((run >= 276831) && (run < 277772))*(1.0/0.947)
	    +((run >= 277772) && (run < 278820))*(1.0/0.942)
	    +((run >= 278820) && (run < 280919))*(1.0/0.906)
	    +((run >= 280919) && (run < 284045))*(1.0/0.950);
	}
	else embedded_stitching_weight = 1.;

	embedded_rate_weight = embedded_trigger_weight * embedded_tracking_weight;

	// add flags for cut categories which correspond to htxs_stage1cats
	htxs_reco_flag_ggh = 0;
	htxs_reco_flag_qqh = 0;
	// reco bins 101 and 102 merged into 2jets categories (sicne hardly populated)
	if(njets==0)                                    htxs_reco_flag_ggh = 103;
	else if((njets==1) && (pt_tt>0)  &&(pt_tt<60))  htxs_reco_flag_ggh = 104;
	else if((njets==1) && (pt_tt>60) &&(pt_tt<120)) htxs_reco_flag_ggh = 105;
	else if((njets==1) && (pt_tt>120)&&(pt_tt<200)) htxs_reco_flag_ggh = 106;
	else if((njets==1) && (pt_tt>200))              htxs_reco_flag_ggh = 107;
	else if((njets>=2) && (pt_tt>0)  &&(pt_tt<60))  htxs_reco_flag_ggh = 108;
	else if((njets>=2) && (pt_tt>60) &&(pt_tt<120)) htxs_reco_flag_ggh = 109;
	else if((njets>=2) && (pt_tt>120)&&(pt_tt<200)) htxs_reco_flag_ggh = 110;
	else if((njets>=2) && (pt_tt>200))              htxs_reco_flag_ggh = 111;

	//stxs categorization currently not used, kept commented
	/*
	if((jpt_1>0)&&(jpt_1<200)&&(njets>=2)&&(mjj>400)&&(jdeta>2.8)&&(pt_ttjj>0)&&(pt_ttjj<25)) htxs_reco_flag_qqh = 201;
	else if((jpt_1>0)&&(jpt_1<200)&&(njets>=2)&&(mjj>400)&&(jdeta>2.8)&&(pt_ttjj>25)) htxs_reco_flag_qqh = 202;
	else if((jpt_1>0)&&(jpt_1<200)&&(njets>=2)&&(mjj>60)&&(mjj<120)) htxs_reco_flag_qqh = 203;
	else if(( (jpt_1>0&&jpt_1<200&&njets<2) || (jpt_1>0&&jpt_1<200&&njets>=2&&mjj>400&&jdeta<2.8) || (jpt_1>0&&jpt_1<200&&njets>=2&&mjj>0&&mjj<60) || (jpt_1>0&&jpt_1<200&&njets>=2&&mjj>120&&mjj<400))) htxs_reco_flag_qqh = 204;
	else if(jpt_1>200) htxs_reco_flag_qqh = 205;*/

	// prefiring weights (from AN-18-255)
	prefiring_weight=1;
	if( sample.first.Contains("TTBar") && era == "2016")      prefiring_weight = 0.989;
	else if( sample.first.Contains("TTBar") && era == "2017") prefiring_weight = 0.984;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 201 && era == "2016") prefiring_weight = 0.972;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 201 && era == "2017") prefiring_weight = 0.950;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 202 && era == "2016") prefiring_weight = 0.972;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 202 && era == "2017") prefiring_weight = 0.950;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 203 && era == "2016") prefiring_weight = 0.972;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 203 && era == "2017") prefiring_weight = 0.950;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 204 && era == "2016") prefiring_weight = 0.983;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 204 && era == "2017") prefiring_weight = 0.970;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 205 && era == "2016") prefiring_weight = 0.920;
	else if( sample.first.Contains("VBFH") && htxs_reco_flag_qqh == 205 && era == "2017") prefiring_weight = 0.850;

	//stxs categorization currently not used, kept commented
	/*
	// Select hadronic and leptonic part of VH sample
	if( subsample.Contains("VH") || subsample.Contains("WplusH") || subsample.Contains("WminusH") ){
	  if( sample.first.Contains("VBFH") && (htxs_stage1cat>206 || htxs_stage1cat<200) ) continue;
	  if( (sample.first.Contains("WH") || sample.first.Contains("ZH")) && htxs_stage1cat<=206 && htxs_stage1cat>=200 ) continue;
	  }*/

	currentTree->Fill();
      }
      cout<<currentTree->GetEntries()<<endl;
      treeList->Add(currentTree);
      inFile->Close();

      //fill xsec_lumi_weight in a histogram for reference, and write to the tree
      TH1F* xsec_lumi_weightH=new TH1F("xsec_lumi_weightH","xsec_lumi_weightH",1,0,1);
      xsec_lumi_weightH->Fill(0.5,xsec_lumi_weight);
      xsec_lumi_weightH->Write();
      delete xsec_lumi_weightH;
      
    }
    cout<<"Before merging trees"<<endl;
    outTree = TTree::MergeTrees(treeList); 
    cout<<"Together : "<<outTree->GetEntries()<<endl;
    outTree  -> Write( "" , TObject::kOverwrite );
    treeList -> Delete();
    cout << outFile << endl;
    outFile  -> Close();
  }
  cout << endl; 
}
