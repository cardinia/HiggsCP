#ifndef SETTINGS_FOR_ERAS_H
#define SETTINGS_FOR_ERAS_H

// **************************************************************************************************
// Define the subsamples that belong to a certain process
// 2018
const vector<TString> SingleMuon_Run2018       = { "SingleMuon_Run2018A",
						   "SingleMuon_Run2018B", 
						   "SingleMuon_Run2018C",
						   "SingleMuon_Run2018D"};

const vector<TString> SingleElectron_Run2018       = { "SingleElectron_Run2018A",
						       "SingleElectron_Run2018B", 
						       "SingleElectron_Run2018C",
						       "SingleElectron_Run2018D"};
const vector<TString> EmbeddedMuTau_2018        = { "EmbeddingMuTau_Run2018A",
						    "EmbeddingMuTau_Run2018B",
						    "EmbeddingMuTau_Run2018C",
						    "EmbeddingMuTau_Run2018D"};
const vector<TString> EmbeddedElTau_2018        = { "EmbeddingElTau_Run2018A",
						    "EmbeddingElTau_Run2018B",
						    "EmbeddingElTau_Run2018C",
						    "EmbeddingElTau_Run2018D"};
const vector<TString> DYJets_2018          = { "DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" ,
                                               "DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" ,
                                               "DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" ,
                                               "DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" ,
                                               "DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"};
const vector<TString> WJets_2018           = { "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" ,
                                               "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" ,
                                               "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" ,
                                               "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" ,
                                               "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" };
const vector<TString> TTbar_2018           = { "TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8" ,
                                               "TTToHadronic_TuneCP5_13TeV_powheg_pythia8" ,
                                               "TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8" };
const vector<TString> SingleTop_2018       = { "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8" ,
                                               "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8" ,
                                               "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8" ,
                                               "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8" };
const vector<TString> Diboson_2018         = { "WW_TuneCP5_13TeV-pythia8" ,
                                               "WZ_TuneCP5_13TeV-pythia8" ,
                                               "ZZ_TuneCP5_13TeV-pythia8" };

// 2017
const vector<TString> SingleMuon_Run2017       = {"SingleMuon_Run2017B",
						  "SingleMuon_Run2017C",
						  "SingleMuon_Run2017D",
						  "SingleMuon_Run2017E",
						  "SingleMuon_Run2017F"};
const vector<TString> SingleElectron_Run2017       = {"SingleElectron_Run2017B",
						      "SingleElectron_Run2017C",
						      "SingleElectron_Run2017D",
						      "SingleElectron_Run2017E",
						      "SingleElectron_Run2017F"};
const vector<TString> EmbeddedMuTau_2017        = { "EmbeddingMuTau_Run2017B",
						    "EmbeddingMuTau_Run2017C",
						    "EmbeddingMuTau_Run2017D",
						    "EmbeddingMuTau_Run2017E",
						    "EmbeddingMuTau_Run2017F"};
const vector<TString> EmbeddedElTau_2017        = { "EmbeddingElTau_Run2017B",
						    "EmbeddingElTau_Run2017C",
						    "EmbeddingElTau_Run2017D",
						    "EmbeddingElTau_Run2017E",
						    "EmbeddingElTau_Run2017F"};
const vector<TString> Embedded_2017        = { "Embedding_Run2017" };
const vector<TString> DYJets_2017          = { "DY1JetsToLL_M-50" ,
					       "DY2JetsToLL_M-50" ,
					       "DY3JetsToLL_M-50" ,
					       "DY4JetsToLL_M-50" ,
					       "DYJetsToLL_M-50"/*,
					       "EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8",
					       "EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8"*/};
const vector<TString> WJets_2017           = { "W1JetsToLNu" ,
					       "W2JetsToLNu" ,
					       "W3JetsToLNu" ,
					       "W4JetsToLNu" ,
					       "WJetsToLNu"
					       /*"WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
					       "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" ,
					       "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8"*/};
const vector<TString> TTbar_2017           = { "TTTo2L2Nu" ,
					       "TTToHadronic" ,
					       "TTToSemiLeptonic" };
const vector<TString> SingleTop_2017       = { "ST_t-channel_antitop_4f" ,
					       "ST_t-channel_top_4f" ,
					       "ST_tW_antitop_5f" ,
					       "ST_tW_top_5f" };
const vector<TString> Diboson_2017         = { "WW" ,
					       "WZ" ,
					       "ZZ" };
const vector<TString> GluGluHToTauTau_2017 = { "GluGluHToTauTau_M125" };
const vector<TString> SUSYGluGluToHToTauTau_2017 = { "SUSYGluGluToHToTauTau" };
const vector<TString> VBFHToTauTau_2017    = { "VBFHToTauTau_M125"
					       /* "ZHToTauTau_M125_13TeV_powheg_pythia8",
					       "WplusHToTauTau_M125_13TeV_powheg_pythia8",
					       "WminusHToTauTau_M125_13TeV_powheg_pythia8"*/}; 

//the 2017 samples below seem not needed for mu-tau
const vector<TString> ZHToTauTau_2017      = { "ZHToTauTau_M125" };//seems absent for mutau
const vector<TString> WHToTauTau_2017      = { "WplusHToTauTau_M125" , "WminusHToTauTau_M125" };//seems absent for mutau
const vector<TString> ggHToWW_2017         = { "GluGluHToWWTo2L2Nu" }; //seems absent for mu-tau
const vector<TString> VBFHToWW_2017        = { "VBFHToWWTo2L2Nu_M125" };//seems absent for mu-tau
const vector<TString> ttH_2017             = { "ttHToTauTau_M125" };//seems absent for mu-tau

// 2016
const vector<TString> SingleMuon_Run2016       = { "SingleMuon_Run2016B" ,
						   "SingleMuon_Run2016C" ,
						   "SingleMuon_Run2016D" ,
						   "SingleMuon_Run2016E" ,
						   "SingleMuon_Run2016F" ,
						   "SingleMuon_Run2016G" ,
						   "SingleMuon_Run2016H" };

const vector<TString> SingleElectron_Run2016       = { "SingleElectron_Run2016B" ,
						       "SingleElectron_Run2016C" ,
						       "SingleElectron_Run2016D" ,
						       "SingleElectron_Run2016E" ,
						       "SingleElectron_Run2016F" ,
						       "SingleElectron_Run2016G" ,
						       "SingleElectron_Run2016H" };

const vector<TString> EmbeddedElTau_2016       = { "EmbeddedElTau_Run2016B" ,
						   "EmbeddedElTau_Run2016C" ,
						   "EmbeddedElTau_Run2016D" ,
						   "EmbeddedElTau_Run2016E" ,
						   "EmbeddedElTau_Run2016F" ,
						   "EmbeddedElTau_Run2016G" ,
						   "EmbeddedElTau_Run2016H" };
const vector<TString> EmbeddedMuTau_2016       = { "EmbeddedMuTau_Run2016B" ,
						   "EmbeddedMuTau_Run2016C" ,
						   "EmbeddedMuTau_Run2016D" ,
						   "EmbeddedMuTau_Run2016E" ,
						   "EmbeddedMuTau_Run2016F" ,
						   "EmbeddedMuTau_Run2016G" ,
						   "EmbeddedMuTau_Run2016H" };
const vector<TString> Embedded_2016        = { "Embedding_Run2016" };
const vector<TString> DYJets_2016          = { "DY1JetsToLL_M-50" , "DY2JetsToLL_M-50" , "DY3JetsToLL_M-50" , "DY4JetsToLL_M-50" , "DYJetsToLL_M-50" , "DYJetsToLL_M-10to50" , "EWKZ2Jets" };
const vector<TString> WJets_2016           = { "W1JetsToLNu" , "W2JetsToLNu" , "W3JetsToLNu" , "W4JetsToLNu" , "WJetsToLNu" , "WGToLNuG" , "WGstarToLNuEE" , "WGstarToLNuMuMu" , "EWKWPlus2Jets" , "EWKWMinus2Jet" };
const vector<TString> TTbar_2016           = { "TTbar" };
const vector<TString> SingleTop_2016       = { "ST_t-channel_antitop" , "ST_t-channel_top" , "ST_tW_antitop" , "ST_tW_top" };
const vector<TString> Diboson_2016         = { "VVTo2L2Nu" , "WZJToLLLNu" , "WZTo1L1Nu2Q" , "WZTo1L3Nu" , "WZTo2L2Q" , "ZZTo2L2Q" , "ZZTo4L" , "WWToLNuQQ" };
const vector<TString> GluGluHToTauTau_2016 = { "GluGluHToTauTau_M125" };
const vector<TString> VBFHToTauTau_2016    = { "VBFHToTauTau_M125" , "ZHToTauTau_M125" , "WplusHToTauTau_M125" , "WminusHToTauTau_M125" };
const vector<TString> ZHToTauTau_2016      = { "ZHToTauTau_M125" };
const vector<TString> WHToTauTau_2016      = { "WplusHToTauTau_M125" , "WminusHToTauTau_M125" };
const vector<TString> ggHToWW_2016         = { "GluGluHToWWTo2L2Nu_M125" };
const vector<TString> VBFHToWW_2016        = { "VBFHToWWTo2L2Nu_M125" };
const vector<TString> ttH_2016             = { "ttHJetToTT_M125" };

// **************************************************************************************************
// **************************************************************************************************
// Cross-section map
// 2018 (taken from 2018, have to be checked!)
const map<TString, double> xsec_map_2018 = {
   { "DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"  , 5765.4*1.079 }, // 5765.4 is the old NNLO xsec - it has been updated to 6225.42 - has been decided to keep the old value for consistency
   { "DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 877.8*1.079 },
   { "DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 304.4*1.079 },
   { "DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 111.5*1.079 },
   { "DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 44.03*1.079 },
   { "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"  , 52940.0*1.162 },
   { "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 8104.0*1.162 },
   { "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 2793.0*1.162 },
   { "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 992.5*1.162 },
   { "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 544.3*1.162 },
   { "TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8"        , 88.29 },
   { "TTToHadronic_TuneCP5_13TeV_powheg_pythia8"     , 377.96 },
   { "TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8" , 365.35 },
   { "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8" , 80.95 },
   { "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"     , 136.02 },
   { "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                  , 35.85 },
   { "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                      , 35.85 },
   { "WW_TuneCP5_13TeV-pythia8" , 75.88 },
   { "WZ_TuneCP5_13TeV-pythia8" , 27.57 },
   { "ZZ_TuneCP5_13TeV-pythia8" , 12.14 }
};
// 2017 (checked ! - reference is https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2017 )

//Merijn: List below is cross checked.
//2019 5 6: There is one DY sample that we have that currently does not appear, plus the susy sample. Rest cross checked, all seems consistent.
const map<TString, double> xsec_map_2017 = {
  { "DYJetsToLL_M-50"  , 5765.4*1.079 }, // 5765.4 is the old NNLO xsec - it has been updated to 6225.42 - has been decided to keep the old value for consistency
  //=> keep
  { "DY1JetsToLL_M-50" , 877.8*1.079 }, //keep
  { "DY2JetsToLL_M-50" , 304.4*1.079 }, //keep
  { "DY3JetsToLL_M-50" , 111.5*1.079 }, //keep
  { "DY4JetsToLL_M-50" , 44.03*1.079 }, //keep
  { "WJetsToLNu"  , 52940.0*1.162 }, //keep
  { "W1JetsToLNu" , 8104.0*1.162 }, //keep
  { "W2JetsToLNu" , 2793.0*1.162 }, //keep
  { "W3JetsToLNu" , 992.5*1.162 },//keep
  { "W4JetsToLNu" , 544.3*1.162 },//keep
  { "TTTo2L2Nu"        , 88.29 }, //keep
  { "TTToHadronic"     , 377.96 },//keep
  { "TTToSemiLeptonic" , 365.35 }, //keep
  { "ST_t-channel_antitop_4f" , 80.95 },//keep
  { "ST_t-channel_top_4f"     , 136.02 },//keep
  { "ST_tW_antitop_5f"                  , 35.85 },//keep
  { "ST_tW_top_5f"                      , 35.85 },//keep
  { "WW" , 75.88 },//keep
  { "WZ" , 27.57 },//keep
  { "ZZ" , 12.14 }, //keep
  { "GluGluHToTauTau_M125" , 48.58*0.0627 },//keep
  { "SUSYGluGluToHToTauTau" , 48.58*0.0627 },//keep
  { "VBFHToTauTau_M125"    , 3.782*0.0627 },//keep. Below are things that seem not needed for mt analysis
  { "EWKWMinus2Jets_WToLNu_M-50" , 23.24 },
  { "EWKWPlus2Jets_WToLNu_M-50" , 29.59 },
  { "WGToLNuG"    , 464.4 },
  { "EWKZ2Jets_ZToLL_M-50" , 4.321 },
  { "EWKZ2Jets_ZToNuNu" , 10.66 },
  { "ZHToTauTau_M125"      , 0.0594 },
  { "WplusHToTauTau_M125"  , 0.0527 },
  { "WminusHToTauTau_M125" , 0.0358 },
  { "GluGluHToWWTo2L2Nu"                        , 48.6*0.02374 },  // from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG#Higgs_cross_sections_and_decay_b and https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR#H_llll_ll
  { "VBFHToWWTo2L2Nu_M125" , 3.78*0.02374 }, // from: see links above
  { "ttHToTauTau_M125"         , 0.5071*0.0627 }  // from: see links above
};


// 2016 (taken from AN2016_355_v10 with minor unrelevant deviations - everything was checked)
const map<TString, double> xsec_map_2016 = {
    { "DYJetsToLL_M-10to50"      , 18610 },
    { "DYJetsToLL_M-50"          , 5765 },
    { "DY1JetsToLL_M-50"         , 1.164*1012.5 },
    { "DY2JetsToLL_M-50"         , 1.164*332.8 },
    { "DY3JetsToLL_M-50"         , 1.164*101.8 },
    { "DY4JetsToLL_M-50"         , 1.164*54.8 },
    { "WJetsToLNu"               , 61526.7 },
    { "W1JetsToLNu"              , 1.221*9644.5 },
    { "W2JetsToLNu"              , 1.221*3144.5 },
    { "W3JetsToLNu"              , 1.221*954.8 },
    { "W4JetsToLNu"              , 1.221*485.6 },
    { "TTbar"                    , 831.76 },
    { "ST_t-channel_antitop"     , 80.95 },
    { "ST_t-channel_top"         , 136.02 },
    { "ST_tW_antitop"            , 35.6 },
    { "ST_tW_top"                , 35.6 },
    { "VVTo2L2Nu"                , 11.95 },
    { "WWToLNuQQ"                , 49.997 },
    { "WZTo2L2Q"                 , 5.595 },
    { "WZTo1L1Nu2Q"              , 10.71 },
    { "WZTo1L3Nu"                , 3.05 },
    { "WZJToLLLNu"               , 5.26 },
    { "ZZTo4L"                   , 1.212 },
    { "ZZTo2L2Q"                 , 3.22 },
    { "WGToLNuG"                 , 178.4 }, // xsdb
    { "WGstarToLNuMuMu"          , 2.793 },
    { "WGstarToLNuEE"            , 3.526 },
    { "EWKWPlus2Jets"            , 25.62 },
    { "EWKWMinus2Jet"            , 20.20 },
    { "EWKZ2Jets"                , 3.987 },
    { "GluGluHToTauTau_M125"     , 48.58*0.0627 },
    { "VBFHToTauTau_M125"        , 3.782*0.0627 },
    { "ZHToTauTau_M125"          , 0.0594 },
    { "WplusHToTauTau_M125"      , 0.0527 },
    { "WminusHToTauTau_M125"     , 0.0358 },
    { "GluGluHToWWTo2L2Nu_M125"  , 48.6*0.02374 }, // see 2017
    { "VBFHToWWTo2L2Nu_M125"     , 3.78*0.02374 }, // see 2017
    { "ttHJetToTT_M125"          , 0.5071*0.0627 } // see 2017
  };
// **************************************************************************************************
// **************************************************************************************************
// maps short name to filenames - needed for stitching
const map<TString , TString> process_map_2018 = {
  { "WJets"   , "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"},
  { "W1Jets"  , "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"},
  { "W2Jets"  , "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"},
  { "W3Jets"  , "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"},
  { "W4Jets"  , "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"},
  { "DYJets"  , "DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"},
  { "DY1Jets" , "DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"},
  { "DY2Jets" , "DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"},
  { "DY3Jets" , "DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"},
  { "DY4Jets" , "DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"},
};


const map<TString , TString> process_map_2017 = {
  { "WJets"   , "WJetsToLNu"}, //all present in dir klundert
  { "W1Jets"  , "W1JetsToLNu"},
  { "W2Jets"  , "W2JetsToLNu"},
  { "W3Jets"  , "W3JetsToLNu"},
  { "W4Jets"  , "W4JetsToLNu"},
  { "DYJets"  , "DYJetsToLL_M-50"},
  { "DY1Jets" , "DY1JetsToLL_M-50"},
  { "DY2Jets" , "DY2JetsToLL_M-50"},
  { "DY3Jets" , "DY3JetsToLL_M-50"},
  { "DY4Jets" , "DY4JetsToLL_M-50"},
};

const map<TString , TString> process_map_2016 = {
  { "WJets"   , "WJetsToLNu"},
  { "W1Jets"  , "W1JetsToLNu"},
  { "W2Jets"  , "W2JetsToLNu"},
  { "W3Jets"  , "W3JetsToLNu"},
  { "W4Jets"  , "W4JetsToLNu"},
  { "DYJets"  , "DYJetsToLL_M-50"},
  { "DY1Jets" , "DY1JetsToLL_M-50"},
  { "DY2Jets" , "DY2JetsToLL_M-50"},
  { "DY3Jets" , "DY3JetsToLL_M-50"},
  { "DY4Jets" , "DY4JetsToLL_M-50"},
};

const map<TString, int> n_events_per_sample = {
  {"WJetsToLNu"  , 74635450},
  {"W1JetsToLNu" , 54988117},
  {"W2JetsToLNu" , 32368249},
  {"W3JetsToLNu" , 19700377},
  {"W4JetsToLNu" , 11333705},
  {"W1JetsToLNu_LHEWpT_100-150" , 61651108},
  {"W1JetsToLNu_LHEWpT_150-250" , 126166905},
  {"W1JetsToLNu_LHEWpT_250-400" , 20841145},
  {"W1JetsToLNu_LHEWpT_400-inf" , 4465538},
  {"W2JetsToLNu_LHEWpT_100-150" , 36719823},
  {"W2JetsToLNu_LHEWpT_150-250" , 170437978},
  {"W2JetsToLNu_LHEWpT_250-400" , 81630276},
  {"W2JetsToLNu_LHEWpT_400-inf" , 30340311},
  {"ZJetsToNuNu_HT-100To200"    , 22737266},
  {"ZJetsToNuNu_HT-200To400"    , 21675916},
  {"ZJetsToNuNu_HT-400To600"    , 9134120},
  {"ZJetsToNuNu_HT-600To800"    , 5697594},
  {"ZJetsToNuNu_HT-800To1200"   , 2058077},
  {"ZJetsToNuNu_HT-1200To2500"  , 338948},
  {"ZJetsToNuNu_HT-2500ToInf"   , 6734},
  {"ZZ" , 1949768},
  {"WW" , 7791498},
  {"WZ" , 3928630},
  {"DYJetsToLL_M-10to50" , 39521230},
  {"DYJetsToLL_M-50"     , 97800939},
  {"DY1JetsToLL_M-50"    , 34859434},
  {"DY2JetsToLL_M-50"    , 9790490},
  {"DY3JetsToLL_M-50"    , 6897933},
  {"DY4JetsToLL_M-50"    , 4346952},
  {"TTTo2L2Nu"        , 50669657999},//960752},
  {"TTToHadronic"     , 41729120},
  {"TTToSemiLeptonic" , 41221873},
  {"ST_t-channel_top_4f"     , 5982064},
  {"ST_t-channel_antitop_4f" , 3675910},
  {"ST_tW_top_5f"            , 7794186},
  {"ST_tW_antitop_5f"        , 7977430},
  {"WToTauNu_M-200"         , 2000000},
  {"WToTauNu_M-200_jesUp"   , 2000000},
  {"WToTauNu_M-200_jesDown" , 2000000},
  {"WToTauNu_M-200_taues_1prong0pizerosUp"       , 2000000},
  {"WToTauNu_M-200_taues_1prong0pizerosDown"     , 2000000},
  {"WToTauNu_M-200_taues_1prongUpTo4pizerosUp"   , 2000000},
  {"WToTauNu_M-200_taues_1prongUpTo4pizerosDown" , 2000000},
  {"WToTauNu_M-200_taues_3prong0pizerosUp"       , 2000000},
  {"WToTauNu_M-200_taues_3prong0pizerosDown"     , 2000000},
  {"WToTauNu_M-200_uesUp"   , 2000000},
  {"WToTauNu_M-200_uesDown" , 2000000},
  {"WToMuNu_M-200"          , 1953601},
  {"WToMuNu_M-200_jesUp"    , 1953601},
  {"WToMuNu_M-200_jesDown"  , 1953601},
  {"WToMuNu_M-200_muUp"     , 1953601},
  {"WToMuNu_M-200_muDown"   , 1953601},
  {"WToMuNu_M-200_uesUp"    , 1953601},
  {"WToMuNu_M-200_uesDown"  , 1953601},
  {"GluGluHToTauTau_M125", 121801800},
  {"VBFHToTauTau_M125", 2977152},
};

template<typename KeyTemp, typename ValTemp>
std::vector<KeyTemp> extract_keys(std::map<KeyTemp, ValTemp> const& input_map) {
  std::vector<KeyTemp> keys;
  for (auto const& element : input_map) {
    keys.push_back(element.first);
  }
  return keys;
}

double getNEventsProcessed(TString input_dir, TString sample)
{
  double nevents;
  vector<TString> keys_list = extract_keys(n_events_per_sample);
  if(std::find(keys_list.begin(), keys_list.end(), sample) != keys_list.end()){
    nevents = n_events_per_sample.at(sample);
  }else if(!sample.Contains("Run")){
    TFile * file = new TFile(input_dir+"/"+sample+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("nWeightedEvents");
    if(!histWeightsH){
      cout << endl << endl << "Histogram nWeightedEvents doesn't exist in the file "<< input_dir+"/"+sample+".root" <<". Quit program." << endl << endl;
      exit(-1);
    }
    nevents = histWeightsH->GetSumOfWeights();
      cout << "WARNING: normalization taken from nWeightedEvents, please check! ";
    file -> Close();
    delete file;
  }else nevents=0;

  return nevents;
}



#endif
