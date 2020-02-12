map<TString, double> xsecs = {
  {"WJetsToLNu"  , 52760*1.166}, // NNLO (1)
  {"W1JetsToLNu" , 1.166*8104.}, // NNLO (2)
  {"W2JetsToLNu" , 1.166*2796.}, // NNLO (3)
  {"W3JetsToLNu" , 1.166*993.5}, // NNLO (4)
  {"W4JetsToLNu" , 1.166*544.4}, // NNLO (5)
  {"W1JetsToLNu_LHEWpT_50-150"  , 1.0181*2661},  // NNLO (6)
  {"W1JetsToLNu_LHEWpT_100-150" , 1.05662*1.0181*286.1}, // NNLO (6a)
  {"W1JetsToLNu_LHEWpT_150-250" , 0.91087*1.0181*71.9},  // NNLO (7)
  {"W1JetsToLNu_LHEWpT_250-400" , 0.69395*1.0181*8.05},  // NNLO (8)
  {"W1JetsToLNu_LHEWpT_400-inf" , 0.64493*1.0181*0.885}, // NNLO (9)
  {"ZZ" , 12.19},  // LO (17) -> could be improved
  {"WW" , 118.7},  // NNLO QCD (18)
  {"WZ" , 27.68},  // LO (19) -> could be improved
  {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
  {"DY1JetsToLL_M-50"      , 878.7*1.079}, // NNLO (20a)
  {"DY2JetsToLL_M-50"      , 304.4*1.079}, // NNLO (20b)
  {"DY3JetsToLL_M-50"      , 111.5*1.079}, // NNLO (20c)
  {"DY4JetsToLL_M-50"      , 44.03*1.079}, // NNLO (20d)
  {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
  {"TTToHadronic"     , 377.96}, // NNLO (22)
  {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
  {"ST_t-channel_top_4f"     , 136.02}, // ? (24) -> could be improved
  {"ST_t-channel_antitop_4f" , 80.95}, // ? (25) -> could be improved
  {"ST_tW_top_5f"            , 35.85}, // ? (26) -> could be improved
  {"ST_tW_antitop_5f"        , 35.85}, // ? (27) -> could be improved
};

map<TString, double> n_events_per_sample = {
  {"WJetsToLNu"  , 74635450},
  {"W1JetsToLNu" , 54988117},
  {"W2JetsToLNu" , 32368249},
  {"W3JetsToLNu" , 19700377},
  {"W4JetsToLNu" , 11333705},
  {"W1JetsToLNu_LHEWpT_100-150" , 61651108},
  {"W1JetsToLNu_LHEWpT_150-250" , 126166905},
  {"W1JetsToLNu_LHEWpT_250-400" , 20841145},
  {"W1JetsToLNu_LHEWpT_400-inf" , 4465538},
  {"ZZ" , 1949768},
  {"WW" , 7791498},
  {"WZ" , 3928630},
  {"TTTo2L2Nu"        , 5.05414e+09},
  {"TTToHadronic"     , 4.09177e+10},
  {"TTToSemiLeptonic" , 4.54904e+10},
  {"ST_t-channel_top_4f"     , 5.98206e+06},
  {"ST_t-channel_antitop_4f" , 3.67591e+06},
  {"ST_tW_top_5f"            , 2.72081e+08},
  {"ST_tW_antitop_5f"        , 2.79005e+08},
  {"DYJetsToLL_M-50"     , 97800939},
  {"DY1JetsToLL_M-50"    , 34859434},
  {"DY2JetsToLL_M-50"    , 9790490},
  {"DY3JetsToLL_M-50"    , 6897933},
  {"DY4JetsToLL_M-50"    , 4346952}
};
