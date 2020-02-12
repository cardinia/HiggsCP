#include "HiggsCP/Inputs/interface/DataCards.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include <iostream>
#include "TString.h"
#include "TMath.h"

using namespace std;

int main(int argc, char * argv[]) {

  Config cfg(argv[1]);

  TString channel="";
  int classIndex=-1;
  if(argc>=3){
    classIndex = std::atoi(argv[2]);
    channel = argv[3];
  }

  const string era = cfg.get<string>("Era");
  TString Era(era);
  const bool embedded = cfg.get<bool>("IsEmbedded");

  const int nbinsPhiCP = cfg.get<int>("NbinsPhiCP");
  const vector<double> DNNbins = cfg.get< vector<double> >("DNNbins");

  const bool useTH2forZtt = cfg.get<bool>("UseTH2ForZtt");
  const bool mvaDM = cfg.get<bool>("mvaDM");
  const bool applyIPcut = cfg.get<bool>("ApplyIPcut");
  const bool runSystematics = cfg.get<bool>("RunSystematics");

  const string input_dir = cfg.get<string>("InputDirectory");
  TString Input_dir(input_dir);

  const string output_dir = cfg.get<string>("OutputDirectory");
  TString Output_dir(output_dir);

  const string output_filename = cfg.get<string>("OutputFileName");
  TString Output_filename(output_filename);

  /*
  for (auto XDNN : DNNbins) 
    cout << " " << XDNN;
  cout << endl;
  exit(1);
  */

  int nbins = nbinsPhiCP;
  double xmin = 0.0;
  double xmax = 2*TMath::Pi();
  vector<double> xDNN = DNNbins;
  
  DataCards * cards = new DataCards(Era,
				    embedded,
				    nbins,
				    xmin,
				    xmax,
				    xDNN,
				    useTH2forZtt,
				    mvaDM,
				    applyIPcut,
				    runSystematics);

  cards->SetInputDirectory(Input_dir);
  cards->SetOutputDirectory(Output_dir);
  cards->SetOutputFileName(Output_filename);
  bool error = cards->Run(classIndex,channel);
  delete cards;

}
