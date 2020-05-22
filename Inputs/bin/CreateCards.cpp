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
  if(argc==4){
    classIndex = std::atoi(argv[2]);
    channel = argv[3];
    cout << "Running on class " << classIndex <<" - channel " << channel<<endl;
  }else if(argc==3){
    classIndex = std::atoi(argv[2]);
    cout << "Running on class " << classIndex <<endl;
  }

  const string era = cfg.get<string>("Era");
  TString Era(era);
  const bool embedded = cfg.get<bool>("IsEmbedded");
  const bool fakeFactor = cfg.get<bool>("FFmethod");

  const int nbinsMuPi = cfg.get<int>("NbinsPhiCPmupi");
  const int nbinsMuRho = cfg.get<int>("NbinsPhiCPmurho");
  const int nbinsMuA1 = cfg.get<int>("NbinsPhiCPmua1");
  const int nbinsMu0A1 = cfg.get<int>("NbinsPhiCPmu0a1");
  const vector<double> DNNbins = cfg.get< vector<double> >("DNNbins");
  const vector<double> DNNbinsZtt = cfg.get< vector<double> >("DNNbinsZtt");
  const vector<double> DNNbinsFakes = cfg.get< vector<double> >("DNNbinsFakes");
  const bool splitBkg= cfg.get<bool>("SplitBkg");

  const bool useTH1forHiggs = cfg.get<bool>("UseTH1ForHiggs");
  const bool useTH2forZtt = cfg.get<bool>("UseTH2ForZtt");
  const bool useTH2forFakes = cfg.get<bool>("UseTH2ForFakes");
  const bool mvaDM = cfg.get<bool>("mvaDM");
  const bool applyIPcut = cfg.get<bool>("ApplyIPcut");
  const bool applyIPcutOnBkg = cfg.get<bool>("ApplyIPcutOnBkg");
  const bool runSystematics = cfg.get<bool>("RunSystematics");

  const string variableCP = cfg.get<string>("CPvariables");
  TString VariableCP(variableCP);

  const string input_dir = cfg.get<string>("InputDirectory");
  TString Input_dir(input_dir);

  const string output_dir = cfg.get<string>("OutputDirectory");
  TString Output_dir(output_dir);

  const string output_filename = cfg.get<string>("OutputFileName");
  TString Output_filename(output_filename);

  const string cutIpMuon = cfg.get<string>("CutIPmuon");
  TString CutIpMuon(cutIpMuon);

  const string cutIpPion = cfg.get<string>("CutIPpion");
  TString CutIpPion(cutIpPion);

  /*
  for (auto XDNN : DNNbins) 
    cout << " " << XDNN;
  cout << endl;
  exit(1);
  */

  map<TString,int> binsperchannel;
  binsperchannel["mupi"]=nbinsMuPi;
  binsperchannel["murho"]=nbinsMuRho;
  binsperchannel["mua1"]=nbinsMuA1;
  binsperchannel["mu0a1"]=nbinsMu0A1;
  double xmin = 0.0;
  double xmax = 2*TMath::Pi();
  vector<double> xDNNSig = DNNbins;
  vector<double> xDNNZtt = DNNbinsZtt;
  vector<double> xDNNFakes = DNNbinsFakes;
  bool SplitBkg = splitBkg;
  bool UseTH1forHiggs = useTH1forHiggs;
  bool UseTH2forZtt = useTH2forZtt;
  bool UseTH2forFakes = useTH2forFakes;
  if(argc==3&&classIndex!=0){
    SplitBkg=false;
    UseTH2forZtt=false;
    UseTH2forFakes=false;
  }else if(argc==4&&classIndex==0) UseTH1forHiggs=false;
  
  if(!SplitBkg&&argc==4){
    cout << "When bkg are not splitted by decay channel, the code cannot run for separate decay channels" << endl << "Please run the code for Signal category as" <<endl << "CreateCards config classIndex decay-channel" << endl << "and Bkg categories as" << endl << "CreateCards config classIndex" <<endl;
    exit(EXIT_FAILURE);
  }

  DataCards * cards = new DataCards(Era,
				    embedded,
				    fakeFactor,
				    VariableCP,
				    binsperchannel,
				    xmin,
				    xmax,
				    xDNNSig,
				    xDNNZtt,
				    xDNNFakes,
				    SplitBkg,
				    UseTH1forHiggs,
				    UseTH2forZtt,
				    UseTH2forFakes,
				    mvaDM,
				    applyIPcut,
				    applyIPcutOnBkg,
				    runSystematics);

  cards->SetCutIPmuon(CutIpMuon);
  cards->SetCutIPpion(CutIpPion);

  cards->SetInputDirectory(Input_dir);
  cards->SetOutputDirectory(Output_dir);
  cards->SetOutputFileName(Output_filename);
  bool error = cards->Run(classIndex,channel);
  delete cards;

  cout << "DONE" << endl;

}
