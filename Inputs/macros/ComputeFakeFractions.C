#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "DesyTauAnalyses/NTupleMaker/test/Plotting.h"
#include "DesyTauAnalyses/NTupleMaker/test/Plotting_Style.h"
#include "DesyTauAnalyses/NTupleMaker/test/HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "HiggsCP/Inputs/settings_for_eras_newNTuples.h"

void ComputeFakeFractions(TString Variable = "m_vis",
			  TString VariableName = "m_vis",
			  TString xtitle = "m_{vis} [GeV]",
			  const int nBins  =   15,
			  float xmin = 0,
			  float xmax = 250,  
			  TString Weight = "puweight*mcweight*effweight*",
			  TString Cuts = "(((trg_singlemuon>0.5&&pt_1>25)||(trg_mutaucross&&pt_1>21&&pt_2>32&&abs(eta_1)<2.1&&abs(eta_2)<2.1))&&iso_1<0.15&&pt_1>21&&pt_2>20&&byTightDeepTau2017v2p1VSmu_2>0.5&&byVVVLooseDeepTau2017v2p1VSe_2>0.5&&extraelec_veto<0.5&&extramuon_veto<0.5&&dilepton_veto<0.5&&byMediumDeepTau2017v2p1VSjet_2<0.5&&puppimt_1<50&&nbtag==0&&abs(eta_1)<2.4&&abs(eta_2)<2.3&&os>0.5)",
			  TString directory = "/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/Jan20/CMSSW_10_2_16/src/DesyTauAnalyses/NTupleMaker/test/mutau/2017/",
			  TString outputDir = "./",
			  TString suffix = "_antiIso",
			  TString era = "2017") {

  TH1::SetDefaultSumw2();
  SetStyle();
  const map<TString, double>  *xsec_map    = 0;
  if (era == "2018") xsec_map = &xsec_map_2018;
  else if (era == "2017") xsec_map = &xsec_map_2017;
  else if (era == "2016") xsec_map = &xsec_map_2016;
  else exit(EXIT_FAILURE);

  float lumi;
  if(era=="2018") lumi = 59970;
  else if(era=="2017") lumi = 41860;
  else if(era=="2016") lumi = 36773;
  else exit(EXIT_FAILURE);

  TString ytitle = "process / data";
  
   
  // These samples names should match the root file names (-> sampleNames[i].root is read in later)

  vector<TString> sampleNames;
 if(era!="2016") sampleNames= {
        "SingleMuon_Run", // data (0)
        "DYJetsToLL_M-50",// (1)Drell-Yan
        "WJetsToLNu",// (2)WJets
        "WW",// (3)WW
        "WZ",// (4)WZ
        "ZZ",// (5)ZZ
        "ST_tW_antitop_5f", // (6) SingleTop tW tbar
        "ST_tW_top_5f", // (7) SingleTop tW t
        "ST_t-channel_antitop_4f",// (8) SingleTop t antitop
        "ST_t-channel_top_4f",// (9) SingleTop t top
        "TTTo2L2Nu",//(10)TTbar leptonic
        "TTToHadronic",//(11) hadronic
        "TTToSemiLeptonic"//(12) semileptonic
  };
  else  sampleNames = {
        "SingleMuon_Run", // data (0)
        "DYJetsToLL_M-50",// (1)Drell-Yan Z->TT
        "WJetsToLNu",// (2)WJets
        "WW",// (3)WW
        "WZ",// (4)WZ
        "ZZ",// (5)ZZ
        "ST_tW_antitop_5f", // (6) SingleTop tW tbar
        "ST_tW_top_5f", // (7) SingleTop tW t
        "ST_t-channel_antitop_4f",// (8) SingleTop t antitop
        "ST_t-channel_top_4f",// (9) SingleTop t top
        "TT"//(10) semileptonic
   };
  sampleNames[0]+=era;
  const int nSamples = sampleNames.size(); //DY is used twice, for Zll and Ztt


  // *******************************
  // ***** Selection Cuts    *******
  // *******************************

  TString cuts[nSamples];
  TString cutstrueT[nSamples];
  TString fakeT="*(gen_match_2==6)";
  TString trueT="*(gen_match_2!=6)";

  // MC specific cuts to select certain type of particle
  
  // Selection cuts applied to all samples
  for (int i=0; i<nSamples; ++i) {
    cuts[i]   = Weight+Cuts+fakeT;
    cutstrueT[i] = Weight+Cuts+trueT;
  };

  cuts[0] = Cuts; //DATA
   
  cutstrueT[0] = Cuts;

  // *******************************
  // ***** Filling Histograms ******
  // *******************************

  TH1D * hist[nSamples];
  TH1D * histtrueT[nSamples];


  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);

  // Draw main selection for all histograms in sampleNames
  for (int i=0; i<nSamples; ++i) {

    cout << endl << sampleNames[i] << ":" << endl;
    // Reading input file
    TFile * file = new TFile( directory + sampleNames[i] + ".root");

    // Get tree and one further important histogram from input file
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents"); //inputEventsH
    TTree * tree = (TTree*)file->Get("TauCheck");

    // Calculate normalization of this sample
    float xsec = 1;
    if(!sampleNames[i].Contains("Run")) xsec = xsec_map->at(sampleNames[i]);
    else cout << cuts[i] <<endl;
    double norm = xsec*lumi/nWeightedEvents->GetSumOfWeights(); 

    cout << "xsec: " << xsec << endl;
    cout << "lumi: " << lumi << endl;
    cout << "norm: " << norm << endl;

    // Name and initialize histograms
    TString histName   = sampleNames[i] + VariableName + "_fake";
    TString histNametrueT = sampleNames[i] + VariableName + "_trueT";
    hist[i]   = new TH1D(histName,"",nBins,xmin,xmax);
    histtrueT[i] = new TH1D(histNametrueT,"",nBins,xmin,xmax);

    cout << "Drawing ..." << endl;
    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histNametrueT,cutstrueT[i]);

    if (i>0) // if sample is MC sample -> Scale to xsec and luminosity
      {
	hist[i]   -> Scale(norm);
	histtrueT[i] -> Scale(norm);
      }

    cout << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }
  cout << endl;
  delete dummyCanv;


  // *****************************************
  // ***** Stitching of W+Jets samples *******
  // *****************************************

  TH1D * histW[10];
  TH1D * histWtrueT[10];
  int nSamplesW = 9;
    
  TString npartonCuts[9] = {"*(gen_noutgoing==0||gen_noutgoing>4)", //cut on inclusive sample
			    "*(gen_noutgoing==1)",//cut on inclusive sample
			    "*(gen_noutgoing==2)",//cut on inclusive sample
			    "*(gen_noutgoing==3)",//cut on inclusive sample
			    "*(gen_noutgoing==4)",//cut on inclusive sample
			    "",
			    "",
			    "",
			    ""
  };

  TString refSamples[6];
  double refEvents[6] = {0,0,0,0,0,0};
  double refXSec[6] = {0,0,0,0,0,0};

  // redefine reference cross sections and reference samples
  refSamples[0] = "WJetsToLNu";
  refSamples[1] = "W1JetsToLNu";
  refSamples[2] = "W2JetsToLNu";
  refSamples[3] = "W3JetsToLNu";
  refSamples[4] = "W4JetsToLNu";
  refEvents[0] = 0;
  refEvents[1] = 0;
  refEvents[2] = 0;
  refEvents[3] = 0;
  refEvents[4] = 0;

  for (int iW=0; iW<5; ++iW) {
    refXSec[iW] = xsec_map->at(refSamples[iW]);
    TFile * file = new TFile(directory+refSamples[iW]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iW] = nWeightedEvents->GetSumOfWeights();   //number of events with amc@NLO weight
  }

  TString wSampleNames[9] = {"WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "WJetsToLNu",
			     "W1JetsToLNu",
			     "W2JetsToLNu",
			     "W3JetsToLNu",
			     "W4JetsToLNu"
  };

  double wNorm[9];
  wNorm[0] = lumi*refXSec[0]/refEvents[0];  //norm for inlcusive with all events
  wNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]); //add inclusive sample (all events with one parton) with W1Jet sample
  wNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  wNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  wNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  wNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  wNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);

  TString cutsW[9];
  TString cutsWtrueT[9];

  for (int iW=0; iW<9; ++iW) {
    cutsW[iW]   = cuts[2]+npartonCuts[iW]; //apply parton cuts
    cutsWtrueT[iW] = cutstrueT[2]+npartonCuts[iW];
  }

  // filling histograms for WJets samples
  for (int i=0; i<nSamplesW; ++i) { // run over W+Jets samples
    cout<< wSampleNames[i] <<endl;
    TFile * file = new TFile(directory+wSampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = wNorm[i];

    TString histNameW   = wSampleNames[i] + VariableName + "_w_fake";
    TString histNameWtrueT = wSampleNames[i] + VariableName + "_w_trueT";
    histW[i]   = new TH1D(histNameW,"",nBins,xmin,xmax);
    histWtrueT[i] = new TH1D(histNameWtrueT,"",nBins,xmin,xmax);

    tree->Draw(Variable+">>"+histNameW,  cutsW[i]); //fill histogram with cuts applied
    tree->Draw(Variable+">>"+histNameWtrueT,cutsWtrueT[i]);

    histW[i]   -> Scale(norm);
    histWtrueT[i] -> Scale(norm);
  }

  hist[2]   = histW[0];
  histtrueT[2] = histWtrueT[0];

  for (int iW=1; iW<9; ++iW)
    {
      hist[2]   -> Add(hist[2],histW[iW]);
      histtrueT[2] -> Add(histtrueT[2],histWtrueT[iW]);
    }


  // ********************************************
  // ***** Stitching of Drell-Yan samples *******
  // ********************************************
  int nSamplesDY = 9;

  TH1D * histZtt[9];
  TH1D * histZtttrueT[9];

  refSamples[0] = "DYJetsToLL_M-50";
  refSamples[1] = "DY1JetsToLL_M-50";
  refSamples[2] = "DY2JetsToLL_M-50";
  refSamples[3] = "DY3JetsToLL_M-50";
  refSamples[4] = "DY4JetsToLL_M-50";

  for (int iDY=0; iDY<5; ++iDY) {
    refXSec[iDY] = xsec_map->at(refSamples[iDY]);
    TFile * file = new TFile(directory+refSamples[iDY]+".root");
    TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
    refEvents[iDY] = nWeightedEvents->GetSumOfWeights();
  }

  

  TString dySampleNames[9] = {"DYJetsToLL_M-50",
			      "DYJetsToLL_M-50",
			      "DYJetsToLL_M-50",
			      "DYJetsToLL_M-50",
			      "DYJetsToLL_M-50",
			      "DY1JetsToLL_M-50",
			      "DY2JetsToLL_M-50",
			       "DY3JetsToLL_M-50",
			      "DY4JetsToLL_M-50"
  };

  double dyNorm[9];
  dyNorm[0] = lumi*refXSec[0]/refEvents[0];
  dyNorm[1] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[2] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[3] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[4] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  dyNorm[5] = lumi/(refEvents[0]/refXSec[0]+refEvents[1]/refXSec[1]);
  dyNorm[6] = lumi/(refEvents[0]/refXSec[0]+refEvents[2]/refXSec[2]);
  dyNorm[7] = lumi/(refEvents[0]/refXSec[0]+refEvents[3]/refXSec[3]);
  dyNorm[8] = lumi/(refEvents[0]/refXSec[0]+refEvents[4]/refXSec[4]);
  //dyNorm[9] = lumi*refXSec[5]/refEvents[5];

  TString cutsZtt[9];
  TString cutsZtttrueT[9];

  for (int iDY=0; iDY<nSamplesDY; ++iDY) {
    cutsZtt[iDY]   = cuts[1]+npartonCuts[iDY];
    cutsZtttrueT[iDY] = cutstrueT[1]+npartonCuts[iDY];
  }


  // filling histograms for DY samples
  for (int i=0; i<nSamplesDY; ++i) { // run over samples
    cout<<dySampleNames[i]<<endl;
    TFile * file = new TFile(directory+dySampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    double norm = dyNorm[i];

    TString histNameZtt   = dySampleNames[i] + VariableName + "_ztt_fake";
    TString histNameZtttrueT = dySampleNames[i] + VariableName + "_ztt_trueT";
    histZtt[i]   = new TH1D(histNameZtt,"",nBins,xmin,xmax);
    histZtttrueT[i] = new TH1D(histNameZtttrueT,"",nBins,xmin,xmax);

    tree -> Draw(Variable+">>"+histNameZtt,  cutsZtt[i]);
    tree -> Draw(Variable+">>"+histNameZtttrueT,cutsZtttrueT[i]);

    histZtt[i]   -> Scale(norm);
    histZtttrueT[i] -> Scale(norm);

    // cout << dySampleNames[i] << " -> ZTT : Entries = " << histZtt[i]->GetEntries()
    // 	 << " : Sum of weights = " << histZtt[i]->GetSumOfWeights()
    // 	 << "    ZLL : Entries = " << histZll[i]->GetEntries() << " : Sum of weights = " << histZll[i]->GetSumOfWeights() << endl;
  }

  hist[1]   = histZtt[0];
  histtrueT[1] = histZtttrueT[0];

  for (int iDY=1; iDY<9; ++iDY) {
    hist[1]  -> Add(hist[1],histZtt[iDY]);
    histtrueT[1]-> Add(histtrueT[1],histZtttrueT[iDY]);
  }

  // ********************************************
  // ***** Adding similar backgrounds     *******
  // ********************************************


  //For FF we need the following groups:
  //V(j->tau) = DY, VV, W
  //tt(j->tau) = top
  //QCD = ss->os

  //Adding DY, W and VV
  for (int iH=2; iH<10; ++iH) hist[1]->Add(hist[1],hist[iH]);

  //Adding top
  if(era!="2016") for (int iH=11; iH<nSamples; ++iH) hist[10]->Add(hist[10],hist[iH]);

  // Adding trueT backgrounds
  for (int iH=2; iH<nSamples; ++iH) histtrueT[1]->Add(histtrueT[1],histtrueT[iH]);


  // ********************************************
  // ***** QCD background estimation      *******
  // ********************************************

  float datatrueT     = histtrueT[0]->GetSumOfWeights();
  float datatrueTfull = histtrueT[0]->Integral(0,nBins+1);
  
  // Subtracting background from trueT
  
  TH1D * nonQCD = (TH1D*)histtrueT[1]  -> Clone("nonQCD");
  nonQCD->Add(histtrueT[1],hist[1]);
  nonQCD->Add(nonQCD,hist[10]);
  histtrueT[0] ->Add(histtrueT[0],nonQCD,1,-1);

  //histtrueT[0]->Add(hist[0],histtrueT[1],1,-1);
  //for (int iH=1; iH<nSamples; ++iH) histtrueT[0]->Add(histtrueT[0],hist[iH],1,-1);
  // Remove negative values for QCD bkg
  for (int iB=1; iB<=nBins; ++iB) {
    double value= histtrueT[0] -> GetBinContent(iB);
    if(value<1e-5) histtrueT[0] -> SetBinContent(iB,0.);
  }

  // ************************************
  // ***** Summarize backgrounds  *******
  // ************************************
  
  TH1D * histData = (TH1D*)hist[0]   -> Clone("data_obs");
  TH1D * QCD      = (TH1D*)histtrueT[0] -> Clone("QCD");
  TH1D * Vfakes      = (TH1D*)hist[1]   -> Clone("Vfakes");
  TH1D * Topfakes      = (TH1D*)hist[10]   -> Clone("Topfakes");
  TH1D * TrueTaus      = (TH1D*)histtrueT[1]  -> Clone("TrueTaus");

 //TH1D * SMH      = (TH1D*)hist[14]  -> Clone("SMH");
  //TH1D * BSMH     = (TH1D*)hist[15]  -> Clone("BSMH");
  for(int i=0;i<nSamples;i++){
    cout << setw(nSamples) << sampleNames[i] << " : Entries = " << hist[i]->GetEntries() << " : Integral = " << hist[i]->Integral(0,nBins+1) << endl;
  }

  cout << endl;
  cout << "QCD : Sum of weights = " << QCD->GetSumOfWeights() << " : Integral = " << QCD->Integral(1,nBins+1) << endl;
  cout << "Vfakes  : Sum of weights = " << Vfakes->GetSumOfWeights()  << " : Integral = " << Vfakes->Integral(1,nBins+1)  << endl;
  cout << "Topfakes  : Sum of weights = " << Topfakes->GetSumOfWeights()  << " : Integral = " << Topfakes->Integral(1,nBins+1)  << endl;
  cout << "TrueTaus : Sum of weights = " << TrueTaus->GetSumOfWeights() << " : Integral = " << TrueTaus->Integral(1,nBins+1) << endl;
  cout << "Data : Sum of weights = " << histData->GetSumOfWeights() << " : Integral = " << histData->Integral(1,nBins+1) << endl;
  
  QCD->Divide(QCD,histData);
  Vfakes->Divide(Vfakes,histData);
  Topfakes->Divide(Topfakes,histData);
  TrueTaus->Divide(TrueTaus,histData);
  histData->Divide(histData,histData);
  

  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************

  TH1D * dummy = (TH1D*)Vfakes->Clone("dummy");
  float errQCD = 0.15; // ad-hoc sys uncertainty of QCD background
  float errW   = 0.15; // ad-hoc sys uncertainty of W+Jets background
  float errTT  = 0.08; // ad-hoc sys uncertainty of TT background

  for (int iB=1; iB<=nBins; ++iB) {   // Add general systematic uncertainties to each bin as error
    float eQCD   = errQCD*QCD->GetBinContent(iB);
    float eW     = errW*Vfakes->GetBinContent(iB);
    float eTT    = errTT*Topfakes->GetBinContent(iB);
    float err2   = eQCD*eQCD+ eW*eW + eTT*eTT;     //TO DO: WAS IST MIT DEM FEHLER AUF DY?
    float errTot = TMath::Sqrt(err2);
    cout << "eQCD: " << eQCD << "  eW: " << eW << "  eTT: " << eTT << "  eTotal: " << errTot << endl;
    dummy -> SetBinError(iB,errTot);
    //SMH   -> SetBinError(iB,0);
    //BSMH  -> SetBinError(iB,0);
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING 
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();
  bool LargeScale = false;  
  TCanvas* canv1;
  if(LargeScale)canv1 = new TCanvas("c1", "c1", 2000,800);
  else canv1 = new TCanvas("c1", "c1", 1200,1000);
  canv1->cd();
  vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);

  vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
  h[0] -> Draw();

  string units="";
  string xtitle_ = (string) xtitle;
  size_t pos = xtitle_.find("[");
  if(pos!=string::npos) {
    units = xtitle_.substr(pos+1, xtitle_.find("]") - pos -1 );
    xtitle_ = xtitle_.substr(0, pos);
  }

  pads[1] -> cd();
  h[1]    -> Draw();
  SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.65, 1.35);
  StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
  h[1] -> GetYaxis()->SetNdivisions(4);
  h[1] -> GetXaxis()->SetTitleOffset(0.95);
  h[1] -> GetXaxis()->SetNdivisions(505);
  h[1] -> GetYaxis()->SetTitleOffset(1.1);
  if(LargeScale)  h[1] -> GetYaxis()->SetTitleOffset(0.9);
  pads[0] -> cd();
  h[0] -> GetYaxis()->SetTitleOffset(1.6);
  if(LargeScale)h[0] -> GetYaxis()->SetTitleOffset(1);
  pads[1] -> SetGrid(0,1);
  pads[0] -> cd();
  gPad->SetTicky();
  // Setup legend
  TLegend *legend = PositionedLegend(0.25, 0.1, 3, 0.03);
  legend -> SetTextFont(42);
  histData -> SetMarkerColor(1);
  histData -> SetLineColor(1);
  histData -> SetFillColor(1);
  histData -> SetFillStyle(0);
  histData -> SetLineWidth(2);
  histData -> SetMarkerStyle(20);
  histData -> SetMarkerSize(1.1);

  InitHist(QCD,TColor::GetColor("#FFCCFF"));
  InitHist(Topfakes,TColor::GetColor("#9999CC"));
  InitHist(Vfakes,TColor::GetColor("#DE5A6A"));
  InitHist(TrueTaus,TColor::GetColor("#FFCC66"));

  legend -> AddEntry(histData, "Observed", "ple");
  legend -> AddEntry(Topfakes,"top+(j->#tau)","f");
  legend -> AddEntry(Vfakes,"V+(j->#tau)","f");
  legend -> AddEntry(QCD,"QCD","f");
  legend -> AddEntry(TrueTaus,"True Taus","f");
  legend -> SetNColumns(2);

  // Add all bkg contributions to one stack plot
  THStack *stack = new THStack("Background","");
  stack -> Add(Topfakes);
  stack -> Add(Vfakes);
  stack -> Add(QCD);
  stack -> Add(TrueTaus);
  stack -> Draw("hsame");

  canv1->Update();
  /*
  InitSignal(SMH ,2);
  InitSignal(BSMH,4);
  if (showSignal)
    {
      legend->AddEntry(SMH,"SM Higgs(125) #times 200","f");
      SMH->Draw("hsame");
      legend->AddEntry(BSMH,"BSM Higgs(125) #times 200","f");
      BSMH->Draw("hsame");
    }
  */
  canv1->Update();

  // Initialize a histogram which adds all error up
  TH1D * bkgdErr = (TH1D*)stack->GetStack()->Last()->Clone("bkgdErr");
  float errLumi = 0.03;
  float errMuon = 0.03;
  float errElectron = 0.04;
  for (int iB=1; iB<=nBins; ++iB) {
    QCD->SetBinError(iB,0);
    Vfakes->SetBinError(iB,0);
    Topfakes->SetBinError(iB,0);
    float eStat =  bkgdErr->GetBinError(iB);
    float X = bkgdErr->GetBinContent(iB);
    float eLumi = errLumi * X;
    float eMuon = errMuon * X;
    float eElectron = errElectron * X;
    float eBkg = dummy->GetBinError(iB);
    float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg+eMuon*eMuon+eElectron*eElectron);
    bkgdErr->SetBinError(iB,Err);
    cout << "eStat = " << eStat << " : eLumi = "<< eLumi <<" : eBkg = " << eBkg << endl;
  }

  bkgdErr -> SetMarkerSize(0);
  int new_idx = 923;//CreateTransparentColor(13,1.0);
  bkgdErr -> SetFillColor(new_idx);
  bkgdErr -> SetFillStyle(3004);
  bkgdErr -> SetLineWidth(1);
  bkgdErr -> Draw("e2same");
  legend  -> AddEntry(bkgdErr, "Bkg. uncertainty" , "F" );
  canv1   -> Update();

  TH1D * ratioH    = (TH1D*)histData -> Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr  -> Clone("ratioErrH");
  ratioH -> Divide((TH1D*)stack->GetStack()->Last()); // Divide by the sum of the THStack

  // Set error of MC bkg correctly in ratio
  for (int iB=1; iB<=nBins; ++iB) {
    ratioErrH -> SetBinContent(iB,1.0);
    ratioErrH -> SetBinError(iB,0.0);
    float xBkg   = bkgdErr -> GetBinContent(iB);
    float errBkg = bkgdErr -> GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);
    }
  }

  pads[1]->cd();
  ratioErrH->Draw("e2same");
  ratioH->Draw("pe0same");

  pads[0]->cd();
  histData->Draw("pesame");

  FixTopRange(pads[0], GetPadYMax(pads[0]), 0.115);
  //DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
  if (lumi== 28686) DrawTitle(pads[0], "28.7 fb^{-1} (13 TeV, 2017)", 3);
  else if (lumi== 13960) DrawTitle(pads[0], "14.0 fb^{-1} (13 TeV, 2017)", 3);
  else if (lumi== 41465) DrawTitle(pads[0], "41.5 fb^{-1} (13 TeV, 2017)", 3);
  else DrawTitle(pads[0], "42.8 fb^{-1} (13 TeV, 2017)", 3);
  DrawTitle(pads[0], "#scale[1.2]{#bf{CMS} Work in progress}", 1);
  FixBoxPadding(pads[0], legend, 0.05);
  legend->Draw();
  FixOverlay();
  canv1->Update();
  pads[0]->GetFrame()->Draw();

  canv1 -> Print( outputDir + "Fakes_" + VariableName +suffix + ".png" );

  TFile * file = new TFile("./FakeFractions_"+VariableName+".root","recreate");

  TH1D * ff_QCD      = new TH1D("ff_QCD","",nBins,xmin,xmax);
  TH1D * ff_W      = new TH1D("ff_W","",nBins,xmin,xmax);
  TH1D * ff_tt      = new TH1D("ff_tt","",nBins,xmin,xmax);
  ff_QCD->Divide(QCD,histData);
  ff_W->Divide(Vfakes,histData);
  ff_tt->Divide(Topfakes,histData);
  // Renormalize to 1 the sum of fractions
  
    for (int iB=1; iB<=nBins; ++iB) {
    double value_QCD=ff_QCD->GetBinContent(iB);
    double value_tt =ff_tt ->GetBinContent(iB);
    double value_W  =ff_W  ->GetBinContent(iB);
    ff_QCD->SetBinContent(iB,value_QCD/(value_QCD+value_tt+value_W));
    ff_tt ->SetBinContent(iB, value_tt/(value_QCD+value_tt+value_W));
    ff_W  ->SetBinContent(iB,  value_W/(value_QCD+value_tt+value_W));
  }
  ff_QCD->Write("ff_QCD");
  ff_W->Write("ff_W");
  ff_tt->Write("ff_tt");
}
