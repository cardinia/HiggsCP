#include <iostream>
#include <fstream>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLatex.h"
#include "HiggsCP/Inputs/interface/HttStylesNew.cc"
#include "HiggsCP/Inputs/interface/CMS_lumi.C"

// embedded = true (EMB), false (DY MC)
// era = 2016, 2017, 2018
// acotautau_refitbs_00, acotautau_refitbs_01
// alphaminus, alpha_IP_1, alpha_IP_1, alpha_plane_2, alpha_plane_2
// trigger = 0 (single-e OR e-tau), 1 (single-e only), 2 (e-tau only) 

struct HistAttr {
  TString variable;
  TString xtitle;
  TString ytitle;
  int nBins;
  float xmin;
  float xmax;
  float yLower;
  float scaleYUpper;
  bool logX;
  bool logY;
  bool plotLegend;
};

using namespace std;

typedef std::map<TString,TTree*> MapOfTrees;

void PlotDNN(TString channel,
	     MapOfTrees mapOfTrees,
	     HistAttr histAttr,
	     bool embedded,
	     TString era,
	     int trigger,
	     int pageType
	     ) {

  
  TString Variable = histAttr.variable;
  TString xtitle = histAttr.xtitle;
  TString ytitle = histAttr.ytitle;
  int nBins  = histAttr.nBins;
  float xmin = histAttr.xmin;
  float xmax = histAttr.xmax;
  float yLower = histAttr.yLower;
  float scaleYUpper = histAttr.scaleYUpper;
  bool logY = histAttr.logY;
  bool logX = histAttr.logX;
  bool plotLegend = histAttr.plotLegend;

  std::cout << std::endl;
  std::cout << "Variable : " << Variable << std::endl;
  std::cout << std::endl;

  double fake_scale = 1.0;
  double zll_scale = 1.0;
  double ztt_scale = 1.0;

  TString pTe_cut("26");
  if (era=="2017")
    pTe_cut = "28";
  if (era=="2018")
    pTe_cut = "33";
  
  TString triggerOption("");
  if (trigger==1)
    triggerOption = "_singleEleOnly";
  else if (trigger==2)
    triggerOption = "_xTrigOnly";
  else if (trigger==3)
    triggerOption = "_resolved";

  if (era=="2016") triggerOption = "";
  
  SetStyle();
  
  TString suffix = "_mc";
  if (embedded) suffix = "_embedded";
  

  lumi_13TeV = "2018, 59.7 fb^{-1}";
  if (era=="2017")
    lumi_13TeV = "2017, 41.0 fb^{-1}";
  if (era=="2016")
    lumi_13TeV = "2016, 35.9 fb^{-1}";

  bool nonEquidistant = false;
  float bins[3] = {0.05,1.5,10.0};
  if (Variable.Contains("IP_signif")&&Variable.Contains("_2")) {
    nBins = 2;
    nonEquidistant = true;
  }

  TString Weight("weight*");
  TString Weight1("weight*");
  TString Cuts("trg_singleelectron>0.5&&pt_1>26&&TMath::Abs(eta_1)<2.1&&pt_2>20&&TMath::Abs(eta_2)<2.3&&os>0.5&&puppimt_1<50");
  TString Cuts1("pt_1>25");

  if (era=="2017"||era=="2018") {
    if (trigger==1) {
      Cuts = "trg_singleelectron>0.5&&pt_1>"+pTe_cut+"&&TMath::Abs(eta_1)<2.1&&pt_2>20&&TMath::Abs(eta_2)<2.3&&os>0.5&&puppimt_1<50";
      Weight = "(weight*trigweight_1/trigweight)*";
    }
    else if (trigger==2) {
      Cuts = "trg_etaucross>0.5&&pt_1>25&&TMath::Abs(eta_1)<2.1&&pt_2>35&&TMath::Abs(eta_2)<2.1&&os>0.5&&puppimt_1<50";
      Weight = "(weight*trigweight_2/trigweight)*";
    }
    else if (trigger==3) {
      Cuts = "trg_singleelectron>0.5&&pt_1>"+pTe_cut+"&&TMath::Abs(eta_1)<2.1&&pt_2>20&&TMath::Abs(eta_2)<2.3&&os>0.5&&puppimt_1<50";
      Cuts1 = "trg_etaucross>0.5&&pt_1>25&&pt_1<="+pTe_cut+"&&TMath::Abs(eta_1)<2.1&&pt_2>35&&TMath::Abs(eta_2)<2.1&&os>0.5&&puppimt_1<50";
      Weight = "(weight*trigweight_1/trigweight)*";
      Weight1 = "(weight*trigweight_2/trigweight)*";
    }
    else
      Cuts = "pt_1>25&&TMath::Abs(eta_1)<2.1&&pt_2>20&&TMath::Abs(eta_2)<2.3&&os>0.5&&puppimt_1<50&&((trg_singleelectron>0.5&&pt_1>"+pTe_cut+")||(trg_etaucross>0.5&&pt_1>25&&pt_2>35&&TMath::Abs(eta_2)<2.1))";
  }
  Cuts += "&&nbtag==0";
  Cuts1 += "&&nbtag==0";

  if (Variable=="jpt_1"||Variable=="jeta_1") {
    Cuts += "&&njets>0";
    Cuts1 += "&&njets>0";
  }
  if (Variable=="jpt_2"||Variable=="jeta_2"||Variable=="mjj"||Variable=="jdeta") {
    Cuts += "&&njets>1";
    Cuts1 += "&&njets>1";
  }
  if (Variable.Contains("IP_signif")&&Variable.Contains("_2")) {
    Cuts  += "&&dmMVA_2==0";
    Cuts1 += "&&dmMVA_2==0";
  }

  TString Cuts_Sig    = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2>0.5");
  TString Cuts_FF     = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5");

  TString Cuts1_Sig    = Cuts1 + TString("&&byMediumDeepTau2017v2p1VSjet_2>0.5");
  TString Cuts1_FF     = Cuts1 + TString("&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5");

  TString CutsZTT_Sig = Cuts_Sig + TString("&&gen_match_1==3&&gen_match_2==5");
  TString CutsZTT_FF  = Cuts_FF  + TString("&&gen_match_1==3&&gen_match_2==5");
  TString Cuts1ZTT_Sig = Cuts1_Sig + TString("&&gen_match_1==3&&gen_match_2==5");
  TString Cuts1ZTT_FF  = Cuts1_FF  + TString("&&gen_match_1==3&&gen_match_2==5");

  TString CutsZLL_Sig  = Cuts_Sig + TString("&&!(gen_match_1==3&&gen_match_2==5)");
  TString CutsZLL_FF   = Cuts_FF  + TString("&&!(gen_match_1==3&&gen_match_2==5)");
  TString Cuts1ZLL_Sig = Cuts1_Sig + TString("&&!(gen_match_1==3&&gen_match_2==5)");
  TString Cuts1ZLL_FF  = Cuts1_FF  + TString("&&!(gen_match_1==3&&gen_match_2==5)");

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  std::vector<TString> sampleNames = {
    "Data",    // data  (0)
    "TTbar",   // TT    (1)
    "Diboson", // VV    (2)
    "WJets",   // W     (3)
    "ZLL",  // Z->ll (4)
    "ZTT"   // Z->tt or EMB (5)
  };

  TString cuts_Sig[10];
  TString cuts_FF[10];
  TString cuts1_Sig[10];
  TString cuts1_FF[10];
  for (int i=0; i<10; ++i) {
    if (embedded) {
      cuts_Sig[i] = Weight+"("+CutsZLL_Sig+"&&gen_match_2!=6)";
      cuts_FF[i]  = Weight+"ff_mva*("+CutsZLL_FF+"&&gen_match_2!=6)";
      cuts1_Sig[i] = Weight1+"("+Cuts1ZLL_Sig+"&&gen_match_2!=6)";
      cuts1_FF[i]  = Weight1+"ff_mva*("+Cuts1ZLL_FF+"&&gen_match_2!=6)";
    }
    else {
      cuts_Sig[i] = Weight+"("+Cuts_Sig+"&&gen_match_2!=6)";
      cuts_FF[i]  = Weight+"ff_mva*("+Cuts_FF+"&&gen_match_2!=6)";
      cuts1_Sig[i] = Weight1+"("+Cuts1_Sig+"&&gen_match_2!=6)";
      cuts1_FF[i]  = Weight1+"ff_mva*("+Cuts1_FF+"&&gen_match_2!=6)";
    }
  }
  cuts_Sig[0] = Cuts_Sig;
  cuts_FF[0] = "ff_mva*("+Cuts_FF+")";
  cuts1_Sig[0] = Cuts1_Sig;
  cuts1_FF[0] = "ff_mva*("+Cuts1_FF+")";

  if (embedded) {
    cuts_Sig[5] = Weight+"("+Cuts_Sig+"&&gen_match_2!=6)";
    cuts_FF[5]  = Weight+"ff_mva*("+Cuts_FF+"&&gen_match_2!=6)";
    cuts1_Sig[5] = Weight1+"("+Cuts1_Sig+"&&gen_match_2!=6)";
    cuts1_FF[5]  = Weight1+"ff_mva*("+Cuts1_FF+"&&gen_match_2!=6)";
  }
  else {
    cuts_Sig[4] = Weight+"("+CutsZLL_Sig+"&&gen_match_2!=6)";
    cuts_FF[4]  = Weight+"ff_mva*("+CutsZLL_FF+"&&gen_match_2!=6)";
    cuts_Sig[5] = Weight+"("+CutsZTT_Sig+"&&gen_match_2!=6)";
    cuts_FF[5]  = Weight+"ff_mva*("+CutsZTT_FF+"&&gen_match_2!=6)";
    cuts1_Sig[4] = Weight1+"("+Cuts1ZLL_Sig+"&&gen_match_2!=6)";
    cuts1_FF[4]  = Weight1+"ff_mva*("+Cuts1ZLL_FF+"&&gen_match_2!=6)";
    cuts1_Sig[5] = Weight1+"("+Cuts1ZTT_Sig+"&&gen_match_2!=6)";
    cuts1_FF[5]  = Weight1+"ff_mva*("+Cuts1ZTT_FF+"&&gen_match_2!=6)";
  }

  TH1D * hist_Sig[10];
  TH1D * hist_FF[10];

  TCanvas * dummy = new TCanvas("dummy","",400,400);

  // filling histograms
  for (unsigned int i=0; i<sampleNames.size(); ++i) {
    TTree * tree = mapOfTrees[sampleNames.at(i)];
    //    if (i==0) tree->MakeClass();
    TString histNameSig = sampleNames[i] + "_sig";
    TString histNameFF  = sampleNames[i] + "_ff";
    if (nonEquidistant) {
      hist_Sig[i] = new TH1D(histNameSig,"",nBins,bins);
      hist_FF[i] = new TH1D(histNameFF,"",nBins,bins);
    }
    else {
      hist_Sig[i] = new TH1D(histNameSig,"",nBins,xmin,xmax);
      hist_FF[i] = new TH1D(histNameFF,"",nBins,xmin,xmax);
    }
    hist_Sig[i]->Sumw2();
    hist_FF[i]->Sumw2();
    tree->Draw(Variable+">>"+histNameSig,cuts_Sig[i]);
    tree->Draw(Variable+">>"+histNameFF,cuts_FF[i]);
    if (trigger==3) {
      TH1D * hist1_Sig;
      TH1D * hist1_FF;
      if (nonEquidistant) {
	hist1_Sig = new TH1D("hist1_Sig","",nBins,bins);
	hist1_FF  = new TH1D("hist1_FF","",nBins,bins);
      }
      else {
	hist1_Sig = new TH1D("hist1_Sig","",nBins,xmin,xmax);
	hist1_FF  = new TH1D("hist1_FF","",nBins,xmin,xmax);
      }
      hist1_Sig->Sumw2();
      hist1_FF->Sumw2();
      tree->Draw(Variable+">>hist1_Sig",cuts1_Sig[i]);
      tree->Draw(Variable+">>hist1_FF",cuts1_FF[i]);
      hist_Sig[i]->Add(hist_Sig[i],hist1_Sig);
      hist_FF[i]->Add(hist_FF[i],hist1_FF);
      delete hist1_Sig;
      delete hist1_FF;
    }
  }
  delete dummy;
  
  // forming Fakes
  for (unsigned int i=1; i<sampleNames.size(); ++i) {
    hist_FF[0]->Add(hist_FF[0],hist_FF[i],1,-1);
  }

  TH1D * histData = (TH1D*)hist_Sig[0]->Clone("data_obs");
  TH1D * TT       = (TH1D*)hist_Sig[1]->Clone("TT");
  TH1D * EWK      = (TH1D*)hist_Sig[2]->Clone("EWK");
  TH1D * ZLL      = (TH1D*)hist_Sig[4]->Clone("ZLL");
  TH1D * ZTT      = (TH1D*)hist_Sig[5]->Clone("ZTT");
  TH1D * QCD      = (TH1D*)hist_FF[0]->Clone("fakes");

  std::cout << "Top : " << TT->GetSum() << std::endl;
  std::cout << "EWK : " << EWK->GetSum() << std::endl;
  std::cout << "QCD : " << QCD->GetSum() << std::endl;
  std::cout << "ZLL : " << ZLL->GetSum() << std::endl;
  std::cout << "ZTT : " << ZTT->GetSum() << std::endl;

  //  adding normalization systematics
  double ZTT_norm = 0.04; //  normalization ZTT :  4% (EMBEDDED)
  double EWK_norm = 0.07; //  normalization EWK :  7%
  double QCD_norm = 0.12; //  normalization Fakes : 12%
  double ZLL_mtau = 0.10; //  mu->tau fake rate ZLL : 10%
  double TT_norm  = 0.07; //  normalization TT  :  7%
  double lep_ID   = 0.03; //  lepton trigger ID :  3%
  double tau_ID   = 0.05; //  tau ID            :  5%

  for (int iB=1; iB<=nBins; ++iB) {

    float ztt  = ZTT->GetBinContent(iB);
    float ztte = ZTT->GetBinError(iB);
    ztte = TMath::Sqrt(ztte*ztte+ztt*ztt*(ZTT_norm*ZTT_norm+lep_ID*lep_ID+tau_ID*tau_ID));
    ZTT->SetBinError(iB,ztte);

    float ewk  = EWK->GetBinContent(iB);
    float ewke = EWK->GetBinError(iB);
    ewke = TMath::Sqrt(ewke*ewke+ewk*ewk*EWK_norm*EWK_norm);
    EWK->SetBinError(iB,ewke);

    float qcd  = QCD->GetBinContent(iB);
    float qcde = QCD->GetBinError(iB);
    qcde = TMath::Sqrt(qcde*qcde+qcd*qcd*QCD_norm*QCD_norm);
    QCD->SetBinError(iB,qcde);

    float tt  = TT->GetBinContent(iB);
    float tte = TT->GetBinError(iB);
    tte = TMath::Sqrt(tte*tte+tt*tt*TT_norm*TT_norm);
    TT->SetBinError(iB,tte);

    float zll  = ZLL->GetBinContent(iB);
    float zlle = ZLL->GetBinError(iB);
    zlle = TMath::Sqrt(zlle*zlle+zll*zll*(lep_ID*lep_ID+ZLL_mtau*ZLL_mtau));
    ZLL->SetBinError(iB,zlle);

  }

  ZTT->Scale(ztt_scale);
  QCD->Scale(fake_scale);
  ZLL->Scale(zll_scale);

  EWK->Add(EWK,TT);
  QCD->Add(QCD,EWK);
  ZLL->Add(ZLL,QCD);
  ZTT->Add(ZTT,ZLL);

  float yieldData   = histData->GetSum();
  float underData   = histData->GetBinContent(0);
  float overData    = histData->GetBinContent(nBins+1);
  float centerData  = histData->GetSumOfWeights();
  int   entriesData = histData->GetEntries();
  float checkData = underData + overData + centerData;

  float yieldMC  = ZTT->GetSum();
  float underMC  = ZTT->GetBinContent(0);
  float overMC   = ZTT->GetBinContent(nBins+1);
  float centerMC = ZTT->GetSumOfWeights();
  float checkMC  = underMC + overMC + centerMC;

  std::cout << std::endl;
  std::cout << "MC   : " << yieldMC << "   " << underMC << " (underflow)   " << overMC << " (overflow)" << std::endl;
  std::cout << "DATA : " << yieldData << "   " << underData << " (underflow)   " << overData << " (overflow)    " << entriesData << " (entries) "<< std::endl;
  std::cout << "---------" << std::endl;
  std::cout << "MC   : " << checkMC << std::endl;
  std::cout << "DATA : " << checkData << std::endl;
  std::cout << std::endl;

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  InitData(histData);

  InitHist(QCD,"","",TColor::GetColor(192,232,100),1001);
  InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(EWK,"","",TColor::GetColor("#DE5A6A"),1001);
  InitHist(ZLL,"","",TColor::GetColor("#4496C8"),1001);
  InitHist(ZTT,"","",TColor::GetColor("#FFCC66"),1001);

  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.6);
  histData->GetYaxis()->SetTitleSize(0.06);
  float yUpper = histData->GetMaximum();
  if (ZTT->GetMaximum()>yUpper)
    yUpper = ZTT->GetMaximum();
  histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);
  ZTT->GetYaxis()->SetRangeUser(0,1.2*ZTT->GetMaximum());
  if (logY) {
    histData->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
    ZTT->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
  }

  histData->SetMarkerSize(1.4);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  
  histData->GetXaxis()->SetLabelSize(0);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.2);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);

  histData->Draw("e1");
  ZTT->Draw("sameh");
  ZLL->Draw("sameh");
  QCD->Draw("sameh");
  EWK->Draw("sameh");
  TT->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  histData->Draw("e1same");
  float chi2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = ZTT->GetBinContent(iB);
    float eMC = bkgdErr->GetBinError(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      float e2 = xData + eMC*eMC;
      chi2 += diff2/e2;
    }
  }
  double prob = TMath::Prob(chi2,double(nBins));
  double chi2norm = chi2/double(nBins);
  std::cout << std::endl;
  std::cout << "chi2/ndof = " << chi2 << "/" << nBins << " = " << chi2norm 
	    << "  ->  p-value = " << prob << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.61,0.47,0.85,0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  leg->AddEntry(histData,"Data","lp");
  if (embedded)
    leg->AddEntry(ZTT,"embedded #tau#tau","f");
  else
    leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
  if (channel=="et") 
    leg->AddEntry(ZLL,"Z#rightarrowee","f");
  else
    leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(QCD,"j#rightarrow#tau misId","f");
  leg->AddEntry(EWK,"electroweak","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  if (plotLegend) leg->Draw();
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);
  if (logX) upper->SetLogx(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.401,1.599);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.7);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  for (int iB=1; iB<=nBins; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = ZTT->GetBinContent(iB);
    ratioErrH->SetBinContent(iB,1.0);
    ratioErrH->SetBinError(iB,0.0);
    float xBkg = bkgdErr->GetBinContent(iB);
    float errBkg = bkgdErr->GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);

    }
    if (x1>0&&x2>0) {
      float e1 = histData->GetBinError(iB);
      float ratio = x1/x2;
      float eratio = e1/x2;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }

  // ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.32);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.2);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);

  ratioH->Draw("e1");
  ratioErrH->Draw("e2same");
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  if (logX) lower->SetLogx(true);
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Update();
  if (pageType<0)
    canv1->Print("figures/"+channel+"_"+era+suffix+triggerOption+".pdf(","pdf");
  else if (pageType>0) 
    canv1->Print("figures/"+channel+"_"+era+suffix+triggerOption+".pdf)","pdf");
  else
    canv1->Print("figures/"+channel+"_"+era+suffix+triggerOption+".pdf","pdf");

  canv1->Print("figures/"+channel+"_"+Variable+"_"+era+suffix+".pdf","pdf");
  canv1->Print("figures/"+channel+"_"+Variable+"_"+era+suffix+".png");

  for (unsigned int i=0; i<sampleNames.size(); ++i) {
    delete hist_Sig[i];
    delete hist_FF[i];
  }
  delete histData;
  delete TT;
  delete EWK;
  delete ZLL;
  delete ZTT;
  delete QCD;
  delete bkgdErr;
  delete ratioH;
  delete ratioErrH;
  delete canv1;
  std::cout << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;
    
}

void warningMessage_exit() {
  std::cout << std::endl;
  std::cout << "Usage : PlotDNN [channel: et, mt] [era: 2016, 2017, 2018] [zttModel: embedded, MC] [trigger: comb, singleLep, lepTau, resolved]" << std::endl;
  exit(-1);
}

int main(int argc, char **argv) {

  if (argc!=5) 
    warningMessage_exit();
  TString channel(argv[1]);
  TString era(argv[2]);
  TString zttModel(argv[3]);
  TString trigger(argv[4]);
  if (channel!="et"&&channel!="mt") 
    warningMessage_exit();
  if (era!="2016"&&era!="2017"&&era!="2018") 
    warningMessage_exit();
  if (zttModel!="embedded"&&zttModel!="MC") 
    warningMessage_exit();
  if (trigger!="comb"&&trigger!="singleLep"&&trigger!="lepTau"&&trigger!="resolved") 
    warningMessage_exit();

  bool isEmbedded = (zttModel=="embedded");
  int itrigger = 0;
  if (trigger=="singleLep") itrigger = 1;
  else if (trigger=="lepTau") itrigger = 2;
  else if (trigger=="resolved") itrigger = 3;

  TString dir = "/nfs/dust/cms/user/filatovo/HTT/CMSSW_10_2_16/src/mlFramework/In_Tuples_"+era+"/et/21July";

  std::vector<TString> sampleNames;
  if (channel=="et") 
    sampleNames.push_back("SingleElectron");
  else
    sampleNames.push_back("SingleMuon");
  sampleNames.push_back("TTbar");
  sampleNames.push_back("Diboson");
  sampleNames.push_back("WJets");
  sampleNames.push_back("DYJets");
  if (isEmbedded) {
    if (channel=="et")
      sampleNames.push_back("EmbeddedElTau");
    else 
      sampleNames.push_back("EmbeddedMuTau");
  }
  else 
    sampleNames.push_back("DYJets");
  std::vector<TString> samplesMap = {
    "Data", "TTbar", "Diboson", "WJets", "ZLL", "ZTT" 
  };
  MapOfTrees trees;
  for (unsigned int i=0; i<samplesMap.size(); ++i) {
    TFile * file = new TFile(dir+"/et-NOMINAL_ntuple_"+sampleNames.at(i)+"_"+era+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    trees[samplesMap.at(i)] = tree;    
  }
  std::ifstream ifs ("variables.txt", std::ifstream::in);
  TString variable;
  TString xtitle;
  TString ytitle;
  int nBins;
  float xmin;
  float xmax;
  float yLower;
  float scaleYUpper;
  bool logY;
  bool logX;
  bool plotLegend;
  bool addGeV;
  std::vector<HistAttr> variablesToPlot;
  while (ifs >> variable >> xtitle >> ytitle >> addGeV >> nBins >> xmin >> xmax >> yLower >> scaleYUpper >> logX >> logY >> plotLegend) {
    HistAttr histAttr;
    histAttr.variable = variable;
    if (addGeV) {
      histAttr.xtitle = xtitle + " [GeV]";
      histAttr.ytitle = ytitle + " GeV";
    }
    else {
      histAttr.xtitle = xtitle;
      histAttr.ytitle = ytitle;
    }
    histAttr.nBins = nBins;
    histAttr.xmin = xmin;
    histAttr.xmax = xmax;
    histAttr.yLower = yLower;
    histAttr.scaleYUpper = scaleYUpper;
    histAttr.logX = logX;
    histAttr.logY = logY;
    histAttr.plotLegend = plotLegend;
    variablesToPlot.push_back(histAttr);
  }

  for (unsigned int iV=0; iV<variablesToPlot.size(); ++iV) {
    int pageType = 0;
    unsigned int nMax = variablesToPlot.size() - 1;
    if (iV==0) pageType = -1;
    if (iV==nMax) pageType = 1;
    PlotDNN(channel,trees,variablesToPlot.at(iV),isEmbedded,era,itrigger,pageType);
  }



}
