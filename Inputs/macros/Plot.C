#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include "Plotting.h"
#include "Plotting_Style.h"
#include "settings.h"
void Plot( bool embedded = true) {

  TString suffix = "";
  if (embedded) suffix = "_embedded";

  TString dir = "/nfs/dust/cms/user/rasp/HiggsCP/2018/";
  double qcd_scale = 1.2;

  TString DataFile = "SingleMuon_Run2018";
  TString Variable = "IP_signif_PV_with_BS_2";
  TString xtitle = "pion IP_{sig}";
  TString ytitle = "Events";
  int nBins  =                50;
  float xmin =                 0;
  float xmax =                10;
  float yLower =       0;
  float scaleYUpper = 10;

  TString Weight("puweight*mcweight*effweight*");
  TString WeightEmb("mcweight*embweight*effweight*");
  TString Cuts("((trg_singlemuon>0.5&&pt_1>25)||(trg_mutaucross&&pt_1>21&&pt_2>32))&&iso_1<0.15&&pt_1>21&&pt_2>20&&byTightDeepTau2017v2p1VSmu_2>0.5&&byVLooseDeepTau2017v2p1VSe_2&&extraelec_veto<0.5&&extramuon_veto<0.5&&dilepton_veto<0.5&&byMediumDeepTau2017v2p1VSjet_2>0.5&&puppimt_1<50");

  TString CutsOS = Cuts + TString("&&os>0.5");
  TString CutsSS = Cuts + TString("&&os<0.5");
  TString CutsZTT_OS  = CutsOS + TString("&&gen_match_1==4&&gen_match_2==5");
  TString CutsZLL_OS  = CutsOS + TString("&&!(gen_match_1==4&&gen_match_2==5)");
  TString CutsZTT_SS  = CutsSS + TString("&&gen_match_1==4&&gen_match_2==5");
  TString CutsZLL_SS  = CutsSS + TString("&&!(gen_match_1==4&&gen_match_2==5)");

  SetStyle();

  bool logY = false;
  bool logX = false;

  double lumi = 59740;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TString sampleNames[26] = {
    DataFile, // data (0)
    "TTTo2L2Nu",        // TTTo2L2Nu        (1)
    "TTToSemiLeptonic", // TTToSemiLeptonic (2)
    "TTToHadronic",     // TTToHadronic     (3)
    "ST_t-channel_antitop_4f",// ST_t_atop  (4)
    "ST_t-channel_top_4f",    // ST_t_top   (5)
    "ST_tW_antitop_5f",       // ST_tW_atop (6)
    "ST_tW_top_5f",           // ST_tW_top  (7)
    "WW",                     // WW         (8) 
    "WZ",                     // WZ         (9)
    "ZZ",                     // ZZ        (10)
    "WJetsToLNu",    // WJets inclusive    (11)
    "W1JetsToLNu",   // W1Jets inclusive   (12)
    "W2JetsToLNu",   // W2Jets inclusive   (13)
    "W3JetsToLNu",   // W3Jets inclusive   (14)
    "W4JetsToLNu",   // W4Jets inclusive   (15)
    "DYJetsToLL_M-50" , // DYJetsToLL      (16)
    "DY1JetsToLL_M-50", // DY1JetsToLL     (17)
    "DY2JetsToLL_M-50", // DY2JetsToLL     (18)
    "DY3JetsToLL_M-50", // DY3JetsToLL     (19)
    "DY4JetsToLL_M-50", // DY4JetsToLL     (20)
    "DYJetsToLL_M-50" , // DYJetsToLL      (21)
    "DY1JetsToLL_M-50", // DY1JetsToLL     (22)
    "DY2JetsToLL_M-50", // DY2JetsToLL     (23)
    "DY3JetsToLL_M-50", // DY3JetsToLL     (24)
    "DY4JetsToLL_M-50"  // DY4JetsToLL     (25)
  };
  if (embedded)
    sampleNames[21] = TString("Embedded_Run2018");

  TString cuts[30];
  TString cutsSS[30];
  for (int i=0; i<30; ++i) {
    if (embedded) {
      cuts[i] = Weight+"("+CutsZLL_OS+")";
      cutsSS[i] = Weight+"("+CutsZLL_SS+")";
      if (i==11) {
	cuts[i] = Weight+"("+CutsZLL_OS+"&&gen_noutgoing==0)";
	cutsSS[i] = Weight+"("+CutsZLL_SS+"&&gen_noutgoing==0)";
      }
    }
    else {
      cuts[i] = Weight+"("+CutsOS+")";
      cutsSS[i] = Weight+"("+CutsSS+")";
      if (i==11) {
	cuts[i] = Weight+"("+CutsOS+"&&gen_noutgoing==0)";
	cutsSS[i] = Weight+"("+CutsSS+"&&gen_noutgoing==0)";
      }
    }
  }
  cuts[0] = CutsOS;
  cutsSS[0] = CutsSS;

  int nSamples = 26;
  for (int i=16; i<=20; ++i) {
    cuts[i] = Weight+"("+CutsZLL_OS+")";
    cutsSS[i] = Weight+"("+CutsZLL_SS+")";
    if (i==16) {
      cuts[i] = Weight+"("+CutsZLL_OS+"&&gen_noutgoing==0)";
      cutsSS[i] = Weight+"("+CutsZLL_SS+"&&gen_noutgoing==0)";
    }
  }
  for (int i=21; i<=25; ++i) {
    cuts[i] = Weight+"("+CutsZTT_OS+")";
    cutsSS[i] = Weight+"("+CutsZTT_SS+")";
    if (i==21) {
      cuts[i] = Weight+"("+CutsZTT_OS+"&&gen_noutgoing==0)";
      cutsSS[i] = Weight+"("+CutsZTT_SS+"&&gen_noutgoing==0)";
    }
  }

  if (embedded) {
    nSamples = 22;
    cuts[21] = WeightEmb+"("+CutsOS+")";
    cutsSS[21] = WeightEmb+"("+CutsSS+")";
  }

  TH1D * hist[30];
  TH1D * histSS[30];

  TCanvas * dummy = new TCanvas("dummy","",400,400);

  // filling histograms
  for (int i=0; i<nSamples; ++i) {
    TFile * file = new TFile(dir+"/"+sampleNames[i]+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    TH1D * histWeightsH = (TH1D*)file->Get("nWeightedEvents");
    double norm = 1; 
    double nevents = 1;
    if (i==21&&embedded) { 
      norm = 1.0;
    }
    else {
      double xsec = xsecs[sampleNames[i]];
      nevents = histWeightsH->GetSumOfWeights();;
      norm = xsec*lumi/nevents;
    }
    TString histName = sampleNames[i];
    TString histNameSS = sampleNames[i] + "_ss";
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    hist[i]->Sumw2();
    histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);
    histSS[i]->Sumw2();
    //    cout << cuts[i] << endl;
    tree->Draw(Variable+">>"+histName,cuts[i]);
    tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
    std::cout << sampleNames[i] 
	      << " : norm = " << norm
      	      << " : events = " << hist[i]->GetSumOfWeights() 
	      << " : initial events from settings = " << nevents << std::endl;
    if (i>0) {
      for (int iB=1; iB<=nBins; ++iB) {

	double x = hist[i]->GetBinContent(iB);
	double e = hist[i]->GetBinError(iB);
    	hist[i]->SetBinContent(iB,norm*x);
    	hist[i]->SetBinError(iB,norm*e);

	x = histSS[i]->GetBinContent(iB);
        e = histSS[i]->GetBinError(iB);
        histSS[i]->SetBinContent(iB,norm*x);
        histSS[i]->SetBinError(iB,norm*e);

      }
    }
  }

  delete dummy;
  
  // adding up SS
  for (int i=1; i<nSamples; ++i) {
    histSS[0]->Add(histSS[0],histSS[i],1,-1);
  }

  for (int iB=1; iB<=nBins; ++iB) {
    histSS[0]->SetBinContent(iB,qcd_scale*histSS[0]->GetBinContent(iB));
    histSS[0]->SetBinError(iB,qcd_scale*histSS[0]->GetBinError(iB));
  }

  // adding top samples
  for (int i=2; i<=3; ++i) {
    hist[1]->Add(hist[1],hist[i]);
  }

  // adding up electroweak samples
  for (int i=5; i<=10; ++i) {
    hist[4]->Add(hist[4],hist[i]);
  }
  // adding W samples 
  for (int i=12; i<=15; ++i) {
    hist[11]->Add(hist[11],hist[i]);
  }
  // adding ZLL
  for (int i=17; i<=20; ++i) {
    hist[16]->Add(hist[16],hist[i]);
  }

  // adding ZTT
  if (!embedded) {
    for (int i=22; i<=25; ++i) {
      hist[21]->Add(hist[21],hist[i]);
    }
  }

  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * W        = (TH1D*)hist[11]->Clone("W");
  TH1D * TT       = (TH1D*)hist[1]->Clone("TT");
  TH1D * EWK      = (TH1D*)hist[4]->Clone("EWK");
  TH1D * ZLL      = (TH1D*)hist[16]->Clone("ZLL");
  TH1D * ZTT      = (TH1D*)hist[21]->Clone("ZTT");
  TH1D * QCD      = (TH1D*)histSS[0]->Clone("QCD");

  std::cout << "Top : " << TT->GetSumOfWeights() << std::endl;
  std::cout << "EWK : " << EWK->GetSumOfWeights() << std::endl;
  std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
  std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << std::endl;

  EWK->Add(EWK,TT);
  W->Add(W,EWK);
  QCD->Add(QCD,W);
  ZLL->Add(ZLL,QCD);
  ZTT->Add(ZTT,ZLL);

  std::cout << std::endl;
  std::cout << "MC : " << ZTT->GetSumOfWeights() << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << std::endl;

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  InitData(histData);

  InitHist(QCD,TColor::GetColor("#FFCCFF"));
  InitHist(TT,TColor::GetColor("#9999CC"));
  InitHist(EWK,TColor::GetColor("#6F2D35"));
  InitHist(W,TColor::GetColor("#DE5A6A"));
  InitHist(ZLL,TColor::GetColor("#4496C8"));
  InitHist(ZTT,TColor::GetColor("#FFCC66"));

  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  float yUpper = histData->GetMaximum();
  if (ZTT->GetMaximum()>yUpper)
    yUpper = ZTT->GetMaximum();
  histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);
  ZTT->GetYaxis()->SetRangeUser(0,1.2*ZTT->GetMaximum());
  if (logY) {
    histData->GetYaxis()->SetRangeUser(yLower,scaleYUpper*yUpper);
  }

  histData->SetMarkerSize(1.4);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
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
  ZTT->Draw("h");
  ZLL->Draw("sameh");
  QCD->Draw("sameh");
  W->Draw("sameh");
  EWK->Draw("sameh");
  TT->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  histData->Draw("e1same");
  float chi2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = W->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.55,0.48,0.8,0.78);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  //  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(ZTT,"embedded Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(QCD,"QCD","f");
  leg->AddEntry(W,"W#rightarrow#mu#nu","f");
  leg->AddEntry(EWK,"electroweak","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  leg->Draw();
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

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
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
  ratioH->GetYaxis()->SetTitleOffset(0.5);
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
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
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
  canv1->Print("figures/plot_"+Variable+suffix+".png");

}
