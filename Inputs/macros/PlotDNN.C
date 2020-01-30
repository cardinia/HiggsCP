#include "HttStylesNew.cc"
#include "CMS_lumi.C"

void PlotDNN( bool embedded = true,
	      TString era = "2018",
	      bool method = true) {



  TString suffix = "";
  if (embedded) suffix = "_embedded";

  TString dir = "/nfs/dust/cms/user/rasp/HiggsCP/2018/DNN/";
  double qcd_scale = 1.2;

  TString DataFile = "SingleMuon";
  TString Variable = "dijetpt";
  TString xtitle = "dijet p_{T} [GeV]";
  TString ytitle = "Events";
  int nBins  =                  20;
  float xmin =                   0;
  float xmax =                 400;
  float yLower =                 0;
  float scaleYUpper =           10;

  TString Cuts("pt_1>21&&os>0.5&&puppimt_1<50&&njets>=2");

  TString Cuts_Sig    = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2>0.5");
  TString Cuts_FF     = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5");

  TString CutsZTT_Sig = Cuts_Sig + TString("&&gen_match_1==4&&gen_match_2==5");
  TString CutsZTT_FF  = Cuts_FF  + TString("&&gen_match_1==4&&gen_match_2==5");

  TString CutsZLL_Sig = Cuts_Sig + TString("&&!(gen_match_1==4&&gen_match_2==5)");
  TString CutsZLL_FF  = Cuts_FF  + TString("&&!(gen_match_1==4&&gen_match_2==5)");

  SetStyle();

  bool logY = false;
  bool logX = false;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TString sampleNames[10] = {
    "SingleMuon",   // data (0)
    "TTbar",        // TTbar     (1)
    "SingleTop",    // SingleTop (2)
    "Diboson",      // Dibosons  (3)
    "WJets",        // WJets     (4)
    "DYJets",       // DYJets    (5)
    "DYJets"        // DYJets    (6)
    "",
    "",
    ""
  };
  if (embedded)
    sampleNames[6] = TString("EmbeddedMuTau");

  TString cuts_Sig[10];
  TString cuts_FF[10];
  for (int i=0; i<10; ++i) {
    if (embedded) {
      cuts_Sig[i] = "weight*("+CutsZLL_Sig+")";
      cuts_FF[i]  = "weight*ff_nom*("+CutsZLL_FF+"&&gen_match_2!=6)";
    }
    else {
      cuts_Sig[i] = "weight*("+Cuts_Sig+")";
      cuts_FF[i]  = "weight*ff_nom*("+Cuts_FF+"&&gen_match_2!=6)";
    }
  }
  cuts_Sig[0] = Cuts_Sig;
  cuts_FF[0] = "ff_nom*("+Cuts_FF+")";

  if (embedded) {
    cuts_Sig[6] = "weight*("+Cuts_Sig+")";
    cuts_FF[6]  = "weight*ff_nom*("+Cuts_FF+"&&gen_match_2!=6)";
  }
  else {
    cuts_Sig[5] = "weight*("+CutsZLL_Sig+")";
    cuts_FF[5]  = "weight*ff_nom*("+CutsZLL_FF+"&&gen_match_2!=6)";
    cuts_Sig[6] = "weight*("+CutsZTT_Sig+")";
    cuts_FF[6]  = "weight*ff_nom*("+CutsZTT_FF+"&&gen_match_2!=6)";
  }

  int nSamples = 7;

  TH1D * hist_Sig[10];
  TH1D * hist_FF[10];

  TCanvas * dummy = new TCanvas("dummy","",400,400);

  // filling histograms
  for (int i=0; i<nSamples; ++i) {
    cout << sampleNames[i] << endl;
    TFile * file = new TFile(dir+"/mt-NOMINAL_ntuple_"+sampleNames[i]+"_"+era+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
    TString histNameSig = sampleNames[i] + "_sig";
    TString histNameFF  = sampleNames[i] + "_ff";
    hist_Sig[i] = new TH1D(histNameSig,"",nBins,xmin,xmax);
    hist_Sig[i]->Sumw2();
    hist_FF[i] = new TH1D(histNameFF,"",nBins,xmin,xmax);
    hist_FF[i]->Sumw2();
    tree->Draw(Variable+">>"+histNameSig,cuts_Sig[i]);
    tree->Draw(Variable+">>"+histNameFF,cuts_FF[i]);
  }
  delete dummy;
  
  // forming Fakes
  for (int i=1; i<nSamples; ++i) {
    hist_FF[0]->Add(hist_FF[0],hist_FF[i],1,-1);
  }

  // summing up dibosons and single-top into ewk
  hist_Sig[2]->Add(hist_Sig[2],hist_Sig[3],1,1);

  TH1D * histData = (TH1D*)hist_Sig[0]->Clone("data_obs");
  TH1D * TT       = (TH1D*)hist_Sig[1]->Clone("TT");
  TH1D * EWK      = (TH1D*)hist_Sig[2]->Clone("EWK");
  TH1D * ZLL      = (TH1D*)hist_Sig[5]->Clone("ZLL");
  TH1D * ZTT      = (TH1D*)hist_Sig[6]->Clone("ZTT");
  TH1D * QCD      = (TH1D*)hist_FF[0]->Clone("fakes");

  std::cout << "Top : " << TT->GetSumOfWeights() << std::endl;
  std::cout << "EWK : " << EWK->GetSumOfWeights() << std::endl;
  std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
  std::cout << "ZLL : " << ZLL->GetSumOfWeights() << std::endl;
  std::cout << "ZTT : " << ZTT->GetSumOfWeights() << std::endl;

  EWK->Add(EWK,TT);
  QCD->Add(QCD,EWK);
  ZLL->Add(ZLL,QCD);
  ZTT->Add(ZTT,ZLL);

  std::cout << std::endl;
  std::cout << "MC   : " << ZTT->GetSumOfWeights() << std::endl;
  std::cout << "DATA : " << histData->GetSumOfWeights() << std::endl;

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

  InitHist(QCD,TColor::GetColor("#FFCCFF"));
  InitHist(TT,TColor::GetColor("#9999CC"));
  InitHist(EWK,TColor::GetColor("#DE5A6A"));
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
  EWK->Draw("sameh");
  TT->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  histData->Draw("e1same");
  float chi2 = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = ZTT->GetBinContent(iB);
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
  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(ZTT,"embedded Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(QCD,"fakes","f");
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