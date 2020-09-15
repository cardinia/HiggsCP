#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include "CorrectZTT.C"
// embedded = true (EMB), false (DY MC)
// era = 2016, 2017, 2018
// acotautau_refitbs_00, acotautau_refitbs_01
// alphaminus, alpha_IP_1, alpha_IP_1, alpha_plane_2, alpha_plane_2
void PlotDNN( bool embedded = true,
	      TString era = "2017"
	      ) {

  TString Variable = "alpha_plane_2";
  TString xtitle = "#alpha_{#rho}";
  //  TString Variable = "alpha_plane_2";
  //  TString xtitle = "#alpha_{#rho}";
  //  TString Variable = "alpha_plane_2";
  //  TString xtitle = "#alpha_{#rho}";
  //  TString Variable = "m_vis";
  //  TString xtitle = "m_{vis} [GeV]";
  TString ytitle = "Events / 0.079";
  int nBins  =                   20;
  float xmin =                    0;
  float xmax =      0.5*TMath::Pi();
  float yLower =                  1;
  float scaleYUpper =            10;
  bool logY = false;
  bool logX = false;
  bool plotLegend = true;

  double fake_scale = 1.0;
  double zll_scale = 1.0;
  double ztt_scale = 1.0;
void PlotDNN( bool embedded = true,
	      TString era = "2017") {

  SetStyle();

  TString suffix = "_mc";
  if (embedded) suffix = "_embedded";


  TString dir = "/nfs/dust/cms/user/filatovo/HTT/CMSSW_10_2_16/src/mlFramework/In_Tuples_"+era+"/et/6July";

  lumi_13TeV = "2018, 59.7 fb^{-1}";
  if (era=="2017")
    lumi_13TeV = "2017, 41.0 fb^{-1}";
  if (era=="2016")
    lumi_13TeV = "2016, 35.9 fb^{-1}";


  bool nonEquidistant = false;
  float bins[3] = {0.0,1.5,10.0};
  if (Variable.Contains("IP_signif")&&Variable.Contains("_2")) {
    nBins = 2;
    nonEquidistant = true;
  }
  TString DataFile = "SingleMuon";
  TString Variable = "m_vis";
  TString xtitle = "m_{vis} [GeV]";
  TString ytitle = "Events";
  int nBins  =                  30;
  float xmin =                   0;
  float xmax =                 300;
  float yLower =                 0;
  float scaleYUpper =           10;
  bool logY = false;
  bool logX = false;


  TString Cuts("pt_1>21&&os>0.5&&puppimt_1<50");

  TString Cuts_Sig    = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2>0.5");
  TString Cuts_FF     = Cuts + TString("&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5");

  TString CutsZTT_Sig = Cuts_Sig + TString("&&gen_match_1==4&&gen_match_2==5");
  TString CutsZTT_FF  = Cuts_FF  + TString("&&gen_match_1==4&&gen_match_2==5");

  TString CutsZLL_Sig = Cuts_Sig + TString("&&!(gen_match_1==4&&gen_match_2==5)");
  TString CutsZLL_FF  = Cuts_FF  + TString("&&!(gen_match_1==4&&gen_match_2==5)");

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TString sampleNames[10] = {
    "SingleElectron",  // data  (0)
    "TTbar",       // TT    (1)
    "Diboson",     // VV    (2)
    "WJets",       // W     (3)
    "DYJets",      // Z->ll (4)
    "DYJets",      // Z->tt or EMB (5)
    "",
    "",
    "",
    ""
  };

  if (embedded)
    sampleNames[5] = TString("EmbeddedElTau");

  TString cuts_Sig[10];
  TString cuts_FF[10];
  for (int i=0; i<10; ++i) {
    if (embedded) {
      cuts_Sig[i] = "weight*("+CutsZLL_Sig+"&&gen_match_2!=6)";
      cuts_FF[i]  = "weight*ff_mva*("+CutsZLL_FF+"&&gen_match_2!=6)";
    }
    else {
      cuts_Sig[i] = "weight*("+Cuts_Sig+"&&gen_match_2!=6)";
      cuts_FF[i]  = "weight*ff_mva*("+Cuts_FF+"&&gen_match_2!=6)";
    }
  }
  cuts_Sig[0] = Cuts_Sig;
  cuts_FF[0] = "ff_mva*("+Cuts_FF+")";

  if (embedded) {
    cuts_Sig[5] = "weight*("+Cuts_Sig+"&&gen_match_2!=6)";
    cuts_FF[5]  = "weight*ff_mva*("+Cuts_FF+"&&gen_match_2!=6)";
  }
  else {
    cuts_Sig[4] = "weight*("+CutsZLL_Sig+"&&gen_match_2!=6)";
    cuts_FF[4]  = "weight*ff_mva*("+CutsZLL_FF+"&&gen_match_2!=6)";
    cuts_Sig[5] = "weight*("+CutsZTT_Sig+"&&gen_match_2!=6)";
    cuts_FF[5]  = "weight*ff_mva*("+CutsZTT_FF+"&&gen_match_2!=6)";
  }

  int nSamples = 6;

  TH1D * hist_Sig[10];
  TH1D * hist_FF[10];

  TCanvas * dummy = new TCanvas("dummy","",400,400);

  // filling histograms
  for (int i=0; i<nSamples; ++i) {
    TFile * file = new TFile(dir+"/mt-NOMINAL_ntuple_"+sampleNames[i]+"_"+era+".root");
    TTree * tree = (TTree*)file->Get("TauCheck");
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
    cout << cuts_Sig[i] << endl;
    tree->Draw(Variable+">>"+histNameFF,cuts_FF[i]);
  }
  delete dummy;
  
  // forming Fakes
  for (int i=1; i<nSamples; ++i) {
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
  double ZTT_norm = 0.06; //  normalization ZTT :  4% (EMBEDDED)
  double EWK_norm = 0.07; //  normalization EWK :  7%
  double QCD_norm = 0.15; //  normalization Fakes : 15%
  double ZLL_mtau = 0.10; //  mu->tau fake rate ZLL : 10%
  double TT_norm  = 0.07; //  normalization TT  :  7%
  double mu_ID    = 0.02; //  muon trigger ID   :  2%
  double tau_ID   = 0.05; //  tau ID            :  5%

  for (int iB=1; iB<=nBins; ++iB) {

    float ztt  = ZTT->GetBinContent(iB);
    float ztte = ZTT->GetBinError(iB);
    ztte = TMath::Sqrt(ztte*ztte+ztt*ztt*(ZTT_norm*ZTT_norm+mu_ID*mu_ID+tau_ID*tau_ID));
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
    zlle = TMath::Sqrt(zlle*zlle+zll*zll*(mu_ID*mu_ID+ZLL_mtau*ZLL_mtau));
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
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.61,0.47,0.85,0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  leg->AddEntry(histData,"Data","lp");
  if (embedded)
    leg->AddEntry(ZTT,"embedded Z#rightarrow#tau#tau","f");
  else
    leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(QCD,"misidentified j#rightarrow#tau","f");
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
  canv1->Print("figures/"+Variable+"_"+era+"_"+DecayMode[decayMode]+alpha[alphaIndex]+suffix+"_"+era+".png");
  canv1->Print("figures/"+Variable+"_"+era+"_"+DecayMode[decayMode]+alpha[alphaIndex]+suffix+"_"+era+".C");

  TFile * fileCanvas = new TFile("figures/"+Variable+"_"+era+"_"+DecayMode[decayMode]+alpha[alphaIndex]+suffix+"_"+era+".root","recreate");
  fileCanvas->cd("");
  canv1->Write("canvas");
  fileCanvas->Close();


}
