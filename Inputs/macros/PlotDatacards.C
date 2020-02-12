#include "HttStylesNew.cc"

// mt_2017_higgs_mupi_IP 
// mt_2017_ztt_mupi_IP
// mt_2017_top_mupi_IP
// mt_2017_fakes_mupi_IP

// mt_2017_higgs_murho_DP
// mt_2017_ztt_murho_DP
// mt_2017_top_murho_DP
// mt_2017_fakes_murho_DP

void PlotDatacards(TString dir = "/.",
		   bool do_fakes = false,
		   TString cat = "mt_murho_sig",
		   bool blindData = false) {

  SetStyle();
  //  gStyle->SetErrorX(0);

  TString xtitle = "bins";
  TString ytitle = "Events";
  TString inputFileName = "output_oleg_ip.root";

  TFile * inputs = new TFile(inputFileName);

  bool logY = false;

  TH1F * histData = (TH1F*)inputs->Get(cat+"/data_obs");
  TH1F * Fakes    = (TH1F*)inputs->Get(cat+"/QCD");
  if (do_fakes)
    Fakes    = (TH1F*)inputs->Get(cat+"/fakes");
  TH1F * W        = (TH1F*)inputs->Get(cat+"/W");
  TH1F * VV       = (TH1F*)inputs->Get(cat+"/VV");
  TH1F * ST       = (TH1F*)inputs->Get(cat+"/ST");
  TH1F * EMB      = (TH1F*)inputs->Get(cat+"/EMB");
  TH1F * Zll      = (TH1F*)inputs->Get(cat+"/ZLL");
  TH1F * TT       = (TH1F*)inputs->Get(cat+"/TT");

  TH1F * hSM      = (TH1F*)inputs->Get(cat+"/ggH_sm_htt125");
  TH1F * hPS      = (TH1F*)inputs->Get(cat+"/ggH_ps_htt125");
  TH1F * hMM      = (TH1F*)inputs->Get(cat+"/ggH_mm_htt125");
  TH1F * qqSM      = (TH1F*)inputs->Get(cat+"/qqH_sm_htt125");
  TH1F * qqPS      = (TH1F*)inputs->Get(cat+"/qqH_ps_htt125");
  TH1F * qqMM      = (TH1F*)inputs->Get(cat+"/qqH_mm_htt125");

  Fakes->Scale(0.95);
  hSM->Add(hSM,qqSM,1.,1.);
  hPS->Add(hPS,qqPS,1.,1.);
  hMM->Add(hMM,qqMM,1.,1.);

  cout << "OK" << endl;

  float nEMB   = EMB->GetSumOfWeights();
  float nFakes = Fakes->GetSumOfWeights();
  float nW     = W->GetSumOfWeights();
  float nZll   = Zll->GetSumOfWeights();
  float nST    = ST->GetSumOfWeights();
  float nVV    = VV->GetSumOfWeights();
  float nTT    = TT->GetSumOfWeights();
  float nData = histData->GetSumOfWeights();
  float nRest = nZll + nST + nVV + nTT;
  if (!do_fakes) nFakes += nW;
  
  float nBkgd = nEMB + nFakes + nRest;
  float nSM   = hSM->GetSumOfWeights();
  float nPS   = hPS->GetSumOfWeights();
  
  std::cout << "W+Jets : " << nW << std::endl;    
  std::cout << "Fakes  : " << nFakes << std::endl;
  std::cout << "ZTT    : " << nEMB << std::endl;
  std::cout << "Rest   : " << nRest << std::endl;
  std::cout << "Bkgd   : " << nBkgd << std::endl;
  std::cout << "Data   : " << nData << std::endl;

  hSM->Scale(20);
  hPS->Scale(20);
  hMM->Scale(20);

  int nBins = histData->GetNbinsX();
  
  double muUnc = 0.02;
  double tauUnc = 0.05;
  double embUnc = 0.06;
  double ttUnc = 0.07;
  double fakeUnc = 0.15;
  double zllUnc = 0.1;
  double wUnc   = 0.20;
  for (int i=1; i<=nBins; ++i) {

    double ttE = TT->GetBinError(i);
    double tt = TT->GetBinContent(i);
    ttE = TMath::Sqrt(ttE*ttE+ttUnc*ttUnc*tt*tt);
    TT->SetBinError(i,ttE);

    double embE = EMB->GetBinError(i);
    double emb = EMB->GetBinContent(i);
    embE = TMath::Sqrt(embE*embE+emb*emb*(muUnc*muUnc+tauUnc*tauUnc+embUnc*embUnc));
    EMB->SetBinError(i,embE);

    double fakeE = Fakes->GetBinError(i);
    double fake = Fakes->GetBinContent(i);
    fakeE = TMath::Sqrt(fakeE*fakeE+fakeUnc*fakeUnc*fake*fake);
    Fakes->SetBinError(i,fakeE);

    double zllE = Zll->GetBinError(i);
    double zll = Zll->GetBinContent(i);
    zllE = TMath::Sqrt(zllE*zllE+zllUnc*zllUnc*zll*zll);
    Zll->SetBinError(i,zllE);

    double wE = W->GetBinError(i);
    double w = W->GetBinContent(i);
    wE = TMath::Sqrt(wE*wE+wUnc*wUnc*w*w);
    W->SetBinError(i,wE);


  }

  TT->Add(TT,Zll);
  TT->Add(TT,ST);
  TT->Add(TT,VV);
  
  Fakes->Add(Fakes,TT);
  if (!do_fakes)
    Fakes->Add(Fakes,W);
  EMB->Add(EMB,Fakes);

  hSM->SetLineColor(2);
  hPS->SetLineColor(4);
  hMM->SetLineColor(3);
  hSM->SetLineWidth(3);
  hPS->SetLineWidth(3);
  hMM->SetLineWidth(3);

  cout << "Tot background after combining : " << EMB->GetSumOfWeights() << endl;

  if (cat.Contains("_ztt")) {
    histData->GetXaxis()->SetRangeUser(0,17.5);
    EMB->GetXaxis()->SetRangeUser(0,17.5);
  }
  TH1F * bkgdErr = (TH1F*)EMB->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  

  float ymax = histData->GetMaximum();
  if (EMB->GetMaximum()>ymax)
    ymax = EMB->GetMaximum();

  histData->GetYaxis()->SetRangeUser(0,1.2*ymax);

  for (int iB=1; iB<=nBins; ++iB) {
  //    hSM->SetBinError(iB,0);
  //    hPS->SetBinError(iB,0);
    TT->SetBinError(iB,0);
    Fakes->SetBinError(iB,0);
    EMB->SetBinError(iB,0);
  }
  InitData(histData);
  InitHist(TT,"","",TColor::GetColor("#6F2D35"),1001);
  InitHist(Fakes,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(EMB,"","",TColor::GetColor("#FFCC66"),1001);
  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);

  histData->SetMarkerSize(1.2);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  //  nData = histData->GetSum();
  //  float nMC   = TT->GetSum();
  //  float eData = TMath::Sqrt(nData);


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

  if (blindData) {
    for (int i=1; i<=nBins; ++i) 
      histData->SetBinContent(i,-100);
  }

  histData->Draw("e1");
  EMB->Draw("sameh");
  Fakes->Draw("sameh");
  TT->Draw("sameh");
  if (blindData) {
    hSM->Draw("sameh");
    hPS->Draw("sameh");
    hMM->Draw("sameh");
  }
  if (!blindData)  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  float chi2 = 0;
  //  int nBins = histData->GetNbinsX();
  for (int iB=1; iB<=nBins; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = EMB->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.65,0.5,0.85,0.88);
  SetLegendStyle(leg);
  leg->SetTextSize(0.047);
  leg->SetHeader(cat);
  leg->AddEntry(histData,"data","lp");
  leg->AddEntry(EMB,"EMB","f");
  leg->AddEntry(Fakes,"fakes","f");
  leg->AddEntry(TT,"rest","f");
  if (blindData) { 
    leg->AddEntry(hSM,"0^{+} (x20)","l");
    leg->AddEntry(hPS,"0^{-} (x20)","l");
    leg->AddEntry(hMM,"max-mix (x20)","l");
  }
  if (blindData)
    leg->Draw();
  //  writeExtraText = true;
  //  extraText = "Preliminary";
  //  CMS_lumi(upper,4,33); 
  //  plotchannel("#tau#nu");

  if (logY) upper->SetLogy(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  cout << histData->GetNbinsX() << endl;
  float low = histData->GetXaxis()->GetBinLowEdge(1);
  float up = histData->GetXaxis()->GetBinLowEdge(nBins+1);
  TH1F * ratioH = (TH1F*)histData->Clone("ratioH");
  TH1F * ratioErrH = (TH1F*)histData->Clone("ratioErrH");
  ratioErrH->SetFillStyle(3013);
  ratioErrH->SetFillColor(1);
  ratioErrH->SetMarkerStyle(21);
  ratioErrH->SetMarkerSize(0);
  
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.501,1.499);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("Obs./Exp.");
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
    float x2 = EMB->GetBinContent(iB);
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
      if (blindData)
	ratioH->SetBinContent(iB,-100);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }

  cout << "nBins = " << ratioH->GetNbinsX() << endl;

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

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Update();
  canv1->Print(cat+".png");
  
}
