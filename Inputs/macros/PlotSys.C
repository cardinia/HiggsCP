#include "HttStylesNew.cc"
// CMS_htt_boson_reso_met_13TeV
// CMS_htt_boson_scale_met_13TeV
// CMS_scale_mu_13TeV
// CMS_scale_t_1prong_13TeV
// CMS_scale_t_1prong1pizero_13TeV
// CMS_scale_t_3prong_13TeV
// CMS_scale_j_FlavorQCD_13TeV
// CMS_scale_j_RelativeBal_13TeV
// CMS_scale_j_HF_13TeV
// CMS_scale_j_BBEC1_13TeV
// CMS_scale_j_EC2_13TeV
// CMS_scale_j_Absolute_13TeV
// CMS_scale_j_Absolute_2018_13TeV
// CMS_scale_j_HF_2018_13TeV
// CMS_scale_j_EC2_2018_13TeV
// CMS_scale_j_RelativeSample_2018_13TeV
// CMS_scale_j_BBEC1_2018_13TeV
// CMS_res_j_13TeV
// CMS_eff_b_13TeV

void PlotSys() {

  TString dirName = "/nfs/dust/cms/user/rasp/CMSSW/CMSSW_10_2_16/src/DesyTauAnalyses/NTupleMaker/test/mutau/2018";
  TString fileName = "TTTo2L2Nu";
  TString varName = "bpt_1";
  TString sysName  = "CMS_eff_b_13TeV";
  TString header = sysName;
  TString cuts = "os>0.5&&byMediumDeepTau2017v2p1VSjet_2>0.5";

  TString ytitle("Events");
  TString xtitle("leading bjet jet p_{T} [GeV]");

  int nBins =       10;
  float xmin =       0;
  float xmax =     200;

  SetStyle();
  gStyle->SetErrorX(0);
  TFile * file = new TFile(dirName+"/"+fileName+".root");

  TTree * treeNominal = (TTree*)file->Get("TauCheck");
  TTree * treeUp = (TTree*)file->Get("TauCheck_"+sysName+"Up");
  TTree * treeDown = (TTree*)file->Get("TauCheck_"+sysName+"Down");

  TH1D * histNominal = new TH1D("histNominal","",nBins,xmin,xmax);
  TH1D * histUp = new TH1D("histUp","",nBins, xmin,xmax);
  TH1D * histDown = new TH1D("histDown","",nBins, xmin,xmax);

  TCanvas * dummy = new TCanvas("dummy","",600,600);
  treeNominal->Draw(varName+">>histNominal",cuts);
  treeUp->Draw(varName+">>histUp",cuts);
  treeDown->Draw(varName+">>histDown",cuts);
  float yMax = histUp->GetMaximum();
  if (histNominal->GetMaximum()>yMax)
    yMax = histNominal->GetMaximum();
  if (histDown->GetMaximum()>yMax)
    yMax = histDown->GetMaximum();

  histUp->GetYaxis()->SetRangeUser(0.01,1.1*yMax);
  histNominal->SetLineColor(1);
  histUp->SetLineColor(2);
  histDown->SetLineColor(4);
  histDown->SetLineStyle(3);
  histNominal->SetMarkerColor(1);
  histUp->SetMarkerColor(2);
  histDown->SetMarkerColor(4);
  histNominal->SetMarkerSize(1.3);
  histNominal->GetYaxis()->SetTitle(ytitle);
  histNominal->GetXaxis()->SetTitle(xtitle);
  histUp->GetYaxis()->SetTitle(ytitle);
  histUp->GetXaxis()->SetTitle(xtitle);
  histDown->GetYaxis()->SetTitle(ytitle);
  histDown->GetXaxis()->SetTitle(xtitle);
  TH1D * ratioUp = (TH1D*)histUp->Clone("ratioUp");
  TH1D * ratioDown = (TH1D*)histDown->Clone("ratioDown");
  TH1D * ratioCentral = (TH1D*)histNominal->Clone("ratioCentral");
  //  ratioCentral->SetFillStyle(3013);
  //  ratioCentral->SetFillColor(1);
  //  ratioCentral->SetMarkerStyle(21);
  //  ratioCentral->SetMarkerSize(0);


  for (int iB=1; iB<=nBins; ++iB) {
    histUp->SetBinError(iB,0); 
    histDown->SetBinError(iB,0); 
    float xUp = histUp->GetBinContent(iB);
    float xDown = histDown->GetBinContent(iB);
    float xCentral = histNominal->GetBinContent(iB);
    float xratioUp = 1;
    float xratioDown = 1;
    if (xCentral>0) {
      xratioUp   = xUp/xCentral;
      xratioDown = xDown/xCentral;
    }
    ratioUp->SetBinContent(iB,xratioUp);
    ratioDown->SetBinContent(iB,xratioDown);
    ratioUp->SetBinError(iB,0);
    ratioDown->SetBinError(iB,0);
    ratioCentral->SetBinContent(iB,1);
    ratioCentral->SetBinError(iB,0);
    if (histNominal->GetBinContent(iB)>0)
      ratioCentral->SetBinError(iB,histNominal->GetBinError(iB)/histNominal->GetBinContent(iB));
  }

  histUp->GetYaxis()->SetTitleOffset(1.4);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
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

  histUp->Draw("h");
  histNominal->Draw("pesame");
  histDown->Draw("hsame");
  TLegend * leg = new TLegend(0.55,0.65,0.92,0.9);
  leg->SetHeader(header);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(histNominal,"Central","l");
  leg->AddEntry(histUp,"Up","l");
  leg->AddEntry(histDown,"Down","l");
  leg->Draw();
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  ratioUp->GetYaxis()->SetRangeUser(0.9,1.1);
  ratioUp->GetYaxis()->SetNdivisions(505);
  ratioUp->GetXaxis()->SetLabelFont(42);
  ratioUp->GetXaxis()->SetLabelOffset(0.04);
  ratioUp->GetXaxis()->SetLabelSize(0.1);
  ratioUp->GetXaxis()->SetTitleSize(0.13);
  ratioUp->GetXaxis()->SetTitleOffset(1.2);
  ratioUp->GetYaxis()->SetTitle("ratio");
  ratioUp->GetYaxis()->SetLabelFont(42);
  ratioUp->GetYaxis()->SetLabelOffset(0.015);
  ratioUp->GetYaxis()->SetLabelSize(0.1);
  ratioUp->GetYaxis()->SetTitleSize(0.14);
  ratioUp->GetYaxis()->SetTitleOffset(0.5);
  ratioUp->GetXaxis()->SetTickLength(0.07);
  ratioUp->GetYaxis()->SetTickLength(0.04);
  ratioUp->GetYaxis()->SetLabelOffset(0.01);

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

  ratioUp->Draw("h");
  ratioDown->Draw("hsame");
  ratioCentral->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();


  canv1->Print("figures/"+varName+"_"+sysName+".png");

} 
