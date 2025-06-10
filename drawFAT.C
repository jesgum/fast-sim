
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

// Usage: in Delphes O2 run `root DetectorK/HistoManager.cxx+ DetectorK/DetectorK.cxx+ drawFAT.C`
// Requires AliPhysics and to be ran here: https://github.com/AliceO2Group/DelphesO2/tree/master/src

void drawFAT()
{
  DetectorK fat = DetectorK("example", "example");
  Double_t x0IB = 0.001;
  Double_t x0OB = 0.01;
  Double_t xrhoIB = 2.3292e-02; // 100 mum Si
  Double_t xrhoOB = 2.3292e-01; // 1000 mum Si
  Double_t resRPhiIB = 0.00025;
  Double_t resZIB = 0.00025;
  Double_t resRPhiOB = 0.00100;
  Double_t resZOB = 0.00100;
  Double_t eff = 0.98;
  fat.AddLayer("vertex", 0.0, 0, 0);                // dummy vertex for matrix calculation
  fat.AddLayer("bpipe0", 0.48, 0.00042, 2.772e-02); // 150 mum Be
  fat.AddLayer("B00", 0.50, x0IB, xrhoIB, resRPhiIB, resZIB, eff);
  fat.AddLayer("B01", 1.20, x0IB, xrhoIB, resRPhiIB, resZIB, eff);
  fat.AddLayer("B02", 2.50, x0IB, xrhoIB, resRPhiIB, resZIB, eff);
  fat.AddLayer("bpipe1", 3.7, 0.0014, 9.24e-02); // 500 mum Be
  fat.AddLayer("B03", 3.75, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B04", 7.00, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B05", 12.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B06", 20.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B07", 30.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B08", 45.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B09", 60.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B10", 80.0, x0OB, xrhoOB, resRPhiOB, resZOB, eff);
  fat.AddLayer("B11", 100., x0OB, xrhoOB, resRPhiOB, resZOB, eff);

  struct binning {
    const float min;
    const float max;
    const int nBins;
    float x(int i) const
    {
      const float w = (width()) / 2;
      return min + (i + 1) * w;
    }
    float width() const
    {
      return (max - min) / nBins;
    }
  };
  const float mass = 0.13957; // pion mass in GeV/c^2
  const float q = 1.0;        // charge in e
  const int nch = 2200;
  fat.SetdNdEtaCent(10); // set the dN/deta for the FAT
  binning ptBinning = {0., 10, 1000};
  const float etaVsPt = 0.0; // pseudorapidity vs pT
  binning etaBinning = {0., 10, 1000};
  const float ptVsEta = 0.1; // pseudorapidity vs eta

  TGraph* gPt = new TGraph();
  gPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gPt->GetYaxis()->SetTitle("Efficiency (%)");
  gPt->SetLineColor(kBlue + 1);
  gPt->SetLineWidth(2);
  gPt->SetMarkerStyle(20);
  gPt->SetMarkerSize(0.7);
  gPt->SetMarkerColor(kBlue + 1);

  TGraph* gPtReso = new TGraph();
  gPtReso->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gPtReso->GetYaxis()->SetTitle("#sigma_{p_{T}} / p_{T} (%)");
  gPtReso->SetLineColor(kRed + 1);
  gPtReso->SetLineWidth(2);
  gPtReso->SetMarkerStyle(21);
  gPtReso->SetMarkerSize(0.7);
  gPtReso->SetMarkerColor(kRed + 1);

  TGraph* gEtaEfficiency = new TGraph();
  gEtaEfficiency->GetXaxis()->SetTitle("#eta");
  gEtaEfficiency->GetYaxis()->SetTitle("Efficiency (%)");
  gEtaEfficiency->SetLineColor(kGreen + 2);
  gEtaEfficiency->SetLineWidth(2);
  gEtaEfficiency->SetMarkerStyle(22);
  gEtaEfficiency->SetMarkerSize(0.7);
  gEtaEfficiency->SetMarkerColor(kGreen + 2);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetLabelSize(0.03, "XY");
  gStyle->SetTitleSize(0.04, "XY");

  // Now we compute the efficiency for a range of pt values and plot the results
  for (int i = 0; i < ptBinning.nBins; ++i) {
    float pt = ptBinning.x(i);
    // solve track
    TrackSol tr(1, pt, etaVsPt, q, mass);
    if (!fat.SolveTrack(tr))
      continue;
    AliExternalTrackParam* trPtr = (AliExternalTrackParam*)tr.fTrackCmb.At(0);
    if (!trPtr)
      continue;
    std::array<float, 15> covm;
    for (int j = 0; j < 15; ++j)
      covm[j] = trPtr->GetCovariance()[j];

    float sigma1Pt2 = covm[14];                // sigma^2 1/pT
    float mom_res = std::sqrt(sigma1Pt2) * pt; // momentum resolution (%)
    gPtReso->AddPoint(pt, mom_res * 100.);     // pt resolution (%)
    // define the efficiency
    float eff = 1;
    for (int j = 1; j < 20; ++j) {
      auto igoodhit = fat.GetGoodHitProb(j);
      if (igoodhit <= 0.)
        continue;
      eff *= igoodhit;
    }
    gPt->AddPoint(pt, eff);
    Printf("pt = %f, eff = %f", pt, eff);
  }

  // Now we compute the efficiency for a range of eta values and plot the results
  for (int i = 0; i < etaBinning.nBins; ++i) {
    float eta = etaBinning.x(i);
    // solve track
    TrackSol tr(1, ptVsEta, eta, q, mass);
    if (!fat.SolveTrack(tr))
      continue;
    AliExternalTrackParam* trPtr = (AliExternalTrackParam*)tr.fTrackCmb.At(0);
    if (!trPtr)
      continue;
    float eff = 1;
    for (int j = 1; j < 20; ++j) {
      auto igoodhit = fat.GetGoodHitProb(j);
      if (igoodhit <= 0.)
        continue;
      eff *= igoodhit;
    }
    gEtaEfficiency->AddPoint(eta, eff * 100.); // efficiency (%)
    Printf("eta = %f, eff = %f", eta, eff);
  }

  // Now we plot the results
  TCanvas* canvPt = new TCanvas("canvPt", "gPt", 800, 600);
  gPt->Draw("ALP");
  canvPt->SaveAs("canvPt.pdf");

  TCanvas* canvPtReso = new TCanvas("canvPtReso", "gPtReso", 800, 600);
  gPtReso->Draw("ALP");
  canvPtReso->SaveAs("canvPtReso.pdf");

  TCanvas* canvEtaEfficiency = new TCanvas("canvEtaEfficiency", "gEtaEfficiency", 800, 600);
  gEtaEfficiency->Draw("ALP");
  canvEtaEfficiency->SaveAs("canvEtaEfficiency.pdf");
}
