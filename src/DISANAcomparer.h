#ifndef DISANA_COMPARER_H
#define DISANA_COMPARER_H

// ROOT headers
#include <TCanvas.h>
#include <TLegend.h>

#include <sys/stat.h>
// STL headers
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project-specific headers
#include <chrono>

#include "DISANAplotter.h"
#include "DrawStyle.h"

namespace fs = std::filesystem;
// color palettes for different models
std::vector<std::tuple<double, double, double>> modelShades = {
    {0.20, 0.30, 0.85},  // Blue
    {0.90, 0.45, 0.10},  // Orange
    {0.00, 0.60, 0.60},  // Teal green
    {0.00, 0.70, 0.00},  // Green
    {0.60, 0.30, 0.80},  // Purple
    {0.85, 0.10, 0.25},  // Red
    {0.40, 0.40, 0.40}   // Gray (fallback)
};

class DISANAcomparer {
 public:
  // Set the bin ranges used for cross-section calculations and plotting
  void SetXBinsRanges(BinManager bins) { fXbins = bins; }

  void NormalizeHistogram(TH1* hist) {
    if (!hist) return;
    double integral = hist->Integral();
    if (integral > 0) hist->Scale(1.0 / integral);
  }
  void AddModelPhi(ROOT::RDF::RNode df, const std::string& label, double beamEnergy) {
    auto plotter = std::make_unique<DISANAplotter>(df, beamEnergy);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV without Pi0 Correction" << std::endl;
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("kMinus");
    plotter->GenerateKinematicHistos("kPlus");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  // Set the output directory for saving plots
  void SetOutputDir(const std::string& outdir) {
    outputDir = outdir;
    if (!fs::exists(outputDir)) {
      fs::create_directories(outputDir);
    }
  }

  // Enable or disable individual variable plotting
  void PlotIndividual(bool plotInd) { plotIndividual = plotInd; }

  // Set plot styles for various plot types
  void SetKinStyle(const DrawStyle& style) { styleKin_ = style; }
  void SetDVCSStyle(const DrawStyle& style) { styleDVCS_ = style; }
  void SetCrossSectionStyle(const DrawStyle& style) { styleCrossSection_ = style; }
  void SetBSAStyle(const DrawStyle& style) { styleBSA_ = style; }
  void UseFittedPhiYields(bool on = true) { useFittedYields_ = on; }
  /// get mean values of Q^2 and x_B
  std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> getMeanQ2xBt(const BinManager& bins, std::unique_ptr<DISANAplotter>& plotter) {
    const auto& xb_bins = bins.GetXBBins();
    const auto& q2_bins = bins.GetQ2Bins();
    const auto& t_bins = bins.GetTBins();

    size_t n_xb = xb_bins.size() - 1;
    size_t n_q2 = q2_bins.size() - 1;
    size_t n_t = t_bins.size() - 1;

    auto rdf = plotter->GetRDF();

    std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
        n_xb, std::vector<std::vector<std::tuple<double, double, double>>>(n_q2, std::vector<std::tuple<double, double, double>>(n_t)));

    for (size_t ix = 0; ix < n_xb; ++ix) {
      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t it = 0; it < n_t; ++it) {
          double xb_lo = xb_bins[ix], xb_hi = xb_bins[ix + 1];
          double q2_lo = q2_bins[iq], q2_hi = q2_bins[iq + 1];
          double t_lo = t_bins[it], t_hi = t_bins[it + 1];

          // Apply filter
          auto rdf_cut = rdf.Filter(Form("xB >= %f && xB < %f", xb_lo, xb_hi)).Filter(Form("Q2 >= %f && Q2 < %f", q2_lo, q2_hi)).Filter(Form("t >= %f && t < %f", t_lo, t_hi));

          // Compute means
          double mean_xB = rdf_cut.Mean("xB").GetValue();
          double mean_Q2 = rdf_cut.Mean("Q2").GetValue();
          double mean_t = rdf_cut.Mean("t").GetValue();

          result[ix][iq][it] = std::make_tuple(mean_xB, mean_Q2, mean_t);
        }
      }
    }

    return result;
  }

  // Plot all basic kinematic distributions (p, theta, phi) for all particle types
  void PlotKinematicComparison_phiAna() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 4);

    std::vector<std::string> types = {"el", "pro", "kPlus", "kMinus"};
    std::vector<std::string> vars = {"p", "theta", "phi"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotVariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "KinematicComparison_phiAna.pdf").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison_phiAna.pdf" << std::endl;
    delete canvas;
  }

  // Plot one specific variable (e.g., p) for a given particle type (e.g., "el")
  void PlotVariableComparison(const std::string& type, const std::string& var, int padIndex, TCanvas* canvas) {
    canvas->cd(padIndex);
    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;
    styleKin_.StylePad((TPad*)gPad);

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }
      NormalizeHistogram(target);
      styleKin_.StyleTH1(target);
      auto [cr, cg, cb] = modelShades[i % modelShades.size()];
      const int colorIdx = 4000 + int(i) * 20;
      if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
      target->SetMarkerColor(colorIdx);
      target->SetLineColorAlpha(colorIdx, 0.8);
      target->SetLineWidth(1);
      target->SetTitle(Form("%s;%s;Count", typeToParticle[type].c_str(), VarName[var].c_str()));

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
  }

  // Save an individual variable comparison plot as PNG
  void PlotSingleVariableComparison(const std::string& type, const std::string& var) {
    TCanvas* canvas = new TCanvas(("c_" + type + "_" + var).c_str(), ("Comparison " + type + " " + var).c_str(), 800, 600);
    gPad->SetGrid();

    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotSingleVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }

      target->SetLineColor(i + 1);

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
    canvas->Update();
    canvas->SaveAs((outputDir + "/compare_" + type + "_" + var + ".pdf").c_str());
    delete canvas;
  }

  void PlotPhiElectroProKinematicsComparison(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1400);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.04);

      bool first = true;
      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto rdf = plotters[i]->GetRDF();
        if (!rdf.HasColumn(var)) {
          std::cerr << "[ERROR] Column " << var << " not found in RDF for model " << labels[i] << "\n";
          continue;
        }

        double min = *(rdf.Min(var));
        double max = *(rdf.Max(var));
        if (min == max) {
          min -= 0.1;
          max += 0.1;
        }
        double margin = std::max(1e-3, 0.05 * (max - min));

        // Get histogram (RResultPtr) and clone it
        auto htmp = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);
        auto h = (TH1D*)htmp->Clone(Form("h_%s_%zu_clone", var.c_str(), i));

        if (!h) continue;  // guard against failed clone

        h->SetDirectory(0);  // prevent ROOT from managing ownership
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);
        auto [cr, cg, cb] = modelShades[i % modelShades.size()];
        const int colorIdx = 4000 + int(i) * 20;
        if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
        h->SetMarkerColor(colorIdx);
        h->SetLineColorAlpha(colorIdx, 0.8);
        h->SetLineWidth(1.0);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }

      if (!histos_to_draw.empty()) {
        legend->Draw();
      }

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi")) {
        PlotSingleVariableComparison("el", var);
      }
    }
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 60, 0, 8.0, 60, 0, 10.0}, "t", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/PhiAna_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVCS kinematics comparison to: " << outputDir + "/PhiAna_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  

  void PlotxBQ2tBin(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, 0, 1.0, 500, 0, 10.0}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, 0, 1.0, 500, 0, 10.0}, "t", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, 0, 1.0, 500, 0, 1.0}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBin.pdf").c_str());
    std::cout << "Saved xBQ2tBin kinematics to: " << outputDir + "/xBQ2tBin.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  // phi analysis
  /// For exclusivity cuts, you can use the following function to select one triplet
  void PlotPhiAnaExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", 0.0, 1.3},
        {"Mx2_epKpKm", "Missing Mass Squared (epK^{+}K^{-})", "MM^{2}(epK^{+}K^{-}) [GeV^{2}]", -0.08, 0.08},
        {"Mx2_eKpKm", "Invariant Mass (eK^{+}K^{-})", "M^{2}(eK^{+}K^{-}) [GeV^{2}]", -0.5, 3},
        {"Mx2_epKp", "Missing Mass Squared (epK^{+})", "MM^{2}(epK^{+}) [GeV^{2}]", -0.5, 1.5},
        {"Mx2_epKm", "Missing Mass Squared (epK^{-})", "MM^{2}(epK^{-}) [GeV^{2}]", -0.5, 1.5},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.5},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
        {"Theta_e_phimeson", "Angle: e-#phi", "#theta(e, #phi) [deg]", 0.0, 60.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 3;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cut = plotters[m]->GetRDF().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColorAlpha(m + 4,0.8);
          auto [cr, cg, cb] = modelShades[m % modelShades.size()];
          const int colorIdx = 4000 + int(m) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h_clone->SetMarkerColor(colorIdx);
          h_clone->SetLineColorAlpha(colorIdx, 0.8);
          h_clone->SetLineWidth(1);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColorAlpha(colorIdx, 0.8);
          line2->SetLineColorAlpha(colorIdx, 0.8);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_Phi_Ana" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  bool file_exists(const char* name) {
    struct stat buffer;
    return (stat(name, &buffer) == 0);
  }

  void dumpHistogram(TH1D* h,
                   double xB,
                   double Q2,
                   double t,
                   double xBmin,
                   double xBmax,
                   double Q2min,
                   double Q2max,
                   double tmin,
                   double tmax,
                   const char* filename="h_data.txt") {
    bool exists = file_exists(filename);

    std::ofstream fout(filename, std::ios::out | std::ios::app);
    if (!fout.is_open()) {
        std::cerr << "cannot open " << filename << " to write!\n";
        return;
    }

    if (!exists) {
        fout << "# xB\tQ2\t-t\tphi\tvalue\terror\txBmin\txBmax\tQ2min\tQ2max\ttmin\ttmax\n";
    }
    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        double phi   = h->GetBinCenter(ibin);
        double value = h->GetBinContent(ibin);
        double err   = h->GetBinError(ibin);
        fout << xB   << "\t"
             << Q2   << "\t"
             << t    << "\t"
             << phi  << "\t"
             << value<< "\t"
             << err  << "\t"
             << xBmin<< "\t"
             << xBmax<< "\t"
             << Q2min<< "\t"
             << Q2max<< "\t"
             << tmin << "\t"
             << tmax << "\n";
    }

    fout.close();
    std::cout << "Data " << h->GetName()
              << " written into " << filename
              << (exists ? " (appended)" : "") << "\n";
  }

  
  void PlotPhiDSigmaDt_FromCache(bool logy = true) {
    if (plotters.empty()) return;

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& t  = fXbins.GetTBins();
    const auto& w  = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_dsdt_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);
        if (logy) gPad->SetLogy();

        // legend in the same spirit as DVCS cross-section plots
        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        // find a sensible common Y-range across models for this (Q2,W) slice
        double yMinPos = std::numeric_limits<double>::infinity();
        double yMaxVal = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& xs3D = plotters[im]->GetPhiDSigmaDt3D();
          if (iq >= xs3D.size() || iw >= xs3D[iq].size()) continue;
          TH1D* h = xs3D[iq][iw];
          if (!h) continue;
          const double I = h->Integral(1, h->GetNbinsX());          // sum of contents
          //if (I > 0)  h->Scale(1.0 / I);
          for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double v = h->GetBinContent(b);
            if (v > 0.0 && v < yMinPos) yMinPos = v;
            if (v > yMaxVal) yMaxVal = v;
          }
        }
        if (!std::isfinite(yMinPos)) yMinPos = 1e-4;
        if (yMaxVal <= 0.0) yMaxVal = 1.0;
        if (logy) { yMinPos *= 0.5; yMaxVal *= 3.0; }

        // header text (bin labels)
        const TString head = hasW
            ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq+1], w[iw], w[iw+1])
            : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq+1]);

        bool first = true;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& xs3D = plotters[im]->GetPhiDSigmaDt3D();
          if (iq >= xs3D.size() || iw >= xs3D[iq].size()) continue;
          TH1D* h = xs3D[iq][iw];
          if (!h) continue;

          // apply cross-section style + consistent palette
          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);

          // axis cosmetics consistent with DVCS cross-section look
          h->SetTitle("");
          h->GetXaxis()->SetTitle("-t [GeV^{2}]");
          h->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          if (logy) h->GetYaxis()->SetRangeUser(yMinPos, yMaxVal);

          if (first) {
            h->Draw("E1X0");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
      
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        // save inside the configured outputDir
        TString out = hasW
            ? Form("%s/phi_dsdt_Q%zu_W%zu.pdf", outputDir.c_str(), iq, iw)
            : Form("%s/phi_dsdt_Q%zu.pdf", outputDir.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }
  // === Add to DISANAcomparer (public): =========================
  void PlotPhiInvMassPerBin_AllModels(const std::string& baseOutDir = "PhiInvMassFits", int nBins = 120, double mMin = 0.98, double mMax = 1.08, bool constrainSigma = true,
                                      double sigmaRef = 0.004, double sigmaFrac = 0.25, double luminosity_rga_fall18 = 1.0, double branching = 1.0, bool IsMC = false) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiInvMassPerBin] no models.\n";
      return;
    }
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing per-bin K^{+}K^{-} mass for model: " << labels[i] << " → " << subdir << std::endl;
      (void)plotters[i]->MakePhiMassFitCanvases3D(fXbins, subdir, nBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, luminosity_rga_fall18,branching,IsMC);
    }
  }


 private:
  BinManager fXbins;
  bool plotIndividual = false;
  bool useFittedYields_ = true;

  DrawStyle style_;              // Default style
  DrawStyle styleKin_;           // Kin plot style
  DrawStyle styleDVCS_;          // DVCS plot style
  DrawStyle styleCrossSection_;  // Cross-section plot style
  DrawStyle styleBSA_;           // BSA plot style

  std::unique_ptr<ROOT::RDF::RNode> rdf;
  std::string outputDir = ".";

  std::vector<std::unique_ptr<DISANAplotter>> plotters;
  std::vector<std::string> labels;

  std::vector<std::string> particleName = {"e", "p", "#gamma"};
  std::map<std::string, std::string> typeToParticle = {{"el", "electron"}, {"pro", "proton"}, {"pho", "#gamma"}, {"kMinus", "K^{-}"}, {"kPlus", "K^{+}"}};
  std::map<std::string, std::string> VarName = {{"p", "p (GeV/#it{c})"}, {"theta", "#theta (rad)"}, {"phi", "#phi(rad)"},{"vz", "v_{z}(cm)"}};
};
#endif  // DISANA_COMPARER_H
