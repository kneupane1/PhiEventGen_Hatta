#ifndef DISANA_PLOTTER_H
#define DISANA_PLOTTER_H

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include <string>
#include <vector>

#include "DISANAMath.h"

class BinManager;

class DISANAplotter {
public:
  DISANAplotter(ROOT::RDF::RNode df_dvcs_data, double beamEnergy)
      : rdf(std::move(df_dvcs_data)), beam_energy(beamEnergy) {}

  ~DISANAplotter() {
    disHistos.clear();
    kinematicHistos.clear();
    // delete cached cross-section histos
    for (auto &row : phi_dsdt_QW_) {
      for (auto *h : row) {
        delete h;
      }
    }
    phi_dsdt_QW_.clear();
  }
  using AccEffProvider =
      std::function<double(double q2lo, double q2hi, double wlo, double whi,
                           double tlo, double thi)>;
  void SetAccEffProvider(AccEffProvider f) { accEffProvider_ = std::move(f); }
  ROOT::RDF::RNode GetRDF() { return rdf; }
  // ---------------- NEW: phi mass-fit results & acceptance provider
  // ---------------
  struct PhiMassFitResult {
    double mu{};    // peak position
    double sigma{}; // gaussian width
    double mLo{}, mHi{};
    ULong64_t Nwin{}; // raw counts in ±3σ
    double Nbkg{};    // fitted background under window
    double Nsig{};    // signal = Nwin - Nbkg
    double dNsig{};   // simple stat ⊕ bkg-fit uncertainty
  };

  struct PhiMassBin {
    TH1D *hM = nullptr;    // invariant-mass histogram (owned, SetDirectory(0))
    TF1 *fTotal = nullptr; // total fit
    TF1 *fSig = nullptr;   // Gaussian signal
    TF1 *fBkg = nullptr;   // background
    double mu = 0, sigma = 0, mLo = 0, mHi = 0;
    double Nwin = 0, Nbkg = 0, Nsig = 0, dNsig = 0;
  };
  // --- in public: add a const accessor ---
  const std::vector<std::vector<TH1D *>> &GetPhiDSigmaDt3D() const {
    return phi_dsdt_QW_;
  }
  std::vector<std::vector<TH1D *>> &GetPhiDSigmaDt3D() {
    return phi_dsdt_QW_;
  } // (optional mutable)

  void GenerateKinematicHistos(const std::string &type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi"}; // for Phi analysis
    // std::vector<std::string> vars = {"p", "theta", "phi"}; // for DVCS
    // analysis
    for (const auto &v : vars) {
      std::string base = "rec" + type + "_" + v;
      // std::cout << " base is " << base << std::endl;
      double histMin, histMax;
      auto it = kinematicAxisRanges.find(base);
      if (it != kinematicAxisRanges.end()) {
        histMin = it->second.first;
        histMax = it->second.second;
      } else {
        auto minVal = rdf.Min(base);
        auto maxVal = rdf.Max(base);
        double lo = *minVal;
        double hi = *maxVal;
        double margin = std::max(1e-3, 0.05 * (hi - lo));
        histMin = lo - margin;
        histMax = hi + margin;
        kinematicAxisRanges[base] = {histMin, histMax}; // store it
      }

      // auto h_all = rdf.Histo1D({(base + "_all").c_str(), "", 100, histMin,
      // histMax}, base);
      auto h_acc =
          rdf.Histo1D({(base).c_str(), "", 100, histMin, histMax}, base);
      kinematicHistos.push_back(h_acc);
    }

    std::map<std::string, std::pair<double, double>> cachedRanges;

    for (const auto &var : disvars) {
      double histMin = axisRanges[var].first;
      double histMax = axisRanges[var].second;
      // auto h_all = rdf.Histo1D({(var + "_all").c_str(), (var + "
      // All").c_str(), 100, histMin, histMax}, var);
      auto h_acc =
          rdf.Histo1D({var.c_str(), var.c_str(), 100, histMin, histMax}, var);
      disHistos.push_back(h_acc);
    }
  }

  std::vector<TH1 *> GetDISHistograms() {
    std::vector<TH1 *> allDIShisto;
    for (auto &h : disHistos)
      allDIShisto.push_back(h.GetPtr());
    return allDIShisto;
  }

  std::vector<TH1 *> GetAllHistograms() {
    std::vector<TH1 *> all;
    for (auto &h : kinematicHistos)
      all.push_back(h.GetPtr());
    for (auto &h : disHistos)
      all.push_back(h.GetPtr());
    return all;
  }

  std::vector<TH1D *> ComputePhiCrossSection(const BinManager &bins,
                                             double luminosity) {
    return kinCalc.ComputePhi_CrossSection(rdf, bins, luminosity);
  }

  // === in public: ==============================================
  struct PhiMassDraw {
    TH1D *h{nullptr};
    TF1 *fTot{nullptr};
    TF1 *fSig{nullptr};
    TF1 *fBkg{nullptr};
    double mu{0}, sigma{0}, mLo{0}, mHi{0};
    TCanvas *c{nullptr};
    std::string name;
  };

  inline std::vector<std::vector<std::vector<PhiMassDraw>>>
  MakePhiMassFitCanvases3D(
      const BinManager &bins, const std::string &outDirPerModel,
      int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
      bool constrainSigma = true, double sigmaRef = 0.004,
      double sigmaFrac = 0.30,
      // --- NEW (optional): compute dσ/dt when luminosity > 0 ---
      double luminosity_nb_inv = -1.0, double branching = 0.492,
      bool IsMC = false) {
    if (luminosity_nb_inv <= 0.0) {
      std::cerr << "[MakePhiMassFitCanvases3D] luminosity<=0 → will NOT "
                   "compute dσ/dt cache.\n";
    }
    // ensure directory exists
    gSystem->mkdir(outDirPerModel.c_str(), /*recursive=*/true);

    const auto &q2Edges = bins.GetQ2Bins();
    const auto &tEdges = bins.GetTBins();
    const auto &wEdges = bins.GetWBins();
    const bool hasW = !wEdges.empty();

    const size_t nQ = q2Edges.size() > 1 ? q2Edges.size() - 1 : 0;
    const size_t nW = hasW ? (wEdges.size() - 1) : 1;
    const size_t nT = tEdges.size() > 1 ? tEdges.size() - 1 : 0;

    std::vector<std::vector<std::vector<PhiMassDraw>>> out;
    if (!nQ || !nT)
      return out;
    out.resize(nQ);

    // ---- NEW: prepare the 3D TH1 cache for dσ/dt (Q,W -> TH1 over t) ----
    phi_dsdt_QW_.assign(nQ, std::vector<TH1D *>(nW, nullptr));
    static std::atomic<unsigned long> uid_dsdt{0};

    for (size_t iq = 0; iq < nQ; ++iq) {
      out[iq].resize(nW);

      const double qLo = q2Edges[iq], qHi = q2Edges[iq + 1];
      auto df_q =
          rdf.Filter([=](double Q2) { return Q2 > qLo && Q2 <= qHi; }, {"Q2"});

      for (size_t iw = 0; iw < nW; ++iw) {
        const double wLo = hasW ? wEdges[iw] : 0.0;
        const double wHi = hasW ? wEdges[iw + 1] : 1e9;
        auto df_qw = df_q.Filter(
            [=](double W) { return (!hasW) || (W > wLo && W <= wHi); }, {"W"});

        out[iq][iw].resize(nT);

        // --- NEW: create the dσ/dt histogram for this (Q2,W) slice (if we will
        // compute) ---
        if (luminosity_nb_inv > 0.0) {
          auto hname =
              Form("phi_dsdt_Q%zu_W%zu_%lu", iq, iw, uid_dsdt.fetch_add(1));
          auto htitle = Form("d#sigma/dt; -t [GeV^{2}]; nb/GeV^{2}");
          TH1D *hDsDt =
              new TH1D(hname, htitle, static_cast<int>(tEdges.size() - 1),
                       tEdges.data());
          hDsDt->SetDirectory(nullptr);
          phi_dsdt_QW_[iq][iw] = hDsDt;
        }

        for (size_t it = 0; it < nT; ++it) {
          const double tLo = tEdges[it], tHi = tEdges[it + 1];
          const double dT = tHi - tLo;

          auto df_bin = df_qw.Filter(
              [=](double t) { return t > tLo && t <= tHi; }, {"t"});

          static std::atomic<unsigned long> uid{0};
          auto hname = Form("hM_%zu_%zu_%zu_%lu", iq, iw, it, uid.fetch_add(1));
          auto hR = df_bin.Histo1D(
              ROOT::RDF::TH1DModel(hname, ";M_{K^{+}K^{-}} [GeV];Counts", 200,
                                   0.8, 1.8),
              "invMass_KpKm");
          hR.GetValue();
          TH1D *h = (TH1D *)hR.GetPtr();
          if (!h || h->GetEntries() <= 20) {
            // still record an empty point in dσ/dt if we're computing it
            if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
              const int b = static_cast<int>(it + 1);
              phi_dsdt_QW_[iq][iw]->SetBinContent(b, 0.0);
              phi_dsdt_QW_[iq][iw]->SetBinError(b, 0.0);
            }
            continue;
          }

          // ---- (existing) pretty clone for drawing ----
          TH1D *hDraw = (TH1D *)h->Clone(Form("%s_draw", h->GetName()));
          hDraw->SetDirectory(0);
          hDraw->SetTitle("");
          hDraw->SetLineColor(kBlue + 1);
          hDraw->SetLineWidth(2);
          hDraw->GetYaxis()->SetTitleOffset(1.2);
          hDraw->SetFillColorAlpha(kBlue - 9, 0.3);
          hDraw->SetMarkerStyle(20);
          hDraw->SetMarkerSize(1.2);
          hDraw->SetMarkerColor(kBlue + 2);
          hDraw->GetXaxis()->SetTitle("M(K^{+}K^{-}) [GeV]");
          hDraw->GetYaxis()->SetTitle("Counts");
          hDraw->GetXaxis()->SetRangeUser(mMin - 0.02, mMax + .020);

          const double bw = hDraw->GetXaxis()->GetBinWidth(1);
          std::cout << Form("[MakePhiMassFitCanvases3D] Fitting %s with %.0f "
                            "bins, M=[%.4f,%.4f], BW=%.4f GeV",
                            hDraw->GetName(), (double)hDraw->GetNbinsX(), mMin,
                            mMax, bw)
                    << std::endl;

          std::string formula = Form(
              // Background: [0]*(x-0.9874)^[1] * exp(-[2]*(x-0.9874))
              "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))"
              " + "
              // Signal (counts/bin): [3] (=Nsig) * bw * Gaussian_pdf(x;
              // [4]=mean, [5]=sigma)
              "[3]*0.398942*%g*TMath::Exp(-0.5*TMath::Power((x-[4])/([5]),2))/"
              "TMath::Abs([5])",
              bw);

          TF1 *fTot =
              new TF1(Form("fitTotal_%s", hname), formula.c_str(), mMin, mMax);

          //[0]*0.398942*bw*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])
          // params: [0]=Abkg(scale), [1]=alpha, [2]=lambda, [3]=Nsig(yield),
          // [4]=mu, [5]=sigma
          fTot->SetParNames("A", "alpha", "lambda", "N", "mu", "sigma");
          fTot->SetParameters(4, 0.9, 2, 10, 1.02, 0.010);
          fTot->SetParLimits(4, 1.005, 1.022);
          fTot->SetParLimits(5, 0.0025, 0.025);
          fTot->SetParLimits(3, 0.000001, 100000.0);
          fTot->SetLineColor(kRed + 1);
          fTot->SetLineWidth(3);
          if (!IsMC)
            hDraw->Fit(fTot, "R0QL");
          if (!IsMC)
            fTot->Draw("SAME C");

          const double A = fTot->GetParameter(0);
          const double alpha = fTot->GetParameter(1);
          const double lambda = fTot->GetParameter(2);
          const double N = fTot->GetParameter(3);
          const double mu = fTot->GetParameter(4);
          const double sigma = fTot->GetParameter(5);
          const double chi2 = fTot->GetChisquare();
          const double ndf = fTot->GetNDF();

          const double mLo = mu - 3.0 * sigma;
          const double mHi = mu + 3.0 * sigma;

          TF1 *fSig = new TF1(Form("fSig_%s", hname),
                              Form("[0]*0.398942*%g*TMath::Exp(-0.5*TMath::"
                                   "Power((x-[1])/([2]),2))/TMath::Abs([2])",
                                   bw),
                              mMin, mMax);
          fSig->SetParameters(N, mu, sigma);
          fSig->SetLineColor(kOrange + 1);
          fSig->SetLineStyle(2);
          fSig->SetLineWidth(2);
          fSig->SetFillColorAlpha(kOrange - 3, 0.3);
          fSig->SetFillStyle(1001);
          if (!IsMC)
            fSig->Draw("SAME FC");

          TF1 *fBkg = new TF1(
              Form("fBkg_%s", hname),
              "[0]*TMath::Power(x-0.9874,[1])*TMath::Exp(-[2]*(x-0.9874))",
              mMin, mMax);
          fBkg->SetParameters(A, alpha, lambda);
          fBkg->SetLineColor(kGreen + 2);
          fBkg->SetLineStyle(3);
          fBkg->SetLineWidth(2);
          fBkg->SetFillColorAlpha(kGreen - 7, 0.3);
          fBkg->SetFillStyle(1001);
          if (!IsMC)
            fBkg->Draw("SAME FC");
          // background under ±3σ (bin-center sum; consistent with fit mode)
          double Nbkg = 0.0;
          double Nsig = 0.0;     // fTot->GetParameter(3);
          double Nsig_err = 0.0; // fTot->GetParError(3);
          if (!IsMC) {

            // 2) Define ±3σ window, clamp to [mMin,mMax]
            double mLo_ = std::max(mMin, 1.0195 - 5.0 * 0.004);
            double mHi_ = std::min(mMax, 1.0195 + 5.0 * 0.004);

            // 3) Integrate counts in [mLo,mHi] with errors from sumw2
            int iLo = hDraw->GetXaxis()->FindBin(mLo_);
            int iHi = hDraw->GetXaxis()->FindBin(mHi_);
            double errInt = 0.0;
            Nsig = hDraw->IntegralAndError(iLo, iHi, errInt);
            Nsig_err = errInt;

          } else {
            Nsig = fTot->GetParameter(3);
            Nsig_err = fTot->GetParError(3);
          }
          // window counts & signal
          auto dfWin = df_bin.Filter(
              [=](float m) { return m > mLo && m < mHi; }, {"invMass_KpKm"});
          const ULong64_t Nwin = *dfWin.Count();

          // ---- NEW: fill dσ/dt if requested ----
          if (luminosity_nb_inv > 0.0 && phi_dsdt_QW_[iq][iw]) {
            const double Aeps = accEffProvider_(
                qLo, qHi, wLo, wHi, tLo, tHi); // user-supplied or default 1
            const int b = static_cast<int>(it + 1);
            double val = 0.0, err = 0.0;
            if (Aeps > 0.0 && branching > 0.0 && dT > 0.0) {
              val = Nsig / (luminosity_nb_inv * Aeps * branching * dT);
              err = Nsig_err / (luminosity_nb_inv * Aeps * branching * dT);
            }
            phi_dsdt_QW_[iq][iw]->SetBinContent(b, val);
            phi_dsdt_QW_[iq][iw]->SetBinError(b, err);
          }

          // ---- (existing) canvas beautification & save ----
          TCanvas *c =
              new TCanvas(Form("c_%s", hname), "K^{+}K^{-} mass", 1200, 1000);
          gStyle->Reset();
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(0);
          gStyle->SetTitleFont(42, "XYZ");
          gStyle->SetLabelFont(42, "XYZ");
          gStyle->SetTitleSize(0.05, "XYZ");
          gStyle->SetLabelSize(0.04, "XYZ");
          c->SetTicks(1, 1);
          c->SetTopMargin(0.04);
          c->SetRightMargin(0.03);
          c->SetBottomMargin(0.11);
          c->SetLeftMargin(0.13);
          gStyle->SetCanvasColor(0);
          gStyle->SetPadColor(0);
          c->SetFillColor(0);
          if (c->GetPad(0))
            c->GetPad(0)->SetFillColor(0);
          c->cd();
          hDraw->SetMinimum(0.0);
          hDraw->Draw("PE");
          fTot->SetLineColor(kRed + 1);
          fBkg->Draw("SAME FC");
          fSig->Draw("SAME FC");
          fTot->Draw("SAME C");

          double yMax = hDraw->GetMaximum();
          TLine *L1 = new TLine(mLo, 0.0, mLo, yMax * 0.75);
          TLine *L2 = new TLine(mHi, 0.0, mHi, yMax * 0.75);
          L1->SetLineColor(kMagenta + 2);
          L2->SetLineColor(kMagenta + 2);
          L1->SetLineStyle(2);
          L2->SetLineStyle(2);
          L1->SetLineWidth(2);
          L2->SetLineWidth(2);
          L1->Draw("SAME");
          L2->Draw("SAME");

          TLegend *leg = new TLegend(0.55, 0.5, 0.85, 0.75);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->AddEntry(hDraw, "K^{+}K^{-} Inv. Mass", "lep");
          leg->AddEntry(fTot, "Total Fit: Exp + Gauss", "l");
          leg->AddEntry(fSig, "Signal (Gauss)", "f");
          leg->AddEntry(fBkg, "Background (Exp)", "f");
          leg->AddEntry(L1, "#pm 3#sigma", "l");
          leg->Draw();

          TLatex latex;
          latex.SetNDC();
          latex.SetTextSize(0.035);
          latex.DrawLatex(0.45, 0.89,
                          Form("Q^{2}[%.2f,%.2f]  %s  -t[%.2f,%.2f]", qLo, qHi,
                               hasW ? Form("W[%.1f,%.1f]", wLo, wHi) : "", tLo,
                               tHi));
          latex.DrawLatex(0.55, 0.85,
                          Form("#mu = %.3f GeV, #sigma = %.3f GeV", mu, sigma));
          const double chi2ndf =
              fTot->GetNDF() > 0 ? fTot->GetChisquare() / fTot->GetNDF() : 0.0;
          latex.DrawLatex(0.55, 0.81, Form("#chi^{2}/ndf = %.2f", chi2ndf));
          latex.DrawLatex(
              0.55, 0.77,
              Form("N_{#phi} (total) = %.1f #pm %.1f", Nsig, Nsig_err));

          auto tagRaw =
              hasW ? Form("Q2_%.1f_%.1f__W_%.1f_%.2f__t_%.1f_%.1f", qLo, qHi,
                          wLo, wHi, tLo, tHi)
                   : Form("Q2_%.1f_%.1f__t_%.1f_%.1f", qLo, qHi, tLo, tHi);
          std::string tag = tagRaw;
          std::replace(tag.begin(), tag.end(), '.', '_');
          c->SaveAs(
              Form("%s/KKmass_%s.pdf", outDirPerModel.c_str(), tag.c_str()));

          // pack draw products
          PhiMassDraw pack;
          pack.h = hDraw;
          pack.fTot = fTot;
          pack.fSig = fSig;
          pack.fBkg = fBkg;
          pack.mu = mu;
          pack.sigma = sigma;
          pack.mLo = mLo;
          pack.mHi = mHi;
          pack.c = c;
          pack.name = tag;
          out[iq][iw][it] = pack;

          // clean transient objects we don't cache (canvas snapshot saved
          // already)
          delete c;
          // delete hDraw;
          if (gPad) {
            gPad->Clear();
            gPad->Update();
            gPad->SetSelected(nullptr);
          }
          delete leg;
          delete L1;
          delete L2;
        } // it (t bins)
      } // iw (W bins)
    } // iq (Q2 bins)
    return out;
  }

private:
  std::vector<std::string> disvars = {"Q2", "xB", "t", "W", "phi"};
  std::map<std::string, std::pair<double, double>> axisRanges = {
      {"Q2", {0.0, 15.0}},
      {"xB", {0.0, 1.0}},
      {"W", {1.0, 10.0}},
      {"t", {0.0, 10.0}},
      {"phi", {-180.0, 180.0}}};
  std::map<std::string, std::pair<double, double>> kinematicAxisRanges = {
      {"recel_p", {-0.05, 13.0}},     {"recel_theta", {-0.01, 1.0}},
      {"recel_phi", {-0.01, 6.1}},    {"recpho_p", {-0.01, 10.0}},
      {"recpho_theta", {-0.01, 1.0}}, {"recpho_phi", {-0.01, 6.1}},
      {"recpro_p", {-0.01, 2.0}},     {"recpro_theta", {-0.01, 2.0}},
      {"recpro_phi", {-0.01, 6.1}}
      // Add more as needed
  };

  DISANAMath kinCalc;
  std::unique_ptr<TFile> predFile;
  std::string predFileName = ".";
  double beam_energy;
  std::string ttreeName;
  std::vector<ROOT::RDF::RResultPtr<TH1>> kinematicHistos, disHistos;
  std::vector<std::vector<TH1D *>> phi_dsdt_QW_;
  ROOT::RDF::RNode rdf;
  AccEffProvider accEffProvider_ = [](double, double, double, double, double,
                                      double) { return 1.0; };
};

#endif // DISANA_PLOTTER_H