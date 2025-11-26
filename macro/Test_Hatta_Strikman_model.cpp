// plot_dsigdt_twopanel.C
// Two-panel ROOT macro to draw theoretical dσ/dt vs |t| using the
// “paper-like” Hatta–Strikman (2021) shape embedded in your generator.
// Panel 1:  W=2.5 GeV, Q2=3.8 GeV^2  (nb/GeV^2)
// Panel 2:  W=2.5 GeV, Q2=20  GeV^2  (pb/GeV^2)
//
// Usage:
//   root -l -b -q 'Test_Hatta_Strikman_model.C("phi_dsigdt_twopanels.pdf")'

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TH1.h"

namespace HS {

static constexpr double Mp = 0.9382720813;  // GeV

struct Opts {
  // Model parameters (paper-consistent defaults)
  double Ds0   = -0.7;   // will be overridden by scan
  double As0   = 0.04;   // A_s(0)
  double mA    = 1.13;   // GeV
  double mD    = 0.76;   // GeV
  double facHS = 2.5;    // higher-spin factor
  double Knorm = 2.21;   // Eq.(22)
  double c0 = 1.0, c1 = 2.0, c2 = 1.0; // not fixed by the paper; kept tunable
  double ms = 0.10;      // GeV (used in (Q^2 + m_s^2)^-4)
};

static inline double dipoleA(double As0, double t, double mA) {
  return As0 / std::pow(1.0 - t/(mA*mA), 2.0);
}
static inline double tripoleD(double Ds0, double t, double mD) {
  return Ds0 / std::pow(1.0 - t/(mD*mD), 3.0);
}

// paper-like kernel; expects Mandelstam t<0
static inline double dsigma_dt_shape(double Q2, double /*W*/, double t, const Opts& o) {
  const double kernel = 1.0 / std::pow(Q2 + o.ms*o.ms, 4);
  const double A = dipoleA(o.As0, t, o.mA);
  const double D = tripoleD(o.Ds0, t, o.mD);
  const double tau = std::max(0.0, -t) / (4.0 * Mp * Mp);
  const double base = o.c0*(A*A) + o.c1*(A*D)*tau + o.c2*(D*D)*(tau*tau);
  const double norm = o.facHS * o.Knorm;
  const double w = norm * kernel * base;
  return (w>0.0 && std::isfinite(w)) ? w : 0.0;
}

} // namespace HS

void draw_panel(int ipad,
                double W, double Q2,
                const char* yunit_label,
                double yscale,
                const std::vector<double>& Ds_scan,
                double tmin_abs = 0.0,   // <— NEW
                double tmax_abs = 3.0,
                int Nt = 600,
                bool logY = true)
{
  using namespace HS;
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.08);
  if (logY) gPad->SetLogy();
  gPad->SetGrid(1,1);

  std::vector<TGraph*> curves;
  curves.reserve(Ds_scan.size());
  double ymin = 1e30, ymax = 0.0;

  int color_idx = 0;
  for (double Ds0 : Ds_scan) {
    Opts o; o.Ds0 = Ds0;
    auto* g = new TGraph();
    g->SetLineWidth(2);
    g->SetLineColor(800 + 20*(++color_idx));

    int ip = 0;
    for (int i = 0; i < Nt; ++i) {
      const double t_abs = tmin_abs + (tmax_abs - tmin_abs) * (Nt==1 ? 0.0 : double(i)/(Nt-1));
      const double t_man = -t_abs; // Mandelstam t
      double w = dsigma_dt_shape(Q2, W, t_man, o) * yscale;
      if (!std::isfinite(w)) continue;
      w = std::max(w, logY ? 1e-16 : 0.0);
      g->SetPoint(ip++, t_abs, w);
      ymin = std::min(ymin, w);
      ymax = std::max(ymax, w);
    }
    curves.push_back(g);
  }

  const double ylo = (logY ? std::max((ymin>0 ? 0.6*ymin : 1e-16), 1e-16) : 0.0);
  const double yhi = (ymax>0 ? 1.25*ymax : 1.0);
  auto* frame = gPad->DrawFrame(tmin_abs, ylo, tmax_abs, yhi);  // <— start at tmin_abs
  frame->GetXaxis()->SetTitle("|t|  [GeV^{2}]");
  frame->GetYaxis()->SetTitle(Form("d#sigma/dt  [%s]", yunit_label));
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->SetTitle("");

  TLatex tl;
  tl.SetTextSize(0.042);
  tl.DrawLatexNDC(0.14, 0.93,
                  Form("W = %.2f GeV,  Q^{2} = %.1f GeV^{2}", W, Q2));

  auto* leg = new TLegend(0.62, 0.58, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.040);
  for (size_t i = 0; i < curves.size(); ++i) {
    curves[i]->Draw("L SAME");
    leg->AddEntry(curves[i], Form("D_{s}(0) = %+0.1f", Ds_scan[i]), "l");
  }
  leg->Draw();
}

void Test_Hatta_Strikman_model(const char* outpdf = "phi_dsigdt_twopanels.png",
                          // panel 1 (JLab-like)
                          double W1 = 2.5, double Q2_1 = 3.8,  double yscale_1 = 1.0,
                          // panel 2 (EIC-like)
                          double W2 = 2.5, double Q2_2 = 20.0, double yscale_2 = 1.0,
                          // NEW: common t-range for both panels
                          double tmin_abs = 0.5, double tmax_abs = 3.0)
{
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_phi_twopanels", "phi d#sigma/dt", 800, 1000);
  c->Divide(1,2);

  std::vector<double> Ds_scan = {+0.4, 0.0, -0.4, -0.7, -1.3};

  c->cd(1);
  draw_panel(1, W1, Q2_1, "nb / GeV^{2}", yscale_1, Ds_scan, tmin_abs, tmax_abs, 600, true);

  c->cd(2);
  draw_panel(2, W2, Q2_2, "pb / GeV^{2}", yscale_2, Ds_scan, tmin_abs, tmax_abs, 600, true);

  c->SaveAs(outpdf);
  ::Info("plot_dsigdt_twopanel", "Wrote %s", outpdf);
}
