// HattaEventGen_fixed.cpp
// Standalone C++ event generator for e p -> e' p K+ K- (via phi production)
// Writes GEMC/CLAS12 LUND ASCII format. Includes optional Hatta–Strikman (2021) weight.
//
// Build:  g++ -O3 -std=c++17 HattaEventGen_fixed.cpp -o epKK_lund
// run ./epKK_lund -n 20000 -E 10.6 --vx 0.0 --vy 0.0 --vz -3.0 -o epKK.lund

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

static constexpr double PI   = 3.14159265358979323846;
static constexpr double Me   = 0.000510999;
static constexpr double Mp   = 0.9382720813;
static constexpr double MK   = 0.493677;
static constexpr double Mphi = 1.019461;

struct FourVec {
  double t;       // E
  double x, y, z; // p
  FourVec(double E = 0, double px = 0, double py = 0, double pz = 0)
      : t(E), x(px), y(py), z(pz) {}
  double E() const { return t; }
  double px() const { return x; }
  double py() const { return y; }
  double pz() const { return z; }
  double m2() const { return t * t - (x * x + y * y + z * z); }
  double M() const { return std::sqrt(std::max(0.0, m2())); }
  double p() const { return std::sqrt(x * x + y * y + z * z); }
  void boost(const double bx, const double by, const double bz) {
    double b2 = bx * bx + by * by + bz * bz;
    if (b2 >= 1.0) return;
    double gamma  = 1.0 / std::sqrt(1.0 - b2);
    double bp     = bx * x + by * y + bz * z;
    double gamma2 = (b2 > 0) ? (gamma - 1.0) / b2 : 0.0;
    double x_ = x + gamma2 * bp * bx + gamma * bx * t;
    double y_ = y + gamma2 * bp * by + gamma * by * t;
    double z_ = z + gamma2 * bp * bz + gamma * bz * t;
    double t_ = gamma * (t + bp);
    x = x_; y = y_; z = z_; t = t_;
  }
};

struct Opts {
  int nEvents = 1000;
  double Ebeam = 10.6;
  double Q2min = 1.0;
  double Q2max = 6.0;
  double ymin  = 0.10;
  double ymax  = 0.85;
  double xmax  = 0.90;
  double vx = 0.0, vy = 0.0, vz = 0.0; // cm
  uint64_t seed = 0ULL;                // 0 -> random_device
  std::string outFile = "epKK.lund";
  double wMargin = 0.02; // GeV above threshold
  bool   flatW = false;
  double Wmin = -1.0, Wmax = -1.0; // used only if flatW
  // HS2021 model switches/params
  bool useHS = true;
  double Ds0 = -0.7;
  double As0 = 0.04;
  double mA = 1.13;
  double mD = 0.76;
  double facHS = 2.5;
  double Knorm = 2.21; // from Eq.(22)
  double c0 = 1.0, c1 = 2.0, c2 = 1.0;
  double longFrac = -1.0; // <0 -> auto via epsilon
  std::string sanityPDF = ""; // if non-empty, generate the two-panel PDF and exit
  double fig1_scale = 1.0;    // multiply Fig1 ds/dt by this (to match nb/GeV^2 if desired)
  double fig2_scale = 1.0;    // multiply Fig2 ds/dt by this (to match pb/GeV^2 if desired)
  bool sanityAlsoRun = false;   // <— new

};

static void usage(const char *prog){
  std::cerr
    << "\nUsage: " << prog << " [options]\n"
    << "  -n <int>          number of events (default 1000)\n"
    << "  -E <GeV>          beam energy (default 10.6)\n"
    << "  --Q2min <GeV2>    min Q^2 (default 1.0)\n"
    << "  --Q2max <GeV2>    max Q^2 (default 6.0)\n"
    << "  --ymin <val>      min y (default 0.10)\n"
    << "  --ymax <val>      max y (default 0.85)\n"
    << "  --xmax <val>      max x_B (default 0.90)\n"
    << "  --vx/--vy/--vz    vertex in cm (defaults 0)\n"
    << "  --seed <u64>      RNG seed (default random)\n"
    << "  -o <file>         output LUND file (default epKK.lund)\n"
    << "  --wMargin <GeV>   extra W threshold margin (default 0.02)\n"
    << "  --flatW [--Wmin --Wmax]\n"
    << "  --useHS           enable Hatta-Strikman 2021 weight\n"
    << "  --Ds0 <val>       D_s(0) (default -0.7)\n"
    << "  --As0 <val>       A_s(0) (default 0.04)\n"
    << "  --mA <GeV>        A dipole mass (default 1.13)\n"
    << "  --mD <GeV>        D tripole mass (default 0.76)\n"
    << "  --facHS <val>     higher-twist factor (default 2.5)\n"
    << "  --Knorm <val>     overall norm (default 2.21)\n"
    << "  --c0/--c1/--c2    A^2, A*D, D^2 coeffs (defaults 1,2,1)\n"
    << "  --longFrac <0-1>  override longitudinal blend (auto if omitted)\n"
    << "  --sanityPDF <file.pdf>  make 2-pad PDF (Fig.1 & Fig.2) and exit\n"
    << "  --fig1scale <val>        scale factor for Fig.1 y-values (default 1)\n"
    << "  --fig2scale <val>        scale factor for Fig.2 y-values (default 1)\n"
    << "  --sanityAlsoRun          after making the PDF, continue to generate events\n";
}

// ROOT 
// ===== Optional plotting support (requires ROOT) =====
#ifdef __has_include
#  if __has_include(<TCanvas.h>)
#    define HATTA_HAVE_ROOT 1
#  else
#    define HATTA_HAVE_ROOT 0
#  endif
#else
#  define HATTA_HAVE_ROOT 0
#endif

#if HATTA_HAVE_ROOT
#  include <TCanvas.h>
#  include <TLegend.h>
#  include <TGraph.h>
#  include <TPad.h>
#  include <TAxis.h>
#  include <TStyle.h>
#  include <TH1.h>  
#endif
// =====================================================

static bool parseArgs(int argc, char **argv, Opts &o) {
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto need = [&](int j, const char *what) {
      if (j >= argc) { std::cerr << "Missing value for " << what << "\n"; std::exit(1); }
      return true;
    };
    if      (a == "-n"       && need(++i, "-n"))        o.nEvents = std::stoi(argv[i]);
    else if (a == "-E"       && need(++i, "-E"))        o.Ebeam   = std::stod(argv[i]);
    else if (a == "--Q2min"  && need(++i, "--Q2min"))   o.Q2min   = std::stod(argv[i]);
    else if (a == "--Q2max"  && need(++i, "--Q2max"))   o.Q2max   = std::stod(argv[i]);
    else if (a == "--ymin"   && need(++i, "--ymin"))    o.ymin    = std::stod(argv[i]);
    else if (a == "--ymax"   && need(++i, "--ymax"))    o.ymax    = std::stod(argv[i]);
    else if (a == "--xmax"   && need(++i, "--xmax"))    o.xmax    = std::stod(argv[i]);
    else if (a == "--vx"     && need(++i, "--vx"))      o.vx      = std::stod(argv[i]);
    else if (a == "--vy"     && need(++i, "--vy"))      o.vy      = std::stod(argv[i]);
    else if (a == "--vz"     && need(++i, "--vz"))      o.vz      = std::stod(argv[i]);
    else if (a == "--seed"   && need(++i, "--seed"))    o.seed    = std::stoull(argv[i]);
    else if (a == "-o"       && need(++i, "-o"))        o.outFile = argv[i];
    else if (a == "--wMargin"&& need(++i, "--wMargin")) o.wMargin = std::stod(argv[i]);
    else if (a == "--flatW")                            o.flatW   = true;
    else if (a == "--Wmin"   && need(++i, "--Wmin"))    o.Wmin    = std::stod(argv[i]);
    else if (a == "--Wmax"   && need(++i, "--Wmax"))    o.Wmax    = std::stod(argv[i]);
    else if (a == "--useHS")                            o.useHS   = true;
    else if (a == "--Ds0"    && need(++i, "--Ds0"))     o.Ds0     = std::stod(argv[i]);
    else if (a == "--As0"    && need(++i, "--As0"))     o.As0     = std::stod(argv[i]);
    else if (a == "--mA"     && need(++i, "--mA"))      o.mA      = std::stod(argv[i]);
    else if (a == "--mD"     && need(++i, "--mD"))      o.mD      = std::stod(argv[i]);
    else if (a == "--facHS"  && need(++i, "--facHS"))   o.facHS   = std::stod(argv[i]);
    else if (a == "--Knorm"  && need(++i, "--Knorm"))   o.Knorm   = std::stod(argv[i]);
    else if (a == "--c0"     && need(++i, "--c0"))      o.c0      = std::stod(argv[i]);
    else if (a == "--c1"     && need(++i, "--c1"))      o.c1      = std::stod(argv[i]);
    else if (a == "--c2"     && need(++i, "--c2"))      o.c2      = std::stod(argv[i]);
    else if (a == "--longFrac"&&need(++i, "--longFrac"))o.longFrac= std::stod(argv[i]);
    else if (a == "--sanityPDF" && need(++i, "--sanityPDF")) o.sanityPDF = argv[i];
    else if (a == "--fig1scale" && need(++i, "--fig1scale")) o.fig1_scale = std::stod(argv[i]);
    else if (a == "--fig2scale" && need(++i, "--fig2scale")) o.fig2_scale = std::stod(argv[i]);
    else if (a == "--sanityAlsoRun") o.sanityAlsoRun = true;   

    else { usage(argv[0]); return false; }
  }
  return true;
}
static inline double hs_weight(double Q2, double W, double t, double xB, double y, const Opts &o);

#if HATTA_HAVE_ROOT
static void MakeSanityTwoPanel(const Opts& base, const char* outpdf) {
  // Paper settings
  const double W_fig = 2.5;          // GeV
  const double Q2_fig1 = 3.8;        // GeV^2
  const double Q2_fig2 = 20.0;       // GeV^2
  const double tmin = 0.0, tmax = 3.0; // |t| range (GeV^2)
  const int    Nt   = 600;
  const double Ds_order[5] = { +0.4, 0.0, -0.4, -0.7, -1.3 };

  // Choose a beam energy to define y and xB for epsilon_LT; this does not change W,Q2
  // Fig.1 was fixed-target-like; Fig.2 mentions sqrt{s}_ep=30 GeV but our weight
  // uses (Q2,W, t, xB, y) — once Q2,W are fixed, nu is fixed, and y depends on E.
  auto y_xB_from_Q2W_E = [](double Q2, double W, double Ebeam){
    const double nu = (W*W + Q2 - Mp*Mp)/(2.0*Mp);
    const double y  = nu / Ebeam;
    const double xB = (Q2 > 0 && nu>0) ? Q2/(2.0*Mp*nu) : 0.0;
    return std::pair<double,double>(y, xB);
  };

  // Beam choices (feel free to tweak if you know the paper’s exact E’s)
  const double E_fig1 = 10.6;  // CLAS12-like
  const double E_fig2 = 15.0;  // arbitrary collider-ish; only enters epsilon via y

  // ROOT canvas
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_phi_sanity", "phi sanity", 720, 960);
  c->Divide(1,2);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.01);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);

  auto draw_panel = [&](int ipad, double Q2, double Ebeam, const char* ttl, const char* yunit, double yscale){
    c->cd(ipad);
    gPad->SetTicks(1,1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);
    gPad->SetMin
    gPad->SetLogy();

    auto [y, xB] = y_xB_from_Q2W_E(Q2, W_fig, Ebeam);
    // Build curves
    std::vector<TGraph*> graphs;
    graphs.reserve(5);
    double ymax = 0.0;

    for (int i = 0; i < 5; ++i) {
      Opts o = base;
      o.Ds0 = Ds_order[i];
      auto gr = new TGraph();
      gr->SetLineWidth(2);
      gr->SetLineColor(800 + i*20); // distinct colors
      int idx = 0;
      for (int it = 0; it < Nt; ++it) {
        const double t_abs = tmin + (tmax - tmin) * (double(it)/(Nt-1));
        const double t_mand = -t_abs; // our hs_weight expects t (Mandelstam); sign-convention: take negative
        double w = hs_weight(Q2, W_fig, t_mand, xB, y, o) * yscale;
        if (!std::isfinite(w)) continue;
        gr->SetPoint(idx++, t_abs, std::max(0.0, w));
        ymax = std::max(ymax, w);
      }
      graphs.push_back(gr);
    }

    // Frame
    auto *frame = gPad->DrawFrame(tmin, 0.0, tmax, (ymax>0 ? 1.15*ymax : 1.0));
    frame->GetXaxis()->SetTitle("|t|  [GeV^{2}]");
    frame->GetYaxis()->SetTitle(Form("d#sigma/dt  [%s]", yunit));
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->SetTitle(ttl);

    // Legend in paper order
    TLegend* leg = new TLegend(0.62, 0.58, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    for (int i = 0; i < 5; ++i) {
      graphs[i]->Draw("L SAME");
      leg->AddEntry(graphs[i], Form("D_{s} = %.1f", Ds_order[i]), "l");
    }
    leg->Draw();
    gPad->SetGrid(1,1);
  };

  draw_panel(1, Q2_fig1, E_fig1,
             Form("W = %.2f GeV, Q^{2} = %.1f GeV^{2}", W_fig, Q2_fig1),
             "nb / GeV^{2}", base.fig1_scale);

  draw_panel(2, Q2_fig2, E_fig2,
             Form("W = %.2f GeV, Q^{2} = %.0f GeV^{2}", W_fig, Q2_fig2),
             "pb / GeV^{2}", base.fig2_scale);

  c->SaveAs(outpdf);
  fprintf(stderr, "[sanity] wrote %s\n", outpdf);
}
#endif

static inline double uniform(std::mt19937_64 & rng, double a, double b) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  return a + (b - a) * U(rng);
}

static inline double twoBodyMomentum(double W, double m1, double m2) {
  double s = W * W;
  double l = (s - (m1 + m2) * (m1 + m2)) * (s - (m1 - m2) * (m1 - m2));
  if (l <= 0) return 0.0;
  return 0.5 / W * std::sqrt(l);
}

// ---------- HS2021 model helpers ----------
static inline double dipoleA(double As0, double t, double mA) {
  return As0 / std::pow(1.0 - t / (mA * mA), 2.0);
}
static inline double tripoleD(double Ds0, double t, double mD) {
  return Ds0 / std::pow(1.0 - t / (mD * mD), 3.0);
}

// epsilon (longitudinal-to-transverse flux ratio), Eq.(6)-style with target-mass gamma
static inline double epsilon_LT(double Q2, double xB, double y) {
  double gamma = 2.0 * xB * Mp / std::sqrt(std::max(1e-12, Q2));
  double num = 1.0 - y - (y * y * gamma * gamma) / 4.0;
  double den = 1.0 - y + (y * y) / 2.0 + (y * y * gamma * gamma) / 4.0;
  if (den == 0) return 0.0;
  return num / den;
}

// Simplified HS weight (shape in t): see chat notes for mapping to paper
/*static inline double hs_weight(double Q2, double W, double t, double xB,
                               double y, const Opts &o) {
  double A = dipoleA(o.As0, t, o.mA);
  double D = tripoleD(o.Ds0, t, o.mD);
  double tau = std::abs(t) / (Mp * Mp);
  double base = o.c0 * A * A + o.c1 * A * D * tau + o.c2 * D * D * tau * tau;
  double kin = 1.0 / std::max(0.2, Q2);
  double eps = epsilon_LT(Q2, xB, y);
  double Lfrac = (o.longFrac >= 0.0 ? o.longFrac : std::min(1.0, std::max(0.0, eps)));
  double longBoost = 1.0 + 0.3 * Lfrac;
  double w = o.Knorm * o.facHS * base * kin * longBoost;
  return std::max(0.0, w);
}*/

// paper-like dσ/dt shape at fixed (Q2,W,t):
static inline double hs_weight(double Q2, double W, double t, double xB,
                               double y, const Opts &o) {
  // use mandelstam t (<0); paper’s form factors use t directly
  const double ms = 0.10;                         // Fig.1 choice in paper
  const double denom = Q2 + ms*ms;
  const double kernel = 1.0 / std::pow(denom, 4); // (Q^2 + m_s^2)^(-4)

  // Dipole / tripole per Eq.(20), As(0)=0.04
  const double A = o.As0 / std::pow(1.0 - t/(o.mA*o.mA), 2.0);
  const double D = o.Ds0 / std::pow(1.0 - t/(o.mD*o.mD), 3.0);

  // D-term enters with Δ factors; use τ = (-t)/(4 m_N^2) from Eq.(19) scaling
  const double tau = std::max(0.0, (-t)) / (4.0 * Mp * Mp);

  // linear combination noted under Figs.1–2: A^2, A*D, D^2
  const double base = o.c0 * (A*A) + o.c1 * (A*D) * tau + o.c2 * (D*D) * (tau*tau);

  // Overall normalization: factor 2.5 (higher-spin) × 2.21 (width), Eq.(22)
  const double norm = o.facHS * o.Knorm;

  // No extra L/T “boost” beyond what ε already does in the exact formula.
  // We keep only a weak kinematic safety clamp; shape is governed by base(t).
  const double w = norm * kernel * base;

  return (w > 0.0 && std::isfinite(w)) ? w : 0.0;
}


int main(int argc, char **argv) {
  long long nGenerated = 0;     // number actually written
  long long nTried     = 0;     // attempts (optional diagnostic)
  long double sumW  = 0.0L;     // Σ w_i
  long double sumW2 = 0.0L;     // Σ w_i^2  (for error estimate)
  Opts opts; if (!parseArgs(argc, argv, opts)) return 1;
#if HATTA_HAVE_ROOT
  if (!opts.sanityPDF.empty()) {
    MakeSanityTwoPanel(opts, opts.sanityPDF.c_str());
    if (!opts.sanityAlsoRun) return 0;   // old behavior unless flag is set
  }
#else
  if (!opts.sanityPDF.empty()) {
    std::cerr << "Rebuild with ROOT to enable --sanityPDF plotting.\n";
    if (!opts.sanityAlsoRun) return 3;   // bail if they only wanted the PDF
    // else fall through to event gen
  }
#endif


  std::mt19937_64 rng;
  if (opts.seed == 0) { std::random_device rd; rng.seed(((uint64_t)rd() << 32) ^ rd()); }
  else rng.seed(opts.seed);

  std::ofstream out(opts.outFile);
  if (!out) { std::cerr << "Cannot open output file: " << opts.outFile << "\n"; return 2; }
  out.setf(std::ios::fixed);

  const double Ethr = Mp + Mphi + opts.wMargin;

  // W-range if flatW
  double Wmin_auto = std::sqrt(std::max(0.0, Mp*Mp + 2*Mp*opts.Ebeam*opts.ymin - opts.Q2max));
  double Wmax_auto = std::sqrt(std::max(0.0, Mp*Mp + 2*Mp*opts.Ebeam*opts.ymax - opts.Q2min));
  if (opts.Wmin < 0) opts.Wmin = std::max(Ethr, Wmin_auto);
  if (opts.Wmax < 0) opts.Wmax = std::max(opts.Wmin + 0.05, Wmax_auto);

  // Lab 4-vectors
  FourVec beam(opts.Ebeam, 0, 0, opts.Ebeam);
  FourVec targ(Mp, 0, 0, 0);

  std::uniform_real_distribution<double> U01(0.0, 1.0);

  int written = 0;
  int attempts = 0;

  while (written < opts.nEvents) {
    attempts++;
    double Q2, y, xB, W2, W;
    if (opts.flatW) {
      Q2 = uniform(rng, opts.Q2min, opts.Q2max);
      W  = uniform(rng, opts.Wmin, opts.Wmax);
      W2 = W * W;
      double nu = (W2 + Q2 - Mp*Mp)/(2*Mp);
      y = nu / opts.Ebeam;
      if (y <= opts.ymin || y >= opts.ymax) continue;
      xB = Q2 / (2*Mp*nu);
    } else {
      Q2 = uniform(rng, opts.Q2min, opts.Q2max);
      y  = uniform(rng, opts.ymin,  opts.ymax);
      double nu = y*opts.Ebeam;
      xB = Q2 / (2*Mp*nu);
      if (xB <= 0 || xB > opts.xmax) continue;
      W2 = Mp*Mp + 2*Mp*nu - Q2;
      if (W2 <= 0) continue;
      W = std::sqrt(W2);
    }
    if (W < Ethr) continue;

    // e' kinematics
    double Eprime = opts.Ebeam - y*opts.Ebeam;
    if (Eprime <= Me) continue;
    double arg = Q2/(4.0*opts.Ebeam*Eprime);
    if (arg > 1.0) continue;
    double theta_e = 2.0 * std::asin(std::sqrt(std::min(1.0, std::max(0.0, arg))));
    double phi_e   = 2*PI*U01(rng);
    double pe      = std::sqrt(std::max(0.0, Eprime*Eprime - Me*Me));
    FourVec eprime(Eprime,
                   pe*std::sin(theta_e)*std::cos(phi_e),
                   pe*std::sin(theta_e)*std::sin(phi_e),
                   pe*std::cos(theta_e));

    // Virtual photon and hadronic system
    FourVec q(beam.E() - eprime.E(), beam.px() - eprime.px(),
              beam.py() - eprime.py(), beam.pz() - eprime.pz());
    double q2 = q.m2();
    if (q2 >= 0) continue;
    FourVec had(targ.E()+q.E(), targ.px()+q.px(), targ.py()+q.py(), targ.pz()+q.pz());
    double What = had.M();
    if (What < (Mp + Mphi)) continue;

    // Boost to hadronic CM
    double bx = had.px()/had.E();
    double by = had.py()/had.E();
    double bz = had.pz()/had.E();

    // had -> p + phi (isotropic in CM)
    double pstar = twoBodyMomentum(What, Mp, Mphi);
    if (pstar <= 0) continue;
    double cth = 2.0*U01(rng) - 1.0;
    double sth = std::sqrt(std::max(0.0, 1.0 - cth*cth));
    double phi = 2*PI*U01(rng);
    double px_cm = pstar*sth*std::cos(phi);
    double py_cm = pstar*sth*std::sin(phi);
    double pz_cm = pstar*cth;
    FourVec p_prot_cm(std::sqrt(Mp*Mp + pstar*pstar), +px_cm, +py_cm, +pz_cm);
    FourVec phi_cm  (std::sqrt(Mphi*Mphi + pstar*pstar), -px_cm, -py_cm, -pz_cm);

    // phi -> K+ K- (isotropic in phi rest frame)
    double kstar = twoBodyMomentum(Mphi, MK, MK);
    double cth2  = 2.0*U01(rng) - 1.0;
    double sth2  = std::sqrt(std::max(0.0, 1.0 - cth2*cth2));
    double phi2  = 2*PI*U01(rng);
    FourVec kp_phi(std::sqrt(MK*MK + kstar*kstar),
                   kstar*sth2*std::cos(phi2),
                   kstar*sth2*std::sin(phi2),
                   kstar*cth2);
    FourVec km_phi(kp_phi.E(), -kp_phi.x, -kp_phi.y, -kp_phi.z);

    // Boost kaons up to CM and then lab
    double bphi_x = phi_cm.px()/phi_cm.E();
    double bphi_y = phi_cm.py()/phi_cm.E();
    double bphi_z = phi_cm.pz()/phi_cm.E();
    kp_phi.boost(bphi_x, bphi_y, bphi_z);
    km_phi.boost(bphi_x, bphi_y, bphi_z);

    FourVec p_prot = p_prot_cm;
    FourVec phi_lab = phi_cm;
    p_prot.boost(bx, by, bz);
    kp_phi.boost(bx, by, bz);
    km_phi.boost(bx, by, bz);
    phi_lab.boost(bx, by, bz);

    if (p_prot.E() < Mp || kp_phi.E() < MK || km_phi.E() < MK) continue;

    // Weight
    double weight = 1.0;
    double nu = y*opts.Ebeam;
    double xB_use = (Q2 > 0 ? Q2/(2*Mp*nu) : 0.0);
    FourVec P_in  = targ;
    FourVec P_out = p_prot;
    double tmand = (P_out.t - P_in.t)*(P_out.t - P_in.t)
                 - (P_out.x - P_in.x)*(P_out.x - P_in.x)
                 - (P_out.y - P_in.y)*(P_out.y - P_in.y)
                 - (P_out.z - P_in.z)*(P_out.z - P_in.z);
    if (opts.useHS) weight = hs_weight(Q2, W, tmand, xB_use, y, opts);
    sumW  += weight;
    sumW2 += weight * weight;
    nGenerated++;
    // LUND record
    const int Npart = 4; // e', p, K+, K-
    int A = 1, Z = 1, beamType = 11, nucleonID = 2212, processID = 7701;

    out << std::setprecision(4);
    out << Npart << "  " << A << "  " << Z << "  " << 0.0 << "  " << 0.0
        << " " << beamType << "  " << std::setprecision(3) << opts.Ebeam
        << "   " << nucleonID << "      " << processID << "      "
        << std::setprecision(7) << weight << "\n";

    auto pline = [&](int idx, int pdg, const FourVec &v, double mass) {
      out << std::setprecision(1) << idx << "  " << 0.0 << "  " << 1 << "  "
          << pdg << "  " << 0 << "  " << 0 << "  " << std::setprecision(4)
          << std::setw(8) << v.px() << "  " << std::setw(8) << v.py()
          << "  " << std::setw(8) << v.pz() << "  " << std::setw(8) << v.E()
          << "  " << std::setw(7) << mass << "  " << std::setprecision(4)
          << std::setw(6) << opts.vx << " " << std::setw(6) << opts.vy
          << " " << std::setw(8) << opts.vz << "\n";
    };

    pline(1,  11, eprime, Me);
    pline(2,2212, p_prot, Mp);
    pline(3, 321, kp_phi, MK);
    pline(4,-321, km_phi, MK);

    written++;
  }

  std::cerr << "Generated " << written << " events to " << opts.outFile
            << " (attempts=" << attempts << ")\n";
   // Final summary & MC-equivalent luminosity
  const long double sigma_tot = sumW;                 // in same units as 'weight'
  const long double sigma_err = std::sqrt(sumW2);     // conservative
  const long double N         = static_cast<long double>(nGenerated);
  const bool hasXsec          = (sigma_tot > 0.0L);

  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(6);
  std::cout << "\n========== Event Generation Summary ==========\n";
  std::cout << "Events generated (written): " << static_cast<long long>(N) << "\n";
  std::cout << "Total cross-section (sum w): " << sigma_tot
            << "  +/- " << sigma_err << "  [units of 'weight']\n";
  if (hasXsec) {
    const long double L_eff = N / sigma_tot;  // [1 / units of 'weight']
    std::cout << "MC-equivalent luminosity N/sigma: " << L_eff
              << "  [1 / units of 'weight']\n";
  } else {
    std::cout << "MC-equivalent luminosity: undefined (total cross-section is 0)\n";
  }
  std::cout << "Generation attempts: " << attempts << "\n";
  std::cout << "==============================================\n";
  return 0;
}
