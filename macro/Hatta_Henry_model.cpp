// Hatta_Henry_model.cpp
// DVMP j=1 (NLO) longitudinal dσ/dt from H,E  + transverse via R(Q), total with ε.
// Two panels: (W,Q^2) chosen by caller; |t| in [max(0.5,|tmin|), min(3.0,|tmax|)].
//
// Example:
//   root -l -b -q 'Hatta_Henry_model.cpp++("phi_Hatta_Henry_total_twopanels.pdf", 10.6, 2.5, 3.8, 2.5, 20.0)'

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>   // std::pair
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"

namespace DVMP {

static constexpr double pi  = 3.14159265358979323846;
static constexpr double M   = 0.9382720813;   // proton [GeV]
static constexpr double Mphi= 1.019461;       // phi [GeV]
static constexpr double alpha_em = 1.0/137.035999084;

// ------- Paper parameters (Sec. IV; Eqs. 47–49) -------
struct Params {
  double mA = 1.6;      // GeV (dipole mass)
  double mD = 1.1;      // GeV (tripole mass)
  double As0 = 0.03;    // A_s(0)
  double Ag0 = 0.42;    // A_g(0)
  double Auc0= 0.55;    // A_{u+d+c}(0)
  double Ds0 = 0.0;     // D_s(0) (user scan)
  double Dg0 = -1.0;    // D_g(0) (user)
  double Duc0= -1.2;    // D_{u+d+c}(0)
  double fphi= 0.221;   // GeV (Eq. 46)
  double es   = -1.0/3.0;
  double CF   = 4.0/3.0;
  double Nc   = 3.0;
  int    nf   = 4;
  double mu   = -1.0;   // default μ=Q
};

// 1-loop αs (Eq. 45), fixed by αs/(2π)=0.0606 at Q^2=2.5 GeV^2
static inline double alpha_s_1loop(double mu2){
  const int nf = 4; const double Nc=3.0;
  const double beta0 = (11.0*Nc - 2.0*nf)/3.0;
  const double Q2_ref=2.5, alpha_ref=2.0*pi*0.0606;
  const double Lref  = 4.0*pi/(beta0*alpha_ref);
  const double Lambda2 = Q2_ref*std::exp(-Lref);
  const double L = std::log(mu2/Lambda2);
  return 4.0*pi/(beta0*L);
}

// Dipole/tripole GFFs (Eq. 47–48)
static inline double A_form(double A0,double t,double mA){return A0/std::pow(1.0 - t/(mA*mA),2);}
static inline double D_form(double D0,double t,double mD){return D0/std::pow(1.0 - t/(mD*mD),3);}

// γ*p CM kinematics (Eqs. 3,5–8) and t-range; (we won't use ξ from here)
struct Kin { double pcm,kcm,tmin,tmax; };
static Kin kin_WQ(double W,double Q2){
  const double W2=W*W;
  const double pcm2=(W2*W2 - 2.0*W2*(M*M - Q2) + (M*M + Q2)*(M*M + Q2))/(4.0*W2);
  const double kcm2=((W2-(Mphi+M)*(Mphi+M))*(W2-(Mphi-M)*(Mphi-M)))/(4.0*W2);
  const double pcm=(pcm2>0? std::sqrt(pcm2):0.0);
  const double kcm=(kcm2>0? std::sqrt(kcm2):0.0);
  const double common = - std::pow((Q2+Mphi*Mphi)/(2.0*W),2);
  const double tmin = common + std::pow(pcm - kcm,2);
  const double tmax = common + std::pow(pcm + kcm,2);
  return {pcm,kcm,tmin,tmax};
}

// κ (Eq. 33): es * CF * fφ / (Nc * Q)
static inline double kappa(double Q,const Params& P){ return P.es*P.CF*P.fphi/(P.Nc*Q); }

// NLO j=1 amplitudes H,E (Eqs. 35–36), with μ=Q
struct Amp{ double H,E; };
static Amp HE_threshold_NLO(double Q2,double t,double xi,const Params& P){
  const double Q=std::sqrt(Q2);
  const double mu2=(P.mu>0? P.mu*P.mu: Q2);
  const double as = alpha_s_1loop(mu2);
  const double LQ = std::log(Q2/mu2);
  const double kapp= kappa(Q,P);
  // GFF combos
  const double As = A_form(P.As0, t,P.mA);
  const double Ag = A_form(P.Ag0, t,P.mA);
  const double Aqs= A_form(P.Auc0, t,P.mA)+As;
  const double Ds = D_form(P.Ds0, t,P.mD);
  const double Dg = D_form(P.Dg0, t,P.mD);
  const double Dqs= D_form(P.Duc0,t,P.mD)+P.Ds0;

  const double nf=P.nf;
  // H coefficients
  const double Hs_LO  = as;
  const double Hs_NLO = (as*as)/(2.0*pi)*(25.7309 - 2.0*nf + (-131.0/18.0 + nf/3.0)*LQ);
  const double Hps    = (as*as)/(2.0*pi)*(-2.3889 + (2.0/3.0)*LQ);
  const double Hg_LO  = (3.0/8.0)*as;
  const double Hg_NLO = (3.0/8.0)*(as*as)/(2.0*pi)*(13.8682 - (83.0/18.0)*LQ);

  const double H_s  = (Hs_LO+Hs_NLO)*(As + xi*xi*Ds);
  const double H_ps = Hps*(Aqs + xi*xi*Dqs);
  const double H_g  = (Hg_LO+Hg_NLO)*(Ag + xi*xi*Dg);

  // Overall factor (Eqs. 33,35): 2 κ / ξ^2 * (15/2)
  const double OF = (15.0 * kapp)/(xi*xi);
  double H   = OF * (H_s + H_ps + H_g);

  // E: B-terms neglected; keep -ξ^2 D pieces
  const double E_s  = (Hs_LO+Hs_NLO)*(-xi*xi*Ds);
  const double E_ps = Hps*(-xi*xi*Dqs);
  const double E_g  = (Hg_LO+Hg_NLO)*(-xi*xi*Dg);
  double Eamp = OF * (E_s + E_ps + E_g);

  if (!std::isfinite(H) || !std::isfinite(Eamp)) return {0.0, 0.0};
  return {H,Eamp};
}

// Longitudinal differential cross section (Eq. 26)
static inline double dsigmaL_dt(double W,double Q2,double t,double xi,const Amp& A,double pcm){
  const double W2=W*W;
  const double pref=(2.0*pi*pi*alpha_em)/((W2 - M*M)*W*pcm);
  double quad = (1.0 - xi*xi)*A.H*A.H
              - ((t/(4.0*M*M)) + xi*xi)*A.E*A.E
              - 2.0*xi*xi*(A.H*A.E);
  if (!(quad >= 0.0) || !std::isfinite(quad)) quad = 0.0;
  return pref*quad;
}

// Bjorken x_B from (W,Q^2)
static inline double xB_from_WQ(double W,double Q2){
  const double W2=W*W;
  return Q2/(W2 - M*M + Q2);
}

// epsilon(y, gamma) with guards; returns {eps, y}
// s_ep: squared ep COM energy
static inline std::pair<double,double> epsilon_guarded(double W,double Q2,double s_ep){
  const double xB = xB_from_WQ(W, Q2);
  const double y  = Q2 / (s_ep * xB);
  const double gamma = 2.0*M*xB/std::sqrt(Q2);
  double num = 1.0 - y - (y*y*gamma*gamma)/4.0;
  double den = 1.0 - y + (y*y)/2.0 + (y*y*gamma*gamma)/4.0;
  double eps = (den!=0.0? num/den : 0.0);
  // If y is out of [0,1], warn and clamp eps into [0,1)
  if (!(y >= 0.0 && y <= 1.0) || !std::isfinite(eps)) {
    ::Warning("epsilon","Unphysical y=%.3f for W=%.2f Q2=%.2f at this beam energy; showing photon-level σ only.", y, W, Q2);
    eps = 0.9; // placeholder if you still want a combined curve; or skip combining
  }
  eps = std::max(0.0, std::min(0.999, eps));
  return {eps, y};
}

} // namespace DVMP

// Robust skewness from x_B; clip to avoid xi=0 or xi>=1
static inline double xi_from_WQ(double W, double Q2){
  const double W2 = W*W;
  const double xB = Q2 / (W2 - DVMP::M*DVMP::M + Q2);
  double xi = xB / (2.0 - xB);
  if (!std::isfinite(xi) || xi <= 0.0) xi = 1e-3;
  return std::max(1e-3, std::min(0.95, xi));
}

// Numeric integrate by trapezoid on [tmin_abs, tmax_abs]
static double trapz(const std::vector<double>& x,const std::vector<double>& y){
  double s=0.0;
  for (size_t i=1;i<x.size();++i) s += 0.5*(y[i]+y[i-1])*(x[i]-x[i-1]);
  return s;
}

static void draw_panel(double W,double Q2,const char* ylab,double Dg0,double Ebeam_GeV,
                       double t_abs_min=0.5,double t_abs_max_cap=3.0,int Nt=500,bool logY=true)
{
  using namespace DVMP;
  gPad->SetTicks(1,1); gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.08);
  if (logY) gPad->SetLogy();

  // Fixed-target s_ep (electron on proton at rest)
  const double s_ep = M*M + 2.0*M*Ebeam_GeV;

  // Kinematics & limits
  const DVMP::Kin K = kin_WQ(W,Q2);
  double tmin_abs = (K.tmin<0 ? std::sqrt(K.tmin*K.tmin) : 0.0);
  double tmax_abs = (K.tmax<0 ? std::sqrt(K.tmax*K.tmax) : 0.0);
  tmin_abs = std::max(tmin_abs, t_abs_min);
  tmax_abs = std::min(tmax_abs, t_abs_max_cap);
  if (tmin_abs>=tmax_abs){ tmax_abs=tmin_abs+0.01; }

  // Photon polarization ε (Eq. 13)
  const auto [eps, y] = epsilon_guarded(W, Q2, s_ep);

  // Ratio R(Q) = σ_L/σ_T (CLAS fit, Eq. 59); override here if desired
  const double R = 0.4*Q2/(DVMP::Mphi*DVMP::Mphi);

  // Scan Ds(0)
  const std::vector<double> Ds_scan = {+0.4, 0.0, -0.4, -0.7, -1.3};

  std::vector<TGraph*> g_tot; g_tot.reserve(Ds_scan.size());
  double ymin=1e30, ymax=0.0;

  // Also collect integrands for totals
  for (size_t k=0;k<Ds_scan.size();++k){
    DVMP::Params P; P.Ds0 = Ds_scan[k]; P.Dg0 = Dg0;

    std::vector<double> vx, vL, vT, vTOT;
    vx.reserve(Nt); vL.reserve(Nt); vT.reserve(Nt); vTOT.reserve(Nt);

    auto* g = new TGraph(); g->SetLineWidth(2); g->SetLineColor(800 + 20*(int)k);

    for (int i=0;i<Nt;++i){
      const double u = (Nt==1?0.0: double(i)/(Nt-1));
      const double t_abs = tmin_abs + (tmax_abs - tmin_abs)*u;
      const double t = -t_abs;

      // robust ξ from x_B (not from tmin)
      const double xi = xi_from_WQ(W, Q2);

      DVMP::Amp A = DVMP::HE_threshold_NLO(Q2,t,xi,P);
      const double dsL = std::max( DVMP::dsigmaL_dt(W,Q2,t,xi,A,K.pcm), 0.0 );
      const double dsT = dsL / std::max(1e-12, R);          // assume weak t-dependence of R
      const double dsTot = dsT + eps*dsL;                   // dσ/dt = dσ_T/dt + ε dσ_L/dt

      vx.push_back(t_abs); vL.push_back(dsL); vT.push_back(dsT); vTOT.push_back(dsTot);

      double yv = std::max(dsTot, logY? 1e-16 : 0.0);
      g->SetPoint(i, t_abs, yv);
      ymin = std::min(ymin, yv);
      ymax = std::max(ymax, yv);
    }

    // Integrate and print
    const double sigL  = trapz(vx, vL);
    const double sigT  = trapz(vx, vT);
    const double sigTot= trapz(vx, vTOT);
    ::Info("panel","W=%.2f Q2=%.2f Ds(0)=%.1f | eps=%.3f R=%.3f  ->  sigma_L=%.4g, sigma_T=%.4g, sigma_tot=%.4g [%s]",
           W,Q2, Ds_scan[k], eps, R, sigL, sigT, sigTot, ylab);

    g_tot.push_back(g);
  }

  // Axes owner and limits
  g_tot.front()->SetTitle("");
  g_tot.front()->GetXaxis()->SetTitle("|t|  [GeV^{2}]");
  g_tot.front()->GetYaxis()->SetTitle(Form("d#sigma/dt  [%s]", ylab));
  g_tot.front()->Draw("AL");
  gPad->Update();
  if (auto* h = g_tot.front()->GetHistogram()){
    h->GetXaxis()->SetRangeUser(tmin_abs, tmax_abs);
    const double ylo=(logY? std::max((ymin>0?0.6*ymin:1e-16),1e-16):0.0);
    const double yhi=(ymax>0? 1.25*ymax : 1.0);
    h->SetMinimum(ylo); h->SetMaximum(yhi);
  }
  for (size_t i=1;i<g_tot.size();++i) g_tot[i]->Draw("L SAME");

  TLatex tl; tl.SetTextSize(0.042);
  tl.DrawLatexNDC(0.14,0.93,Form("W = %.2f GeV,  Q^{2} = %.1f GeV^{2}", W,Q2));

  auto* leg = new TLegend(0.58,0.56,0.90,0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
  for (size_t i=0;i<Ds_scan.size();++i)
    leg->AddEntry(g_tot[i], Form("D_{s}(0) = %+0.1f", Ds_scan[i]), "l");
  leg->AddEntry((TObject*)nullptr, Form("#varepsilon = %.3f", eps), "");
  leg->AddEntry((TObject*)nullptr, Form("R(Q)=%.3f", 0.4*Q2/(DVMP::Mphi*DVMP::Mphi)), "");
  leg->Draw();
}

void Hatta_Henry_model(const char* outpdf="phi_Hatta_Henry_total_twopanels.pdf",
                       double Ebeam_GeV = 10.6,       // electron beam energy (fixed-target ep)
                       // panel 1 (nb): JLab-like
                       double W1=2.5, double Q2_1=3.8,
                       // panel 2 (pb): EIC-like
                       double W2=2.5, double Q2_2=20.0,
                       // gluon D-term choice
                       double Dg0 = -1.0,
                       // |t| range
                       double tmin_abs=0.5, double tmax_abs=4.0)
{
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_phi_total", "phi total d#sigma/dt", 850, 1000);
  c->Divide(1,2);

  c->cd(1);
  draw_panel(W1, Q2_1, "nb / GeV^{2}", Dg0, Ebeam_GeV, tmin_abs, tmax_abs, 500, true);

  c->cd(2);
  draw_panel(W2, Q2_2, "pb / GeV^{2}", Dg0, Ebeam_GeV, tmin_abs, tmax_abs, 500, true);

  c->SaveAs(outpdf);
}
