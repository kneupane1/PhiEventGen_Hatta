#ifndef DISANAMATH_H
#define DISANAMATH_H

// ROOT + std
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TVector2.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
constexpr double pi = 3.14159265358979323846;
const double m_e = 0.000511;       // GeV
//const double m_e = 0.0;       // GeV
const double m_p = 0.938272;       // GeV
//const double m_p = 0.938;       // GeV
const double m_kMinus = 0.493677;  // GeV
const double m_kPlus = 0.493677;   // GeV

// -----------------------------------------------------------------------------
// (p,theta,phi) helpers
// -----------------------------------------------------------------------------
namespace {
TVector3 SphericalToCartesian(double p, double theta, double phi) {
  const double px = p * std::sin(theta) * std::cos(phi);
  const double py = p * std::sin(theta) * std::sin(phi);
  const double pz = p * std::cos(theta);
  return TVector3(px, py, pz);
}
TLorentzVector Build4Vector(double p, double theta, double phi, double mass) {
  TVector3 v = SphericalToCartesian(p, theta, phi);
  const double E = std::sqrt(v.Mag2() + mass * mass);
  return TLorentzVector(v, E);
}
}  // namespace

// -----------------------------------------------------------------------------
// Bin manager 
// -----------------------------------------------------------------------------
class BinManager {
 public:
  BinManager() {
    q2_bins_ = {1.0, 2.0, 4.0, 6.0};
    t_bins_ = {0.1, 0.3, 0.6, 1.0};
    xb_bins_ = {0.1, 0.2, 0.4, 0.6};
    W_bins_ = {0.1, 10.0};
  }
  const std::vector<double> &GetQ2Bins() const { return q2_bins_; }
  const std::vector<double> &GetTBins() const { return t_bins_; }
  const std::vector<double> &GetXBBins() const { return xb_bins_; }
  const std::vector<double> &GetWBins() const { return W_bins_; }
  void SetQ2Bins(const std::vector<double> &v) { q2_bins_ = v; }
  void SetTBins(const std::vector<double> &v) { t_bins_ = v; }
  void SetXBBins(const std::vector<double> &v) { xb_bins_ = v; }
  void SetWBins(const std::vector<double> &v) { W_bins_ = v; }

 private:
  std::vector<double> q2_bins_, t_bins_, xb_bins_, W_bins_;
};

// -----------------------------------------------------------------------------
// DISANAMath
// -----------------------------------------------------------------------------
class DISANAMath {
 private:
  // Kinematics
  double Q2_{}, xB_{}, t_{}, phi_deg_{}, W_{}, nu_{}, y_{};

  // Exclusivity
  double mx2_ep_{};
  double emiss_{};
  double ptmiss_{};
  double mx2_epg_{};
  double mx2_epKpKm_{};
  double delta_phi_{};
  double theta_gg_{};
  double theta_gphi_{};
  double mx2_egamma_{};
  double mx2_eKpKm_{};
  double mx2_epKp_{-1.};
  double mx2_epKm_{};
  double Theta_e_gamma_{};
  double Theta_e_phimeson_{};
  double Cone_Kp_{};
  double Cone_Km_{};
  double Cone_p_{};
  double coplanarity_had_normals_deg_{};

  double DeltaE_{};

 public:
  DISANAMath() = default;
  // φ→K+K− analysis (both kaons given)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double kMinus_p, double kMinus_theta,
             double kMinus_phi, double kPlus_p, double kPlus_theta, double kPlus_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector kMinus = Build4Vector(kMinus_p, kMinus_theta, kMinus_phi, m_kMinus);
    TLorentzVector kPlus = Build4Vector(kPlus_p, kPlus_theta, kPlus_phi, m_kPlus);
    ComputeKinematics(e_in, e_out, p_in, p_out, kPlus, kMinus);
  }
  // Accessors
  double GetQ2() const { return Q2_; }
  double GetxB() const { return xB_; }
  double GetT() const { return t_; }
  double GetPhi() const { return phi_deg_; }
  double GetW() const { return W_; }
  double GetNu() const { return nu_; }
  double Gety() const { return y_; }

  double GetMx2_ep() const { return mx2_ep_; }
  double GetEmiss() const { return emiss_; }
  double GetPTmiss() const { return ptmiss_; }
  double GetMx2_epKpKm() const { return mx2_epKpKm_; }
  double GetDeltaPhi() const { return delta_phi_; }
  double GetTheta_gamma_gamma() const { return theta_gg_; }
  double GetMx2_egamma() const { return mx2_egamma_; }
  double GetMx2_eKpKm() const { return mx2_eKpKm_; }
  double GetMx2_epKp() const { return mx2_epKp_; }
  double GetMx2_epKm() const { return mx2_epKm_; }
  double GetTheta_e_gamma() const { return Theta_e_gamma_; }
  double GetTheta_g_phimeson() const { return theta_gphi_; }
  double GetTheta_e_phimeson() const { return Theta_e_phimeson_; }
  double GetDeltaE() const { return DeltaE_; }
  double GetMx_epKp() const { return (mx2_epKp_ > 0) ? std::sqrt(mx2_epKp_) : -999.0; }
  double GetCone_Kp() const { return Cone_Kp_; }
  double GetCone_Km() const { return Cone_Km_; }
  double GetCone_p() const { return Cone_p_; }

  double GetCoplanarity_had_normals_deg() const { return coplanarity_had_normals_deg_; }

  // Helpers
  double ComputePhiH(const TVector3 &q1_, const TVector3 &k1_, const TVector3 &q2_) const {
    const double t1 = ((q1_.Cross(k1_)).Dot(q2_)) / std::abs((q1_.Cross(k1_)).Dot(q2_));
    const TVector3 t2 = q1_.Cross(k1_);
    const TVector3 t3 = q1_.Cross(q2_);
    const double n2 = t2.Mag(), n3 = t3.Mag();
    const double t4 = (t2.Dot(t3)) / (n2 * n3);
    return t1 * std::acos(t4) * 180. / pi + 180.;
  }

  // Core (phi)
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &kPlus, const TLorentzVector &kMinus) {
    TLorentzVector q = electron_in - electron_out;
    TLorentzVector phi = kPlus + kMinus;

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());

    TVector3 nL = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 nH = q.Vect().Cross(proton_out.Vect()).Unit();
    const double cos_phi = nL.Dot(nH);
    const double sin_phi = (nL.Cross(nH)).Dot(q.Vect().Unit());
    phi_deg_ = (std::atan2(sin_phi, cos_phi) + pi) * 180. / pi;

    const TLorentzVector totI = electron_in + proton_in;
    const TLorentzVector totF = electron_out + proton_out + phi;
    const TLorentzVector miss = totI - totF;

    mx2_ep_ = (totI - electron_out - proton_out).Mag2();
    emiss_ = miss.E();
    ptmiss_ = miss.Vect().Perp();
    mx2_epKpKm_ = miss.Mag2();

    TVector3 qv = q.Vect();
    TVector3 ev = electron_in.Vect();
    TVector3 phiv = phi.Vect();
    TVector3 pv = proton_out.Vect();
    TVector3 kpv = kPlus.Vect();
    TVector3 kmv = kMinus.Vect();
    TVector3 n_qphi = qv.Cross(phiv);
    TVector3 n_pphi = pv.Cross(phiv);

    delta_phi_ = std::abs(ComputePhiH(qv, ev, phiv) - ComputePhiH(qv, ev, -pv));
    theta_gphi_ = phiv.Angle((totI - (electron_out + proton_out)).Vect()) * 180. / pi;
    Cone_p_ = (electron_in + proton_in - electron_out - phi).Angle(pv) * 180. / pi;  // p vs K+
    Cone_Km_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Angle(kmv) * 180. / pi;  // p vs K−
    Cone_Kp_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Angle(kpv) * 180. / pi;  // p vs K−

    mx2_eKpKm_ = (electron_in + proton_in - electron_out - phi).Mag2();
    mx2_epKp_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Mag2();
    mx2_epKm_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Mag2();

    Theta_e_phimeson_ = electron_out.Angle(phi.Vect()) * 180. / pi;
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + phi.E());
    double coplanarity_had_normals_deg_ = std::numeric_limits<double>::quiet_NaN();
    if (n_qphi.Mag() > 0 && n_pphi.Mag() > 0) {
      n_qphi = n_qphi.Unit();
      n_pphi = n_pphi.Unit();

      double c = n_qphi.Dot(n_pphi);
      // numerical safety
      c = std::max(-1.0, std::min(1.0, c));
      double ang = std::acos(c) * 180. / pi; // in [0, 180]

      // treat ±normals as equivalent (optional but common)
      coplanarity_had_normals_deg_ = std::min(ang, 180.0 - ang);
    }
  }
  // φ dσ/dt 
  std::vector<TH1D *> ComputePhi_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const size_t n_q2 = q2_bins.size() - 1, n_t = t_bins.size() - 1;
    std::vector<TH1D *> hist(n_q2, nullptr);

    for (size_t iq = 0; iq < n_q2; ++iq) {
      const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
      auto hname = Form("phi_cs_q2bin%zu", iq);
      auto htitle = Form("d#sigma/dt for Q^{2}=[%.2f, %.2f]", qmin, qmax);
      hist[iq] = new TH1D(hname, htitle, n_t, &t_bins[0]);
      hist[iq]->SetDirectory(nullptr);
      hist[iq]->GetXaxis()->SetTitle("-t [GeV^{2}]");
      hist[iq]->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    }

    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };
    auto fill = [&](double Q2, double t) {
      int iq = findBin(Q2, q2_bins), it = findBin(t, t_bins);
      if (iq >= 0 && it >= 0) hist[iq]->Fill(t);
    };
    df.Foreach(fill, {"Q2", "t"});

    for (auto *h : hist) {
      if (!h) continue;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b), N = h->GetBinContent(b);
        const double ds = N / (luminosity * w), es = std::sqrt(std::max(0.0, N)) / (luminosity * w);
        h->SetBinContent(b, ds);
        h->SetBinError(b, es);
      }
    }
    timer.Stop();
    std::cout << "[ComputePhi_CrossSection] Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
  }
};

#endif  // DISANAMATH_H
