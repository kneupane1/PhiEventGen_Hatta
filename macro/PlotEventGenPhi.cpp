#include <TApplication.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TString.h>

#include <ROOT/RDataFrame.hxx>

#include "../src/DISANAMath.h"
#include "../src/DISANAcomparer.h"
#include "../src/DrawStyle.h"

// --- local helpers (needed in this TU) ---
struct Part {
  int pdg = 0;
  double px = 0, py = 0, pz = 0, E = 0, m = 0;
};
static inline double rad2deg(double r) {
  return r * 180.0 / 3.14159265358979323846;
}
static inline double Deg2Rad(double d) { return d * M_PI / 180.0; }

// Build the final columns requested, starting from p/th/phi (in degrees) that
// your init() wrote.
ROOT::RDF::RNode AddFinalEPKKColumns(ROOT::RDF::RNode df, double EbeamGeV) {
  // physical masses (GeV)
  constexpr double m_e = 0.000510999;
  constexpr double m_p = 0.9382720813;
  constexpr double m_k = 0.493677;

  // === 1) components from (p,theta,phi) in degrees ===
  auto d = df
               // electron components
               .Define("ele_px",
                       [](double p, double th_deg, double ph_deg) {
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::cos(ph);
                       },
                       {"p_e", "th_e", "ph_e"})
               .Define("ele_py",
                       [](double p, double th_deg, double ph_deg) {
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::sin(ph);
                       },
                       {"p_e", "th_e", "ph_e"})
               .Define("ele_pz",
                       [](double p, double th_deg) {
                         double th = Deg2Rad(th_deg);
                         return p * std::cos(th);
                       },
                       {"p_e", "th_e"})

               // proton components
               .Define("pro_px",
                       [](double p, double th_deg, double ph_deg) {
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::cos(ph);
                       },
                       {"p_p", "th_p", "ph_p"})
               .Define("pro_py",
                       [](double p, double th_deg, double ph_deg) {
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::sin(ph);
                       },
                       {"p_p", "th_p", "ph_p"})
               .Define("pro_pz",
                       [](double p, double th_deg) {
                         double th = Deg2Rad(th_deg);
                         return p * std::cos(th);
                       },
                       {"p_p", "th_p"})

               // K+ components
               .Define("kPlus_px",
                       [](double p, double th_deg, double ph_deg) {
                         if (p < 0)
                           return 0.0; // handle missing
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::cos(ph);
                       },
                       {"p_kp", "th_kp", "ph_kp"})
               .Define("kPlus_py",
                       [](double p, double th_deg, double ph_deg) {
                         if (p < 0)
                           return 0.0;
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::sin(ph);
                       },
                       {"p_kp", "th_kp", "ph_kp"})
               .Define("kPlus_pz",
                       [](double p, double th_deg) {
                         if (p < 0)
                           return 0.0;
                         double th = Deg2Rad(th_deg);
                         return p * std::cos(th);
                       },
                       {"p_kp", "th_kp"})

               // K- components (measured; may be missing in some events)
               .Define("kMinus_px",
                       [](double p, double th_deg, double ph_deg) {
                         if (p < 0)
                           return 0.0;
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::cos(ph);
                       },
                       {"p_km", "th_km", "ph_km"})
               .Define("kMinus_py",
                       [](double p, double th_deg, double ph_deg) {
                         if (p < 0)
                           return 0.0;
                         double th = Deg2Rad(th_deg), ph = Deg2Rad(ph_deg);
                         return p * std::sin(th) * std::sin(ph);
                       },
                       {"p_km", "th_km", "ph_km"})
               .Define("kMinus_pz",
                       [](double p, double th_deg) {
                         if (p < 0)
                           return 0.0;
                         double th = Deg2Rad(th_deg);
                         return p * std::cos(th);
                       },
                       {"p_km", "th_km"})
               .Define("recel_p",
                       [](double px, double py, double pz) {
                         return std::sqrt(px * px + py * py + pz * pz);
                       },
                       {"ele_px", "ele_py", "ele_pz"})
               .Define("recel_theta",
                       [](double px, double py, double pz) {
                         return std::acos(
                             pz / std::sqrt(px * px + py * py + pz * pz));
                       },
                       {"ele_px", "ele_py", "ele_pz"})
               .Define("recel_phi",
                       [](double px, double py) {
                         double ph = std::atan2(py, px);
                         return ph < 0 ? ph + 2 * M_PI : ph;
                       },
                       {"ele_px", "ele_py"})

               .Define("recpro_p",
                       [](double px, double py, double pz) {
                         return std::sqrt(px * px + py * py + pz * pz);
                       },
                       {"pro_px", "pro_py", "pro_pz"})
               .Define("recpro_theta",
                       [](double px, double py, double pz) {
                         return std::acos(
                             pz / std::sqrt(px * px + py * py + pz * pz));
                       },
                       {"pro_px", "pro_py", "pro_pz"})
               .Define("recpro_phi",
                       [](double px, double py) {
                         double ph = std::atan2(py, px);
                         return ph < 0 ? ph + 2 * M_PI : ph;
                       },
                       {"pro_px", "pro_py"})

               .Define("reckPlus_p",
                       [](double px, double py, double pz) {
                         return std::sqrt(px * px + py * py + pz * pz);
                       },
                       {"kPlus_px", "kPlus_py", "kPlus_pz"})
               .Define("reckPlus_theta",
                       [](double px, double py, double pz) {
                         return std::acos(
                             pz / std::sqrt(px * px + py * py + pz * pz));
                       },
                       {"kPlus_px", "kPlus_py", "kPlus_pz"})
               .Define("reckPlus_phi",
                       [](double px, double py) {
                         double ph = std::atan2(py, px);
                         return ph < 0 ? ph + 2 * M_PI : ph;
                       },
                       {"kPlus_px", "kPlus_py"})

               .Define("reckMinus_p",
                       [](double px, double py, double pz) {
                         return std::sqrt(px * px + py * py + pz * pz);
                       },
                       {"kMinus_px", "kMinus_py", "kMinus_pz"})
               .Define("reckMinus_theta",
                       [](double px, double py, double pz) {
                         return std::acos(
                             pz / std::sqrt(px * px + py * py + pz * pz));
                       },
                       {"kMinus_px", "kMinus_py", "kMinus_pz"})
               .Define("reckMinus_phi",
                       [](double px, double py) {
                         double ph = std::atan2(py, px);
                         return ph < 0 ? ph + 2 * M_PI : ph;
                       },
                       {"kMinus_px", "kMinus_py"})
                .Define("invMass_KpKm",
                     [](double px1, double py1, double pz1, double px2, double py2, double pz2) -> float {
                       constexpr float mK = 0.493677;  // Kaon mass in GeV/c²
                       float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mK * mK);
                       float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mK * mK);
                       float px = px1 + px2;
                       float py = py1 + py2;
                       float pz = pz1 + pz2;
                       float E = E1 + E2;
                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_px", "kMinus_py", "kMinus_pz"});

  return d;
}
// Holders so the returned RNode stays valid after the function returns
struct _DFKeepAlive {
  std::string outFilePath;
  std::shared_ptr<TFile> outFile;   // keep the file open
};
static _DFKeepAlive g_hold;
// -----------------------------------------------------------------------------
// Forward declare your final selections (implemented below)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df);
// convenience 3-vector helpers
static double MomentumFunc(float px, float py, float pz) {
  return std::sqrt(px * px + py * py + pz * pz);
}
static double ThetaFunc(float px, float py, float pz) {
  return std::acos(pz / std::sqrt(px * px + py * py + pz * pz));
}
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}
template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string &name,
                               const Method method, float beam_energy) {
  return node.Define(
      name,
      [method, beam_energy](
          double recel_p, double recel_theta, double recel_phi, double recpro_p,
          double recpro_theta, double recpro_phi, double reckMinus_p,
          double reckMinus_theta, double reckMinus_phi, double reckPlus_p,
          double reckPlus_theta, double reckPlus_phi) {
        return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi,
                           recpro_p, recpro_theta, recpro_phi, reckMinus_p,
                           reckMinus_theta, reckMinus_phi, reckPlus_p,
                           reckPlus_theta, reckPlus_phi).*
                method)();
      },
      {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta",
       "recpro_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi",
       "reckPlus_p", "reckPlus_theta", "reckPlus_phi"});
}
// -----------------------------------------------------------------------------
// Plot styling (unchanged)
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2); // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13,
                    0.06); // For DVCS plots
DrawStyle csStyle(0.05, 0.05, .95, 1.1, 42, 5, 510, 0.12, 0.03, 0.12,
                  0.02); // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, .8, .8, 42, 5, 510, 0.15, 0.07, 0.16,
                   0.06); // For BSA

// For exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pro_det_region == 2", "CD"},
    {"pro_det_region == 1", "FD"},
};

ROOT::RDF::RNode init(const std::string &lund_path);
// -----------------------------------------------------------------------------
// Main plotter with toggles
void PlotEventGenPhi() // subset toggle inside missing-mass
{
  //ROOT::EnableImplicitMT();
  float beam_energy_fall2018 = 6.000f;

  // -----------------------------
  // Input locations
  // Exclusive reconstruction K+K-
  std::string input_path_LundFiles =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/dvcsgen/"
      "PhiEventGen/Phi_Hatta_model/outputs/";
  std::string filename_Lund_inputs =
      Form("%s/epKK.lund", input_path_LundFiles.c_str());
  auto df_clas6_all = init(filename_Lund_inputs);
  auto df_clas6 = AddFinalEPKKColumns(df_clas6_all, /*Ebeam=*/10.6);

  // -----------------------------
  // Comparer setup (we’ll add models conditionally below)
  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);
  comparer.PlotIndividual(false);

  // Binning (unchanged)
  BinManager xBins;
  xBins.SetQ2Bins({1.4,  3.8});
  xBins.SetXBBins({0, 0.99});
  xBins.SetWBins({2.0, 3.0});
  xBins.SetTBins({0.2, .3, 0.4, 0.5, .7, .8, 1.0, 1.4, 1.8, 2.5,3.0, 3.5, 4.0, 5.0, 8.0 });
  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);

  // Some global constants you had
  double luminosity_clas6 = 1248.495759;
  double polarisation = 0.85;
  double branching = 0.49;

  // Apply your “final” DVEP selections and then pick exclusive phi event
  auto df_clas6_final = ApplyFinalDVEPSelections(df_clas6);
  comparer.AddModelPhi(df_clas6_final, "Clas 6", beam_energy_fall2018);

  // -----------------------------
  // Shared summary plots
  std::cout << "Applying further cuts and plotting…" << std::endl;
  comparer.PlotPhiElectroProKinematicsComparison();
  comparer.PlotKinematicComparison_phiAna();
  // comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);
  comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", 40, 0.988, 1.15,
                                          true, luminosity_clas6, branching);
  comparer.PlotPhiDSigmaDt_FromCache();
  gApplication->Terminate(0);
}

// -----------------------------------------------------------------------------
// Final DVEP selections (kept, just fixed labels to match values)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("W > 2.0", "Cut: W > 2.0 GeV")
      .Filter("recel_p  > 1.5", "Cut: recel_p < 1.5 GeV")
      //.Filter("bestEle_idx > 0", "Cut: reckPlus_p < 3.5 GeV")
      .Filter("reckPlus_p  < 7.5", "Cut: reckPlus_p < 3.5 GeV")
      .Filter("reckMinus_p < 7.5", "Cut: reckMinus_p < 3.5 GeV")

      // 9. Missing energy / exclusivity
      .Filter("Mx2_eKpKm > 0.8*0.8 && Mx2_eKpKm < 1.08*1.08",
              "Cut: Proton Missing Mass Squared in [0.8,1.08] GeV^2")
      .Filter("Mx2_epKm > .08 && Mx2_epKm < 0.48",
              "Cut: Kaon Missing Mass Squared in [0.08,.48] GeV^2")
      .Filter("Mx2_epKp > .08 && Mx2_epKp < 0.48",
              "Cut: Kaon Missing Mass Squared in [0.08,.48] GeV^2")
      .Filter("Cone_Kp < 6.0", "Cut: cone Angle  K+ < 6°")
      .Filter("Cone_Km < 6.0", "Cut: cone Angle  K- < 6°")
      .Filter("Cone_p < 6.0", "Cut:  cone between p < 6°")
      .Filter("Coplanarity_had_normals_deg <15", "Cut: Coplanarity angle 15°")

      .Filter("PTmiss < 0.120", "Cut: Total Missing PTmiss < .120 GeV")
      .Filter("Mx2_epKpKm < 0.0075",
              "Cut: Total Missing Mass squared < 0.0074 GeV")
      .Filter("Emiss <0.32 && Emiss> -0.175",
              "Cut: Missing energy Emiss < 0.2 GeV");
}

ROOT::RDF::RNode init(const std::string &lund_path) {
  // If we fail to read, return an empty RDataFrame node with 0 entries.
  auto returnEmptyNode = []() -> ROOT::RDF::RNode {
    static ROOT::RDataFrame emptyDF(0);
    return emptyDF;
  };

  // Build an in-memory tree first
  auto tree = std::make_unique<TTree>("lund", "LUND ep->epK+K-");

  // ----------------------------
  // Base per-particle branches
  // ----------------------------
  double p_e = -999, th_e = -999, ph_e = -999;
  double p_p = -999, th_p = -999, ph_p = -999;
  double p_kp = -999, th_kp = -999, ph_kp = -999;
  double p_km = -999, th_km = -999, ph_km = -999;

  double Q2 = 0, Tabs = 0, xB = 0, W = 0, Phi = 0, Wgt = 1.0;

  tree->Branch("p_e", &p_e);
  tree->Branch("th_e", &th_e);
  tree->Branch("ph_e", &ph_e);
  tree->Branch("p_p", &p_p);
  tree->Branch("th_p", &th_p);
  tree->Branch("ph_p", &ph_p);
  tree->Branch("p_kp", &p_kp);
  tree->Branch("th_kp", &th_kp);
  tree->Branch("ph_kp", &ph_kp);
  tree->Branch("p_km", &p_km);
  tree->Branch("th_km", &th_km);
  tree->Branch("ph_km", &ph_km);

  tree->Branch("Q2", &Q2);
  tree->Branch("t", &Tabs);      // |t| = -t
  tree->Branch("xB", &xB);
  tree->Branch("W", &W);
  tree->Branch("phi", &Phi);     // 0..360
  tree->Branch("w", &Wgt);

  // ------------------------------------
  // Exclusivity / geometry branches
  // ------------------------------------
  double mx2_ep = -999, emiss = -999, ptmiss = -999, mx2_epg = -999, mx2_epKpKm = -999;
  double delta_phi = -999, theta_gg = -999, theta_gphi = -999;
  double mx2_egamma = -999, mx2_eKpKm = -999, mx2_epKp = -999, mx2_epKm = -999;
  double theta_e_gamma = -999, theta_e_phi = -999;
  double cone_kp = -999, cone_km = -999, cone_p = -999;
  double coplanarity_had_normals_deg = -999;
  double DeltaE = -999, Mx_epKp = -999;

  tree->Branch("Mx2_ep", &mx2_ep);
  tree->Branch("Emiss", &emiss);
  tree->Branch("PTmiss", &ptmiss);
  tree->Branch("Mx2_epg", &mx2_epg); // fill only if available
  tree->Branch("Mx2_epKpKm", &mx2_epKpKm);

  tree->Branch("DeltaPhi", &delta_phi);
  tree->Branch("Theta_gg", &theta_gg);
  tree->Branch("Theta_g_phimeson", &theta_gphi);

  tree->Branch("Mx2_egamma", &mx2_egamma);
  tree->Branch("Mx2_eKpKm", &mx2_eKpKm);
  tree->Branch("Mx2_epKp", &mx2_epKp);
  tree->Branch("Mx2_epKm", &mx2_epKm);

  tree->Branch("Theta_e_gamma", &theta_e_gamma);
  tree->Branch("Theta_e_phimeson", &theta_e_phi);

  tree->Branch("Cone_Kp", &cone_kp);
  tree->Branch("Cone_Km", &cone_km);
  tree->Branch("Cone_p", &cone_p);

  tree->Branch("Coplanarity_had_normals_deg", &coplanarity_had_normals_deg);
  tree->Branch("DeltaE", &DeltaE);
  tree->Branch("Mx_epKp", &Mx_epKp);

  // ---------- Parse LUND ----------
  std::ifstream in(lund_path);
  if (!in) {
    std::cerr << "[init] ERROR: cannot open LUND file: " << lund_path << "\n";
    return returnEmptyNode();
  }

  struct Part { int pdg=0; double px=0,py=0,pz=0,E=0,m=0; };
  auto rad2deg = [](double r){ return r*180.0/3.14159265358979323846; };

  std::size_t nEv = 0, nFilled = 0;
  while (true) {
    std::string header;
    if (!std::getline(in, header)) break;
    if (header.empty()) continue;

    int Npart=0, A=0, Z=0, beamType=0, nucleonID=0, processID=0;
    double Tpol=0, spinZ=0, Ebeam=0, weight=1.0;
    {
      std::istringstream hs(header);
      hs >> Npart >> A >> Z >> Tpol >> spinZ >> beamType >> Ebeam >> nucleonID >> processID >> weight;
      if (!hs) {
        for (int i=0;i<Npart;i++) { std::string dummy; std::getline(in,dummy); }
        continue;
      }
    }

    Part pe{}, pp{}, pkp{}, pkm{};
    bool gotE=false, gotP=false, gotKp=false, gotKm=false;

    for (int i=0;i<Npart;i++) {
      std::string pline; if (!std::getline(in, pline)) break;
      if (pline.empty()) { --i; continue; }
      std::istringstream ps(pline);
      int idx=0, type=0, pdg=0, parent=0, fd=0; double lifetime=0;
      double px=0, py=0, pz=0, E=0, m=0, vx=0, vy=0, vz=0;
      ps >> idx >> lifetime >> type >> pdg >> parent >> fd >> px >> py >> pz >> E >> m >> vx >> vy >> vz;
      if (!ps) continue;

      if (pdg==11)       { pe  = Part{pdg,px,py,pz,E,m}; gotE=true; }
      else if (pdg==2212){ pp  = Part{pdg,px,py,pz,E,m}; gotP=true; }
      else if (pdg==321) { pkp = Part{pdg,px,py,pz,E,m}; gotKp=true; }
      else if (pdg==-321){ pkm = Part{pdg,px,py,pz,E,m}; gotKm=true; }
    }

    ++nEv;
    if (!(gotE && gotP)) continue;

    auto to_ptp_rad = [](const Part& pr, double& p, double& th, double& ph){
      p = std::sqrt(pr.px*pr.px + pr.py*pr.py + pr.pz*pr.pz);
      double c = (p>0)? pr.pz/p : 1.0;
      if (c> 1.0) c= 1.0; if (c<-1.0) c=-1.0;
      th = std::acos(c);
      ph = std::atan2(pr.py, pr.px);
    };

    double pe_p, pe_th, pe_ph; to_ptp_rad(pe,  pe_p, pe_th, pe_ph);
    double pp_p, pp_th, pp_ph; to_ptp_rad(pp,  pp_p, pp_th, pp_ph);

    double kp_p=0, kp_th=0, kp_ph=0, km_p=0, km_th=0, km_ph=0;
    if (gotKp) to_ptp_rad(pkp, kp_p, kp_th, kp_ph);
    if (gotKm) to_ptp_rad(pkm, km_p, km_th, km_ph);

    // Build kinematics object (φ channel)
    DISANAMath kin(Ebeam, pe_p, pe_th, pe_ph,
                         pp_p, pp_th, pp_ph,
                         km_p, km_th, km_ph,
                         kp_p, kp_th, kp_ph);

    auto set_deg = [&](double th_rad, double ph_rad, double& th_deg, double& ph_deg){
      th_deg = rad2deg(th_rad);
      double phd = rad2deg(ph_rad);
      if (phd>180) phd -= 360;
      if (phd<-180) phd += 360;
      ph_deg = phd;
    };

    // Fill base per-particle cols (degrees)
    p_e = pe_p; set_deg(pe_th, pe_ph, th_e, ph_e);
    p_p = pp_p; set_deg(pp_th, pp_ph, th_p, ph_p);
    if (gotKp){ p_kp=kp_p; set_deg(kp_th, kp_ph, th_kp, ph_kp); } else { p_kp=-999; th_kp=-999; ph_kp=-999; }
    if (gotKm){ p_km=km_p; set_deg(km_th, km_ph, th_km, ph_km); } else { p_km=-999; th_km=-999; ph_km=-999; }

    // DVEP
    Q2   = kin.GetQ2();
    xB   = kin.GetxB();
    W    = kin.GetW();
    Tabs = kin.GetT();      // |t|
    Phi  = kin.GetPhi();    // 0..360
    Wgt  = weight;

    // Exclusivity / geometry
    mx2_ep   = kin.GetMx2_ep();
    emiss    = kin.GetEmiss();
    ptmiss   = kin.GetPTmiss();
    // mx2_epg = kin.GetMx2_epg(); // uncomment if available
    mx2_egamma = kin.GetMx2_egamma();
    theta_e_gamma = kin.GetTheta_e_gamma();
    DeltaE   = kin.GetDeltaE();

    mx2_epKpKm = kin.GetMx2_epKpKm();
    delta_phi  = kin.GetDeltaPhi();
    theta_gg   = kin.GetTheta_gamma_gamma();
    theta_gphi = kin.GetTheta_g_phimeson();
    mx2_eKpKm  = kin.GetMx2_eKpKm();
    mx2_epKp   = kin.GetMx2_epKp();
    mx2_epKm   = kin.GetMx2_epKm();
    theta_e_phi= kin.GetTheta_e_phimeson();
    cone_kp    = kin.GetCone_Kp();
    cone_km    = kin.GetCone_Km();
    cone_p     = kin.GetCone_p();
    coplanarity_had_normals_deg = kin.GetCoplanarity_had_normals_deg();
    Mx_epKp    = kin.GetMx_epKp();

    tree->Fill();
    ++nFilled;
  }

  std::cerr << "[init] Parsed events: " << nEv << ", filled: " << nFilled << "\n";

  // ---------- Persist tree to file and build RDataFrame from file ----------
  g_hold.outFilePath = "epKK_from_lund.root";          // choose your path
  {
    std::unique_ptr<TFile> fout(TFile::Open(g_hold.outFilePath.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      std::cerr << "[init] ERROR: cannot create ROOT file " << g_hold.outFilePath << "\n";
      return returnEmptyNode();
    }
    fout->cd();
    tree->Write();  // writes as "lund"
    fout->Write();
    // close here; we will reopen read-only below
  }

  // Keep the file open while the dataframe is used
  g_hold.outFile.reset(TFile::Open(g_hold.outFilePath.c_str(), "READ"));
  if (!g_hold.outFile || g_hold.outFile->IsZombie()) {
    std::cerr << "[init] ERROR: cannot reopen ROOT file " << g_hold.outFilePath << " for reading\n";
    return returnEmptyNode();
  }

  // Construct RDataFrame from the file-backed tree (works in all ROOT builds)
  static ROOT::RDataFrame df_from_file("lund", g_hold.outFilePath.c_str());
  return df_from_file;
}
