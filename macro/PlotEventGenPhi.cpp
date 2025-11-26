#include <TApplication.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h> // Needed for reading ROOT files
#include <TTree.h> // Needed for reading TTree
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <ROOT/RDataFrame.hxx>

// Include your custom files (assuming they are still needed)
#include "../src/DISANAMath.h"
#include "../src/DISANAcomparer.h"
#include "../src/DrawStyle.h"

// --- Global KeepAlive structure for ROOT file ownership ---
// The TFile needs to be kept open for the RDataFrame to read the TTree data.
struct _DFKeepAlive {
  std::string inFilePath;
  std::shared_ptr<TFile> inFile; // keep the file open
};
static _DFKeepAlive g_hold;


// --- local helpers (needed in this TU) ---
// Note: struct Part and get_ptp_for_pdg removed as we read vectors now.

static inline double rad2deg(double r) {
  return r * 180.0 / 3.14159265358979323846;
}
static inline double Deg2Rad(double d) { return d * M_PI / 180.0; }

/**
 * Helper to define p, theta, phi (in degrees) for a specific particle.
 * This is the new core logic to extract single particle kinematics from
 * the array columns (pid, px, py, pz) read from the ROOT file.
 */
static void get_ptp_deg_from_arrays(const std::vector<int>    &pid,
                                    const std::vector<double> &px,
                                    const std::vector<double> &py,
                                    const std::vector<double> &pz,
                                    int targetPdg,
                                    double &p, double &th_deg, double &ph_deg) {
  p      = -999.0;
  th_deg = -999.0;
  ph_deg = -999.0;
  for (std::size_t i = 0; i < pid.size(); ++i) {
    if (pid[i] != targetPdg) continue;

    double px_i = px[i], py_i = py[i], pz_i = pz[i];
    double pp   = std::sqrt(px_i*px_i + py_i*py_i + pz_i*pz_i);
    double c    = (pp > 0.0) ? pz_i / pp : 1.0;
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;

    double th_rad = std::acos(c);
    double ph_rad = std::atan2(py_i, px_i);

    p = pp;
    th_deg = rad2deg(th_rad);
    
    // Convert phi to -180 to 180 range
    double phd = rad2deg(ph_rad);
    if (phd > 180)  phd -= 360;
    if (phd < -180) phd += 360;
    ph_deg = phd;

    // Found the particle, stop searching
    return;
  }
}

// Build the final columns requested, starting from p/th/phi (in degrees) that
// your init() wrote. (This function is largely unchanged, just the source of p/th/phi changes)
ROOT::RDF::RNode AddFinalEPKKColumns(ROOT::RDF::RNode df, double EbeamGeV) {
  // physical masses (GeV)
  constexpr double m_e = 0.000510999;
  constexpr double m_p = 0.9382720813;
  constexpr double m_k = 0.493677;

  // === 1) components from (p,theta,phi) in degrees ===
  // Note: All columns are defined using the 'new' p_e, th_e, etc. columns
  // which are defined in the updated init() function.
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
               // The rest of the definitions remain the same, calculating P, Theta, Phi from Px, Py, Pz
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

// -----------------------------------------------------------------------------
// Forward declare your final selections (implemented below)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df);

// convenience 3-vector helpers (These can be simplified or removed, 
// as the primary kinematics columns are now created in init() and AddFinalEPKKColumns())
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

// Function signature changed to take root_path instead of lund_path
ROOT::RDF::RNode init(const std::string &root_path);
// -----------------------------------------------------------------------------
// Main plotter with toggles
void PlotEventGenPhi() // subset toggle inside missing-mass
{
  //ROOT::EnableImplicitMT();
  float beam_energy_fall2018 = 10.600f;
  ROOT::EnableImplicitMT(32);
  // -----------------------------
  // Input locations
  // Exclusive reconstruction K+K-
  // NOTE: You must provide the path to the *ROOT* file generated by
  // ConvertOneLundToRootBare/ConvertLundDirToRootBare
  std::string input_path_RootFiles =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/analysis/outputs/test_1/input_rootfiles/";
  
  // *** EXAMPLE: REPLACE THIS WITH YOUR ACTUAL GENERATED ROOT FILE NAME ***
  std::string filename_Root_inputs =
      Form("%s/merged_lund_data.root", input_path_RootFiles.c_str());
      
  auto df_clas6_all = init(filename_Root_inputs);
  
  // Ensure we check if the DataFrame is valid (i.e., if the ROOT file was opened)
  if (df_clas6_all.GetNRuns() == 0) {
      std::cerr << "ERROR: Failed to initialize RDataFrame from ROOT file. Exiting." << std::endl;
      gApplication->Terminate(1);
      return;
  }
  
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
  xBins.SetQ2Bins({0.4,  10.8});
  xBins.SetXBBins({0, 0.99});
  xBins.SetWBins({1.8, 10.6});
  xBins.SetTBins({0.2, .3, 0.4, 0.5, .7, .8, 1.0, 1.4, 1.8, 2.5,3.0, 3.5, 4.0, 5.0, 8.0 });
  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);

  // Some global constants you had
  double luminosity_clas6 = 1263.356542;
  double polarisation = 0.85;
  double branching = 0.49;

  // Apply your “final” DVEP selections and then pick exclusive phi event
  auto df_clas6_final = ApplyFinalDVEPSelections(df_clas6);
  comparer.AddModelPhi(df_clas6_final, "MC 10.6 GeV", beam_energy_fall2018);

  // -----------------------------
  // Shared summary plots
  std::cout << "Applying further cuts and plotting…" << std::endl;
  comparer.PlotPhiElectroProKinematicsComparison();
  comparer.PlotKinematicComparison_phiAna();
  //comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);
  comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", 40, 0.988, 1.15, true, luminosity_clas6, branching);
  comparer.PlotPhiDSigmaDt_FromCache();
  gApplication->Terminate(0);
}

// -----------------------------------------------------------------------------
// Final DVEP selections (kept, just fixed labels to match values)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df) {
  // NOTE: You may need to uncomment the lines below if your DISANAMath
  // produces the kinematics columns (Q2, W, Mx2_eKpKm, etc.)
  // If your ROOT files already contain these columns from a previous step,
  // this is okay, but typically you need to define them via DISANAMath
  // after getting the particle P, Theta, Phi.

  // Since you had them in your original LUND parser, I'll keep them as Filters,
  // assuming DISANAMath calls *will* be used later or were already defined.

  // Define DISANAMath columns (if not already defined)
  // The 'define_DISCAT' template is used to call DISANAMath methods
  ROOT::RDF::RNode d = df;
  float EbeamGeV = 10.600f; // Assuming 10.6 GeV

  // Redefine DIS Kinematics
  d = define_DISCAT(d, "Q2", &DISANAMath::GetQ2, EbeamGeV);
  d = define_DISCAT(d, "W",  &DISANAMath::GetW,  EbeamGeV);
  d = define_DISCAT(d, "xB", &DISANAMath::GetxB, EbeamGeV);
  d = define_DISCAT(d, "t",  &DISANAMath::GetT,  EbeamGeV); // |t|
  d = define_DISCAT(d, "phi", &DISANAMath::GetPhi, EbeamGeV); // 0..360

  // Redefine Exclusivity Kinematics (needed for the commented-out filters below)
  d = define_DISCAT(d, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, EbeamGeV);
  d = define_DISCAT(d, "Mx2_epKp", &DISANAMath::GetMx2_epKp, EbeamGeV);
  d = define_DISCAT(d, "Mx2_epKm", &DISANAMath::GetMx2_epKm, EbeamGeV);
  d = define_DISCAT(d, "Cone_Kp", &DISANAMath::GetCone_Kp, EbeamGeV);
  d = define_DISCAT(d, "Cone_Km", &DISANAMath::GetCone_Km, EbeamGeV);
  d = define_DISCAT(d, "Cone_p", &DISANAMath::GetCone_p, EbeamGeV);
  d = define_DISCAT(d, "Coplanarity_had_normals_deg", &DISANAMath::GetCoplanarity_had_normals_deg, EbeamGeV);
  d = define_DISCAT(d, "PTmiss", &DISANAMath::GetPTmiss, EbeamGeV);
  d = define_DISCAT(d, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, EbeamGeV);
  d = define_DISCAT(d, "Emiss", &DISANAMath::GetEmiss, EbeamGeV);


  return d
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("W > 2.0", "Cut: W > 2.0 GeV")
      .Filter("recel_p  > 1.5", "Cut: recel_p > 1.5 GeV")
      .Filter("reckPlus_p  < 7.5", "Cut: reckPlus_p < 7.5 GeV")
      .Filter("reckMinus_p < 7.5", "Cut: reckMinus_p < 7.5 GeV");

      // Uncomment these lines if you want the exclusivity cuts applied:
      /*.Filter("Mx2_eKpKm > 0.8*0.8 && Mx2_eKpKm < 1.08*1.08",
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
              "Cut: Missing energy Emiss < 0.2 GeV");*/
}

/**
 * Initializes RDataFrame from a ROOT file produced by ConvertOneLundToRootBare.
 * It uses the vector columns (pid, px, py, pz) to define the single-particle
 * columns (p_e, th_e, ph_e, etc.) needed by the rest of the analysis.
 */
ROOT::RDF::RNode init(const std::string &root_path) {
  // If we fail to read, return an empty RDataFrame node with 0 entries.
  auto returnEmptyNode = []() -> ROOT::RDF::RNode {
    static ROOT::RDataFrame emptyDF(0);
    return emptyDF;
  };

  // Open the ROOT file and store the pointer in the global "keep alive" struct.
  // We use TFile::Open to handle network paths if needed.
  g_hold.inFile = std::shared_ptr<TFile>(TFile::Open(root_path.c_str()));
  if (!g_hold.inFile || g_hold.inFile->IsZombie()) {
    std::cerr << "[init] ERROR: cannot open ROOT file: " << root_path << "\n";
    g_hold.inFile.reset(); // clear the failed file pointer
    return returnEmptyNode();
  }

  // Get the TTree
  TTree *tree = nullptr;
  g_hold.inFile->GetObject("lund", tree);
  if (!tree) {
    std::cerr << "[init] ERROR: TTree 'lund' not found in file: " << root_path << "\n";
    g_hold.inFile.reset();
    return returnEmptyNode();
  }

  // Create RDataFrame from the TTree
  ROOT::RDataFrame df(*tree);
  std::cout << "[init] Opened ROOT file and TTree 'lund' with " << df.Count().GetValue() << " entries.\n";

  // --- Particle Data Group (PDG) codes ---
  constexpr int PDG_E  = 11;    // Electron
  constexpr int PDG_P  = 2212;  // Proton
  constexpr int PDG_KP = 321;   // K+
  constexpr int PDG_KM = -321;  // K-
  
  // Lambda function to call the helper for RDataFrame::Define
  auto get_particle_kin = [](const std::vector<int>& pid,
                             const std::vector<double>& px,
                             const std::vector<double>& py,
                             const std::vector<double>& pz,
                             int targetPdg) -> std::array<double, 3> {
      double p, th, ph;
      get_ptp_deg_from_arrays(pid, px, py, pz, targetPdg, p, th, ph);
      return {p, th, ph};
  };

  // Define a new column that holds {p, th_deg, ph_deg} for each particle
  auto df_with_kin = df
    // Electron
    .Define("e_kin", [=](const std::vector<int>& pid, const std::vector<double>& px, 
                         const std::vector<double>& py, const std::vector<double>& pz) {
             return get_particle_kin(pid, px, py, pz, PDG_E);
           }, {"pid", "px", "py", "pz"})
    // Proton
    .Define("p_kin", [=](const std::vector<int>& pid, const std::vector<double>& px, 
                         const std::vector<double>& py, const std::vector<double>& pz) {
             return get_particle_kin(pid, px, py, pz, PDG_P);
           }, {"pid", "px", "py", "pz"})
    // K+
    .Define("kp_kin", [=](const std::vector<int>& pid, const std::vector<double>& px, 
                          const std::vector<double>& py, const std::vector<double>& pz) {
             return get_particle_kin(pid, px, py, pz, PDG_KP);
           }, {"pid", "px", "py", "pz"})
    // K-
    .Define("km_kin", [=](const std::vector<int>& pid, const std::vector<double>& px, 
                          const std::vector<double>& py, const std::vector<double>& pz) {
             return get_particle_kin(pid, px, py, pz, PDG_KM);
           }, {"pid", "px", "py", "pz"});
  
  // --- Split the {p, th_deg, ph_deg} array into the original column names ---
  auto df_final_cols = df_with_kin
    // Electron
    .Define("p_e",  "e_kin[0]")
    .Define("th_e", "e_kin[1]")
    .Define("ph_e", "e_kin[2]")
    // Proton
    .Define("p_p",  "p_kin[0]")
    .Define("th_p", "p_kin[1]")
    .Define("ph_p", "p_kin[2]")
    // K+
    .Define("p_kp",  "kp_kin[0]")
    .Define("th_kp", "kp_kin[1]")
    .Define("ph_kp", "kp_kin[2]")
    // K-
    .Define("p_km",  "km_kin[0]")
    .Define("th_km", "km_kin[1]")
    .Define("ph_km", "km_kin[2]");

  // The original ROOT file contains the event-level LUND header fields directly
  // which can be renamed/kept as they are. The converter code uses the names:
  // Q2, t, xB, W, phi, w. The ROOT file produced by the converter *doesn't*
  // have these, but the LUND parser *did* calculate them via DISANAMath.
  // Since we've replaced the LUND parser with a ROOT reader, these columns
  // must be redefined later (which is done in ApplyFinalDVEPSelections).
  // The columns present in the ROOT file are: Npart, A, Z, Tpol, SpinZ, Ebeam, etc.
  // We keep the essential ones:
  return df_final_cols
    //.Define("Ebeam", "(double)Ebeam")
    .Define("Wgt", "(double)Weight");
}