#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TSystem.h>
#include <TList.h>
#include <TCollection.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ========================
// Single-file converter
// ========================
bool ConvertOneLundToRootBare(const std::string &lund_path,
                              const std::string &root_path) {
  std::ifstream in(lund_path);
  if (!in) {
    std::cerr << "[ConvertOneLundToRootBare] ERROR: cannot open "
              << lund_path << "\n";
    return false;
  }

  // --- Output file and tree
  TFile out(root_path.c_str(), "RECREATE");
  if (out.IsZombie()) {
    std::cerr << "[ConvertOneLundToRootBare] ERROR: cannot create "
              << root_path << "\n";
    return false;
  }

  TTree tree("lund", "Bare LUND content");

  // -------- Event-level LUND header fields ----------
  Int_t   Npart = 0, A = 0, Z = 0;
  Int_t   BeamType = 0, NucleonID = 0, ProcessID = 0;
  Double_t Tpol = 0, SpinZ = 0, Ebeam = 0, Weight = 1.0;

  tree.Branch("Npart",      &Npart,      "Npart/I");
  tree.Branch("A",          &A,          "A/I");
  tree.Branch("Z",          &Z,          "Z/I");
  tree.Branch("Tpol",       &Tpol,       "Tpol/D");
  tree.Branch("SpinZ",      &SpinZ,      "SpinZ/D");
  tree.Branch("BeamType",   &BeamType,   "BeamType/I");
  tree.Branch("Ebeam",      &Ebeam,      "Ebeam/D");
  tree.Branch("NucleonID",  &NucleonID,  "NucleonID/I");
  tree.Branch("ProcessID",  &ProcessID,  "ProcessID/I");
  tree.Branch("Weight",     &Weight,     "Weight/D");

  // -------- Per-particle arrays ----------
  std::vector<int>    pid, type, parent, daughter;
  std::vector<double> px, py, pz, E, m, vx, vy, vz;

  tree.Branch("pid",      &pid);
  tree.Branch("type",     &type);
  tree.Branch("parent",   &parent);
  tree.Branch("daughter", &daughter);
  tree.Branch("px",       &px);
  tree.Branch("py",       &py);
  tree.Branch("pz",       &pz);
  tree.Branch("E",        &E);
  tree.Branch("m",        &m);
  tree.Branch("vx",       &vx);
  tree.Branch("vy",       &vy);
  tree.Branch("vz",       &vz);

  std::size_t nEv = 0;

  while (true) {
    std::string header;
    if (!std::getline(in, header)) break;
    if (header.empty()) continue;

    // Parse header:
    // Npart A Z Tpol SpinZ BeamType Ebeam NucleonID ProcessID Weight
    {
      std::istringstream hs(header);
      hs >> Npart >> A >> Z >> Tpol >> SpinZ >> BeamType
         >> Ebeam >> NucleonID >> ProcessID >> Weight;
      if (!hs) {
        // malformed header -> skip supposed body lines
        for (int i = 0; i < Npart; ++i) {
          std::string dummy;
          std::getline(in, dummy);
        }
        continue;
      }
    }

    // reset particle vectors
    pid.clear(); type.clear(); parent.clear(); daughter.clear();
    px.clear(); py.clear(); pz.clear(); E.clear(); m.clear();
    vx.clear(); vy.clear(); vz.clear();

    // Read Npart particle lines
    for (int i = 0; i < Npart; ++i) {
      std::string pline;
      if (!std::getline(in, pline)) break;
      if (pline.empty()) {
        --i;
        continue;
      }

      std::istringstream ps(pline);
      int idx = 0, t = 0, pdg = 0, par = 0, fd = 0;
      double lifetime = 0;
      double px_i = 0, py_i = 0, pz_i = 0, E_i = 0, m_i = 0;
      double vx_i = 0, vy_i = 0, vz_i = 0;

      ps >> idx >> lifetime >> t >> pdg >> par >> fd
         >> px_i >> py_i >> pz_i >> E_i >> m_i >> vx_i >> vy_i >> vz_i;
      if (!ps) continue;

      pid.push_back(pdg);
      type.push_back(t);
      parent.push_back(par);
      daughter.push_back(fd);
      px.push_back(px_i);
      py.push_back(py_i);
      pz.push_back(pz_i);
      E.push_back(E_i);
      m.push_back(m_i);
      vx.push_back(vx_i);
      vy.push_back(vy_i);
      vz.push_back(vz_i);
    }

    tree.Fill();
    ++nEv;
  }

  std::cout << "[ConvertOneLundToRootBare] " << lund_path
            << " -> " << root_path
            << " (events: " << nEv << ")\n";

  out.cd();
  tree.Write();
  out.Write();
  out.Close();
  return true;
}

//========================
// Directory-level driver
// ========================
void ConvertLundDirToRootBare() {
  // Hard-coded input/output directories
  std::string inDir  =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/"
      "PhiEventGen/analysis/lund_outputs/lager_runs";
  std::string outDir =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/"
      "PhiEventGen/analysis/outputs/test_1/input_rootfiles/";
  
  // Define the name for your merged output file
  std::string mergedFileName = "merged_lund_data.root";

  // No ExpandPathName – we stay with std::string
  std::string inDirExp  = inDir;
  std::string outDirExp = outDir; // always use the explicit outDir here

  gSystem->mkdir(outDirExp.c_str(), true);

  TSystemDirectory dir(inDirExp.c_str(), inDirExp.c_str());
  TList *files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "[ConvertLundDirToRootBare] No files in " << inDirExp << "\n";
    return;
  }

  std::cout << "[ConvertLundDirToRootBare] Input dir:  " << inDirExp  << "\n"
            << "[ConvertLundDirToRootBare] Output dir: " << outDirExp << "\n";

  TIter next(files);
  while (auto *file = (TSystemFile *)next()) {
    std::string fname = file->GetName();
    if (file->IsDirectory()) continue;
    if (fname.size() < 4)   continue;                 // need at least ".txt"

    // Only process *.txt (your cleaned LUND-ish files)
    if (fname.substr(fname.size() - 4) != ".txt") continue;

    std::string inPath  = inDirExp  + "/" + fname;
    std::string base    = fname.substr(0, fname.size() - 4); // strip ".txt"
    std::string outPath = outDirExp + "/" + base + ".root";

    // 1. Convert the single LUND file to a single ROOT file
    ConvertOneLundToRootBare(inPath, outPath);
  }

  std::cout << "[ConvertLundDirToRootBare] Conversion done.\n";

  // --- 2. Add the hadd command here ---
  TString haddCommand = TString::Format(
    "hadd -f %s/%s %s/*.root", 
    outDirExp.c_str(), 
    mergedFileName.c_str(), 
    outDirExp.c_str()
  );

  std::cout << "[ConvertLundDirToRootBare] Merging files with command:\n"
            << "  " << haddCommand << "\n";

  int status = gSystem->Exec(haddCommand.Data());
  
  if (status == 0) {
      std::cout << "[ConvertLundDirToRootBare] Successfully merged files into " 
                << outDirExp << "/" << mergedFileName << ".\n";
  } else {
      std::cerr << "[ConvertLundDirToRootBare] ERROR: hadd command failed with status " 
                << status << ". Check command path and directory permissions.\n";
  }
  
  // --- End of hadd command ---

  std::cout << "[ConvertLundDirToRootBare] Done.\n";
}