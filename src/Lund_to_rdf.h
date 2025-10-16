#include "DISPhiMath.h"  // your header that defines DISANAMath

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TTree.h>

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <map>

struct Part { int pdg=0; double px=0,py=0,pz=0,E=0,m=0; };
static inline double rad2deg(double r){ return r*180.0/3.14159265358979323846; }

static int ColElectron() { return TColor::GetColor("#1f77b4"); } // blue
static int ColProton()   { return TColor::GetColor("#d62728"); } // red
static int ColKplus()    { return TColor::GetColor("#2ca02c"); } // green
static int ColKminus()   { return TColor::GetColor("#9467bd"); } // purple

static void StyleHist(TH1D* h, int col){
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetFillColorAlpha(col, 0.15);
  h->SetTitleFont(42, "");
  h->SetTitleSize(0.050, "");
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.038);
  h->GetYaxis()->SetLabelSize(0.038);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.15);
}

// Build a TTree with one row per event and columns for kinematics.
// Columns (all doubles): p_e,th_e,ph_e, p_p,th_p,ph_p, p_kp,th_kp,ph_kp, p_km,th_km,ph_km,
// Q2, t, xb, W, phi, w
static TTree* LundToTree(const char* fname){
  auto *tree = new TTree("lund", "LUND rows");

  // Branch vars
  double p_e=-999, th_e=-999, ph_e=-999;
  double p_p=-999, th_p=-999, ph_p=-999;
  double p_kp=-999, th_kp=-999, ph_kp=-999;
  double p_km=-999, th_km=-999, ph_km=-999;
  double Q2=0, Tabs=0, xB=0, W=0, Phi=0, Wgt=1.0;

  tree->Branch("p_e", &p_e);   tree->Branch("th_e",&th_e); tree->Branch("ph_e",&ph_e);
  tree->Branch("p_p", &p_p);   tree->Branch("th_p",&th_p); tree->Branch("ph_p",&ph_p);
  tree->Branch("p_kp",&p_kp);  tree->Branch("th_kp",&th_kp); tree->Branch("ph_kp",&ph_kp);
  tree->Branch("p_km",&p_km);  tree->Branch("th_km",&th_km); tree->Branch("ph_km",&ph_km);
  tree->Branch("Q2",&Q2); tree->Branch("t",&Tabs); tree->Branch("xB",&xB);
  tree->Branch("W",&W);  tree->Branch("phi",&Phi); tree->Branch("w",&Wgt);

  std::ifstream in(fname);
  if(!in){ std::cerr << "ERROR: cannot open " << fname << std::endl; return tree; }

  while(true){
    std::string header;
    if(!std::getline(in, header)) break;
    if(header.empty()) continue;

    int Npart=0, A=0, Z=0, beamType=0, nucleonID=0, processID=0;
    double Tpol=0, spinZ=0, Ebeam=0, weight=1.0;
    {
      std::istringstream hs(header);
      hs >> Npart >> A >> Z >> Tpol >> spinZ >> beamType >> Ebeam >> nucleonID >> processID >> weight;
      if(!hs) continue;
    }

    Part pe{}, pp{}, pkp{}, pkm{};
    bool gotE=false, gotP=false, gotKp=false, gotKm=false;

    for(int i=0;i<Npart;i++){
      std::string pline; if(!std::getline(in, pline)) break;
      std::istringstream ps(pline);
      int idx=0, type=0, pdg=0, parent=0, fd=0; double lifetime=0;
      double px=0, py=0, pz=0, E=0, m=0, vx=0, vy=0, vz=0;
      ps >> idx >> lifetime >> type >> pdg >> parent >> fd >> px >> py >> pz >> E >> m >> vx >> vy >> vz;
      if(!ps) continue;

      if(pdg==11)       { pe  = Part{pdg,px,py,pz,E,m}; gotE=true; }
      else if(pdg==2212){ pp  = Part{pdg,px,py,pz,E,m}; gotP=true; }
      else if(pdg==321) { pkp = Part{pdg,px,py,pz,E,m}; gotKp=true; }
      else if(pdg==-321){ pkm = Part{pdg,px,py,pz,E,m}; gotKm=true; }
    }

    if(!(gotE && gotP)) continue;

    // helper: (p,theta,phi) in radians + degrees
    auto to_ptp = [](const Part& pr){
      const double p  = std::sqrt(pr.px*pr.px + pr.py*pr.py + pr.pz*pr.pz);
      const double c  = (p>0)? pr.pz/p : 1.0;
      const double th = std::acos(std::max(-1.0,std::min(1.0,c)));
      const double ph = std::atan2(pr.py, pr.px);
      return std::make_tuple(p, th, ph);
    };

    double pe_p, pe_th, pe_ph; std::tie(pe_p,pe_th,pe_ph) = to_ptp(pe);
    double pp_p, pp_th, pp_ph; std::tie(pp_p,pp_th,pp_ph) = to_ptp(pp);
    double kp_p=0, kp_th=0, kp_ph=0;
    double km_p=0, km_th=0, km_ph=0;

    if(gotKp) std::tie(kp_p,kp_th,kp_ph) = to_ptp(pkp);
    if(gotKm) std::tie(km_p,km_th,km_ph) = to_ptp(pkm);

    // Build kinematics object (phi channel if both kaons exist; fallback DVCS-like)
    // DISANAMath(e_in_E, e_out_p, e_out_th, e_out_ph, p_out_p, p_out_th, p_out_ph,
    //            kMinus_p, kMinus_th, kMinus_ph, kPlus_p, kPlus_th, kPlus_ph)
    DISANAMath kin;
    if(gotKp && gotKm){
      kin = DISANAMath(Ebeam, pe_p, pe_th, pe_ph, pp_p, pp_th, pp_ph,
                       km_p, km_th, km_ph, kp_p, kp_th, kp_ph);
    } else {
      // If your header has the DVCS-like ctor with the 'gamma' (photon/proxy), use it:
      // DISANAMath(Ebeam, pe_p, pe_th, pe_ph, pp_p, pp_th, pp_ph, g_p, g_th, g_ph, /*IsMissing*/true)
      // As a simple, safe proxy here, reuse the proton direction for the "hadron":
      kin = DISANAMath(Ebeam, pe_p, pe_th, pe_ph, pp_p, pp_th, pp_ph,
                       pp_p, pp_th, pp_ph, /*IsMissing*/ true);
    }

    // Fill branch variables (angles stored in degrees for convenience)
    auto fill_deg = [](double th_rad, double ph_rad, double &th_deg, double &ph_deg){
      th_deg = rad2deg(th_rad);
      double phd = rad2deg(ph_rad);
      if (phd>180) phd-=360; if (phd<-180) phd+=360;
      ph_deg = phd;
    };

    p_e = pe_p; fill_deg(pe_th, pe_ph, th_e, ph_e);
    p_p = pp_p; fill_deg(pp_th, pp_ph, th_p, ph_p);

    if(gotKp){ p_kp=kp_p; fill_deg(kp_th, kp_ph, th_kp, ph_kp); } else { p_kp=-999; th_kp=-999; ph_kp=-999; }
    if(gotKm){ p_km=km_p; fill_deg(km_th, km_ph, th_km, ph_km); } else { p_km=-999; th_km=-999; ph_km=-999; }

    Q2   = kin.GetQ2();
    xB   = kin.GetxB();
    W    = kin.GetW();
    Tabs = kin.GetT();   // header code defines GetT() ≡ |t| (i.e. -t >= 0)
    Phi  = kin.GetPhi(); // 0..360 deg, lepton–hadron plane
    Wgt  = weight;

    tree->Fill();
  }

  return tree;
}

void LundPlot_RDF(const char* fname="epKK.lund", double lumi=1.0){
  // style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetTitleFont(42, "");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetTickLength(0.02, "XYZ");

  // 1) Build a TTree
  std::unique_ptr<TTree> tree(LundToTree(fname));

  // 2) DataFrame over the tree
  ROOT::RDataFrame df(*tree);

  // 3) Book all histograms from df (weights via "w" when needed)

  // (A) species p/θ/φ
  auto hP_e   = df.Filter("p_e>0").Histo1D({"hP_e","e' momentum; p [GeV/c]; Entries",120,0,12}, "p_e");
  auto hTh_e  = df.Filter("th_e>-990").Histo1D({"hTh_e","e' polar angle; #theta [deg]; Entries",180,0,180}, "th_e");
  auto hPh_e  = df.Filter("ph_e>-990").Histo1D({"hPh_e","e' azimuth; #phi [deg]; Entries",180,-180,180}, "ph_e");

  auto hP_p   = df.Filter("p_p>0").Histo1D({"hP_p","p momentum; p [GeV/c]; Entries",120,0,12}, "p_p");
  auto hTh_p  = df.Filter("th_p>-990").Histo1D({"hTh_p","p polar angle; #theta [deg]; Entries",180,0,180}, "th_p");
  auto hPh_p  = df.Filter("ph_p>-990").Histo1D({"hPh_p","p azimuth; #phi [deg]; Entries",180,-180,180}, "ph_p");

  auto hP_kp  = df.Filter("p_kp>0").Histo1D({"hP_kp","K^{+} momentum; p [GeV/c]; Entries",120,0,12}, "p_kp");
  auto hTh_kp = df.Filter("th_kp>-990").Histo1D({"hTh_kp","K^{+} polar angle; #theta [deg]; Entries",180,0,180}, "th_kp");
  auto hPh_kp = df.Filter("ph_kp>-990").Histo1D({"hPh_kp","K^{+} azimuth; #phi [deg]; Entries",180,-180,180}, "ph_kp");

  auto hP_km  = df.Filter("p_km>0").Histo1D({"hP_km","K^{-} momentum; p [GeV/c]; Entries",120,0,12}, "p_km");
  auto hTh_km = df.Filter("th_km>-990").Histo1D({"hTh_km","K^{-} polar angle; #theta [deg]; Entries",180,0,180}, "th_km");
  auto hPh_km = df.Filter("ph_km>-990").Histo1D({"hPh_km","K^{-} azimuth; #phi [deg]; Entries",180,-180,180}, "ph_km");

  // (B) DVEP kinematics
  auto hQ2    = df.Histo1D({"hQ2","Q^{2}; Q^{2} [GeV^{2}]; Entries",120,0,12}, "Q2");
  auto hTneg  = df.Histo1D({"hT","-t; -t [GeV^{2}]; Entries",120,0,10.0}, "t");
  auto hXB    = df.Histo1D({"hXB","x_{B}; x_{B}; Entries",120,0,1.0}, "xB");
  auto hW     = df.Histo1D({"hW","W; W [GeV]; Entries",120,1.5,10.0}, "W");
  auto hPhi   = df.Histo1D({"hPhi","Lepton–hadron #phi; #phi [deg]; Entries",72,0,360}, "phi");
  auto hQ2vT  = df.Histo2D({"hQ2vT","Q^{2} vs -t; -t [GeV^{2}]; Q^{2} [GeV^{2}]",100,0,10.5,100,0,12}, "t", "Q2");

  // (C) weighted -t spectrum (arb. normalization)
  auto hDsDt  = df.Histo1D({"hDsDt","Weighted spectrum vs -t; -t [GeV^{2}]; arb. units",60,0,7.5}, "t", "w");

  // Force creation now (so we can style)
  auto *hP_eH  = (TH1D*)hP_e.GetPtr();   auto *hTh_eH = (TH1D*)hTh_e.GetPtr(); auto *hPh_eH=(TH1D*)hPh_e.GetPtr();
  auto *hP_pH  = (TH1D*)hP_p.GetPtr();   auto *hTh_pH = (TH1D*)hTh_p.GetPtr(); auto *hPh_pH=(TH1D*)hPh_p.GetPtr();
  auto *hP_kpH = (TH1D*)hP_kp.GetPtr();  auto *hTh_kpH= (TH1D*)hTh_kp.GetPtr();auto *hPh_kpH=(TH1D*)hPh_kp.GetPtr();
  auto *hP_kmH = (TH1D*)hP_km.GetPtr();  auto *hTh_kmH= (TH1D*)hTh_km.GetPtr();auto *hPh_kmH=(TH1D*)hPh_km.GetPtr();
  auto *hQ2H   = (TH1D*)hQ2.GetPtr();    auto *hTnegH = (TH1D*)hTneg.GetPtr();
  auto *hXBH   = (TH1D*)hXB.GetPtr();    auto *hWH    = (TH1D*)hW.GetPtr();
  auto *hPhiH  = (TH1D*)hPhi.GetPtr();   auto *hDsDtH = (TH1D*)hDsDt.GetPtr();
  auto *hQ2vTH = (TH2D*)hQ2vT.GetPtr();

  // Style
  StyleHist(hP_eH,  ColElectron()); StyleHist(hTh_eH, ColElectron()); StyleHist(hPh_eH, ColElectron());
  StyleHist(hP_pH,  ColProton());   StyleHist(hTh_pH, ColProton());   StyleHist(hPh_pH, ColProton());
  StyleHist(hP_kpH, ColKplus());    StyleHist(hTh_kpH,ColKplus());    StyleHist(hPh_kpH,ColKplus());
  StyleHist(hP_kmH, ColKminus());   StyleHist(hTh_kmH,ColKminus());   StyleHist(hPh_kmH,ColKminus());

  // ---- (A) canvas: species
  TCanvas* cA = new TCanvas("c_species", "Species distributions", 1800, 1300);
  cA->Divide(3,4, 0.005, 0.005);
  for (int i=1;i<=12;i++){ cA->cd(i)->SetGrid(1,1); cA->cd(i)->SetTicks(1,1);
    cA->cd(i)->SetLeftMargin(0.12); cA->cd(i)->SetRightMargin(0.04);
    cA->cd(i)->SetTopMargin(0.06);  cA->cd(i)->SetBottomMargin(0.12); }
  cA->cd(1);  hP_eH->Draw("HIST");   cA->cd(2);  hTh_eH->Draw("HIST");   cA->cd(3);  hPh_eH->Draw("HIST");
  cA->cd(4);  hP_pH->Draw("HIST");   cA->cd(5);  hTh_pH->Draw("HIST");   cA->cd(6);  hPh_pH->Draw("HIST");
  cA->cd(7);  hP_kpH->Draw("HIST");  cA->cd(8);  hTh_kpH->Draw("HIST");  cA->cd(9);  hPh_kpH->Draw("HIST");
  cA->cd(10); hP_kmH->Draw("HIST");  cA->cd(11); hTh_kmH->Draw("HIST");  cA->cd(12); hPh_kmH->Draw("HIST");
  cA->SaveAs("lund_species.pdf"); cA->SaveAs("lund_species.png");

  // ---- (B) canvas: DVEP
  TCanvas* cB = new TCanvas("c_dvep", "DVEP kinematics", 1600, 900);
  cB->Divide(3,2, 0.01, 0.01);
  for (int i=1;i<=6;i++){ cB->cd(i)->SetGrid(1,1); cB->cd(i)->SetTicks(1,1);
    cB->cd(i)->SetLeftMargin(0.12); cB->cd(i)->SetRightMargin(0.05);
    cB->cd(i)->SetTopMargin(0.06);  cB->cd(i)->SetBottomMargin(0.12); }
  cB->cd(1); hQ2H->SetLineWidth(2);   hQ2H->Draw("HIST");
  cB->cd(2); hTnegH->SetLineWidth(2); hTnegH->Draw("HIST");
  cB->cd(3); hXBH->SetLineWidth(2);   hXBH->Draw("HIST");
  cB->cd(4); hWH->SetLineWidth(2);    hWH->Draw("HIST");
  cB->cd(5); gPad->SetRightMargin(0.12); hQ2vTH->SetContour(80); hQ2vTH->Draw("COLZ");
  cB->cd(6); hPhiH->SetLineWidth(2);  hPhiH->Draw("HIST");
  cB->SaveAs("lund_dvep.pdf"); cB->SaveAs("lund_dvep.png");

  // ---- (C) canvas: weighted dσ/dt proxy (arb. norm; weights from LUND header)
  TCanvas* cC = new TCanvas("c_xs", "Weighted d#sigma/dt (arb. units)", 900, 650);
  cC->SetGrid(1,1); cC->SetTicks(1,1);
  cC->SetLeftMargin(0.12); cC->SetRightMargin(0.04);
  cC->SetTopMargin(0.06);  cC->SetBottomMargin(0.12);
  hDsDtH->SetLineColor(TColor::GetColor("#444444"));
  hDsDtH->SetLineWidth(2);
  hDsDtH->SetMarkerStyle(20);
  hDsDtH->SetMarkerSize(0.9);
  hDsDtH->SetMarkerColor(hDsDtH->GetLineColor());
  hDsDtH->Draw("HIST E1");
  TLatex l; l.SetTextFont(42); l.SetTextSize(0.04);
  l.DrawLatexNDC(0.16,0.93,"Weighted d#sigma/dt (header weights; arb. units)");
  cC->SaveAs("lund_dsigma_dt.pdf"); cC->SaveAs("lund_dsigma_dt.png");

  // --------------------------------------------------------------------------
  // Optional: Use your DISANAMath cross-section helpers directly on the df
  // --------------------------------------------------------------------------
  // BinManager bins;             // default bin edges (edit in DISPhiMath.h)
  // DISANAMath math;
  //
  // // 3D dσ/dφ(Q2, t, xB): expects df to have columns "Q2","t","xB","phi"
  // auto ds_dphi = math.ComputeDVCS_CrossSection(df, bins, /*luminosity=*/lumi);
  //
  // // dσ/dt in Q2 bins: expects "Q2","t"
  // auto ds_dt   = math.ComputePhi_CrossSection(df, bins, /*luminosity=*/lumi);
  //
  // // You can draw ds_dt[iq] etc. here, or serialize to a ROOT file for later.
  // // Example:
  // // TFile fout("xs_out.root","RECREATE");
  // // for (auto *h : ds_dt) if(h) h->Write();
  // // fout.Close();

  std::cout << "Wrote: lund_species.*, lund_dvep.*, lund_dsigma_dt.*" << std::endl;
}

// Entry point convenience: root -l -q 'LundPlot_RDF.C("epKK.lund")'
void LundPlot(const char* fname="epKK.lund", double lumi=1.0){
  LundPlot_RDF(fname, lumi);
}
