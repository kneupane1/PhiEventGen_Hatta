// LundPlot.C — uses DISANAMath 
// (A) 4x3: p, theta, phi for e, p, K+, K-
// (B) 2x3: Q^2, -t(=|t|), x_B, W, Q^2 vs -t, lepton–hadron φ
// (C) weighted spectrum vs -t (uses LUND header weight)

#include "DISPhiMath.h"   // your header that defines DISANAMath

#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <map>

struct Part { int pdg=0; double px=0,py=0,pz=0,E=0,m=0; };
static inline double rad2deg(double r){ return r*180.0/3.14159265358979323846; }
static inline double deg2rad(double d){ return d*3.14159265358979323846/180.0; }

// colors
static int ColElectron() { return TColor::GetColor("#1f77b4"); }
static int ColProton()   { return TColor::GetColor("#d62728"); }
static int ColKplus()    { return TColor::GetColor("#2ca02c"); }
static int ColKminus()   { return TColor::GetColor("#9467bd"); }

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

void LundPlot(const char* fname="epKK.lund") {
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

  // (A) species hists
  std::map<int, TH1D*> hP, hTh, hPh;
  auto mk1 = [&](const char* n, const char* t, int nb, double lo, double hi){
    TH1D* h = new TH1D(n, t, nb, lo, hi); h->Sumw2(false); return h; };
  hP[11]     = mk1("hP_e",   "e' momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[11]    = mk1("hTh_e",  "e' polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[11]    = mk1("hPh_e",  "e' azimuth; #phi [deg]; Entries", 180, -180, 180);
  hP[2212]   = mk1("hP_p",   "p momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[2212]  = mk1("hTh_p",  "p polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[2212]  = mk1("hPh_p",  "p azimuth; #phi [deg]; Entries", 180, -180, 180);
  hP[321]    = mk1("hP_kp",  "K^{+} momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[321]   = mk1("hTh_kp", "K^{+} polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[321]   = mk1("hPh_kp", "K^{+} azimuth; #phi [deg]; Entries", 180, -180, 180);
  hP[-321]   = mk1("hP_km",  "K^{-} momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[-321]  = mk1("hTh_km", "K^{-} polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[-321]  = mk1("hPh_km", "K^{-} azimuth; #phi [deg]; Entries", 180, -180, 180);

  // (B) DVEP kinematics
  TH1D* hQ2   = new TH1D("hQ2",  "Q^{2}; Q^{2} [GeV^{2}]; Entries", 120, 0, 12);
  TH1D* hTneg = new TH1D("hT",   "-t; -t [GeV^{2}]; Entries",       120, 0, 10.0);
  TH1D* hXB   = new TH1D("hXB",  "x_{B}; x_{B}; Entries",           120, 0, 1.0);
  TH1D* hW    = new TH1D("hW",   "W; W [GeV]; Entries",             120, 1.5, 10.0);
  TH2D* hQ2vT = new TH2D("hQ2vT","Q^{2} vs -t; -t [GeV^{2}]; Q^{2} [GeV^{2}]", 100, 0, 10.5, 100, 0, 12);
  TH1D* hPhi  = new TH1D("hPhi","Lepton–hadron #phi; #phi [deg]; Entries", 72, 0, 360);

  // (C) weighted dσ/dt proxy
  TH1D* hDsDt = new TH1D("hDsDt","Weighted spectrum vs -t; -t [GeV^{2}]; arb. units", 60, 0, 7.5);
  hDsDt->Sumw2();

  // read LUND
  std::ifstream in(fname);
  if (!in) { std::cerr << "ERROR: cannot open " << fname << std::endl; return; }

  while (true) {
    std::string header;
    if (!std::getline(in, header)) break;
    if (header.empty()) continue;

    // Header: Npart A Z Tpol spinZ beamType Ebeam nucleonID processID weight
    int Npart=0, A=0, Z=0, beamType=0, nucleonID=0, processID=0;
    double Tpol=0, spinZ=0, Ebeam=0, weight=1.0;
    {
      std::istringstream hs(header);
      hs >> Npart >> A >> Z >> Tpol >> spinZ >> beamType >> Ebeam >> nucleonID >> processID >> weight;
      if (!hs) continue;
    }

    Part pe{}, pp{}, pkp{}, pkm{};
    bool gotE=false, gotP=false, gotKp=false, gotKm=false;

    for (int i=0;i<Npart;i++) {
      std::string pline; if (!std::getline(in, pline)) break;
      std::istringstream ps(pline);
      int idx=0, type=0, pdg=0, parent=0, fd=0; double lifetime=0;
      double px=0, py=0, pz=0, E=0, m=0, vx=0, vy=0, vz=0;
      ps >> idx >> lifetime >> type >> pdg >> parent >> fd >> px >> py >> pz >> E >> m >> vx >> vy >> vz;
      if (!ps) continue;

      // fill species
      const double p  = std::sqrt(px*px + py*py + pz*pz);
      double ct = (p>0) ? (pz/p) : 1.0; ct = std::max(-1.0,std::min(1.0,ct));
      const double th = rad2deg(std::acos(ct));
      double ph = rad2deg(std::atan2(py,px)); if (ph>180) ph-=360; if (ph<-180) ph+=360;

      if (pdg==11 || pdg==2212 || pdg==321 || pdg==-321) {
        hP[pdg]->Fill(p);
        hTh[pdg]->Fill(th);
        hPh[pdg]->Fill(ph);
      }

      if (pdg==11)       { pe  = Part{pdg,px,py,pz,E,m}; gotE=true; }
      else if (pdg==2212){ pp  = Part{pdg,px,py,pz,E,m}; gotP=true; }
      else if (pdg==321) { pkp = Part{pdg,px,py,pz,E,m}; gotKp=true; }
      else if (pdg==-321){ pkm = Part{pdg,px,py,pz,E,m}; gotKm=true; }
    }

    if (!(gotE && gotP)) continue;

    // convert cartesian to (p,theta,phi) in RADIANS (DISANAMath expects radians)
    auto cart_to_ptp = [&](const Part& pr){
      const double p  = std::sqrt(pr.px*pr.px + pr.py*pr.py + pr.pz*pr.pz);
      const double c  = (p>0)? pr.pz/p : 1.0;
      const double th = std::acos(std::max(-1.0,std::min(1.0,c)));
      const double ph = std::atan2(pr.py, pr.px);
      return std::make_tuple(p, th, ph);
    };

    // inputs for DISANAMath ctor
    auto [pe_p,  pe_th,  pe_ph]  = cart_to_ptp(pe);
    auto [pp_p,  pp_th,  pp_ph]  = cart_to_ptp(pp);
    auto [kp_p,  kp_th,  kp_ph]  = cart_to_ptp(pkp);
    auto [km_p,  km_th,  km_ph]  = cart_to_ptp(pkm);

    // Build kinematics object (phi channel); beam along +z with energy Ebeam
    // DISANAMath(e_in_E, e_out_p, e_out_th, e_out_ph, p_out_p, p_out_th, p_out_ph,
    //            kMinus_p, kMinus_th, kMinus_ph, kPlus_p, kPlus_th, kPlus_ph)
    DISANAMath kin(Ebeam,
                   pe_p, pe_th, pe_ph,
                   pp_p, pp_th, pp_ph,
                   km_p, km_th, km_ph,
                   kp_p, kp_th, kp_ph);

    // Fill DVEP
    const double Q2   = kin.GetQ2();
    const double xB   = kin.GetxB();
    const double W    = kin.GetW();
    const double tneg = kin.GetT();         // |t| (i.e., -t) from header logic
    const double phi  = kin.GetPhi();       // 0..360 deg (lepton–hadron plane)

    hQ2->Fill(Q2);
    hXB->Fill(xB);
    hW->Fill(W);
    hTneg->Fill(tneg); 
    hQ2vT->Fill(tneg, Q2); 
    hPhi->Fill(phi);

    // weighted -t spectrum
    hDsDt->Fill(tneg, weight);
  }

  // style species
  StyleHist(hP[11],    ColElectron()); StyleHist(hTh[11],   ColElectron()); StyleHist(hPh[11],   ColElectron());
  StyleHist(hP[2212],  ColProton());   StyleHist(hTh[2212], ColProton());   StyleHist(hPh[2212], ColProton());
  StyleHist(hP[321],   ColKplus());    StyleHist(hTh[321],  ColKplus());    StyleHist(hPh[321],  ColKplus());
  StyleHist(hP[-321],  ColKminus());   StyleHist(hTh[-321], ColKminus());   StyleHist(hPh[-321], ColKminus());

  // (A) canvas
  TCanvas* cA = new TCanvas("c_species", "Species distributions", 1800, 1300);
  cA->Divide(3,4, 0.005, 0.005);
  for (int i=1;i<=12;i++){ cA->cd(i)->SetGrid(1,1); cA->cd(i)->SetTicks(1,1);
    cA->cd(i)->SetLeftMargin(0.12); cA->cd(i)->SetRightMargin(0.04);
    cA->cd(i)->SetTopMargin(0.06);  cA->cd(i)->SetBottomMargin(0.12); }
  cA->cd(1);  hP[11]->Draw("HIST");   cA->cd(2);  hTh[11]->Draw("HIST");   cA->cd(3);  hPh[11]->Draw("HIST");
  cA->cd(4);  hP[2212]->Draw("HIST"); cA->cd(5);  hTh[2212]->Draw("HIST"); cA->cd(6);  hPh[2212]->Draw("HIST");
  cA->cd(7);  hP[321]->Draw("HIST");  cA->cd(8);  hTh[321]->Draw("HIST");  cA->cd(9);  hPh[321]->Draw("HIST");
  cA->cd(10); hP[-321]->Draw("HIST"); cA->cd(11); hTh[-321]->Draw("HIST"); cA->cd(12); hPh[-321]->Draw("HIST");
  cA->SaveAs("lund_species.pdf"); cA->SaveAs("lund_species.png");

  // (B) canvas
  TCanvas* cB = new TCanvas("c_dvep", "DVEP kinematics", 1600, 900);
  cB->Divide(3,2, 0.01, 0.01);
  for (int i=1;i<=6;i++){ cB->cd(i)->SetGrid(1,1); cB->cd(i)->SetTicks(1,1);
    cB->cd(i)->SetLeftMargin(0.12); cB->cd(i)->SetRightMargin(0.05);
    cB->cd(i)->SetTopMargin(0.06);  cB->cd(i)->SetBottomMargin(0.12); }
  cB->cd(1); hQ2->SetLineWidth(1);   hQ2->Draw("HIST");
  cB->cd(2); hTneg->SetLineWidth(1); hTneg->Draw("HIST");
  cB->cd(3); hXB->SetLineWidth(1);   hXB->Draw("HIST");
  cB->cd(4); hW->SetLineWidth(1);    hW->Draw("HIST");
  cB->cd(5); gPad->SetRightMargin(0.12); hQ2vT->SetContour(80); hQ2vT->Draw("COLZ");
  cB->cd(6); hPhi->SetLineWidth(1);  hPhi->Draw("HIST"); // lepton–hadron φ (0..360)
  cB->SaveAs("lund_dvep.pdf"); cB->SaveAs("lund_dvep.png");

  // (C) canvas
  TCanvas* cC = new TCanvas("c_xs", "Weighted d#sigma/dt (arb. units)", 900, 650);
  cC->SetGrid(1,1); cC->SetTicks(1,1);
  cC->SetLeftMargin(0.12); cC->SetRightMargin(0.04);
  cC->SetTopMargin(0.06);  cC->SetBottomMargin(0.12);
  hDsDt->SetLineColor(TColor::GetColor("#444444"));
  hDsDt->SetLineWidth(1);
  hDsDt->SetMarkerStyle(20);
  hDsDt->SetMarkerSize(0.9);
  hDsDt->SetMarkerColor(hDsDt->GetLineColor());
  hDsDt->Draw("HIST E1");
  TLatex l; l.SetTextFont(42); l.SetTextSize(0.04);
  l.DrawLatexNDC(0.16,0.93,"Weighted d#sigma/dt (header weights; arb. units)");
  cC->SaveAs("lund_dsigma_dt.pdf"); cC->SaveAs("lund_dsigma_dt.png");

  std::cout << "Wrote: lund_species.*, lund_dvep.*, lund_dsigma_dt.*" << std::endl;
}
