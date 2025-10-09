// LundPlot.C
// ROOT macro to plot momentum, theta, phi from a LUND file
// 4 rows (e, p, K+, K-), 3 cols (p, theta, phi), with improved styling.
//
// Usage:
//   root -l -q 'LundPlot_beautified.C("epKK.lund")'

#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TColor.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TLine.h>

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <map>

struct Part { int pdg=0; double px=0,py=0,pz=0,E=0,m=0; };
static double rad2deg(double r){ return r*180.0/3.14159265358979323846; }

// Nice color picks (ROOT color indices or custom)
static int ColElectron() { return TColor::GetColor("#1f77b4"); } // blue
static int ColProton()   { return TColor::GetColor("#d62728"); } // red
static int ColKplus()    { return TColor::GetColor("#2ca02c"); } // green
static int ColKminus()   { return TColor::GetColor("#9467bd"); } // purple

// Style helper for 1D hists
static void StyleHist(TH1D* h, int col){
  h->SetLineColor(col);
  h->SetLineWidth(3);
  h->SetFillColorAlpha(col, 0.15);
  h->SetTitleFont(42, "");   // overall title font
  h->SetTitleSize(0.07, ""); // overall title size
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleOffset(0.75);
  h->GetYaxis()->SetTitleOffset(0.90);
}

void LundPlot(const char* fname="epKK.lund") {
  // Global style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetTitleFont(42, ""); // Helvetica-like
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.2, "XYZ");
  gStyle->SetLabelSize(0.2, "XYZ");
  gStyle->SetTickLength(0.03, "XYZ");

  // Hist containers
  std::map<int, TH1D*> hP, hTh, hPh;

  auto mk = [&](const char* n, const char* t, int nb, double lo, double hi){
    TH1D* h = new TH1D(n, t, nb, lo, hi);
    h->Sumw2(false);
    return h;
  };

  // Binning ranges: tweakable
  hP[11]     = mk("hP_e",   "e' momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[11]    = mk("hTh_e",  "e' polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[11]    = mk("hPh_e",  "e' azimuth; #phi [deg]; Entries", 180, -180, 180);

  hP[2212]   = mk("hP_p",   "p momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[2212]  = mk("hTh_p",  "p polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[2212]  = mk("hPh_p",  "p azimuth; #phi [deg]; Entries", 180, -180, 180);

  hP[321]    = mk("hP_kp",  "K^{+} momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[321]   = mk("hTh_kp", "K^{+} polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[321]   = mk("hPh_kp", "K^{+} azimuth; #phi [deg]; Entries", 180, -180, 180);

  hP[-321]   = mk("hP_km",  "K^{-} momentum; p [GeV/c]; Entries", 120, 0, 12);
  hTh[-321]  = mk("hTh_km", "K^{-} polar angle; #theta [deg]; Entries", 180, 0, 180);
  hPh[-321]  = mk("hPh_km", "K^{-} azimuth; #phi [deg]; Entries", 180, -180, 180);

  // Open file
  std::ifstream in(fname);
  if(!in){
    std::cerr << "ERROR: cannot open " << fname << std::endl;
    return;
  }

  int lineNo = 0;
  while (true) {
    std::string header;
    if(!std::getline(in, header)) break;
    lineNo++;
    if(header.empty()) continue;

    std::istringstream hs(header);
    int Npart = 0; hs >> Npart;
    if(!hs) continue;

    for(int i=0;i<Npart;i++){
      std::string pline;
      if(!std::getline(in, pline)) break;
      std::istringstream ps(pline);

      int idx=0, type=0, pdg=0, parent=0, fd=0;
      double lifetime=0;
      double px=0, py=0, pz=0, E=0, m=0;
      double vx=0, vy=0, vz=0;

      ps >> idx >> lifetime >> type >> pdg >> parent >> fd >> px >> py >> pz >> E >> m >> vx >> vy >> vz;
      if(!ps) continue;

      if(pdg==11 || pdg==2212 || pdg==321 || pdg==-321){
        double p  = std::sqrt(px*px + py*py + pz*pz);
        double ct = (p>0) ? (pz/p) : 1.0;
        if(ct> 1) ct =  1;
        if(ct<-1) ct = -1;
        double th = rad2deg(std::acos(ct));
        double ph = rad2deg(std::atan2(py,px));
        if(ph>180) ph -= 360;
        if(ph<-180) ph += 360;

        hP[pdg]->Fill(p);
        hTh[pdg]->Fill(th);
        hPh[pdg]->Fill(ph);
      }
    }
  }

  // Apply consistent styling after filling (titles auto-size to content)
  StyleHist(hP[11],    ColElectron());
  StyleHist(hTh[11],   ColElectron());
  StyleHist(hPh[11],   ColElectron());

  StyleHist(hP[2212],  ColProton());
  StyleHist(hTh[2212], ColProton());
  StyleHist(hPh[2212], ColProton());

  StyleHist(hP[321],   ColKplus());
  StyleHist(hTh[321],  ColKplus());
  StyleHist(hPh[321],  ColKplus());

  StyleHist(hP[-321],  ColKminus());
  StyleHist(hTh[-321], ColKminus());
  StyleHist(hPh[-321], ColKminus());

  // Canvas with generous margins and tick marks on both axes
  TCanvas* c = new TCanvas("c", "ep -> e' p K^{+} K^{-}  |  LUND distributions", 1800, 1300);
  c->Divide(3,4, 0.005, 0.005); // small gaps
  for(int i=1;i<=12;i++){
    c->cd(i)->SetGrid(1,1);
    c->cd(i)->SetTicks(1,1);
    c->cd(i)->SetLeftMargin(0.12);
    c->cd(i)->SetRightMargin(0.04);
    c->cd(i)->SetTopMargin(0.06);
    c->cd(i)->SetBottomMargin(0.12);
  }

  auto annotateRow = [&](int padIndex, const char* rowTitle){
    c->cd(padIndex);
    TPaveText* pt = new TPaveText(0.12,0.83,0.95,0.93,"NDC");
    pt->SetFillColorAlpha(0,0.0); pt->SetTextFont(62); pt->SetTextSize(0.05);
    pt->SetTextAlign(13); // left-middle
    pt->AddText(rowTitle);
    pt->Draw();
  };

  // Draw by rows
  // Row 1: electron
  c->cd(1);  hP[11]->Draw("HIST");  //annotateRow(1, "electron (PDG 11)");
  c->cd(2);  hTh[11]->Draw("HIST");
  c->cd(3);  hPh[11]->Draw("HIST");

  // Row 2: proton
  c->cd(4);  hP[2212]->Draw("HIST"); //annotateRow(4, "proton (PDG 2212)");
  c->cd(5);  hTh[2212]->Draw("HIST");
  c->cd(6);  hPh[2212]->Draw("HIST");

  // Row 3: K+
  c->cd(7);  hP[321]->Draw("HIST");  //annotateRow(7, "K^{+} (PDG 321)");
  c->cd(8);  hTh[321]->Draw("HIST");
  c->cd(9);  hPh[321]->Draw("HIST");

  // Row 4: K-
  c->cd(10); hP[-321]->Draw("HIST"); //annotateRow(10, "K^{-} (PDG -321)");
  c->cd(11); hTh[-321]->Draw("HIST");
  c->cd(12); hPh[-321]->Draw("HIST");

  // Global super-title
  c->cd(1);
  TLatex latex;
  latex.SetTextFont(62);
  latex.SetTextSize(0.050);
  latex.DrawLatexNDC(0.12, 0.965, "ep #rightarrow e' p K^{+} K^{-}  —  LUND distributions");

  c->Update();
  c->SaveAs("lund_plots_beautified.pdf");
  c->SaveAs("lund_plots_beautified.png");
  std::cout << "Saved: lund_plots_beautified.pdf / .png" << std::endl;
}
