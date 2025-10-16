#ifndef DRAWSTYLE_H
#define DRAWSTYLE_H

#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TGaxis.h>

class DrawStyle {
 public:
  DrawStyle(double titleSize = 0.05, double labelSize = 0.04,
            double xTitleOffset = 1.1, double yTitleOffset = 1.6,
            int font = 42, int maxDigits = 5, int nDivisions = 510,
            double leftMargin = 0.16, double rightMargin = 0.07,
            double bottomMargin = 0.13, double topMargin = 0.06)
      : titleSize_(titleSize),
        labelSize_(labelSize),
        xTitleOffset_(xTitleOffset),
        yTitleOffset_(yTitleOffset),
        font_(font),
        maxDigits_(maxDigits),
        nDivisions_(nDivisions),
        leftMargin_(leftMargin),
        rightMargin_(rightMargin),
        bottomMargin_(bottomMargin),
        topMargin_(topMargin) {}

  void ApplyGlobalStyle() const {
    gStyle->SetTitleSize(titleSize_, "XYZ");
    gStyle->SetLabelSize(labelSize_, "XYZ");
    gStyle->SetTitleFont(font_, "XYZ");
    gStyle->SetLabelFont(font_, "XYZ");
    gStyle->SetTitleOffset(xTitleOffset_, "X");
    gStyle->SetTitleOffset(yTitleOffset_, "Y");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPalette(kRainBow);

    TGaxis::SetMaxDigits(maxDigits_);
  }

  void NormalizeHistogram(TH1* hist) const {
    if (!hist) return;
    double integral = hist->Integral();
    if (integral > 0) hist->Scale(1.0 / integral);
  }

  void StylePad(TPad* pad,
                double left = -1, double right = -1,
                double bottom = -1, double top = -1) const {
    if (!pad) return;
    pad->SetLeftMargin(left >= 0 ? left : leftMargin_);
    pad->SetRightMargin(right >= 0 ? right : rightMargin_);
    pad->SetBottomMargin(bottom >= 0 ? bottom : bottomMargin_);
    pad->SetTopMargin(top >= 0 ? top : topMargin_);
    pad->SetTicks(1, 1);
    pad->SetFillStyle(4000);
  }

  void StyleTH1(TH1* hist,
              const char* title = "",
              double titleSize = -1, double labelSize = -1,
              double xOffset = -1, double yOffset = -1) const {
  if (!hist) return;

  hist->SetStats(0);
  hist->SetLineWidth(2);
  hist->SetTitle(title);
  hist->SetTitleFont(font_, "");
  hist->SetTitleSize(titleSize > 0 ? titleSize : titleSize_, "");

  hist->GetXaxis()->SetTitleSize(titleSize > 0 ? titleSize : titleSize_);
  hist->GetXaxis()->SetLabelSize(labelSize > 0 ? labelSize : labelSize_);
  hist->GetXaxis()->SetTitleOffset(xOffset > 0 ? xOffset : xTitleOffset_);
  hist->GetXaxis()->SetTitleFont(font_);
  hist->GetXaxis()->SetLabelFont(font_);
  hist->GetXaxis()->SetNdivisions(nDivisions_);

  hist->GetYaxis()->SetTitleSize(titleSize > 0 ? titleSize : titleSize_);
  hist->GetYaxis()->SetLabelSize(labelSize > 0 ? labelSize : labelSize_);
  hist->GetYaxis()->SetTitleOffset(yOffset > 0 ? yOffset : yTitleOffset_);
  hist->GetYaxis()->SetTitleFont(font_);
  hist->GetYaxis()->SetLabelFont(font_);
  hist->GetYaxis()->SetNdivisions(nDivisions_);
}


 void StyleTH2(TH2* hist,
              const char* title = "",
              double titleSize = -1, double labelSize = -1,
              double xOffset = -1, double yOffset = -1,
              double zLabelSize = -1.0) const {
  if (!hist) return;

  hist->SetStats(0);
  hist->SetTitle(title);
  hist->SetTitleFont(font_, "");
  hist->SetTitleSize(titleSize > 0 ? titleSize : titleSize_, "");

  hist->GetXaxis()->SetTitleSize(titleSize > 0 ? titleSize : titleSize_);
  hist->GetXaxis()->SetLabelSize(labelSize > 0 ? labelSize : labelSize_);
  hist->GetXaxis()->SetTitleOffset(xOffset > 0 ? xOffset : xTitleOffset_);
  hist->GetXaxis()->SetTitleFont(font_);
  hist->GetXaxis()->SetLabelFont(font_);
  hist->GetXaxis()->SetNdivisions(nDivisions_);

  hist->GetYaxis()->SetTitleSize(titleSize > 0 ? titleSize : titleSize_);
  hist->GetYaxis()->SetLabelSize(labelSize > 0 ? labelSize : labelSize_);
  hist->GetYaxis()->SetTitleOffset(yOffset > 0 ? yOffset : yTitleOffset_);
  hist->GetYaxis()->SetTitleFont(font_);
  hist->GetYaxis()->SetLabelFont(font_);
  hist->GetYaxis()->SetNdivisions(nDivisions_);

  hist->GetZaxis()->SetLabelSize(zLabelSize > 0 ? zLabelSize : (labelSize > 0 ? labelSize : labelSize_));
  hist->GetZaxis()->SetTitleFont(font_);
  hist->GetZaxis()->SetNdivisions(nDivisions_);
}


 private:
  double titleSize_;
  double labelSize_;
  double xTitleOffset_;
  double yTitleOffset_;
  int font_;
  int maxDigits_;
  int nDivisions_;
  double leftMargin_, rightMargin_, bottomMargin_, topMargin_;
};

#endif  // DRAWSTYLE_H
