#include "Lcore.h"
int upsilon_compare()
{
  const int N = 100;
  double W1 = 20.0;
  double W2 = W1+1.0;
  double W3 = W1+2.0;
  TGraph *g20 = new TGraph(N);
  TGraph *g21 = new TGraph(N);
  TGraph *g22 = new TGraph(N);
  g20->SetLineColor(kBlue);
  g21->SetLineColor(kGreen);
  g22->SetLineColor(kBlack);
  g20->SetLineWidth(2);
  g21->SetLineWidth(2);
  g22->SetLineWidth(2);


  double tmin = -2.5;
  double tmax = 0.0;
  double tstep = (tmax-tmin)/(1.0*N);
  UPSILONMODEL::SetModel("v1");


  for(int i = 0 ; i < N ; i++)
    {
      double t = tmin + i*tstep;
      g20->SetPoint(i, -t, 3.89*pow(10,5)*UPSILONMODEL::dSigmaY1S(W1,t));
      g21->SetPoint(i, -t, 3.89*pow(10,5)*UPSILONMODEL::dSigmaY1S(W2,t));
      g22->SetPoint(i, -t, 3.89*pow(10,5)*UPSILONMODEL::dSigmaY1S(W3,t));
    }

  TCanvas *c = new TCanvas("c","c",800,600);
  TLegend *l = new TLegend(0.6,0.7,0.8,0.85);
  l->AddEntry(g20,Form("W=%.0f",W1),"l");
  l->AddEntry(g21,Form("W=%.0f",W2),"l");
  l->AddEntry(g22,Form("W=%.0f",W3),"l");
  gPad->SetLogy();
  g20->GetXaxis()->SetTitle("-t [GeV^{2}]");
  g20->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
  g20->GetYaxis()->SetRangeUser(0.00005,1.05);
  g20->Draw("AC");
  g21->Draw("C same");
  g22->Draw("C same");
  l->Draw("same");
 
  /*
  TGraph *g = new TGraph(N);
  double Wmin = 14.0;
  double Wmax = 22.0;
  double Wstep = (Wmax - Wmin)/(1.0*N);
  for(int i = 0 ; i < N ; i++)
    {
      double W = Wmin + i*Wstep;
      g->SetPoint(i,W,3.89*pow(10,5)*UPSILONMODEL::dSigmaY1S(W,0.0));
    }
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  gPad->SetLogy();
  g->GetXaxis()->SetTitle("W [GeV]");
  g->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
  g->GetYaxis()->SetRangeUser(0.0001,1);
  g->Draw("AC");



  */
  /*
  TGraph *g2_Re = new TGraph(N);
  TGraph *g2_Im = new TGraph(N);
  g2_Re->SetTitle("Complex Scattering Amplitude T");
  double Wmin2 = 14.0;
  double Wmax2 = 50.0;
  double Wstep2 = (Wmax2 - Wmin2)/(1.0*N);
  for(int i = 0 ; i < N ; i++)
    {
      double s = pow(Wmin2 + i*Wstep2,2.0);
      g2_Re->SetPoint(i,sqrt(s),abs(UPSILONMODEL::RealT(s)));
      g2_Im->SetPoint(i,sqrt(s),abs(UPSILONMODEL::ImaginaryT(s)));
    }
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  //  gPad->SetLogy();
  gPad->SetLeftMargin(0.15);
  g2_Re->GetXaxis()->SetTitle("W [GeV]");
  g2_Re->GetYaxis()->SetTitle("T_{#Upsilon_{p}}");
  g2_Re->GetYaxis()->SetRangeUser(0,1400);
  g2_Re->SetLineColor(kViolet+3);
  g2_Re->SetLineWidth(2);
  g2_Im->SetLineColor(kOrange);
  g2_Im->SetLineWidth(2);
  TLegend *l2 = new TLegend(0.2,0.7,0.375,0.8);
  l2->AddEntry(g2_Re,"Re T (T_{0}=20.5)","l");
  l2->AddEntry(g2_Im,"Im T","l");
  g2_Re->Draw("AC");
  g2_Im->Draw("C same");
  l2->Draw("same");
  */
return 0;
}
