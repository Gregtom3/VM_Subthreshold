#include "Lcore.h"
int plot_dodt()
{
  //-------------------------------------------------
  //  Parameters for determining simulation to plot
  //-------------------------------------------------
  string nuc = "D";
  int nevents = 1000000;
  int batches = 200;
  double Ebeam = 10;
  double Hbeam = 100;
  double ymin  = 0.01;
  double ymax  = 0.80;
  double Q2min = 0.00;
  double Q2max = 1.00;
  double tmin  = -2.5;
  double tmax  = 0.0;
  //-------------------------------------------------
  //  Luminosity for error bars
  //-------------------------------------------------  
  int lumi = 100;
  double nuclear_factor = 0.5; // 1/2 for deuteron, 1 for proton
  const double BR = 2.38e-2;
  //-------------------------------------------------
  //  Open the rootfile and TTree
  //-------------------------------------------------  
  TFile *f = new TFile(Form("%s-%d-%d_%.0fx%.0f_Y_%.2f_%.2f_Q2_%.2f_%.2f/combinedTree.root",nuc.c_str(),nevents,batches,Ebeam,Hbeam,ymin,ymax,Q2min,Q2max),"READ");
  TTree *t = (TTree*)f->Get("tree");
  //-------------------------------------------------
  //  Calculate the %-success rate of the simulation
  //-------------------------------------------------  
  const int tEntries = t->GetEntries();
  const double acceptance = (1.0 * tEntries)/(1.0 * nevents);
  double phase_space = 0.0;
  TH1D *h_phase_space = new TH1D("h_phase_space","",1,0,100);
  t->Draw("event.phase_space>>h_phase_space","","goff",1);
  phase_space = h_phase_space->GetBinContent(1);
  const double phase_space_acceptance = acceptance * phase_space; 
  //-------------------------------------------------
  //  Construct the histograms for filling
  //-------------------------------------------------  
  const int nbins = 25;      // # data points in final plot
  const double Wmin = 40;    
  const double Wmax = 42;
  const double tstep = (tmax - tmin) / ( 1.0 * nbins );
  TH1F *h_dodt = new TH1F("h_dodt","Upsilon Photoproduction Diff. Cross Section;-t [GeV^{2}];d#sigma/dt(#gamma+p-->V+p') [nb/GeV^{2}]",nbins,-tmax,-tmin);
  //-------------------------------------------------
  //  Fill the Histogram
  //-------------------------------------------------  
  t->Draw("-event.t>>h_dodt",Form("weight/factor*(event.W>%f&&event.W<%f)",Wmin,Wmax),"goff");
  h_dodt->Scale(1.0 / (h_dodt->GetEntries()*BR*tstep));
  //-------------------------------------------------
  //  Create the plot
  //-------------------------------------------------  
  TCanvas *c = new TCanvas("c","c",800,600);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  h_dodt->GetXaxis()->SetRangeUser(-tmax,-tmin);
  h_dodt->GetYaxis()->SetRangeUser(0.0001,10);
  h_dodt->Draw("hist");
  //-------------------------------------------------
  //  Create the truth dodt plot
  //-------------------------------------------------  
  double W = (Wmax+Wmin)/2.0;
  UPSILONMODEL::SetModel("v1");
  UPSILONMODEL::set_B(9.46030,0.938272,W);
  TGraph *tg = new TGraph(nbins);
  for(int i = 1 ; i <= nbins ; i++)
    {
      double t = tmin + (i + 0.5) * tstep;
      tg->SetPoint(i, -t, UPSILONMODEL::dSigmaY1S(W,t));
    }
  tg->SetMarkerStyle(20);
  tg->SetMarkerColor(kRed);
  tg->SetLineWidth(2);
  tg->SetLineColor(kRed);
  tg->Draw("PC same");
  return 0;


}
