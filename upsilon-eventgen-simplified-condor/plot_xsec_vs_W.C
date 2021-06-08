#include "Lcore.h"
int plot_xsec_vs_W()
{
  //-------------------------------------------------
  //  Parameters for determining simulation to plot
  //-------------------------------------------------
  string nuc = "D";
  int nevents = 10000000;
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
  //  Y1S -> e+e- branching ratio
  //-------------------------------------------------  
  const double BR = 2.38e-2;
  //-------------------------------------------------
  //  Open the rootfile and TTree
  //-------------------------------------------------  
  TFile *f = new TFile(Form("%s-%d-%d_%.0fx%.0f_Y_%.2f_%.2f_Q2_%.2f_%.2f/combinedTree.root",nuc.c_str(),nevents,batches,Ebeam,Hbeam,ymin,ymax,Q2min,Q2max),"READ");
  TTree *t = (TTree*)f->Get("tree");
  //-------------------------------------------------
  //  Create the sigma vs. W plot
  //-------------------------------------------------  
  const int n_Wpoints = 20;
  const double Wmin = 10;
  const double Wmax = 100;
  const double Wstep = (Wmax - Wmin) / (1.0 * n_Wpoints);
  TGraphErrors *tge_W = new TGraphErrors(n_Wpoints);
  tge_W->SetTitle("Upsilon Photoproduction Cross Section; W [GeV]; #sigma(#gamma+p-->V+p') [nb]");
  tge_W->SetMarkerColor(kBlack);
  tge_W->SetMarkerStyle(20);
  tge_W->SetLineColor(kBlack);
  tge_W->SetLineWidth(2);
  //-------------------------------------------------
  //  Create the default dodt vs. t plot
  //-------------------------------------------------  
  const int nbins = 25;
  const double tstep = (tmax - tmin) / ( 1.0 * nbins );
  TH1F *h_dodt = new TH1F("h_dodt","",nbins,-tmax,-tmin);
  //-------------------------------------------------
  //  Fill the dodt vs. t histogram for each bin of W
  //-------------------------------------------------  


  for(int i = 0; i < n_Wpoints; i++) 
    {


      //-------------------------------------------------
      //  Clear the histogram, 
      //  find the Wmin(W1), Wmid, and Wmax(W2) for this point
      //-------------------------------------------------  
      h_dodt->Reset();
      double W1 = Wmin + i * Wstep;
      double W2 = Wmin + (i + 1) * Wstep;
      double Wmid = (W2+W1)/2.0;
      //-------------------------------------------------
      // Fill the dodt histogram within W range
      // (weight) is equal to d^3 o / dQ2 dy d(exp_b0t) x [ [exp_b0t]max - [exp_b0t]min ]
      // (factor) is equal to d^2 o / dQ2 dy
      // (weight)/(factor) = sigma(gamma+p-->V+p') = do/d(exp_b0t) x [ [exp_b0t]max - [exp_b0t]min ]
      //-------------------------------------------------  
      t->Draw("-event.t>>h_dodt",Form("weight/factor*(event.W>%f&&event.W<%f)",W1,W2),"goff");
      //-------------------------------------------------
      //  Scrap event if there are no entries (unphysical W)
      //-------------------------------------------------  
      if(h_dodt->GetEntries()<=0.0) continue;
      //-------------------------------------------------
      //  Scale by branching ratio, total histogram entries, and delta t
      //-------------------------------------------------  
      h_dodt->Scale(1.0 / (h_dodt->GetEntries()*BR*tstep));      
      //-------------------------------------------------
      //  Calculate the integral of the histogram (to extract sigma)
      //-------------------------------------------------  
      double integral = h_dodt->Integral(0,nbins,"width");
      //-------------------------------------------------
      //  Set the data point for the xsec vs. W plot
      //-------------------------------------------------  
      tge_W->SetPoint(i,Wmid,integral);
      tge_W->SetPointError(i,Wstep/2.0,0); 
      cout << " W = " << Wmid << " | Sigma = " << integral << endl;
    }

  //-------------------------------------------------
  //  Create the plot
  //-------------------------------------------------  
  TCanvas *c = new TCanvas("c","c",800,600);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  tge_W->GetXaxis()->SetRangeUser(Wmin-1,Wmax+1);
  tge_W->GetYaxis()->SetRangeUser(0.0001,10);
  tge_W->Draw("AP");
  
  //-------------------------------------------------
  //  Create the truth sigma vs. W plot
  //-------------------------------------------------  
  cout << " ---------- Truth ------------ " << endl;
  UPSILONMODEL::SetModel("v1");
  TGraph *tg = new TGraph(nbins);
  for(int i = 0 ; i < n_Wpoints ; i++)
    {
      double W1 = Wmin + i * Wstep;
      double W2 = Wmin + (i + 1) * Wstep;
      double Wmid = (W2+W1)/2.0;
      UPSILONMODEL::set_B(9.46030,0.938272,Wmid);
      double tmin = UPSILONMODEL::tmin(9.46030,0.938272,Wmid);
      double tmax = UPSILONMODEL::tmax(9.46030,0.938272,Wmid);
      double A = UPSILONMODEL::eventA;
      double B = UPSILONMODEL::eventB;
      if(A<0||B<0) continue;
      double integral = A/B * (exp(B*tmax) - exp(B*(-2.5)));
      tg->SetPoint(i, Wmid, integral);
      cout << " W = " << Wmid << " | Sigma = " << integral << endl;
    }
  tg->SetMarkerStyle(20);
  tg->SetMarkerColor(kRed);
  tg->SetLineWidth(2);
  tg->SetLineColor(kRed);
  tg->Draw("PC same");
  return 0;


}
