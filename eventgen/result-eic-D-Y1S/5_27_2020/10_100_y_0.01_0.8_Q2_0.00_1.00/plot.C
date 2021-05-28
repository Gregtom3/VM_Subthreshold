int plot()
{

  double BR = 0.0238;
  TFile *f2 = new TFile("eD_noFermiMotion.root","OPEN");
  TTree *t2 = (TTree*)f2->Get("tree");
  const int N = t2->GetEntries();
  int plot = 2;
  // Plot 1
  if(plot==1)
    {
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      gStyle->SetOptStat(0);
      const int nbins_1 = 100;
      const double Wmin_1 = 10.0;
      const double Wmax_1 = 60.0;
      const double Wstep_1 = (Wmax_1-Wmin_1)/(1.*nbins_1);
      TH1F *h1 = new TH1F("h1","Upsilon 1S photoproduction cross section;W [GeV];#sigma^{#gamma + p --> #Upsilon + p'}[nb]",nbins_1,Wmin_1,Wmax_1);
      h1->SetLineColor(kBlack);
      h1->SetLineWidth(2);
      t2->Draw("event.W>>h1","weight_gamma_p","goff");
      gPad->SetLogy();
      h1->Scale(1.0/N/BR/Wstep_1);
      h1->Draw("hist");
    }
  

  // Plot 2
  if(plot==2)
    {
      double Wmin_2 = 12;
      double Wmax_2 = 14;
      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      gStyle->SetOptStat(0);
      const int nbins_2 = 50;
      const double tmin_2 = 0.0;
      const double tmax_2 = 2.5;
      const double tstep_2 = (tmax_2-tmin_2)/(1.*nbins_2);
      TH1F *h2 = new TH1F("h2","Upsilon 1S photoproduction diff. cross section;-t [GeV^{2}];d#sigma^{#gamma + p --> #Upsilon + p'}/dt [nb/GeV^{2}]",nbins_2,tmin_2,tmax_2);
      h2->SetLineColor(kBlack);
      h2->SetLineWidth(2);
      t2->Draw("-event.t>>h2",Form("weight_gamma_p*(event.W>%.2f&&event.W<%.2f)",Wmin_2,Wmax_2),"goff");
      h2->GetYaxis()->SetRangeUser(0.00005,1);
      gPad->SetLogy();
      h2->Scale(1.0/(h2->GetEntries()*BR*tstep_2));
      h2->GetYaxis()->SetRangeUser(0.0001,1);
      h2->Draw("hist");
      TF1 *f2 = new TF1("f2","expo",1.0,2.5);
      f2->SetLineColor(kRed);
      f2->SetLineWidth(2);
      h2->Fit(f2,"NQR");
      f2->SetRange(0,100);
      double integral = f2->Integral(0,100);
      TLegend *l2 = new TLegend(0.4,0.65,0.8,0.825);
      gStyle->SetLegendBorderSize(0);
      l2->AddEntry(h2,Form("Simulation | %.0f < W < %.0f",Wmin_2,Wmax_2),"l");
      l2->AddEntry(f2,Form("Exp. Fit | Integral = #sigma [nb] = %.4f",integral),"l");
      l2->Draw("same");
      f2->Draw("same");
    }
  return 0;
}
