int plot_proton()
{
  TFile *f1 = new TFile("../result-photo-psi2S/p_solid_photo_17.2GeV.root","READ");
  TFile *f2 = new TFile("../result-electro-psi2S/p_solid_electro_17.2GeV.root","READ");

  TTree *t1 = (TTree*)f1->Get("tree");
  TTree *t2 = (TTree*)f2->Get("tree");

  double Emin = 13;
  double Emax = 13.5;
  
  TH2F *h1 = new TH2F("h1","#psi(2S) photoproduction from 17.2GeV beam; Proton Theta [deg]; Proton Momentum [GeV]",100,0,35,100,0,10);
  TH2F *h2 = new TH2F("h2","#psi(2S) electroproduction from 17.2GeV beam; Proton Theta [deg]; Proton Momentum [GeV]",100,0,35,100,0,10);

  t1->Draw("pOut->P():pOut->Theta()*180/3.14159265>>h1",Form("Eg_true>%f&&Eg_true<%f",Emin,Emax),"goff");
  t2->Draw("pOut->P():pOut->Theta()*180/3.14159265>>h2",Form("Eg_true>%f&&Eg_true<%f",Emin,Emax),"goff");

  TLatex latex;
  latex.SetTextSize(0.05);

  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);
  h1->Draw("colz");
  latex.DrawLatexNDC(0.5,0.8,Form("%.2f < E_{#gamma} < %.2f [GeV]",Emin,Emax));
  TLine l1;
  l1.SetLineStyle(9);
  l1.SetLineWidth(2);
  l1.SetLineColor(kRed);
  l1.DrawLine(18,0,18,2);  l1.DrawLine(28,0,28,2);  l1.DrawLine(18,2,28,2);  l1.DrawLine(18,0,28,0);
  l1.DrawLine(8,0,8,4.5);  l1.DrawLine(18,0,18,4.5);  l1.DrawLine(8,4.5,18,4.5);  l1.DrawLine(8,0,18,0);

  TCanvas *c2 = new TCanvas();
  gStyle->SetOptStat(0);
  h2->Draw("colz");
  l1.DrawLine(18,0,18,2);  l1.DrawLine(28,0,28,2);  l1.DrawLine(18,2,28,2);  l1.DrawLine(18,0,28,0);
  l1.DrawLine(8,0,8,4.5);  l1.DrawLine(18,0,18,4.5);  l1.DrawLine(8,4.5,18,4.5);  l1.DrawLine(8,0,18,0);
  latex.DrawLatexNDC(0.5,0.8,Form("%.2f < E_{#gamma} < %.2f [GeV]",Emin,Emax));
  
  return 0;
}
