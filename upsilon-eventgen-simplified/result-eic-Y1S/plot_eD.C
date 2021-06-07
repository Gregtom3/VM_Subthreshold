int plot_eD()
{
  int lumi = 100;
  double nuclear_factor = 0.5; // 1/2 for deuteron, 1 for proton

  TFile *f = new TFile("eD.root","READ");
  TTree *t = (TTree*)f->Get("tree");

  TH1D *h_phase_space = (TH1D*)f->Get("h_phase_space");
  double phase_space_volume = h_phase_space->GetBinContent(1);

  const int nbins = 10;
  const double Emin = 50;
  const double Emax = 75;
  const double Estep = (Emax - Emin) / ( 1.0 * nbins );
  TH1F *h_Egamma = new TH1F("h_Egamma","",nbins,Emin,Emax);
  TH1F *h_Egamma_counts = new TH1F("h_Egamma_counts","",nbins,Emin,Emax);
  TGraphErrors *tge_Egamma = new TGraphErrors(nbins);
  tge_Egamma->SetTitle(Form("Upsilon Photoproduction;E_{gamma} [GeV];Counts per %dfb^{-1}",lumi));
  tge_Egamma->SetMarkerStyle(21);
  tge_Egamma->SetMarkerColor(kBlue);
  tge_Egamma->SetLineColor(kBlue);
  tge_Egamma->SetLineWidth(2);

  t->Draw("gamma_boost.E()>>h_Egamma","weight","goff");
  t->Draw("gamma_boost.E()>>h_Egamma_counts","","goff");
 
  h_Egamma->Scale(phase_space_volume);

  for(int i = 1 ; i <= nbins ; i++)
    {
      double Emid = Emin + (i+0.5)*Estep;
      double entries = h_Egamma_counts->GetBinContent(i);
      if(entries<=0.0) continue;
      double val = h_Egamma->GetBinContent(i);
      double xsec = val/entries * pow(10,6); // fb
      double nevents = xsec*lumi*nuclear_factor*Estep;
      tge_Egamma->SetPoint(i,Emid,nevents);
      tge_Egamma->SetPointError(i,Estep/2,sqrt(nevents));
    }
  
  TCanvas *c = new TCanvas("c","c",800,600);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  tge_Egamma->GetXaxis()->SetRangeUser(Emin,Emax);
  tge_Egamma->GetYaxis()->SetRangeUser(0.01,10);
  tge_Egamma->Draw("AP");
  
  return 0;
}
