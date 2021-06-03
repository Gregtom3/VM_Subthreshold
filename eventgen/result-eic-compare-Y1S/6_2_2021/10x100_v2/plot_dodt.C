int plot_dodt()
{
  TFile *f1 = new TFile("ep.root","READ");
  TFile *f2 = new TFile("eD.root","READ");
  TH1F *h1[10];
  TH1F *h2[10];
  double W1 = 10;
  double W2 = 15;
  double Wstep = 0.5;
  for(int i = 0 ; i < 10 ; i++)
    {
      h1[i]=(TH1F*)f1->Get(Form("h_W_%d",i+1));
      h2[i]=(TH1F*)f2->Get(Form("h_W_%d",i+1));
      h2[i]->GetYaxis()->SetRangeUser(0.0001,1);
      h2[i]->SetTitle(Form("%.1f < W < %.1f",W1+i*Wstep,W1+(i+1)*Wstep));
      h1[i]->SetLineColor(kBlack);
      h2[i]->SetLineColor(kBlue);
      h1[i]->SetMarkerColor(kBlack);
      h2[i]->SetMarkerColor(kBlue);
      h1[i]->SetMarkerStyle(20);
      h2[i]->SetMarkerStyle(20);
      h1[i]->SetLineWidth(2);
      h2[i]->SetLineWidth(2);
    }
  TCanvas *c[10];
  for(int j = 5 ; j < 10 ; j++)
    {
      c[j]=new TCanvas(Form("c%d",j),Form("c%d",j),800,600);
      gPad->SetLogy();
      gStyle->SetOptStat(0);
      TLegend *l1 = new TLegend(0.72,0.65,0.85,0.8);
      l1->AddEntry(h1[j],"e+p","l");
      l1->AddEntry(h2[j],"e+D","l");
      l1->SetTextSize(0.05);
      h2[j]->Draw("hist");
      if(h1[j]->GetEntries()>0)
	{
	  h1[j]->Draw("hist same");
	}
      l1->Draw("same");
    }
  return 0; 
}
