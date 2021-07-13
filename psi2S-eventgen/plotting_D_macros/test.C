int test()
{
  TFile *fIn = new TFile("../result-photo-psi2S/p_solid_photo_17.2GeV.root","READ");
  TTree *tIn = (TTree*)fIn->Get("tree");
  TH1F *h = new TH1F("h","; |t| [GeV^{2}]; d#sigma/dt^{#gamma+p} [nb / GeV^{2}]",100,0,8);
  TH1F *hCounts = new TH1F("hCounts","; |t| [GeV^{2}]; d#sigma/dt^{#gamma+p} [nb / GeV^{2}]",100,0,8);

  tIn->Draw("-event.t>>h","weight_dodt * (abs(Eg_true - 12.25) < 0.25)","goff");
  tIn->Draw("-event.t>>hCounts","(abs(Eg_true - 12.25) < 0.25)","goff");  
  for(int i = 1; i <= h->GetNbinsX() ; i++)
    {
      if(hCounts->GetBinContent(i)==0) continue;
      h->SetBinContent(i, h->GetBinContent(i) / hCounts->GetBinContent(i));
    }
  h->Draw("hist");
  return 0;
}
