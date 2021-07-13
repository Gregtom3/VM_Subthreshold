int test2()
{
  double beamE = 17.2;
  double Emin = 13.0;
  double Emax = 13.5;
  
  TFile *fIn = new TFile(Form("../result-photo-psi2S/p_solid_photo_%.1fGeV.root",beamE),"READ");
  TTree *tIn = (TTree*)fIn->Get("tree");
  const int N = 25;
  const double tmin = 0;
  const double tmax = 8;
  TH1F *h = new TH1F("h","; |t| [GeV^{2}]; d#sigma/dt^{(#gamma+p)} [nb / GeV^{2}]",N,tmin,tmax);
  TH1F *hCounts = new TH1F("hCounts","; |t| [GeV^{2}]; Events Detected",N,tmin,tmax);
  TH1F *hEvents = new TH1F("hEvents","; |t| [GeV^{2}]; Events",N,tmin,tmax);
  TH1F *hEvents_Acc = new TH1F("hEvents_Acc","; |t| [GeV^{2}]; Events Detected",N,tmin,tmax);
  TGraphErrors *tge = new TGraphErrors(N);
  tge->SetLineWidth(2);  tge->SetLineColor(kBlue); tge->SetMarkerStyle(20); tge->SetMarkerColor(kBlue);
  tge->SetTitle("; |t| [GeV^{2}] ; d#sigma/dt^{#gamma+p} [nb / GeV^{2}]");

  tIn->Draw("-event.t>>h",Form("weight_dodt * (Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");
  tIn->Draw("-event.t>>hCounts",Form("(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
  tIn->Draw("-event.t>>hEvents",Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
  tIn->Draw("-event.t>>hEvents_Acc",Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)*is_accept_pOut*is_accept_ePlusOut*is_accept_eMinusOut",Emin,Emax),"goff");  

  for(int i = 1; i <= h->GetNbinsX() ; i++)
    {
      if(hEvents_Acc->GetBinContent(i)<1) continue;
      double x = h->GetBinCenter(i);
      double errx = 0.5 * h->GetBinWidth(i);
      
      double total_events = hEvents->GetBinContent(i);
      double accepted_events = hEvents_Acc->GetBinContent(i);
      double acceptance = accepted_events/total_events;
      
      double y = h->GetBinContent(i) / hCounts->GetBinContent(i);
      double erry = sqrt(pow(y/sqrt(accepted_events),2)+pow(0.1*y,2));

      tge->SetPoint(i,x,y);
      tge->SetPointError(i,errx,erry);

      hEvents_Acc->SetBinError(i,sqrt(hEvents_Acc->GetBinContent(i)));
      hEvents->SetBinError(i,sqrt(hEvents->GetBinContent(i)));
    }

  // ************************************************************* //
  
  TStyle *myStyle = new TStyle("Style","My Style");
  myStyle->SetOptStat(0);
  myStyle->SetNdivisions(505);
  myStyle->SetHistLineWidth(2);
  myStyle->SetPadGridY(true);
  myStyle->SetOptLogy();
  myStyle->SetTitleXSize(0.045); 
  myStyle->SetTitleYSize(0.045);
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleYOffset(0.92);
  myStyle->SetPadRightMargin(0.05);
  myStyle->SetHistLineWidth(2);
  gROOT->SetStyle("Style");

  // ************************************************************* //
  
  TCanvas *c = new TCanvas("c","c",1000,500);
  c->Divide(2,1);
  TLatex latex;
  latex.SetTextSize(0.03);
  latex.SetTextFont(42);
  // ----------------- Pad 1 ------------------ //
  c->cd(1);
  c->SetHighLightColor(0);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  tge->Draw("AP");
  tge->GetYaxis()->SetRangeUser(0.00001,10);
  tge->GetYaxis()->SetTitleOffset(1.2);
  tge->GetXaxis()->SetLimits(0,8);
  latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
  latex.DrawLatexNDC(0.18,0.82,"50 days at L=1.2e37 cm^{-2}s^{-1}");
  latex.DrawLatexNDC(0.18,0.78,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
  // ----------------- Pad 2 ------------------ //
  c->cd(2);
  c->SetHighLightColor(0);
  gPad->SetLogy();
  TLegend * legend = new TLegend(0.72,0.69,0.91,0.77);
  hEvents_Acc->SetLineColor(kYellow+2);
  hEvents_Acc->SetFillColor(kYellow+2);
  hEvents_Acc->SetFillStyle(3001);
  hEvents_Acc->SetLineWidth(2);
  hEvents->SetLineColor(kBlack);
  hEvents->SetFillColor(kGray+1);
  hEvents->SetFillStyle(3001);
  hEvents->Draw("hist");
  hEvents->Draw("E1 same");
  hEvents->GetXaxis()->SetLimits(0,8);
  hEvents->GetYaxis()->SetRangeUser(.1,50000);
  hEvents_Acc->Draw("hist same");
  hEvents_Acc->Draw("E1 same");
  legend->AddEntry(hEvents,"Total","f");
  legend->AddEntry(hEvents_Acc,"Accepted","f");
  legend->SetTextFont(42);
  legend->Draw("same");
  latex.DrawLatexNDC(0.18,0.855,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
  latex.DrawLatexNDC(0.18,0.815,"50 days at L=1.2e37 cm^{-2}s^{-1}");
  latex.DrawLatexNDC(0.18,0.775,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
  return 0;
}
