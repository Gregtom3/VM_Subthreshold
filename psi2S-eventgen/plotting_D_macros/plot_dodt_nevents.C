int plot_dodt_nevents()
{

  double beamE = 17.2;
  double Emin = 13.0;
  double Emax = 13.5;
  
  TFile *fIn[2];
  fIn[0] = new TFile(Form("../result-photo-psi2S/p_solid_photo_%.1fGeV.root",beamE),"READ");
  fIn[1] = new TFile(Form("../result-electro-psi2S/p_solid_electro_%.1fGeV.root",beamE),"READ");
  TTree *tIn[2];
  tIn[0] = (TTree*)fIn[0]->Get("tree");
  tIn[1] = (TTree*)fIn[1]->Get("tree");
  const int N[2] = {25, 15}; //photo, electro

  const double tmin = 0;
  const double tmax = 8;

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

  TH1F *h[2];
  TH1F *hCounts[2];
  TH1F *hEvents[2];
  TH1F *hEvents_Acc[2];
  TGraphErrors *tge[2];
  TCanvas *c[2];
  TLegend *legend[2];

  int colors[2] = {596, 629};
  for(int m = 0 ; m < 2 ; m++)
    {
      h[m] = new TH1F(Form("h%d",m),"; |t| [GeV^{2}]; d#sigma/dt^{(#gamma+p)} [nb / GeV^{2}]",N[m],tmin,tmax);
      hCounts[m] = new TH1F(Form("hCounts%d",m),"; |t| [GeV^{2}]; Events Detected",N[m],tmin,tmax);
      hEvents[m] = new TH1F(Form("hEvents%d",m),"; |t| [GeV^{2}]; Events",N[m],tmin,tmax);
      hEvents_Acc[m] = new TH1F(Form("hEvents_Acc%d",m),"; |t| [GeV^{2}]; Events Detected",N[m],tmin,tmax);
      tge[m] = new TGraphErrors(N[m]);
      tge[m]->SetLineWidth(2);  tge[m]->SetLineColor(colors[m]); tge[m]->SetMarkerStyle(20); tge[m]->SetMarkerColor(colors[m]);
      tge[m]->SetTitle("; |t| [GeV^{2}] ; d#sigma/dt^{#gamma+p} [nb / GeV^{2}]");
      
      tIn[m]->Draw(Form("-event.t>>h%d",m),Form("weight_dodt * (Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");
      tIn[m]->Draw(Form("-event.t>>hCounts%d",m),Form("(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
      tIn[m]->Draw(Form("-event.t>>hEvents%d",m),Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
      tIn[m]->Draw(Form("-event.t>>hEvents_Acc%d",m),Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)*is_accept_pOut*is_accept_ePlusOut*is_accept_eMinusOut",Emin,Emax),"goff");  
      
      for(int i = 1; i <= h[m]->GetNbinsX() ; i++)
	{
	  if(hEvents_Acc[m]->GetBinContent(i)<1) continue;
	  double x = h[m]->GetBinCenter(i);
	  double errx = 0.5 * h[m]->GetBinWidth(i);
	  
	  double total_events = hEvents[m]->GetBinContent(i);
	  double accepted_events = hEvents_Acc[m]->GetBinContent(i);
	  double acceptance = accepted_events/total_events;
	  
	  double y = h[m]->GetBinContent(i) / hCounts[m]->GetBinContent(i);
	  double erry = sqrt(pow(y/sqrt(accepted_events),2)+pow(0.1*y,2));
	  
	  tge[m]->SetPoint(i,x,y);
	  tge[m]->SetPointError(i,errx,erry);
	  
	  hEvents_Acc[m]->SetBinError(i,sqrt(hEvents_Acc[m]->GetBinContent(i)));
	  hEvents[m]->SetBinError(i,sqrt(hEvents[m]->GetBinContent(i)));
	}
  
      // ************************************************************* //
  
      c[m] = new TCanvas(Form("c%d",m),"c",1000,500);
      c[m]->Divide(2,1);
      TLatex latex;
      latex.SetTextSize(0.03);
      latex.SetTextFont(42);
      // ----------------- Pad 1 ------------------ //
      c[m]->cd(1);
      c[m]->SetHighLightColor(0);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      gStyle->SetOptStat(0);
      tge[m]->Draw("AP");
      tge[m]->GetYaxis()->SetRangeUser(0.00001,10);
      tge[m]->GetYaxis()->SetTitleOffset(1.2);
      tge[m]->GetXaxis()->SetLimits(0,8);
      if(m==0) // photo
	latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
      else if(m==1) // electro
	latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[2]{electroproduction}} , %.1f GeV Beam",beamE));
      latex.DrawLatexNDC(0.3,0.82,"50 days at L=1.2e37 cm^{-2}s^{-1}");
      latex.DrawLatexNDC(0.3,0.78,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
      // ----------------- Pad 2 ------------------ //
      c[m]->cd(2);
      c[m]->SetHighLightColor(0);
      gPad->SetLogy();
      legend[m] = new TLegend(0.57,0.69,0.91,0.77);
      hEvents_Acc[m]->SetLineColor(colors[m]);
      hEvents_Acc[m]->SetFillColor(colors[m]);
      hEvents_Acc[m]->SetFillStyle(3001);
      hEvents_Acc[m]->SetLineWidth(2);
      hEvents[m]->SetLineColor(kBlack);
      hEvents[m]->SetFillColor(kGray+1);
      hEvents[m]->SetFillStyle(3001);
      hEvents[m]->Draw("hist");
      hEvents[m]->Draw("E1 same");
      hEvents[m]->GetXaxis()->SetLimits(0,8);
      hEvents[m]->GetYaxis()->SetRangeUser(.1,50000);
      hEvents_Acc[m]->Draw("hist same");
      hEvents_Acc[m]->Draw("E1 same");
      legend[m]->AddEntry(hEvents[m],Form("Total (%.0f events)",hEvents[m]->Integral()),"f");
      legend[m]->AddEntry(hEvents_Acc[m],Form("Accepted (%.0f events)", hEvents_Acc[m]->Integral()),"f");
      legend[m]->SetTextFont(42);
      legend[m]->Draw("same");
      if(m==0) // photo
	latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
      else if(m==1) // electro
	latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[2]{electroproduction}} , %.1f GeV Beam",beamE));
      latex.DrawLatexNDC(0.3,0.815,"50 days at L=1.2e37 cm^{-2}s^{-1}");
      latex.DrawLatexNDC(0.3,0.775,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
    }
  return 0;
}
