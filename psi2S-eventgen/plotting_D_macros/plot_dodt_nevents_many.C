int plot_dodt_nevents_many()
{
  gROOT->SetBatch(true);
  
  // ****************************************************** //
  int A = 1;
  double beamE = 17.2;
  double Elower = 9.5;
  double Ehigher = 17;
  // ****************************************************** //
  // Number of Egamma bins
  int nE[2] = {5,3}; // photo, electro
                       // ALWAYS ASSUME nE_photo > nE_electro
  
  // Number of data points in dsigma/dt graphs
  const int N[2] = {25, 15}; //photo, electro
  // ****************************************************** //
    
  TFile *fIn[2];
  TString nuc = "";
  if(A==1)
    nuc = "p";
  else if(A==2)
    nuc = "D";

  TString production[2];
  production[0]="photoproduction";
  production[1]="electroproduction";

  fIn[0] = new TFile(Form("../result-photo-psi2S/p_solid_photo_%.1fGeV.root",beamE),"READ");
  fIn[1] = new TFile(Form("../result-electro-psi2S/p_solid_electro_%.1fGeV.root",beamE),"READ");
  TTree *tIn[2];
  tIn[0] = (TTree*)fIn[0]->Get("tree");
  tIn[1] = (TTree*)fIn[1]->Get("tree");

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

  TH1F *h[2][nE[0]];
  TH1F *hCounts[2][nE[0]];
  TH1F *hEvents[2][nE[0]];
  TH1F *hEvents_Acc[2][nE[0]];
  TGraphErrors *tge[2][nE[0]];
  TCanvas *c[2][nE[0]];
  TLegend *legend[2][nE[0]];

  int colors[2] = {596, 629};

  
  

  for(int m = 0 ; m < 2 ; m++)
    {
      for(int M = 0 ; M < nE[m] ; M++)
	{
	  double Emin = Elower + M * (Ehigher-Elower)/(nE[m]);
	  double Emax = Elower + (M+1) * (Ehigher-Elower)/(nE[m]);
	  h[m][M] = new TH1F(Form("h%d%d",m,M),"; |t| [GeV^{2}]; d#sigma/dt^{(#gamma+p)} [nb / GeV^{2}]",N[m],tmin,tmax);
	  hCounts[m][M] = new TH1F(Form("hCounts%d%d",m,M),"; |t| [GeV^{2}]; Events Detected",N[m],tmin,tmax);
	  hEvents[m][M] = new TH1F(Form("hEvents%d%d",m,M),"; |t| [GeV^{2}]; Events",N[m],tmin,tmax);
	  hEvents_Acc[m][M] = new TH1F(Form("hEvents_Acc%d%d",m,M),"; |t| [GeV^{2}]; Events Detected",N[m],tmin,tmax);
	  tge[m][M] = new TGraphErrors(N[m]);
	  tge[m][M]->SetLineWidth(2);  tge[m][M]->SetLineColor(colors[m]); tge[m][M]->SetMarkerStyle(20); tge[m][M]->SetMarkerColor(colors[m]);
	  tge[m][M]->SetTitle("; |t| [GeV^{2}] ; d#sigma/dt^{#gamma+p} [nb / GeV^{2}]");
      
	  tIn[m]->Draw(Form("-event.t>>h%d%d",m,M),Form("weight_dodt * (Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");
	  tIn[m]->Draw(Form("-event.t>>hCounts%d%d",m,M),Form("(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
	  tIn[m]->Draw(Form("-event.t>>hEvents%d%d",m,M),Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)",Emin,Emax),"goff");  
	  tIn[m]->Draw(Form("-event.t>>hEvents_Acc%d%d",m,M),Form("weight_total*lumi*time/Nsim*(Eg_true > %f && Eg_true < %f)*is_accept_pOut*is_accept_ePlusOut*is_accept_eMinusOut",Emin,Emax),"goff");  
      
	  for(int i = 1; i <= h[m][M]->GetNbinsX() ; i++)
	    {
	      if(hEvents_Acc[m][M]->GetBinContent(i)<1) continue;
	      double x = h[m][M]->GetBinCenter(i);
	      double errx = 0.5 * h[m][M]->GetBinWidth(i);
	  
	      double total_events = hEvents[m][M]->GetBinContent(i);
	      double accepted_events = hEvents_Acc[m][M]->GetBinContent(i);
	      double acceptance = accepted_events/total_events;
	  
	      double y = h[m][M]->GetBinContent(i) / hCounts[m][M]->GetBinContent(i);
	      double erry = sqrt(pow(y/sqrt(accepted_events),2)+pow(0.1*y,2));
	  
	      tge[m][M]->SetPoint(i,x,y);
	      tge[m][M]->SetPointError(i,errx,erry);
	  
	      hEvents_Acc[m][M]->SetBinError(i,sqrt(hEvents_Acc[m][M]->GetBinContent(i)));
	      hEvents[m][M]->SetBinError(i,sqrt(hEvents[m][M]->GetBinContent(i)));
	    }
  
	  // ************************************************************* //
  
	  c[m][M] = new TCanvas(Form("c%d%d",m,M),"c",1000,500);
	  c[m][M]->Divide(2,1);
	  TLatex latex;
	  latex.SetTextSize(0.03);
	  latex.SetTextFont(42);
	  // ----------------- Pad 1 ------------------ //
	  c[m][M]->cd(1);
	  c[m][M]->SetHighLightColor(0);
	  gPad->SetLogy();
	  gPad->SetLeftMargin(0.15);
	  gStyle->SetOptStat(0);
	  tge[m][M]->Draw("AP");
	  tge[m][M]->GetYaxis()->SetRangeUser(0.00001,10);
	  tge[m][M]->GetYaxis()->SetTitleOffset(1.2);
	  tge[m][M]->GetXaxis()->SetLimits(0,8);
	  if(m==0) // photo
	    latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
	  else if(m==1) // electro
	    latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[2]{electroproduction}} , %.1f GeV Beam",beamE));
	  latex.DrawLatexNDC(0.18,0.82,"50 days at L=1.2e37 cm^{-2}s^{-1}");
	  latex.DrawLatexNDC(0.18,0.78,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
	  // ----------------- Pad 2 ------------------ //
	  c[m][M]->cd(2);
	  c[m][M]->SetHighLightColor(0);
	  gPad->SetLogy();
	  legend[m][M] = new TLegend(0.72,0.69,0.91,0.77);
	  hEvents_Acc[m][M]->SetLineColor(colors[m]);
	  hEvents_Acc[m][M]->SetFillColor(colors[m]);
	  hEvents_Acc[m][M]->SetFillStyle(3001);
	  hEvents_Acc[m][M]->SetLineWidth(2);
	  hEvents[m][M]->SetLineColor(kBlack);
	  hEvents[m][M]->SetFillColor(kGray+1);
	  hEvents[m][M]->SetFillStyle(3001);
	  hEvents[m][M]->Draw("hist");
	  hEvents[m][M]->Draw("E1 same");
	  hEvents[m][M]->GetXaxis()->SetLimits(0,8);
	  hEvents[m][M]->GetYaxis()->SetRangeUser(.1,50000);
	  hEvents_Acc[m][M]->Draw("hist same");
	  hEvents_Acc[m][M]->Draw("E1 same");
	  legend[m][M]->AddEntry(hEvents[m][M],"Total","f");
	  legend[m][M]->AddEntry(hEvents_Acc[m][M],"Accepted","f");
	  legend[m][M]->SetTextFont(42);
	  legend[m][M]->Draw("same");
	  if(m==0) // photo
	    latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[4]{photoproduction}} , %.1f GeV Beam",beamE));
	  else if(m==1) // electro
	    latex.DrawLatexNDC(0.18,0.86,Form("#bf{SoLID Simulation} #psi(2S) #bf{#color[2]{electroproduction}} , %.1f GeV Beam",beamE));
	  latex.DrawLatexNDC(0.18,0.815,"50 days at L=1.2e37 cm^{-2}s^{-1}");
	  latex.DrawLatexNDC(0.18,0.775,Form("%.2f < E_{#gamma} < %.2f",Emin,Emax));
	  c[m][M]->SaveAs(Form("PLOTS_DODT_NEVENTS_MANY/%s_%s_%.1f_GeV_beam_%.2f_Egamma_%.2f.pdf",nuc.Data(),production[m].Data(),beamE,Emin,Emax));
	}
    }
      return 0;
    }
