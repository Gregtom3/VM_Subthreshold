int plot_nevents_vs_Egamma()
{
  int A = 1;
  //  int type = 0; // 0 --> photo , 1 --> electro
  double eIn_E = 17.0;
  
  // ************************************************************* //

  TString stringA = (A==1) ? "p" : "D";
  //  TString stringtype = (type==0) ? "photo" : "electro";
  TFile *fIn1 = new TFile(Form("../result-photo-psi2S/%s_solid_photo_%.1fGeV.root",stringA.Data(),eIn_E),"READ");
  TTree *tIn1 = (TTree*)fIn1->Get("tree"); 

  TFile *fIn2 = new TFile(Form("../result-electro-psi2S/%s_solid_electro_%.1fGeV.root",stringA.Data(),eIn_E),"READ");
  TTree *tIn2 = (TTree*)fIn2->Get("tree"); 

  // ************************************************************* //
  
  TStyle *myStyle = new TStyle("Style","My Style");
  myStyle->SetOptStat(0);
  myStyle->SetNdivisions(505);
  myStyle->SetHistLineColor(kRed);
  myStyle->SetHistLineWidth(2);
  myStyle->SetPadGridY(true);
  myStyle->SetOptLogy();
  myStyle->SetTitleXSize(0.045); 
  myStyle->SetTitleYSize(0.045);
  myStyle->SetPadBottomMargin(0.15);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  gROOT->SetStyle("Style");

  // ************************************************************* //
  
  const int nbins_x = 100;
  const double xmin = 9.5;
  const double xmax = eIn_E+0.5;
  TH1F *h_1 = new TH1F("h_1",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);
  TH1F *h_2 = new TH1F("h_2",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);
  TH1F *h_3 = new TH1F("h_3",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);

  TH1F *h_sub_1 = new TH1F("h_sub_1",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);
  TH1F *h_sub_2 = new TH1F("h_sub_2",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);
  TH1F *h_sub_3 = new TH1F("h_sub_3",";E_{#gamma} [GeV];Nevents",nbins_x,xmin,xmax);

  const double histmin = 1;
  const double histmax = 9.5*pow(10,4);
  h_1->GetYaxis()->SetRangeUser(histmin,histmax);
  h_2->GetYaxis()->SetRangeUser(histmin,histmax);
  h_1->GetYaxis()->SetTitleOffset(0.95);

  h_sub_1->SetFillColor(kRed);
  h_sub_2->SetFillColor(kRed);
  h_sub_3->SetFillColor(kRed);

  h_sub_1->SetFillStyle(3001);
  h_sub_2->SetFillStyle(3001);
  h_sub_3->SetFillStyle(3001);

  // ************************************************************* //

  TLatex latex;
  latex.SetTextSize(0.032);
  latex.SetTextFont(42);

  const double text_xmin = 0.12;
  const double text_xmax = 0.7;
  const double text_ymin = 0.83;
  const double text_ymax = 0.9;

  TCanvas *c = new TCanvas("c","c",1800,600);
  c->Divide(3,1,0);

  c->cd(1);
  tIn1->Draw("gamma->E()>>h_1","weight_smear2","hist");
  latex.DrawLatexNDC(0.15,0.95,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(0.15,0.9,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #bf{Photoprodution}",stringA.Data()));
  latex.DrawLatexNDC(0.15,0.85,Form("50 days at 1.2e37 cm^{-2}s^{-1}A^{-1} #rightarrow %.0f events",h_1->Integral()));
  if(A!=1)
    {
      tIn1->Draw("gamma->E()>>h_sub_1","weight_smear2*is_sub","hist same");
      latex.DrawLatexNDC(0.15,0.8,Form("#color[2]{%.0f events subthreshold}",h_sub_1->Integral()));
    }
  gPad->RedrawAxis();

  c->cd(2);
  tIn2->Draw("gamma->E()>>h_2","weight_smear2","hist");
  latex.DrawLatexNDC(0.1,0.95,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(0.1,0.9,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #bf{Electroprodution}",stringA.Data()));
  latex.DrawLatexNDC(0.1,0.85,Form("50 days at 1.2e37 cm^{-2}s^{-1}A^{-1} #rightarrow %.0f events",h_2->Integral()));
  if(A!=1)
    {
      tIn2->Draw("gamma->E()>>h_sub_2","weight_smear2*is_sub","hist same");
      latex.DrawLatexNDC(0.1,0.8,Form("#color[2]{%.0f events subthreshold}",h_sub_2->Integral()));
    }
  gPad->RedrawAxis();
  
  c->cd(3);
  //  gPad->SetRightMargin(0.1);
  h_3->Add(h_1,h_2);
  h_3->Draw("hist");
  h_3->GetYaxis()->SetRangeUser(histmin,histmax);  
  latex.DrawLatexNDC(0.1,0.95,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(0.1,0.9,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #bf{Photo+electro prodution}",stringA.Data()));
  latex.DrawLatexNDC(0.1,0.85,Form("50 days at 1.2e37 cm^{-2}s^{-1}A^{-1} #rightarrow %.0f events",h_3->Integral()));
  h_sub_3->Add(h_sub_1,h_sub_2);
  if(A!=1)
    {
      h_sub_3->Draw("hist same");
      latex.DrawLatexNDC(0.1,0.8,Form("#color[2]{%.0f events subthreshold}",h_sub_3->Integral()));
    }
  gPad->RedrawAxis();

  return 0;
}
