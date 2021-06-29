#include "../Lcore.h"
double x(double Egamma);
void tlimits_W(double W, double * t);
void tlimits_Eg(double Eg, double * t);
double psi2S_dodt_23g(double * xx, double * p)
{
  double X = x(xx[0]);
  double t = xx[1];
  double N2g = 1818.89;
  double N3g = 809.949;
  double v = 1.0 / (16.0 * M_PI);
  double R = 1.0;
  double M = 3.686097;//GeV
  double result = N2g * v * pow(1.0 - X, 2) * exp(1.13 * t) / (R * R * M * M) + N3g * v * exp(1.13 * t) / pow(R * M, 4);//nb GeV^-2
  return result;//ds/dt in unit GeV^-4
}

double psi2S_o_23g(double * xx, double * p)
{
  double X = x(xx[0]);
  double tlimits[2];
  tlimits_Eg(xx[0], tlimits);
  double tmin = tlimits[0];
  double tmax = tlimits[1];

  double N2g = 1818.89;
  double N3g = 809.949;
  double v = 1.0 / (16.0 * M_PI);
  double R = 1.0;
  double M = 3.686097;//GeV
  double result_lower_bound = N2g * v * pow(1.0 - X, 2) * exp(1.13 * tmin) /1.13 / (R * R * M * M) + N3g * v * exp(1.13 * tmin)/ 1.13 / pow(R * M, 4);//nb GeV^-2
  double result_upper_bound = N2g * v * pow(1.0 - X, 2) * exp(1.13 * tmax) /1.13 / (R * R * M * M) + N3g * v * exp(1.13 * tmax)/ 1.13 / pow(R * M, 4);//nb GeV^-2
  return result_upper_bound-result_lower_bound;
}

TF2 *f_dodt_23g = new TF2("psi2S_dodt_23g_function",psi2S_dodt_23g,10,20,-10,0);
TF1 *f_o_23g = new TF2("psi2S_o_23g_function",psi2S_o_23g,10,20);

int plot_dodt_vs_Egamma_project()
{

  int A = 1;
  double eIn_E = 17.0;
  
  // ************************************************************* //
  PSI2SMODEL::SetModel("23g");
  const double BR = 0.00793;
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
  myStyle->SetHistLineWidth(2);
  myStyle->SetPadGridY(true);
  myStyle->SetOptLogy();
  myStyle->SetTitleXSize(0.045); 
  myStyle->SetTitleYSize(0.045);
  myStyle->SetPadBottomMargin(0.15);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleYOffset(0.92);
  gROOT->SetStyle("Style");

  // ************************************************************* //
  
  const int nbins_x = 15;
  const double xmin = 0;
  const double xmax = 6;
  const double xstep = (xmax-xmin)/(1.0*nbins_x);
  TH1F *h_1 = new TH1F("h_1",";|t| [GeV^{2}];Counts",nbins_x,xmin,xmax);
  TH1F *h_2 = new TH1F("h_2",";|t| [GeV^{2}];Counts",nbins_x,xmin,xmax);
  TGraphErrors *tge[3];
  tge[0]= new TGraphErrors(nbins_x);
  tge[1]= new TGraphErrors(nbins_x);
  tge[2]= new TGraphErrors(nbins_x);

  tge[0]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
  tge[1]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
  tge[2]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");

  tge[0]->SetLineWidth(2);
  tge[0]->SetLineColor(kBlue);
  tge[0]->SetMarkerColor(kBlue); 
  tge[0]->SetMarkerStyle(20);

  tge[1]->SetLineWidth(2);
  tge[1]->SetLineColor(kRed);
  tge[1]->SetMarkerColor(kRed);
  tge[1]->SetMarkerStyle(20);

  tge[2]->SetLineWidth(2);
  tge[2]->SetLineColor(kViolet+2);
  tge[2]->SetMarkerColor(kViolet+2);
  tge[2]->SetMarkerStyle(20);

  // ************************************************************* //

  double Egamma_min = 12.2;
  double Egamma_max = 12.8;
  double Egamma_mid = (Egamma_max+Egamma_min)/2.0;

  tIn1->Project("h_1","-(event.t-tmin)",Form("weight_smear2*(Eg_true>%f&&Eg_true<%f)",Egamma_min,Egamma_max));
  tIn2->Project("h_2","-(event.t-tmin)",Form("weight_smear2*is_accept_eOut*(Eg_true>%f&&Eg_true<%f)",Egamma_min,Egamma_max));

  // ************************************************************* //
  // Get the nEvents for photo and electro per bin
  // Calculate the statistical and systematic error
  // Calculate the true dsigma/dt(Egamma, t)
  // Motivate the error in dsigma/dt by taking the errors in quadrature
  // Calculate the total error

  double tlimits_Egamma_mid[2];
  tlimits_Eg(Egamma_mid, tlimits_Egamma_mid);
  double tmin_Egamma_mid = -tlimits_Egamma_mid[1];
  double tmax_Egamma_mid = -tlimits_Egamma_mid[0];

  for(int i = 0 ; i < nbins_x ; i++)
    {
      double binCount_photo = h_1->GetBinContent(i+1);
      double binCount_electro = h_2->GetBinContent(i+1);
      double binCount_total = binCount_photo + binCount_electro;
      
      double statError_photo = 1.0/sqrt(binCount_photo);
      double statError_electro = 1.0/sqrt(binCount_electro);
      double statError_total = 1.0/sqrt(binCount_total);

      double sysError_photo = 0.1;
      double sysError_electro = 0.1;
      double sysError_total = 0.1;

      double dodt_true = f_dodt_23g->Eval(Egamma_mid,-(h_1->GetBinCenter(i+1)+tmin_Egamma_mid)); 
      
      double err_photo = dodt_true * sqrt( pow(statError_photo,2) + pow(sysError_photo,2) );
      double err_electro = dodt_true * sqrt( pow(statError_electro,2) + pow(sysError_electro,2) );
      double err_total = dodt_true * sqrt( pow(statError_total,2) + pow(sysError_total,2) );
      
      if(binCount_photo>0)
	{
	  tge[0]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	  tge[0]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_photo);
	}
      if(binCount_electro>0)
	{
	  tge[1]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	  tge[1]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_electro);
	}
      if(binCount_total>0)
	{ 
	  tge[2]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	  tge[2]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_total);
	}
    }

  // ************************************************************* //

  TLatex latex;
  latex.SetTextSize(0.032);
  latex.SetTextFont(42);
  double text_left = 0.33;
  double text_top  = 0.86;
  double text_step = 0.05;
  TCanvas *c0 = new TCanvas("c0","c0",800,600);
  tge[0]->Draw("AP");
  tge[0]->GetXaxis()->SetLimits(xmin,xmax);
  tge[0]->GetYaxis()->SetRangeUser(0.00001,1);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[4]{#bf{Photoprodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  tge[1]->Draw("AP");
  tge[1]->GetXaxis()->SetLimits(xmin,xmax);
  tge[1]->GetYaxis()->SetRangeUser(0.00001,1);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[2]{#bf{Electroprodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  tge[2]->Draw("AP");
  tge[2]->GetXaxis()->SetLimits(xmin,xmax);
  tge[2]->GetYaxis()->SetRangeUser(0.00001,1);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[882]{#bf{Photo+Electro prodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));
  c2->SetHighLightColor(0);
  return 0;
}

double x(double Egamma)
{
  const double Mp = 0.938272;
  const double Mv = 3.686;
  double X = (2.0 * Mp * Mv + Mv * Mv)/(2.0 * Mp * Egamma);
  return X;
}

void tlimits_W(double W, double * t)
{

  const double m1_2 = 0.00;
  const double m2_2 = pow(0.938272,2);
  const double m3_2 = pow(3.686,2);
  const double m4_2 = pow(0.938272,2);
  
  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  t[0] = m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) - 2 * p1 * p2;
  t[1] = m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) + 2 * p1 * p2;
  
  return;
}

void tlimits_Eg(double Eg, double * t)
{
  const double m1_2 = 0.00;
  const double m2_2 = pow(0.938272,2);
  const double m3_2 = pow(3.686,2);
  const double m4_2 = pow(0.938272,2);

  const double W = sqrt(2  * sqrt(m2_2) * Eg + m2_2);

  double p1 = sqrt(pow(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
  double p2 = sqrt(pow(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);

  t[0] = m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) - 2 * p1 * p2;
  t[1] = m2_2 + m4_2 - 2 * sqrt(m2_2 + p1*p1) * sqrt(m4_2 + p2*p2) + 2 * p1 * p2;
  
  return;
}
