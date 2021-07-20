#include "../Lcore.h"
double x(double Egamma);
void tlimits_W(double W, double * t);
void tlimits_Eg(double Eg, double * t);
double psi2S_dodt_23g(double * xx, double * p)
{
  double X = x(xx[0]);
  double t = xx[1];
  double N2g = 1569.43;
  double N3g = 698.868;
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

  double N2g = 1569.43;
  double N3g = 698.868;
  double v = 1.0 / (16.0 * M_PI);
  double R = 1.0;
  double M = 3.686097;//GeV
  double result_lower_bound = N2g * v * pow(1.0 - X, 2) * exp(1.13 * tmin) /1.13 / (R * R * M * M) + N3g * v * exp(1.13 * tmin)/ 1.13 / pow(R * M, 4);//nb GeV^-2
  double result_upper_bound = N2g * v * pow(1.0 - X, 2) * exp(1.13 * tmax) /1.13 / (R * R * M * M) + N3g * v * exp(1.13 * tmax)/ 1.13 / pow(R * M, 4);//nb GeV^-2
  return result_upper_bound-result_lower_bound;
}

TF2 *f_dodt_23g = new TF2("psi2S_dodt_23g_function",psi2S_dodt_23g,10,20,-10,0);
TF1 *f_o_23g = new TF1("psi2S_o_23g_function",psi2S_o_23g,10.93,20);

int plot_dodt_vs_Egamma_many()
{
  gROOT->SetBatch(true);
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
  
  double Egamma_bottom = 11;
  double Egamma_top = 17;
  int n_Egamma = 8;
  double xVals[n_Egamma];
  double sigma[n_Egamma][3];
  double sigma_err[n_Egamma][3];
  double Egamma_step = (Egamma_top - Egamma_bottom)/n_Egamma;
  TCanvas *c0[n_Egamma];
  TCanvas *c1[n_Egamma];
  TCanvas *c2[n_Egamma];

  // ************************************************************* //

  const int nbins_x = 20;
  const double xmin = 0;
  const double xmax = 6;
  const double xstep = (xmax-xmin)/(1.0*nbins_x);

  // ************************************************************* //

  TH1F *h_1 = new TH1F("h_1",";|t| [GeV^{2}];Counts",nbins_x,xmin,xmax);
  TH1F *h_2 = new TH1F("h_2",";|t| [GeV^{2}];Counts",nbins_x,xmin,xmax);
  TGraphErrors *tge[3][n_Egamma];
  for(int k = 0 ; k < n_Egamma; k++)
    {
      tge[0][k]= new TGraphErrors(nbins_x);
      tge[1][k]= new TGraphErrors(nbins_x);
      tge[2][k]= new TGraphErrors(nbins_x);
      
      tge[0][k]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
      tge[1][k]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
      tge[2][k]->SetTitle(";|t-t_{min}| [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
      
      tge[0][k]->SetLineWidth(2);
      tge[0][k]->SetLineColor(kBlue);
      tge[0][k]->SetMarkerColor(kBlue); 
      tge[0][k]->SetMarkerStyle(20);
      
      tge[1][k]->SetLineWidth(2);
      tge[1][k]->SetLineColor(kRed);
      tge[1][k]->SetMarkerColor(kRed);
      tge[1][k]->SetMarkerStyle(20);
      
      tge[2][k]->SetLineWidth(2);
      tge[2][k]->SetLineColor(kViolet+2);
      tge[2][k]->SetMarkerColor(kViolet+2);
      tge[2][k]->SetMarkerStyle(20);
    }
  
  // ************************************************************* //

  TF1 *fit_photo = new TF1("fit_photo","[0]*exp([1]*x)",0,100);
  TF1 *fit_electro = new TF1("fit_electro","[0]*exp([1]*x)",0,100);
  TF1 *fit_total = new TF1("fit_total","[0]*exp([1]*x)",0,100);


  fit_photo->SetLineStyle(2);
  fit_electro->SetLineStyle(2);
  fit_total->SetLineStyle(2);

  fit_photo->SetLineWidth(1);
  fit_electro->SetLineWidth(1);
  fit_total->SetLineWidth(1);

  // ************************************************************* //

  for(int j = 0 ; j < n_Egamma; j++)
    {
      double Egamma_min = Egamma_bottom + (j)*Egamma_step;
      double Egamma_max = Egamma_min + Egamma_step;
      double Egamma_mid = (Egamma_max+Egamma_min)/2.0;
      xVals[j] = Egamma_mid;
      h_1->Reset();
      h_2->Reset();

      tIn1->Draw("-(event.t-tmin)>>h_1",Form("weight_smear2*(Eg_true>%f&&Eg_true<%f)",Egamma_min,Egamma_max),"goff");
      tIn2->Draw("-(event.t-tmin)>>h_2",Form("weight_smear2*is_accept_eOut*(Eg_true>%f&&Eg_true<%f)",Egamma_min,Egamma_max),"goff");
      
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
	      tge[0][j]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	      tge[0][j]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_photo);
	    }
	  if(binCount_electro>0)
	    {
	      tge[1][j]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	      tge[1][j]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_electro);
	    }
	  if(binCount_total>0)
	    { 
	      tge[2][j]->SetPoint(i, h_1->GetBinCenter(i+1) , dodt_true);
	      tge[2][j]->SetPointError(i, 0.5 * h_1->GetBinWidth(i+1) , err_total);
	    }
	}

      // ************************************************************* //

      fit_photo->SetParameters(0.05,-1.13);
      fit_electro->SetParameters(0.05,-1.13);
      fit_total->SetParameters(0.05,-1.13);
      
      fit_photo->SetParLimits(0,0.0001,5);
      fit_electro->SetParLimits(0,0.0001,5);
      fit_total->SetParLimits(0,0.0001,5);

      TFitResultPtr photo_r = tge[0][j]->Fit(fit_photo,"RNQS");
      TFitResultPtr electro_r = tge[1][j]->Fit(fit_electro,"RNQS");
      TFitResultPtr total_r = tge[2][j]->Fit(fit_total,"RNQS");

      double photo_A = fit_photo->GetParameter(0);
      double photo_B = fit_photo->GetParameter(1);
      double photo_A_err = fit_photo->GetParError(0);
      double photo_B_err = fit_photo->GetParError(1);
      double photo_sigma = fit_photo->Integral(0, tmax_Egamma_mid - tmin_Egamma_mid);
      double photo_sigma_err = fit_photo->IntegralError(0, tmax_Egamma_mid - tmin_Egamma_mid,photo_r->GetParams(), photo_r->GetCovarianceMatrix().GetMatrixArray());

      double electro_A = fit_electro->GetParameter(0);
      double electro_B = fit_electro->GetParameter(1);
      double electro_A_err = fit_electro->GetParError(0);
      double electro_B_err = fit_electro->GetParError(1);
      double electro_sigma = fit_electro->Integral(0, tmax_Egamma_mid - tmin_Egamma_mid);
      double electro_sigma_err = fit_electro->IntegralError(0, tmax_Egamma_mid - tmin_Egamma_mid,electro_r->GetParams(), electro_r->GetCovarianceMatrix().GetMatrixArray());

      double total_A = fit_total->GetParameter(0);
      double total_B = fit_total->GetParameter(1);
      double total_A_err = fit_total->GetParError(0);
      double total_B_err = fit_total->GetParError(1);
      double total_sigma = fit_total->Integral(0, tmax_Egamma_mid - tmin_Egamma_mid);
      double total_sigma_err = fit_total->IntegralError(0, tmax_Egamma_mid - tmin_Egamma_mid,total_r->GetParams(), total_r->GetCovarianceMatrix().GetMatrixArray());

      // ************************************************************* //

      sigma[j][0] = photo_sigma;
      sigma[j][1] = electro_sigma;
      sigma[j][2] = total_sigma;
      sigma_err[j][0] = photo_sigma_err;
      sigma_err[j][1] = electro_sigma_err;
      sigma_err[j][2] = total_sigma_err;

      // ************************************************************* //
      
      TLatex latex;
      latex.SetTextSize(0.032);
      latex.SetTextFont(42);
      double text_left = 0.33;
      double text_top  = 0.86;
      double text_step = 0.05;
      c0[j] = new TCanvas(Form("c0%d",j),Form("c0%d",j),800,600);
      tge[0][j]->Draw("AP");
      fit_photo->Draw("C same");
      tge[0][j]->Draw("P same");
      tge[0][j]->GetXaxis()->SetLimits(xmin,xmax);
      tge[0][j]->GetYaxis()->SetRangeUser(0.00001,1);
      latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
      latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[4]{#bf{Photoprodution}}",stringA.Data()));
      latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));
      latex.DrawLatexNDC(text_left+text_step,text_top-3*text_step,Form("d#sigma/dt = A exp(B t) , #sigma = %.5f +/- %.5f [nb]", photo_sigma, photo_sigma_err));
      latex.DrawLatexNDC(text_left+text_step,text_top-4*text_step,Form("A = %.3f +/- %.3f  , B = %.3f +/- %.3f",photo_A,photo_A_err,photo_B,photo_B_err));
      c0[j]->SaveAs(Form("dodt_vs_Egamma_plots/%s-photo-%.2f-Egamma-%.2f.pdf",stringA.Data(),Egamma_min,Egamma_max));
      
      c1[j] = new TCanvas(Form("c1%d",j),Form("c1%d",j),800,600);
      tge[1][j]->Draw("AP");
      fit_electro->Draw("C same");
      tge[1][j]->Draw("P same");
      tge[1][j]->GetXaxis()->SetLimits(xmin,xmax);
      tge[1][j]->GetYaxis()->SetRangeUser(0.00001,1);
      latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
      latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[2]{#bf{Electroprodution}}",stringA.Data()));
      latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));
      latex.DrawLatexNDC(text_left+text_step,text_top-3*text_step,Form("d#sigma/dt = A exp(B t) , #sigma = %.5f +/- %.5f [nb]", electro_sigma, electro_sigma_err));
      latex.DrawLatexNDC(text_left+text_step,text_top-4*text_step,Form("A = %.3f +/- %.3f  , B = %.3f +/- %.3f",electro_A,electro_A_err,electro_B,electro_B_err));
      c1[j]->SaveAs(Form("dodt_vs_Egamma_plots/%s-electro-%.2f-Egamma-%.2f.pdf",stringA.Data(),Egamma_min,Egamma_max));

      c2[j] = new TCanvas(Form("c2%d",j),Form("c2%d",j),800,600);
      tge[2][j]->Draw("AP");
      fit_total->Draw("C same");
      tge[2][j]->Draw("P same");
      tge[2][j]->GetXaxis()->SetLimits(xmin,xmax);
      tge[2][j]->GetYaxis()->SetRangeUser(0.00001,1);
      latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
      latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[882]{#bf{Photo+Electro prodution}}",stringA.Data()));
      latex.DrawLatexNDC(text_left,text_top-2*text_step,Form("%.2f GeV < E_{#gamma} < %.2f GeV",Egamma_min,Egamma_max));
      latex.DrawLatexNDC(text_left+text_step,text_top-3*text_step,Form("d#sigma/dt = A exp(B t) , #sigma = %.5f +/- %.5f [nb]", total_sigma, total_sigma_err));
      latex.DrawLatexNDC(text_left+text_step,text_top-4*text_step,Form("A = %.3f +/- %.3f  , B = %.3f +/- %.3f",total_A,total_A_err,total_B,total_B_err));
      c2[j]->SetHighLightColor(0);
      c2[j]->SaveAs(Form("dodt_vs_Egamma_plots/%s-total-%.2f-Egamma-%.2f.pdf",stringA.Data(),Egamma_min,Egamma_max));


      // ************************************************************* //
      

    }

  f_o_23g->SetLineStyle(9);
  f_o_23g->SetLineColor(kBlack);

  TGraphErrors *tge_sigma[3];
  for(int i = 0 ; i < 3 ; i++)
    {
      tge_sigma[i]=new TGraphErrors(n_Egamma);
      tge_sigma[i]->SetMarkerStyle(20);
      tge_sigma[i]->SetLineWidth(2);
      tge_sigma[i]->SetTitle(";E_{#gamma}[GeV];#sigma (#gamma+p#rightarrow #psi(2S)+p') [nb]");
  }

  tge_sigma[0]->SetLineColor(kBlue);
  tge_sigma[1]->SetLineColor(kRed);
  tge_sigma[2]->SetLineColor(882);

  tge_sigma[0]->SetMarkerColor(kBlue);
  tge_sigma[1]->SetMarkerColor(kRed);
  tge_sigma[2]->SetMarkerColor(882);

  for(int i = 0 ; i < n_Egamma; i++)
    {
      tge_sigma[0]->SetPoint(i,xVals[i],sigma[i][0]);
      tge_sigma[1]->SetPoint(i,xVals[i],sigma[i][1]);
      tge_sigma[2]->SetPoint(i,xVals[i],sigma[i][2]);

      tge_sigma[0]->SetPointError(i,Egamma_step/2.0,sigma_err[i][0]);
      tge_sigma[1]->SetPointError(i,Egamma_step/2.0,sigma_err[i][1]);
      tge_sigma[2]->SetPointError(i,Egamma_step/2.0,sigma_err[i][2]);
    }
  
  TLatex latex;
  latex.SetTextSize(0.032);
  latex.SetTextFont(42);
  double text_left = 0.33;
  double text_top  = 0.86;
  double text_step = 0.05;
  
  TCanvas *c_photo = new TCanvas("c_photo","c_photo",800,600);
  gPad->SetLogy();
  tge_sigma[0]->Draw("AP");
  tge_sigma[0]->GetXaxis()->SetLimits(10,20);
  tge_sigma[0]->GetYaxis()->SetRangeUser(0.001,10);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[4]{#bf{Photoprodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,"50 days at 1.2e37/cm^{2}s");
  
  f_o_23g->Draw("C same");
  c_photo->SaveAs("dodt_vs_Egamma_plots/c_photo.pdf");

  TCanvas *c_electro = new TCanvas("c_electro","c_electro",800,600);
  gPad->SetLogy();
  tge_sigma[1]->Draw("AP");
  tge_sigma[1]->GetXaxis()->SetLimits(10,20);
  tge_sigma[1]->GetYaxis()->SetRangeUser(0.001,10);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[2]{#bf{Electroprodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,"50 days at 1.2e37/cm^{2}s");
  f_o_23g->Draw("C same");
  c_electro->SaveAs("dodt_vs_Egamma_plots/c_electro.pdf");

  TCanvas *c_total = new TCanvas("c_total","c_total",800,600);
  gPad->SetLogy();
  tge_sigma[2]->Draw("AP");
  tge_sigma[2]->GetXaxis()->SetLimits(10,20);
  tge_sigma[2]->GetYaxis()->SetRangeUser(0.001,10);
  latex.DrawLatexNDC(text_left,text_top,Form("#bf{SoLID} Simulation (%.1f GeV beam)",eIn_E));
  latex.DrawLatexNDC(text_left,text_top-text_step,Form("e+%s #rightarrow e'+p'+#psi(2S) (e^{-}e^{+}) , #color[882]{#bf{Photo+Electro prodution}}",stringA.Data()));
  latex.DrawLatexNDC(text_left,text_top-2*text_step,"50 days at 1.2e37/cm^{2}s");
  f_o_23g->Draw("C same");
  c_total->SaveAs("dodt_vs_Egamma_plots/c_total.pdf"); 
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
