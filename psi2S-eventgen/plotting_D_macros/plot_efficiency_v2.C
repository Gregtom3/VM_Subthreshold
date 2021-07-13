int plot_efficiency_v2()
{
  //  gROOT->SetBatch(true);
  
  // ****************************************************** //
  int A = 1;
  double beamE = 17;
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

  fIn[0] = new TFile(Form("../data/%s-photo-500000-%.0f/%s_solid_photo_%.1fGeV.root",nuc.Data(),beamE,nuc.Data(),beamE),"READ");
  fIn[1] = new TFile(Form("../data/%s-electro-10000000-%.0f/%s_solid_electro_%.1fGeV.root",nuc.Data(),beamE,nuc.Data(),beamE),"READ");
  TTree *tIn[2];
  tIn[0] = (TTree*)fIn[0]->Get("tree");
  tIn[1] = (TTree*)fIn[1]->Get("tree");
  
  TStyle *myStyle = new TStyle("Style","My Style");
  myStyle->SetOptStat(0);
  myStyle->SetNdivisions(505);
  myStyle->SetHistLineWidth(2);
  myStyle->SetPadGridY(true);
  myStyle->SetTitleXSize(0.045); 
  myStyle->SetTitleYSize(0.045);
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleYOffset(0.92);
  myStyle->SetPadRightMargin(0.05);
  myStyle->SetHistLineWidth(2);
  gROOT->SetStyle("Style");

  TH1F *hEvents[2];
  TH1F *hEvents_Acc_2fold[2];
  TH1F *hEvents_Acc_3fold_e[2];
  TH1F *hEvents_Acc_3fold_p[2];
  TH1F *hEvents_Acc_4fold[2];
  
  TCanvas *c[2];
  double Emin = 9.5;
  double Emax = 17;
  int N = 25;

 TLatex latex;
 latex.SetTextSize(0.03);
 latex.SetTextFont(42);
 TLegend *legend[2];

  for(int m = 0 ; m < 2 ; m++)
    {
      hEvents[m] = new TH1F(Form("hEvents%d",m),"",N,Emin,Emax);
      hEvents_Acc_2fold[m] = new TH1F(Form("hEvents_Acc_2fold%d",m),"",N,Emin,Emax);
      hEvents_Acc_3fold_e[m] = new TH1F(Form("hEvents_Acc_3fold_e%d",m),"",N,Emin,Emax);
      hEvents_Acc_3fold_p[m] = new TH1F(Form("hEvents_Acc_3fold_p%d",m),"",N,Emin,Emax);
      hEvents_Acc_4fold[m] = new TH1F(Form("hEvents_Acc_4fold%d",m),"",N,Emin,Emax);

      if(m==0) // photo
	{
	  hEvents_Acc_2fold[m]->SetLineColor(860);
	  hEvents_Acc_3fold_p[m]->SetLineColor(863);
	}
      else if(m==1) // electro
	{
	  hEvents_Acc_2fold[m]->SetLineColor(627);
	  hEvents_Acc_3fold_p[m]->SetLineColor(625);
	  hEvents_Acc_3fold_e[m]->SetLineColor(632);
	  hEvents_Acc_4fold[m]->SetLineColor(635);
	}
      tIn[m]->Draw(Form("Eg_true>>hEvents%d",m),"weight_total * time * lumi / Nsim","goff");
      tIn[m]->Draw(Form("Eg_true>>hEvents_Acc_2fold%d",m),"weight_total * time * lumi / Nsim * is_accept_ePlusOut * is_accept_eMinusOut","goff");
      tIn[m]->Draw(Form("Eg_true>>hEvents_Acc_3fold_p%d",m),"weight_total * time * lumi / Nsim * is_accept_ePlusOut * is_accept_eMinusOut * is_accept_pOut","goff");
      if(m==1)
	{
	  tIn[m]->Draw(Form("Eg_true>>hEvents_Acc_3fold_e%d",m),"weight_total * time * lumi / Nsim * is_accept_ePlusOut * is_accept_eMinusOut * is_accept_eOut","goff");
	  tIn[m]->Draw(Form("Eg_true>>hEvents_Acc_4fold%d",m),"weight_total * time * lumi / Nsim * is_accept_ePlusOut * is_accept_eMinusOut * is_accept_pOut * is_accept_eOut","goff");
	}

      hEvents_Acc_2fold[m]->Divide(hEvents[m]);
      hEvents_Acc_3fold_p[m]->Divide(hEvents[m]);
      if(m==1)
	{
	  hEvents_Acc_3fold_e[m]->Divide(hEvents[m]);
	  hEvents_Acc_4fold[m]->Divide(hEvents[m]);
	}

      
      c[m] = new TCanvas(Form("c%d",m),"c",1200,600);
      hEvents_Acc_2fold[m]->Draw("][");
      hEvents_Acc_2fold[m]->GetYaxis()->SetRangeUser(0,2);
      hEvents_Acc_2fold[m]->GetXaxis()->SetRangeUser(Emin,Emax);
      hEvents_Acc_3fold_p[m]->Draw("][ same");
      legend[m]=new TLegend(0.567,0.669,0.927,0.794);
      legend[m]->AddEntry(hEvents_Acc_2fold[m],"2-fold (e^{-}e^{+})","l");
      legend[m]->AddEntry(hEvents_Acc_3fold_p[m],"3-fold p'(e^{-}e^{+})","l");
      if(m==1)
	{
	  hEvents_Acc_3fold_e[m]->Draw("][ same");
	  hEvents_Acc_4fold[m]->Draw("][ same");
	  legend[m]->AddEntry(hEvents_Acc_3fold_e[m],"3-fold e'(e^{-}e^{+})","l");
	  legend[m]->AddEntry(hEvents_Acc_4fold[m],"4-fold e'p'(e^{-}e^{+})","l");
	}
      legend[m]->Draw("same");
    }


  
  
  return 0;
}
