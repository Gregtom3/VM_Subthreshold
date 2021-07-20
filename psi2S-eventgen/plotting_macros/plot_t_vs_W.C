int plot_t_vs_W()
{
  int A = 1;
  int type = 0; // 0 --> photo , 1 --> electro
  double eIn_E = 17.0;
  
  // ************************************************************* //

  TString stringA = (A==1) ? "p" : "D";
  TString stringtype = (type==0) ? "photo" : "electro";
  TFile *fIn = new TFile(Form("../result-%s-psi2S/%s_solid_%s_%.1fGeV.root",stringtype.Data(),stringA.Data(),stringtype.Data(),eIn_E),"READ");
  TTree *tIn = (TTree*)fIn->Get("tree"); 
  
  // ************************************************************* //
  
  const int nbins_x = 50;
  const int nbins_y = 50;
  const double xmin = 4.4;
  const double xmax = 5.8;
  const double ymin = 0;
  const double ymax = 8;
  TH2F *h_1 = new TH2F("h_1",Form("e(%.1fGeV)+%s , Psi2S %s-production , No Acceptance;W [GeV];|t-t_{min}| [GeV^{2}]",eIn_E,stringA.Data(),stringtype.Data()),nbins_x,xmin,xmax,nbins_y,ymin,ymax);
  TH2F *h_2 = new TH2F("h_2",Form("e(%.1fGeV)+%s , Psi2S %s-production , Decay Accepted;W [GeV];|t-t_{min}| [GeV^{2}]",eIn_E,stringA.Data(),stringtype.Data()),nbins_x,xmin,xmax,nbins_y,ymin,ymax);
  TH2F *h_3 = new TH2F("h_3",Form("e(%.1fGeV)+%s , Psi2S %s-production , Decay+pOut Accepted ;W [GeV];|t-t_{min}| [GeV^{2}]",eIn_E,stringA.Data(),stringtype.Data()),nbins_x,xmin,xmax,nbins_y,ymin,ymax);
  
  // ************************************************************* //
  
  const double text_xmin = 0.35;
  const double text_xmax = 0.835;
  const double text_ymin = 0.75;
  const double text_ymax = 0.825;
  TCanvas *c = new TCanvas("c","c",2000,600);
  gStyle->SetOptStat(0);
  c->Divide(3,1,0.001,0.001);
  
  c->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);
  tIn->Draw("-(event.t-tmin):event.W>>h_1","weight_lumi","colz");
  TPaveText *pt_1 = new TPaveText(text_xmin,text_ymin,text_xmax,text_ymax,"NDC NB");
  pt_1->SetFillColor(kWhite);
  pt_1->AddText(Form("Nevents per 1.2e37 cm^{-2}s^{-1}A^{-1} = %.0f" , h_1->Integral()));
  pt_1->Draw("same");

  c->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);
  tIn->Draw("-(event.t-tmin):event.W>>h_2","weight_lumi*is_accept_ePlusOut*is_accept_eMinusOut","colz");
  TPaveText *pt_2=new TPaveText(text_xmin,text_ymin,text_xmax,text_ymax,"NDC NB");
  pt_2->SetFillColor(kWhite);
  pt_2->AddText(Form("Nevents per 1.2e37 cm^{-2}s^{-1}A^{-1} = %.0f" , h_2->Integral()));
  pt_2->Draw("same");

  c->cd(3);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);
  tIn->Draw("-(event.t-tmin):event.W>>h_3","weight_smear2","colz");
  TPaveText *pt_3 = new TPaveText(text_xmin,text_ymin,text_xmax,text_ymax,"NDC NB");
  pt_3->SetFillColor(kWhite);
  pt_3->AddText(Form("Nevents per 1.2e37 cm^{-2}s^{-1}A^{-1} = %.0f" , h_3->Integral()));
  pt_3->Draw("same");

  return 0;
}
