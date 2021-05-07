#include "Lcore.h"

int SetStyle(TH1D * hB){
  hB->GetXaxis()->CenterTitle(true);
  hB->GetXaxis()->SetTitleSize(0.06);
  hB->GetXaxis()->SetTitleOffset(1.15);
  hB->GetXaxis()->SetLabelSize(0.06);
  hB->GetXaxis()->SetNdivisions(6, 5, 0);
  hB->GetYaxis()->CenterTitle(true);
  hB->GetYaxis()->SetTitleSize(0.06);
  hB->GetYaxis()->SetTitleOffset(1.15);
  hB->GetYaxis()->SetLabelSize(0.06);
  hB->GetYaxis()->SetNdivisions(6, 5, 0);
  hB->SetStats(0);
  return 0;
}

int SetStyle(TH2D * hB){
  hB->GetXaxis()->CenterTitle(true);
  hB->GetXaxis()->SetTitleSize(0.06);
  hB->GetXaxis()->SetTitleOffset(1.15);
  hB->GetXaxis()->SetLabelSize(0.06);
  hB->GetXaxis()->SetNdivisions(6, 5, 0);
  hB->GetYaxis()->CenterTitle(true);
  hB->GetYaxis()->SetTitleSize(0.06);
  hB->GetYaxis()->SetTitleOffset(1.15);
  hB->GetYaxis()->SetLabelSize(0.06);
  hB->GetYaxis()->SetNdivisions(6, 5, 0);
  hB->SetStats(0);
  return 0;
}
  

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./plotit <loadfile> <savefile> <opt>" << endl;
    return 1;
  }

  TString loadfile = argv[1];
  TString savefile = argv[2];
  const int opt = atoi(argv[3]);

  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);

  if (opt == 1){
    TFile * fs = new TFile(loadfile.Data(), "r");
    TH1D * hMJpsi = (TH1D *) fs->Get("Mass_e+e-_Jpsi");
    TH2D * hMomentum = (TH2D *) fs->Get("Pe+Pe-_Jpsi");
    TH2D * hThetaPelectron = (TH2D *) fs->Get("ThetaP_e-_Jpsi");
    TH2D * hThetaPpositron = (TH2D *) fs->Get("ThetaP_e+_Jpsi");
    TH2D * hThetaPproton = (TH2D *) fs->Get("ThetaP_proton_Jpsi");
    TH2D * hAngleP = (TH2D *) fs->Get("AngleP_Jpsi");
    TH1D * hFermiP = (TH1D *) fs->Get("FermiP_Jpsi");
    TH1D * hFermiPz = (TH1D *) fs->Get("FermiPz_Jpsi");
    TH1D * hPJpsi = (TH1D *) fs->Get("PJpsi");
    TH2D * hFermiPPJpsi = (TH2D *) fs->Get("FermiPPJpsi");
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(hMJpsi);
    hMJpsi->DrawClone("");
    c0->Print(savefile + "(", "pdf");

    SetStyle(hMomentum);
    hMomentum->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPelectron);
    hThetaPelectron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPpositron);
    hThetaPpositron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPproton);
    hThetaPproton->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hAngleP);
    hAngleP->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiP);
    hFermiP->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiPz);
    hFermiPz->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hPJpsi);
    hPJpsi->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiPPJpsi);
    hFermiPPJpsi->DrawClone("colz");     
    c0->Print(savefile + ")", "pdf");
   
  }

  if (opt == 2){
    TFile * fs = new TFile(loadfile.Data(), "r");
    TH1D * hMPhi = (TH1D *) fs->Get("Mass_e+e-_Phi");
    TH2D * hMomentum = (TH2D *) fs->Get("Pe+Pe-_Phi");
    TH2D * hThetaPelectron = (TH2D *) fs->Get("ThetaP_e-_Phi");
    TH2D * hThetaPpositron = (TH2D *) fs->Get("ThetaP_e+_Phi");
    TH2D * hThetaPproton = (TH2D *) fs->Get("ThetaP_proton_Phi");
    TH2D * hAngleP = (TH2D *) fs->Get("AngleP_Phi");
    TH1D * hFermiP = (TH1D *) fs->Get("FermiP_Phi");
    TH1D * hFermiPz = (TH1D *) fs->Get("FermiPz_Phi");
    TH1D * hPPhi = (TH1D *) fs->Get("PPhi");
    TH2D * hFermiPPPhi = (TH2D *) fs->Get("FermiPPPhi");
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(hMPhi);
    hMPhi->DrawClone("");
    c0->Print(savefile + "(", "pdf");

    SetStyle(hMomentum);
    hMomentum->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPelectron);
    hThetaPelectron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPpositron);
    hThetaPpositron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPproton);
    hThetaPproton->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hAngleP);
    hAngleP->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiP);
    hFermiP->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiPz);
    hFermiPz->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hPPhi);
    hPPhi->DrawClone("");
    c0->Print(savefile, "pdf");

    SetStyle(hFermiPPPhi);
    hFermiPPPhi->DrawClone("colz");     
    c0->Print(savefile + ")", "pdf");
   
  }


  if (opt == 11){
    TFile * fs = new TFile(loadfile.Data(), "r");
    TH1D * hMJpsi = (TH1D *) fs->Get("Mass_e+e-_Jpsi");
    TH2D * hMomentum = (TH2D *) fs->Get("Pe+Pe-_Jpsi");
    TH2D * hThetaPelectron = (TH2D *) fs->Get("ThetaP_e-_Jpsi");
    TH2D * hThetaPpositron = (TH2D *) fs->Get("ThetaP_e+_Jpsi");
    TH2D * hThetaPJpsi = (TH2D *) fs->Get("ThetaPJpsi");
    TH1D * hPmin = (TH1D *) fs->Get("Pmin_Jpsi");
        
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(hMJpsi);
    hMJpsi->DrawClone("");
    c0->Print(savefile + "(", "pdf");

    SetStyle(hMomentum);
    hMomentum->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPelectron);
    hThetaPelectron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPpositron);
    hThetaPpositron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPJpsi);
    hThetaPJpsi->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hPmin);
    hPmin->DrawClone("colz");
    c0->Print(savefile + ")", "pdf");

  }

  if (opt == 111){
    TFile * fs = new TFile(loadfile.Data(), "r");
    TH1D * hMJpsi = (TH1D *) fs->Get("Mass_e+e-_Jpsi");
    TH2D * hMomentum = (TH2D *) fs->Get("Pe+Pe-_Jpsi");
    TH2D * hThetaPelectron = (TH2D *) fs->Get("ThetaP_e-_Jpsi");
    TH2D * hThetaPpositron = (TH2D *) fs->Get("ThetaP_e+_Jpsi");
    TH2D * hThetaPJpsi = (TH2D *) fs->Get("ThetaPJpsi");
    //TH1D * hPmin = (TH1D *) fs->Get("Pmin_Jpsi");
        
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(hMJpsi);
    hMJpsi->DrawClone("");
    c0->Print(savefile + "(", "pdf");

    SetStyle(hMomentum);
    hMomentum->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPelectron);
    hThetaPelectron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPpositron);
    hThetaPpositron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPJpsi);
    hThetaPJpsi->DrawClone("colz");
    c0->Print(savefile + ")", "pdf");

    //SetStyle(hPmin);
    //hPmin->DrawClone("colz");
    //c0->Print(savefile + ")", "pdf");

  }

  if (opt == 22){
    TFile * fs = new TFile(loadfile.Data(), "r");
    TH1D * hMPhi = (TH1D *) fs->Get("Mass_e+e-_Phi");
    TH2D * hMomentum = (TH2D *) fs->Get("Pe+Pe-_Phi");
    TH2D * hThetaPelectron = (TH2D *) fs->Get("ThetaP_e-_Phi");
    TH2D * hThetaPpositron = (TH2D *) fs->Get("ThetaP_e+_Phi");
    TH2D * hThetaPphi = (TH2D *) fs->Get("ThetaPPhi");
      
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(hMPhi);
    hMPhi->DrawClone("");
    c0->Print(savefile + "(", "pdf");

    SetStyle(hMomentum);
    hMomentum->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPelectron);
    hThetaPelectron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPpositron);
    hThetaPpositron->DrawClone("colz");
    c0->Print(savefile, "pdf");

    SetStyle(hThetaPphi);
    hThetaPphi->DrawClone("colz");
    c0->Print(savefile + ")", "pdf");
   
  }

  return 0;
}

  
