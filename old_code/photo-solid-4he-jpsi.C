#include "Lcore.h"

int CheckAcceptance(const TLorentzVector P, const double thmin, const double thmax){
  if (P.Theta() < thmin * M_PI / 180.0) return 0;
  if (P.Theta() > thmax * M_PI / 180.0) return 0;
  if (P.P() < 0.3) return 0;
  return 1;
}

int main(const int argc, const char * argv[]){

  gErrorIgnoreLevel = 5000;
  // Set simulation
  gRandom->SetSeed(0);

  if (argc < 3){
    cout << "./photo-solid-4he-jpsi <model> <Nsim>" << endl;
    return 0;
  }
  Long64_t Nsim = atoi(argv[2]);

  // Electron beam energy and luminosity
  double Ebeam = 7.2;//GeV
  double lumi = 1.2e37 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
    
  // Set Jpsi production model
  JPSIHE4::SetModel(argv[1]);

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 6.2;
  double kmax = 7.2;
  
  TString app = argv[1];
  TFile * fraw = new TFile("result-photon/He4raw-"+app+".root", "RECREATE");
  TH1D * hMJpsi1 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum1 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPJpsi1 = new TH2D("ThetaPJpsi", ";#theta[J/#psi] (deg);P[J/#psi] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hPmin1 = new TH1D("Pmin_Jpsi", ";P[p]_{min} (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hPmax1 = new TH1D("Pmax_Jpsi", ";P[p]_{max} (GeV);Events / hour", 100, 0.0, 1.0);
  hMJpsi1->SetDirectory(fraw);
  hMomentum1->SetDirectory(fraw);
  hThetaPelectron1->SetDirectory(fraw);
  hThetaPpositron1->SetDirectory(fraw);
  hThetaPJpsi1->SetDirectory(fraw);
  hPmin1->SetDirectory(fraw);
  hPmax1->SetDirectory(fraw);

  // detected subthreshold 1
  TFile * fsub = new TFile("result-photon/He4sub-"+app+".root", "RECREATE");
  TH1D * hMJpsi2 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum2 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron2 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron2 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPJpsi2 = new TH2D("ThetaPJpsi", ";#theta[J/#psi] (deg);P[J/#psi] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hPmin2 = new TH1D("Pmin_Jpsi", ";P[p]_{min} (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hPmax2 = new TH1D("Pmax_Jpsi", ";P[p]_{max} (GeV);Events / hour", 100, 0.0, 1.0);
  hMJpsi2->SetDirectory(fsub);
  hMomentum2->SetDirectory(fsub);
  hThetaPelectron2->SetDirectory(fsub);
  hThetaPpositron2->SetDirectory(fsub);
  hThetaPJpsi2->SetDirectory(fsub);
  hPmin2->SetDirectory(fsub);
  hPmax2->SetDirectory(fsub);

  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam) * 1.99 / 2;//15cm LHe4 target
    weight *= GENERATE::Event_g4He2ee_Jpsi(ki, kf);

    weight *= CheckAcceptance(kf[0], 6.0, 29.0) * CheckAcceptance(kf[1], 6.0, 29.0);

    if (weight > 0.0){

      hMJpsi1->Fill( (kf[0]+kf[1]).M(), weight);
      hMomentum1->Fill( kf[0].P(), kf[1].P(), weight);
      hThetaPelectron1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
      hThetaPpositron1->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
      hThetaPJpsi1->Fill( (kf[0]+kf[1]).Theta()/M_PI*180, (kf[0]+kf[1]).P(), weight);
      //hPmin1->Fill( JPSIHE4::hpmin->GetBinContent(JPSIHE4::hpmin->FindBin(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0)), weight);
      //hPmax1->Fill( JPSIHE4::hpmax->GetBinContent(JPSIHE4::hpmax->FindBin(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0)), weight);
      hPmin1->Fill( JPSIHE4::hpmin->Interpolate(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0), weight);
      hPmax1->Fill( JPSIHE4::hpmax->Interpolate(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0), weight);
      if (ki[0].E() < 8.2){
	hMJpsi2->Fill( (kf[0]+kf[1]).M(), weight);
	hMomentum2->Fill( kf[0].P(), kf[1].P(), weight);
	hThetaPelectron2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPpositron2->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hThetaPJpsi2->Fill( (kf[0]+kf[1]).Theta()/M_PI*180, (kf[0]+kf[1]).P(), weight);
	//hPmin2->Fill( JPSIHE4::hpmin->GetBinContent(JPSIHE4::hpmin->FindBin(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0)), weight);
	//hPmax2->Fill( JPSIHE4::hpmax->GetBinContent(JPSIHE4::hpmax->FindBin(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0)), weight);
        hPmin2->Fill( JPSIHE4::hpmin->Interpolate(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0), weight);
        hPmax2->Fill( JPSIHE4::hpmax->Interpolate(ki[0].E(), (kf[0]+kf[1]).P(), (kf[0]+kf[1]).Angle(ki[0].Vect()) / M_PI * 180.0), weight);
      }
    }
    
  }

  hMJpsi1->Scale(lumi*time/Nsim);
  hMomentum1->Scale(lumi*time/Nsim);
  hThetaPelectron1->Scale(lumi*time/Nsim);
  hThetaPpositron1->Scale(lumi*time/Nsim);
  hThetaPJpsi1->Scale(lumi*time/Nsim);
  hPmin1->Scale(lumi*time/Nsim);
  hPmax1->Scale(lumi*time/Nsim);
  fraw->Write();
  fraw->Close();

  hMJpsi2->Scale(lumi*time/Nsim);
  hMomentum2->Scale(lumi*time/Nsim);
  hThetaPelectron2->Scale(lumi*time/Nsim);
  hThetaPpositron2->Scale(lumi*time/Nsim);
  hThetaPJpsi2->Scale(lumi*time/Nsim);
  hPmin2->Scale(lumi*time/Nsim);
  hPmax2->Scale(lumi*time/Nsim);
  fsub->Write();
  fsub->Close();

  return 0;
}
