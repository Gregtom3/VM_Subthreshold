#include "Lcore.h"

int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 100000000;

  // Electron beam energy and luminosity
  double Ebeam = 10.6;//GeV
  double lumi = 2.0e35 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
  
  // Set nuclear
  NUCLEAR::SetNuclear("p");  
  //GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  //GENERATE::TF_fMomentum->SetNpx(1000);
    
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set scattered electron range
  GENERATE::cthrange[0] = cos(4.5/180.0 * M_PI);
  GENERATE::cthrange[1] = cos(2.5/180.0 * M_PI);
  GENERATE::perange[0] = 0.5;//GeV
  GENERATE::perange[1] = 5.0;//GeV

  // raw
  TFile * fraw = new TFile("result/praw.root", "RECREATE");  
  TH1D * hMJpsi = new TH1D("Mass_e+e-_Jpsi", "", 100, 2.6, 3.6);
  TH2D * hMomentum = new TH2D("Pe+Pe-_Jpsi", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron = new TH2D("ThetaP_e-_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron = new TH2D("ThetaP_e+_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP = new TH2D("AngleP_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  hMJpsi->SetDirectory(fraw);
  hMomentum->SetDirectory(fraw);
  hThetaPelectron->SetDirectory(fraw);
  hThetaPpositron->SetDirectory(fraw);
  hAngleP->SetDirectory(fraw);

  // subthreshold
  TFile * fsub = new TFile("result/psub.root", "RECREATE");
  TH1D * hMJpsi1 = new TH1D("Mass_e+e-_Jpsi", "", 100, 2.6, 3.6);
  TH2D * hMomentum1 = new TH2D("Pe+Pe-_Jpsi", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_e-_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_e+_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP1 = new TH2D("AngleP_Jpsi", "", 90, 0.0, 180.0, 100, 0.0, 10.0);
  hMJpsi1->SetDirectory(fsub);
  hMomentum1->SetDirectory(fsub);
  hThetaPelectron1->SetDirectory(fsub);
  hThetaPpositron1->SetDirectory(fsub);
  hAngleP1->SetDirectory(fsub);
  
  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight > 0.0){
      hMJpsi->Fill( (kf[2]+kf[3]).M(), weight);
      hMomentum->Fill( kf[2].P(), kf[3].P(), weight);
      hThetaPelectron->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
      hThetaPpositron->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
      hAngleP->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);

      if (kf[0].E() > Ebeam - 8.2){
	hMJpsi1->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum1->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron1->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hAngleP1->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
      }
    }
    
  }

  hMJpsi->Scale(lumi*time/Nsim);
  hMomentum->Scale(lumi*time/Nsim);
  hThetaPelectron->Scale(lumi*time/Nsim);
  hThetaPpositron->Scale(lumi*time/Nsim);
  hAngleP->Scale(lumi*time/Nsim);
  fraw->Write();

  hMJpsi1->Scale(lumi*time/Nsim);
  hMomentum1->Scale(lumi*time/Nsim);
  hThetaPelectron1->Scale(lumi*time/Nsim);
  hThetaPpositron1->Scale(lumi*time/Nsim);
  hAngleP1->Scale(lumi*time/Nsim);
  fsub->Write();

  return 0;
}
