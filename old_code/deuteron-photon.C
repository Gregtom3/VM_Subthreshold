#include "Lcore.h"

int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 100000000;

  if (argc > 1) Nsim = atoi(argv[1]);

  // Electron beam energy and luminosity
  double Ebeam = 10.6;//GeV
  double lumi = 2.0e35 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Set Acceptance
  DETECTOR::SetDetector("CLAS12");
  
  // Set nuclear
  NUCLEAR::SetNuclear("D");  
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
    
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 5.0;
  double kmax = 10.59;

  // raw
  TFile * fraw = new TFile("result-photon/Draw.root", "RECREATE");  
  TH1D * hMJpsi1 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum1 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton1 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP1 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP1 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz1 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi1 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi1 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi1->SetDirectory(fraw);
  hMomentum1->SetDirectory(fraw);
  hThetaPelectron1->SetDirectory(fraw);
  hThetaPpositron1->SetDirectory(fraw);
  hThetaPproton1->SetDirectory(fraw);
  hAngleP1->SetDirectory(fraw);
  hFermiP1->SetDirectory(fraw);
  hFermiPz1->SetDirectory(fraw);
  hPJpsi1->SetDirectory(fraw);
  hFermiPPJpsi1->SetDirectory(fraw);
  
  // subthreshold
  TFile * fsub = new TFile("result-photon/Dsub.root", "RECREATE");
  TH1D * hMJpsi2 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum2 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron2 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron2 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton2 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP2 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP2 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz2 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi2 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi2 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi2->SetDirectory(fsub);
  hMomentum2->SetDirectory(fsub);
  hThetaPelectron2->SetDirectory(fsub);
  hThetaPpositron2->SetDirectory(fsub);
  hThetaPproton2->SetDirectory(fsub);
  hAngleP2->SetDirectory(fsub);
  hFermiP2->SetDirectory(fsub);
  hFermiPz2->SetDirectory(fsub);
  hPJpsi2->SetDirectory(fsub);
  hFermiPPJpsi2->SetDirectory(fsub);

  // detected
  TFile * fdetected = new TFile("result-photon/Ddetected.root", "RECREATE");
  TH1D * hMJpsi3 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum3 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron3 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron3 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton3 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP3 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP3 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz3 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi3 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi3 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi3->SetDirectory(fdetected);
  hMomentum3->SetDirectory(fdetected);
  hThetaPelectron3->SetDirectory(fdetected);
  hThetaPpositron3->SetDirectory(fdetected);
  hThetaPproton3->SetDirectory(fdetected);
  hAngleP3->SetDirectory(fdetected);
  hFermiP3->SetDirectory(fdetected);
  hFermiPz3->SetDirectory(fdetected);
  hPJpsi3->SetDirectory(fdetected);
  hFermiPPJpsi3->SetDirectory(fdetected);

  // detected subthreshold
  TFile * fdetectedsub = new TFile("result-photon/Ddetectedsub.root", "RECREATE");
  TH1D * hMJpsi4 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum4 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron4 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron4 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton4 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP4 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP4 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz4 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi4 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi4 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi4->SetDirectory(fdetectedsub);
  hMomentum4->SetDirectory(fdetectedsub);
  hThetaPelectron4->SetDirectory(fdetectedsub);
  hThetaPpositron4->SetDirectory(fdetectedsub);
  hThetaPproton4->SetDirectory(fdetectedsub);
  hAngleP4->SetDirectory(fdetectedsub);
  hFermiP4->SetDirectory(fdetectedsub);
  hFermiPz4->SetDirectory(fdetectedsub);
  hPJpsi4->SetDirectory(fdetectedsub);
  hFermiPPJpsi4->SetDirectory(fdetectedsub);
  
  
  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam);
    weight *= GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_gN2Nee_Jpsi(ki, kf); 

    if (weight > 0.0){
      acceptance = DETECTOR::AcceptanceCLAS12(kf[1], "e+") * DETECTOR::AcceptanceCLAS12(kf[2], "e-");
      
      hMJpsi1->Fill( (kf[1]+kf[2]).M(), weight);
      hMomentum1->Fill( kf[1].P(), kf[2].P(), weight);
      hThetaPelectron1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
      hThetaPpositron1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
      hThetaPproton1->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
      hAngleP1->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
      hFermiP1->Fill( ki[1].P(), weight);
      hFermiPz1->Fill( ki[1].Pz(), weight);
      hPJpsi1->Fill( (kf[1]+kf[2]).P(), weight);
      hFermiPPJpsi1->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);

      hMJpsi3->Fill( (kf[1]+kf[2]).M(), weight * acceptance);
      hMomentum3->Fill( kf[1].P(), kf[2].P(), weight * acceptance);
      hThetaPelectron3->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight * acceptance);
      hThetaPpositron3->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight * acceptance);
      hThetaPproton3->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight * acceptance);
      hAngleP3->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight * acceptance);
      hFermiP3->Fill( ki[1].P(), weight * acceptance);
      hFermiPz3->Fill( ki[1].Pz(), weight * acceptance);
      hPJpsi3->Fill( (kf[1]+kf[2]).P(), weight * acceptance);
      hFermiPPJpsi3->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight * acceptance);
      
      if (ki[0].E() < 8.2){
	hMJpsi2->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum2->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton2->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP2->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP2->Fill( ki[1].P(), weight);
	hFermiPz2->Fill( ki[1].Pz(), weight);
	hPJpsi2->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPJpsi2->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);

	hMJpsi4->Fill( (kf[1]+kf[2]).M(), weight * acceptance);
	hMomentum4->Fill( kf[1].P(), kf[2].P(), weight * acceptance);
	hThetaPelectron4->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight * acceptance);
	hThetaPpositron4->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight * acceptance);
	hThetaPproton4->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight * acceptance);
	hAngleP4->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight * acceptance);
	hFermiP4->Fill( ki[1].P(), weight * acceptance);
	hFermiPz4->Fill( ki[1].Pz(), weight * acceptance);
	hPJpsi4->Fill( (kf[1]+kf[2]).P(), weight * acceptance);
	hFermiPPJpsi4->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight * acceptance);
      }
    }
    
  }

  hMJpsi1->Scale(lumi*time/Nsim);
  hMomentum1->Scale(lumi*time/Nsim);
  hThetaPelectron1->Scale(lumi*time/Nsim);
  hThetaPpositron1->Scale(lumi*time/Nsim);
  hThetaPproton1->Scale(lumi*time/Nsim);
  hAngleP1->Scale(lumi*time/Nsim);
  hFermiP1->Scale(lumi*time/Nsim);
  hFermiPz1->Scale(lumi*time/Nsim);
  hPJpsi1->Scale(lumi*time/Nsim);
  hFermiPPJpsi1->Scale(lumi*time/Nsim);
  fraw->Write();
  fraw->Close();

  hMJpsi2->Scale(lumi*time/Nsim);
  hMomentum2->Scale(lumi*time/Nsim);
  hThetaPelectron2->Scale(lumi*time/Nsim);
  hThetaPpositron2->Scale(lumi*time/Nsim);
  hThetaPproton2->Scale(lumi*time/Nsim);
  hAngleP2->Scale(lumi*time/Nsim);
  hFermiP2->Scale(lumi*time/Nsim);
  hFermiPz2->Scale(lumi*time/Nsim);
  hPJpsi2->Scale(lumi*time/Nsim);
  hFermiPPJpsi2->Scale(lumi*time/Nsim);
  fsub->Write();
  fsub->Close();

  hMJpsi3->Scale(lumi*time/Nsim);
  hMomentum3->Scale(lumi*time/Nsim);
  hThetaPelectron3->Scale(lumi*time/Nsim);
  hThetaPpositron3->Scale(lumi*time/Nsim);
  hThetaPproton3->Scale(lumi*time/Nsim);
  hAngleP3->Scale(lumi*time/Nsim);
  hFermiP3->Scale(lumi*time/Nsim);
  hFermiPz3->Scale(lumi*time/Nsim);
  hPJpsi3->Scale(lumi*time/Nsim);
  hFermiPPJpsi3->Scale(lumi*time/Nsim);
  fdetected->Write();
  fdetected->Close();

  hMJpsi4->Scale(lumi*time/Nsim);
  hMomentum4->Scale(lumi*time/Nsim);
  hThetaPelectron4->Scale(lumi*time/Nsim);
  hThetaPpositron4->Scale(lumi*time/Nsim);
  hThetaPproton4->Scale(lumi*time/Nsim);
  hAngleP4->Scale(lumi*time/Nsim);
  hFermiP4->Scale(lumi*time/Nsim);
  hFermiPz4->Scale(lumi*time/Nsim);
  hPJpsi4->Scale(lumi*time/Nsim);
  hFermiPPJpsi4->Scale(lumi*time/Nsim);
  fdetectedsub->Write();
  fdetectedsub->Close();

  

  return 0;
}
