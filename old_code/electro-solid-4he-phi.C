#include "Lcore.h"

int CheckAcceptance(const TLorentzVector P, const double thmin, const double thmax){
  if (P.Theta() < thmin * M_PI / 180.0) return 0;
  if (P.Theta() > thmax * M_PI / 180.0) return 0;
  if (P.P() < 0.3) return 0;
  return 1;
}

int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 100000000;

  if (argc > 1) Nsim = atoi(argv[1]);
  else {
    cout << "./electro-solid-study-phi <Nsim>" << endl;
    return 0;
  }

  // Electron beam energy and luminosity
  double Ebeam = 2.2;//GeV
  double lumi = 1.2e37 / 2 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Detector
  DETECTOR::SetDetector("CLAS12");
  
  // Set phi production model
  PHIHE4::SetModel();
  
  // Set bremsstrahlung photon
  double degtorad = M_PI / 180.0;
  GENERATE::cthrange[0] = cos(4.5 * degtorad);
  GENERATE::cthrange[1] = cos(2.5 * degtorad);
  GENERATE::perange[0] = 0.5;//GeV
  GENERATE::perange[1] = 0.9;//GeV

  TFile * fraw = new TFile("result-electro-phi/He4raw.root", "RECREATE");
  TH1D * hMPhi1 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum1 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * hThetaPphi1 = new TH2D("ThetaPPhi", ";#theta[#phi] (deg);P[#phi] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  hMPhi1->SetDirectory(fraw);
  hMomentum1->SetDirectory(fraw);
  hThetaPelectron1->SetDirectory(fraw);
  hThetaPpositron1->SetDirectory(fraw);
  hThetaPphi1->SetDirectory(fraw);
 
  // detected subthreshold 1
  TFile * fsub = new TFile("result-electro-phi/He4sub.root", "RECREATE");
  TH1D * hMPhi2 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum2 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * hThetaPelectron2 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * hThetaPpositron2 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * hThetaPphi2 = new TH2D("ThetaPPhi", ";#theta[#phi] (deg);P[#phi] (GeV)", 90, 0.0, 180.0, 100, 0.0, 2.0);
  hMPhi2->SetDirectory(fsub);
  hMomentum2->SetDirectory(fsub);
  hThetaPelectron2->SetDirectory(fsub);
  hThetaPpositron2->SetDirectory(fsub);
  hThetaPphi2->SetDirectory(fsub);
  
  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;
  double Mphi = 1.019;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_e4He2eee_Phi(ki, kf);

    weight *= (DETECTOR::AcceptanceCLAS12FD(kf[1], "e+") + DETECTOR::AcceptanceCLAS12CD(kf[1], "e+"))
      * (DETECTOR::AcceptanceCLAS12FD(kf[2], "e-") + DETECTOR::AcceptanceCLAS12CD(kf[2], "e-"));

    if (weight > 0.0){
      q = ki[0] - kf[0];

      hMPhi1->Fill( (kf[1]+kf[2]).M(), weight);
      hMomentum1->Fill( kf[1].P(), kf[2].P(), weight);
      hThetaPelectron1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
      hThetaPpositron1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
      hThetaPphi1->Fill( (kf[1]+kf[2]).Theta()/M_PI*180, (kf[1]+kf[2]).P(), weight);
	
      if (q.E() < Mphi + (Mphi * Mphi - q * q) / (2.0 * Mp)){
	hMPhi2->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum2->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPphi2->Fill( (kf[1]+kf[2]).Theta()/M_PI*180, (kf[1]+kf[2]).P(), weight);      
      }
    }
    
  }

  hMPhi1->Scale(lumi*time/Nsim);
  hMomentum1->Scale(lumi*time/Nsim);
  hThetaPelectron1->Scale(lumi*time/Nsim);
  hThetaPpositron1->Scale(lumi*time/Nsim);
  hThetaPphi1->Scale(lumi*time/Nsim);
  fraw->Write();
  fraw->Close();

  hMPhi2->Scale(lumi*time/Nsim);
  hMomentum2->Scale(lumi*time/Nsim);
  hThetaPelectron2->Scale(lumi*time/Nsim);
  hThetaPpositron2->Scale(lumi*time/Nsim);
  hThetaPphi2->Scale(lumi*time/Nsim);
  fsub->Write();
  fsub->Close();

  return 0;
}
