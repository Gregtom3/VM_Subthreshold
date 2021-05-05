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
  Long64_t Nsim = 10000;

  if (argc > 1) Nsim = atoi(argv[1]);
  else {
    cout << "./electro-solid-study <Nsim>" << endl;
    return 0;
  }

  // Electron beam energy and luminosity
  double Ebeam = 11.0;//GeV
  double lumi = 1.2e37 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Set nuclear
  NUCLEAR::SetNuclear("D");
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
  GENERATE::Set_psi_2S();
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set scattered electron range
  double degtorad = M_PI / 180.0;
  GENERATE::cthrange[0] = cos(29.0 * degtorad);
  GENERATE::cthrange[1] = cos(6.0 * degtorad);
  GENERATE::perange[0] = 0.3;//GeV
  GENERATE::perange[1] = 6.0;//GeV

  // detected 1
  TFile * fall1 = new TFile("result-electro-psi2S/Dsolid1.root", "RECREATE");
  TH1D * hMJpsi1 = new TH1D("Mass_ee_Psi2S", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 3.0, 4.0);
  TH2D * hMomentum1 = new TH2D("PePe_Psi2S", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_eminus_Psi2S", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_eplus_Psi2S", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton1 = new TH2D("ThetaP_proton_Psi2S", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP1 = new TH2D("AngleP_Psi2S", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP1 = new TH1D("FermiP_Psi2S", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz1 = new TH1D("FermiPz_Psi2S", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi1 = new TH1D("PPsi2S", ";P[#psi (2S)] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi1 = new TH2D("FermiPPPsi2S", ";P[N] (GeV);P[#psi (2S)] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi1->SetDirectory(fall1);
  hMomentum1->SetDirectory(fall1);
  hThetaPelectron1->SetDirectory(fall1);
  hThetaPpositron1->SetDirectory(fall1);
  hThetaPproton1->SetDirectory(fall1);
  hAngleP1->SetDirectory(fall1);
  hFermiP1->SetDirectory(fall1);
  hFermiPz1->SetDirectory(fall1);
  hPJpsi1->SetDirectory(fall1);
  hFermiPPJpsi1->SetDirectory(fall1);

  // detected subthreshold 1
  TFile * fsub1 = new TFile("result-electro-psi2S/Dsolidsub1.root", "RECREATE");
  TH1D * hMJpsisub1 = new TH1D("Mass_ee_Psi2S", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 3.0, 4.0);
  TH2D * hMomentumsub1 = new TH2D("PePe_Psi2S", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub1 = new TH2D("ThetaP_eminus_Psi2S", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub1 = new TH2D("ThetaP_eplus_Psi2S", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub1 = new TH2D("ThetaP_proton_Psi2S", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub1 = new TH2D("AngleP_Psi2S", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub1 = new TH1D("FermiP_Psi2S", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub1 = new TH1D("FermiPz_Psi2S", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub1 = new TH1D("PPsi2S", ";P[#psi (2S)] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub1 = new TH2D("FermiPPPsi2S", ";P[N] (GeV);P[#psi (2S)] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub1->SetDirectory(fsub1);
  hMomentumsub1->SetDirectory(fsub1);
  hThetaPelectronsub1->SetDirectory(fsub1);
  hThetaPpositronsub1->SetDirectory(fsub1);
  hThetaPprotonsub1->SetDirectory(fsub1);
  hAnglePsub1->SetDirectory(fsub1);
  hFermiPsub1->SetDirectory(fsub1);
  hFermiPzsub1->SetDirectory(fsub1);
  hPJpsisub1->SetDirectory(fsub1);
  hFermiPPJpsisub1->SetDirectory(fsub1);

  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;
  double Mjpsi = 3.68609;
  int count = 0;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight > 0.0){
      q = ki[0] - kf[0];
      //case 1

      if (CheckAcceptance(kf[0], 6.0, 22.0) * CheckAcceptance(kf[2], 6.0, 22.0) * CheckAcceptance(kf[3], 6.0, 22.0)){
	count++;
	hMJpsi1->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum1->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron1->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP1->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP1->Fill( ki[1].P(), weight);
	hFermiPz1->Fill( ki[1].Pz(), weight);
	hPJpsi1->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi1->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub1->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub1->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub1->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub1->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub1->Fill( ki[1].P(), weight);
	  hFermiPzsub1->Fill( ki[1].Pz(), weight);
	  hPJpsisub1->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub1->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
	
      }
    }
  }
  cout << count << endl;
  //1
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
  fall1->Write();
  fall1->Close();

  hMJpsisub1->Scale(lumi*time/Nsim);
  hMomentumsub1->Scale(lumi*time/Nsim);
  hThetaPelectronsub1->Scale(lumi*time/Nsim);
  hThetaPpositronsub1->Scale(lumi*time/Nsim);
  hThetaPprotonsub1->Scale(lumi*time/Nsim);
  hAnglePsub1->Scale(lumi*time/Nsim);
  hFermiPsub1->Scale(lumi*time/Nsim);
  hFermiPzsub1->Scale(lumi*time/Nsim);
  hPJpsisub1->Scale(lumi*time/Nsim);
  hFermiPPJpsisub1->Scale(lumi*time/Nsim);
  fsub1->Write();
  fsub1->Close();
  

  return 0;
}
