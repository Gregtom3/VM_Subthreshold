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
    cout << "./photo-solid-study-phi <Nsim>" << endl;
    return 0;
  }

  // Electron beam energy and luminosity
  double Ebeam = 11.0;//GeV
  double lumi = 1.2e37 / 2 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
    
  // Set nuclear
  NUCLEAR::SetNuclear("D");  
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
    
  // Set Phi production model
  PHIMODEL::SetModel("fit");

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 1.0;
  double kmax = 10.99;

  // detected 1
  TFile * fall1 = new TFile("result-photon-phi/Dsolid1.root", "RECREATE");
  TH1D * hMPhi1 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum1 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron1 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron1 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton1 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP1 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP1 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz1 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi1 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi1 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi1->SetDirectory(fall1);
  hMomentum1->SetDirectory(fall1);
  hThetaPelectron1->SetDirectory(fall1);
  hThetaPpositron1->SetDirectory(fall1);
  hThetaPproton1->SetDirectory(fall1);
  hAngleP1->SetDirectory(fall1);
  hFermiP1->SetDirectory(fall1);
  hFermiPz1->SetDirectory(fall1);
  hPPhi1->SetDirectory(fall1);
  hFermiPPPhi1->SetDirectory(fall1);

  // detected subthreshold 1
  TFile * fsub1 = new TFile("result-photon-phi/Dsolidsub1.root", "RECREATE");
  TH1D * hMPhisub1 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub1 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub1 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub1 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub1 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub1 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub1 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub1 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub1 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub1 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub1->SetDirectory(fsub1);
  hMomentumsub1->SetDirectory(fsub1);
  hThetaPelectronsub1->SetDirectory(fsub1);
  hThetaPpositronsub1->SetDirectory(fsub1);
  hThetaPprotonsub1->SetDirectory(fsub1);
  hAnglePsub1->SetDirectory(fsub1);
  hFermiPsub1->SetDirectory(fsub1);
  hFermiPzsub1->SetDirectory(fsub1);
  hPPhisub1->SetDirectory(fsub1);
  hFermiPPPhisub1->SetDirectory(fsub1);

  // detected 2
  TFile * fall2 = new TFile("result-photon-phi/Dsolid2.root", "RECREATE");
  TH1D * hMPhi2 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum2 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron2 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron2 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton2 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP2 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP2 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz2 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi2 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi2 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi2->SetDirectory(fall2);
  hMomentum2->SetDirectory(fall2);
  hThetaPelectron2->SetDirectory(fall2);
  hThetaPpositron2->SetDirectory(fall2);
  hThetaPproton2->SetDirectory(fall2);
  hAngleP2->SetDirectory(fall2);
  hFermiP2->SetDirectory(fall2);
  hFermiPz2->SetDirectory(fall2);
  hPPhi2->SetDirectory(fall2);
  hFermiPPPhi2->SetDirectory(fall2);

  // detected subthreshold 2
  TFile * fsub2 = new TFile("result-photon-phi/Dsolidsub2.root", "RECREATE");
  TH1D * hMPhisub2 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub2 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub2 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub2 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub2 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub2 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub2 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub2 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub2 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub2 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub2->SetDirectory(fsub2);
  hMomentumsub2->SetDirectory(fsub2);
  hThetaPelectronsub2->SetDirectory(fsub2);
  hThetaPpositronsub2->SetDirectory(fsub2);
  hThetaPprotonsub2->SetDirectory(fsub2);
  hAnglePsub2->SetDirectory(fsub2);
  hFermiPsub2->SetDirectory(fsub2);
  hFermiPzsub2->SetDirectory(fsub2);
  hPPhisub2->SetDirectory(fsub2);
  hFermiPPPhisub2->SetDirectory(fsub2);

  // detected 3
  TFile * fall3 = new TFile("result-photon-phi/Dsolid3.root", "RECREATE");
  TH1D * hMPhi3 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum3 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron3 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron3 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton3 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP3 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP3 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz3 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi3 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi3 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi3->SetDirectory(fall3);
  hMomentum3->SetDirectory(fall3);
  hThetaPelectron3->SetDirectory(fall3);
  hThetaPpositron3->SetDirectory(fall3);
  hThetaPproton3->SetDirectory(fall3);
  hAngleP3->SetDirectory(fall3);
  hFermiP3->SetDirectory(fall3);
  hFermiPz3->SetDirectory(fall3);
  hPPhi3->SetDirectory(fall3);
  hFermiPPPhi3->SetDirectory(fall3);

  // detected subthreshold 3
  TFile * fsub3 = new TFile("result-photon-phi/Dsolidsub3.root", "RECREATE");
  TH1D * hMPhisub3 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub3 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub3 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub3 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub3 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub3 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub3 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub3 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub3 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub3 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub3->SetDirectory(fsub3);
  hMomentumsub3->SetDirectory(fsub3);
  hThetaPelectronsub3->SetDirectory(fsub3);
  hThetaPpositronsub3->SetDirectory(fsub3);
  hThetaPprotonsub3->SetDirectory(fsub3);
  hAnglePsub3->SetDirectory(fsub3);
  hFermiPsub3->SetDirectory(fsub3);
  hFermiPzsub3->SetDirectory(fsub3);
  hPPhisub3->SetDirectory(fsub3);
  hFermiPPPhisub3->SetDirectory(fsub3);

  // detected 4
  TFile * fall4 = new TFile("result-photon-phi/Dsolid4.root", "RECREATE");
  TH1D * hMPhi4 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum4 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron4 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron4 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton4 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP4 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP4 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz4 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi4 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi4 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi4->SetDirectory(fall4);
  hMomentum4->SetDirectory(fall4);
  hThetaPelectron4->SetDirectory(fall4);
  hThetaPpositron4->SetDirectory(fall4);
  hThetaPproton4->SetDirectory(fall4);
  hAngleP4->SetDirectory(fall4);
  hFermiP4->SetDirectory(fall4);
  hFermiPz4->SetDirectory(fall4);
  hPPhi4->SetDirectory(fall4);
  hFermiPPPhi4->SetDirectory(fall4);

  // detected subthreshold 4
  TFile * fsub4 = new TFile("result-photon-phi/Dsolidsub4.root", "RECREATE");
  TH1D * hMPhisub4 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub4 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub4 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub4 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub4 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub4 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub4 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub4 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub4 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub4 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub4->SetDirectory(fsub4);
  hMomentumsub4->SetDirectory(fsub4);
  hThetaPelectronsub4->SetDirectory(fsub4);
  hThetaPpositronsub4->SetDirectory(fsub4);
  hThetaPprotonsub4->SetDirectory(fsub4);
  hAnglePsub4->SetDirectory(fsub4);
  hFermiPsub4->SetDirectory(fsub4);
  hFermiPzsub4->SetDirectory(fsub4);
  hPPhisub4->SetDirectory(fsub4);
  hFermiPPPhisub4->SetDirectory(fsub4);

  // detected 5
  TFile * fall5 = new TFile("result-photon-phi/Dsolid5.root", "RECREATE");
  TH1D * hMPhi5 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum5 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron5 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron5 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton5 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP5 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP5 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz5 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi5 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi5 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi5->SetDirectory(fall5);
  hMomentum5->SetDirectory(fall5);
  hThetaPelectron5->SetDirectory(fall5);
  hThetaPpositron5->SetDirectory(fall5);
  hThetaPproton5->SetDirectory(fall5);
  hAngleP5->SetDirectory(fall5);
  hFermiP5->SetDirectory(fall5);
  hFermiPz5->SetDirectory(fall5);
  hPPhi5->SetDirectory(fall5);
  hFermiPPPhi5->SetDirectory(fall5);

  // detected subthreshold 5
  TFile * fsub5 = new TFile("result-photon-phi/Dsolidsub5.root", "RECREATE");
  TH1D * hMPhisub5 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub5 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub5 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub5 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub5 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub5 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub5 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub5 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub5 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub5 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub5->SetDirectory(fsub5);
  hMomentumsub5->SetDirectory(fsub5);
  hThetaPelectronsub5->SetDirectory(fsub5);
  hThetaPpositronsub5->SetDirectory(fsub5);
  hThetaPprotonsub5->SetDirectory(fsub5);
  hAnglePsub5->SetDirectory(fsub5);
  hFermiPsub5->SetDirectory(fsub5);
  hFermiPzsub5->SetDirectory(fsub5);
  hPPhisub5->SetDirectory(fsub5);
  hFermiPPPhisub5->SetDirectory(fsub5);

  // detected 6
  TFile * fall6 = new TFile("result-photon-phi/Dsolid6.root", "RECREATE");
  TH1D * hMPhi6 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentum6 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron6 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron6 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton6 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP6 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP6 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz6 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhi6 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhi6 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhi6->SetDirectory(fall6);
  hMomentum6->SetDirectory(fall6);
  hThetaPelectron6->SetDirectory(fall6);
  hThetaPpositron6->SetDirectory(fall6);
  hThetaPproton6->SetDirectory(fall6);
  hAngleP6->SetDirectory(fall6);
  hFermiP6->SetDirectory(fall6);
  hFermiPz6->SetDirectory(fall6);
  hPPhi6->SetDirectory(fall6);
  hFermiPPPhi6->SetDirectory(fall6);

  // detected subthreshold 6
  TFile * fsub6 = new TFile("result-photon-phi/Dsolidsub6.root", "RECREATE");
  TH1D * hMPhisub6 = new TH1D("Mass_e+e-_Phi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 0.0, 2.0);
  TH2D * hMomentumsub6 = new TH2D("Pe+Pe-_Phi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub6 = new TH2D("ThetaP_e-_Phi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub6 = new TH2D("ThetaP_e+_Phi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub6 = new TH2D("ThetaP_proton_Phi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub6 = new TH2D("AngleP_Phi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub6 = new TH1D("FermiP_Phi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub6 = new TH1D("FermiPz_Phi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPPhisub6 = new TH1D("PPhi", ";P[#phi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPPhisub6 = new TH2D("FermiPPPhi", ";P[N] (GeV);P[#phi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMPhisub6->SetDirectory(fsub6);
  hMomentumsub6->SetDirectory(fsub6);
  hThetaPelectronsub6->SetDirectory(fsub6);
  hThetaPpositronsub6->SetDirectory(fsub6);
  hThetaPprotonsub6->SetDirectory(fsub6);
  hAnglePsub6->SetDirectory(fsub6);
  hFermiPsub6->SetDirectory(fsub6);
  hFermiPzsub6->SetDirectory(fsub6);
  hPPhisub6->SetDirectory(fsub6);
  hFermiPPPhisub6->SetDirectory(fsub6);

  
  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam) * 1.95 / 2;//15cm LD2 target
    weight *= GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_gN2Nee_Phi(ki, kf); 

    if (weight > 0.0){
      //case 1
      if (CheckAcceptance(kf[1], 6.0, 22.0) * CheckAcceptance(kf[2], 6.0, 22.0) * CheckAcceptance(kf[0], 6.0, 22.0)){
	hMPhi1->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum1->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton1->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP1->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP1->Fill( ki[1].P(), weight);
	hFermiPz1->Fill( ki[1].Pz(), weight);
	hPPhi1->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi1->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub1->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub1->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub1->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub1->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub1->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub1->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub1->Fill( ki[1].P(), weight);
	  hFermiPzsub1->Fill( ki[1].Pz(), weight);
	  hPPhisub1->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub1->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }

      //case 2
      if (CheckAcceptance(kf[1], 7.0, 23.0) * CheckAcceptance(kf[2], 7.0, 23.0) * CheckAcceptance(kf[0], 7.0, 23.0)){
	hMPhi2->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum2->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton2->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP2->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP2->Fill( ki[1].P(), weight);
	hFermiPz2->Fill( ki[1].Pz(), weight);
	hPPhi2->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi2->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub2->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub2->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub2->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub2->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub2->Fill( ki[1].P(), weight);
	  hFermiPzsub2->Fill( ki[1].Pz(), weight);
	  hPPhisub2->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub2->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }

      //case 3
      if (CheckAcceptance(kf[1], 8.0, 24.0) * CheckAcceptance(kf[2], 8.0, 24.0) * CheckAcceptance(kf[0], 8.0, 24.0)){
	hMPhi3->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum3->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron3->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron3->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton3->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP3->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP3->Fill( ki[1].P(), weight);
	hFermiPz3->Fill( ki[1].Pz(), weight);
	hPPhi3->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi3->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub3->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub3->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub3->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub3->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub3->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub3->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub3->Fill( ki[1].P(), weight);
	  hFermiPzsub3->Fill( ki[1].Pz(), weight);
	  hPPhisub3->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub3->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }

      //case 4
      if (CheckAcceptance(kf[1], 9.0, 26.0) * CheckAcceptance(kf[2], 9.0, 26.0) * CheckAcceptance(kf[0], 9.0, 26.0)){
	hMPhi4->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum4->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron4->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron4->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton4->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP4->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP4->Fill( ki[1].P(), weight);
	hFermiPz4->Fill( ki[1].Pz(), weight);
	hPPhi4->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi4->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub4->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub4->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub4->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub4->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub4->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub4->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub4->Fill( ki[1].P(), weight);
	  hFermiPzsub4->Fill( ki[1].Pz(), weight);
	  hPPhisub4->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub4->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }

      //case 5
      if (CheckAcceptance(kf[1], 10.0, 28.0) * CheckAcceptance(kf[2], 10.0, 28.0) * CheckAcceptance(kf[0], 10.0, 28.0)){
	hMPhi5->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum5->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron5->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron5->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton5->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP5->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP5->Fill( ki[1].P(), weight);
	hFermiPz5->Fill( ki[1].Pz(), weight);
	hPPhi5->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi5->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub5->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub5->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub5->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub5->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub5->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub5->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub5->Fill( ki[1].P(), weight);
	  hFermiPzsub5->Fill( ki[1].Pz(), weight);
	  hPPhisub5->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub5->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }

      //case 6
      if (CheckAcceptance(kf[1], 11.0, 29.0) * CheckAcceptance(kf[2], 11.0, 29.0) * CheckAcceptance(kf[0], 11.0, 29.0)){
	hMPhi6->Fill( (kf[1]+kf[2]).M(), weight);
	hMomentum6->Fill( kf[1].P(), kf[2].P(), weight);
	hThetaPelectron6->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPpositron6->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPproton6->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hAngleP6->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	hFermiP6->Fill( ki[1].P(), weight);
	hFermiPz6->Fill( ki[1].Pz(), weight);
	hPPhi6->Fill( (kf[1]+kf[2]).P(), weight);
	hFermiPPPhi6->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	
	if (ki[0].E() < 1.57){
	  hMPhisub6->Fill( (kf[1]+kf[2]).M(), weight);
	  hMomentumsub6->Fill( kf[1].P(), kf[2].P(), weight);
	  hThetaPelectronsub6->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPpositronsub6->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hThetaPprotonsub6->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	  hAnglePsub6->Fill( kf[2].Angle(kf[1].Vect())/M_PI*180, kf[1].P(), weight);
	  hFermiPsub6->Fill( ki[1].P(), weight);
	  hFermiPzsub6->Fill( ki[1].Pz(), weight);
	  hPPhisub6->Fill( (kf[1]+kf[2]).P(), weight);
	  hFermiPPPhisub6->Fill( ki[1].P(), (kf[1]+kf[2]).P(), weight);
	}
      }
    }
    
  }

  //1
  hMPhi1->Scale(lumi*time/Nsim);
  hMomentum1->Scale(lumi*time/Nsim);
  hThetaPelectron1->Scale(lumi*time/Nsim);
  hThetaPpositron1->Scale(lumi*time/Nsim);
  hThetaPproton1->Scale(lumi*time/Nsim);
  hAngleP1->Scale(lumi*time/Nsim);
  hFermiP1->Scale(lumi*time/Nsim);
  hFermiPz1->Scale(lumi*time/Nsim);
  hPPhi1->Scale(lumi*time/Nsim);
  hFermiPPPhi1->Scale(lumi*time/Nsim);
  fall1->Write();
  fall1->Close();

  hMPhisub1->Scale(lumi*time/Nsim);
  hMomentumsub1->Scale(lumi*time/Nsim);
  hThetaPelectronsub1->Scale(lumi*time/Nsim);
  hThetaPpositronsub1->Scale(lumi*time/Nsim);
  hThetaPprotonsub1->Scale(lumi*time/Nsim);
  hAnglePsub1->Scale(lumi*time/Nsim);
  hFermiPsub1->Scale(lumi*time/Nsim);
  hFermiPzsub1->Scale(lumi*time/Nsim);
  hPPhisub1->Scale(lumi*time/Nsim);
  hFermiPPPhisub1->Scale(lumi*time/Nsim);
  fsub1->Write();
  fsub1->Close();

  //2
  hMPhi2->Scale(lumi*time/Nsim);
  hMomentum2->Scale(lumi*time/Nsim);
  hThetaPelectron2->Scale(lumi*time/Nsim);
  hThetaPpositron2->Scale(lumi*time/Nsim);
  hThetaPproton2->Scale(lumi*time/Nsim);
  hAngleP2->Scale(lumi*time/Nsim);
  hFermiP2->Scale(lumi*time/Nsim);
  hFermiPz2->Scale(lumi*time/Nsim);
  hPPhi2->Scale(lumi*time/Nsim);
  hFermiPPPhi2->Scale(lumi*time/Nsim);
  fall2->Write();
  fall2->Close();

  hMPhisub2->Scale(lumi*time/Nsim);
  hMomentumsub2->Scale(lumi*time/Nsim);
  hThetaPelectronsub2->Scale(lumi*time/Nsim);
  hThetaPpositronsub2->Scale(lumi*time/Nsim);
  hThetaPprotonsub2->Scale(lumi*time/Nsim);
  hAnglePsub2->Scale(lumi*time/Nsim);
  hFermiPsub2->Scale(lumi*time/Nsim);
  hFermiPzsub2->Scale(lumi*time/Nsim);
  hPPhisub2->Scale(lumi*time/Nsim);
  hFermiPPPhisub2->Scale(lumi*time/Nsim);
  fsub2->Write();
  fsub2->Close();

  //3
  hMPhi3->Scale(lumi*time/Nsim);
  hMomentum3->Scale(lumi*time/Nsim);
  hThetaPelectron3->Scale(lumi*time/Nsim);
  hThetaPpositron3->Scale(lumi*time/Nsim);
  hThetaPproton3->Scale(lumi*time/Nsim);
  hAngleP3->Scale(lumi*time/Nsim);
  hFermiP3->Scale(lumi*time/Nsim);
  hFermiPz3->Scale(lumi*time/Nsim);
  hPPhi3->Scale(lumi*time/Nsim);
  hFermiPPPhi3->Scale(lumi*time/Nsim);
  fall3->Write();
  fall3->Close();

  hMPhisub3->Scale(lumi*time/Nsim);
  hMomentumsub3->Scale(lumi*time/Nsim);
  hThetaPelectronsub3->Scale(lumi*time/Nsim);
  hThetaPpositronsub3->Scale(lumi*time/Nsim);
  hThetaPprotonsub3->Scale(lumi*time/Nsim);
  hAnglePsub3->Scale(lumi*time/Nsim);
  hFermiPsub3->Scale(lumi*time/Nsim);
  hFermiPzsub3->Scale(lumi*time/Nsim);
  hPPhisub3->Scale(lumi*time/Nsim);
  hFermiPPPhisub3->Scale(lumi*time/Nsim);
  fsub3->Write();
  fsub3->Close();

  //4
  hMPhi4->Scale(lumi*time/Nsim);
  hMomentum4->Scale(lumi*time/Nsim);
  hThetaPelectron4->Scale(lumi*time/Nsim);
  hThetaPpositron4->Scale(lumi*time/Nsim);
  hThetaPproton4->Scale(lumi*time/Nsim);
  hAngleP4->Scale(lumi*time/Nsim);
  hFermiP4->Scale(lumi*time/Nsim);
  hFermiPz4->Scale(lumi*time/Nsim);
  hPPhi4->Scale(lumi*time/Nsim);
  hFermiPPPhi4->Scale(lumi*time/Nsim);
  fall4->Write();
  fall4->Close();

  hMPhisub4->Scale(lumi*time/Nsim);
  hMomentumsub4->Scale(lumi*time/Nsim);
  hThetaPelectronsub4->Scale(lumi*time/Nsim);
  hThetaPpositronsub4->Scale(lumi*time/Nsim);
  hThetaPprotonsub4->Scale(lumi*time/Nsim);
  hAnglePsub4->Scale(lumi*time/Nsim);
  hFermiPsub4->Scale(lumi*time/Nsim);
  hFermiPzsub4->Scale(lumi*time/Nsim);
  hPPhisub4->Scale(lumi*time/Nsim);
  hFermiPPPhisub4->Scale(lumi*time/Nsim);
  fsub4->Write();
  fsub4->Close();

  //5
  hMPhi5->Scale(lumi*time/Nsim);
  hMomentum5->Scale(lumi*time/Nsim);
  hThetaPelectron5->Scale(lumi*time/Nsim);
  hThetaPpositron5->Scale(lumi*time/Nsim);
  hThetaPproton5->Scale(lumi*time/Nsim);
  hAngleP5->Scale(lumi*time/Nsim);
  hFermiP5->Scale(lumi*time/Nsim);
  hFermiPz5->Scale(lumi*time/Nsim);
  hPPhi5->Scale(lumi*time/Nsim);
  hFermiPPPhi5->Scale(lumi*time/Nsim);
  fall5->Write();
  fall5->Close();

  hMPhisub5->Scale(lumi*time/Nsim);
  hMomentumsub5->Scale(lumi*time/Nsim);
  hThetaPelectronsub5->Scale(lumi*time/Nsim);
  hThetaPpositronsub5->Scale(lumi*time/Nsim);
  hThetaPprotonsub5->Scale(lumi*time/Nsim);
  hAnglePsub5->Scale(lumi*time/Nsim);
  hFermiPsub5->Scale(lumi*time/Nsim);
  hFermiPzsub5->Scale(lumi*time/Nsim);
  hPPhisub5->Scale(lumi*time/Nsim);
  hFermiPPPhisub5->Scale(lumi*time/Nsim);
  fsub5->Write();
  fsub5->Close();

  //6
  hMPhi6->Scale(lumi*time/Nsim);
  hMomentum6->Scale(lumi*time/Nsim);
  hThetaPelectron6->Scale(lumi*time/Nsim);
  hThetaPpositron6->Scale(lumi*time/Nsim);
  hThetaPproton6->Scale(lumi*time/Nsim);
  hAngleP6->Scale(lumi*time/Nsim);
  hFermiP6->Scale(lumi*time/Nsim);
  hFermiPz6->Scale(lumi*time/Nsim);
  hPPhi6->Scale(lumi*time/Nsim);
  hFermiPPPhi6->Scale(lumi*time/Nsim);
  fall6->Write();
  fall6->Close();

  hMPhisub6->Scale(lumi*time/Nsim);
  hMomentumsub6->Scale(lumi*time/Nsim);
  hThetaPelectronsub6->Scale(lumi*time/Nsim);
  hThetaPpositronsub6->Scale(lumi*time/Nsim);
  hThetaPprotonsub6->Scale(lumi*time/Nsim);
  hAnglePsub6->Scale(lumi*time/Nsim);
  hFermiPsub6->Scale(lumi*time/Nsim);
  hFermiPzsub6->Scale(lumi*time/Nsim);
  hPPhisub6->Scale(lumi*time/Nsim);
  hFermiPPPhisub6->Scale(lumi*time/Nsim);
  fsub6->Write();
  fsub6->Close();

  return 0;
}
