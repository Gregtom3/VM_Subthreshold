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
  TFile * fall1 = new TFile("result-electro/Dsolid1.root", "RECREATE");
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
  TFile * fsub1 = new TFile("result-electro/Dsolidsub1.root", "RECREATE");
  TH1D * hMJpsisub1 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub1 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub1 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub1 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub1 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub1 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub1 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub1 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub1 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub1 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
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

  // detected 2
  TFile * fall2 = new TFile("result-electro/Dsolid2.root", "RECREATE");
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
  hMJpsi2->SetDirectory(fall2);
  hMomentum2->SetDirectory(fall2);
  hThetaPelectron2->SetDirectory(fall2);
  hThetaPpositron2->SetDirectory(fall2);
  hThetaPproton2->SetDirectory(fall2);
  hAngleP2->SetDirectory(fall2);
  hFermiP2->SetDirectory(fall2);
  hFermiPz2->SetDirectory(fall2);
  hPJpsi2->SetDirectory(fall2);
  hFermiPPJpsi2->SetDirectory(fall2);

  // detected subthreshold 2
  TFile * fsub2 = new TFile("result-electro/Dsolidsub2.root", "RECREATE");
  TH1D * hMJpsisub2 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub2 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub2 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub2 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub2 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub2 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub2 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub2 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub2 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub2 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub2->SetDirectory(fsub2);
  hMomentumsub2->SetDirectory(fsub2);
  hThetaPelectronsub2->SetDirectory(fsub2);
  hThetaPpositronsub2->SetDirectory(fsub2);
  hThetaPprotonsub2->SetDirectory(fsub2);
  hAnglePsub2->SetDirectory(fsub2);
  hFermiPsub2->SetDirectory(fsub2);
  hFermiPzsub2->SetDirectory(fsub2);
  hPJpsisub2->SetDirectory(fsub2);
  hFermiPPJpsisub2->SetDirectory(fsub2);

  // detected 3
  TFile * fall3 = new TFile("result-electro/Dsolid3.root", "RECREATE");
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
  hMJpsi3->SetDirectory(fall3);
  hMomentum3->SetDirectory(fall3);
  hThetaPelectron3->SetDirectory(fall3);
  hThetaPpositron3->SetDirectory(fall3);
  hThetaPproton3->SetDirectory(fall3);
  hAngleP3->SetDirectory(fall3);
  hFermiP3->SetDirectory(fall3);
  hFermiPz3->SetDirectory(fall3);
  hPJpsi3->SetDirectory(fall3);
  hFermiPPJpsi3->SetDirectory(fall3);

  // detected subthreshold 3
  TFile * fsub3 = new TFile("result-electro/Dsolidsub3.root", "RECREATE");
  TH1D * hMJpsisub3 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub3 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub3 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub3 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub3 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub3 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub3 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub3 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub3 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub3 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub3->SetDirectory(fsub3);
  hMomentumsub3->SetDirectory(fsub3);
  hThetaPelectronsub3->SetDirectory(fsub3);
  hThetaPpositronsub3->SetDirectory(fsub3);
  hThetaPprotonsub3->SetDirectory(fsub3);
  hAnglePsub3->SetDirectory(fsub3);
  hFermiPsub3->SetDirectory(fsub3);
  hFermiPzsub3->SetDirectory(fsub3);
  hPJpsisub3->SetDirectory(fsub3);
  hFermiPPJpsisub3->SetDirectory(fsub3);

  // detected 4
  TFile * fall4 = new TFile("result-electro/Dsolid4.root", "RECREATE");
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
  hMJpsi4->SetDirectory(fall4);
  hMomentum4->SetDirectory(fall4);
  hThetaPelectron4->SetDirectory(fall4);
  hThetaPpositron4->SetDirectory(fall4);
  hThetaPproton4->SetDirectory(fall4);
  hAngleP4->SetDirectory(fall4);
  hFermiP4->SetDirectory(fall4);
  hFermiPz4->SetDirectory(fall4);
  hPJpsi4->SetDirectory(fall4);
  hFermiPPJpsi4->SetDirectory(fall4);

  // detected subthreshold 4
  TFile * fsub4 = new TFile("result-electro/Dsolidsub4.root", "RECREATE");
  TH1D * hMJpsisub4 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub4 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub4 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub4 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub4 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub4 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub4 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub4 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub4 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub4 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub4->SetDirectory(fsub4);
  hMomentumsub4->SetDirectory(fsub4);
  hThetaPelectronsub4->SetDirectory(fsub4);
  hThetaPpositronsub4->SetDirectory(fsub4);
  hThetaPprotonsub4->SetDirectory(fsub4);
  hAnglePsub4->SetDirectory(fsub4);
  hFermiPsub4->SetDirectory(fsub4);
  hFermiPzsub4->SetDirectory(fsub4);
  hPJpsisub4->SetDirectory(fsub4);
  hFermiPPJpsisub4->SetDirectory(fsub4);

  // detected 5
  TFile * fall5 = new TFile("result-electro/Dsolid5.root", "RECREATE");
  TH1D * hMJpsi5 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum5 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron5 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron5 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton5 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP5 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP5 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz5 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi5 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi5 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi5->SetDirectory(fall5);
  hMomentum5->SetDirectory(fall5);
  hThetaPelectron5->SetDirectory(fall5);
  hThetaPpositron5->SetDirectory(fall5);
  hThetaPproton5->SetDirectory(fall5);
  hAngleP5->SetDirectory(fall5);
  hFermiP5->SetDirectory(fall5);
  hFermiPz5->SetDirectory(fall5);
  hPJpsi5->SetDirectory(fall5);
  hFermiPPJpsi5->SetDirectory(fall5);

  // detected subthreshold 5
  TFile * fsub5 = new TFile("result-electro/Dsolidsub5.root", "RECREATE");
  TH1D * hMJpsisub5 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub5 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub5 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub5 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub5 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub5 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub5 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub5 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub5 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub5 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub5->SetDirectory(fsub5);
  hMomentumsub5->SetDirectory(fsub5);
  hThetaPelectronsub5->SetDirectory(fsub5);
  hThetaPpositronsub5->SetDirectory(fsub5);
  hThetaPprotonsub5->SetDirectory(fsub5);
  hAnglePsub5->SetDirectory(fsub5);
  hFermiPsub5->SetDirectory(fsub5);
  hFermiPzsub5->SetDirectory(fsub5);
  hPJpsisub5->SetDirectory(fsub5);
  hFermiPPJpsisub5->SetDirectory(fsub5);

  // detected 6
  TFile * fall6 = new TFile("result-electro/Dsolid6.root", "RECREATE");
  TH1D * hMJpsi6 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum6 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron6 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron6 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPproton6 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAngleP6 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiP6 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPz6 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsi6 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsi6 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsi6->SetDirectory(fall6);
  hMomentum6->SetDirectory(fall6);
  hThetaPelectron6->SetDirectory(fall6);
  hThetaPpositron6->SetDirectory(fall6);
  hThetaPproton6->SetDirectory(fall6);
  hAngleP6->SetDirectory(fall6);
  hFermiP6->SetDirectory(fall6);
  hFermiPz6->SetDirectory(fall6);
  hPJpsi6->SetDirectory(fall6);
  hFermiPPJpsi6->SetDirectory(fall6);

  // detected subthreshold 6
  TFile * fsub6 = new TFile("result-electro/Dsolidsub6.root", "RECREATE");
  TH1D * hMJpsisub6 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentumsub6 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectronsub6 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositronsub6 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPprotonsub6 = new TH2D("ThetaP_proton_Jpsi", ";#theta[p] (deg);P[p] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hAnglePsub6 = new TH2D("AngleP_Jpsi", ";#alpha<e^{+}e^{-}> (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hFermiPsub6 = new TH1D("FermiP_Jpsi", ";P[N] (GeV);Events / hour", 100, 0.0, 1.0);
  TH1D * hFermiPzsub6 = new TH1D("FermiPz_Jpsi", ";Pz[N] (GeV);Events / hour", 100, -1.0, 1.0);
  TH1D * hPJpsisub6 = new TH1D("PJpsi", ";P[J/#psi] (GeV);Events / hour", 100, 0.0, 10.0);
  TH2D * hFermiPPJpsisub6 = new TH2D("FermiPPJpsi", ";P[N] (GeV);P[J/#psi] (GeV)", 100, 0.0, 1.0, 100, 0.0, 10.0);
  hMJpsisub6->SetDirectory(fsub6);
  hMomentumsub6->SetDirectory(fsub6);
  hThetaPelectronsub6->SetDirectory(fsub6);
  hThetaPpositronsub6->SetDirectory(fsub6);
  hThetaPprotonsub6->SetDirectory(fsub6);
  hAnglePsub6->SetDirectory(fsub6);
  hFermiPsub6->SetDirectory(fsub6);
  hFermiPzsub6->SetDirectory(fsub6);
  hPJpsisub6->SetDirectory(fsub6);
  hFermiPPJpsisub6->SetDirectory(fsub6);

  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;
  double Mjpsi = 3.68609;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight > 0.0){
      q = ki[0] - kf[0];
      //case 1
      if (CheckAcceptance(kf[0], 6.0, 22.0) * CheckAcceptance(kf[2], 6.0, 22.0) * CheckAcceptance(kf[3], 6.0, 22.0)){
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
		
      //case 2
      if (CheckAcceptance(kf[0], 7.0, 23.0) * CheckAcceptance(kf[2], 7.0, 23.0) * CheckAcceptance(kf[3], 7.0, 23.0)){
	hMJpsi2->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum2->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron2->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP2->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP2->Fill( ki[1].P(), weight);
	hFermiPz2->Fill( ki[1].Pz(), weight);
	hPJpsi2->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi2->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub2->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub2->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub2->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub2->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub2->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub2->Fill( ki[1].P(), weight);
	  hFermiPzsub2->Fill( ki[1].Pz(), weight);
	  hPJpsisub2->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub2->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
      }

      //case 3
      if (CheckAcceptance(kf[0], 8.0, 24.0) * CheckAcceptance(kf[2], 8.0, 24.0) * CheckAcceptance(kf[3], 8.0, 24.0)){
	hMJpsi3->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum3->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron3->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron3->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton3->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP3->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP3->Fill( ki[1].P(), weight);
	hFermiPz3->Fill( ki[1].Pz(), weight);
	hPJpsi3->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi3->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub3->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub3->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub3->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub3->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub3->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub3->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub3->Fill( ki[1].P(), weight);
	  hFermiPzsub3->Fill( ki[1].Pz(), weight);
	  hPJpsisub3->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub3->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
      }

      //case 4
      if (CheckAcceptance(kf[0], 9.0, 26.0) * CheckAcceptance(kf[2], 9.0, 26.0) * CheckAcceptance(kf[3], 9.0, 26.0)){
	hMJpsi4->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum4->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron4->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron4->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton4->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP4->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP4->Fill( ki[1].P(), weight);
	hFermiPz4->Fill( ki[1].Pz(), weight);
	hPJpsi4->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi4->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub4->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub4->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub4->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub4->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub4->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub4->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub4->Fill( ki[1].P(), weight);
	  hFermiPzsub4->Fill( ki[1].Pz(), weight);
	  hPJpsisub4->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub4->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
      }

      //case 5
      if (CheckAcceptance(kf[0], 10.0, 28.0) * CheckAcceptance(kf[2], 10.0, 28.0) * CheckAcceptance(kf[3], 10.0, 28.0)){
	hMJpsi5->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum5->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron5->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron5->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton5->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP5->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP5->Fill( ki[1].P(), weight);
	hFermiPz5->Fill( ki[1].Pz(), weight);
	hPJpsi5->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi5->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub5->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub5->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub5->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub5->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub5->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub5->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub5->Fill( ki[1].P(), weight);
	  hFermiPzsub5->Fill( ki[1].Pz(), weight);
	  hPJpsisub5->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub5->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
      }

      //case 6
      if (CheckAcceptance(kf[0], 11.0, 29.0) * CheckAcceptance(kf[2], 11.0, 29.0) * CheckAcceptance(kf[3], 11.0, 29.0)){
	hMJpsi6->Fill( (kf[2]+kf[3]).M(), weight);
	hMomentum6->Fill( kf[2].P(), kf[3].P(), weight);
	hThetaPelectron6->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	hThetaPpositron6->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	hThetaPproton6->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hAngleP6->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	hFermiP6->Fill( ki[1].P(), weight);
	hFermiPz6->Fill( ki[1].Pz(), weight);
	hPJpsi6->Fill( (kf[2]+kf[3]).P(), weight);
	hFermiPPJpsi6->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	
	if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	  hMJpsisub6->Fill( (kf[2]+kf[3]).M(), weight);
	  hMomentumsub6->Fill( kf[2].P(), kf[3].P(), weight);
	  hThetaPelectronsub6->Fill( kf[3].Theta()/M_PI*180, kf[3].P(), weight);
	  hThetaPpositronsub6->Fill( kf[2].Theta()/M_PI*180, kf[2].P(), weight);
	  hThetaPprotonsub6->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	  hAnglePsub6->Fill( kf[3].Angle(kf[2].Vect())/M_PI*180, kf[2].P(), weight);
	  hFermiPsub6->Fill( ki[1].P(), weight);
	  hFermiPzsub6->Fill( ki[1].Pz(), weight);
	  hPJpsisub6->Fill( (kf[2]+kf[3]).P(), weight);
	  hFermiPPJpsisub6->Fill( ki[1].P(), (kf[2]+kf[3]).P(), weight);
	}
      }
    }
  }

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

  //2
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
  fall2->Write();
  fall2->Close();

  hMJpsisub2->Scale(lumi*time/Nsim);
  hMomentumsub2->Scale(lumi*time/Nsim);
  hThetaPelectronsub2->Scale(lumi*time/Nsim);
  hThetaPpositronsub2->Scale(lumi*time/Nsim);
  hThetaPprotonsub2->Scale(lumi*time/Nsim);
  hAnglePsub2->Scale(lumi*time/Nsim);
  hFermiPsub2->Scale(lumi*time/Nsim);
  hFermiPzsub2->Scale(lumi*time/Nsim);
  hPJpsisub2->Scale(lumi*time/Nsim);
  hFermiPPJpsisub2->Scale(lumi*time/Nsim);
  fsub2->Write();
  fsub2->Close();

  //3
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
  fall3->Write();
  fall3->Close();

  hMJpsisub3->Scale(lumi*time/Nsim);
  hMomentumsub3->Scale(lumi*time/Nsim);
  hThetaPelectronsub3->Scale(lumi*time/Nsim);
  hThetaPpositronsub3->Scale(lumi*time/Nsim);
  hThetaPprotonsub3->Scale(lumi*time/Nsim);
  hAnglePsub3->Scale(lumi*time/Nsim);
  hFermiPsub3->Scale(lumi*time/Nsim);
  hFermiPzsub3->Scale(lumi*time/Nsim);
  hPJpsisub3->Scale(lumi*time/Nsim);
  hFermiPPJpsisub3->Scale(lumi*time/Nsim);
  fsub3->Write();
  fsub3->Close();

  //4
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
  fall4->Write();
  fall4->Close();

  hMJpsisub4->Scale(lumi*time/Nsim);
  hMomentumsub4->Scale(lumi*time/Nsim);
  hThetaPelectronsub4->Scale(lumi*time/Nsim);
  hThetaPpositronsub4->Scale(lumi*time/Nsim);
  hThetaPprotonsub4->Scale(lumi*time/Nsim);
  hAnglePsub4->Scale(lumi*time/Nsim);
  hFermiPsub4->Scale(lumi*time/Nsim);
  hFermiPzsub4->Scale(lumi*time/Nsim);
  hPJpsisub4->Scale(lumi*time/Nsim);
  hFermiPPJpsisub4->Scale(lumi*time/Nsim);
  fsub4->Write();
  fsub4->Close();

  //5
  hMJpsi5->Scale(lumi*time/Nsim);
  hMomentum5->Scale(lumi*time/Nsim);
  hThetaPelectron5->Scale(lumi*time/Nsim);
  hThetaPpositron5->Scale(lumi*time/Nsim);
  hThetaPproton5->Scale(lumi*time/Nsim);
  hAngleP5->Scale(lumi*time/Nsim);
  hFermiP5->Scale(lumi*time/Nsim);
  hFermiPz5->Scale(lumi*time/Nsim);
  hPJpsi5->Scale(lumi*time/Nsim);
  hFermiPPJpsi5->Scale(lumi*time/Nsim);
  fall5->Write();
  fall5->Close();

  hMJpsisub5->Scale(lumi*time/Nsim);
  hMomentumsub5->Scale(lumi*time/Nsim);
  hThetaPelectronsub5->Scale(lumi*time/Nsim);
  hThetaPpositronsub5->Scale(lumi*time/Nsim);
  hThetaPprotonsub5->Scale(lumi*time/Nsim);
  hAnglePsub5->Scale(lumi*time/Nsim);
  hFermiPsub5->Scale(lumi*time/Nsim);
  hFermiPzsub5->Scale(lumi*time/Nsim);
  hPJpsisub5->Scale(lumi*time/Nsim);
  hFermiPPJpsisub5->Scale(lumi*time/Nsim);
  fsub5->Write();
  fsub5->Close();

  //6
  hMJpsi6->Scale(lumi*time/Nsim);
  hMomentum6->Scale(lumi*time/Nsim);
  hThetaPelectron6->Scale(lumi*time/Nsim);
  hThetaPpositron6->Scale(lumi*time/Nsim);
  hThetaPproton6->Scale(lumi*time/Nsim);
  hAngleP6->Scale(lumi*time/Nsim);
  hFermiP6->Scale(lumi*time/Nsim);
  hFermiPz6->Scale(lumi*time/Nsim);
  hPJpsi6->Scale(lumi*time/Nsim);
  hFermiPPJpsi6->Scale(lumi*time/Nsim);
  fall6->Write();
  fall6->Close();

  hMJpsisub6->Scale(lumi*time/Nsim);
  hMomentumsub6->Scale(lumi*time/Nsim);
  hThetaPelectronsub6->Scale(lumi*time/Nsim);
  hThetaPpositronsub6->Scale(lumi*time/Nsim);
  hThetaPprotonsub6->Scale(lumi*time/Nsim);
  hAnglePsub6->Scale(lumi*time/Nsim);
  hFermiPsub6->Scale(lumi*time/Nsim);
  hFermiPzsub6->Scale(lumi*time/Nsim);
  hPJpsisub6->Scale(lumi*time/Nsim);
  hFermiPPJpsisub6->Scale(lumi*time/Nsim);
  fsub6->Write();
  fsub6->Close();
  

  return 0;
}
