#include "Lcore.h"

int MisPID(TLorentzVector * k);
int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km);

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./twopinew <pemin> <savepath> <filename> <Nsim>" << endl;
    return 0;
  }

  const double pemin = atof(argv[1]);
  const double pemax = 4.0;

  TString path = argv[2];
  TString filename = path + "twopinew.root";
  
  double Nsim = atof(argv[4]);

  Initialize();

  double convert = 1e-4 / pow(Phys::hbar, 2) * 79.0;//convert unit to GeV^-2 per gold

  TLorentzVector ki[2], kf[5];
  TLorentzVector kgold;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, Mp * 197.0);

  TFile * fs = new TFile(filename.Data(), "RECREATE");

  TH1D * h0 = new TH1D("MpPiPi", "", 2200, 1.18, 1.62);
  TH1D * h1 = new TH1D("MpKK_PiPi", "", 2200, 1.88, 2.32);

  TH1D * h0a = new TH1D("MpPip", "", 2200, 1.03, 1.47);
  TH1D * h1a = new TH1D("MpKp_PiPi", "", 2200, 1.38, 1.82);

  TH1D * h0b = new TH1D("MpPim", "", 2200, 1.03, 1.47);
  TH1D * h1b = new TH1D("MpKm_PiPi", "", 2200, 1.38, 1.82);

  TH1D * h0c = new TH1D("MPiPi", "", 2200, 0.26, 1.70);
  TH1D * h1c = new TH1D("MKK_PiPi", "", 2200, 0.98, 1.42);

  TH1D * h0d = new TH1D("Missing", "", 2200, 183.6, 186.0);
  TH1D * h1d = new TH1D("Missing_PiPi", "", 2200, 183.6, 186.0);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);

  h0a->SetDirectory(fs);
  h1a->SetDirectory(fs);

  h0b->SetDirectory(fs);
  h1b->SetDirectory(fs);

  h0c->SetDirectory(fs);
  h1c->SetDirectory(fs);

  h0d->SetDirectory(fs);
  h1d->SetDirectory(fs);

  TH2D * d0a = new TH2D("Momentum_p_Pip", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0b = new TH2D("Momentum_p_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0c = new TH2D("Momentum_Pip_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);

  d0a->SetDirectory(fs);
  d0b->SetDirectory(fs);
  d0c->SetDirectory(fs);

  TH2D * d1a = new TH2D("Momentum_p_Pip_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1b = new TH2D("Momentum_p_Pim_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1c = new TH2D("Momentum_Pip_Pim_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);

  d1a->SetDirectory(fs);
  d1b->SetDirectory(fs);
  d1c->SetDirectory(fs);

  TLorentzVector PP, Pa, Pb, Pc, Pd;

  ifstream infile(argv[3]);
  double temp, px, py, pz, E, weight;
  while (infile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> weight){
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[0].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[1].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[2].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[3].SetXYZT(px, py, pz, E);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      if (kf[0].Theta() / M_PI * 180.0 > 4.5) continue;
      if (kf[0].Theta() / M_PI * 180.0 < 2.5) continue;
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "pi+") * DETECTOR::Smear(&kf[3], "pi-");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      Pd = ki[0] + kgold - kf[0] - PP;
      h0->Fill(PP.M(), weight);
      h0a->Fill(Pa.M(), weight);
      h0b->Fill(Pb.M(), weight);
      h0c->Fill(Pc.M(), weight);
      h0d->Fill(Pd.M(), weight);
      d0a->Fill(kf[2].P(), kf[1].P(), weight);
      d0b->Fill(kf[3].P(), kf[1].P(), weight);
      d0c->Fill(kf[2].P(), kf[3].P(), weight);
      MisPID(&kf[2]);
      MisPID(&kf[3]);
      weight *= Cut(kf[1], kf[2], kf[3]);
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      Pd = ki[0] + kgold - kf[0] - PP;
      h1->Fill(PP.M(), weight);
      h1a->Fill(Pa.M(), weight);
      h1b->Fill(Pb.M(), weight);
      h1c->Fill(Pc.M(), weight);
      h1d->Fill(Pd.M(), weight);
      d1a->Fill(kf[2].P(), kf[1].P(), weight);
      d1b->Fill(kf[3].P(), kf[1].P(), weight);
      d1c->Fill(kf[2].P(), kf[3].P(), weight);
    }
  }

  h0->Scale(convert/Nsim);
  h1->Scale(convert/Nsim);
 
  h0a->Scale(convert/Nsim);
  h1a->Scale(convert/Nsim);
 
  h0b->Scale(convert/Nsim);
  h1b->Scale(convert/Nsim);
 
  h0c->Scale(convert/Nsim);
  h1c->Scale(convert/Nsim);

  h0d->Scale(convert/Nsim);
  h1d->Scale(convert/Nsim);
 
  d0a->Scale(convert/Nsim);
  d0b->Scale(convert/Nsim);
  d0c->Scale(convert/Nsim);

  d1a->Scale(convert/Nsim);
  d1b->Scale(convert/Nsim);
  d1c->Scale(convert/Nsim);
 
  infile.close();
  fs->Write();

  return 0;
}

int MisPID(TLorentzVector * k){
  k->SetE(sqrt(PARTICLE::K.M() * PARTICLE::K.M() + k->P() * k->P()));
  return 0;
}

int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km){
  TLorentzVector kk = p + Kp;
  if (kk.M() > 1.48) return 0;
  kk = p + Km;
  if (kk.M() > 1.48) return 0;
  kk = Kp + Km;
  if (kk.M() > 1.04) return 0;
  if (p.P() > 0.8) return 0;
  if (Kp.P() > 0.5) return 0;
  if (Km.P() > 0.5) return 0;
  return 1;
}

