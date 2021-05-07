#include "Lcore.h"

int main(const int argc, const char * argv[]){
  
  if (argc < 2){
    return 0;
  }

  TString path = argv[1];
  TString filename = path + "distri.root";

  Long64_t Nsim;
  if (argc < 3) Nsim = 10000000;
  else Nsim = atoi(argv[2]);

  Initialize();

  TLorentzVector ki[2], kf[5];
  TLorentzVector kgold;
  double weight = 0;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, Mp * 197.0);

  TFile * fs = new TFile(filename.Data(), "RECREATE");

  TH1D * hQ2 = new TH1D("Q2_vertex", "", 1000, 1e-5, 2e-1);
  TH1D * hW = new TH1D("W_vertex", "", 100, 1.8, 2.2);
  TH2D * hWQ2 = new TH2D("WQ2_vertex", "", 100, 1.8, 2.2, 1000, 1e-5, 2e-1);

  TH1D * hQ2d = new TH1D("Q2_detected", "", 1000, 1e-5, 2e-1);
  TH1D * hWd = new TH1D("W_detected", "", 100, 1.8, 2.2);
  TH2D * hWQ2d = new TH2D("WQ2_detected", "", 100, 1.8, 2.2, 1000, 1e-5, 2e-1);
 
  TH1D * hQ2s = new TH1D("Q2_smeared", "", 1000, 1e-5, 2e-1);
  TH1D * hWs = new TH1D("W_smeared", "", 100, 1.8, 2.2);
  TH2D * hWQ2s = new TH2D("WQ2_smeared", "", 100, 1.8, 2.2, 1000, 1e-5, 2e-1);
  
  hQ2->SetDirectory(fs);
  hW->SetDirectory(fs);
  hWQ2->SetDirectory(fs);

  hQ2d->SetDirectory(fs);
  hWd->SetDirectory(fs);
  hWQ2d->SetDirectory(fs);

  hQ2s->SetDirectory(fs);
  hWs->SetDirectory(fs);
  hWQ2s->SetDirectory(fs);

  TLorentzVector PP, q;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
    weight = GENERATE::Event_eNKKN_BoundState(ki, kf);
    if (weight > 0){
      PP = kf[2] + kf[3] + kf[4];
      q = ki[0] - kf[0];
      hQ2->Fill(-q*q, weight);
      hW->Fill(PP.M(), weight);
      hWQ2->Fill(PP.M(), -q*q, weight);

      weight *= DETECTOR::Acceptance(kf[4], "p") * DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      hQ2d->Fill(-q*q, weight);
      hWd->Fill(PP.M(), weight);
      hWQ2d->Fill(PP.M(), -q*q, weight);

      DETECTOR::Smear(&kf[4], "p");
      DETECTOR::Smear(&kf[2], "K+");
      DETECTOR::Smear(&kf[3], "K-");
      PP = kf[2] + kf[3] + kf[4];
      hQ2s->Fill(-q*q, weight);
      hWs->Fill(PP.M(), weight);
      hWQ2s->Fill(PP.M(), -q*q, weight);
    }
  }

  hQ2->Scale(1.0/Nsim);
  hW->Scale(1.0/Nsim);
  hWQ2->Scale(1.0/Nsim);

  hQ2d->Scale(1.0/Nsim);
  hWd->Scale(1.0/Nsim);
  hWQ2d->Scale(1.0/Nsim);

  hQ2s->Scale(1.0/Nsim);
  hWs->Scale(1.0/Nsim);
  hWQ2s->Scale(1.0/Nsim);
  
  fs->Write();

  return 0;
}
