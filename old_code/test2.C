#include "Lcore.h"

int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 1000000;

  // Electron beam energy and luminosity
  double Ebeam = 12.0;//GeV
  double lumi = 1.3e31 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
  
  // Set nuclear
  NUCLEAR::SetNuclear("p");  
   
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 8.2;
  double kmax = 11.8;
  
  TH1D * hMJpsi1 = new TH1D("Mass_e+e-_Jpsi_1", "", 100, 2.6, 3.6);
  TH1D * hMJpsi2 = new TH1D("Mass_e+e-_Jpsi_2", "", 100, 2.6, 3.6);
  TH1D * hMJpsi3 = new TH1D("Mass_e+e-_Jpsi_3", "", 100, 2.6, 3.6);
  
  TLorentzVector ki[2], kf[4];
  double weight = 0.0;
  double acceptance = 1.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam);
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_gN2Nee_Jpsi(ki, kf); 

    if (weight > 0.0){     
      hMJpsi1->Fill( (kf[1]+kf[2]).M(), weight * acceptance);

      if (kf[1].Theta() / M_PI * 180.0 < 2.0) continue;
      if (kf[1].Theta() / M_PI * 180.0 > 120.0) continue;
      if (kf[1].P() < 0.4) continue;
      if (kf[2].Theta() / M_PI * 180.0 < 2.0) continue;
      if (kf[2].Theta() / M_PI * 180.0 > 120.0) continue;
      if (kf[2].P() < 0.4) continue;

      hMJpsi2->Fill( (kf[1]+kf[2]).M(), weight * acceptance);

      if (kf[0].Theta() / M_PI * 180.0 < 2.0) continue;
      if (kf[0].Theta() / M_PI * 180.0 > 120.0) continue;
      if (kf[0].P() < 0.4) continue;

      hMJpsi3->Fill( (kf[1]+kf[2]).M(), weight * acceptance);


      
    }
    
  }

  hMJpsi1->Scale(lumi*time/Nsim);
  hMJpsi2->Scale(lumi*time/Nsim);
  hMJpsi3->Scale(lumi*time/Nsim);

  cout << hMJpsi1->Integral() << endl;
  cout << hMJpsi2->Integral() << endl;
  cout << hMJpsi3->Integral() << endl;
  
  
  return 0;
}
