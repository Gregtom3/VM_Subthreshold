#include "Lcore.h"

int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 10000000;

  // Electron beam energy and luminosity
  double Ebeam = 11.0;//GeV
  double lumi = 1.2e37 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
  
  // Set nuclear
  NUCLEAR::SetNuclear("p");  
   
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  PHIMODEL::SetModel("fit");

  // Set scattered electron range
  
  TH1D * hMJpsi = new TH1D("Mass_e+e-_Jpsi", "", 100, 0.0, 2.0);
  
  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZT(0, 0, Ebeam, Ebeam);
  double weight = 0.0;
  double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_gN2Nee_Phi(ki, kf); 

    if (weight > 0.0){
      acceptance = 1.0;
      
      hMJpsi->Fill( (kf[1]+kf[2]).M(), weight * acceptance);
      //cout << weight << "  " << (kf[1] + kf[2]).M() << endl;
    }
    
  }

  hMJpsi->Scale(lumi*time/Nsim);

  cout << hMJpsi->Integral() << endl;

  
  return 0;
}
