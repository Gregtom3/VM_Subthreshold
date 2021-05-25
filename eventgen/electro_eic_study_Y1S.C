#include "Lcore.h"
#include "TTree.h"

struct event{
  double x;
  double y;
  double Q2;
  double W;
  double W2;
  double t;
  double weight;
  double weight_lumi; 
  int is_sub;
};

event myEvent;

int electro_eic_study_Y1S(){

  gRandom->SetSeed(0);
  Long64_t Nsim = 1000000;
  double Ebeam = 10; // GeV
  double Hbeam = 100;

  // Electron beam energy and luminosity
  double lumi = pow(10,34); // cm^-2 s^-1
  lumi *= pow(10,-24); // cm^-2 = 10-24 barn^-1
  lumi *= pow(10,-9);  // barn^-1 = 10^-9 nb^-1
  lumi *= 1.0;         // Nuclear factor
  // lumi is now in units of nb^-1 s^-1
  //* pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Set nuclear
  NUCLEAR::SetNuclear("p");
  //  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  //  GENERATE::TF_fMomentum->SetNpx(1000);
  // Set Upsilon1S production model
  UPSILONMODEL::SetModel("v1");

  // Set Collider Mode
  COLLIDER::SetCollider("D");
  //  double beta = COLLIDER::GetBeta();
  
  // Set kinematic range
  double degtorad = M_PI / 180.0;
  GENERATE::Q2range[0] = 0.0;//GeV^2
  GENERATE::Q2range[1] = 100.0;//GeV^2
  GENERATE::Wrange[0]  = 12.0;//GeV
  GENERATE::Wrange[1]  = 14.0;//GeV
  GENERATE::trange[0]  = -1.0;//GeV^2
  GENERATE::trange[1]  = -0.0;//GeV^2

  // Set detector
  //  DETECTOR::SetDetector("SoLID");
  

  // TTree Info
  TFile * fall2 = new TFile(Form("result-eic-Y1S/eic_electro_%.1fGeV_%.1fGeV.root",Ebeam,Hbeam), "RECREATE");
  TTree * tree = new TTree("tree","");
  tree->SetDirectory(fall2);

  TLorentzVector ki[2], kf[4], q;
  
  double weight = 0.0;
  
  double MY1S = 9.46030;
  int count = 0;
  int is_sub;
  
  TLorentzVector *eIn = new TLorentzVector();
  TLorentzVector *pIn = new TLorentzVector();
  TLorentzVector *eOut = new TLorentzVector();
  TLorentzVector *pOut = new TLorentzVector();
  TLorentzVector *ePlusOut = new TLorentzVector();
  TLorentzVector *eMinusOut = new TLorentzVector();
  TLorentzVector vm;
  TLorentzVector gamma;
  
  // Add TTree Branches
  tree->Branch("event", &myEvent.x,
	       "x/D:y/D:Q2/D:W/D:W2/D:t/D:weight/D:weight_lumi/D:is_sub/I");
  tree->Branch("eIn","TLorentzVector",&eIn);
  tree->Branch("pIn","TLorentzVector",&pIn);
  tree->Branch("eOut","TLorentzVector",&eOut);
  tree->Branch("pOut","TLorentzVector",&pOut);
  tree->Branch("ePlusOut","TLorentzVector",&ePlusOut);
  tree->Branch("eMinusOut","TLorentzVector",&eMinusOut);
  tree->Branch("vm","TLorentzVector",&vm);
  tree->Branch("gamma","TLorentzVector",&q);
  tree->Branch("gamma_boost","TLorentzVector",&gamma);
  tree->Branch("Nsim",&Nsim,"Sims/I");
  
  double beta = 0.0; // beta of the incoming proton

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    // Create lab beam electron and proton
    ki[0].SetXYZM(0, 0, -Ebeam, PARTICLE::e.M());
    ki[1].SetXYZM(0, 0, sqrt(Hbeam*Hbeam-PARTICLE::proton.M()*PARTICLE::proton.M()), PARTICLE::proton.M());
    // Get the beta of the incoming proton
    beta = ki[1].Beta();
    // Boost into the proton rest frame
    ki[1].Boost(0,0,-beta);
    // Change the momentum of this proton in that frame
    //    weight = GENERATE::GetNucleon(&ki[1]);
    // Boost back into the lab frame
    ki[1].Boost(0,0,beta);
    // Run the event production
    weight = GENERATE::Event_eN2eNee_Upsilon1S(ki, kf); 

    if (weight > 0.0){
      count++;
      q = ki[0] - kf[0];
      gamma = q;
      gamma.Boost(0,0,-beta);
      // Store event info
      myEvent.x = - (q*q) / (2.0 * ki[1] * q);
      myEvent.y = (ki[1] * q ) / (ki[1] * ki[0]);
      myEvent.Q2 = - (q * q);
      myEvent.W = sqrt( (ki[1] + q) * (ki[1] + q) );
      myEvent.W2 = (ki[1] + q) * (ki[1] + q);
      myEvent.t = (ki[1] - kf[1]) * (ki[1] - kf[1]);
      
      eIn=(TLorentzVector*)ki[0].Clone();
      pIn=(TLorentzVector*)ki[1].Clone();
      eOut=(TLorentzVector*)kf[0].Clone();
      pOut=(TLorentzVector*)kf[1].Clone();
      ePlusOut=(TLorentzVector*)kf[2].Clone();
      eMinusOut=(TLorentzVector*)kf[3].Clone();
      // Calculate if the process was "subthreshold"
      myEvent.is_sub = 0;
      if (gamma.E() < MY1S + (MY1S * MY1S - gamma * gamma) / (2.0 * Mp)){
	myEvent.is_sub = 1;
      }

      myEvent.weight = weight;
      myEvent.weight_lumi = weight * lumi * time / Nsim; // Evts/Hour
      vm=kf[2]+kf[3];
      tree->Fill();
    }
  }
  cout << count << endl;
  GENERATE::print_fail();
 

  //2
  fall2->Write();
  fall2->Close();


  return 0;
}
