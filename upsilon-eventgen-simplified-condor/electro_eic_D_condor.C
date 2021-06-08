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
  double factor;
  double weight_gamma_p;
  double phase_space;
  int is_sub;
  int is_accept;
};

event myEvent;

int electro_eic_D_condor(Long64_t Nsim, int Ntotalbatch, int Nbatch, double Ebeam, double Hbeam, 
			 double ymin, double ymax, double Q2min, double Q2max){


  bool doFermiMotion = true;
  gRandom->SetSeed(0);

  double tmin = -2.5;
  double tmax = 0.0;
  // Electron beam energy and luminosity
  double lumi = pow(10,34); // cm^-2 s^-1
  lumi *= pow(10,-24); // cm^-2 = 10-24 barn^-1
  lumi *= pow(10,-9);  // barn^-1 = 10^-9 nb^-1
  lumi *= 0.5;         // Nuclear factor
  // lumi is now in units of nb^-1 s^-1

  double time = 3600.0 * 24 * 116; // 116 days

  // Set nuclear
  NUCLEAR::SetNuclear("D");
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
  // Set Upsilon1S production model
  UPSILONMODEL::SetModel("v2");
  
  // Set kinematic range
  double degtorad = M_PI / 180.0;
  GENERATE::Q2range[0] = Q2min;//GeV^2
  GENERATE::Q2range[1] = Q2max;//GeV^2
  GENERATE::yrange[0]  = ymin;
  GENERATE::yrange[1]  = ymax;
  GENERATE::trange[0]  = tmin;//GeV^2
  GENERATE::trange[1]  = tmax;//GeV^2

  // TTree Info
  TFile * fall2 = new TFile(Form("D-%lld-%d_%.0fx%.0f_Y_%.2f_%.2f_Q2_%.2f_%.2f/run%d.root",Nsim,Ntotalbatch,Ebeam,Hbeam,ymin,ymax,Q2min,Q2max,Nbatch),"RECREATE");
  TTree * tree = new TTree("tree","");
  tree->SetDirectory(fall2);

  TLorentzVector ki[2], kf[4], q;  
  double weight = 0.0;
  double MY1S = 9.46030;
  int is_sub;
  
  TLorentzVector *eIn = new TLorentzVector();
  TLorentzVector *pIn = new TLorentzVector();
  TLorentzVector *eOut = new TLorentzVector();
  TLorentzVector *pOut = new TLorentzVector();
  TLorentzVector *ePlusOut = new TLorentzVector();
  TLorentzVector *eMinusOut = new TLorentzVector();
  TLorentzVector vm;
  TLorentzVector gamma;
  TLorentzVector target_proton;
  // Add TTree Branches
  tree->Branch("event", &myEvent.x,
	       "x/D:y/D:Q2/D:W/D:W2/D:t/D:weight/D:factor/D:weight_gamma_p/D:phase_space/D:is_sub/I:is_accept/I");
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

  double BR = 2.38e-2;
  
  double beta = 0.0; // beta of the incoming proton
  for (Long64_t i = 0; i < Nsim/Ntotalbatch; i++){
    if (i % (Nsim/Ntotalbatch/10) == 0) cout << i/(Nsim/Ntotalbatch/10)*10 << "%" << endl;

    weight = 1.0;

    // Step 1.) Create lab beam electron and proton

    ki[0].SetXYZM(0, 0, -Ebeam, PARTICLE::e.M());
    ki[1].SetXYZM(0, 0, sqrt(Hbeam*Hbeam-PARTICLE::proton.M()*PARTICLE::proton.M()), PARTICLE::proton.M());
    target_proton=ki[1];

    // Step 2.) Get the beta of the incoming proton

    beta = ki[1].Beta();

    // Step 3.) Boost into the proton rest frame

    ki[1].Boost(0,0,-beta);

    // Step 4.) Change the momentum of this proton in that frame

    if(doFermiMotion)
      weight *= GENERATE::GetNucleon(&ki[1]);
    
    // Step 5.) Boost back into the lab frame

    ki[1].Boost(0,0,beta);

    // Step 6.) Run the event production

    weight *= GENERATE::Event_eN2eNee_Upsilon1S(ki, kf); 

    if (weight > 0.0){
      q = ki[0] - kf[0];
      gamma = q;
      gamma.Boost(0,0,-beta);

      // Step 6a.) Store event info

      myEvent.x = - (q*q) / (2.0 * ki[1] * q);
      myEvent.y = (ki[1] * q ) / (ki[1] * ki[0]);
      myEvent.Q2 = - (q * q);
      myEvent.W = sqrt( (ki[1] + q) * (ki[1] + q) );
      myEvent.W2 = (ki[1] + q) * (ki[1] + q);
      myEvent.t = (ki[1] - kf[1]) * (ki[1] - kf[1]);
      myEvent.phase_space = (Q2max-Q2min)*(ymax-ymin);
      eIn=(TLorentzVector*)ki[0].Clone();
      pIn=(TLorentzVector*)ki[1].Clone();
      eOut=(TLorentzVector*)kf[0].Clone();
      pOut=(TLorentzVector*)kf[1].Clone();
      ePlusOut=(TLorentzVector*)kf[2].Clone();
      eMinusOut=(TLorentzVector*)kf[3].Clone();
      vm=kf[2]+kf[3];
      // Step 6b.) Calculate if the process was "subthreshold"
      myEvent.is_sub = 0;
      if (gamma.E() < MY1S + (MY1S * MY1S - gamma * gamma) / (2.0 * Mp)){
	myEvent.is_sub = 1;
      }
      // Step 6c.) Calculate miscellaneous variables
      double gy = ki[1].M() * sqrt(myEvent.Q2) / (ki[1]*ki[0]);
      double epsilon = (1.0 - myEvent.y - 0.25 * gy * gy) / (1.0 - myEvent.y + 0.5 * myEvent.y * myEvent.y + 0.25 * gy * gy);
      double dipole = pow((MY1S*MY1S)/(myEvent.Q2+MY1S*MY1S),2.575);
      double R = pow((2.164*MY1S*MY1S + myEvent.Q2)/(2.164*MY1S*MY1S),2.131)-1.0;
      double alpha_em = 1.0 / 137.0;
      double gammaT = alpha_em/(2*M_PI)*(1.0+pow(1.0-myEvent.y,2))/(myEvent.y*myEvent.Q2);
      double factor = gammaT*(1.0+epsilon*R)*dipole;
      myEvent.factor = factor;
      // Step 6d.) Calculate event weights
      myEvent.weight = weight; // raw
      myEvent.weight_gamma_p = weight; // sigma (gamma + p --> V + p')
      // Step 6e.) Calculate acceptance of final state particles at EIC
      bool accept = (abs(eOut->Eta())<=5.0)&&(abs(ePlusOut->Eta())<=5.0)&&(abs(eMinusOut->Eta())<=5.0)&&(abs(pOut->Theta())>=0.002)&&(myEvent.y>0.01)&&(myEvent.y<0.8);
      myEvent.is_accept = accept;
      // Step 6f.) Fill tree
      tree->Fill();
    }
  }

  // Write the tree
  fall2->Write();
  fall2->Close();


  return 0;
}
