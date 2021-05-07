#include "Lcore.h"

int main(const int argc, const char * argv[]){

  //number of events to generate
  Long64_t Nsim = 100000000;

  //initialization
  Initialize();

  //incoming and outgoing particle 4-vectors
  TLorentzVector ki[2], kf[5];

  //weighting factor
  double weight = 0;

  //incoming electron
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);


  //generate event with a bound-state production
  //outgoing particles kf: e', N', [K+, K-, p]
  weight = GENERATE::Event_eNKKN_BoundState(ki, kf);

  //generate event with phi production
  //outgoing particles kf: e', p, [K+, K-]
  weight = GENERATE::Event_eNKK_Phi(ki, kf);

  //generate event with direct KK production
  //outgoing particles kf: e', p, K+, K-
  weight = GENERATE::Event_eNKK_KK(ki, kf);

  //generate event with Lambda(1520) production
  //outgoing particles kf: e', p, K+, K-
  weight = GENERATE::Event_eNKK_L1520(ki, kf);


  return 0;
}

  
