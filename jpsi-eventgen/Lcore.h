#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include "TROOT.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TFile.h"
//#include "TTree.h"
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/SpecFuncMathCore.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "TLegend.h"

#include "Lparticle.h"

using namespace std;
using namespace ROOT;
const double Mp = PARTICLE::proton.M();

namespace COLLIDER{

  double A = 1.0;
  double ml = 0.000511; // GeV^2
  double mn = 0.938272013; // GeV^2
  
  double gamma = 1.0;
  double beta  = 1.0;

  double GetBeta(){
    return beta;
  }
  
  int SetCollider(const char * Astring = "p",
		  const double Enucleon = 1.0){

    if (strcmp(Astring, "D") == 0){
      A = 2.0;
    }

    if (strcmp(Astring, "Au") == 0){
      A = 197.0;
    }
    
    mn = (0.938272013 * A); // Total Mass of Incoming Nucleus

    if(Enucleon<0.938272013){ // User-entered Enucleon < Physically possible
      beta = 0.;
      gamma = 1.;
    }

    else{
      double En_total = A * Enucleon; // Total Energy of Incoming Nucleus
      double Pn_total = sqrt(En_total*En_total - mn*mn); // Total Momentum of Incoming Nucleus
      beta = -Pn_total/sqrt(mn*mn+Pn_total*Pn_total);
      gamma = 1./sqrt(1-beta*beta);
    }
    return 0;
  }
}
namespace NUCLEAR{

  int flag = 0;

  double A = 1.0;
  double Z = 1.0;

  double (* fMomentum)(const double * p0, const double * par);
  double (* fEnergy)(const double * E0, const double * par);

  double Momentum_D(const double * p0, const double * par){//non-normalized
    double p = p0[0];
    double a = 0.0456;
    double b = 0.2719;
    double result = pow(1.0 / (p * p + a * a) - 1.0 / (p * p + b * b), 2);//C. Weiss 2014
    return p * p * result;
  }

  double Momentum_C12(const double * p0, const double * par){//non-normalized
    double p = p0[0];//nucleon momentum in C12 in unit of GeV
    if (p < 0.0){
      std::cerr << "Unphysical momentum value in CarbonMomentum!" << std::endl;
      return -1.0;
    }
    double A0 = 10.6552;
    double A2 = 32.3105;
    double B2 = 8.43226;
    double result = (A0 + pow(A2 * p, 2)) * exp(-pow(B2 * p, 2));
    return p * p * result / 0.0241519;//Normalized momentum distribution
  }
  
  double Energy_C12(const double * E0, const double * par){
    double E = E0[0];//nucleon missing energy in C12 in unit of GeV
    if (E <= 0.0){
      std::cerr << "Unphysical energy value in CarbonEnergy!" << std::endl;
      return -1.0;
    }
    double A1 = 0.569161;
    double a1 = 0.0175;
    double b1 = 0.00510586;
    double A2 = 0.0683607;
    double a2 = 0.0369238;
    double b2 = 0.0173035;
    double result = A1 * exp(-pow((E - a1) / b1, 2)) + A2 * exp(-pow((E - a2) / b2, 2));
    return result / 0.00724478;//Normalized missing energy distribution
  }

  double Momentum_Au197(const double * p0, const double * par = 0){//non-normalized
    double p = p0[0];//nucleon momentum in Au197 in unit of GeV
    if (p < 0.0){
      std::cerr << "Unphysical momentum value in GoldMomentum!" << std::endl;
      return -1.0;
    }
    double A0 = 58.3382;
    double A2 = 69.2938;
    double B2 = 7.82756;
    double result = (A0 + pow(A2 * p, 2)) * exp(-pow(B2 * p, 2));
    return p * p * result / 0.162508;//Normalized momentum distribution
  }
  
  double Energy_Au197(const double * E0, const double * par = 0){
    double E = E0[0];//nucleon missing energy in Au197 in unit of GeV
    if (E <= 0.0){
      std::cerr << "Unphysical energy value in GoldEnergy!" << std::endl;
      return -1.0;
    }
    double A1 = 1.73622;
    double a1 = 3.07375;
    double b1 = 0.645561;
    double A2 = 14.1433;
    double a2 = 0.795058;
    double result = A1 * atan(A2 * pow(E/0.01, a1)) * exp(-b1 * pow(E/0.01, a2));
    return result / 0.0433967;//Normalized missing energy distribution
  }

  int SetNuclear(const char * nuclear = "p"){
    flag = 1;
    if (strcmp(nuclear, "D") == 0){
      flag = 2;
      fMomentum = &Momentum_D;
    }
    else if (strcmp(nuclear, "C12") == 0){
      fMomentum = &Momentum_C12;
      fEnergy = &Energy_C12;
    }
    else if (strcmp(nuclear, "Au197") == 0){
      fMomentum = &Momentum_Au197;
      fEnergy = &Energy_Au197;
    }
    else {
      cout << "No matching nuclear! Set to proton!" << endl;
      flag = 0;
    }
    return 0;
  }

  
}

namespace PSI2SMODEL{//Model of Psi2S production (based on 2+3g model with a (x0.16) Psi2S / J/Psi  production ratio

  double (*dSigmaPsi2S)(const double, const double);

  double dSigmaPsi2S_2g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 1569.43; // <---- Needs to be changed
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.686097;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / pow(R * M, 2);//nb GeV^-2
    return result;//ds/dt in unit GeV^-4
  }
  
  double dSigmaPsi2S_23g(const double x, const double t){//Brodsky et al. PLB498 
    double N2g = 1569.43;
    double N3g = 698.868;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.686097;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / (R * R * M * M) + N3g * v * exp(1.13 * t) / pow(R * M, 4);//nb GeV^-2
    return result;//ds/dt in unit GeV^-4
  }

  int SetModel(const char * model = "2g"){
    if (strcmp(model, "2g") == 0)
      dSigmaPsi2S = &dSigmaPsi2S_2g;
    else if (strcmp(model, "23g") == 0)
      dSigmaPsi2S = &dSigmaPsi2S_23g;
    else {
      cout << "No matching model! Set to 2g model!" << endl;
      dSigmaPsi2S = dSigmaPsi2S_2g;
    }
    return 0;
  }
}

namespace JPSIMODEL{//Model of J/psi production

  double (*dSigmaJpsi)(const double, const double);
  
  double dSigmaJpsi_2g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 7.5671e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / pow(R * M, 2);//nb GeV^-2
    return result;
  }

  double dSigmaJpsi_23g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 6.499e3;
    double N3g = 2.894e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / (R * R * M * M) + N3g * v * exp(1.13 * t) / pow(R * M, 4);//nb GeV^-2
    return result;
  }


  int SetModel(const char * model = "2g"){
    if (strcmp(model, "2g") == 0)
      dSigmaJpsi = &dSigmaJpsi_2g;
    else if (strcmp(model, "23g") == 0)
      dSigmaJpsi = &dSigmaJpsi_23g;
    else {
      cout << "No matching model! Set to 2g model!" << endl;
      dSigmaJpsi = dSigmaJpsi_2g;
    }
    return 0;
  }
}
namespace JPSIPomLQCD{//Pomeron-LQCD model of J/psi production from the proton

  double cth[142];
  double ds[21][142];
  double W[21];
  ROOT::Math::Interpolator * dsigma[21];

  double dSigmaJpsi(const double W0, const double cth){
    int idx = (int) floor((W0 - 4.036) / 0.05);
    if (idx < 0) return 0;
    if (idx < 21){
      double ds1 = dsigma[idx]->Eval(cth);
      double ds2 = dsigma[idx+1]->Eval(cth);
      double res = ds1 + (ds2 - ds1) * (W0 - W[idx]) / 0.05;
      return res * 1.0e-7 / pow(0.197327, 2);//in unit GeV^-2
    }
    if (idx >= 21){
      return dsigma[21]->Eval(cth) * 1.0e-7 / pow(0.197327, 2);
    }
    return 0;
  }       

  int SetModel(){
    ifstream infile("harrymodel/pomeron-lqcd/diff.dat");
    char tmp[100];
    infile.getline(tmp, 100);
    double dum, th;
    cth[0] = -1.0;
    cth[141] = 1.0;
    for (int i = 0; i < 21; i++){
      for (int j = 140; j >= 1; j--){
	infile >> dum >> th >> ds[i][j];
	cth[j] = cos(th * M_PI / 180.0);
      }
      W[i] = dum;
      ds[i][0] = ds[i][1];
      ds[i][141] = ds[i][140];
    }
    infile.close();
    for (int i = 0; i < 21; i++){
      dsigma[i] = new ROOT::Math::Interpolator(142, ROOT::Math::Interpolation::kCSPLINE);
      dsigma[i]->SetData(142, cth, ds[i]);
    }
    return 0;
  }

}

namespace GENERATE{

  TRandom3 random(0);
  TGenPhaseSpace GenPhase;
  double Weight = 0.0;

  TF1 * TF_fBremsstrahlung;
  TF1 * TF_fMomentum;
  TF1 * TF_fEnergy;

  bool do_PomLQCD = false;
  /* Bremsstrahlung photon */

  int fail = 0;

  // Weights
  double WEIGHT_JACOBIAN = 0.0;
  double WEIGHT_NUCLEON_FLUX = 0.0;
  double WEIGHT_PHOTON_FLUX = 0.0;
  double WEIGHT_DODT = 0.0;
  double WEIGHT_DECAY = 0.0;
  
  double Bremsstrahlung(const double * y, const double * par){//ds/dy approximate expression
    //E0: electron beam energy; k: photon energy
    if (y[0] < 0.01) {// Infrared cut
      std::cerr << "Out of range in Bremsstrahlung!" << std::endl;
      return -1.0;
    }
    double result = (4.0 / 3.0 - 4.0 / 3.0 * y[0] + y[0] * y[0]) / y[0];
    return result;
  }
  
  double BremsstrahlungPhoton(TLorentzVector * q, const double kmin, const double kmax, const double E){//Generate a Bremsstrahlung photon ! Assume d / X0 = 0.01 radiator!
    //q: photon; E: electron beam energy; [kmin, kmax]: photon energy range
    double ymin = kmin / E;
    double ymax = kmax / E;
    double y = TF_fBremsstrahlung->GetRandom(ymin, ymax);
    q[0].SetXYZT(0.0, 0.0, y * E, y * E);
    WEIGHT_PHOTON_FLUX = 0.01 * (4.0 / 3.0 * log(ymax / ymin) - 4.0 / 3.0 * (ymax - ymin) + 1.0 / 2.0 * (ymax * ymax - ymin * ymin));
    return WEIGHT_PHOTON_FLUX;
  }

  int SetBremsstrahlung(){
    TF_fBremsstrahlung = new TF1("fBremsstrahlung", Bremsstrahlung, 0.01, 1.0, 0);
    TF_fBremsstrahlung->SetNpx(1000);
    return 0;
  }

  int print_fail(){ //debugger
    std::cout << "The specific weight calculation has failed " << fail << " times." << std::endl;
    return 0;
  }

  void SetPomLQCD(){
    do_PomLQCD = true;
  }
  /* Nucleon from a nuclear target */
  
  double GetNucleon(TLorentzVector * P){
    if (NUCLEAR::flag > 0){
      double p = TF_fMomentum->GetRandom();
      double cth = random.Uniform(-1.0, 1.0);
      double phi = random.Uniform(-M_PI, M_PI);
      double sth = sqrt(1.0 - cth * cth);
      if (NUCLEAR::flag == 1){
	double dE = TF_fEnergy->GetRandom();
	P->SetXYZT(p * sth * cos(phi), p * sth * sin(phi), p * cth, sqrt(p * p + Mp * Mp) - dE);
      }
      else {//flag = 2, deuteron
	double E = 2.0 * Mp - sqrt(Mp * Mp + p * p);
	P->SetXYZT(p * sth * cos(phi), p * sth * sin(phi), p * cth, E);
      }
    }
    else //flag = 0, proton
      P->SetXYZT(0.0, 0.0, 0.0, Mp);
    return 1.0;
  }

  double cthrange[2] = {-1.0, 1.0};
  double perange[2] = {0.0, 10.0};


  /* Added 5/19/2021 for Upsilon Production Model based on Slyvester */
  double Q2range[2] = {0.0, 10.0};
  double Wrange[2] = {0.0, 9999.0};
  double trange[2] = {-100.0, 0.0};
  double yrange[2] = {0.1,0.8};

  double VirtualPhoton(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', gamma
    double m = PARTICLE::e.M();
    double Pe = random.Uniform(perange[0], perange[1]);
    double cth = random.Uniform(cthrange[0], cthrange[1]);
    double sth = sqrt(1.0 - cth * cth);
    double phi = random.Uniform(-M_PI, M_PI);
    kf[0].SetXYZM(Pe * sth * cos(phi), Pe * sth * sin(phi), Pe * cth, m);//e'
    kf[1] = ki[0] - kf[0];//virtual photon
    double W2 = (kf[1] + ki[1]) * (kf[1] + ki[1]);//W^2 = (P + q)^2
    if (W2 < Mp * Mp)
      {
	return 0;//below the lowest state
      }
    double Q2 = - kf[1] * kf[1];//Q^2 = -q^2
    double alpha_em = 1.0 / 137.0;
    double couple = 4.0 * M_PI * alpha_em;
    double flux = sqrt(pow(ki[1] * kf[1], 2) + Q2 * Mp * Mp) / sqrt(pow(ki[0] * ki[1], 2) - m * m * Mp * Mp);
    double amp = (2.0 * Q2 - 4.0 * m * m) / (Q2 * Q2);
    double phase = kf[0].P() * kf[0].P() / (2.0 * kf[0].E() * pow(2.0 * M_PI, 3));
    double volume = 2.0 * M_PI * abs(perange[1] - perange[0]) * abs(cthrange[1] - cthrange[0]);
    double y = (ki[1] * kf[1]) / (ki[1] * ki[0]);
    double gy = ki[1].M() * sqrt(Q2) / (ki[1] * ki[0]);
    double epsilon = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    WEIGHT_PHOTON_FLUX = couple * flux * amp * phase * volume / (1.0 - epsilon);
    return WEIGHT_PHOTON_FLUX;
  }

  /* Psi2S productions */
   double Psi2SPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: Psi2S, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    double W = Pout.M();
    double Mpsi2S = PARTICLE::psi2S.RandomM();
    if (W < Mpsi2S + Mp)
      {
	fail++;
	return 0;//below the threshold
      }
    double mass[2] = {Mpsi2S, Mp};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();//uniform generate in 4pi solid angle
    kf[0] = *GenPhase.GetDecay(0);//Psi2S
    kf[1] = *GenPhase.GetDecay(1);//N'    
    double t = (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double x = (2.0 * Mp * Mpsi2S + Mpsi2S * Mpsi2S) / (W * W - Mp * Mp);
    double k = sqrt(pow(W * W - kf[0] * kf[0] - kf[1] * kf[1], 2) - 4.0 * (kf[0] * kf[0]) * (kf[1] * kf[1])) / (2.0 * W);//final state c.m. momentum
    double q = sqrt(pow(W * W - ki[0] * ki[0] - ki[1] * ki[1], 2) - 4.0 * (ki[0] * ki[0]) * (ki[1] * ki[1])) / (2.0 * W);//initial state c.m. momentum
    double volume = 4.0 * M_PI;
    double Jac = 2.0 * k * q / (2.0 * M_PI);
    double flux = sqrt(pow(ki[0] * ki[1], 2) - (ki[0] * ki[0]) * (ki[1] * ki[1])) / (Mp * ki[0].P());
    double cth = (sqrt(ki[0] * ki[0] + q * q) * sqrt(kf[0] * kf[0] + k * k) - ki[0] * kf[0]) / (q * k);
    WEIGHT_DODT =  PSI2SMODEL::dSigmaPsi2S(x,t);
    WEIGHT_JACOBIAN = Jac * volume;
    WEIGHT_NUCLEON_FLUX = flux;
    double weight = WEIGHT_DODT * WEIGHT_JACOBIAN * WEIGHT_NUCLEON_FLUX;
    return weight;//GeV^-2
   }
   

   double Psi2SElectroproduction(TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', Psi2S, N'
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector ki2[2] = {kf[1], ki[1]};//Initial state: virtual photon N
    double weight2 = Psi2SPhotoproduction(ki2, &kf[1]);//Generate Psi2S N' from virtual photon production
    return weight1 * weight2;
  }
  
  double Event_gN2Nee_Psi2S(TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: N', [e+, e-]
    TLorentzVector kf1[2];//Psi2S, N'
    double weight = Psi2SPhotoproduction(ki, kf1);
    kf[0] = kf1[1];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[0], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double Mpsi = kf1[0].M();
    double Ep = kf1[0] * kf1[1] / Mpsi;//recoil proton energy in Psi2S rest frame
    double p = sqrt(Ep * Ep - Mp * Mp);//recoil proton momentum in Psi2S rest frame
    double l = sqrt(pow(Mpsi * Mpsi - kf[1] * kf[1] - kf[2] * kf[2], 2) - 4.0 * (kf[1] * kf[1]) * (kf[2] * kf[2])) / (2.0 * Mpsi);//decayed lepton momentum in Psi2S rest frame
    double cth = (Ep * Mpsi / 2.0 - kf[1] * kf[0]) / (p * l);//cos(theta) between final lepton and final proton in Psi2S rest frame
    double r = 0.0;
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 7.93e-3;//Branch ratio to e+e- (Psi_2S)
    WEIGHT_DECAY = wth * branch;
    return weight * WEIGHT_DECAY;
  }
  
  double Event_eN2eNee_Psi2S(TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', N', [e+, e-]
    TLorentzVector kf1[3];//e', Psi2S, N'
    double weight = Psi2SElectroproduction(ki, kf1);
    kf[0] = kf1[0];//e'
    kf[1] = kf1[2];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//e+
    kf[3] = *GenPhase.GetDecay(1);//e-
    double Mpsi = kf1[1].M();
    double Ep = kf1[2] * kf1[1] / Mpsi;//recoil proton energy in Psi2S rest frame
    double p = sqrt(Ep * Ep - Mp * Mp);//recoil proton momentum in Psi2S rest frame
    double l = sqrt(pow(Mpsi * Mpsi - kf[2] * kf[2] - kf[3] * kf[3], 2) - 4.0 * (kf[2] * kf[2]) * (kf[3] * kf[3])) / (2.0 * Mpsi);//decayed lepton momentum in Psi2S rest frame
    double cth = (Ep * Mpsi / 2.0 - kf[2] * kf[1]) / (p * l);//cos(theta) between final lepton and final proton in Psi2S rest frame
    double y = (ki[0].E() - kf[0].E()) / ki[0].E();
    double Q2 = - (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double gy = sqrt(Q2) / ki[0].E();
    double epsilon = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    double R = pow(1.0 + Q2 / 2.164 / pow(Mpsi,2), 2.131) - 1.0;
    double r = epsilon * R / (1.0 + epsilon * R);
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 7.93e-3;//Branch ratio to e+e- (Psi_2S)
    WEIGHT_DECAY = wth * branch;
    return weight * WEIGHT_DECAY;
  }


  /* J/psi productions */
  
  /* modified 2020-09-09 to Harry's pomeron-LQCD model, arXiv:2004.13934 */
  double JpsiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: Jpsi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    double W = Pout.M();
    double MJpsi = PARTICLE::Jpsi.RandomM();
    if (W < MJpsi + Mp)
      {
	fail++;
	return 0;//below the threshold
      }
    double mass[2] = {MJpsi, Mp};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();//uniform generate in 4pi solid angle
    kf[0] = *GenPhase.GetDecay(0);//Jpsi
    kf[1] = *GenPhase.GetDecay(1);//N'
    double t = (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double x = (2.0 * Mp * MJpsi + MJpsi * MJpsi) / (W * W - Mp * Mp);
    double k = sqrt(pow(W * W - kf[0] * kf[0] - kf[1] * kf[1], 2) - 4.0 * (kf[0] * kf[0]) * (kf[1] * kf[1])) / (2.0 * W);//final state c.m. momentum
    double q = sqrt(pow(W * W - ki[0] * ki[0] - ki[1] * ki[1], 2) - 4.0 * (ki[0] * ki[0]) * (ki[1] * ki[1])) / (2.0 * W);//initial state c.m. momentum
    double volume = 4.0 * M_PI;
    double Jac = 2.0 * k * q / (2.0 * M_PI);
    double flux = sqrt(pow(ki[0] * ki[1], 2) - (ki[0] * ki[0]) * (ki[1] * ki[1])) / (Mp * ki[0].P());
    WEIGHT_DODT = JPSIMODEL::dSigmaJpsi(x,t);
    WEIGHT_JACOBIAN = Jac * volume;
    WEIGHT_NUCLEON_FLUX = flux;
    double weight = WEIGHT_DODT * WEIGHT_JACOBIAN * WEIGHT_NUCLEON_FLUX;
    return weight;//GeV^-2
  }

  double JpsiElectroproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', Jpsi, N'
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector ki2[2] = {kf[1], ki[1]};//Initial state: virtual photon N
    double weight2 = JpsiPhotoproduction(ki2, &kf[1]);//Generate Jpsi N' from virtual photon production
    return weight1 * weight2;
  }

  double Event_gN2Nee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: N', [e+, e-]
    TLorentzVector kf1[2];//Jpsi, N'
    double weight = JpsiPhotoproduction(ki, kf1);
    kf[0] = kf1[1];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[0], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double Mj = kf1[0].M();
    double Ep = kf1[0] * kf1[1] / Mj;//recoil proton energy in Jpsi rest frame
    double p = sqrt(Ep * Ep - Mp * Mp);//recoil proton momentum in Jpsi rest frame
    double l = sqrt(pow(Mj * Mj - kf[1] * kf[1] - kf[2] * kf[2], 2) - 4.0 * (kf[1] * kf[1]) * (kf[2] * kf[2])) / (2.0 * Mj);//decayed lepton momentum in Jpsi rest frame
    double cth = (Ep * Mj / 2.0 - kf[1] * kf[0]) / (p * l);//cos(theta) between final lepton and final proton in Jpsi rest frame
    double r = 0.0;
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 5.971e-2;//Branch ratio to e+e- (J/Psi)
    WEIGHT_DECAY = wth * branch;
    return weight * WEIGHT_DECAY;
  }
  
  double Event_eN2eNee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', N', [e+, e-]
    TLorentzVector kf1[3];//e', Jpsi, N'
    double weight = JpsiElectroproduction(ki, kf1);
    kf[0] = kf1[0];//e'
    kf[1] = kf1[2];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//e+
    kf[3] = *GenPhase.GetDecay(1);//e-
    double Mj = kf1[1].M();
    double Ep = kf1[2] * kf1[1] / Mj;//recoil proton energy in Jpsi rest frame
    double p = sqrt(Ep * Ep - Mp * Mp);//recoil proton momentum in Jpsi rest frame
    double l = sqrt(pow(Mj * Mj - kf[2] * kf[2] - kf[3] * kf[3], 2) - 4.0 * (kf[2] * kf[2]) * (kf[3] * kf[3])) / (2.0 * Mj);//decayed lepton momentum in Jpsi rest frame
    double cth = (Ep * Mj / 2.0 - kf[2] * kf[1]) / (p * l);//cos(theta) between final lepton and final proton in Jpsi rest frame
    double y = (ki[0].E() - kf[0].E()) / ki[0].E();
    double Q2 = - (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double gy = sqrt(Q2) / ki[0].E();
    double epsilon = (1.0 - y - 0.24 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    double R = pow(1.0 + Q2 / 2.164 / pow(Mj,2), 2.131) - 1.0;
    double r = epsilon * R / (1.0 + epsilon * R);
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 5.971e-2;//Branch ratio to e+e- (J/Psi)
    WEIGHT_DECAY = wth * branch;
    return weight * wth * branch;
  }
}


namespace DETECTOR{

  TRandom3 random(0);

  TFile * facc_clas, * facc_solid, * facc_alert;
  TFile * fres_clas, * fres_electron_solid, * fres_positron_solid, * fres_proton_solid, * fres_alert1, * fres_alert2;
  TH3F * acc_ele_clas, * acc_pos_clas, * acc_pip_clas, * acc_pim_clas, * acc_Kp_clas, * acc_Km_clas, * acc_proton_clas;
  TH2D * acc_ele_solid, * acc_pos_solid, * acc_proton_solid;
  TH2D * acc_proton_alert, * acc_Kp_alert, * acc_pip_alert, * acc_Km_alert, * acc_pim_alert;
  TH2D * res_Kp_alert_p, * res_Kp_alert_theta, * res_Kp_alert_phi;
  TH2D * res_Km_alert_p, * res_Km_alert_theta, * res_Km_alert_phi;
  TH2D * res_pip_alert_p, * res_pip_alert_theta, * res_pip_alert_phi;
  TH2D * res_pim_alert_p, * res_pim_alert_theta, * res_pim_alert_phi;
  TH2D * res_proton_alert_p, * res_proton_alert_theta, * res_proton_alert_phi;
  TH2D * res_electron_solid_p, * res_electron_solid_theta, * res_electron_solid_phi;
  TH2D * res_positron_solid_p, * res_positron_solid_theta, * res_positron_solid_phi;
  TH2D * res_proton_solid_p, * res_proton_solid_theta, * res_proton_solid_phi;
  

  int SetDetector(const char * detector = 0){
    if (strcmp(detector, "CLAS12") == 0){
      facc_clas = new TFile("acceptance/clasev_acceptance.root", "r");
      acc_pip_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pip");
      acc_pim_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pim");
      acc_ele_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_ele");
      acc_pos_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pip");//!!!!
      acc_proton_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pip");//!!!!
      acc_Kp_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pip");//!!!!
      acc_Km_clas = (TH3F *) facc_clas->Get("acceptance_PThetaPhi_pim");//!!!!
    }
    else if (strcmp(detector, "SoLID") == 0){
      facc_solid = new TFile("acceptance/acceptance_solid_JPsi_electron_target315_1e7_output_final.root", "r");
      acc_ele_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
      acc_pos_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
      acc_proton_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
      fres_electron_solid = new TFile("acceptance/JPsi_electron_resolution_2d.root", "r");
      fres_positron_solid = new TFile("acceptance/JPsi_electron_resolution_2d.root", "r");//!!!!
      fres_proton_solid = new TFile("acceptance/JPsi_proton_resolution_2d.root", "r");
      res_electron_solid_p = (TH2D *) fres_electron_solid->Get("p_resolution");
      res_electron_solid_theta = (TH2D *) fres_electron_solid->Get("theta_resolution");
      res_electron_solid_phi = (TH2D *) fres_electron_solid->Get("phi_resolution");
      res_positron_solid_p = (TH2D *) fres_positron_solid->Get("p_resolution");
      res_positron_solid_theta = (TH2D *) fres_positron_solid->Get("theta_resolution");
      res_positron_solid_phi = (TH2D *) fres_positron_solid->Get("phi_resolution");
      res_proton_solid_p = (TH2D *) fres_proton_solid->Get("p_resolution");
      res_proton_solid_theta = (TH2D *) fres_proton_solid->Get("theta_resolution");
      res_proton_solid_phi = (TH2D *) fres_proton_solid->Get("phi_resolution");
      res_electron_solid_p->Scale(1.5);
      res_electron_solid_theta->Scale(1.5);
      res_electron_solid_phi->Scale(1.5);
      res_positron_solid_p->Scale(1.5);
      res_positron_solid_theta->Scale(1.5);
      res_positron_solid_phi->Scale(1.5);
      res_proton_solid_p->Scale(1.5);
      res_proton_solid_theta->Scale(1.5);
      res_proton_solid_phi->Scale(1.5);
    }
    else if (strcmp(detector, "ALERT") == 0){
      facc_alert = new TFile("acceptance/acc_alert_20190427.root", "r");
      fres_alert1 = new TFile("acceptance/res_kp_20190429.root", "r");
      fres_alert2 = new TFile("acceptance/res_proton_20190429.root", "r");    
      acc_proton_alert = (TH2D *) facc_alert->Get("h0");
      acc_Kp_alert = (TH2D *) facc_alert->Get("h1");
      acc_pip_alert = (TH2D *) facc_alert->Get("h2");
      acc_Km_alert = (TH2D *) facc_alert->Get("h3");
      acc_pim_alert = (TH2D *) facc_alert->Get("h4");
      res_Kp_alert_p = (TH2D *) fres_alert1->Get("h1");
      res_Kp_alert_theta = (TH2D *) fres_alert1->Get("h2");
      res_Kp_alert_phi = (TH2D *) fres_alert1->Get("h3");
      res_Km_alert_p = (TH2D *) fres_alert1->Get("h1");
      res_Km_alert_theta = (TH2D *) fres_alert1->Get("h2");
      res_Km_alert_phi = (TH2D *) fres_alert1->Get("h3");
      res_pip_alert_p = (TH2D *) fres_alert1->Get("h1");
      res_pip_alert_theta = (TH2D *) fres_alert1->Get("h2");
      res_pip_alert_phi = (TH2D *) fres_alert1->Get("h3");
      res_pim_alert_p = (TH2D *) fres_alert1->Get("h1");
      res_pim_alert_theta = (TH2D *) fres_alert1->Get("h2");
      res_pim_alert_phi = (TH2D *) fres_alert1->Get("h3");
      res_proton_alert_p = (TH2D *) fres_alert2->Get("h1");
      res_proton_alert_theta = (TH2D *) fres_alert2->Get("h2");
      res_proton_alert_phi = (TH2D *) fres_alert2->Get("h3");
    } 
    return 0;
  }
  
  double AcceptanceCLAS12FD(const TLorentzVector P, const char * part){
    double p = P.P();
    double theta = P.Theta() * 180.0 / M_PI;
    if (theta > 35.0) return 0;//only forward detector
    double phi = P.Phi() * 180.0 / M_PI;
    if (phi < 0) phi = phi + 360.0;
    TH3F * acc;
    if (strcmp(part, "e") == 0 || strcmp(part, "e-") == 0) acc = acc_ele_clas;
    else if (strcmp(part, "e+") == 0) acc = acc_pos_clas;
    else if (strcmp(part, "p") == 0) acc = acc_proton_clas;
    else if (strcmp(part, "K+") == 0) acc = acc_Kp_clas;
    else if (strcmp(part, "K-") == 0) acc = acc_Km_clas;
    else if (strcmp(part, "pi+") == 0) acc = acc_pip_clas;
    else if (strcmp(part, "pi-") == 0) acc = acc_pim_clas;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(phi);
    int biny = acc->GetYaxis()->FindBin(theta);
    int binz = acc->GetZaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny, binz);
    if (strcmp(part, "K+") == 0 || strcmp(part, "K-") == 0) result *= exp(-6.5 / Phys::c / PARTICLE::K.Tau() / P.Beta() / P.Gamma());//kaon decay
    return result;
  }

  double AcceptanceCLAS12CD(const TLorentzVector P, const char * part){
    double p = P.P();
    double theta = P.Theta() * 180.0 / M_PI;
    if (theta < 35.0 || theta > 125.0) return 0;
    if (p > 0.3) return 1;
    return 0;
  } 

  double AcceptanceSoLID(const TLorentzVector P, const char * part){
    double p = P.P();
    double theta = P.Theta() * 180.0 / M_PI;
    TH2D * acc;
    if (strcmp(part, "e") == 0 || strcmp(part, "e-") == 0) acc = acc_ele_solid;
    else if (strcmp(part, "e+") == 0) acc = acc_pos_solid;
    else if (strcmp(part, "p") == 0) {
      if (theta > 8.0 && theta < 15.0 && p < 4.5)
	acc = acc_proton_solid;
      else if (theta > 15.0 && theta < 26.0 && p < 2.0)
	acc = acc_proton_solid;
      else return 0;
    }
    else return 0;
    int binx = acc->GetXaxis()->FindBin(theta);
    int biny = acc->GetYaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny);
    return result;
  }

  double SmearSoLID(TLorentzVector &P, const char * part){
    double acc = AcceptanceSoLID(P, part);
    if (acc == 0) return acc;
    double dp, dtheta, dphi;
    double p = P.P();
    double theta = P.Theta();
    double phi = P.Phi();
    double m = P.M();
    if (strcmp(part, "e") == 0 || strcmp(part, "e-") == 0){
      dp = res_electron_solid_p->GetBinContent(res_electron_solid_p->FindBin(p, theta * 180.0 / M_PI)) / 100.0;//dp / p
      dtheta = res_electron_solid_theta->GetBinContent(res_electron_solid_theta->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
      dphi = res_electron_solid_phi->GetBinContent(res_electron_solid_phi->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
    }
    else if (strcmp(part, "e+") == 0){
      dp = res_positron_solid_p->GetBinContent(res_positron_solid_p->FindBin(p, theta * 180.0 / M_PI)) / 100.0;//dp / p
      dtheta = res_positron_solid_theta->GetBinContent(res_positron_solid_theta->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
      dphi = res_positron_solid_phi->GetBinContent(res_positron_solid_phi->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
    }
    else if (strcmp(part, "p") == 0){
      dp = res_proton_solid_p->GetBinContent(res_proton_solid_p->FindBin(p, theta * 180.0 / M_PI)) / 100.0;//dp / p
      dtheta = res_proton_solid_theta->GetBinContent(res_proton_solid_theta->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
      dphi = res_proton_solid_phi->GetBinContent(res_proton_solid_phi->FindBin(p, theta * 180.0 / M_PI)) / 1000.0;
    }
    else 
      return acc;
    p = p * random.Gaus(1.0, dp);
    theta = random.Gaus(theta, dtheta);
    phi = random.Gaus(phi, dphi);
    P.SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), m);
    return acc;
  }
    

  double AcceptanceALERT(const TLorentzVector P, const char * part){
    double p = P.P() * 1000.0;//MeV
    if (p > 350.0) return 0;//sharp cut on momentum
    double theta = P.Theta() * 180.0 / M_PI;
    TH2D * acc;
    if (strcmp(part, "p") == 0) acc = acc_pip_alert;
    else if (strcmp(part, "K+") == 0) acc = acc_Kp_alert;
    else if (strcmp(part, "K-") == 0) acc = acc_Km_alert;
    else if (strcmp(part, "pi+") == 0) acc = acc_pip_alert;
    else if (strcmp(part, "pi-") == 0) acc = acc_pim_alert;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(theta);
    int biny = acc->GetYaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny);
    return result;
  }

  double Acceptance(const TLorentzVector P, const char * part){
    double acc_clas = AcceptanceCLAS12FD(P, part);
    double acc_alert = AcceptanceALERT(P, part);
    return 1.0 - (1.0 - acc_clas) * (1.0 - acc_alert);
  }

  double Smear(TLorentzVector * P, const char * part){
    double m = P->M();
    double p = P->P();
    double theta = P->Theta();
    double phi = P->Phi();
    double res[3];
    double acc = AcceptanceALERT(P[0], part);
    TH2D * resp, * restheta, * resphi;
    if (acc > 0){
      if (strcmp(part, "p") == 0){
	resp = res_proton_alert_p;
	restheta = res_proton_alert_theta;
	resphi = res_proton_alert_phi;
      }
      else if (strcmp(part, "K+") == 0){
	resp = res_Kp_alert_p;
	restheta = res_Kp_alert_theta;
	resphi = res_Kp_alert_phi;
      }
      else if (strcmp(part, "K-") == 0){
	resp = res_Km_alert_p;
	restheta = res_Km_alert_theta;
	resphi = res_Km_alert_phi;
      }
      else if (strcmp(part, "pi+") == 0){
        resp = res_pip_alert_p;
        restheta = res_pip_alert_theta;
        resphi = res_pip_alert_phi;
      }
      else if (strcmp(part, "pi-") == 0){
        resp = res_pim_alert_p;
        restheta = res_pim_alert_theta;
        resphi = res_pim_alert_phi;
      }
      else
	return 0;
      res[0] = resp->GetBinContent(resp->GetXaxis()->FindBin(theta / M_PI * 180.0), resp->GetYaxis()->FindBin(p * 1000.0)) / 100.0;
      res[1] = restheta->GetBinContent(restheta->GetXaxis()->FindBin(theta / M_PI * 180.0), restheta->GetYaxis()->FindBin(p * 1000.0));
      res[2] = resphi->GetBinContent(resphi->GetXaxis()->FindBin(theta / M_PI * 180.0), resphi->GetYaxis()->FindBin(p * 1000.0));
      p = p * abs(random.Gaus(1, res[0]));
      theta = random.Gaus(theta, res[1]);
      phi = random.Gaus(phi, res[2]);
      P->SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), m);
      return acc;
    }
    acc = AcceptanceCLAS12FD(P[0], part);
    if (acc > 0){
      res[0] = 0.01;
      res[1] = 0.001;
      res[2] = 0.004;
      p = p * abs(random.Gaus(1, res[0]));
      theta = random.Gaus(theta, res[1]);
      phi = random.Gaus(phi, res[2]);
      P->SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), m);
      return acc;
    }
    return 0;
  }


}


#endif
