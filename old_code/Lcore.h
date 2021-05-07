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
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"

#include "Lparticle.h"

using namespace std;

const double Mp = PARTICLE::proton.M();

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

namespace JPSIMODEL{//Model of J/psi production

  double (*dSigmaJpsi)(const double, const double);
  
  double dSigmaJpsi_2g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 7.5671e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / pow(R * M, 2);//nb GeV^-2
    return result * 1e-7 / pow(0.197327, 2);//ds/dt in unit GeV^-4
  }

  double dSigmaJpsi_23g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 6.499e3;
    double N3g = 2.894e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / (R * R * M * M) + N3g * v * exp(1.13 * t) / pow(R * M, 4);//nb GeV^-2
    return result * 1e-7 / pow(0.197327, 2);//ds/dt in unit GeV^-4
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

namespace PHIMODEL{//Model of phi production

  double (* dSigmaPhi)(const double, const double);

  double dSigmaPhi_fit(const double dW, const double cth){
    if (dW <= 0) return 0;
    double par[5] = {0.232612, 1.95038, 4.02454, 1.52884, 0.525636};
    double b1 = par[3] * pow(dW * (dW + 2.0 * (Mp + PARTICLE::phi.M())), par[4]);
    double a2 = par[2] * dW;
    double a0 = par[0] * atan(par[1] * par[1] * (dW * (dW + 2.0 * (Mp + PARTICLE::phi.M()))));
    double r0 = exp(b1 * cth) * b1 / (2.0 * sinh(b1));
    double r2 = cth * cth * exp(b1 * cth) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
    double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
    return ds;//GeV^-2
  }

  int SetModel(const char * model = "fit"){
    if (strcmp(model, "fit") == 0)
      dSigmaPhi = &dSigmaPhi_fit;
    else {
      cout << "No matching model! Set to fit model!" << endl;
      dSigmaPhi = &dSigmaPhi_fit;
    }
    return 0;
  }

}


namespace JPSIHE4{//Harry's model for J/psi production from He4

  TRandom3 random(0);

  TFile * fJpsiHe4 = new TFile("harrymodel/harrymodel-jpsi-cut.root", "r");
  TH2D * hds[11];
  TH3D * hpmin, * hpmax;
  double Egarray[11] = {6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2};
  double logsarray[11];
  ROOT::Math::Interpolator Eglogs(11, ROOT::Math::Interpolation::kCSPLINE);
  
  double sigma(const double Eg){
    double result = 1.0e-7 / pow(0.197367,2) * exp(Eglogs.Eval(Eg));
    return result;// in unit GeV^-2
  }

  double GetJpsi(const double Eg, TLorentzVector * kj){
    double k, theta, phi;
    if (Eg < 6.2)
      return 0;
    else if (Eg > 6.2 && Eg < 6.3)
      hds[0]->GetRandom2(k, theta);
    else if (Eg > 6.3 && Eg < 6.5)
      hds[1]->GetRandom2(k, theta);
    else if (Eg > 6.5 && Eg < 6.7)
      hds[2]->GetRandom2(k, theta);
    else if (Eg > 6.7 && Eg < 6.9)
      hds[3]->GetRandom2(k, theta);
    else if (Eg > 6.9 && Eg < 7.1)
      hds[4]->GetRandom2(k, theta);
    else if (Eg > 7.1 && Eg < 7.3)
      hds[5]->GetRandom2(k, theta);
    else if (Eg > 7.3 && Eg < 7.5)
      hds[6]->GetRandom2(k, theta);
    else if (Eg > 7.5 && Eg < 7.7)
      hds[7]->GetRandom2(k, theta);
    else if (Eg > 7.7 && Eg < 7.9)
      hds[8]->GetRandom2(k, theta);
    else if (Eg > 7.9 && Eg < 8.1)
      hds[9]->GetRandom2(k, theta);
    else if (Eg > 8.1)
      hds[10]->GetRandom2(k, theta);
    else
      return 0;
    theta = theta / 180.0 * M_PI;
    phi = random.Uniform(-M_PI, M_PI);
    kj->SetXYZM(k * sin(theta) * cos(phi), k * sin(theta) * sin(phi), k * cos(theta), PARTICLE::Jpsi.RandomM());
    return sigma(Eg);
  }
    

  int SetModel(const char * model = "cut1"){
    if (strcmp(model, "cut1") == 0){
      hds[0] = (TH2D *) fJpsiHe4->Get("cut1_E=6.2");
      hds[1] = (TH2D *) fJpsiHe4->Get("cut1_E=6.4");
      hds[2] = (TH2D *) fJpsiHe4->Get("cut1_E=6.6");
      hds[3] = (TH2D *) fJpsiHe4->Get("cut1_E=6.8");
      hds[4] = (TH2D *) fJpsiHe4->Get("cut1_E=7.0");
      hds[5] = (TH2D *) fJpsiHe4->Get("cut1_E=7.2");
      hds[6] = (TH2D *) fJpsiHe4->Get("cut1_E=7.4");
      hds[7] = (TH2D *) fJpsiHe4->Get("cut1_E=7.6");
      hds[8] = (TH2D *) fJpsiHe4->Get("cut1_E=7.8");
      hds[9] = (TH2D *) fJpsiHe4->Get("cut1_E=8.0");
      hds[10] = (TH2D *) fJpsiHe4->Get("cut1_E=8.2");
      hpmin = (TH3D *) fJpsiHe4->Get("cut1_pmin");
      hpmax = (TH3D *) fJpsiHe4->Get("cut1_pmax");
    }
    else if (strcmp(model, "cut2") == 0){
      hds[0] = (TH2D *) fJpsiHe4->Get("cut2_E=6.2");
      hds[1] = (TH2D *) fJpsiHe4->Get("cut2_E=6.4");
      hds[2] = (TH2D *) fJpsiHe4->Get("cut2_E=6.6");
      hds[3] = (TH2D *) fJpsiHe4->Get("cut2_E=6.8");
      hds[4] = (TH2D *) fJpsiHe4->Get("cut2_E=7.0");
      hds[5] = (TH2D *) fJpsiHe4->Get("cut2_E=7.2");
      hds[6] = (TH2D *) fJpsiHe4->Get("cut2_E=7.4");
      hds[7] = (TH2D *) fJpsiHe4->Get("cut2_E=7.6");
      hds[8] = (TH2D *) fJpsiHe4->Get("cut2_E=7.8");
      hds[9] = (TH2D *) fJpsiHe4->Get("cut2_E=8.0");
      hds[10] = (TH2D *) fJpsiHe4->Get("cut2_E=8.2");
      hpmin = (TH3D *) fJpsiHe4->Get("cut2_pmin");
      hpmax = (TH3D *) fJpsiHe4->Get("cut2_pmax");
    }
    else
      return 0;
    for (int i = 0; i < 11; i++)
      logsarray[i] = log(hds[i]->Integral() * 0.01 * 2.0 / 180.0 * M_PI + 1e-50);
    Eglogs.SetData(11, Egarray, logsarray);
    return 0;
  }

}

namespace JPSID{//Harry's model for J/psi production from Deuteron

  TRandom3 random(0);

  TFile * fJpsiD = new TFile("harrymodel/harrymodel-2h-jpsi-cut.root", "r");
  TH2D * hds[6];
  TH3D * hpmin, * hpmax;
  double Egarray[6] = {7.2, 7.4, 7.6, 7.8, 8.0, 8.2};
  double logsarray[6];
  ROOT::Math::Interpolator Eglogs(6, ROOT::Math::Interpolation::kCSPLINE);
  
  double sigma(const double Eg){
    double result = 1.0e-7 / pow(0.197367,2) * exp(Eglogs.Eval(Eg));
    return result;// in unit GeV^-2
  }

  double GetJpsi(const double Eg, TLorentzVector * kj){
    double k, theta, phi;
    if (Eg < 7.2)
      return 0;
    else if (Eg > 7.2 && Eg < 7.3)
      hds[0]->GetRandom2(k, theta);
    else if (Eg > 7.3 && Eg < 7.5)
      hds[1]->GetRandom2(k, theta);
    else if (Eg > 7.5 && Eg < 7.7)
      hds[2]->GetRandom2(k, theta);
    else if (Eg > 7.7 && Eg < 7.9)
      hds[3]->GetRandom2(k, theta);
    else if (Eg > 7.9 && Eg < 8.1)
      hds[4]->GetRandom2(k, theta);
    else if (Eg > 8.1)
      hds[5]->GetRandom2(k, theta);
    else
      return 0;
    theta = theta / 180.0 * M_PI;
    phi = random.Uniform(-M_PI, M_PI);
    kj->SetXYZM(k * sin(theta) * cos(phi), k * sin(theta) * sin(phi), k * cos(theta), PARTICLE::Jpsi.RandomM());
    return sigma(Eg);
  }
    

  int SetModel(const char * model = "cut1"){
    if (strcmp(model, "cut1") == 0){
      hds[0] = (TH2D *) fJpsiD->Get("cut1_E=7.2");
      hds[1] = (TH2D *) fJpsiD->Get("cut1_E=7.4");
      hds[2] = (TH2D *) fJpsiD->Get("cut1_E=7.6");
      hds[3] = (TH2D *) fJpsiD->Get("cut1_E=7.8");
      hds[4] = (TH2D *) fJpsiD->Get("cut1_E=8.0");
      hds[5] = (TH2D *) fJpsiD->Get("cut1_E=8.2");
    }
    else if (strcmp(model, "cut2") == 0){
      hds[0] = (TH2D *) fJpsiD->Get("cut2_E=7.2");
      hds[1] = (TH2D *) fJpsiD->Get("cut2_E=7.4");
      hds[2] = (TH2D *) fJpsiD->Get("cut2_E=7.6");
      hds[3] = (TH2D *) fJpsiD->Get("cut2_E=7.8");
      hds[4] = (TH2D *) fJpsiD->Get("cut2_E=8.0");
      hds[5] = (TH2D *) fJpsiD->Get("cut2_E=8.2");
    }
    else
      return 0;
    for (int i = 0; i < 6; i++)
      logsarray[i] = log(hds[i]->Integral() * 0.01 * 2.0 / 180.0 * M_PI + 1e-50);
    Eglogs.SetData(6, Egarray, logsarray);
    return 0;
  }

}

namespace JPSIHE4_old{//Harry's model for J/psi production from He4

  TRandom3 random(0);

  TFile * fJpsiHe4 = new TFile("harrymodel/harrymodel-jpsi.root", "r");
  TH2D * hds[11];
  TH3D * hpmin, * hpmax;
  double Egarray[11] = {6.2, 6.5, 6.8, 7.1, 7.4, 7.7, 8.0, 8.3, 8.6, 8.9, 9.2};
  double logsarray[11];
  ROOT::Math::Interpolator Eglogs(11, ROOT::Math::Interpolation::kCSPLINE);
  
  double sigma(const double Eg){
    double result = 1.0e-7 / pow(0.197367,2) * exp(Eglogs.Eval(Eg));
    return result;// in unit GeV^-2
  }

  double GetJpsi(const double Eg, TLorentzVector * kj){
    double k, theta, phi;
    if (Eg < 6.2)
      return 0;
    else if (Eg > 6.2 && Eg < 6.35)
      hds[0]->GetRandom2(k, theta);
    else if (Eg > 6.35 && Eg < 6.65)
      hds[1]->GetRandom2(k, theta);
    else if (Eg > 6.65 && Eg < 6.95)
      hds[2]->GetRandom2(k, theta);
    else if (Eg > 6.95 && Eg < 7.25)
      hds[3]->GetRandom2(k, theta);
    else if (Eg > 7.25 && Eg < 7.55)
      hds[4]->GetRandom2(k, theta);
    else if (Eg > 7.55 && Eg < 7.85)
      hds[5]->GetRandom2(k, theta);
    else if (Eg > 7.85 && Eg < 8.15)
      hds[6]->GetRandom2(k, theta);
    else if (Eg > 8.15 && Eg < 8.45)
      hds[7]->GetRandom2(k, theta);
    else if (Eg > 8.45 && Eg < 8.75)
      hds[8]->GetRandom2(k, theta);
    else if (Eg > 8.75 && Eg < 9.05)
      hds[9]->GetRandom2(k, theta);
    else if (Eg > 9.05)
      hds[10]->GetRandom2(k, theta);
    else
      return 0;
    theta = theta / 180.0 * M_PI;
    phi = random.Uniform(-M_PI, M_PI);
    kj->SetXYZM(k * sin(theta) * cos(phi), k * sin(theta) * sin(phi), k * cos(theta), PARTICLE::Jpsi.RandomM());
    return sigma(Eg);
  }
    

  int SetModel(const char * model = "all"){
    if (strcmp(model, "all") == 0){
      hds[0] = (TH2D *) fJpsiHe4->Get("ds_E=6.2");
      hds[1] = (TH2D *) fJpsiHe4->Get("ds_E=6.5");
      hds[2] = (TH2D *) fJpsiHe4->Get("ds_E=6.8");
      hds[3] = (TH2D *) fJpsiHe4->Get("ds_E=7.1");
      hds[4] = (TH2D *) fJpsiHe4->Get("ds_E=7.4");
      hds[5] = (TH2D *) fJpsiHe4->Get("ds_E=7.7");
      hds[6] = (TH2D *) fJpsiHe4->Get("ds_E=8.0");
      hds[7] = (TH2D *) fJpsiHe4->Get("ds_E=8.3");
      hds[8] = (TH2D *) fJpsiHe4->Get("ds_E=8.6");
      hds[9] = (TH2D *) fJpsiHe4->Get("ds_E=8.9");
      hds[10] = (TH2D *) fJpsiHe4->Get("ds_E=9.2");
      hpmin = (TH3D *) fJpsiHe4->Get("pmin");
      hpmax = (TH3D *) fJpsiHe4->Get("pmax");
    }
    else if (strcmp(model, "c300") == 0){
      hds[0] = (TH2D *) fJpsiHe4->Get("ds_E=6.2_c300");
      hds[1] = (TH2D *) fJpsiHe4->Get("ds_E=6.5_c300");
      hds[2] = (TH2D *) fJpsiHe4->Get("ds_E=6.8_c300");
      hds[3] = (TH2D *) fJpsiHe4->Get("ds_E=7.1_c300");
      hds[4] = (TH2D *) fJpsiHe4->Get("ds_E=7.4_c300");
      hds[5] = (TH2D *) fJpsiHe4->Get("ds_E=7.7_c300");
      hds[6] = (TH2D *) fJpsiHe4->Get("ds_E=8.0_c300");
      hds[7] = (TH2D *) fJpsiHe4->Get("ds_E=8.3_c300");
      hds[8] = (TH2D *) fJpsiHe4->Get("ds_E=8.6_c300");
      hds[9] = (TH2D *) fJpsiHe4->Get("ds_E=8.9_c300");
      hds[10] = (TH2D *) fJpsiHe4->Get("ds_E=9.2_c300");
      hpmin = (TH3D *) fJpsiHe4->Get("pmin_c300");
      hpmax = (TH3D *) fJpsiHe4->Get("pmax_c300");
    }
    else
      return 0;
    for (int i = 0; i < 11; i++)
      logsarray[i] = log(hds[i]->Integral() * 0.01 * 2.0 / 180.0 * M_PI + 1e-50);
    Eglogs.SetData(11, Egarray, logsarray);
    return 0;
  }

}

namespace PHIHE4{//Harry's model for phi production form He4

  TFile * fphiHe4 = new TFile("harrymodel/harrymodel-phi.root", "r");
  TH2D * hE13 = (TH2D *) fphiHe4->Get("E=1.3");
  TH2D * hE14 = (TH2D *) fphiHe4->Get("E=1.4");
  TH2D * hE15 = (TH2D *) fphiHe4->Get("E=1.5");
  TH2D * hE16 = (TH2D *) fphiHe4->Get("E=1.6");
  TH2D * hE17 = (TH2D *) fphiHe4->Get("E=1.7");
  TH2D * hE18 = (TH2D *) fphiHe4->Get("E=1.8");
  TH2D * hE19 = (TH2D *) fphiHe4->Get("E=1.9");
  TH2D * hE20 = (TH2D *) fphiHe4->Get("E=2.0");

  double sE13 = hE13->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE14 = hE14->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE15 = hE15->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE16 = hE16->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE17 = hE17->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE18 = hE18->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE19 = hE19->Integral() * 0.005 * 2.0 / 180.0 * M_PI;
  double sE20 = hE20->Integral() * 0.005 * 2.0 / 180.0 * M_PI;

  double Egarray[8] = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
  double logsarray[8] = {log(sE13), log(sE14), log(sE15), log(sE16), log(sE17), log(sE18), log(sE19), log(sE20)};
  ROOT::Math::Interpolator Eglogs(8, ROOT::Math::Interpolation::kCSPLINE);

  double sigma(const double Eg){
    double result = 1.0e-7 / pow(0.197367,2) * exp(Eglogs.Eval(Eg));
    return result;// in unit GeV^-2
  }
  
  int SetModel(){
    Eglogs.SetData(8, Egarray, logsarray);
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

  /* Bremsstrahlung photon */
  
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
    return 0.01 * (4.0 / 3.0 * log(ymax / ymin) - 4.0 / 3.0 * (ymax - ymin) + 1.0 / 2.0 * (ymax * ymax - ymin * ymin));
  }

  int SetBremsstrahlung(){
    TF_fBremsstrahlung = new TF1("fBremsstrahlung", Bremsstrahlung, 0.01, 1.0, 0);
    TF_fBremsstrahlung->SetNpx(1000);
    return 0;
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

  /* J/psi productions */
  
  double JpsiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: Jpsi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    double W = Pout.M();
    double MJpsi = PARTICLE::Jpsi.RandomM();
    if (W < MJpsi + Mp) return 0;//below the threshold
    double x = (2.0 * Mp * MJpsi + MJpsi * MJpsi) / (W * W - Mp * Mp);
    double mass[2] = {MJpsi, Mp};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//Jpsi
    kf[1] = *GenPhase.GetDecay(1);//N'
    double t = (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double k = sqrt(pow(W * W - kf[0] * kf[0] - kf[1] * kf[1], 2) - 4.0 * (kf[0] * kf[0]) * (kf[1] * kf[1])) / (2.0 * W);
    double q = sqrt(pow(W * W - ki[0] * ki[0] - ki[1] * ki[1], 2) - 4.0 * (ki[0] * ki[0]) * (ki[1] * ki[1])) / (2.0 * W);
    double Jac = 2.0 * k * q / (2.0 * M_PI);
    double flux = sqrt(pow(ki[0] * ki[1], 2) - (ki[0] * ki[0]) * (ki[1] * ki[1])) / (Mp * ki[0].P());
    double volume = 4.0 * M_PI;
    double weight = JPSIMODEL::dSigmaJpsi(x, t) * Jac * flux * volume;//ds/dOmega * volumn(4pi)
    return weight;//GeV^-2
  }

  double cthrange[2] = {-1.0, 1.0};
  double perange[2] = {0.0, 10.0};
  double VirtualPhoton(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', gamma
    double Pe = random.Uniform(perange[0], perange[1]);
    double cth = random.Uniform(cthrange[0], cthrange[1]);
    double sth = sqrt(1.0 - cth * cth);
    double phi = random.Uniform(-M_PI, M_PI);
    kf[0].SetXYZT(Pe * sth * cos(phi), Pe * sth * sin(phi), Pe * cth, sqrt(PARTICLE::e.M() * PARTICLE::e.M() + Pe * Pe));//e'
    kf[1] = ki[0] - kf[0];//photon  
    double W2 = (kf[1] + ki[1]) * (kf[1] + ki[1]);
    if (W2 < Mp * Mp) return 0;//below the lowest state
    double Q2 = - kf[1] * kf[1];
    if (Q2 < 1e-5) return 0;//small Q2 cut
    double alpha_em = 1.0 / 137.0;
    double s = (ki[0] + ki[1]) * (ki[0] + ki[1]);//(P + l)^2
    if (s < W2) return 0;//energy conservation
    double y = (ki[1] * kf[1]) / (ki[1] * ki[0]);//(P * q) / (P * l)
    if (y < 0.0 || y > 1.0) return 0;//unphysical transfer
    double M2 = ki[1] * ki[1];
    double flux = sqrt(pow(W2 - M2 + Q2, 2) + 4.0 * M2 * Q2) / (s - M2);
    double coup = alpha_em / (M_PI * M_PI * Q2);
    double tran = (1.0 - y + 0.5 * y * y + M2 * Q2 / pow(s - M2, 2)) / (y * y + 4.0 * M2 * Q2 / pow(s - M2, 2));
    double phas = kf[0].P();
    double volume = 2.0 * M_PI * abs(perange[1] - perange[0]) * abs(cthrange[1] - cthrange[0]);
    cout << "Virtual Photon Output = " << flux * coup * tran * phas * volume << endl;
    return flux * coup * tran * phas * volume;
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
    TLorentzVector Pc = kf1[0] + kf1[1];//total V+N
    TLorentzVector PV = kf1[0];//V
    TLorentzVector P1 = kf[1];//e+
    PV.Boost(-Pc.BoostVector());
    P1.Boost(-Pc.BoostVector());
    P1.Boost(-PV.BoostVector());
    double theta = P1.Angle(PV.Vect());
    double r = 0.29;
    double wth = 3.0 / 8.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cos(theta),2));
    double branch = 5.971e-2;//Branch ratio to e+e-
    return weight * wth * branch;
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
    TLorentzVector Pc = kf1[1] + kf1[2];//total V+N
    TLorentzVector PV = kf1[1];//V
    TLorentzVector P1 = kf[2];//e+
    PV.Boost(-Pc.BoostVector());
    P1.Boost(-Pc.BoostVector());
    P1.Boost(-PV.BoostVector());
    double theta = P1.Angle(PV.Vect());
    double r = 0.29;
    double wth = 3.0 / 8.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cos(theta),2));
    double branch = 5.971e-2;//Branch ratio to e+e-
    return weight * wth * branch;
  }

  /* phi productions */

  double PhiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: phi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    double W = Pout.M();
    double Mphi = PARTICLE::phi.RandomM();
    if (W < Mphi + Mp) return 0;//below the threshold
    double mass[2] = {Mphi, Mp};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//phi
    kf[1] = *GenPhase.GetDecay(1);//N'
    double dW = W - Mphi - Mp;
    double k = sqrt(pow(W * W - kf[0] * kf[0] - kf[1] * kf[1], 2) - 4.0 * (kf[0] * kf[0]) * (kf[1] * kf[1])) / (2.0 * W);
    double q = sqrt(pow(W * W - ki[0] * ki[0] - ki[1] * ki[1], 2) - 4.0 * (ki[0] * ki[0]) * (ki[1] * ki[1])) / (2.0 * W);
    double cth = (sqrt(ki[0] * ki[0] + q * q) * sqrt(kf[0] * kf[0] + k * k) - ki[0] * kf[0]) / (q * k);
    double flux = sqrt(pow(ki[0] * ki[1], 2) - (ki[0] * ki[0]) * (ki[1] * ki[1])) / (Mp * ki[0].P());
    double volume = 4.0 * M_PI;
    double weight = PHIMODEL::dSigmaPhi(dW, cth) * flux * volume;//ds/dOmega * volume(4pi)
    return weight;//GeV^-2
  }

  double PhiElectroproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', phi, N'
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector ki2[2] = {kf[1], ki[1]};//Initial state: virtual photon N
    double weight2 = PhiPhotoproduction(ki2, &kf[1]);//Generate phi N' from virtual photon production
    return weight1 * weight2;
  }

  double Event_gN2Nee_Phi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: N', [e+, e-]
    TLorentzVector kf1[2];//phi, N'
    double weight = PhiPhotoproduction(ki, kf1);
    kf[0] = kf1[1];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[0], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double branch = 2.973e-4;//Branch ratio to e+e-
    return weight * branch;
  }

  double Event_eN2eNee_Phi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', N', [e+, e-]
    TLorentzVector kf1[3];//e', phi, N'
    double weight = PhiElectroproduction(ki, kf1);
    kf[0] = kf1[0];//e'
    kf[1] = kf1[2];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//e+
    kf[3] = *GenPhase.GetDecay(1);//e-
    double branch = 2.973e-4;
    return weight * branch;
  }

  double Event_g4He2ee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: g; kf: [e+, e-]
    double Eg = ki[0].E();
    TLorentzVector j;
    double weight = JPSIHE4::GetJpsi(Eg, &j);
    if (weight == 0) return 0;
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//e+
    kf[1] = *GenPhase.GetDecay(1);//e-
    double branch = 5.971e-2;
    return weight * branch;
  }

  double Event_gD2ee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: g; kf: [e+, e-]
    double Eg = ki[0].E();
    TLorentzVector j;
    double weight = JPSID::GetJpsi(Eg, &j);
    if (weight == 0) return 0;
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//e+
    kf[1] = *GenPhase.GetDecay(1);//e-
    double branch = 5.971e-2;
    return weight * branch;
  }

  double Event_e4He2eee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e; kf: e', [e+, e-]
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector q = ki[0] - kf[0];
    //double s = pow(q.E() + 4.0 * Mp, 2) - pow(q.P(), 2);
    //double Eg = (s - pow(4.0 * Mp, 2)) / (2.0 * 4.0 * Mp);
    double Eg = q.E();
    TLorentzVector j;
    double weight2 = JPSIHE4::GetJpsi(Eg, &j);
    if (weight2 == 0) return 0;
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double branch = 5.971e-2;
    return weight1 * weight2 * branch;
  }

  double Event_eD2eee_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e; kf: e', [e+, e-]
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector q = ki[0] - kf[0];
    //double s = pow(q.E() + 4.0 * Mp, 2) - pow(q.P(), 2);
    //double Eg = (s - pow(4.0 * Mp, 2)) / (2.0 * 4.0 * Mp);
    double Eg = q.E();
    TLorentzVector j;
    double weight2 = JPSID::GetJpsi(Eg, &j);
    if (weight2 == 0) return 0;
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double branch = 5.971e-2;
    return weight1 * weight2 * branch;
  }

  double Event_g4He2ee_Phi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: g; kf: [e+, e-]
    double Eg = ki[0].E();
    double k, theta, phi;
    if (Eg > 1.95){
      PHIHE4::hE20->GetRandom2(k, theta);
    }
    else if (Eg > 1.85){
      PHIHE4::hE19->GetRandom2(k, theta);
    }
    else if (Eg > 1.75){
      PHIHE4::hE18->GetRandom2(k, theta);
    }
    else if (Eg > 1.65){
      PHIHE4::hE17->GetRandom2(k, theta);
    }
    else if (Eg > 1.55){
      PHIHE4::hE16->GetRandom2(k, theta);
    }
    else if (Eg > 1.45){
      PHIHE4::hE15->GetRandom2(k, theta);
    }
    else if (Eg > 1.35){
      PHIHE4::hE14->GetRandom2(k, theta);
    }
    else if (Eg > 1.3){
      PHIHE4::hE13->GetRandom2(k, theta);
    }
    else
      return 0;
    theta = theta / 180.0 * M_PI;
    phi = random.Uniform(-M_PI, M_PI);
    TLorentzVector j;
    j.SetXYZM(k * sin(theta) * cos(phi), k * sin(theta) * sin(phi), k * cos(theta), PARTICLE::phi.RandomM());
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//e+
    kf[1] = *GenPhase.GetDecay(1);//e-
    double branch = 2.973e-4;
    double weight = PHIHE4::sigma(Eg);
    return weight * branch;
  }

  double Event_e4He2eee_Phi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e; kf: e', [e+, e-]
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector q = ki[0] - kf[0];
    double s = pow(q.E() + 4.0 * Mp, 2) - pow(q.P(), 2);
    double Eg = (s - pow(4.0 * Mp, 2)) / (2.0 * 4.0 * Mp);
    double k, theta, phi;
    if (Eg > 1.95){
      PHIHE4::hE20->GetRandom2(k, theta);
    }
    else if (Eg > 1.85){
      PHIHE4::hE19->GetRandom2(k, theta);
    }
    else if (Eg > 1.75){
      PHIHE4::hE18->GetRandom2(k, theta);
    }
    else if (Eg > 1.65){
      PHIHE4::hE17->GetRandom2(k, theta);
    }
    else if (Eg > 1.55){
      PHIHE4::hE16->GetRandom2(k, theta);
    }
    else if (Eg > 1.45){
      PHIHE4::hE15->GetRandom2(k, theta);
    }
    else if (Eg > 1.35){
      PHIHE4::hE14->GetRandom2(k, theta);
    }
    else if (Eg > 1.3){
      PHIHE4::hE13->GetRandom2(k, theta);
    }
    else
      return 0;
    theta = theta / 180.0 * M_PI;
    phi = random.Uniform(-M_PI, M_PI);
    TLorentzVector j;
    j.SetXYZM(k * sin(theta) * cos(phi), k * sin(theta) * sin(phi), k * cos(theta), PARTICLE::phi.RandomM());
    j.RotateY(ki[0].Theta());
    j.RotateZ(ki[0].Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double branch = 2.973e-4;
    double weight2 = PHIHE4::sigma(Eg);
    return weight1 * weight2 * branch;
  }

    


}


namespace DETECTOR{

  TRandom3 random(0);

  TFile * facc_clas, * facc_solid, * facc_alert;
  TFile * fres_clas, * fres_solid, * fres_alert1, * fres_alert2;
  TH3F * acc_ele_clas, * acc_pos_clas, * acc_pip_clas, * acc_pim_clas, * acc_Kp_clas, * acc_Km_clas, * acc_proton_clas;
  TH2D * acc_ele_solid, * acc_pos_solid, * acc_proton_solid;
  TH2D * acc_proton_alert, * acc_Kp_alert, * acc_pip_alert, * acc_Km_alert, * acc_pim_alert;
  TH2D * res_Kp_alert_p, * res_Kp_alert_theta, * res_Kp_alert_phi;
  TH2D * res_Km_alert_p, * res_Km_alert_theta, * res_Km_alert_phi;
  TH2D * res_pip_alert_p, * res_pip_alert_theta, * res_pip_alert_phi;
  TH2D * res_pim_alert_p, * res_pim_alert_theta, * res_pim_alert_phi;
  TH2D * res_proton_alert_p, * res_proton_alert_theta, * res_proton_alert_phi;

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
      facc_solid = new TFile("acceptance/acceptance_solid_JPsi_electron_target315_output.root", "r");
      acc_ele_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
      acc_pos_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
      acc_proton_solid = (TH2D *) facc_solid->Get("acceptance_ThetaP_overall");
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
    else if (strcmp(part, "p") == 0) acc = acc_proton_solid;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(theta);
    int biny = acc->GetYaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny);
    return result;
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
