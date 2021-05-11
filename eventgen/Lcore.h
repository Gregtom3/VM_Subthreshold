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

   double dSigmaJpsi_2S_2g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 7.5671e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / pow(R * M, 2);//nb GeV^-2
    return 0.16*result * 1e-7 / pow(0.197327, 2);//ds/dt in unit GeV^-4
  }

  double dSigmaJpsi_2S_23g(const double x, const double t){//Brodsky et al. PLB498 (2001) 23-28
    double N2g = 6.499e3;
    double N3g = 2.894e3;
    double v = 1.0 / (16.0 * M_PI);
    double R = 1.0;
    double M = 3.0969;//GeV
    double result = N2g * v * pow(1.0 - x, 2) * exp(1.13 * t) / (R * R * M * M) + N3g * v * exp(1.13 * t) / pow(R * M, 4);//nb GeV^-2
    return 0.16*result * 1e-7 / pow(0.197327, 2);//ds/dt in unit GeV^-4
  }

  int SetModel(const char * model = "2g"){
    if (strcmp(model, "2g") == 0)
      dSigmaJpsi = &dSigmaJpsi_2g;
    else if (strcmp(model, "23g") == 0)
      dSigmaJpsi = &dSigmaJpsi_23g;
    else if(strcmp(model, "2S_2g") == 0)
      dSigmaJpsi = &dSigmaJpsi_2S_2g;
    else if(strcmp(model, "2S_23g") == 0)
      dSigmaJpsi = &dSigmaJpsi_2S_23g;
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

namespace JPSID4d{//Harry's model for J/psi production from Deuteron 4-dim diff cross section 2020-06-14

  TRandom3 random(0);
  const double Md = 1.8756;
  const double Mn = 0.93957;
  
  TFile * fJpsiD;
  TH3D * hds72[10], * hds73[10], * hds74[10], * hds75[10], * hds76[10], * hds77[10], * hds78[10], * hds79[10], * hds80[10], * hds81[10], * hds82[10], * hds83[10], * hds84[10], * hds85[10];
  TH3D * hp72[10], * hp73[10], * hp74[10], * hp75[10], * hp76[10], * hp77[10], * hp78[10], * hp79[10], * hp80[10], * hp81[10], * hp82[10], * hp83[10], * hp84[10], * hp85[10];
  double ds72[10], ds73[10], ds74[10], ds75[10], ds76[10], ds77[10], ds78[10], ds79[10], ds80[10], ds81[10], ds82[10], ds83[10], ds84[10], ds85[10];

  double GetJpsip(const TLorentzVector q, TLorentzVector * kj){
    //q: gamma(*); kj: J/psi, p
    double k, theta_k, phi_k, p, theta_p, phi_p, ds;
    double deg = M_PI / 180.0;
    double Eg = q.E() + q * q / (2.0 * Md);
    phi_p = random.Uniform(-180.0, 180.0);
    int idx = (int) ((abs(phi_p) + 10) / 20);
    if (Eg < 7.2)
      return 0;
    else if (Eg > 7.2 && Eg < 7.25){
      hds72[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp72[idx]->GetBinContent(hp72[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.25 && Eg < 7.35){
      hds73[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp73[idx]->GetBinContent(hp73[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.35 && Eg < 7.45){
      hds74[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp74[idx]->GetBinContent(hp74[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.45 && Eg < 7.55){
      hds75[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp75[idx]->GetBinContent(hp75[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.55 && Eg < 7.65){
      hds76[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp76[idx]->GetBinContent(hp76[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.65 && Eg < 7.75){
      hds77[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp77[idx]->GetBinContent(hp77[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.75 && Eg < 7.85){
      hds78[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp78[idx]->GetBinContent(hp78[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.85 && Eg < 7.95){
      hds79[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp79[idx]->GetBinContent(hp79[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.95 && Eg < 8.05){
      hds80[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp80[idx]->GetBinContent(hp80[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 8.05 && Eg < 8.15){
      hds81[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp81[idx]->GetBinContent(hp81[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 8.15 && Eg < 8.25){
      hds82[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp82[idx]->GetBinContent(hp82[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 8.25 && Eg < 8.35){
      hds83[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp83[idx]->GetBinContent(hp83[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 8.35 && Eg < 8.45){
      hds84[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp84[idx]->GetBinContent(hp84[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 8.45 && Eg < 8.5){
      hds85[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp85[idx]->GetBinContent(hp85[idx]->FindBin(k, theta_k, theta_p));
    }
    else
      return 0;
    phi_p = phi_p * deg;
    theta_k = theta_k * deg;
    theta_p = theta_p * deg;   
    kj[0].SetXYZM(k * sin(theta_k), 0.0, k * cos(theta_k), PARTICLE::Jpsi.RandomM());
    double A = (pow(q.E() + Md - kj[0].E(), 2) + pow(Mp, 2) - pow(Mn, 2) - pow(q.P(), 2) - pow(kj[0].P(), 2) + 2.0 * q.P() * k * cos(theta_k)) / 2.0;
    double cth = cos(theta_k) * cos(theta_p) + sin(theta_k) * sin(theta_p) * cos(phi_p);
    double a = pow(q.E() + Md - kj[0].E(), 2) - pow(q.P() * cos(theta_p) - k * cth, 2);
    double b = -2.0 * A * (q.P() * cos(theta_p) - k * cth);
    double c = pow(Mp, 2) * pow(q.E() + Md - kj[0].E(), 2) - pow(A, 2);
    if (b * b - 4.0 * a * c <= 0){
      p = - b / (2.0 * a);
    }
    else {
      p = (- b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    }
    kj[1].SetXYZM(p * sin(theta_p) * cos(phi_p), p * sin(theta_p) * sin(phi_p), p * cos(theta_p), Mp);
    phi_k = random.Uniform(-M_PI, M_PI);
    kj[0].RotateZ(phi_k);
    kj[1].RotateZ(phi_k);
    kj[0].RotateY(q.Theta());
    kj[0].RotateZ(q.Phi());
    kj[1].RotateY(q.Theta());
    kj[1].RotateZ(q.Phi());
    if (Eg < 7.3){
      ds = exp((Eg - 7.2) / (7.3 - 7.2) * (log(ds73[idx]) - log(ds72[idx])) + log(ds72[idx]));
    }
    else if (Eg < 7.4){
      ds = exp((Eg - 7.3) / (7.4 - 7.3) * (log(ds74[idx]) - log(ds73[idx])) + log(ds73[idx]));
    }
    else if (Eg < 7.5){
      ds = exp((Eg - 7.4) / (7.5 - 7.4) * (log(ds75[idx]) - log(ds74[idx])) + log(ds74[idx]));
    }
    else if (Eg < 7.6){
      ds = exp((Eg - 7.5) / (7.6 - 7.5) * (log(ds76[idx]) - log(ds75[idx])) + log(ds75[idx]));
    }
    else if (Eg < 7.7){
      ds = exp((Eg - 7.6) / (7.7 - 7.6) * (log(ds77[idx]) - log(ds76[idx])) + log(ds76[idx]));
    }
    else if (Eg < 7.8){
      ds = exp((Eg - 7.7) / (7.8 - 7.7) * (log(ds78[idx]) - log(ds77[idx])) + log(ds77[idx]));
    }
    else if (Eg < 7.9){
      ds = exp((Eg - 7.8) / (7.9 - 7.8) * (log(ds79[idx]) - log(ds78[idx])) + log(ds78[idx]));
    }
    else if (Eg < 8.0){
      ds = exp((Eg - 7.9) / (8.0 - 7.9) * (log(ds80[idx]) - log(ds79[idx])) + log(ds79[idx]));
    }
    else if (Eg < 8.1){
      ds = exp((Eg - 8.0) / (8.1 - 8.0) * (log(ds81[idx]) - log(ds80[idx])) + log(ds80[idx]));
    }
    else if (Eg >= 8.1){
      ds = exp((Eg - 8.1) / (8.2 - 8.1) * (log(ds82[idx]) - log(ds81[idx])) + log(ds81[idx]));
    }
    else if (Eg >= 8.2){
      ds = exp((Eg - 8.2) / (8.3 - 8.2) * (log(ds83[idx]) - log(ds82[idx])) + log(ds82[idx]));
    }
    else if (Eg >= 8.3){
      ds = exp((Eg - 8.3) / (8.4 - 8.3) * (log(ds84[idx]) - log(ds83[idx])) + log(ds83[idx]));
    }
    else if (Eg >= 8.4){
      ds = exp((Eg - 8.4) / (8.5 - 8.4) * (log(ds85[idx]) - log(ds84[idx])) + log(ds84[idx]));
    }
    else
      return 0;
    return ds * 1.0e-7 / pow(0.197367,2);
  }
    

  int SetModel(const char * model = "v18"){
    if (strcmp(model, "v18") == 0)
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-4D-v18.root", "r");
    else {
      cout << "model not matching!" << endl;
      cout << "available options: v18" << endl;
      return 1;
    }
    for (int i = 1; i <= 10; i++){
      hds72[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.2_idx%d",i));
      hp72[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.2_idx%d",i));
      ds72[i-1]  = hds72[i-1]->Integral();
      hds73[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.3_idx%d",i));
      hp73[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.3_idx%d",i));
      ds73[i-1]  = hds73[i-1]->Integral();
      hds74[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.4_idx%d",i));
      hp74[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.4_idx%d",i));
      ds74[i-1]  = hds74[i-1]->Integral();
      hds75[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.5_idx%d",i));
      hp75[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.5_idx%d",i));
      ds75[i-1]  = hds75[i-1]->Integral();
      hds76[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.6_idx%d",i));
      hp76[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.6_idx%d",i));
      ds76[i-1]  = hds76[i-1]->Integral();
      hds77[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.7_idx%d",i));
      hp77[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.7_idx%d",i));
      ds77[i-1]  = hds77[i-1]->Integral();
      hds78[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.8_idx%d",i));
      hp78[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.8_idx%d",i));
      ds78[i-1]  = hds78[i-1]->Integral();
      hds79[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.9_idx%d",i));
      hp79[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.9_idx%d",i));
      ds79[i-1]  = hds79[i-1]->Integral();
      hds80[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.0_idx%d",i));
      hp80[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.0_idx%d",i));
      ds80[i-1]  = hds80[i-1]->Integral();
      hds81[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.1_idx%d",i));
      hp81[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.1_idx%d",i));
      ds81[i-1]  = hds81[i-1]->Integral();
      hds82[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.2_idx%d",i));
      hp82[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.2_idx%d",i));
      ds82[i-1]  = hds82[i-1]->Integral();
      hds83[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.3_idx%d",i));
      hp83[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.3_idx%d",i));
      ds83[i-1]  = hds83[i-1]->Integral();
      hds84[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.4_idx%d",i));
      hp84[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.4_idx%d",i));
      ds84[i-1]  = hds84[i-1]->Integral();
      hds85[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.5_idx%d",i));
      hp85[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.5_idx%d",i));
      ds85[i-1]  = hds85[i-1]->Integral();
    }
    return 0;
  }

}

namespace JPSID4d_old{//Harry's model for J/psi production from Deuteron 4-dim diff cross section 2020-05-14

  TRandom3 random(0);
  const double Md = 1.8756;
  const double Mn = 0.93957;
  
  TFile * fJpsiD;
  TH3D * hds72[10];
  TH3D * hds77[10];
  TH3D * hds82[10];
  TH3D * hp72[10];
  TH3D * hp77[10];
  TH3D * hp82[10];

  double GetJpsip(const TLorentzVector q, TLorentzVector * kj){
    //q: gamma(*); kj: J/psi, p
    double k, theta_k, phi_k, p, theta_p, phi_p, ds, ds72, ds77, ds82;
    double deg = M_PI / 180.0;
    double Eg = q.E() + q * q / (2.0 * Md);
    phi_p = random.Uniform(-180.0, 180.0);
    int idx = (int) ((abs(phi_p) + 10) / 20);
    if (Eg < 7.2)
      return 0;
    else if (Eg > 7.2 && Eg < 7.45){
      hds72[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp72[idx]->GetBinContent(hp72[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.45 && Eg < 7.95){
      hds77[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp77[idx]->GetBinContent(hp82[idx]->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.95 && Eg < 8.2){
      hds82[idx]->GetRandom3(k, theta_k, theta_p);
      //p = hp82[idx]->GetBinContent(hp82[idx]->FindBin(k, theta_k, theta_p));
    }
    else
      return 0;
    phi_p = phi_p * deg;
    theta_k = theta_k * deg;
    theta_p = theta_p * deg;   
    kj[0].SetXYZM(k * sin(theta_k), 0.0, k * cos(theta_k), PARTICLE::Jpsi.RandomM());
    double A = (pow(q.E() + Md - kj[0].E(), 2) + pow(Mp, 2) - pow(Mn, 2) - pow(q.P(), 2) - pow(kj[0].P(), 2) + 2.0 * q.P() * k * cos(theta_k)) / 2.0;
    double cth = cos(theta_k) * cos(theta_p) + sin(theta_k) * sin(theta_p) * cos(phi_p);
    double a = pow(q.E() + Md - kj[0].E(), 2) - pow(q.P() * cos(theta_p) - k * cth, 2);
    double b = -2.0 * A * (q.P() * cos(theta_p) - k * cth);
    double c = pow(Mp, 2) * pow(q.E() + Md - kj[0].E(), 2) - pow(A, 2);
    if (b * b - 4.0 * a * c <= 0){
      p = - b / (2.0 * a);
    }
    else {
      p = (- b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    }
    kj[1].SetXYZM(p * sin(theta_p) * cos(phi_p), p * sin(theta_p) * sin(phi_p), p * cos(theta_p), Mp);
    phi_k = random.Uniform(-M_PI, M_PI);
    kj[0].RotateZ(phi_k);
    kj[1].RotateZ(phi_k);
    kj[0].RotateY(q.Theta());
    kj[0].RotateZ(q.Phi());
    kj[1].RotateY(q.Theta());
    kj[1].RotateZ(q.Phi());
    ds72 = hds72[idx]->Integral();
    ds77 = hds77[idx]->Integral();
    ds82 = hds82[idx]->Integral();
    if (Eg < 7.7){
      ds = exp((Eg - 7.2) / (7.7 - 7.2) * (log(ds77) - log(ds72)) + log(ds72));
    }
    else if (Eg >= 7.7){
      ds = exp((Eg - 7.7) / (8.2 - 7.7) * (log(ds82) - log(ds77)) + log(ds77));
    }
    else
      return 0;
    return ds * 1.0e-7 / pow(0.197367,2);
  }
    

  int SetModel(const char * model = "default"){
    if (strcmp(model, "default") == 0)
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-4D-default.root", "r");
    else if (strcmp(model, "nv2Ia") == 0)
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-4D-nv2Ia.root", "r");
    else {
      cout << "model not matching!" << endl;
      cout << "available options: default, nv2Ia" << endl;
      return 1;
    }
    for (int i = 1; i <= 10; i++){
      hds72[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.2_idx%d",i));
      hp72[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.2_idx%d",i));
      hds77[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E7.7_idx%d",i));
      hp77[i-1] = (TH3D *) fJpsiD->Get(Form("p_E7.7_idx%d",i));
      hds82[i-1] = (TH3D *) fJpsiD->Get(Form("ds_E8.2_idx%d",i));
      hp82[i-1] = (TH3D *) fJpsiD->Get(Form("p_E8.2_idx%d",i));
    }
    return 0;
  }

}

namespace JPSID3d{//Harry's model for J/psi production from Deuteron 3-dim diff cross section 2020-06-07

  TRandom3 random(0);
  const double Md = 1.8756;
  const double Mn = 0.93957;

  TFile * fJpsiD;
  TH3D * hds72, * hds77, * hds82;
  TH3D * hp72, * hp77, * hp82;

  double GetJpsip(const TLorentzVector q, TLorentzVector * kj){
    //q: gamma(*); kj: J/psi, p
    double k, theta_k, phi_k, p, theta_p, phi_p, ds, ds72, ds77, ds82;
    double deg = M_PI / 180.0;
    double Eg = q.E() + q * q / (2.0 * Md);
    if (Eg < 7.2)
      return 0;
    else if (Eg > 7.2 && Eg < 7.45){
      hds72->GetRandom3(k, theta_k, theta_p);
      p = hp72->GetBinContent(hp72->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.45 && Eg < 7.95){
      hds77->GetRandom3(k, theta_k, theta_p);
      p = hp77->GetBinContent(hp77->FindBin(k, theta_k, theta_p));
    }
    else if (Eg > 7.95 && Eg < 8.2){
      hds82->GetRandom3(k, theta_k, theta_p);
      p = hp82->GetBinContent(hp82->FindBin(k, theta_k, theta_p));
    }
    else
      return 0;
    theta_k = theta_k * deg;
    theta_p = theta_p * deg;
    kj[0].SetXYZM(k * sin(theta_k), 0.0, k * cos(theta_k), PARTICLE::Jpsi.RandomM());
    kj[1].SetXYZM(p * sin(theta_p), 0.0, p * cos(theta_p), Mp);
    double cphip = (pow(q.E() + Md - kj[0].E() - kj[1].E(), 2) - pow(Mn, 2) - pow(q.P(), 2) - pow(kj[0].P(), 2) - pow(kj[1].P(), 2) + 2.0 * q.P() * kj[0].Pz() + 2.0 * q.P() * kj[1].Pz() - 2.0 * kj[0].Pz() * kj[1].Pz()) / (2.0 * kj[0].Pt() * kj[1].Pt());
    if (cphip > 1.0) cphip = 1.0;
    if (cphip < -1.0) cphip = -1.0;
    phi_p = acos(cphip);
    phi_k = random.Uniform(-M_PI, M_PI);
    kj[0].RotateZ(phi_k);
    kj[1].RotateZ(phi_k + phi_p);
    kj[0].RotateY(q.Theta());
    kj[0].RotateZ(q.Phi());
    kj[1].RotateY(q.Theta());
    kj[1].RotateZ(q.Phi());
    ds72 = hds72->Integral();
    ds77 = hds77->Integral();
    ds82 = hds82->Integral();
    if (Eg < 7.7){
      ds = exp((Eg - 7.2) / (7.7 - 7.2) * (log(ds77) - log(ds72)) + log(ds72));
    }
    else if (Eg >= 7.7){
      ds = exp((Eg - 7.7) / (8.2 - 7.7) * (log(ds82) - log(ds77)) + log(ds77));
    }
    else
      return 0;
    return ds * 1.0e-7 / pow(0.197367,2);
  }
    

  int SetModel(const char * model = "default"){
    if (strcmp(model, "default") == 0){
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-3D-default.root", "r");
    }
    else if (strcmp(model, "nv2Ia") == 0){
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-3D-nv2Ia.root", "r");
    }
    else if (strcmp(model, "v18") == 0){
      fJpsiD = new TFile("harrymodel/harrymodel-jpsi-d-3D-v18.root", "r");
    }
    else {
      cout << "model not matching!" << endl;
      cout << "available options: default, nv2Ia, v18" << endl;
      return 1;
    }
    hds72 = (TH3D *) fJpsiD->Get("ds_E7.2");
    hds77 = (TH3D *) fJpsiD->Get("ds_E7.7");
    hds82 = (TH3D *) fJpsiD->Get("ds_E8.2");
    hp72 = (TH3D *) fJpsiD->Get("p_E7.2");
    hp77 = (TH3D *) fJpsiD->Get("p_E7.7");
    hp82 = (TH3D *) fJpsiD->Get("p_E8.2");
    cout << "Model: " << model << endl;
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

  bool psi_2S = false;
  /* Bremsstrahlung photon */

  int fail = 0;
  
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

  int Set_psi_2S(){
    psi_2S = true;
    return 0;
  }

  int print_fail(){ //debugger
    std::cout << "The specific weight calculation has failed " << fail << " times." << std::endl;
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

  /* modified 2020-09-09 to Harry's pomeron-LQCD model, arXiv:2004.13934 */
  /* modified 2020-09-09 to Harry's pomeron-LQCD model, arXiv:2004.13934 */
  double JpsiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: gamma, N; kf: Jpsi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    double W = Pout.M();
    double MJpsi = 0.0;
    if(psi_2S==true)
      {
	MJpsi = PARTICLE::psi2S.RandomM();
      }
    else
      {
	MJpsi = PARTICLE::Jpsi.RandomM();
      }
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
    double cth = (sqrt(ki[0] * ki[0] + q * q) * sqrt(kf[0] * kf[0] + k * k) - ki[0] * kf[0]) / (q * k);
    double volume = 4.0 * M_PI;
    //  double weight = JPSIPomLQCD::dSigmaJpsi(W, cth) * volume;//ds/dOmega * volumn(4pi)
    double Jac = 2.0 * k * q / (2.0 * M_PI);
    double flux = sqrt(pow(ki[0] * ki[1], 2) - (ki[0] * ki[0]) * (ki[1] * ki[1])) / (Mp * ki[0].P());
    double weight = JPSIMODEL::dSigmaJpsi(x,t) * Jac * flux * volume;
    return weight;//GeV^-2
  }

  double cthrange[2] = {-1.0, 1.0};
  double perange[2] = {0.0, 10.0};
  double VirtualPhoton_old(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', gamma
    double Pe = random.Uniform(perange[0], perange[1]);
    double cth = random.Uniform(cthrange[0], cthrange[1]);
    double sth = sqrt(1.0 - cth * cth);
    double phi = random.Uniform(-M_PI, M_PI);
    double m = PARTICLE::e.M();
    kf[0].SetXYZT(Pe * sth * cos(phi), Pe * sth * sin(phi), Pe * cth, sqrt(m * m + Pe * Pe));//e'
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
    //    cout << "Virtual Photon Output (OLD) = " << flux * coup * tran * phas * volume << endl;
    return flux * coup * tran * phas * volume;
  }

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
    double epsilon = (1.0 - y - 0.24 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    //    double epsilon2 = (1.0)/(1.0 + 2.0 * pow(abs(kf[1].P()),2) / Q2 * pow(tan(kf[0].Theta()/2),2));
    //    cout << epsilon << " compared with " << epsilon2 << endl;
    //    cout << "Virtual Photon Output = " << couple * flux * amp * phase * volume / (1.0 - epsilon) << endl;
    return couple * flux * amp * phase * volume / (1.0 - epsilon);
  }



  double VirtualPhoton_Compare(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', gamma
    double m = PARTICLE::e.M();
    double Pe = random.Uniform(perange[0], perange[1]);
    double cth = random.Uniform(cthrange[0], cthrange[1]);
    double sth = sqrt(1.0 - cth * cth);
    double phi = random.Uniform(-M_PI, M_PI);
    if(true)
      {
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
	double epsilon = (1.0 - y - 0.24 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
	//    double epsilon2 = (1.0)/(1.0 + 2.0 * pow(abs(kf[1].P()),2) / Q2 * pow(tan(kf[0].Theta()/2),2));
	//    cout << epsilon << " compared with " << epsilon2 << endl;
	cout << "Virtual Photon Output = " << couple * flux * amp * phase * volume / (1.0 - epsilon) << endl;
      }
    if(true)
      {
	//double m = PARTICLE::e.M();
	kf[0].SetXYZT(Pe * sth * cos(phi), Pe * sth * sin(phi), Pe * cth, sqrt(m * m + Pe * Pe));//e'
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
	cout << "Virtual Photon Output (OLD) = " << flux * coup * tran * phas * volume << endl;

      }
    return 1;
    //	return couple * flux * amp * phase * volume / (1.0 - epsilon);
  }

  double JpsiElectroproduction(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', Jpsi, N'
    //    double weight1 = VirtualPhoton_Compare(ki, kf);//Generate scattered electron
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
    //    cout << "r = " << r << endl;
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 0.0;
    if(psi_2S==true)
      {
	branch = 7.73e-3;//Branch ratio to e+e- (Psi_2S)
      }
    else
      {
	branch = 5.971e-2;//Branch ratio to e+e- (J/Psi)
      }
    
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

  /* Harry's model */

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

  /* Harry's model 4d 2020-06-07, 4d 2020-06-12 */
  double Event_gD2eep_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: g; kf: [e+, e-], p
    TLorentzVector kj[2];
    double weight = JPSID4d::GetJpsip(ki[0], kj);
    if (weight == 0) return 0;
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kj[0], 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//e+
    kf[1] = *GenPhase.GetDecay(1);//e-
    kf[2] = kj[1];//p
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
    double s = pow(q.E() + 4.0 * Mp, 2) - pow(q.P(), 2);
    double Eg = (s - pow(4.0 * Mp, 2)) / (2.0 * 4.0 * Mp);
    //double Eg = q.E();
    TLorentzVector j;
    double weight2 = JPSID::GetJpsi(Eg, &j);
    if (weight2 == 0) return 0;
    j.RotateY(q.Theta());
    j.RotateZ(q.Phi());
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(j, 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    double branch = 5.971e-2;
    return weight1 * weight2 * branch;
  }

   /* Harry's model 3d 2020-06-07, 4d 2020-06-12 */
  double Event_eD2eeep_Jpsi(const TLorentzVector * ki, TLorentzVector * kf){
    //ki: e; kf: e', [e+, e-], p
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    TLorentzVector q = ki[0] - kf[0];
    TLorentzVector kj[2];
    double weight2 = JPSID4d::GetJpsip(q, kj);
    if (weight2 == 0) return 0;
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kj[0], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//e+
    kf[2] = *GenPhase.GetDecay(1);//e-
    kf[3] = kj[1];//p
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
      facc_solid = new TFile("acceptance/acceptance_solid_JPsi_electron_target315_output.root", "r");
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
