#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_integration.h>
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

  int SetNuclear(const char * nuclear = "p"){
    flag = 1;
    if (strcmp(nuclear, "D") == 0){
      flag = 2;
      fMomentum = &Momentum_D;
    }
    else {
      cout << "No matching nuclear! Set to proton!" << endl;
      flag = 0;
    }
    return 0;
  }

  
}
namespace UPSILONMODEL{//Model of Upsilon (1S) production 
                       //Heavily inspired by the source code of the lAger event generator (see https://eicweb.phy.anl.gov/monte_carlo/lager/-/blob/master/src/lager/gen/lA/oleksii_2vmp.cc)
                       // Relevant source paper https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.014016
  int vflag = 1;
  const double Mv = 9.46030; //https://pdg.lbl.gov/2014/listings/rpp2014-list-upsilon-1S.pdf
  const double Mp = 0.938272;

  const double C_el = 13.8e-3;
  const double nu_el = 8.88;
  const double b_el = 1.27;
  const double a_el = 1.38;

  const double C_inel = 18.7;
  const double nu_inel = 20.90;
  const double b_inel = 3.53;
  const double a_inel = 1.2;

  const double f = 0.238;
  const double alpha = 1./137.;
  const double e = sqrt(4.0*M_PI*alpha); // electric charge

  
  const double epsilon = 0.0001;
  const double numax = 2000000.0; // maximum bound for infinite integration 

  const double khbar = 4.135667662e-15 / TMath::TwoPi();
  // Speed of light (m/s)
  const double kc = 299792458.;
  // hhbarc (fm * GeV)
  const double khbarc = khbar * kc * 1e6;
  // khbarc2 (nb * GeV^2, using 1fm^2 = 10mb^2 =1e7nb^2)
  const double khbarc2 = khbarc * khbarc * 1e7;
  
  TF1 *_f1 = 0;
  
  // gsl workspace for integration of Dispersion integral
  gsl_integration_workspace * _w = gsl_integration_workspace_alloc(2048);
  gsl_function _F;


  // Event variables
  double eventb0=0.0;
  double eventB=0.0;
  double eventA=0.0;
  
  double ImaginaryT(const double s){
    double nu = 0.5*(s-Mp*Mp-Mv*Mv);
    double Disc_el = 0.0;
    double Disc_inel = 0.0;
    if(nu>nu_el)
      Disc_el = C_el*pow(1.0-nu_el/nu,b_el)*pow(nu/nu_el,a_el);
    if(nu>nu_inel)
      Disc_inel = C_inel*pow(1.0-nu_inel/nu,b_inel)*pow(nu/nu_inel,a_inel);
    return Disc_el+Disc_inel;
  }

  double ImaginaryT_nu(const double nu){
    double Disc_el = 0.0;
    double Disc_inel = 0.0;
    if(nu>nu_el)
      Disc_el = C_el*pow(1.0-nu_el/nu,b_el)*pow(nu/nu_el,a_el);
    if(nu>nu_inel)
      Disc_inel = C_inel*pow(1.0-nu_inel/nu,b_inel)*pow(nu/nu_inel,a_inel);
    return Disc_el+Disc_inel;
  }

  double dispersion_integral(double x, void *p)
  {
    double nuPRIME = x;
    double nu = *(double *) p;
    double s = 2*nu + Mp*Mp + Mv*Mv;
    return (ImaginaryT_nu(nuPRIME)/nuPRIME - ImaginaryT_nu(nu)/nu)/((nuPRIME*nuPRIME-nu*nu));
  }

  double RealT(const double s){
    double nu = 0.5*(s-Mp*Mp-Mv*Mv);
    double T0 = 20.5; // Can be set to 0, 20.5, or 87, see paper for details
    _F.params = &nu;
    double result, error;
    gsl_integration_qagiu(&_F, nu_el, 0, 1e-8, 2048, _w, &result, &error);
    result+= ImaginaryT(nu)/nu * std::log(std::fabs((nu_el+nu)/(nu_el-nu))) / (2 * nu);
    return T0 + 2.0/M_PI * nu*nu * result; 
  }

  double B(double *xx, double *par)
  {
    double b = xx[0];
    double tmin = par[0];
    double tmax = par[1];
    double sigma = par[2];
    double A = par[3];
    return (b - A*exp(b*tmax)/sigma + A*exp(b*tmin)/sigma);
  }
  double (*dSigmaY1S)(const double, const double); // W , t

  double tmax(const double Mv, const double Mp, const double W)
  {
    double s = W*W;
    double qvp = sqrt((0.25/s)*(s-pow(Mv+Mp,2))*(s-pow(Mv-Mp,2)));
    double qgp = (s-Mp*Mp)/(2.0*sqrt(s));
    return Mv * Mv -2.0 * qgp * ( sqrt ( qvp * qvp + Mv * Mv ) - qvp ); // less negative
  }

  double tmin(const double Mv, const double Mp, const double W)
  {
    double s = W*W;
    double qvp = sqrt((0.25/s)*(s-pow(Mv+Mp,2))*(s-pow(Mv-Mp,2)));
    double qgp = (s-Mp*Mp)/(2.0*sqrt(s));
    return Mv * Mv -2.0 * qgp * ( sqrt ( qvp * qvp + Mv * Mv ) + qvp ); // more negative
  }
    
  void set_B(const double Mv, const double Mp, const double W)
  {
    double s = W*W;
    double nu = 0.5*(s-Mp*Mp-Mv*Mv);
    double qvp = sqrt((0.25/s)*(s-pow(Mv+Mp,2))*(s-pow(Mv-Mp,2)));
    double qgp = (s-Mp*Mp)/(2.0*sqrt(s));
    double coeff = pow(e*f/Mv,2)/(64.0*M_PI*s*qgp*qgp);
    double ReT = RealT(s);
    double ImT = ImaginaryT(s);
    double dodt_t0 = coeff * (ReT*ReT + ImT*ImT) * 3.89e5; // units nb/GeV^2
    double A = dodt_t0;
    double tmin_ = tmin(Mv,Mp,W);
    double tmax_ = tmax(Mv,Mp,W);
   
    double sigma_el = (std::pow(e * f / Mv, 2) * 1 / (2 * W * qgp) *
		       (qvp / qgp) * C_el *
		       std::pow(1 - nu_el / nu, b_el) *
		       std::pow(nu / nu_el, a_el) * khbarc2); // units nb

    _f1->SetParameters(tmin_,tmax_,sigma_el,A);

    double B = _f1->GetX(0.0,0.5,10); // units 1/GeV^2
    // Set the event A and event B
    eventA = A;
    eventB = B;
  }

  void set_b0()
  {
    double W = (Mv+Mp)*(1.05);
    double s = W*W;
    double nu = 0.5*(s-Mp*Mp-Mv*Mv);
    double qvp = sqrt((0.25/s)*(s-pow(Mv+Mp,2))*(s-pow(Mv-Mp,2)));
    double qgp = (s-Mp*Mp)/(2.0*sqrt(s));
    double coeff = pow(e*f/Mv,2)/(64.0*M_PI*s*qgp*qgp);
    double ReT = RealT(s);
    double ImT = ImaginaryT(s);
    double dodt_t0 = coeff * (ReT*ReT + ImT*ImT) * 3.89e5; // units nb/GeV^2
    double A = dodt_t0;
    double tmin_ = tmin(Mv,Mp,W);
    double tmax_ = tmax(Mv,Mp,W);
   
    double sigma_el = (std::pow(e * f / Mv, 2) * 1 / (2 * W * qgp) *
		       (qvp / qgp) * C_el *
		       std::pow(1 - nu_el / nu, b_el) *
		       std::pow(nu / nu_el, a_el) * khbarc2); // units nb
    _f1->SetParameters(tmin_,tmax_,sigma_el,A);
    double b0 = _f1->GetX(0.0,0.5,10); // units 1/GeV^2

    // Set the event b0
    eventb0 = b0;
  }


  double dSigmaY1S_v1(const double W, const double t){
    set_B(Mv,Mp,W);
    return eventA * exp(eventB*t); // units nb/GeV^2
  }
  double dSigmaY1S_v2(const double W, const double t){
    double jac = exp(-eventb0*t)/eventb0;
    return eventA*exp(eventB*t) * jac; // dsigma_dexp_b0t 
  }
  int SetModel(const char * model = "v1"){
    if (strcmp(model, "v1") == 0)
      {
	vflag = 1;
	dSigmaY1S = &dSigmaY1S_v1;
	_f1 = new TF1("B_func",B,0,10,4);
	_F.function = &dispersion_integral;
      }
    else if (strcmp(model, "v2") == 0)
      {
	vflag = 2;
	dSigmaY1S = &dSigmaY1S_v2;
	set_b0();
	_f1 = new TF1("B_func",B,0,10,4);
	_F.function = &dispersion_integral;
      }
    else {
      cout << "No matching model! Set to v1 model!" << endl;
      dSigmaY1S = &dSigmaY1S_v1;
      vflag = 1;
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
  
  double VirtualPhoton(TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', gamma
    double m = PARTICLE::e.M();
    double mp = ki[1].M();
    double mY = PARTICLE::upsilon1S.M();

    // Step 1.) Select event Q2, y
    double Q2 = random.Uniform(Q2range[0],Q2range[1]);
    //  double W2 = random.Uniform(pow(Wrange[0],2),pow(Wrange[1],2));
    double y = random.Uniform(yrange[0],yrange[1]);
    // Step 2.) Boost both e & N into "N" rest frame where the kinematics are easier
    TVector3 beta = ki[1].BoostVector();
    ki[0].Boost(-beta);
    ki[1].Boost(-beta);
    // Step 2.5) Rotate the vectors such that the z direction points along ki[0]'s momentum
    TVector3 direction(1.0,0,0);
    direction.SetPhi(ki[0].Phi());
    direction.SetTheta(ki[0].Theta()/2.0);
    ki[0].Rotate(M_PI,direction);
    // Step 3.) Calculate the energy, momentum, etc. of the scattered e- and gamma*
    double _Ee = ki[0].E();
    double W2 = (2 * ki[1].M() * _Ee * y - Q2 + ki[1].M()*ki[1].M());
    double _Eg = (W2 - mp * mp + Q2) / ( 2.0 * mp);
    double _Eeprime = _Ee - _Eg;
    if(_Eeprime < 0)
      {
	return 0; // impossible event
      }
    if ((W2 < mp * mp) || (sqrt(W2) < Wrange[0]) || (sqrt(W2) > Wrange[1]))
      {
	return 0;//below the lowest state
      }
    double _Pe = sqrt(_Ee*_Ee - m * m);
    double _Peprime = sqrt(_Eeprime*_Eeprime - m * m);
    double _th = M_PI-std::acos((-Q2 - 2 * m * m + 2 * _Ee * _Eeprime) / (2 * _Pe * _Peprime));
    if(isnan(_th))
      {
	return 0;
      }
    double _cth= cos(_th);
    double _sth = sqrt(1.0 - _cth * _cth);
    double _phi = random.Uniform(-M_PI, M_PI);
    kf[0].SetXYZM(_Peprime * _sth * cos(_phi) , _Peprime * _sth * sin(_phi) , -_Peprime * _cth, m);//e' in "N" rest frame
    kf[1] = ki[0]-kf[0]; // virtual photon in "N" rest frame
    // Step 3.5) Unrotate the vectors
    ki[0].Rotate(M_PI,direction);
    kf[0].Rotate(M_PI,direction);
    kf[1].Rotate(M_PI,direction);
    // Step 4.) Boost back into originial frame
    ki[0].Boost(beta);
    ki[1].Boost(beta);
    kf[0].Boost(beta);
    kf[1].Boost(beta);
    double alpha_em = 1.0 / 137.0;
    //    double volume = 2.0 * M_PI * abs(Q2range[1] - Q2range[0]) * abs(yrange[1] - yrange[0]);
    double gy = ki[1].M() * sqrt(Q2) / (ki[1] * ki[0]);
    double epsilon = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    double dipole = pow((mY*mY)/(Q2+mY*mY),2.575); // Equation A4
    double R = pow((2.164*mY*mY + Q2)/(2.164*mY*mY),2.131) - 1.0;
    double gammaT = alpha_em/(2*M_PI)*(1.0+pow(1.0-y,2))/(y*Q2);
    return (1.0 + epsilon * R) * dipole * gammaT;
  }

  /* Upsilon1S productions */
  double Upsilon1SElectroproduction(TLorentzVector * ki, TLorentzVector *kf){
    //ki: e, N; kf: e', Psi2S, N'
    double weight1 = VirtualPhoton(ki, kf);//Generate scattered electron
    if (weight1 == 0) return 0;
    double mp = ki[1].M();
    TLorentzVector Pout = kf[1] + ki[1]; // q + N
    double W = Pout.M();
    double W2 = W * W;
    double Q2 = -(kf[1]*kf[1]);
    double Mup = PARTICLE::upsilon1S.RandomM();
    if (W < Mup + Mp)
      {
	return 0; //below the threshold
      }
    // Step 1.) Set & Get the b0 & B for the event
    UPSILONMODEL::set_B(Mup,mp,W);
    double B = UPSILONMODEL::eventB;
    double b0 = UPSILONMODEL::eventb0;
    // Step 2.) Generate a t value within the range from a distribution
    // v1 --> Uniform
    // v2 --> Expontential
    double t = 0.0;
    if(UPSILONMODEL::vflag == 1)
      t = random.Uniform(trange[0],trange[1]);
    else if(UPSILONMODEL::vflag == 2)
      t = (1.0/b0)*std::log(random.Uniform(exp(b0*trange[0]),exp(b0*trange[1])));
    // Step 3.) Calculate beta s.t. we boost into the q + N rest C.O.M frame
    TVector3 beta = Pout.BoostVector();
    // Step 3.5) Calculate axis of rotation such that final state particles are created in a frame where "p" lies along the z-axis
    ki[1].Boost(-beta);
    kf[1].Boost(-beta);
    TVector3 direction(1.0,0,0);
    direction.SetPhi(ki[1].Phi());
    direction.SetTheta(ki[1].Theta()/2.0);
    kf[1].Boost(beta);
    ki[1].Boost(beta);
    // Step 4.) Calculate the final state particles in the C.O.M frame
    // From lAger's code
    // t --> "target"
    // r --> "recoil"
    // v --> "vector meson"
    const double Et_cm = (W2 + Q2 + mp*mp) / (2. * W);
    const double Pt_cm = sqrt(Et_cm * Et_cm - mp*mp);
    const double Er_cm = (W2 - Mup*Mup + Mp*Mp) / (2. * W);
    const double Pr_cm = sqrt(Er_cm * Er_cm - Mp*Mp);
    const double Ev_cm = (W2 + Mup*Mup - Mp*Mp) / (2. * W);
    const double Pv_cm = sqrt(Ev_cm * Ev_cm - Mup*Mup);
    // Step 5.) From the generated "t", get the "theta" of the scattered p in this frame
    const double ctheta_cm =
      (t + 2 * Et_cm * Er_cm - mp * mp - Mp * Mp) / (2 * Pt_cm * Pr_cm);
    if(ctheta_cm>1.0||ctheta_cm<-1.0)
      {
	return 0;
      }
    const double theta_cm = std::acos(ctheta_cm);
    const double phi_cm = random.Uniform(0, TMath::TwoPi());
    const double stheta_cm = sqrt(1.0 - ctheta_cm * ctheta_cm);
    const double theta_cm2 = TMath::Pi() - theta_cm;
    const double ctheta_cm2 = std::cos(theta_cm2);
    const double stheta_cm2 = std::sin(theta_cm2);
    // Step 6.) Set the VM and p' in this C.O.M frame
    kf[1].SetXYZM(Pv_cm * stheta_cm2 * cos(phi_cm), Pv_cm * stheta_cm2 * sin(phi_cm), Pv_cm * ctheta_cm2, Mup);
    kf[2].SetXYZM(-Pr_cm * stheta_cm * cos(phi_cm), Pr_cm * stheta_cm * sin(-phi_cm), Pr_cm * ctheta_cm, Mp);
    // Step 6.5) Rotate the VM and p' such that p is on the z-axis
    kf[1].Rotate(M_PI,direction);
    kf[2].Rotate(M_PI,direction);
    // Step 7.) Boost these particles back into the original frame
    kf[1].Boost(beta);
    kf[2].Boost(beta);
    //    double volume = exp(b0*trange[1])-exp(b0*trange[0]);
    // Need to include flux factor because of moving target in this frame? //
    double weight2 = UPSILONMODEL::dSigmaY1S(W,t);
    return weight1*weight2;
  }


  double Event_eN2eNee_Upsilon1S(TLorentzVector * ki, TLorentzVector * kf){
    //ki: e, N; kf: e', N', [e+, e-]
    TLorentzVector kf1[3];//e', Y1S, N'
    double weight = Upsilon1SElectroproduction(ki, kf1);
    kf[0] = kf1[0];//e'
    kf[1] = kf1[2];//N'
    double mass[2] = {PARTICLE::e.M(), PARTICLE::e.M()};
    GenPhase.SetDecay(kf1[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//e+
    kf[3] = *GenPhase.GetDecay(1);//e-
    double Mup = kf1[1].M();
    double Ep = kf1[2] * kf1[1] / Mup;//recoil proton energy in Upsilon1S rest frame
    double p = sqrt(Ep * Ep - Mp * Mp);//recoil proton momentum in Upsilon1S rest frame
    double l = sqrt(pow(Mup * Mup - kf[2] * kf[2] - kf[3] * kf[3], 2) - 4.0 * (kf[2] * kf[2]) * (kf[3] * kf[3])) / (2.0 * Mup);//decayed lepton momentum in Upsilon1S rest frame
    double cth = (Ep * Mup / 2.0 - kf[2] * kf[1]) / (p * l);//cos(theta) between final lepton and final proton in Upsilon1S rest frame
    double y = (ki[0].E() - kf[0].E()) / ki[0].E();
    double Q2 = - (ki[0] - kf[0]) * (ki[0] - kf[0]);
    double gy = sqrt(Q2) / ki[0].E();
    double epsilon = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy);
    double R = pow(1.0 + Q2 / 2.164 / pow(Mup,2), 2.131) - 1.0;
    double r = epsilon * R / (1.0 + epsilon * R);
    double wth = 3.0 / 4.0 * (1.0 + r + (1.0 - 3.0 * r) * pow(cth,2));
    double branch = 2.38e-2;//Branch ratio to e+e- (Y1S)
    return weight * wth * branch;
  }
}
#endif
