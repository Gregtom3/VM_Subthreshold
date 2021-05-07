#ifndef _LPHASE_H_
#define _LPHASE_H_

#include <cmath>
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"

namespace LPHASE{
  double VolumePHS2(const double E, const double mass[2]);
  double VolumePHS3(const double E, const double mass[3]);
  double E1E2Range3(const double * x, const double * par);
  double VPHS2(const double E, const double mass[2]);
  double dVPHS3(const double * x, const double * par);
  double VPHS3(const double E, const double mass[3]);
  double dVPHS4(const double * x, const double * par);
  double VPHS4(const double E, const double mass[4]);
  double dVPHS5(const double * x, const double * par);
  double VPHS5(const double E, const double mass[5]);
  double dVPHS6(const double * x, const double * par);
  double VPHS6(const double E, const double mass[6]);
  double dVPHS7(const double * x, const double * par);
  double VPHS7(const double E, const double mass[6]);
  double VPHS(const double E, const double * mass, const int Nf);
}

double LPHASE::VolumePHS2(const double E, const double mass[2]){
  if (E <= mass[0] + mass[1])
    return 0;
  double p = sqrt( (E * E - pow(mass[0] + mass[1], 2)) * (E * E - pow(mass[0] - mass[1], 2))) / (2.0 * E);
  return M_PI * p / E;
}

double LPHASE::E1E2Range3(const double * x, const double * par){
  double E1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  if (E - E1 < m2 + m3 || E1 < m1)
    return 0;
  double M23 = sqrt(E * E - 2.0 * E * E1 + m1 * m1);
  double k = sqrt( (M23 * M23 - pow(m2 + m3, 2)) * (M23 * M23 - pow(m2 - m3, 2))) / (2.0 * M23);
  double gamma = (E - E1) / M23;
  double beta = sqrt(E1 * E1 - m1 * m1) / (E - E1);
  double E2max = gamma * (sqrt(m2 * m2 + k * k) + beta * k);
  double E2min = gamma * (sqrt(m2 * m2 + k * k) - beta * k);
  return E2max - E2min;
}

double LPHASE::VolumePHS3(const double E, const double mass[3]){
  if (E <= mass[0] + mass[1] + mass[2])
    return 0;
  double par[4] = {E, mass[0], mass[1], mass[2]};
  double E1min = mass[0];
  double E1max = (E * E + mass[0] * mass[0] - pow(mass[1] + mass[2], 2)) / (2.0 * E);
  TF1 f0("E1 integrand", LPHASE::E1E2Range3, E1min, E1max, 4);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(E1min, E1max);
  return M_PI * M_PI * result;
}

//recursive

double LPHASE::VPHS2(const double E, const double mass[2]){
  double m1 = mass[0];
  double m2 = mass[1];
  if (E <= m1 + m2)
    return 0;
  double p = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  return M_PI * p / E;
}

double LPHASE::dVPHS3(const double * x, const double * par){
  double p1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  if (E <= sqrt(p1 * p1 + m1 * m1) + sqrt(p1 * p1 + pow(m2 + m3, 2)))
    return 0;
  double M = sqrt(E * E - 2.0 * E * sqrt(p1 * p1 + m1 * m1) + m1 * m1);
  double mass[2] = {m2, m3};
  double result = 2.0 * M_PI * p1 * p1 / sqrt(p1 * p1 + m1 * m1) * LPHASE::VPHS2(M, mass);
  return result;
}

double LPHASE::VPHS3(const double E, const double mass[3]){
  double m1 = mass[0];
  double m2 = mass[1] + mass[2];
  if (E <= m1 + m2)
    return 0;
  double par[4] = {E, mass[0], mass[1], mass[2]};
  double pmax = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  TF1 f0("p1 integrand", LPHASE::dVPHS3, 0, pmax, 4);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0, pmax);
  return result;
}

double LPHASE::dVPHS4(const double * x, const double * par){
  double p1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  double m4 = par[4];
  if (E <= sqrt(p1 * p1 + m1 * m1) + sqrt(p1 * p1 + pow(m2 + m3 + m4, 2)))
    return 0;
  double M = sqrt(E * E - 2.0 * E * sqrt(p1 * p1 + m1 * m1) + m1 * m1);
  double mass[3] = {m2, m3, m4};
  double result = 2.0 * M_PI * p1 * p1 / sqrt(p1 * p1 + m1 * m1) * LPHASE::VPHS3(M, mass);
  return result;
}

double LPHASE::VPHS4(const double E, const double mass[4]){
  double m1 = mass[0];
  double m2 = mass[1] + mass[2] + mass[3];
  if (E <= m1 + m2)
    return 0;
  double par[5] = {E, mass[0], mass[1], mass[2], mass[3]};
  double pmax = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  TF1 f0("p1 integrand", LPHASE::dVPHS4, 0, pmax, 5);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0, pmax);
  return result;
}

double LPHASE::dVPHS5(const double * x, const double * par){
  double p1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  double m4 = par[4];
  double m5 = par[5];
  if (E <= sqrt(p1 * p1 + m1 * m1) + sqrt(p1 * p1 + pow(m2 + m3 + m4 + m5, 2)))
    return 0;
  double M = sqrt(E * E - 2.0 * E * sqrt(p1 * p1 + m1 * m1) + m1 * m1);
  double mass[4] = {m2, m3, m4, m5};
  double result = 2.0 * M_PI * p1 * p1 / sqrt(p1 * p1 + m1 * m1) * LPHASE::VPHS4(M, mass);
  return result;
}

double LPHASE::VPHS5(const double E, const double mass[5]){
  double m1 = mass[0];
  double m2 = mass[1] + mass[2] + mass[3] + mass[4];
  if (E <= m1 + m2)
    return 0;
  double par[6] = {E, mass[0], mass[1], mass[2], mass[3], mass[4]};
  double pmax = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  TF1 f0("p1 integrand", LPHASE::dVPHS5, 0, pmax, 6);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0, pmax);
  return result;
}

double LPHASE::dVPHS6(const double * x, const double * par){
  double p1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  double m4 = par[4];
  double m5 = par[5];
  double m6 = par[6];
  if (E <= sqrt(p1 * p1 + m1 * m1) + sqrt(p1 * p1 + pow(m2 + m3 + m4 + m5 + m6, 2)))
    return 0;
  double M = sqrt(E * E - 2.0 * E * sqrt(p1 * p1 + m1 * m1) + m1 * m1);
  double mass[5] = {m2, m3, m4, m5, m6};
  double result = 2.0 * M_PI * p1 * p1 / sqrt(p1 * p1 + m1 * m1) * LPHASE::VPHS5(M, mass);
  return result;
}

double LPHASE::VPHS6(const double E, const double mass[5]){
  double m1 = mass[0];
  double m2 = mass[1] + mass[2] + mass[3] + mass[4] + mass[5];
  if (E <= m1 + m2)
    return 0;
  double par[7] = {E, mass[0], mass[1], mass[2], mass[3], mass[4], mass[5]};
  double pmax = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  TF1 f0("p1 integrand", LPHASE::dVPHS6, 0, pmax, 7);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0, pmax);
  return result;
}

double LPHASE::dVPHS7(const double * x, const double * par){
  double p1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  double m4 = par[4];
  double m5 = par[5];
  double m6 = par[6];
  double m7 = par[7];
  if (E <= sqrt(p1 * p1 + m1 * m1) + sqrt(p1 * p1 + pow(m2 + m3 + m4 + m5 + m6 + m7, 2)))
    return 0;
  double M = sqrt(E * E - 2.0 * E * sqrt(p1 * p1 + m1 * m1) + m1 * m1);
  double mass[6] = {m2, m3, m4, m5, m6, m7};
  double result = 2.0 * M_PI * p1 * p1 / sqrt(p1 * p1 + m1 * m1) * LPHASE::VPHS6(M, mass);
  return result;
}

double LPHASE::VPHS7(const double E, const double mass[5]){
  double m1 = mass[0];
  double m2 = mass[1] + mass[2] + mass[3] + mass[4] + mass[5] + mass[6];
  if (E <= m1 + m2)
    return 0;
  double par[8] = {E, mass[0], mass[1], mass[2], mass[3], mass[4], mass[5], mass[6]};
  double pmax = sqrt((E * E - pow(m1 + m2, 2)) * (E * E - pow(m1 - m2, 2))) / (2.0 * E);
  TF1 f0("p1 integrand", LPHASE::dVPHS6, 0, pmax, 8);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0, pmax);
  return result;
}

double LPHASE::VPHS(const double E, const double * mass, const int Nf){
  if (Nf == 2)
    return LPHASE::VPHS2(E, mass);
  if (Nf == 3)
    return LPHASE::VPHS3(E, mass);
  if (Nf == 4)
    return LPHASE::VPHS4(E, mass);
  if (Nf == 5)
    return LPHASE::VPHS5(E, mass);
  if (Nf == 6)
    return LPHASE::VPHS6(E, mass);
  if (Nf == 7)
    return LPHASE::VPHS7(E, mass);

  return 0;
}


#endif
