int test()
{
  double Mv = 3.686;
  double Mp = 0.938272;
  double t = -0.15;
  double Px = gRandom->Rndm();
  double Py = gRandom->Rndm();
  double Pz = gRandom->Rndm();

  TLorentzVector pIn;
  pIn.SetPxPyPzE(Px,Py,Pz,sqrt(Px*Px+Py*Py+Pz*Pz+Mp*Mp));
  
  TLorentzVector pOut;
  double cth = (t-2*Mp*Mp+2*
  return 0;
}
