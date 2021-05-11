#include "Lcore.h"
#include "TTree.h"
int CheckAcceptance(const TLorentzVector P, const double thmin, const double thmax){
  if (P.Theta() < thmin * M_PI / 180.0) return 0;
  if (P.Theta() > thmax * M_PI / 180.0) return 0;
  if (P.P() < 0.3) return 0;
  return 1;
}
double CalcEg(const TLorentzVector P){
  double Md = 1.8756;
  double Mn = 0.93957;
  double Eg = (pow(P.E() - Md, 2) - pow(P.P(), 2) - Mn * Mn) / (2.0 * (P.E() - P.Pz() - Md));
  return Eg;
}
int main(const int argc, const char * argv[]){

  // Set simulation
  gRandom->SetSeed(0);
  Long64_t Nsim = 10000;

  if (argc > 1) Nsim = atoi(argv[1]);
  else {
    cout << "./electro-solid-study <Nsim>" << endl;
    return 0;
  }

  // Electron beam energy and luminosity
  double Ebeam = 8.5;//GeV
  double lumi = 1.2e37 * 0.5 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Set nuclear
  NUCLEAR::SetNuclear("D");
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
  //  GENERATE::Set_psi_2S();
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set scattered electron range
  double degtorad = M_PI / 180.0;
  GENERATE::cthrange[0] = -1;
  GENERATE::cthrange[1] = 1;
  GENERATE::perange[0] = 0.0;//GeV
  GENERATE::perange[1] = 10.0;//GeV

  // Set detector
  DETECTOR::SetDetector("SoLID");
  
  // detected
  TFile * fall1 = new TFile("fast-result-electro-jpsi/Dsolid.root", "RECREATE");
  TH1F *h_total = new TH1F("h_total","J/#psi electroproduction at SoLID;M[e^{+}e^{-}];Events / hour",40,2.9,3.3);
  TH1F *h_above = new TH1F("h_above","h_above",40,2.9,3.3);
  TH1F *h_below = new TH1F("h_below","h_below",40,2.9,3.3);
  h_total->SetLineColor(kBlack);
  h_above->SetLineColor(kRed);
  h_below->SetLineColor(kBlue);

  
  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  double weight_smear = 0.0;
  double weight_smear2 = 0.0; //will include lumi, time, nsim
  //double acceptance = 0.0;
  double Mjpsi = 3.097;
  int count = 0;
  int is_sub;

  double Eg=0.0;
  double Eg_smear=0.0;
  
  TLorentzVector _eOutSmear;
  TLorentzVector _pOutSmear;
  TLorentzVector _ePlusOutSmear;
  TLorentzVector _eMinusOutSmear;
 
  TLorentzVector vm;
  TLorentzVector vmSmear;
  
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;
    
    weight = GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight > 0.0){
      count++;
      q = ki[0] - kf[0];
 
      _eOutSmear=kf[0];
      _pOutSmear=kf[1];
      _ePlusOutSmear=kf[2];
      _eMinusOutSmear=kf[3];
      
      // Calculate if the process was "subthreshold"
      is_sub = 0;
      if (q.E() < Mjpsi + (Mjpsi * Mjpsi - q * q) / (2.0 * Mp)){
	is_sub = 1;
      }
      
      // Calculate unsmeared quantities
      Eg = CalcEg(kf[2]+kf[3]+kf[1]);
      vm=kf[2]+kf[3];

      // Calculate smeared quantities
      DETECTOR::SmearSoLID(_eOutSmear, "e-"); // smear final state electron, but its detection is unimportant for event reco
      weight_smear = weight*DETECTOR::SmearSoLID(_pOutSmear, "p")*DETECTOR::SmearSoLID(_ePlusOutSmear, "e+")*DETECTOR::SmearSoLID(_eMinusOutSmear, "e-");
      weight_smear2 = weight_smear*lumi*time/Nsim;
      Eg_smear = CalcEg(_ePlusOutSmear+_eMinusOutSmear+_pOutSmear);
      vmSmear=_ePlusOutSmear+_eMinusOutSmear;


      // Fill in histograms
      h_total->Fill(vmSmear.M(),weight_smear2);
      h_above->Fill(vmSmear.M(),weight_smear2*(is_sub==0));
      h_below->Fill(vmSmear.M(),weight_smear2*(is_sub==1));
      
    }
  }
  cout << count << endl;
  GENERATE::print_fail();
  //1
  
  fall1->Write();
  fall1->Close();



  return 0;
}
