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
  double Ebeam = 11.0;//GeV
  double lumi = 1.2e37 * 0.5 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s

  // Set nuclear
  NUCLEAR::SetNuclear("D");
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
  GENERATE::Set_psi_2S();
  // Set Jpsi production model
  JPSIMODEL::SetModel("2S_23g");
  GENERATE::SetBremsstrahlung();
  double kmin = 9.5;
  double kmax = Ebeam;
  
  
  // Set detector
  DETECTOR::SetDetector("SoLID");
 
  // TTree Info
  TFile * fall2 = new TFile("result-photo-psi2S/Dsolid_photo.root", "RECREATE");
  TTree * tree = new TTree("tree","");
  tree->SetDirectory(fall2);

  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  double weight_smear = 0.0;
  double weight_smear2 = 0.0; //will include lumi, time, nsim
  //double acceptance = 0.0;
  double Mjpsi = 3.686097;
  int count = 0;
  int is_sub;

  double Eg=0.0;
  double Eg_smear=0.0;
  
  TLorentzVector *eIn = new TLorentzVector();
  TLorentzVector *pIn = new TLorentzVector();
  TLorentzVector *pOut = new TLorentzVector();
  TLorentzVector *ePlusOut = new TLorentzVector();
  TLorentzVector *eMinusOut = new TLorentzVector();
  TLorentzVector _pOutSmear;
  TLorentzVector _ePlusOutSmear;
  TLorentzVector _eMinusOutSmear;
  TLorentzVector *pOutSmear = new TLorentzVector();
  TLorentzVector *ePlusOutSmear = new TLorentzVector();
  TLorentzVector *eMinusOutSmear = new TLorentzVector();
  TLorentzVector vm;
  TLorentzVector vmSmear;
  
  // Add TTree Branches
  tree->Branch("eIn","TLorentzVector",&eIn);
  tree->Branch("pIn","TLorentzVector",&pIn);
  tree->Branch("pOut","TLorentzVector",&pOut);
  tree->Branch("ePlusOut","TLorentzVector",&ePlusOut);
  tree->Branch("eMinusOut","TLorentzVector",&eMinusOut);
  tree->Branch("pOutSmear","TLorentzVector",&pOutSmear);
  tree->Branch("ePlusOutSmear","TLorentzVector",&ePlusOutSmear);
  tree->Branch("eMinusOutSmear","TLorentzVector",&eMinusOutSmear);
  tree->Branch("vm","TLorentzVector",&vm);
  tree->Branch("vmSmear","TLorentzVector",&vmSmear);
  tree->Branch("weight",&weight,"Double/D");
  tree->Branch("weight_smear",&weight_smear,"Double/D");
  tree->Branch("weight_smear2",&weight_smear2,"Double/D");
  tree->Branch("Eg",&Eg,"Double/D");
  tree->Branch("Eg_smear",&Eg_smear,"Double/D");
  tree->Branch("is_sub",&is_sub,"Sub/I");
  tree->Branch("Nsim",&Nsim,"Sims/I");
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/10) == 0) cout << i/(Nsim/10)*10 << "%" << endl;

    weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam) * 1.95 / 2; // 15cm LD2 target
    weight *= GENERATE::GetNucleon(&ki[1]);
    weight *= GENERATE::Event_gN2Nee_Jpsi(ki, kf); 

    if (weight > 0.0){
      count++;
      //case 1
      eIn=(TLorentzVector*)ki[0].Clone();
      pIn=(TLorentzVector*)ki[1].Clone();
      pOut=(TLorentzVector*)kf[0].Clone();
      ePlusOut=(TLorentzVector*)kf[1].Clone();
      eMinusOut=(TLorentzVector*)kf[2].Clone();
      _pOutSmear=kf[0];
      _ePlusOutSmear=kf[1];
      _eMinusOutSmear=kf[2];

      pOutSmear=(TLorentzVector*)_pOutSmear.Clone();
      ePlusOutSmear=(TLorentzVector*)_ePlusOutSmear.Clone();
      eMinusOutSmear=(TLorentzVector*)_eMinusOutSmear.Clone();
     
      // Filling in TTree

      // Calculate if the process was "subthreshold"
      is_sub = 0;
      if (ki[0].E() < 10.9267)
	{
	  is_sub = 1;
	}
      
     
      // Calculate unsmeared quantities
      Eg = CalcEg(kf[0]+kf[1]+kf[2]);
      vm=kf[1]+kf[2];

      // Calculate smeared quantities

      
      weight_smear = weight*DETECTOR::SmearSoLID(_pOutSmear, "p")*DETECTOR::SmearSoLID(_ePlusOutSmear, "e+")*DETECTOR::SmearSoLID(_eMinusOutSmear, "e-");
      weight_smear2 = weight_smear*lumi*time/Nsim;
      Eg_smear = CalcEg(_ePlusOutSmear+_eMinusOutSmear+_pOutSmear);
      vmSmear=_ePlusOutSmear+_eMinusOutSmear;

      if(weight_smear2>0.0)
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
