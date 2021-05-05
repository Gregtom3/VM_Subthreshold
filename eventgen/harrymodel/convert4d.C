#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH3D.h"
#include "TString.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./convert4d <opt>" << endl;
    return 0;
  }

  const int opt = atoi(argv[1]);

  double eps = 1e-10;
  double zero = 1e-250;

  if (opt == 0){//jpsi from d, 4-fold
    double q, theta_k, k, phi_p, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-4D-v18.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    //double vol = 0.01 * (2.0 * deg) * (2.0 * M_PI) * (1.0 * deg) * (2.0 * M_PI);
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p2-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p3-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.3_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.3_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p4-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.4_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.4_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p5-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.5_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.5_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p6-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.6_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.6_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p7-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p8-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.8_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.8_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p9-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.9_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.9_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p0-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.0_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.0_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p1-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.1_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.1_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p2-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p3-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.3_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.3_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p4-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.4_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.4_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p5-phip-v18/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.5_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.5_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    fs->Close();
  }

  if (opt == 1){//jpsi from d, 4-fold
    double q, theta_k, k, phi_p, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-4D-default.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    //double vol = 0.01 * (2.0 * deg) * (2.0 * M_PI) * (1.0 * deg) * (2.0 * M_PI);
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p2-phip/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i <<": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p7-phip/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p2-phip/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    fs->Close();
  }

  if (opt == 2){//jpsi from d, 4-fold
    double q, theta_k, k, phi_p, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-4D-nv2Ia.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    //double vol = 0.01 * (2.0 * deg) * (2.0 * M_PI) * (1.0 * deg) * (2.0 * M_PI);
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p2-phip-nv2IIa/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p7-phip-nv2IIa/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E7.7_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p2-phip-nv2IIa/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      hp = new TH3D(Form("p_E8.2_idx%d", i), "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg)) * (2.0 * M_PI);
	hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      fs->Write();
      cout << i << ": " << hs->Integral() << endl;
      hs->Delete();
      hp->Delete();
    }
    fs->Close();
  }

  return 0;
}
    
    
