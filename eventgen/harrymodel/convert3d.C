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
  double zero = 1e-200;

  if (opt == 1){//jpsi from d, 3-fold
    double q, theta_k, k, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-3D-default.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    ifstream infile("duke-6-7-7p2/jpsi-7p2.out-default");
    hs = new TH3D("ds_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-7p7/jpsi-7p7.out-default");
    hs = new TH3D("ds_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-8p2/jpsi-8p2.out-default");
    hs = new TH3D("ds_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();
    

    fs->Close();
  }

  if (opt == 2){//jpsi from d, 3-fold
    double q, theta_k, k, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-3D-nv2Ia.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    ifstream infile("duke-6-7-7p2/jpsi-7p2.out-nv2Ia");
    hs = new TH3D("ds_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-7p7/jpsi-7p7.out-nv2Ia");
    hs = new TH3D("ds_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-8p2/jpsi-8p2.out-nv2Ia");
    hs = new TH3D("ds_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();
    

    fs->Close();
  }

  if (opt == 3){//jpsi from d, 3-fold
    double q, theta_k, k, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-3D-v18.root", "RECREATE");
    TH3D * hs, * hp;
    double vol;
    ifstream infile("duke-6-7-7p2/jpsi-7p2.out-v18");
    hs = new TH3D("ds_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-7p7/jpsi-7p7.out-v18");
    hs = new TH3D("ds_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E7.7", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();

    infile.open("duke-6-7-8p2/jpsi-8p2.out-v18");
    hs = new TH3D("ds_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    hp = new TH3D("p_E8.2", "", 900, 1.0, 10.0, 6, 0.0, 12.0, 30, 5.5, 35.5);
    infile.getline(dum, 200);
    while (infile >> q >> theta_k >> k >> theta_p >> p >> crst >> tmp >> tmp){
      vol = 0.01 * (cos(theta_k * deg) - cos((theta_k + 2.0) * deg)) * (2.0 * M_PI) * (cos((theta_p - 0.5) * deg) - cos((theta_p + 0.5) * deg));
      hs->Fill(k, theta_k + eps, theta_p, crst * vol + zero);
      hp->Fill(k, theta_k + eps, theta_p, p);
    }
    infile.close();
    cout << hs->Integral() << endl;
    fs->Write();
    hs->Delete();
    hp->Delete();
    

    fs->Close();
  }

  return 0;
}
    
    
