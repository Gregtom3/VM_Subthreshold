#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2) {
    cout << "./convert <opt>" << endl;
    return 0;
  }

  const int opt = atoi(argv[1]);

  double eps=1e-10;
  double zero = 1e-200;
  
  if (opt == 1){//jpsi
    double q, k, theta, ds, pmin, pmax, tmp;
    char dum[300];
    
    TFile * fs = new TFile("harrymodel-jpsi.root", "RECREATE");
    
    TH2D * hE62 = new TH2D("ds_E=6.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE62->SetDirectory(fs);
    TH2D * hE65 = new TH2D("ds_E=6.5","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE65->SetDirectory(fs);
    TH2D * hE68 = new TH2D("ds_E=6.8","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE68->SetDirectory(fs);
    TH2D * hE71 = new TH2D("ds_E=7.1","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE71->SetDirectory(fs);
    TH2D * hE74 = new TH2D("ds_E=7.4","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE74->SetDirectory(fs);
    TH2D * hE77 = new TH2D("ds_E=7.7","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE77->SetDirectory(fs);
    TH2D * hE80 = new TH2D("ds_E=8.0","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE80->SetDirectory(fs);
    TH2D * hE83 = new TH2D("ds_E=8.3","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE83->SetDirectory(fs);
    TH2D * hE86 = new TH2D("ds_E=8.6","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE86->SetDirectory(fs);
    TH2D * hE89 = new TH2D("ds_E=8.9","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE89->SetDirectory(fs);
    TH2D * hE92 = new TH2D("ds_E=9.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    hE92->SetDirectory(fs);

    TH3D * hpmin = new TH3D("pmin", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    hpmin->SetDirectory(fs);
    TH3D * hpmax = new TH3D("pmax", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    hpmax->SetDirectory(fs);
    
    ifstream fhE62("file-jpsi-pmx-6.2.dat");
    fhE62.getline(dum, 300);
    while(fhE62 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE62->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE62.close();

    ifstream fhE65("file-jpsi-pmx-6.5.dat");
    fhE65.getline(dum, 300);
    while(fhE65 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE65->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE65.close();

    ifstream fhE68("file-jpsi-pmx-6.8.dat");
    fhE68.getline(dum, 300);
    while(fhE68 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE68->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE68.close();

    ifstream fhE71("file-jpsi-pmx-7.1.dat");
    fhE71.getline(dum, 300);
    while(fhE71 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE71->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE71.close();

    ifstream fhE74("file-jpsi-pmx-7.4.dat");
    fhE74.getline(dum, 300);
    while(fhE74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE74.close();

    ifstream fhE77("file-jpsi-pmx-7.7.dat");
    fhE77.getline(dum, 300);
    while(fhE77 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE77->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE77.close();

    ifstream fhE80("file-jpsi-pmx-8.0.dat");
    fhE80.getline(dum, 300);
    while(fhE80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE80.close();

    ifstream fhE83("file-jpsi-pmx-8.3.dat");
    fhE83.getline(dum, 300);
    while(fhE83 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE83->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE83.close();
    
    ifstream fhE86("file-jpsi-pmx-8.6.dat");
    fhE86.getline(dum, 300);
    while(fhE86 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE86->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE86.close();

    ifstream fhE89("file-jpsi-pmx-8.9.dat");
    fhE89.getline(dum, 300);
    while(fhE89 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE89->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE89.close();

    ifstream fhE92("file-jpsi-pmx-9.2.dat");
    fhE92.getline(dum, 300);
    while(fhE92 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      hE92->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      hpmin->Fill(q, k, theta+eps, pmin);
      hpmax->Fill(q, k, theta+eps, pmax);
    }
    fhE92.close();

    TH2D * cE62 = new TH2D("ds_E=6.2_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE62->SetDirectory(fs);
    TH2D * cE65 = new TH2D("ds_E=6.5_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE65->SetDirectory(fs);
    TH2D * cE68 = new TH2D("ds_E=6.8_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE68->SetDirectory(fs);
    TH2D * cE71 = new TH2D("ds_E=7.1_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE71->SetDirectory(fs);
    TH2D * cE74 = new TH2D("ds_E=7.4_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE74->SetDirectory(fs);
    TH2D * cE77 = new TH2D("ds_E=7.7_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE77->SetDirectory(fs);
    TH2D * cE80 = new TH2D("ds_E=8.0_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE80->SetDirectory(fs);
    TH2D * cE83 = new TH2D("ds_E=8.3_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE83->SetDirectory(fs);
    TH2D * cE86 = new TH2D("ds_E=8.6_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE86->SetDirectory(fs);
    TH2D * cE89 = new TH2D("ds_E=8.9_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE89->SetDirectory(fs);
    TH2D * cE92 = new TH2D("ds_E=9.2_c300","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    cE92->SetDirectory(fs);

    TH3D * cpmin = new TH3D("pmin_c300", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    cpmin->SetDirectory(fs);
    TH3D * cpmax = new TH3D("pmax_c300", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    cpmax->SetDirectory(fs);
    
    ifstream fcE62("file-jpsi-pmx-c300-6.2.dat");
    fcE62.getline(dum, 300);
    while(fcE62 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE62->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE62.close();

    ifstream fcE65("file-jpsi-pmx-c300-6.5.dat");
    fcE65.getline(dum, 300);
    while(fcE65 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE65->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE65.close();

    ifstream fcE68("file-jpsi-pmx-c300-6.8.dat");
    fcE68.getline(dum, 300);
    while(fcE68 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE68->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE68.close();

    ifstream fcE71("file-jpsi-pmx-c300-7.1.dat");
    fcE71.getline(dum, 300);
    while(fcE71 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE71->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE71.close();

    ifstream fcE74("file-jpsi-pmx-c300-7.4.dat");
    fcE74.getline(dum, 300);
    while(fcE74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE74.close();

    ifstream fcE77("file-jpsi-pmx-c300-7.7.dat");
    fcE77.getline(dum, 300);
    while(fcE77 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE77->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE77.close();

    ifstream fcE80("file-jpsi-pmx-c300-8.0.dat");
    fcE80.getline(dum, 300);
    while(fcE80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE80.close();

    ifstream fcE83("file-jpsi-pmx-c300-8.3.dat");
    fcE83.getline(dum, 300);
    while(fcE83 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE83->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE83.close();
    
    ifstream fcE86("file-jpsi-pmx-c300-8.6.dat");
    fcE86.getline(dum, 300);
    while(fcE86 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE86->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE86.close();

    ifstream fcE89("file-jpsi-pmx-c300-8.9.dat");
    fcE89.getline(dum, 300);
    while(fcE89 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE89->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE89.close();

    ifstream fcE92("file-jpsi-pmx-c300-9.2.dat");
    fcE92.getline(dum, 300);
    while(fcE92 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      cE92->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      cpmin->Fill(q, k, theta+eps, pmin);
      cpmax->Fill(q, k, theta+eps, pmax);
    }
    fcE92.close();

    fs->Write();

  }

  if (opt == 2){//jpsi cut
    double q, k, theta, ds, pmin, pmax, tmp;
    char dum[300];
    
    TFile * fs = new TFile("harrymodel-jpsi-cut.root", "RECREATE");
    
    TH2D * c1E62 = new TH2D("cut1_E=6.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E62->SetDirectory(fs);
    TH2D * c1E64 = new TH2D("cut1_E=6.4","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E64->SetDirectory(fs);
    TH2D * c1E66 = new TH2D("cut1_E=6.6","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E66->SetDirectory(fs);
    TH2D * c1E68 = new TH2D("cut1_E=6.8","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E68->SetDirectory(fs);
    TH2D * c1E70 = new TH2D("cut1_E=7.0","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E70->SetDirectory(fs);
    TH2D * c1E72 = new TH2D("cut1_E=7.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E72->SetDirectory(fs);
    TH2D * c1E74 = new TH2D("cut1_E=7.4","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E74->SetDirectory(fs);
    TH2D * c1E76 = new TH2D("cut1_E=7.6","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E76->SetDirectory(fs);
    TH2D * c1E78 = new TH2D("cut1_E=7.8","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E78->SetDirectory(fs);
    TH2D * c1E80 = new TH2D("cut1_E=8.0","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E80->SetDirectory(fs);
    TH2D * c1E82 = new TH2D("cut1_E=8.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1E82->SetDirectory(fs);

    TH3D * c1pmin = new TH3D("cut1_pmin", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1pmin->SetDirectory(fs);
    TH3D * c1pmax = new TH3D("cut1_pmax", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    c1pmax->SetDirectory(fs);

    ifstream f1E62("cut-1-jpsi-6.2.dat");
    f1E62.getline(dum, 300);
    while(f1E62 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E62->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E62.close();

    ifstream f1E64("cut-1-jpsi-6.4.dat");
    f1E64.getline(dum, 300);
    while(f1E64 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E64->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E64.close();

    ifstream f1E66("cut-1-jpsi-6.6.dat");
    f1E66.getline(dum, 300);
    while(f1E66 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E66->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E66.close();

    ifstream f1E68("cut-1-jpsi-6.8.dat");
    f1E68.getline(dum, 300);
    while(f1E68 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E68->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E68.close();

    ifstream f1E70("cut-1-jpsi-7.0.dat");
    f1E70.getline(dum, 300);
    while(f1E70 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E70->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E70.close();

    ifstream f1E72("cut-1-jpsi-7.2.dat");
    f1E72.getline(dum, 300);
    while(f1E72 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E72->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E72.close();
    
    ifstream f1E74("cut-1-jpsi-7.4.dat");
    f1E74.getline(dum, 300);
    while(f1E74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E74.close();

    ifstream f1E76("cut-1-jpsi-7.6.dat");
    f1E76.getline(dum, 300);
    while(f1E76 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E76->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E76.close();

    ifstream f1E78("cut-1-jpsi-7.8.dat");
    f1E78.getline(dum, 300);
    while(f1E78 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E78->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E78.close();

    ifstream f1E80("cut-1-jpsi-8.0.dat");
    f1E80.getline(dum, 300);
    while(f1E80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E80.close();

    ifstream f1E82("cut-1-jpsi-8.2.dat");
    f1E82.getline(dum, 300);
    while(f1E82 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E82->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c1pmin->Fill(q, k, theta+eps, pmin);
      c1pmax->Fill(q, k, theta+eps, pmax);
    }
    f1E82.close();



    TH2D * c2E62 = new TH2D("cut2_E=6.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E62->SetDirectory(fs);
    TH2D * c2E64 = new TH2D("cut2_E=6.4","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E64->SetDirectory(fs);
    TH2D * c2E66 = new TH2D("cut2_E=6.6","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E66->SetDirectory(fs);
    TH2D * c2E68 = new TH2D("cut2_E=6.8","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E68->SetDirectory(fs);
    TH2D * c2E70 = new TH2D("cut2_E=7.0","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E70->SetDirectory(fs);
    TH2D * c2E72 = new TH2D("cut2_E=7.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E72->SetDirectory(fs);
    TH2D * c2E74 = new TH2D("cut2_E=7.4","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E74->SetDirectory(fs);
    TH2D * c2E76 = new TH2D("cut2_E=7.6","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E76->SetDirectory(fs);
    TH2D * c2E78 = new TH2D("cut2_E=7.8","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E78->SetDirectory(fs);
    TH2D * c2E80 = new TH2D("cut2_E=8.0","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E80->SetDirectory(fs);
    TH2D * c2E82 = new TH2D("cut2_E=8.2","", 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2E82->SetDirectory(fs);

    TH3D * c2pmin = new TH3D("cut2_pmin", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2pmin->SetDirectory(fs);
    TH3D * c2pmax = new TH3D("cut2_pmax", "", 11, 6.05, 9.35, 700, 3.0, 10.0, 11, 0.0, 22.0);
    c2pmax->SetDirectory(fs);

    ifstream f2E62("cut-2-jpsi-6.2.dat");
    f2E62.getline(dum, 300);
    while(f2E62 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E62->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E62.close();

    ifstream f2E64("cut-2-jpsi-6.4.dat");
    f2E64.getline(dum, 300);
    while(f2E64 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E64->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E64.close();

    ifstream f2E66("cut-2-jpsi-6.6.dat");
    f2E66.getline(dum, 300);
    while(f2E66 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E66->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E66.close();

    ifstream f2E68("cut-2-jpsi-6.8.dat");
    f2E68.getline(dum, 300);
    while(f2E68 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E68->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E68.close();

    ifstream f2E70("cut-2-jpsi-7.0.dat");
    f2E70.getline(dum, 300);
    while(f2E70 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E70->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E70.close();

    ifstream f2E72("cut-2-jpsi-7.2.dat");
    f2E72.getline(dum, 300);
    while(f2E72 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E72->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E72.close();
    
    ifstream f2E74("cut-2-jpsi-7.4.dat");
    f2E74.getline(dum, 300);
    while(f2E74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E74.close();

    ifstream f2E76("cut-2-jpsi-7.6.dat");
    f2E76.getline(dum, 300);
    while(f2E76 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E76->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E76.close();

    ifstream f2E78("cut-2-jpsi-7.8.dat");
    f2E78.getline(dum, 300);
    while(f2E78 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E78->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E78.close();

    ifstream f2E80("cut-2-jpsi-8.0.dat");
    f2E80.getline(dum, 300);
    while(f2E80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E80.close();

    ifstream f2E82("cut-2-jpsi-8.2.dat");
    f2E82.getline(dum, 300);
    while(f2E82 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E82->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
      c2pmin->Fill(q, k, theta+eps, pmin);
      c2pmax->Fill(q, k, theta+eps, pmax);
    }
    f2E82.close();


    

    fs->Write();

  }

  if (opt == 3){//jpsi cut deuteron
    double q, k, theta, ds, pmin, pmax, tmp;
    char dum[300];
    
    TFile * fs = new TFile("harrymodel-2h-jpsi-cut.root", "RECREATE");
    
    TH2D * c1E72 = new TH2D("cut1_E=7.2","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E72->SetDirectory(fs);
    TH2D * c1E74 = new TH2D("cut1_E=7.4","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E74->SetDirectory(fs);
    TH2D * c1E76 = new TH2D("cut1_E=7.6","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E76->SetDirectory(fs);
    TH2D * c1E78 = new TH2D("cut1_E=7.8","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E78->SetDirectory(fs);
    TH2D * c1E80 = new TH2D("cut1_E=8.0","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E80->SetDirectory(fs);
    TH2D * c1E82 = new TH2D("cut1_E=8.2","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c1E82->SetDirectory(fs);

    ifstream f1E72("cut-1-2h-jpsi-7.2.dat");
    f1E72.getline(dum, 300);
    while(f1E72 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E72->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E72.close();
    
    ifstream f1E74("cut-1-2h-jpsi-7.4.dat");
    f1E74.getline(dum, 300);
    while(f1E74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E74.close();

    ifstream f1E76("cut-1-2h-jpsi-7.6.dat");
    f1E76.getline(dum, 300);
    while(f1E76 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E76->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E76.close();

    ifstream f1E78("cut-1-2h-jpsi-7.8.dat");
    f1E78.getline(dum, 300);
    while(f1E78 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E78->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E78.close();

    ifstream f1E80("cut-1-2h-jpsi-8.0.dat");
    f1E80.getline(dum, 300);
    while(f1E80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E80.close();

    ifstream f1E82("cut-1-2h-jpsi-8.2.dat");
    f1E82.getline(dum, 300);
    while(f1E82 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c1E82->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f1E82.close();


    TH2D * c2E72 = new TH2D("cut2_E=7.2","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E72->SetDirectory(fs);
    TH2D * c2E74 = new TH2D("cut2_E=7.4","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E74->SetDirectory(fs);
    TH2D * c2E76 = new TH2D("cut2_E=7.6","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E76->SetDirectory(fs);
    TH2D * c2E78 = new TH2D("cut2_E=7.8","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E78->SetDirectory(fs);
    TH2D * c2E80 = new TH2D("cut2_E=8.0","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E80->SetDirectory(fs);
    TH2D * c2E82 = new TH2D("cut2_E=8.2","", 700, 3.0, 10.0, 6, 0.0, 12.0);
    c2E82->SetDirectory(fs);

    ifstream f2E72("cut-2-2h-jpsi-7.2.dat");
    f2E72.getline(dum, 300);
    while(f2E72 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E72->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E72.close();
    
    ifstream f2E74("cut-2-2h-jpsi-7.4.dat");
    f2E74.getline(dum, 300);
    while(f2E74 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E74->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E74.close();

    ifstream f2E76("cut-2-2h-jpsi-7.6.dat");
    f2E76.getline(dum, 300);
    while(f2E76 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E76->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E76.close();

    ifstream f2E78("cut-2-2h-jpsi-7.8.dat");
    f2E78.getline(dum, 300);
    while(f2E78 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E78->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E78.close();

    ifstream f2E80("cut-2-2h-jpsi-8.0.dat");
    f2E80.getline(dum, 300);
    while(f2E80 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E80->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E80.close();

    ifstream f2E82("cut-2-2h-jpsi-8.2.dat");
    f2E82.getline(dum, 300);
    while(f2E82 >> k >> ds >> q >> theta >> pmin >> pmax >> tmp >> tmp){
      c2E82->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI + zero);
    }
    f2E82.close();

    fs->Write();
  }


  if (false){
    double q, k, theta, ds, pmin, pmax, tmp;
    char dum[300];
    
    TFile * fs = new TFile("harrymodel-phi.root", "RECREATE");
    
    TH2D * hE13 = new TH2D("E=1.3","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE13->SetDirectory(fs);
    ifstream file13("file-phi-1.3.dat");
    file13.getline(dum, 300);
    while(file13 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE13->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file13.close();
    
    TH2D * hE14 = new TH2D("E=1.4","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE14->SetDirectory(fs);
    ifstream file14("file-phi-1.4.dat");
    file14.getline(dum, 300);
    while(file14 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE14->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file14.close();
    
    TH2D * hE15 = new TH2D("E=1.5","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE15->SetDirectory(fs);
    ifstream file15("file-phi-1.5.dat");
    file15.getline(dum, 300);
    while(file15 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE15->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file15.close();
    
    TH2D * hE16 = new TH2D("E=1.6","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE16->SetDirectory(fs);
    ifstream file16("file-phi-1.6.dat");
    file16.getline(dum, 300);
    while(file16 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE16->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file16.close();
    
    TH2D * hE17 = new TH2D("E=1.7","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE17->SetDirectory(fs);
    ifstream file17("file-phi-1.7.dat");
    file17.getline(dum, 300);
    while(file17 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE17->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file17.close();
    
    TH2D * hE18 = new TH2D("E=1.8","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE18->SetDirectory(fs);
    ifstream file18("file-phi-1.8.dat");
    file18.getline(dum, 300);
    while(file18 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE18->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file18.close();
    
    TH2D * hE19 = new TH2D("E=1.9","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE19->SetDirectory(fs);
    ifstream file19("file-phi-1.9.dat");
    file19.getline(dum, 300);
    while(file19 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE19->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file19.close();
    
    TH2D * hE20 = new TH2D("E=2.0","", 400, 0.0+eps, 2.0+eps, 11, 0.0, 22.0);
    hE20->SetDirectory(fs);
    ifstream file20("file-phi-2.0.dat");
    file20.getline(dum, 300);
    while(file20 >> k >> ds >> tmp >> theta >> tmp >> tmp){
      hE20->Fill(k, theta+eps, ds * sin((theta+1.0) / 180.0 * M_PI) * 2.0 * M_PI);
    }
    file20.close();
    
    
    fs->Write();
  }
  
  return 0;
}
