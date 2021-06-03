#include <cstdarg>
#include <iostream>

using namespace std;

int combineTree(double nbatches, TString outputDirectory, TString inputFileInitial)
{
  // There are these pesky, large GEANT4 Node trees files labelled T1 which get outputted by the DST_Reader. It would be optimal to delete these to save on storage
  bool verbose = true;
  
  TString outputFileName = "";
  outputFileName = outputDirectory + TString("/combinedTree.root");
  
  TChain *eictree = new TChain("tree");
  for(int i=0;i<nbatches;i++)
    {
      TString locationOfEICTree = TString(inputFileInitial+TString(Form("/run%d",i))+TString(".root"));
      if(verbose)
	cout << "Adding 'EICTree' " << i+1 << " of " << nbatches << " ... ";
      eictree->Add(locationOfEICTree.Data());
      if(verbose)
	cout << " Added | ";
      if(verbose)
	cout << "\n";
    }
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  eictree->CloneTree(-1,"fast");
  outputFile->Write();
  outputFile->Close();
  return 0;
}
