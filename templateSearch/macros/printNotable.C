/*

  Copying and pasting into spreadsheets is stoops

 */
#include "AnitaEventSummary.h"


void printNotable(string fileNameBase="candidates_oct13") {

  stringstream name;
  name.str("");
  name << fileNameBase << ".root";

  TFile *inFile = TFile::Open(name.str().c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");

  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " to print out" << endl;

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("template",&tempSum);

  name.str("");
  name << fileNameBase << "_values.txt";
  ofstream outfile(name.str());

  outfile << "candidateNum eventNumber templateCorr mapPeak mapSNR mapTheta cohLinPolFrac deconvLinPolFrac cohTotPolFrac deconvTotPolFrac cohHilbPeak deconvHilbPeak" << endl;
  
  for (int i=0; i<lenEntries; i++) {
    summaryTree->GetEntry(i);
    
    AnitaEventSummary::WaveformInfo cohrnt = evSum->coherent_filtered[0][0];
    AnitaEventSummary::WaveformInfo deconv = evSum->deconvolved_filtered[0][0];

    outfile << i << " " << evSum->eventNumber << " " << tempSum->coherent[0][0].cRay[4] << " ";
    outfile << evSum->peak[0][0].value << " " << evSum->peak[0][0].snr << " " << evSum->peak[0][0].theta << " ";
    outfile << cohrnt.linearPolFrac() << " " << deconv.linearPolFrac() << " ";
    outfile << cohrnt.totalPolFrac() << " " << deconv.totalPolFrac() << " ";
    outfile << cohrnt.peakHilbert << " " << deconv.peakHilbert << endl;

  }

  outfile.close();


}
