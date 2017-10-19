/*

  Copying and pasting into spreadsheets is stoops

 */
#include "AnitaEventSummary.h"


void printNotable(string fileNameBase="cutsClust_oct14",double clusterThreshold=40) {

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
  double clusterValue;
  summaryTree->SetBranchAddress("clusterValue",&clusterValue);
  clusterValue = 999;

  name.str("");
  name << fileNameBase << "_values.csv";
  ofstream outfile(name.str());

  outfile << "candidateNum eventNumber coherTemplateCorr deconvTemplateCorr mapPeak mapSNR mapTheta cohLinPolFrac deconvLinPolFrac cohTotPolFrac deconvTotPolFrac cohHilbPeak deconvHilbPeak" << endl;
  
  int candNum=0;
  for (int i=0; i<lenEntries; i++) {
    summaryTree->GetEntry(i);
    
    if ( (clusterThreshold>0 && (clusterValue > clusterThreshold)) || 
	  (evSum->eventNumber==27142546 || evSum->eventNumber==39599205) ) {

      AnitaEventSummary::WaveformInfo cohrnt = evSum->coherent_filtered[0][0];
      AnitaEventSummary::WaveformInfo deconv = evSum->deconvolved_filtered[0][0];
    
      outfile << candNum << " " << evSum->eventNumber << " ";
      outfile << tempSum->coherent[0][0].cRay[4] << " " << tempSum->deconvolved[0][0].cRay[4] << " ";
      outfile << evSum->peak[0][0].value << " " << evSum->peak[0][0].snr << " " << evSum->peak[0][0].theta << " ";
      outfile << cohrnt.linearPolFrac() << " " << deconv.linearPolFrac() << " ";
      outfile << cohrnt.totalPolFrac() << " " << deconv.totalPolFrac() << " ";
      outfile << cohrnt.peakHilbert << " " << deconv.peakHilbert << endl;
      candNum++;
    }
  }
  cout << "Found " << candNum << " within threshold" << endl;
  
  outfile.close();


}
