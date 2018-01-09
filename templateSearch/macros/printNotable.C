/*

  Copying and pasting into spreadsheets is stoops

 */
#include "AnitaEventSummary.h"

void printNotable(string inFileName) {

  stringstream name;

  TFile *inFile = TFile::Open(inFileName.c_str());
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
  //  TString *labelString = NULL;
  //  summaryTree->SetBranchAddress("label",&labelString);
  AnitaEventSummary::WaveformInfo *newStokes = NULL;
  summaryTree->SetBranchAddress("newStokes",&newStokes);


  size_t pos = inFileName.find(".root");
  string baseName = inFileName.substr(0,pos);
  string outFileName = baseName+"_value.csv";
  ofstream outfile(outFileName);
  
  //  outfile << "candidateNum eventNumber realTime sourceLat sourceLon sourceAlt anitaLat anitaLon anitaAlt elevation trueAz coherTemplateCorr deconvTemplateCorr mapPeak mapSNR mapTheta cohLinPolFrac deconvLinPolFrac cohTotPolFrac deconvTotPolFrac cohHilbPeak deconvHilbPeak" << endl;
  outfile << "candidateNum eventNumber realTime sourceLat sourceLon sourceAlt anitaLat anitaLon anitaAlt elevation trueAz newPolPlane" << endl;
  int candNum=0;
  for (int i=0; i<lenEntries; i++) {
    summaryTree->GetEntry(i);
    
    //    if ( (clusterThreshold>0 && (clusterValue > clusterThreshold)) || 
    //	  (evSum->eventNumber==27142546 || evSum->eventNumber==39599205) ) {

    //    if (!(strstr(labelString->Data(),"Clustered Passing"))) continue;
    

    AnitaEventSummary::WaveformInfo cohrnt = evSum->coherent_filtered[0][0];
    AnitaEventSummary::WaveformInfo deconv = evSum->deconvolved_filtered[0][0];
    
    AnitaEventSummary::PointingHypothesis peak = evSum->peak[0][0];
    AnitaEventSummary::PayloadLocation anita = evSum->anitaLocation;

    
    outfile << candNum << " " << evSum->eventNumber << " " << evSum->realTime << " ";
    outfile << peak.latitude << " " << peak.longitude << " " << peak.altitude << " ";
    outfile << anita.latitude << " " << anita.longitude << " " << anita.altitude << " ";
    outfile << peak.theta << " " << FFTtools::wrap(anita.heading-peak.phi,360,0)+180 << " "; //true azimuth. 0-360 CW from north
    outfile << newStokes->linearPolAngle() << " ";
    /*
    outfile << tempSum->coherent[0][0].cRay[4] << " " << tempSum->deconvolved[0][0].cRay[4] << " ";
    outfile << peak.value << " " << peak.snr << " " << peak.theta << " ";
    outfile << cohrnt.linearPolFrac() << " " << deconv.linearPolFrac() << " ";
    outfile << cohrnt.totalPolFrac() << " " << deconv.totalPolFrac() << " ";
    outfile << cohrnt.peakHilbert << " " << deconv.peakHilbert;
    */
    outfile << endl;
    candNum++;
    
  }
  cout << "Found " << candNum << " within threshold" << endl;
  
  outfile.close();


}


void printGeoAssociatedSets() {
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add("trueCandidates_oct14_remasked.root");
  int lenEntries = summaryTree->GetEntries();

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);

  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);

    string filename = "geoAssociated/geoAssociated_ev" + to_string(evSum->eventNumber) + "_newStokes.root";
    cout << filename << endl;

    printNotable(filename);
  }

  return;
}

void printNotable() {
  cout << "loaded printNotable.C" << endl;
  return;
}
