#include "Adu5Pat.h"

/*

  Super simple clustering

 */


bool debug=false;

bool notNotable(AnitaEventSummary *summary) {
  /*
    in case I included ones that aren't interesting

    Double Negative!
    not interesting (skip) = true
    interesting (cluster) = false
  */

  // if it doesn't land on the continent (and therefore can't cluster)
  if (summary->peak[0][0].latitude <= -999) return true;

  return false;
}
  


TChain* importGPS() {

  TChain* outChain = new TChain("adu5PatTree","adu5PatTree");

  char* baseDir = getenv("ANITA_ROOT_DATA");
  stringstream name;
  for (int i=130; i<440; i++) {
    name.str("");
    name << baseDir << "/run" << i << "/gpsEvent" << i << ".root";
    outChain->Add(name.str().c_str());
  }

  cout << "Building gps index" << endl;
  outChain->BuildIndex("eventNumber");

  return outChain;
}


double calcClusterDistance(AnitaEventSummary *eventA, UsefulAdu5Pat *gpsA, AnitaEventSummary *eventB, UsefulAdu5Pat *gpsB) {
  /*

    The function that actual does the "log likelyhood" angle distance calculation.

    It is really just a geometric average angular difference calculator between two events

   */

  
  //where event a was seen (a)
  double thetaA = eventA->peak[0][0].theta;
  double phiA   = eventA->peak[0][0].phi;
  double latA   = eventA->peak[0][0].latitude;
  double lonA   = eventA->peak[0][0].longitude;
  double altA   = eventA->peak[0][0].altitude;
    
  //if event a is further than 1000km from location B, return -9999
  double distance = gpsB->getDistanceFromSource(latA,lonA,altA);
  if (distance > 1000e3) return -9999;

  //where event b was seen (b)
  double thetaB = eventB->peak[0][0].theta;
  double phiB   = eventB->peak[0][0].phi;
  double latB   = eventB->peak[0][0].latitude;
  double lonB   = eventB->peak[0][0].longitude;
  double altB   = eventB->peak[0][0].altitude;      
    
  //where event a is seen from B's location
  double thetaBA,phiBA; 
  gpsB->getThetaAndPhiWave(lonA,latA,altA,thetaBA,phiBA);
  thetaBA *= TMath::RadToDeg();
  phiBA *= TMath::RadToDeg();
  //difference between B->b and B->a
  double diffBAtheta = TMath::Abs(thetaBA - thetaA);
  double diffBAphi = TMath::Abs(FFTtools::wrap(phiBA - phiB,360,0));
  if (diffBAphi > 180) diffBAphi = 360 - diffBAphi;
    
  //where event b is seen from A's location
  double thetaAB,phiAB; 
  gpsA->getThetaAndPhiWave(lonB,latB,altB,thetaAB,phiAB);
  thetaAB *= TMath::RadToDeg();
  phiAB   *= TMath::RadToDeg();
  //difference between A->a and A->b
  double diffABtheta = TMath::Abs(thetaAB - thetaB);
  double diffABphi = TMath::Abs(FFTtools::wrap(phiAB - phiA,360,0));
  if (diffABphi > 180) diffABphi = 360 - diffABphi;      
    

  const double sigmaTheta = 0.2193;
  const double sigmaPhi = 0.5429;
      
  double diffTheta = TMath::Sqrt(pow(diffABtheta,2) + pow(diffBAtheta,2))/sigmaTheta;
  double diffPhi = TMath::Sqrt(pow(diffABphi,2) + pow(diffBAphi,2))/sigmaPhi;
    
  double diff = TMath::Sqrt(pow(diffTheta,2) + pow(diffPhi,2));
      
  //	cout << entryA << " vs " << entryB << " = " << diff << endl;
  
  return diff;

}








void clusterEventList() {
  /*
    Takes in a list of events and compares their source locations to determine whether they cluster

    Most useful thing it outputs (into cluster.root) are a histogram and TGraph that stores the "closest" event to each input event

    Takes forever since it is a O(N^2) process
   */


  TFile *inFile = TFile::Open("cuts.root");

  TTree *summaryTree = (TTree*)inFile->Get("eventSummary");

  if (summaryTree == NULL) {
    cout << "Couldn't find cutSummary in file! Quitting" << endl;
    return;
  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);


  int lenEntries = summaryTree->GetEntries();
  cout << "found " << lenEntries << " events" << endl;

  TH2D* hCluster = new TH2D("hCluster","hCluster",lenEntries,-0.5,lenEntries-0.5,1000,0,1000);

  TGraph *gClosest = new TGraph();
  gClosest->SetName("gClosest");

  TH1D *hClosest = new TH1D("hClosest","closest event in sigma",5000,0,500);


  int eventCountA=0;
  for (int entryA=0; entryA<lenEntries; entryA++) {

    double closest = 999;

    summaryTree->GetEntry(entryA);

    if (notNotable(summary)) continue;

    eventCountA++;
    int eventNumberA = summary->eventNumber;

    //where event a was captured from (A)
    UsefulAdu5Pat *usefulGPSA = new UsefulAdu5Pat(gps);
    
    //where event a was seen (a)
    double thetaA = summary->peak[0][0].theta;
    double phiA   = summary->peak[0][0].phi;
    double latA = summary->peak[0][0].latitude;
    double lonA = summary->peak[0][0].longitude;
    double altA = summary->peak[0][0].altitude;
    
    

    int eventCountB = 0;
    for (int entryB=0; entryB<lenEntries; entryB++) {
      
      //will always be zero
      if (entryA == entryB) continue;
      
      summaryTree->GetEntry(entryB);

      if (notNotable(summary)) continue;

      eventCountB++;
      int eventNumberB = summary->eventNumber;

      //where event b was captured from (B)
      UsefulAdu5Pat *usefulGPSB = new UsefulAdu5Pat(gps);


      //if event a is further than 1000km from location B, just skip it
      double distance = usefulGPSB->getDistanceFromSource(latA,lonA,altA);
      if (distance > 1000e3) continue;

      //where event b was seen (b)
      double thetaB = summary->peak[0][0].theta;
      double phiB   = summary->peak[0][0].phi;
      double latB = summary->peak[0][0].latitude;
      double lonB = summary->peak[0][0].longitude;
      double altB = summary->peak[0][0].altitude;      
      
      //where event a is seen from B's location
      double thetaBA,phiBA; 
      usefulGPSB->getThetaAndPhiWave(lonA,latA,altA,thetaBA,phiBA);
      thetaBA *= TMath::RadToDeg();
      phiBA *= TMath::RadToDeg();
      //difference between B->b and B->a
      double diffBAtheta = TMath::Abs(thetaBA - thetaA);
      double diffBAphi = TMath::Abs(FFTtools::wrap(phiBA - phiB,360,0));
      if (diffBAphi > 180) diffBAphi = 360 - diffBAphi;
      
      //where event b is seen from A's location
      double thetaAB,phiAB; 
      usefulGPSA->getThetaAndPhiWave(lonB,latB,altB,thetaAB,phiAB);
      thetaAB *= TMath::RadToDeg();
      phiAB   *= TMath::RadToDeg();
      //difference between A->a and A->b
      double diffABtheta = TMath::Abs(thetaAB - thetaB);
      double diffABphi = TMath::Abs(FFTtools::wrap(phiAB - phiA,360,0));
      if (diffABphi > 180) diffABphi = 360 - diffABphi;      


      const double sigmaTheta = 0.2193;
      const double sigmaPhi = 0.5429;
      
      double diffTheta = TMath::Sqrt(pow(diffABtheta,2) + pow(diffBAtheta,2))/sigmaTheta;
      double diffPhi = TMath::Sqrt(pow(diffABphi,2) + pow(diffBAphi,2))/sigmaPhi;
      
      double diff = TMath::Sqrt(pow(diffTheta,2) + pow(diffPhi,2));
      
      //	cout << entryA << " vs " << entryB << " = " << diff << endl;

      hCluster->Fill(eventCountA,diff);
    
      if (diff < closest) closest = diff;
  
      delete usefulGPSB;
      

      if (debug) {
	cout << "( " << eventCountA << " ) " << eventNumberA << " , " << eventNumberB << " ( " << eventCountB << ") | ";
	cout << phiBA << " - " << phiB << " , " << phiAB << " - " << phiA << " =  " << diffPhi << " | ";
	cout << thetaA << " - " << thetaB << " , " << diffBAtheta << " - " << diffABtheta << " = " << diffTheta << " | ";
	cout << distance << endl;
      }
    }

    delete usefulGPSA;
    
    gClosest->SetPoint(gClosest->GetN(),eventNumberA,closest);
    hClosest->Fill(closest);

    cout << "eventCountA:" << eventCountA << " eventNumber " << eventNumberA << " closest:" << closest << endl;
    
  }
  
  hCluster->Draw("colz");


  TFile* outFile = new TFile("cluster.root","recreate");
  hCluster->Write();
  gClosest->Write();
  hClosest->Write();
  outFile->Close();
  
  return;
}


void printNotable(TH2D *hist,double threshold=10) {

  int numBinsX = hist->GetNbinsX();
  int numBinsY = hist->GetNbinsY();

  TH1D *minDists = new TH1D("minDists","closest event angular sigma",50,0,50);

  for (int evBin=0; evBin<numBinsX; evBin++) {
    
    TH1D *temp = hist->ProjectionY("temp",evBin+1,evBin+1);

    double minDist = temp->GetMinimumBin();
    minDists->Fill(minDist);

    delete temp;
  }


  
      
    minDists->Draw();
    

}



void saveEventsNearCandidates(double threshold) {
  /*
    For my anthropogenic background estimate, I want to see what the distrubtions of all the events that _would_ have 
    clustered with any given "candidate" event.  So this makes a root file with all those events.  Hopefully not too many!

    Lets pull in the cluster.root file with all the clustered events, find the ones that pass the threshold, then
    build them all at once

   */
  stringstream name,title;


  //get ALL the events
  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  int lenEntries = summaryTree->GetEntries();
  cout << lenEntries << " total entries found" << endl;

  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  AnitaNoiseSummary *noiseSummary = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);

  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);


  //get the "candidate" events
  TFile *clusterFile = TFile::Open("cluster.root");
  if (!clusterFile->IsOpen()) {
    cout << "Couldn't find cluster.root, quitting " << endl;
    return;
  }
  TGraph *clusteredEvs = (TGraph*)clusterFile->Get("gClosest");
  vector<int> candidateEvs;
  for (int ev=0; ev<clusteredEvs->GetN(); ev++) {
    if (clusteredEvs->GetY()[ev] >= threshold) candidateEvs.push_back(clusteredEvs->GetX()[ev]);
  }
  const int numCandidateEvs = candidateEvs.size();
  clusterFile->Close();

  cout << "Found " << numCandidateEvs << " Candidate Events:" << endl;

  //get the candidates (from cuts.root which has way less entries in it)
  AnitaEventSummary *candidateSummaries[numCandidateEvs];
  UsefulAdu5Pat *candidateGPS[numCandidateEvs];

  TFile *cutFile = TFile::Open("cuts.root");
  if (!cutFile->IsOpen()) {
    cout << "Couldn't find cuts.root, quitting " << endl;
    return;
  }
  TTree *cutTree = (TTree*)cutFile->Get("eventSummary");
  cutTree->BuildIndex("eventNumber");

  AnitaEventSummary *cutSum = NULL;
  cutTree->SetBranchAddress("eventSummary",&cutSum);

  Adu5Pat *gpstemp = NULL;
  cutTree->SetBranchAddress("gpsEvent",&gpstemp);

  for (int ev=0; ev<numCandidateEvs; ev++){
    int entry = cutTree->GetEntryNumberWithIndex(candidateEvs[ev]);
    cout << ev << ": event number " << candidateEvs[ev] << " is entry " << entry << endl;
    cutTree->GetEntry(entry);
    candidateSummaries[ev] = (AnitaEventSummary*)cutSum->Clone();
    candidateGPS[ev] = new UsefulAdu5Pat(gpstemp);
  }
 

  //make an output file and trees
  TFile *outFile = TFile::Open("candidateClustering.root","recreate");
  TTree *outTree[numCandidateEvs];
  for (int ev=0; ev<numCandidateEvs; ev++) {
    name.str("");
    name << "ev" << candidateEvs[ev] << "Tree";
    outTree[ev] = new TTree(name.str().c_str(),name.str().c_str());
    outTree[ev]->Branch("eventSummary",&eventSummary);
    outTree[ev]->Branch("template",&templateSummary);
    outTree[ev]->Branch("noiseSummary",&noiseSummary);
    outTree[ev]->Branch("gpsEvent",&gps);
  }
 

  //histograms to record distance distributions past threshold
  TH1D *hCluster[numCandidateEvs];
  for (int ev=0; ev<numCandidateEvs; ev++) {
    name.str("");
    name << "hCluster_" << candidateEvs[ev];
    title.str("");
    title << "Event Cluster Distance to Ev " << candidateEvs[ev];
    hCluster[ev] = new TH1D(name.str().c_str(),title.str().c_str(),1000,0,threshold*4);
  }


  int close[numCandidateEvs];
  for (int ev=0; ev<numCandidateEvs; ev++) {
    close[ev] = 0;
  }

  for (int entryB=0; entryB<lenEntries; entryB++) {

    if (entryB%10000 == 0) {
      cout << entryB << "/" << lenEntries << " ( ";
      for (int ev=0; ev<numCandidateEvs; ev++) {
	cout << close[ev] << " ";
      }
      cout << ")" << endl; }
    
    summaryTree->GetEntry(entryB);
    
    if (notNotable(eventSummary)) continue;

    //where event b was captured from (B)
    UsefulAdu5Pat *usefulGPSB = new UsefulAdu5Pat(gps);

    for (int ev=0; ev<numCandidateEvs; ev++) {
      double dist = calcClusterDistance(candidateSummaries[ev],candidateGPS[ev],eventSummary,usefulGPSB);
      hCluster[ev]->Fill(dist);
      if (dist <= threshold && dist != -9999) {
	outTree[ev]->Fill();
	close[ev]++;
      }
    }

    delete usefulGPSB;
  }      

  cout << "Finished!" << endl;

  for (int ev=0; ev<numCandidateEvs; ev++) {
    hCluster[ev]->Write();
    outTree[ev]->Write();
  }

  outFile->Close();

}
