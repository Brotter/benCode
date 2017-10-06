#include "Adu5Pat.h"
#include "BaseList.h"
#include "AntarcticaGeometry.h"
#include "baseCluster.C"

//needed to import things
#include "loadAll.C"

string dataDateString = "09.27.17_19h";


/*

  Super simple clustering

  Also now it does a billion things

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
  if (distance > 1000e3) {
    //    cout << "distance: " << distance << " ";
    //    cout << "A: " << latA << "," << lonA << "," << altA << " ";
    //    cout << "B: " << gpsB->latitude << "," << gpsB->longitude << "," << gpsB->altitude << endl;
    return -9999;}

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
      
  //  cout << "diff: " << diff << endl;
  
  return diff;

}



/*--- Major Code Piece!  Determines which events do not cluster! ---*/
void clusterEvents(string inFileName="cuts.root",string outFileName="clusterEvents.root") {
  /*
    Takes in a list of events and compares their source locations to determine whether they cluster

    Most useful thing it outputs (into cluster.root) are:
    hClosest: a 2d histogram that stores how close all other events are in log-likilhood units
    gClosest: a TGraph that stores the "closest" event to each input event

    Takes forever since it is a O(N^2) process
   */

  cout << "Using input file: " << inFileName << endl;
  TFile *inFile = TFile::Open(inFileName.c_str());

  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");

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
    AnitaEventSummary *eventSummaryA = new AnitaEventSummary(*summary);
    //where event a was captured from (A)
    UsefulAdu5Pat *usefulGPSA = new UsefulAdu5Pat(gps);
    
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

      double distance = calcClusterDistance(eventSummaryA,usefulGPSA,summary,usefulGPSB);
      if (distance == -9999) continue;

      hCluster->Fill(eventCountA,distance);
    
      if (distance < closest) closest = distance;
  
      delete usefulGPSB;
    }

    delete usefulGPSA;
    delete eventSummaryA;
    gClosest->SetPoint(gClosest->GetN(),eventNumberA,closest);
    hClosest->Fill(closest);

    cout << "eventCountA:" << eventCountA << " eventNumber " << eventNumberA << " closest:" << closest << endl;
    
  }
  
  hCluster->Draw("colz");


  TFile* outFile = new TFile(outFileName.c_str(),"recreate");
  hCluster->Write();
  gClosest->Write();
  hClosest->Write();
  outFile->Close();
  
  return;
}
/* --------------------- */


void categorizeClusters(double threshold) {

  /*

    more complicated way to cluster that categorizes things instead of just comparing them vs each other


    UNFINISHED
   */

  //grab the file with all the events that survive cuts
  TFile *inFile = TFile::Open("cuts.root");
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  if (summaryTree == NULL) {
    cout << "Couldn't find summaryTree in cuts.root! Quitting" << endl;
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

  //initialize the vector of entry numbers (remove them when they've been categorized
  vector<int> entryNumbers;
  for (int entry=0; entry<lenEntries; entry++) {
    entryNumbers.push_back(entry);
  }

  //initialize a vector that will store the clusters (variable size)
  vector< vector <int> > clusterVector;


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



void saveEventsNearCandidates(double threshold, string outFileName="") {
  /*
    For my anthropogenic background estimate, I want to see what the distrubtions of all the events that _would_ have 
    clustered with any given "candidate" event.  So this makes a root file with all those events.  Hopefully not too many!

    Lets pull in the cluster.root file with all the clustered events, find the ones that pass the threshold, then
    build them all at once

   */
  stringstream name,title;


  //get ALL the events
  TChain *summaryTree = (TChain*)loadAll(dataDateString,false);

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
  AnitaTemplateSummary *candidateTemplates[numCandidateEvs];
  UsefulAdu5Pat *candidateGPS[numCandidateEvs];

  TFile *cutFile = TFile::Open("cuts.root");
  if (!cutFile->IsOpen()) {
    cout << "Couldn't find cuts.root, quitting " << endl;
    return;
  }
  TTree *cutTree = (TTree*)cutFile->Get("summaryTree");
  if (cutTree == NULL) {
    cout << "couldn't find summaryTree in cuts.root" << endl;
    return;
  }
  cutTree->BuildIndex("eventNumber");
  AnitaEventSummary *cutSum = NULL;
  cutTree->SetBranchAddress("eventSummary",&cutSum);
  AnitaTemplateSummary *tempSum = NULL;
  cutTree->SetBranchAddress("template",&tempSum);
  Adu5Pat *gpstemp = NULL;
  cutTree->SetBranchAddress("gpsEvent",&gpstemp);

  for (int ev=0; ev<numCandidateEvs; ev++){
    int entry = cutTree->GetEntryNumberWithIndex(candidateEvs[ev]);
    cout << ev << ": event number " << candidateEvs[ev] << " is entry " << entry << endl;
    cutTree->GetEntry(entry);
    candidateSummaries[ev] = (AnitaEventSummary*)cutSum->Clone();
    candidateTemplates[ev] = new AnitaTemplateSummary(*tempSum);
    candidateGPS[ev] = new UsefulAdu5Pat(gpstemp);

  }
  cutFile->Close();

  //make an output file and trees, and save all the candidate info
  if (outFileName == "") outFileName = "candidateClustering.root";
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  TTree *outTree[numCandidateEvs];
  for (int ev=0; ev<numCandidateEvs; ev++) {
    name.str("");
    name << "ev" << candidateEvs[ev] << "Tree";
    outTree[ev] = new TTree(name.str().c_str(),name.str().c_str());
    outTree[ev]->Branch("eventSummary",&eventSummary);
    outTree[ev]->Branch("template",&templateSummary);
    outTree[ev]->Branch("noiseSummary",&noiseSummary);
    outTree[ev]->Branch("gpsEvent",&gps);

    name.str("");
    name << "ev" << candidateEvs[ev] << "Summary";
    outFile->WriteObject(candidateSummaries[ev],name.str().c_str());
    name.str("");
    name << "ev" << candidateEvs[ev] << "Template";
    outFile->WriteObject(candidateTemplates[ev],name.str().c_str());
    name.str("");
    name << "ev" << candidateEvs[ev] << "UsefulGps";
    outFile->WriteObject(candidateGPS[ev],name.str().c_str());

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


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;

  for (int entryB=0; entryB<lenEntries; entryB++) {

    if (entryB%10000 == 0) {
      cout << entryB << "/" << lenEntries << " ( ";
      for (int ev=0; ev<numCandidateEvs; ev++) {
	cout << close[ev] << " ";
      }
      cout << ") ";
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entryB)/totalTimeSec;
      double remaining = (float(lenEntries-entryB)/rate)/60.;
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left" << endl;
      watch.Start();
    }
    
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



/*-------------------
  Candidate Clustering Stuff

 */



void drawCandidateClusters(double threshold,bool save=false) {
  /*

    Draws and saves a bunch of images for each "candidate" (determined by `threshold`) that show their distributions
    compared to a pre-made list of events that cluster with that candidate.  So you have to run saveEventsNearCandidates() first

   */


  stringstream name;
  
  //open the results file
  TFile *inFile = TFile::Open("candidateCluster_40.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't open file" << endl;
    return;
  }

  
  //parsing the results file is a mess I don't want to deal with, so just get the "candidate" events again
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


  //get the "background histograms"
  TFile *backgroundFile = TFile::Open("minbiasBackgrounds.root");
  if (!backgroundFile->IsOpen()) {
    cout << "background file not found, run makeMinbiasBackgroundHist() !" << endl;
    return;
  }
  TH1D* backMapPeak  = (TH1D*)backgroundFile->Get("histMapPeak");
  TH1D* backMapSNR   = (TH1D*)backgroundFile->Get("histMapSNR");
  TH1D* backTemplate = (TH1D*)backgroundFile->Get("histTemplate");
  TH1D* backHilbert  = (TH1D*)backgroundFile->Get("histHilbert");


  //make canvas
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);

  //look through all candidates
  for (int candNum=0; candNum<numCandidateEvs; candNum++) {

    name.str("");
    name << "ev" << candidateEvs[candNum] << "Tree";
    TTree* currTree = (TTree*)inFile->Get(name.str().c_str());

    c1->Clear();
    c1->Divide(2,2);
    c1->SetLogy();

    TVirtualPad *pad1 = c1->cd(1);
    pad1->SetLogy();
    TH1D* frontMapPeak  = new TH1D("frontMapPeak","Map Peak; Map Peak; Count",350,0,0.35);
    backMapPeak->SetTitle(frontMapPeak->GetTitle());
    backMapPeak->Draw();
    frontMapPeak->SetLineColor(kRed);
    currTree->Draw("peak[0][0].value >> frontMapPeak","","same");
    name.str("");
    name << "eventNumber == " << candidateEvs[candNum];
    int n = currTree->Draw("peak[0][0].value",name.str().c_str(),"goff");
    TGraph *gCandidate1 = new TGraph();
    gCandidate1->SetPoint(0,currTree->GetV1()[0],10);
    gCandidate1->SetMarkerStyle(29);
    gCandidate1->SetMarkerColor(kRed);
    gCandidate1->Draw("psame");
    
    
    TVirtualPad *pad2 = c1->cd(2);
    pad2->SetLogy();
    TH1D* frontMapSNR   = new TH1D("frontMapSNR","Map SNR; Map SNR; Count",450,0,45);
    backMapSNR->SetTitle(frontMapSNR->GetTitle());
    backMapSNR->Draw();
    frontMapSNR->SetLineColor(kRed);
    currTree->Draw("peak[0][0].snr >> frontMapSNR","","same");
    n = currTree->Draw("peak[0][0].snr",name.str().c_str(),"goff");
    TGraph *gCandidate2 = new TGraph();
    gCandidate2->SetPoint(0,currTree->GetV1()[0],10);
    gCandidate2->SetMarkerStyle(29);
    gCandidate2->SetMarkerColor(kRed);
    gCandidate2->Draw("psame");
    
    TVirtualPad *pad3 = c1->cd(3);
    pad3->SetLogy();
    TH1D* frontTemplate = new TH1D("frontTemplate","Template Correlation; Template Correlation; Count",100,0,1);
    backTemplate->SetTitle(frontTemplate->GetTitle());
    backTemplate->Draw();
    frontTemplate->SetLineColor(kRed);
    currTree->Draw("template.coherent[0][0].cRay[4] >> frontTemplate","","same");
    n = currTree->Draw("template.coherent[0][0].cRay[4]",name.str().c_str(),"goff");
    TGraph *gCandidate3 = new TGraph();
    gCandidate3->SetPoint(0,currTree->GetV1()[0],10);
    gCandidate3->SetMarkerStyle(29);
    gCandidate3->SetMarkerColor(kRed);
    gCandidate3->Draw("psame");
    
    TVirtualPad *pad4 = c1->cd(4);
    pad4->SetLogy();
    TH1D* frontHilbert  = new TH1D("frontHilbert","Deconvolved Hilbert Peak; Hilbert Peak; Count",800,0,800);
    backHilbert->SetTitle(frontHilbert->GetTitle());
    backHilbert->Draw();
    frontHilbert->SetLineColor(kRed);
    currTree->Draw("deconvolved_filtered[0][0].peakHilbert >> frontHilbert","","same");
    n = currTree->Draw("deconvolved_filtered[0][0].peakHilbert",name.str().c_str(),"goff");
    TGraph *gCandidate4 = new TGraph();
    gCandidate4->SetPoint(0,currTree->GetV1()[0],10);
    gCandidate4->SetMarkerStyle(29);
    gCandidate4->SetMarkerColor(kRed);
    gCandidate4->Draw("psame");
    
    if (save) {
      name.str("");
      name << "plots/ev" << candidateEvs[candNum] << "_four.png";
      c1->SaveAs(name.str().c_str());
    }      
      
    delete frontMapPeak;
    delete frontMapSNR;
    delete frontTemplate;
    delete frontHilbert;
    delete gCandidate1;
    delete gCandidate2;
    delete gCandidate3;
    delete gCandidate4;

  }

  return;

}



void draw2DCandidateClusters(double threshold,bool save=false) {
  /*
    Same as the 1d counterpart
   */

  stringstream name;
  
  //open the results file
  TFile *inFile = TFile::Open("candidateCluster_40.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't open file" << endl;
    return;
  }

  
  //parsing the results file is a mess I don't want to deal with, so just get the "candidate" events again
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


  //get the "background histograms"
  TFile *backgroundFile = TFile::Open("minbiasBackgrounds.root");
  if (!backgroundFile->IsOpen()) {
    cout << "background file not found, run makeMinbiasBackgroundHist() !" << endl;
    return;
  }
  TH2D* back1 = (TH2D*)backgroundFile->Get("histMapPeakVsTemplate");
  TH2D* back2 = (TH2D*)backgroundFile->Get("histMapSNRNVsHilbert");

  //make canvas
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);

  //look through all candidates
  for (int candNum=0; candNum<numCandidateEvs; candNum++) {
    
    c1->Clear();

    //get candidate info from input file
    name.str("");
    name << "hCluster_" << candidateEvs[candNum];
    TH1D* currHist = (TH1D*)inFile->Get(name.str().c_str());

    name.str("");
    name << "ev" << candidateEvs[candNum] << "Tree";
    TTree* currTree = (TTree*)inFile->Get(name.str().c_str());

    //make graphs for plot one
    TH2D* front1 = new TH2D("front1","Map Peak Value vs cRay +4 Template; Template Correlation Value; Map Peak",100,0,1,350,0,0.35);

    currTree->Draw("peak[0][0].value:template.coherent[0][0].cRay[4] >> front1","","goff");
    name.str("");
    name << "eventNumber == " << candidateEvs[candNum];
    int n = currTree->Draw("peak[0][0].value:template.coherent[0][0].cRay[4]",name.str().c_str(),"");
    TGraph *gCandidate = new TGraph(n,currTree->GetV2(),currTree->GetV1());
    cout << gCandidate->GetX()[0] << " " << gCandidate->GetY()[0] << endl;

    //plot graph one
    TExec *ex1 = new TExec("ex1","gStyle->SetPalette(kInvertedDarkBodyRadiator);");
    back1->SetTitle(front1->GetTitle());
    back1->Draw("colz");
    ex1->Draw();
    back1->Draw("col same");
    

    TExec *ex2 = new TExec("ex2","gStyle->SetPalette(kTemperatureMap)");
    ex2->Draw();
    front1->Draw("col Same");
    gCandidate->SetMarkerStyle(29);
    gCandidate->SetMarkerColor(kRed);

    gCandidate->Draw("pSame");


    if (save) {
      name.str("");
      name << "plots/ev" << candidateEvs[candNum] << "_1.png";
      c1->SaveAs(name.str().c_str());
    }

    delete gCandidate;
    delete front1;

    c1->Clear();

    //make graphs for plot two
    TH2D* front2 = new TH2D("front2","Map SNR vs Deconvolved Hilbert Peak;Deconvolved Hilbert Peak; Map SNR",800,0,800,450,0,45);    

    currTree->Draw("peak[0][0].snr:deconvolved_filtered[0][0].peakHilbert >> front2","","goff");

    name.str("");
    name << "eventNumber == " << candidateEvs[candNum];
    n = currTree->Draw("peak[0][0].snr:deconvolved_filtered[0][0].peakHilbert",name.str().c_str(),"");
    gCandidate = new TGraph(n,currTree->GetV2(),currTree->GetV1());
    cout << gCandidate->GetX()[0] << " " << gCandidate->GetY()[0] << endl;

    //plot graph two
    back2->SetTitle(front2->GetTitle());
    back2->Draw("colz");
    ex1->Draw();
    back2->Draw("col same");

    ex2->Draw();
    front2->Draw("col same");
    
    gCandidate->SetMarkerStyle(29);
    gCandidate->SetMarkerColor(kRed);
    gCandidate->Draw("p same");
    
    if (save) {
      name.str("");
      name << "plots/ev" << candidateEvs[candNum] << "_2.png";
      c1->SaveAs(name.str().c_str());
    }

    delete gCandidate;
    delete front2;

 
  }
  

  return;
}


void saveCandidates(double threshold) {
  /*
    
    Find all the candidates further than threshold from another event, and save their summaries to a separate file

    Uses gClosest generated by cluster.C::clusterEvents() and saved in cluster.root

   */

  //get the closest-distance graph
  TFile *clusterFile = TFile::Open("cluster.root");
  if (!clusterFile->IsOpen()) {
    cout << "Couldn't find cluster.root, run clusterEvents() first " << endl;
    return;
  }
  TGraph *clusteredEvs = (TGraph*)clusterFile->Get("gClosest");
 

  //get the summaries for all the "candidates" that passed cuts
  TFile *inFile = TFile::Open("cuts.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't find cuts.root, run separateNotable_fromFile() first " << endl;
    return;
  }
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  if (summaryTree == NULL) {
    cout << "Couldn't find summaryTree in cuts.root" << endl;
    return;
  }
  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *templateSummary = NULL;
  AnitaNoiseSummary *noiseSummary = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&templateSummary);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  summaryTree->BuildIndex("eventNumber");


  //make output tree and output file
  TFile *outFile = TFile::Open("candidates.root","recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");
  outTree->Branch("eventSummary",&evSum); 
  outTree->Branch("eventSummary",&evSum);
  outTree->Branch("template",&templateSummary);
  outTree->Branch("noiseSummary",&noiseSummary);
  outTree->Branch("gpsEvent",&gps);

  
  //loop through clusteredEvents graph and save the ones that are above threshold

  for (int i=0; i<clusteredEvs->GetN(); i++) {
    double closest = clusteredEvs->GetY()[i];
    if (closest > threshold) {
      int evNum = clusteredEvs->GetX()[i];
      int entry = summaryTree->GetEntryNumberWithBestIndex(evNum);
      if (entry < 0) {
	cout << "event not found: " << evNum << " -> " << entry << endl;
	break;
      }
      summaryTree->GetEntry(entry);
      outFile->cd();
      outTree->Fill();
    }
  }
    
  outTree->Write();
  outFile->Close();

  return;
}


void clusterBackground(double threshold=40.,int numSplits=1,int split=0, string outFileName="pseudoBaseEvents.root") {
  /*
    Does the opposite of saveCandidates, instead of picking out things that don't cluster, it finds all the impulsive
    events that DID cluster, then finds ALL the events (impulsive or no) that cluster with those.

    Outputs a summaryTree into pseudoBaseEvents.root with all the "anthropogenic background" events, events near bases
    that produced a few impulsive signals that passed cuts.
    
    Uses gClosest generated by cluster.C::clusterEvents() and saved in cluster.root

    needs two things:
    1) cuts_final.root : the events that have passed all the impulsivity cuts
    2) all the data

   */



  //get the summaries for all the "candidates" that passed cuts
  TFile *inFile = TFile::Open("cuts_final.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't find cuts.root, run separateNotable_fromFile() first " << endl;
    return; }
  TTree *impulsiveTree = (TTree*)inFile->Get("summaryTree");
  if (impulsiveTree == NULL) {
    cout << "Couldn't find summaryTree in cuts.root" << endl;
    return; }
  AnitaEventSummary *evSum_imp = NULL;
  Adu5Pat *gps_imp = NULL;
  impulsiveTree->SetBranchAddress("eventSummary",&evSum_imp);
  impulsiveTree->SetBranchAddress("gpsEvent",&gps_imp);
  
  //lets store all the ones that pass in a vector for easy access later
  vector<UsefulAdu5Pat*> vImpulsiveUsefulGps;
  vector<AnitaEventSummary*> vImpulsiveEvSum;
  for (int i=0; i<impulsiveTree->GetEntries(); i++) {
    //for all the events that got clustered, if the closest event is within the clustering threshold, save it to the vector
    impulsiveTree->GetEntry(i);
    UsefulAdu5Pat *currGPS = new UsefulAdu5Pat(gps_imp);
    vImpulsiveUsefulGps.push_back(currGPS);
    AnitaEventSummary *currEvSum = (AnitaEventSummary*)evSum_imp->Clone();
    vImpulsiveEvSum.push_back(currEvSum);
  }

  int lenImp = vImpulsiveEvSum.size();
  cout << "Found " << lenImp << " events that cluster with each other" << endl;

  //also open up ALL of the events :)
  TChain *summaryTree = loadAll(dataDateString,false);
  int lenEntries = summaryTree->GetEntries();
  cout << "Opened up all the data, found " << lenEntries << " entries" << endl;
  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  

  //make output tree and output file
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");
  outTree->Branch("eventSummary",&evSum);
  outTree->Branch("template",&tempSum);
  outTree->Branch("noiseSummary",&noiseSum);
  outTree->Branch("gpsEvent",&gps);
  //also something to save the event number it clustered with just for fun
  // The "seed" is the point on the ice that accumulates all the events nearby
  // seedIndexNumber is just the entry number in the staring cuts_final.root tree
  int seedEventNumber,seedIndexNumber;
  outTree->Branch("seedEventNumber",&seedEventNumber);
  outTree->Branch("seedIndexNumber",&seedIndexNumber);
  //also the cluster "value"
  double clusterValue;
  outTree->Branch("clusterValue",&clusterValue);


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;

  //I'll need to split this up onto the servers, because it takes forever
  int startEntry,stopEntry;
  if (numSplits == 1) {
    startEntry=0;
    stopEntry=lenEntries;
  }
  else {
    lenEntries /= numSplits;
    startEntry = split*lenEntries;
    stopEntry = (split+1)*lenEntries;
    cout << "Splitting into " << numSplits << " sections, which means " << lenEntries << " events per section" << endl;
    cout << "Doing section: " << split << ", starting at entry " << startEntry << " and stopping at " << stopEntry << endl;
  }

  //then loop through ALL the events, and if it clusters with one of the things in the vector, save it
  int savedCount=0;
  for (int entry=startEntry; entry<stopEntry; entry++) {
    //printing stuff
    if (entry%10000 == 0 && entry>0) {
      int printEntry = entry-startEntry;
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(printEntry)/totalTimeSec;
      double remaining = (float(stopEntry-entry)/rate)/60.;
      watch.Start();
      cout << printEntry << "/" << stopEntry-startEntry << " (" << savedCount << ") ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;

    }
    //get entry
    summaryTree->GetEntry(entry);
    //get UsefulAdu5Pat
    UsefulAdu5Pat *currGPS = new UsefulAdu5Pat(gps);
    
    //default to no event, -999
    seedEventNumber = -999;
    seedIndexNumber = -999;
    clusterValue = -999;

    //and a temporary location to store the values
    double tempClusterValue = -999;

    //try clustering with the major bases first
    //wais (seed == -1)
    tempClusterValue = calcBaseDistance(evSum,currGPS,
				    AnitaLocations::getWaisLatitude(),
				    AnitaLocations::getWaisLongitude(),
				    AnitaLocations::getWaisAltitude());
    if (tempClusterValue > 0) {
      seedEventNumber = -1;
      seedIndexNumber = -1;
      clusterValue = tempClusterValue;
    }

    //ldb (seed == -2)
    tempClusterValue = calcBaseDistance(evSum,currGPS,
				    AnitaLocations::LATITUDE_LDB,
				    AnitaLocations::LONGITUDE_LDB,
				    AnitaLocations::ALTITUDE_LDB);
    if ((tempClusterValue > 0) && (tempClusterValue < clusterValue)) {
      //      cout << "ldb:" << clusterValue << endl;
      seedEventNumber = -2;
      seedIndexNumber = -2;
      clusterValue = tempClusterValue;
    }
    
		
    //cluster with all the impulsive events if it doesn't match a base
    // this should also find the base that it clusters to BEST, not just first
    for (int imp=0; imp<lenImp; imp++) {
      tempClusterValue = calcClusterDistance(evSum,currGPS,vImpulsiveEvSum[imp],vImpulsiveUsefulGps[imp]);
      //if: value below lowest clustered event && value not -9999
      if ((tempClusterValue < clusterValue) && (tempClusterValue > 0)) {
	seedEventNumber = vImpulsiveEvSum[imp]->eventNumber;
	seedIndexNumber = imp;
	clusterValue = tempClusterValue;
      }
    }

    if (clusterValue < threshold && clusterValue > 0) {
      outFile->cd();
      outTree->Fill();
      savedCount++;
    }

    delete currGPS;

  }


  cout << "Done :)" << endl;
  cout << " Just have to delete some stuff and save" << endl;

  outFile->cd();
  outTree->Write();
  outFile->Close();


  for (int imp=0; imp<lenImp; imp++) {
    delete vImpulsiveUsefulGps[imp];
    delete vImpulsiveEvSum[imp];
  }

  
  return;
}





void drawEventErrorEllipse() {
  /*
    I want a beter visualization for what something being "within the pointing uncertainty" means.

    Also I can show the horizon on it for looking at events
   */


  return;
}


/*********************************************************************************
In case you want to call from the command line you gotta edit this because macros are dumb!
*/
void cluster() {
  cout << "loaded cluster.C" << endl;
}


void cluster(int numSplits, int split,string baseDir) {

  stringstream name;
  name << baseDir << "/pseudoBaseCluster_" << split << ".root";

  //set the first one to some huge number so it includes all events
  clusterBackground(40,numSplits,split,name.str());

  return;
}
