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


  TFile* outFile = new TFile("cluster.root","recreate");
  hCluster->Write();
  gClosest->Write();
  hClosest->Write();
  outFile->Close();
  
  return;
}



void categorizeClusters(double threshold) {

  /*

    more complicated way to cluster that categorizes things instead of just comparing them vs each other

   */

  //grab the file with all the events that survive cuts
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
  AnitaTemplateSummary *candidateTemplates[numCandidateEvs];
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



void makeMinbiasBackgroundHist() {
  /*
    drawCanidateClusters needs background plots, gotta run this to make 'em
  */

  TFile *inFile = TFile::Open("07.28.17_17h_decimated.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't open file" << endl;
    return;
  }
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");

  int lenEntries = summaryTree->GetEntries();
  cout << "got " << lenEntries << " entries" << endl;

  TFile *outFile = TFile::Open("minbiasBackgrounds.root","recreate");

  TH2D* histMapPeakVsTemplate = new TH2D("histMapPeakVsTemplate","minbias",100,0,1,350,0,0.35);
  summaryTree->Draw("peak[0][0].value:template.coherent[0][0].cRay[4] >> histMapPeakVsTemplate","!flags.isRF","colz");
  histMapPeakVsTemplate->SetStats(0);
  histMapPeakVsTemplate->Write();

  TH2D* histMapSNRVsHilbert = new TH2D("histMapSNRVsHilbert","minbias",800,0,800,450,0,45);
  summaryTree->Draw("peak[0][0].snr:deconvolved_filtered[0][0].peakHilbert >> histMapSNRVsHilbert","!flags.isRF","colz");
  histMapSNRVsHilbert->SetStats(0);
  histMapSNRVsHilbert->Write();

  TH1D* histMapPeak = new TH1D("histMapPeak","minbias;Map Peak; Count",350,0,0.35);
  summaryTree->Draw("peak[0][0].value >> histMapPeak","!flags.isRF");
  histMapPeak->SetStats(0);
  histMapPeak->Write();

  TH1D* histMapSNR = new TH1D("histMapSNR","minbias;Map SNR; Count",450,0,45);
  summaryTree->Draw("peak[0][0].snr >> histMapSNR","!flags.isRF");
  histMapSNR->SetStats(0);
  histMapSNR->Write();

  TH1D* histTemplate = new TH1D("histTemplate","minbias; cRay +4 Template Correlation; Count",100,0,1);
  summaryTree->Draw("template.coherent[0][0].cRay[4] >> histTemplate","!flags.isRF");
  histTemplate->SetStats(0);
  histTemplate->Write();

  TH1D* histHilbert = new TH1D("histHilbert","minbias; Deconvolved Hilbert Peak; Count",800,0,800);
  summaryTree->Draw("deconvolved_filtered[0][0].peakHilbert >> histHilbert","!flags.isRF");
  histHilbert->SetStats(0);
  histHilbert->Write();

  outFile->Close();

  return;
}

void drawCandidateClusters(double threshold,bool save=false) {
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







/*
  Clustering with bases

*/



TTree *getBaseTree() {

  stringstream name;

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");

  return baseCampTree;
}


TGraph* getBaseGraph() {
  /*
    First get the bases!
  */


  stringstream name;

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  TTree *awsCampTree = (TTree*)fBaseList->Get("awsTree");
  double fullLat,fullLon;
  double fullLataws,fullLongaws;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLon);

  awsCampTree->SetBranchAddress("fullLat",&fullLataws);
  awsCampTree->SetBranchAddress("fullLong",&fullLongaws);

  TGraph *gBaseList = new TGraph();
  gBaseList->SetName("gBaseList_latlon");

  cout << "baseCampTree->GetEntries() = " << baseCampTree->GetEntries() << endl;
  for (int entry=0; entry<baseCampTree->GetEntries(); entry++) {
    baseCampTree->GetEntry(entry);
    gBaseList->SetPoint(entry,fullLat,fullLon);
  }


  cout << "awsCampTree->GetEntries() = " << awsCampTree->GetEntries() << endl;
  for (int entry=0; entry<awsCampTree->GetEntries(); entry++) {
    awsCampTree->GetEntry(entry);
    gBaseList->SetPoint(gBaseList->GetN(),fullLat,fullLon);
  }

  return gBaseList;
}



double calcBaseDistance(AnitaEventSummary *event, UsefulAdu5Pat *gps,double lat, double lon,double alt) {
  /*

    For calculating the clustering to bases which only have one definate known location

   */

  if (gps->getDistanceFromSource(lon,lat,alt) > 1000e3) {
    return -9999;
  }


  //where the base is
  double thetaBase,phiBase;
  gps->getThetaAndPhiWave(lon,lat,alt,thetaBase,phiBase);
  thetaBase *= TMath::RadToDeg();
  phiBase   *= TMath::RadToDeg();

  //where the event is
  double phi =  event->peak[0][0].phi;
  double theta = event->peak[0][0].theta;
  
  //difference between A->a and A->b
  double diffTheta = TMath::Abs(theta - thetaBase);
  double diffPhi = TMath::Abs(FFTtools::wrap(phi - phiBase,360,0));
  if (diffPhi > 180) diffPhi = 360 - diffPhi;
    

  const double sigmaTheta = 0.2193;
  const double sigmaPhi = 0.5429;
      
  double diff = TMath::Sqrt(pow(diffTheta/sigmaTheta,2) + pow(diffPhi/sigmaPhi,2));
      
  //  cout << "diff: " << diff << endl;
  
  return diff;

}



void saveEventsNearBases(double threshold=40.,int numSplits=1,int split=0, string outFileName="") {
/*
  This will take up some space, but otherwise histograms take forever to generate

  Save every event that passes a clustering threshold, default was 40 but probably should be different, near a recorded
  base from jruss's list, or on one of the events that passed all the other cuts (which I also should probably tune!!)

  Basically identical to saveEventsNearCandidates

 */
  stringstream name,title;


  //get ALL the events
  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C(\"07.28.17_17h/\",false)");

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


  //get the bases
  TTree* baseTree = getBaseTree();
  double lat,lon,alt;
  baseTree->SetBranchAddress("fullLat",&lat);
  baseTree->SetBranchAddress("fullLong",&lon);
  baseTree->SetBranchAddress("alt",&alt);

  int numBases = baseTree->GetEntries();

  cout << "Found " << numBases << " Bases~" << endl;

  //make an output file and trees, and save all the candidate info
  if (outFileName == "") outFileName = "baseClustering.root";
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  TTree *outTree = new TTree("eventSummary","eventSummary");
  outTree->Branch("eventSummary",&eventSummary);
  outTree->Branch("template",&templateSummary);
  outTree->Branch("noiseSummary",&noiseSummary);
  outTree->Branch("gpsEvent",&gps);
 

  //histograms to record distance distributions past threshold
  TH1D *hCluster[numBases];
  for (int base=0; base<numBases; base++) {
    name.str("");
    name << "hCluster_" << base;
    title.str("");
    title << "Event Cluster Distance to Base Num " << base;
    hCluster[base] = new TH1D(name.str().c_str(),title.str().c_str(),1000,0,threshold*4);
  }


  int close[numBases];
  for (int base=0; base<numBases; base++) {
    close[base] = 0;
  }


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;



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

  
  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10000 == 0) {
      int printEntry = entry-startEntry;
      cout << printEntry << "/" << lenEntries << " ( ";	
      for (int base=0; base<numBases; base++) {
	cout << close[base] << " ";
      }
      cout << ") ";
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(printEntry)/totalTimeSec;
      double remaining = (float(lenEntries-printEntry)/rate)/60.;
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left" << endl;
      watch.Start();
    }
    
    summaryTree->GetEntry(entry);
    
    if (notNotable(eventSummary)) continue;

    //where event b was captured from (B)
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);


    bool filled=false;
    
    for (int base=0; base<numBases; base++) {
      baseTree->GetEntry(base);
      double tempAlt = alt;
      if (tempAlt<0) tempAlt=0;
      double dist = calcBaseDistance(eventSummary,usefulGPS,lat,lon,tempAlt);
      hCluster[base]->Fill(dist);
      if (dist <= threshold && dist != -9999) {
	close[base]++;
	if (!filled) {
	  outTree->Fill();
	  filled = true;
	}
      }
    }

    delete usefulGPS;
  }      

  cout << "Finished!" << endl;

  for (int base=0; base<numBases; base++) {
    hCluster[base]->Write();  
  }
  outTree->Write();
  outFile->Close();

  return;

}




void mergeClusterHistograms(int numCores=32,int numBases=104,string date="08.23.17_18h") {
  /*
    I do clustering with the cluster servers so I gotta merge the results by hand
   */


  stringstream name;

  TChain *eventSummary = new TChain("eventSummary","eventSummary");


  TH1D *hCluster[numBases];
  TList *histList[numBases];
  for (int base=0; base<numBases; base++) {
    histList[base] = new TList;
  }

  char* basedir = getenv("ANITA3_RESULTSDIR");

  for (int i=0; i<numCores; i++) {
    name.str("");
    name << basedir << "cluster/" << date << "/baseCluster_" << i << ".root";
    cout << "loading: " << name.str() << endl;

    eventSummary->Add(name.str().c_str());

    TFile *inFile = TFile::Open(name.str().c_str());

    for (int base=0; base<numBases; base++) {
      name.str("");
      name << "hCluster_" << base;
      TH1D *currHist = (TH1D*)inFile->Get(name.str().c_str());
      if (i==0) hCluster[base] = (TH1D*)currHist->Clone();
      cout << name.str() << " " << currHist->GetEntries() << endl;
      histList[base]->Add(currHist);
    }

  }

  cout << "merging" << endl;
  for (int base=0; base<numBases; base++) {
    cout << "base:" << base << endl;
    hCluster[base]->Reset();
    hCluster[base]->Merge(histList[base]);

  }
  
    
  TFile *outFile = TFile::Open("mergeBaseClusters.root","recreate");
  for (int base=0; base<numBases; base++) {
    hCluster[base]->Write();
  }
  
  cout << "Copying summary" << endl;
  eventSummary->CloneTree(-1,"fast");
  outFile->Write();
  

  outFile->Close();

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
  name << baseDir << "/baseCluster_" << split << ".root";

  saveEventsNearBases(40,numSplits,split,name.str());

  return;
}
