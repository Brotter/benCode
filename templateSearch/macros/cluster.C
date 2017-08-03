#include "Adu5Pat.h"

/*

  Super simple clustering

 */


bool notNotable(AnitaEventSummary *summary, AnitaTemplateSummary *templateSummary) {
  /*
    in case I included ones that aren't interesting

    Double Negative!
    not interesting (skip) = true
    interesting (cluster) = false
  */

  if (summary->peak[0][0].latitude <= -999) return true;
  if (templateSummary->deconvolved[0][0].impulse < 0.6) return true;

  return false;
}
  


void cluster() {


  TFile *inFile = TFile::Open("notableEvents.root");

  TTree *summaryTree = (TTree*)inFile->Get("cutSummary");

  if (summaryTree == NULL) {
    cout << "Couldn't find cutSummary in file! Quitting" << endl;
    return;
  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  int lenEntries = summaryTree->GetEntries();
  cout << "found " << lenEntries << " events" << endl;

  TH2D* hCluster = new TH2D("hCluster","hCluster",lenEntries,-0.5,lenEntries-0.5,1000,0,1000);

  TGraph *singlets = new TGraph();
  singlets->SetName("singlets");

  for (int entryA=0; entryA<lenEntries; entryA++) {

    summaryTree->GetEntry(entryA);
    if (notNotable(summary,templateSummary)) continue;

    //where event a was captured from (A)
    Adu5Pat *gpsA = new Adu5Pat();
    gpsA->latitude  = summary->anitaLocation.latitude;
    gpsA->longitude = summary->anitaLocation.longitude;
    gpsA->altitude  = summary->anitaLocation.altitude;
    gpsA->heading    = summary->anitaLocation.heading;
    UsefulAdu5Pat *usefulGPSA = new UsefulAdu5Pat(gpsA);
    
    //where event a was seen (a)
    double thetaA = summary->peak[0][0].theta;
    double phiA   = summary->peak[0][0].phi;
    double latA = summary->peak[0][0].latitude;
    double lonA = summary->peak[0][0].longitude;
    double altA = summary->peak[0][0].altitude;
    
    
    int close = 0;
    int eventNum = summary->eventNumber;
    for (int entryB=0; entryB<lenEntries; entryB++) {
      
      //will always be zero
      if (entryA == entryB) continue;
      
      summaryTree->GetEntry(entryB);
      if (notNotable(summary,templateSummary)) continue;

      //where event b was captured from (B)
      Adu5Pat *gpsB = new Adu5Pat();
      gpsB->latitude  = summary->anitaLocation.latitude;
      gpsB->longitude = summary->anitaLocation.longitude;
      gpsB->altitude  = summary->anitaLocation.altitude;
      gpsB->heading    = summary->anitaLocation.heading;
      UsefulAdu5Pat *usefulGPSB = new UsefulAdu5Pat(gpsB);
      
      //if event a is further than 1000km from location B, just skip it
      double distance = usefulGPSB->getDistanceFromSource(latA,lonA,altA);
      if (distance > 1000e3) continue;

      //where event b was seen (b)
      double thetaB = summary->peak[0][0].theta;
      double phiB   = summary->peak[0][0].phi;
      double latB = summary->peak[0][0].latitude;
      double lonB = summary->peak[0][0].longitude;
      double altB = summary->peak[0][0].altitude;      
      
      //where event b is seen from A's location
      double thetaBA,phiBA; 
      usefulGPSB->getThetaAndPhiWave(lonA,latA,altA,thetaBA,phiBA);
      thetaBA *= TMath::RadToDeg();
      phiBA *= TMath::RadToDeg();
      //difference between A->b and A->a
      double diffBAtheta = TMath::Abs(thetaBA - thetaA);
      double diffBAphi = FFTtools::wrap(TMath::Abs(phiBA - phiA));
      if (diffBAphi > 180) diffBAphi = 360 - diffBAphi;
      
      //where event a is seen from B's location
      double thetaAB,phiAB; 
      usefulGPSA->getThetaAndPhiWave(lonB,latB,altB,thetaAB,phiAB);
      thetaAB *= TMath::RadToDeg();
      phiAB   *= TMath::RadToDeg();
      //difference between B->a and B->b
      double diffABtheta = TMath::Abs(thetaAB - thetaB);
      double diffABphi = FFTtools::wrap(TMath::Abs(phiAB - phiB));
      if (diffABphi > 180) diffABphi = 360 - diffABphi;      


      const double sigmaTheta = 0.2193;
      const double sigmaPhi = 0.5429;
      
      double diffTheta = TMath::Sqrt(pow(diffABtheta,2) + pow(diffBAtheta,2))/sigmaTheta;
      double diffPhi = TMath::Sqrt(pow(diffABphi,2) + pow(diffBAphi,2))/sigmaPhi;
      
      double diff = TMath::Sqrt(pow(diffTheta,2) + pow(diffPhi,2));
      
      //	cout << entryA << " vs " << entryB << " = " << diff << endl;

      if (diff < 50) close++;
      
      hCluster->Fill(entryA,diff);
      
      delete gpsB;
      delete usefulGPSB;
      
      cout << entryA << "," << entryB << " | ";
      cout << phiAB << "," << phiB << " - " << phiBA << "," << phiA << " " << diffPhi << " | ";
      cout << thetaA << "," << thetaB << " - " << diffBAtheta << " " << diffABtheta << " " << diffTheta << " | ";
      cout << distance << endl;
      
    }

    delete gpsA;
    delete usefulGPSA;
    
    cout << "entryA:" << entryA << " eventNumber " << eventNum << " close:" << close << endl;

    if (close == 0) {
      cout << "Close! " << singlets->GetN() << " so far " << endl;
      singlets->SetPoint(singlets->GetN(),singlets->GetN(),eventNum);
    }
  }
  
  hCluster->Draw("colz");

  
  cout << "Close Event Numbers:" << endl;
  for (int i=0; i<singlets->GetN(); i++) {
    cout << singlets->GetY()[i] << endl;
  }

  TFile* outFile = new TFile("cluster.root","recreate");
  hCluster->Write();
  singlets->Write();
  outFile->Close();
  
  return;
}
