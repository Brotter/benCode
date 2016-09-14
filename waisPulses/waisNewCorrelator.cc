#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <random>
//root
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TF1.h"
//anita
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAnitaEvent.h"
#include "Adu5Pat.h"
#include "AntarcticaMapPlotter.h"
#include "CrossCorrelator.h"
#include "UsefulAdu5Pat.h"

using namespace std;

/*
  Ben Rotter September 2016

  Ben Strutt and Cosmin created a nice correlator, so how about I try to use that and see how it compares to
  Abby's.  Abby doesn't really reconstruct WAIS pulses very well, so either it is a filtering issue on her end,
  or there is something much more pathalogical with her code?  Either way, if this works I can integrate it in
  to her code and re-run all the events with a more "modern" correlator

 */


int main(int argc, char** argv) {

  cout << "Hello!  Let us do some physics mate" << endl;

  stringstream name;
  //okay lets start by grabbing the wais header files
  TChain *headTree = new TChain("headTree","headTree");  
  name.str("");
  name << "/Users/brotter/Science/ANITA/ANITA3/anita16/benPrograms/waisPulses/waisHeadFile.root";
  headTree->Add(name.str().c_str());
  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);

  int numEntries = headTree->GetEntries();
  cout << "I found " << numEntries << " wais pulser entries using file:" << endl;
  cout << name.str() << endl;
    

  //then we can grab all the event files for the WAIS range (runs 330 to 360 gets them all, though is wide)
  //Also grab the gps stuff for pointing
  TChain *eventTree = new TChain("eventTree","eventTree");  
  TChain *gpsTree = new TChain("adu5PatTree","adu5PatTree");
  for (int i=330; i<354; i++) {
    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << i << "/eventFile" << i << ".root";
    eventTree->Add(name.str().c_str());
    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << i << "/gpsEvent" << i << ".root";
    gpsTree->Add(name.str().c_str());
  }
  RawAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);
  Adu5Pat *gps = NULL;
  gpsTree->SetBranchAddress("pat",&gps);

  
  //Figure out how many events we are dealing with
  cout << "There are " << eventTree->GetEntries() << " events imported too." << endl;
  cout << "There are " << gpsTree->GetEntries() << " gps events imported too." << endl;

  //I want to sort these by eventNumber, since most of the events aren't WAIS pulses
  eventTree->BuildIndex("eventNumber");
  cout << "Event index built" << endl;


  //create the correlator object
  CrossCorrelator *correlator = new CrossCorrelator();

  //also we are working with horizontal polarization, so this is for ease
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  TH2D *pointingMap = new TH2D("pointingMap","pointingMap;lat;long",
					     200,-82,-78,200,-119,-105);

  for (int entry=0; entry<1000; entry++) {
    headTree->GetEntry(entry);
    int eventEntry = eventTree->GetEntryNumberWithBestIndex(head->eventNumber);
    eventTree->GetEntry(eventEntry);
    gpsTree->GetEntry(eventEntry);

    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(event,WaveCalType::kFull,head);
    correlator->reconstructEvent(usefulEvent);
    
    
    Double_t peakValue,peakPhiDeg,peakThetaDeg,lat,lon,alt,theta_adjustment_required;
    TH2D *mapHist = correlator->getMap(AnitaPol::kHorizontal,peakValue,peakPhiDeg,peakThetaDeg);
    delete mapHist;
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    int returnValue = usefulGPS->traceBackToContinent(peakPhiDeg*TMath::DegToRad(),
						      peakThetaDeg*TMath::DegToRad(),
						      &lat,&lon,&alt,&theta_adjustment_required);
    delete usefulGPS;

    cout << peakPhiDeg << " " << peakThetaDeg << " " << lat << " " << lon << " " << returnValue << endl;


    if (returnValue==1) {
      pointingMap->Fill(lat,lon); }
  }
  
  TFile *outFile = TFile::Open("waisNewCorrelator.root","recreate");

  pointingMap->Write();

  outFile->Close();

  //Thats all folks!
  cout << "Okay I quit because I did everything you told me, goodbye!" << endl;


  return 1;
}
  

