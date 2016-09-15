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
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAnitaEvent.h"
#include "Adu5Pat.h"
#include "AntarcticaMapPlotter.h"
#include "CrossCorrelator.h"
#include "UsefulAdu5Pat.h"
#include "AnitaConventions.h"

using namespace std;

/*
  Ben Rotter September 2016

  Ben Strutt and Cosmin created a nice correlator, so how about I try to use that and see how it compares to
  Abby's.  Abby doesn't really reconstruct WAIS pulses very well, so either it is a filtering issue on her end,
  or there is something much more pathalogical with her code?  Either way, if this works I can integrate it in
  to her code and re-run all the events with a more "modern" correlator

  This takes FOREVER to run.  So what I want to do is set it up so it can do it run by run on the big servers.
  That way I can get through all the events without having to wait as long at least.

  So lets do:
  argc=1:
   argv[0] = scriptName (like it always will be)
             otherwise will process all events in runs 330 to 360
  argc=2:
   argv[0] = same as above
   argv[1] = single runNumber to work with (all events in that run)
  argc=3:
   argv[0] = save as above
   argv[1] = runNumber to START at (all events in that run)
   argv[2] = runNumber to STOP at (all events in that run too)
  argc=4:
   argv[0:2] = same as above
   argv[3] =  number of entries to look at (over several runs!)
  argc=5:
   argv[0:2] = same as above
   argv[3] = TChain entry number to start with
   argv[4] = TChain entry number to stop with
  argc=6:
   argv[0:4] = same as above
   argv[5] = special post-fix to add to output file name (since I have to batch run this)

 */


int main(int argc, char** argv) {

  cout << "Hello!  Let us do some physics mate" << endl;

  int startRun=330;
  int stopRun =360;
  int startWaisEntry= 0;
  int stopWaisEntry = -1;
  string postFix = "";

  if (argc == 1) {
    cout << "Using all events in runs 330 through 360" << endl;
  }
  else if (argc == 2) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[1])+1; 
    cout << "Using only run " << startRun << endl;
  }
  
  else if (argc == 3) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    cout << "Using runs " << startRun << " through " << stopRun << endl;
  }
  else if (argc == 4) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    stopWaisEntry = atoi(argv[3]);
    cout << "Using first " << stopWaisEntry << " in runs " << startRun << " through " << stopRun << endl;
  }
  else if (argc == 5) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    startWaisEntry = atoi(argv[3]);
    stopWaisEntry = atoi(argv[4]);
    cout << "Using event " << startWaisEntry << " through " << stopWaisEntry << " in runs " << startRun << " through " << stopRun << endl;
  }
  else if (argc == 6) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    startWaisEntry = atoi(argv[3]);
    stopWaisEntry = atoi(argv[4]);
    postFix.assign(argv[5]);
    postFix.insert(0,"_");
    cout << "Using event " << startWaisEntry << " through " << stopWaisEntry << " in runs " << startRun << " through " << stopRun << " and naming the file with " << postFix << " at the end" << endl;
  }
  else {
    cout << "Somehow you didn't input the correct number of parameters, which is hard" << endl;
  }



  stringstream name;
  //okay lets start by grabbing the wais header files
  TChain *waisHeadTree = new TChain("headTree","headTree");  
  name.str("");
  name << "waisHeadFile.root";
  waisHeadTree->Add(name.str().c_str());
  RawAnitaHeader *waisHead = NULL;
  waisHeadTree->SetBranchAddress("header",&waisHead);

  int numWaisEntries = waisHeadTree->GetEntries();
  cout << "I found " << numWaisEntries << " wais pulser entries using file:" << endl;
  cout << name.str() << endl;
    

  //then we can grab all the event files for the WAIS range (runs 330 to 360 gets them all, though is wide)
  //Also grab the gps stuff for pointing
  TChain *eventTree = new TChain("eventTree","eventTree");  
  TChain *headTree = new TChain("headTree","headTree");
  TChain *gpsTree = new TChain("adu5PatTree","adu5PatTree");
  
  char* dataDir = getenv("ANITA3_DATA");
  for (int i=startRun; i<stopRun; i++) {
    name.str("");
    name << dataDir << "run" << i << "/calEventFile" << i << ".root";
    eventTree->Add(name.str().c_str());
    name.str("");
    name << dataDir << "run" << i << "/headFile" << i << ".root";
    headTree->Add(name.str().c_str());
    name.str("");
    name << dataDir << "run" << i << "/gpsEvent" << i << ".root";
    gpsTree->Add(name.str().c_str());
  }
  int numEventEntries = eventTree->GetEntries();
  //  RawAnitaEvent *event = NULL;
  CalibratedAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);
  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);
  Adu5Pat *gps = NULL;
  gpsTree->SetBranchAddress("pat",&gps);

  
  //Figure out how many events we are dealing with
  cout << "There are " << numEventEntries << " events imported too." << endl;
  cout << "There are " << gpsTree->GetEntries() << " gps events imported too." << endl;

  //I want to sort these by eventNumber, since most of the events aren't WAIS pulses
  cout << "Building event index (this might take awhile)..." << endl;
  headTree->BuildIndex("eventNumber");
  waisHeadTree->BuildIndex("eventNumber");
  cout << "Event index built" << endl;


  //create the correlator object
  CrossCorrelator *correlator = new CrossCorrelator();

  //notch it where I know the sats are I guess?
  CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite Notch", 260 - 26, 260 + 26);
  correlator->addNotch(notch260);

  CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch", 370 - 26, 370 + 26);
  correlator->addNotch(notch370);

  //also we are working with horizontal polarization, so this is for ease
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  //open the output file before making the tree
  name.str("");
  name << "rootFiles/waisNewCorrelator" << postFix << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");

  //save all the output data in a TTree to plot later
  TTree *newCorrelatorTree = new TTree("newCorrelatorTree","newCorrelatorTree");
  Double_t peakValue,peakPhiDeg,peakThetaDeg,lat,lon,alt,theta_adjustment_required,heading,waisTheta,waisPhi;
  Int_t eventNumber,returnValue;
  newCorrelatorTree->Branch("peakValue",&peakValue);
  newCorrelatorTree->Branch("peakPhiDeg",&peakPhiDeg);
  newCorrelatorTree->Branch("peakThetaDeg",&peakThetaDeg);
  newCorrelatorTree->Branch("lat",&lat);
  newCorrelatorTree->Branch("lon",&lon);
  newCorrelatorTree->Branch("alt",&alt);
  newCorrelatorTree->Branch("theta_adjustment_required",&theta_adjustment_required);
  newCorrelatorTree->Branch("eventNumber",&eventNumber);
  newCorrelatorTree->Branch("heading",&heading);
  newCorrelatorTree->Branch("returnValue",&returnValue);
  newCorrelatorTree->Branch("waisTheta",&waisTheta);
  newCorrelatorTree->Branch("waisPhi",&waisPhi);

  waisHeadTree->GetEntry(startWaisEntry);
  int firstWaisEventNumber = waisHead->eventNumber;
  int startEntry = headTree->GetEntryNumberWithIndex(firstWaisEventNumber);
  waisHeadTree->GetEntry(stopWaisEntry);
  int lastWaisEventNumber = waisHead->eventNumber;
  int stopEntry = headTree->GetEntryNumberWithIndex(lastWaisEventNumber);

  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10 == 0) {
      cout << entry << " / " << stopEntry-startEntry << "\r";
      fflush(stdout);
    }

    eventTree->GetEntry(entry);
    headTree->GetEntry(entry);
    eventNumber = head->eventNumber;
    //calibrating requires doing ALL the events IN ORDER, so I need to do this even though I don't use most
    //once I generate all the CalibratedAnitaEvent.root files I don't have to do this
    //    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(event,WaveCalType::kFull,head);
    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(event);

    //if the eventNumber isn't in the waisTree, just move on
    int waisEntry = waisHeadTree->GetEntryNumberWithIndex(eventNumber);
    if (waisEntry==-1) {
      delete usefulEvent;
      continue;
    }

    //Otherwise do what you were doing before!
    waisHeadTree->GetEntry(waisEntry);
    gpsTree->GetEntry(entry);

    correlator->reconstructEvent(usefulEvent,1,1);
    
    //Get the actual map (I don't need to do this! it only returns the course map and I want the fine map)
    //    TH2D *mapHist = correlator->getMap(AnitaPol::kHorizontal,peakValue,peakPhiDeg,peakThetaDeg);
    //    delete mapHist;

    //get the peak phi bin
    peakValue = correlator->fineMapPeakValues[0][0];
    peakPhiDeg = correlator->fineMapPeakPhiDegs[0][0];
    peakThetaDeg = correlator->fineMapPeakThetaDegs[0][0];

    heading = gps->heading;
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    returnValue = usefulGPS->traceBackToContinent(peakPhiDeg*TMath::DegToRad(),
						  -1*peakThetaDeg*TMath::DegToRad(),
						  &lat,&lon,&alt,&theta_adjustment_required);
    usefulGPS->getThetaAndPhiWaveWaisDivide(waisTheta,waisPhi);
    delete usefulGPS;

    //these are in radians!
    waisTheta *= TMath::RadToDeg();
    waisPhi *= TMath::RadToDeg();

    //try to get them to wrap correctly?
    if (waisTheta < -180) waisTheta += 360;
    if (waisPhi < -180) waisPhi += 360;
    if (peakThetaDeg < -180) peakThetaDeg += 360;
    if (peakPhiDeg < -180) peakPhiDeg += 360;

    if (waisTheta > 180) waisTheta -= 360;
    if (waisPhi > 180) waisPhi -= 360;
    if (peakThetaDeg > 180) peakThetaDeg -= 360;
    if (peakPhiDeg > 180) peakPhiDeg -= 360;

    outFile->cd();
    newCorrelatorTree->Fill();
    
  }

  cout << "stored " << newCorrelatorTree->GetEntries() << " entries" << endl;
  outFile->cd();
  newCorrelatorTree->Write();

  outFile->Close();

  //Thats all folks!
  cout << "Okay I quit because I did everything you told me, goodbye!" << endl;


  return 1;
}
  

