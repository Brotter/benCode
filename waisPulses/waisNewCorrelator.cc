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
  int startEntry= 0;
  int stopEntry = -1;
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
    stopEntry = atoi(argv[3]);
    cout << "Using first " << stopEntry << " in runs " << startRun << " through " << stopRun << endl;
  }
  else if (argc == 5) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    startEntry = atoi(argv[3]);
    stopEntry = atoi(argv[4]);
    cout << "Using event " << startEntry << " through " << stopEntry << " in runs " << startRun << " through " << stopRun << endl;
  }
  else if (argc == 6) {
    startRun = atoi(argv[1]);
    stopRun  = atoi(argv[2]); 
    startEntry = atoi(argv[3]);
    stopEntry = atoi(argv[4]);
    postFix.assign(argv[5]);
    postFix.insert(0,"_");
    cout << "Using event " << startEntry << " through " << stopEntry << " in runs " << startRun << " through " << stopRun << " and naming the file with " << postFix << " at the end" << endl;
  }
  else {
    cout << "Somehow you didn't input the correct number of parameters, which is hard" << endl;
  }



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
  for (int i=startRun; i<stopRun; i++) {
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
  cout << "Building event index (this might take awhile)..." << endl;
  eventTree->BuildIndex("eventNumber");
  cout << "Event index built" << endl;


  //create the correlator object
  CrossCorrelator *correlator = new CrossCorrelator();

  //also we are working with horizontal polarization, so this is for ease
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  //open the output file before making the tree
  name.str("");
  name << "waisNewCorrelator" << postFix << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");

  //save all the output data in a TTree to plot later
  TTree *newCorrelatorTree = new TTree("newCorrelatorTree","newCorrelatorTree");
  Double_t peakValue,peakPhiDeg,peakThetaDeg,lat,lon,alt,theta_adjustment_required,heading;
  Int_t eventNumber;
  newCorrelatorTree->Branch("peakValue",&peakValue);
  newCorrelatorTree->Branch("peakPhiDeg",&peakPhiDeg);
  newCorrelatorTree->Branch("peakThetaDeg",&peakThetaDeg);
  newCorrelatorTree->Branch("lat",&lat);
  newCorrelatorTree->Branch("lon",&lon);
  newCorrelatorTree->Branch("alt",&alt);
  newCorrelatorTree->Branch("theta_adjustment_required",&theta_adjustment_required);
  newCorrelatorTree->Branch("eventNumber",&eventNumber);
  newCorrelatorTree->Branch("heading",&heading);


  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10 == 0) {
      cout << entry << " / " << stopEntry-startEntry << "/r";
      fflush(stdout);
    }

    headTree->GetEntry(entry);
    if ((head->run < startRun) || (head->run > stopRun)) {
      continue;
    }
    eventNumber = head->eventNumber;
    int eventEntry = eventTree->GetEntryNumberWithBestIndex(eventNumber);
    eventTree->GetEntry(eventEntry);
    gpsTree->GetEntry(eventEntry);

    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(event,WaveCalType::kFull,head);
    correlator->reconstructEvent(usefulEvent);
    
    TH2D *mapHist = correlator->getMap(AnitaPol::kHorizontal,peakValue,peakPhiDeg,peakThetaDeg);
    delete mapHist;
    heading = gps->heading;
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    int returnValue = usefulGPS->traceBackToContinent(peakPhiDeg*TMath::DegToRad(),
						      peakThetaDeg*TMath::DegToRad(),
						      &lat,&lon,&alt,&theta_adjustment_required);
    delete usefulGPS;
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
  

