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
#include "TGraph.h"
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
#include "FFTtools.h"

using namespace std;

/*
  Ben Rotter October 2016

  So now that I've got a good transfer function for each antenna and channel, it would be nice
  to compare them vs the WAIS pulses.

  so, for now, this code will take each wais pulse, find the phi that points at wais, and averages those
  waveforms together to find the impulse response.

  no arguments: do all the events
  two arguments: first=startWaisEntry second=lastWaisEntry
  three arguments: first=startWaisEntry second=lastWaisEntry third=fileName

 */



int main(int argc, char** argv) {

  cout << "Hello!  Let us do some physics mate" << endl;

  int startEntry = 0;
  int stopEntry = -1;

  if (argc==1) {
    cout << "Using all WAIS Pulser events" << endl;
  }
  else if (argc==3) {
    startEntry = atoi(argv[1]);
    stopEntry = atoi(argv[2]);
    cout << "Using subset of WAIS Pulse events (" << startEntry << " to " << stopEntry << ")"  << endl;
  }
  else if (argc==4) {
    startEntry = atoi(argv[1]);
    stopEntry = atoi(argv[2]);
    cout << "Using subset of WAIS Pulse events (" << startEntry << " to " << stopEntry << ")  with output file name " << argv[3] << endl;
  }
  else {
    cout << "Who knows what the hell you are doing" << endl;
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

  if (stopEntry == -1 || stopEntry > numWaisEntries) {
    stopEntry = numWaisEntries;
  }
    

  //I made calibrated event files so this is a good place to use them!
  TChain *eventTree = new TChain("eventTree","eventTree");  
  TChain *gpsTree = new TChain("adu5PatTree","adu5PatTree");
  

  int startRun = 330;
  int stopRun = 360;

  char* dataDir = getenv("ANITA3_DATA");
  for (int i=startRun; i<stopRun; i++) {
    name.str("");
    name << dataDir << "run" << i << "/calEventFile" << i << ".root";
    eventTree->Add(name.str().c_str());
    name.str("");
    name << dataDir << "run" << i << "/gpsEvent" << i << ".root";
    gpsTree->Add(name.str().c_str());
  }
  int numEventEntries = eventTree->GetEntries();
  //  RawAnitaEvent *event = NULL;
  CalibratedAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);
  Adu5Pat *gps = NULL;
  gpsTree->SetBranchAddress("pat",&gps);

  
  //Figure out how many events we are dealing with
  cout << "There are " << numEventEntries << " events imported too." << endl;
  cout << "There are " << gpsTree->GetEntries() << " gps events imported too." << endl;

  //I want to sort these by eventNumber, since most of the events aren't WAIS pulses
  cout << "Building gps/event index (this might take awhile)..." << endl;
  gpsTree->BuildIndex("eventNumber");
  eventTree->BuildIndex("eventNumber");
  cout << "GPS/event index built" << endl;


  //array of TGraphs to store all the impulse responses?
  TGraph *impulseResponses[NUM_PHI];


  double waisTheta,waisPhi;


  //To check if I've been in that phi sector already (which is stupid to initialize for some reason)
  bool firstGraph[16];
  for (int phii=0; phii<16; phii++) {
    firstGraph[phii] = true;
    cout << firstGraph[phii];
  }
  cout << endl;

  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10 == 0) {
      cout << entry << " / " << numWaisEntries << "\r";
      fflush(stdout);
    }

    //get the wais header and its event number
    waisHeadTree->GetEntry(entry);
    int eventNumber = waisHead->eventNumber;
    int eventEntry = eventTree->GetEntryNumberWithIndex(eventNumber);

    //find the actual events with that info
    eventTree->GetEntry(eventEntry);
    gpsTree->GetEntry(eventEntry);

    //get the calibrated event
    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(event);

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    usefulGPS->getThetaAndPhiWaveWaisDivide(waisTheta,waisPhi);

    int phi = int((waisPhi*TMath::RadToDeg())/22.5);
    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      
      TGraph *currGraph = usefulEvent->getGraph(ring,phi,AnitaPol::kHorizontal);

      if (firstGraph[phi]) {
	cout << "FirstGraph! " << entry << " " << waisPhi << " " << phi << endl;
	impulseResponses[phi] = new TGraph(*currGraph);
	firstGraph[phi] = false;
      }
      else {
	TGraph *grToCorrelate[2];
	grToCorrelate[0] = new TGraph(*currGraph);
	grToCorrelate[1] = new TGraph(*impulseResponses[phi]);
	TGraph *correlated = FFTtools::correlateAndAverage(2,grToCorrelate);
	impulseResponses[phi] = new TGraph(*correlated);
	delete correlated;
	delete grToCorrelate[0];
	delete grToCorrelate[1];
      }


      delete currGraph;

    }

    //Got what I wanted, done with those classes
    delete usefulGPS;
    delete usefulEvent;
    

    
  }

  
  name.str("");
  if (argc!=4) {
    name << "waisImpulseResponse.root";
  }
  else {
    name << argv[3];
  }
      


  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  for (int phi=0; phi<16; phi++) {
    name.str("");
    name << "phi" << phi;
    impulseResponses[phi]->SetName(name.str().c_str());
    impulseResponses[phi]->Write();
  }
  outFile->Close();

  //Thats all folks!
  cout << "Okay I quit because I did everything you told me, goodbye!" << endl;

  for (int phi=0; phi<16; phi++) {
    delete impulseResponses[phi];
  }

  return 1;
}
  

