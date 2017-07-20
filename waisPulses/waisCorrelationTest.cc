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

  This "test" code saves all the correlated graphs because they AREN"T CORRELATING WELL ><

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

  TChain *calEventTree = new TChain("eventTree","eventTree");
  name.str("");
  name << "calEventFileWais.root";
  calEventTree->Add(name.str().c_str());
  CalibratedAnitaEvent *calEvent = NULL;
  calEventTree->SetBranchAddress("event",&calEvent);

  int numWaisEntries = waisHeadTree->GetEntries();
  cout << "I found " << numWaisEntries << " wais pulser entries using file:" << endl;
  cout << name.str() << endl;

  if (stopEntry == -1 || stopEntry > numWaisEntries) {
    stopEntry = numWaisEntries;
  }
    

  //I made calibrated event files so this is a good place to use them!
  //  TChain *eventTree = new TChain("eventTree","eventTree");  
  TChain *gpsTree = new TChain("adu5PatTree","adu5PatTree");
  

  int startRun = 330;
  int stopRun = 360;

  char* dataDir = getenv("ANITA3_DATA");
  for (int i=startRun; i<stopRun; i++) {
    //    name.str("");
    //    name << dataDir << "run" << i << "/calEventFile" << i << ".root";
    //    eventTree->Add(name.str().c_str());
    name.str("");
    name << dataDir << "run" << i << "/gpsEvent" << i << ".root";
    gpsTree->Add(name.str().c_str());
  }
  int numEventEntries = calEventTree->GetEntries();
  //  RawAnitaEvent *event = NULL;
  //  CalibratedAnitaEvent *event = NULL;
  //  eventTree->SetBranchAddress("event",&event);
  Adu5Pat *gps = NULL;
  gpsTree->SetBranchAddress("pat",&gps);

  
  //Figure out how many events we are dealing with
  cout << "There are " << numEventEntries << " events imported too." << endl;
  cout << "There are " << gpsTree->GetEntries() << " gps events imported too." << endl;

  if (numEventEntries==0 || gpsTree->GetEntries()==0 || numWaisEntries==0) {
    cout << "My input root data files are missing!  I didn't find anything to work with! Exiting..." << endl;
    return -1;
  }


  //I want to sort these by eventNumber, since most of the events aren't WAIS pulses
  cout << "Building gps/event index (this might take awhile)..." << endl;
  gpsTree->BuildIndex("eventNumber");
  //  eventTree->BuildIndex("eventNumber");
  cout << "GPS/event index built" << endl;


  //one TGraph to average things in (for phi0)
  TGraph *impulseResponse;

  
  //Lets get an impulse response template for things to correlate to...
  TGraph *templateResponse = new TGraph("~/Science/ANITA/ANITA3/benCode/analysis/impulseResponse/integratedTF/autoPlots/A3ImpulseResponse/09BV.txt");
  //and make that the default
  impulseResponse = new TGraph(*templateResponse);


  
  double waisTheta,waisPhi;

  TGraph *savedGraphs[48][5000];

  //lets make a counter for the phi sector
  int counter[48];
  for (int i=0; i<48; i++) {
    counter[i]=0;
  }

  //make a "test" output file
  TFile *outFile = TFile::Open("waisCorrelationTest.root","recreate");
    
  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10 == 0) {
      cout << entry << " / " << numWaisEntries << "\r";
      fflush(stdout);
    }

    //get the wais header and its event number
    waisHeadTree->GetEntry(entry);
    calEventTree->GetEntry(entry);
    int eventNumber = waisHead->eventNumber;
    //    int eventEntry = eventTree->GetEntryNumberWithIndex(eventNumber);
    int gpsEntry = gpsTree->GetEntryNumberWithIndex(eventNumber);

    //find the actual events with that info
    //    eventTree->GetEntry(eventEntry);
    gpsTree->GetEntry(gpsEntry);

    //get the calibrated event
    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(calEvent);

    //and the gps
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    usefulGPS->getThetaAndPhiWaveWaisDivide(waisTheta,waisPhi);
    delete usefulGPS;
    
    //which phi sector is it in?  well the gps "front" is phi sector 2 so maybe like this
    int phi = int((waisPhi*TMath::RadToDeg())/22.5)+2;
    if (phi>=16) phi -= 16;    

    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      int index = phi*3 + ringi;
      //I want to get the interpolated graph (impulse response is at 0.1ps)
      TGraph *currRawGraph = usefulEvent->getGraph(ring,phi,AnitaPol::kHorizontal);
      savedGraphs[index][counter[index]] = FFTtools::getInterpolatedGraph(currRawGraph,0.1);
      delete currRawGraph;
	
	counter[index]++;
	
    } //end ring loop
    
    
    
    //Got what I wanted, done with those classes
    delete usefulEvent;
    
    
  } //end entry loop
  
  
  for (int phi=0; phi<16; phi++) {
    for (int ring=0; ring<3; ring++) {
      int index = phi*3 + ring;
      TGraph *correlated = FFTtools::correlateAndAverage(counter[index],savedGraphs[index]);
      name.str("");
      name << "Phi" << phi << "Ring" << ring;
      correlated->SetName(name.str().c_str());
      outFile->cd();
      correlated->Write();
      delete correlated;
    }
  }


  outFile->Close();

  //Thats all folks!
  cout << "Okay I quit because I did everything you told me, goodbye!" << endl;
  for (int index=0; index<48; index++) {
    for (int i=0; i<counter[index]; i++) {
      delete savedGraphs[index][i];
    }
  }

  return 1;
}
  

