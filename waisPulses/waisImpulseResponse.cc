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


  FFTtools::loadWisdom("/Users/brotter/macros/fftWisdom.dat");


  //start stuff
  cout << "Hello!  Let us do some physics mate" << endl;

  int startEntry = 0;
  int stopEntry = -1;

  if (argc==1) {
    cout << "Using all WAIS Pulser events" << endl;
  }
  else if (argc==2) {
    startEntry = 0;
    stopEntry = atoi(argv[1]);
    cout << "Using first " << stopEntry << " WAIS Pulse events" << endl;
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
    cout << "Who knows what the hell you are doing, figure it out son" << endl;
    return -1;
  }




  stringstream name;
  //okay lets start by grabbing the wais header files
  TChain *waisHeadTree = new TChain("headTree","headTree");  
  name.str("");
  name << "/Users/brotter/Science/ANITA/ANITA3/anita16/rootFilesLocal/waisHeadFile.root";
  waisHeadTree->Add(name.str().c_str());
  RawAnitaHeader *waisHead = NULL;
  waisHeadTree->SetBranchAddress("header",&waisHead);

  TChain *calEventTree = new TChain("eventTree","eventTree");
  name.str("");
  name << "/Users/brotter/Science/ANITA/ANITA3/anita16/rootFilesLocal/calEventFileWais.root";
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


  //anita geomtool is key too
  AnitaGeomTool *geom = AnitaGeomTool::Instance(3);



  const int numChans= 12*9;

  //array of TGraphs to store all the impulse responses?
  TGraph *impulseResponsesTrig[numChans];
  //lets do it as a function of angle from boresight too! From -60 to +60, in steps of 2 (-1 to 1, 1 to 3, 3 to 5, etc)
  const double startAngle = -60;
  const double stopAngle = 60;
  const int numAngles = 60;

  const double angleStep = (stopAngle - startAngle + 1)/numAngles;
  TGraph *impulseResponsesAng[numChans][numAngles];
  
  //Lets get an impulse response template for things to correlate to...
  //  TGraph *templateResponse = new TGraph("~/Science/ANITA/ANITA3/benCode/analysis/impulseResponse/integratedTF/autoPlots/A3ImpulseResponse/09BV.txt");
  //and make that the default for all the graphs (I'll subtract it at the end I guess?)


  
  //lets make a counter for each antenna too
  int counterTrig[numChans];
  int counterAng[numChans][numAngles];
  for (int chanIndex=0; chanIndex<numChans; chanIndex++) {
    counterTrig[chanIndex] = 0;
    for (int angle=0; angle<numAngles; angle++) {
      counterAng[chanIndex][angle] = 0;
    }
  }


  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%1000 == 0) {
      cout << entry << " / " << numWaisEntries << "\r";
      fflush(stdout);
    }

    //get the wais header and its event number
    waisHeadTree->GetEntry(entry);
    calEventTree->GetEntry(entry);

    //and gps
    int eventNumber = waisHead->eventNumber;
    int gpsEntry = gpsTree->GetEntryNumberWithIndex(eventNumber);
    gpsTree->GetEntry(gpsEntry);

    //get the calibrated event
    UsefulAnitaEvent *usefulEvent = new UsefulAnitaEvent(calEvent);

    //and gps
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
    //direction to wais (phi2 == 0deg)
    double waisTheta,waisPhi;
    usefulGPS->getThetaAndPhiWaveWaisDivide(waisTheta,waisPhi);
    //its in radians ugh
    waisPhi *= (180. / TMath::Pi());


    //define this pointer structure here
    TGraph *grToCorrelate[2];

    for (int phi=0; phi<16; phi++) {

      //is this a triggered phi sector?
      int trig = waisHead->isInL3Pattern(phi,AnitaPol::kHorizontal);

      //how far away from that phi sector is it?
      double angDiff = waisPhi - (phi-2)*22.5;
      if (angDiff > 180) angDiff -= 360;
      if (angDiff <- 180) angDiff += 360;
      //      if (phi < 2) cout << "phi: " << phi << " waisPhi: " << waisPhi << " angDiff: " << angDiff << endl;

      int angBin = -1;
      for (int bin=0; bin<numAngles; bin++) {
	double binLow = startAngle + bin*angleStep;
	double binHigh = startAngle + (bin+1)*angleStep;
	if ((angDiff > binLow) && (angDiff < binHigh)) {
	  angBin = bin;
	  break;
	}
      }


      //loop through the rings in that phi sector (I can't believe this is how you have to do it)
      for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
	AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
	
	int chanIndex = geom->getChanIndexFromRingPhiPol(ring,phi,AnitaPol::kHorizontal);


	//I want to get the interpolated graph (impulse response is at 0.1ps)
	TGraph *currRawGraph = usefulEvent->getGraph(ring,phi,AnitaPol::kHorizontal);
	TGraph *currGraph = FFTtools::getInterpolatedGraph(currRawGraph,1/2.6);
	delete currRawGraph;
	TGraph *currGraphInterp = FFTtools::getInterpolatedGraphFreqDom(currGraph,0.1);
	delete currGraph;
	//Also, we have to pad it because sometimes it starts nearer to the end
	TGraph *goodGraph = FFTtools::padWave(currGraphInterp,3);
	delete currGraphInterp;


	//Do a correlation with the angle bin that corresponds to in that phi sector
	if (counterAng[chanIndex][angBin] == 0 && angBin != -1) {
	  impulseResponsesAng[chanIndex][angBin] = (TGraph*)goodGraph->Clone();
	  cout << "new angle!" << chanIndex << " " << angBin << endl;
	  counterAng[chanIndex][angBin]++;
	}
	else if (angBin != -1) {
	  //Do it for whatever angle bin that corresponds to
	  grToCorrelate[0] = impulseResponsesAng[chanIndex][angBin];
	  grToCorrelate[1] = goodGraph;
	
	  //then average them and store it in the big array of things I'm lugging around
	  impulseResponsesAng[chanIndex][angBin] = FFTtools::correlateAndAverage(2,grToCorrelate);
	
	  //correlateAndAverage divides the end result by the first argument, so lets undo that to get the sum of the two
	  for (int i=0; i<impulseResponsesAng[chanIndex][angBin]->GetN(); i++) {
	    impulseResponsesAng[chanIndex][angBin]->GetY()[i] *= 2;
	  }
	
	  //and delete the old array (helpfully stored as this pointer)
	  delete grToCorrelate[0];
     
	  counterAng[chanIndex][angBin]++;
	}

	//Do it for Trigger 
	//if this is the first one, then just store it
	if (trig) {
	  if (counterTrig[chanIndex] == 0) {
	    impulseResponsesTrig[chanIndex] = (TGraph*)goodGraph->Clone();
	    cout << "new trig! " << chanIndex << endl;
	    counterTrig[chanIndex]++;
	  }
	  //otherwise correlate and average it
	  else {
	    //correlateAndAverage needs it like this I think?
	    grToCorrelate[0] = impulseResponsesTrig[chanIndex];
	    grToCorrelate[1] = goodGraph;
	    
	    //then average them and store it in the big array of things I'm lugging around
	    impulseResponsesTrig[chanIndex] = FFTtools::correlateAndAverage(2,grToCorrelate);
	    
	    //correlateAndAverage divides the end result by the first argument, so lets undo that to get the sum of the two
	    for (int i=0; i<impulseResponsesTrig[chanIndex]->GetN(); i++) {
	      impulseResponsesTrig[chanIndex]->GetY()[i] *= 2;
	    }
	    
	    //and delete the old array (helpfully stored as this pointer)
	    delete grToCorrelate[0];
	    
	    //and incriment the counter for that channel
	    counterTrig[chanIndex]++;
	  }	    
	}//end else for correlating and averaging triggered phi sector
	
	//done with the interpolated waveform
	delete goodGraph;
	
	
      } //end ring loop

    }//end phi loop      
      
      
    //Got what I wanted from that event, done with these classes
    delete usefulGPS;
    delete usefulEvent;
    
    
  } //end entry loop
  

  //now that we have the coherent sum of them all, we have to divide by however many we saw
  for (int chanIndex=0; chanIndex<numChans; chanIndex++) {
    cout << chanIndex << " saw " << counterTrig[chanIndex] << " pulses | ";
    if (counterTrig[chanIndex]!=0) {
      for (int pt=0; pt<impulseResponsesTrig[chanIndex]->GetN(); pt++){
	impulseResponsesTrig[chanIndex]->GetY()[pt] /= counterTrig[chanIndex];
      }
    }
    for (int angBin=0; angBin<numAngles; angBin++) {
      cout << counterAng[chanIndex][angBin] << " ";
      if (counterAng[chanIndex][angBin]!=0) {
	for (int pt=0; pt<impulseResponsesAng[chanIndex][angBin]->GetN(); pt++) {
	  impulseResponsesAng[chanIndex][angBin]->GetY()[pt] /= counterAng[chanIndex][angBin];
	}
      }
    }
  cout << endl;
  }


  
  //get output file
  name.str("");
  if (argc!=4) {
    name << "waisImpulseResponse_wAngs.root";
  }
  else {
    name << argv[3];
  }
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  
  
  int surf,chan,phi,ant;
  AnitaRing::AnitaRing_t ring;
  AnitaPol::AnitaPol_t pol;
  for (int chanIndex=0; chanIndex<numChans; chanIndex++) {
    if (counterTrig[chanIndex] != 0 ) {
      geom->getSurfChanFromChanIndex(chanIndex,surf,chan);
      geom->getRingAntPolPhiFromSurfChan(surf,chan,ring,ant,pol,phi);
      name.str("");
      name << "wais";
      if (phi<9) name << "0";
      name << phi+1 << AnitaRing::ringAsChar(ring) << AnitaPol::polAsChar(pol);
      impulseResponsesTrig[chanIndex]->SetName(name.str().c_str());
      impulseResponsesTrig[chanIndex]->SetTitle(name.str().c_str());
      impulseResponsesTrig[chanIndex]->Write();
      delete impulseResponsesTrig[chanIndex];
    }
    for (int angBin=0; angBin<numAngles; angBin++) {
      if (counterAng[chanIndex][angBin] != 0) {
	name.str("");
	name << "wais";
	if (phi<9) name << "0";
	name << phi+1 << AnitaRing::ringAsChar(ring) << AnitaPol::polAsChar(pol) << "_" << angBin;
	impulseResponsesAng[chanIndex][angBin]->SetName(name.str().c_str());
	impulseResponsesAng[chanIndex][angBin]->SetTitle(name.str().c_str());
	impulseResponsesAng[chanIndex][angBin]->Write();
	delete impulseResponsesAng[chanIndex][angBin];
      }
    }
  }
  outFile->Close(); 
  
  cout << "trying to save wisdom" << endl;
  FFTtools::saveWisdom("/Users/brotter/macros/fftWisdom.dat");
  
  
  cout << "Done!  Just have to return now" << endl;
  return 1;
}
  

