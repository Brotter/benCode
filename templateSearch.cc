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
#include "TH2.h"
#include "TH2D.h"
#include "TEllipse.h" 
#include "TMarker.h" 
#include "TStyle.h" 
#include "TCanvas.h"
//anita
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "Adu5Pat.h"
#include "CalibratedAnitaEvent.h"
#include "UsefulAnitaEvent.h"

//cosmin's stuff
#include "AnitaEventSummary.h"
#include "AnalysisConfig.h" 
#include "UCFilters.h" 
#include "PeakFinder.h" 
#include "FilterStrategy.h" 
#include "Correlator.h" 
#include "Analyzer.h" 
#include "WaveformCombiner.h"
#include "AnitaDataset.h"
#include "BlindDataset.h"

using namespace std;

/*

  Ben Rotter - March 2017 - University of Hawaii at Manoa

  FINAL STRETCH

  I have a solid impulse response for the system now, and a good template for CR events.  So lets search based on that

  Also switch over to the DataSet type format, which makes everything go real quickly.

 */


int main(int argc, char** argv) {


  int runNum;
  string outFileName;
  int lenEntries = -1;

  if (argc==3) {
    runNum = atoi(argv[1]);
    outFileName = argv[2];
    cout << "Hello!  Let us do some physics mate" << endl;
    cout << "Using run " << runNum << " and outfile " << outFileName << endl;
  }
  else if (argc==4) {
    runNum = atoi(argv[1]);
    outFileName = argv[2];
    lenEntries = atoi(argv[3]);
    cout << "Hello!  Let us do some physics mate" << endl;
    cout << "Using run " << runNum << " and outfile " << outFileName << endl;
    cout << "Only doing " << lenEntries << " of the first entries" << endl;
  }
  else {
    cout << "Usage: " << argv[0] << " [run] [output base filename] [opt: num entries]" << endl;
    return -1;
  }


  //  AnitaDataset *data = new AnitaDataset(runNum);
  BlindDataset *data = new BlindDataset(runNum,true);

  int numEntries = data->N();
  cout << "number of entries in run:" << numEntries << endl;


  stringstream name;
  //also make a Tree of event headers that pass cuts
  name.str("");
  name << outFileName << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  outFile->cd();
  TTree *outTree = new TTree("summaryTree","summaryTree");

  //Lets make the summary object that I can shove into the output tree
  AnitaEventSummary *eventSummary = new AnitaEventSummary; 


  outTree->Branch("eventSummary",&eventSummary);

  //Make a filter strategy
  //  with a debug file
  //  name.str("");
  //  name << outFileName << "_filtOutFile.root";
  //  TFile *filterOutFile = TFile::Open(name.str().c_str(),"recreate"); 
  //  FilterStrategy strategy(filterOutFile);
  //  without a debug file
  FilterStrategy *strategy = new FilterStrategy();
  //Add the actual Filters
  //  with the sine subtract alghorithm
  FilterOperation *sineSub = new UCorrelator::SineSubtractFilter();

  strategy->addOperation(sineSub);
  //  with abby's list of filtering
  //  UCorrelator::applyAbbysFilterStrategy(&strategy);



  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  //set the response to my "single" response
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  config.response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;

  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); ;
  
  
  //and get the "averaged" impulse response as the template"
  char* templateDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << templateDir << "/share/UCorrelator/responses/SingleBRotter/all.imp";
  TGraph *grTemplate = new TGraph(name.str().c_str());

  //and add a thing to store whatever it finds
  Double_t templateValue = 0;
  outTree->Branch("templateValue",&templateValue);


  //**loop through entries
  //option to have less entries! (-1 is the default, so in case you don't specify)
  if (lenEntries == -1) lenEntries = numEntries;

  int entryIndex=0;

  cout << "templateSearch(): starting event loop" << endl;

  for (int entry=0; entry<lenEntries; entry++) {
    //    if (entry%1==0) {
    cout << entry << "/" << lenEntries << "(" << entryIndex << ")" << endl;
    fflush(stdout);
    //    }
    //get all the pointers set right
    data->getEntry(entry);
    
    //1) calibrate and then filter the event and get a FilteredAnitaEvent back
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    cout << filteredEvent << endl;
    analyzer->analyze(filteredEvent, eventSummary); 
    delete filteredEvent;

    //Lets figure out which was the trigger (H=0, V=1, also defaults to H)
    AnitaPol::AnitaPol_t whichTrig =  (AnitaPol::AnitaPol_t)eventSummary->flags.isVPolTrigger;

    const TGraphAligned *coherentAligned = analyzer->getCoherent(whichTrig,0)->even();
    TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
    TGraph *grCorr = FFTtools::getCorrelationGraph(grTemplate,coherent);
    templateValue = TMath::Abs(TMath::MaxElement(grCorr->GetN(),grCorr->GetY()));

    delete coherent;
    delete grCorr;


    outFile->cd();
    outTree->Fill();

    analyzer->clearInteractiveMemory();
  }


  outFile->cd();
  outTree->Write();
  outFile->Close();

  delete eventSummary;

  //  cout << "Physics complete!  See ya later buddy :)" << endl;

  return 1;
  
}


  
