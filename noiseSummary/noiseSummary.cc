/*

  template search has become a monster.

  I want to just generate the noise summary simply without all the insane things I've tacked on

  Also it might help me simplify all my nonsense code

 */

#include <sstream>
#include <iostream>

#include "AnitaDataset.h"
#include "AnitaEventSummary.h"
#include "AnitaNoiseSummary.h"
#include "AnitaNoiseMachine.h"
#include "AnalysisConfig.h"
#include "Analyzer.h"
#include "FilterStrategy.h"
#include "FilteredAnitaEvent.h"
#include "SpectrumAverage.h"
#include "UCFilters.h"
#include "RawAnitaHeader.h"

#include "TFile.h"
#include "TTree.h"


using namespace std;

int main(int argc,char** argv) {

  int run,numEntries;
  string outFileName;
  stringstream name;

  if (argc == 3) {
    run = atoi(argv[1]);
    outFileName = argv[2];
    numEntries = -1;
  }
  else if (argc == 4) {
    run = atoi(argv[1]);
    outFileName = argv[2];
    numEntries = atoi(argv[3]);
  }
  else {
    cout << "Usage: " << argv[0] << " [run number] [output file base name]" << endl;
    return -1;
  }

  cout << " Hello! :D" << endl;

  //open data
  AnitaDataset *data = new AnitaDataset(run);

  //determine number of events
  int lenEntries = data->N();
  if (numEntries == -1 || numEntries > lenEntries) {
    numEntries = lenEntries;
  }
  cout << "Found " << lenEntries << " events in run " << run << ".  Using " << numEntries << " of them." << endl;


  //open output file
  name.str(""); name << outFileName << ".root";
  cout << "Making output file " << name.str() << " and output tree and filling it with stuff" << endl;
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  //make output tree
  TTree* summaryTree = new TTree("summaryTree","summaryTree");
  //make output things
  AnitaEventSummary *eventSummary = new AnitaEventSummary();
  summaryTree->Branch("eventSummary",&eventSummary);

  AnitaNoiseSummary *noiseSummary = new AnitaNoiseSummary();
  summaryTree->Branch("noiseSummary",&noiseSummary);

  //make the noise analyzer
  AnitaNoiseMachine *noiseMachine = new AnitaNoiseMachine();
  noiseMachine->fillMap = true;
  noiseMachine->quiet = false;
  
  //make the analyzer  
  cout << "Making the Analyzer" << endl;
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true);

  
  //make a filtering plan
  FilterStrategy *strategy = new FilterStrategy();
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);
  UCorrelator::AdaptiveBrickWallFilter *brickWall = new UCorrelator::AdaptiveBrickWallFilter(specAvgLoader,2,false);
  strategy->addOperation(brickWall);

  int cnt=0;
  cout << "Okay lets go" << endl;
  for (int entry=0; entry<numEntries; entry++) {
    data->getEntry(entry);

    //only want minbias stuff, skip it otherwise
    int trigType = data->header()->trigType&0x0F;
    if (trigType == 1) continue;
    cnt++;

    cout << entry << " / " << numEntries << " (" << cnt << ")" << endl;
      
   
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    analyzer->analyze(filteredEvent, eventSummary); 

    //do the noise analysis (only update it if it is a min bias though
    if (trigType == 1) {
      noiseSummary->isMinBias = false;
    }
    else {
      cout << "updating" << endl;
      noiseSummary->isMinBias = true;
      noiseMachine->updateMachine(analyzer,filteredEvent);
    }

    noiseMachine->fillNoiseSummary(noiseSummary);
    noiseMachine->fillEventSummary(eventSummary);
    summaryTree->Fill();

    delete filteredEvent;
    analyzer->clearInteractiveMemory();
    
  }

  outFile->cd();
  summaryTree->Write();
  outFile->Close();

  cout << "Done!  Thanks :) " << endl;


  return 1;
}
