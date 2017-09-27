#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <random>
#include <iomanip>
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
#include "TStopwatch.h"
#include "Compression.h"
//anita
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "Adu5Pat.h"
#include "CalibratedAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "FFTtools.h"

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
#include "SpectrumAverage.h"
#include "ProgressBar.h"
#include "FilteredAnitaEvent.h"
#include "SystemResponse.h"
#include "UCUtil.h"
#include "GeoMagnetic.h"


//my stuff
//#include "windowing.h"
//#include "templates.h"
#include "AnitaTemplates.h"
#include "AnitaNoiseSummary.h"
#include "AnitaNoiseMachine.h"

//lets handle SIGINT correctly so I can close the files
#include "signal.h"


using namespace std;

TFile* outFile = NULL; //the pointer has to be global otherwise the inturrupts can't call it
void emergencyClose(int sig) {
  
  if (sig==15) cout << endl << "emergencyClose(): SIGTERM Signal caught! Quitting!" << endl;
  else if (sig==3) cout << endl << "emergencyClose(): SIGQUIT Signal caught! Qutting!" << endl;
  else cout << endl << "emergencyClose(): A signal was caught!" << endl;
  cout << endl << "Okay I guess you want me to stop :(  Just kidding I'm not sad!  Let me save first...  " << endl;

  if (outFile != NULL) {
    outFile->cd();
    outFile->Write(); //I think this causes a segfault occasionally...
    outFile->Close();
    cout << "Saved successfully!" << endl;
  }
  else {
    cout << "Warning!  I don't have a file open so I didn't actually do anthing right there" << endl;
  }

  cout << "Goodnight! :D" << endl;

  exit(EXIT_SUCCESS);
}

void systemPause() {
  cout << endl << "systemPause(): SIGINT Signal caught! Pausing!" << endl;
  cout << endl << "Writing data and pausing so that you can read what I've done without ROOT errors" << endl;

  if (outFile != NULL) {
    outFile->cd();
    outFile->Write();
    cout << "Written successfully!" << endl;
  }
  else {
    cout << "Warning!  I don't have a file open so I didn't actually do anthing right there" << endl;
  }

}

volatile sig_atomic_t paused = false;

void pauseSwitch(int sig) {

  cout << endl << "^^^ pauseSwitch(): SIGINT Signal caught!" << endl;

  if (!paused) {
    systemPause();
    cout << "Send another SIGINT signal to continue" << endl;
    sig=0;
    paused = true;
  }
  else {
    cout << "Thanks! Resuming!" << endl;
    paused = false;
  }

}


using namespace std;

/*

  Ben Rotter - March 2017 - University of Hawaii at Manoa

  FINAL STRETCH

  I have a solid impulse response for the system now, and a good template for CR events.  So lets search based on that

  Also switch over to the DataSet type format, which makes everything go real quickly.

 */



void entryToRun(int entry, int &runOut, int &startEntry) {

  ifstream entriesPerRun("entriesPerRun.txt");
  int run,entryLow,entryHigh,numEntries;
  while (entriesPerRun >> run >> entryLow >> entryHigh >> numEntries) {
    if (entry >= entryLow && entry<= entryHigh) {
      runOut = run;
      startEntry = entryLow;
      break;
    }
  }


  entriesPerRun.close();


}


void usage() {
  cout << "Usage: templateSearch [option]" << endl;
  cout << "Available options:" << endl;
  cout << " --entry [output base filename] [start entry] [end entry]   : default, by entry" << endl;
  cout << " --min   [output base filename] [start entry] [end entry]   : same as --entry, but with minbias cut" << endl;
  cout << " --wais  [output base filename] [numSplits]   [splitNumber] : for wais, divides for you, does a wais cut" << endl;

}


int main(int argc, char* argv[]) {
  cout << "♪┏(°.°)┛┗(°.°)┓┗(°.°)┛┏(°.°)┓ ♪" << endl;
  //handle inturrupt signals for real
  signal(SIGTERM, emergencyClose); 
  signal(SIGQUIT, emergencyClose); 
  signal(SIGINT, pauseSwitch);

  /*
    figure out what you're doing
  */ 
  string outFileName;
  int startEntry,endEntry,totalEntriesToDo;
  bool waisFlag = false;
  bool minFlag = false;
  bool entryFlag = false;

  // If you're lost
  if (argc <= 1 ){
    usage();
    return -1;
  }
 else if (strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-h") == 0 ) { 
   cout << argv[1] << " help selected (argc" << argc << ")" << endl;
   usage();
   return -1;
 }
  // Original way, by "entry"
 else if (strcmp("--entry",argv[1]) == 0 || strcmp("--min",argv[1]) == 0) {
   if (argc == 5) {
     if (strcmp(argv[1],"--min") == 0) minFlag = true;
      else                         entryFlag = true;
      outFileName = argv[2];
      startEntry = atoi(argv[3]);
      endEntry = atoi(argv[4]);
      cout << "Hello!  Lets do some physics using entries!" << endl;
      if (minFlag) cout << "---------->Only doing minbias events though!<------------" << endl;
      cout << "Doing entry " << startEntry << " to " << endEntry << endl;
      totalEntriesToDo = endEntry - startEntry;
      cout << "This will be a total of " << totalEntriesToDo << " events to process" << endl;
    }
    else {
      cout << "Usage: templateSearch (--entry || --min) [output base filename] [start entry] [end entry]" << endl;    
      return -1;
    }
  }

  else if (strcmp("--wais",argv[1]) == 0) {
    if (argc == 5) {
      outFileName = argv[2];
      waisFlag = true;
      //a bunch of constants I figured out by looking
      const int waisFirstEventNum = 55576917;
      const int waisLastEventNum = 61702834;
      const int waisFirstEntry = 49759012;
      const int waisLastEntry  = 55860967;
      const int waisEntryRange = waisLastEntry - waisFirstEntry;

      int numSplits = atoi(argv[3]);
      int splitNum  = atoi(argv[4]);
      if (numSplits == 1) splitNum = 0;
      int numEntriesPerSplit = waisEntryRange/numSplits;
      startEntry = waisFirstEntry + (numEntriesPerSplit*(splitNum));
      endEntry   = waisFirstEntry + (numEntriesPerSplit*(splitNum+1));
      totalEntriesToDo = endEntry - startEntry;

      cout << "Hello!  Lets do some physics on wais pulses!" << endl;
      cout << "Splitting the relevent runs into " << numSplits << " splits, and doing split " << splitNum << endl;
      cout << "That means I'm doing entry " << startEntry << " through " << endEntry << endl;
      cout << "This will be a total of " << totalEntriesToDo << " events to process" << endl;
    }
    else {
      cout << "Usage: templateSearch --wais  [output base filename] [numSplits] [splitNumber]" << endl;
      return -1;
    }
  }
  else {
    usage();
    return -1;
  }


  //print out some info
  cout << "Default AnitaVersion is " << AnitaVersion::get() << ", using 3 though" << endl;
  AnitaVersion::set(3);

  char serverName[64];
  gethostname(serverName,64);
  cout << "For reference, this process is running on " << serverName << endl;


  //load wisdom crap that maybe makes the first part go faster (but it doesn't really)
  char* homeDir = getenv("HOME");
  stringstream wisdomDir;
  wisdomDir.str("");
  wisdomDir << homeDir << "/macros/fftWisdom.dat";
  FFTtools::loadWisdom(wisdomDir.str().c_str());


  //all rejoyce: the glorious stringstream declaration
  stringstream name;




  //open up the output file
  name.str("");
  name << outFileName << ".root";
  cout << "Using " << name.str() << " as output file" << endl;
  outFile = TFile::Open(name.str().c_str(),"recreate");
  //LZMA is the "fancy" compression, and 9 is the strongest
  outFile->SetCompressionLevel(9); 
  outFile->SetCompressionAlgorithm(ROOT::kLZMA);


  outFile->cd();
  TTree *outTree = new TTree("summaryTree","summaryTree");

  //Lets make the summary object that I can shove into the output tree
  AnitaEventSummary *eventSummary = new AnitaEventSummary; 
  outTree->Branch("eventSummary",&eventSummary);

  //and add gps data too, don't forget to fill this later!
  Adu5Pat *gpsEvent = NULL;
  outTree->Branch("gpsEvent",&gpsEvent);



  /*  Making the filters */
  cout << "Making the filters" << endl;
  //the spectrum average is used for a couple of filters to make them sort of "adaptive"

  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);

  //Make a filter strategy
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

  //Make (or use) the sine subtract cache
  //  UCorrelator::SineSubtractCache


  /** doing it by construction instead of by key which doesn't work for whatever reason */
  //  with a debug file
  //  name.str("");
  //  name << outFileName << "_filtOutFile.root";
  //  TFile *filterOutFile = TFile::Open(name.str().c_str(),"recreate"); 
  //  FilterStrategy strategy(filterOutFile);
  //  without a debug file
  //  FilterStrategy *fStrat = new FilterStrategy();
  //
  //Sine subtract alghorithm (this is the complicated way to do it)
  //  UCorrelator::SineSubtractFilter *sineSub = new UCorrelator::SineSubtractFilter(10,3);
  //  sineSub->makeAdaptive(specAvgLoader,2);
  //  fStrat->addOperation(sineSub);
  // This seems like it should work and is easier
  //add "adsinsub_2_5_13" (default in MagicDisplay)
  //  UCorrelator::fillStrategyWithKey(fStrat,"sinsub_05_1_ad_1");
  //
  //Brick wall filter, should be way faster
  // UCorrelator::AdaptiveBrickWallFilter(const UCorrelator::SpectrumAverageLoader * spec, double thresh=2, bool fillNotch = true);  
  // Don't fill in the noise because whats the point of that really
  //  UCorrelator::AdaptiveBrickWallFilter *brickWall = new UCorrelator::AdaptiveBrickWallFilter(specAvgLoader,2,false);
  //  fStrat->addOperation(brickWall);
  //
  //  with abby's list of filtering
  //  UCorrelator::applyAbbysFilterStrategy(&fStrat);



  /* CONFIGURATION FOR THE ANALYSIS!!!!!! */
  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 

  //set the response to my "single" response
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;

  //set the stokes windowing to use the entire waveform
  config->windowStokes = false;
  
  //lets try to do only the first peak (I never use the others)
  config->nmaxima = 2;



  //and create an analyzer object
  cout << "Making the Analyzer" << endl;
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); //interactive needs to be true




  cout << "Adding more stuff to output tree:" << endl;
  //and add a thing to store template results
  cout << "+template results" << endl;
  AnitaTemplateSummary *templateSummary = new AnitaTemplateSummary();
  outTree->Branch("template",&templateSummary);
  //and one for the filtered I guess (if I ever fix it)
  AnitaTemplateSummary *templateSummary_filtered = new AnitaTemplateSummary();
  outTree->Branch("template_filtered",&templateSummary_filtered);
  //and the noise summary too
  cout << "+noise summary" << endl;
  AnitaNoiseSummary *noiseSummary = new AnitaNoiseSummary();
  outTree->Branch("noiseSummary",&noiseSummary);

  //the thing that calculates the summary is persistant between events
  cout << "creating noise machine" << endl;
  AnitaNoiseMachine *noiseMachine = new AnitaNoiseMachine();
  noiseMachine->fillMap = false;


  //also lets add the geomagnetic estimates, one for upgoing, one for "reflected"
  cout << "+geomagnetic" << endl;
  double geomagExpect,geomagExpectUp;
  outTree->Branch("geomagExpect",&geomagExpect);
  outTree->Branch("geomagExpectUp",&geomagExpectUp);


  /*=================
  //get the templates
  */
  //The length every waveform I look at from now on will be length, until I window it I guess.
  const int length = 2048;

  AnitaTemplateMachine *templateMachine = new AnitaTemplateMachine(length);
  templateMachine->loadTemplates();
  templateMachine->deconvolveTemplates(apd);

  //should be able to write it to file too
  outFile->cd();
  cout << "Writing templates to file" << endl;
  templateMachine->writeTemplatesToFile(outFile);
  
  /*
    --------*/


  
  //find which run startEntry refers to
  int startRun,startEntryInRun;
  entryToRun(startEntry,startRun,startEntryInRun);
  int entryToStartAt = startEntry - startEntryInRun;
  cout << "startEntry " << startEntry << " starts " << entryToStartAt << " entries into run " << startRun << endl;
  int runToGet = startRun;

  int entriesInCurrRun; //number of entries in current open run file
  int completedRunEvs; //number of entries from previous runs

  //make data storage object pointer for later
  AnitaDataset *data = NULL;



  cout << "templateSearch(): starting event loop" << endl;

  int skippedEvs=0;
  int skippedEvsInst=0;
  int processedEvs=0;
  bool printRate = false;

  Acclaim::ProgressBar p(totalEntriesToDo);
  TStopwatch watch; //!< ROOT's stopwatch class, used to time the progress since object construction
  watch.Start(kTRUE);
  int totalTimeSec;
  

  /*==============================
    Event Loop Begins Here      */
  for (Long64_t entry=0; entry<totalEntriesToDo; entry++) {
    while(paused);

    //entry - tracks how far you are in the requested range
    //entryToGet - which entry you're calling from the run
    int entryToGet = entry + entryToStartAt - completedRunEvs;

    //little bit longer progress bar that normal (Acclaim's is good too but I want my OWN)
    const int refreshRate = 100;
 
    if (processedEvs%refreshRate==0 && processedEvs > 0) {
      if (printRate) {
	int timeElapsed = watch.RealTime(); //+1 to prevent divide by zero error
	totalTimeSec += timeElapsed;
	double totalTimeMin = float(totalTimeSec)/60.;

	double instRateSec = float(refreshRate)/timeElapsed;
	double totalRateSec = float(entry)/totalTimeSec;	
	double totalRateMin = float(entry)/totalTimeMin;

	double procRateSec = float(processedEvs)/totalTimeSec;
	double procRateMin = float(processedEvs)/totalTimeMin;

	double secondsLeft =float(totalEntriesToDo-entry)/totalRateSec;
	double minutesLeft =float(totalEntriesToDo-entry)/totalRateMin;
	
	cout << std::setprecision(3);
	cout << entry << "/" << totalEntriesToDo << " evs in " << totalTimeMin << " mins,";
	cout << " (" << processedEvs << " processed) ";
	if (timeElapsed != 0) {
	  cout << procRateSec << "ev/sec (" << instRateSec << " ev/sec instant) | ";
	}
	cout << secondsLeft << " seconds (" << minutesLeft << " mins) remain";
	cout << " {" << entryToGet << "/" << entriesInCurrRun << " in run}";
	
	cout << endl;
	fflush(stdout);
	watch.Start();
	skippedEvsInst=0;
	printRate = false;
      }
    }
    else printRate = true;
    //end of progress bar


    /* Check to see if you need to load a new AnitaDataset since this can span between runs */
    if (entry == 0 || entryToGet > entriesInCurrRun-1) {
      if (data != NULL) delete data;

      if ( (runToGet >= 257) && (runToGet <= 263) ) runToGet = 264; //dead runs that kill processes
      if (  runToGet == 441) break; //440 doesn't exist, if you get there then break the event loop

      if (entry != 0) entryToStartAt = 0; //startEntryInRun becomes zero after the first runswitch

      data = new AnitaDataset(runToGet,false);
      //set the polarity blinding to be on, unless you're looking at wais pulses specifically maybe
      if (!waisFlag) {data->setStrategy(AnitaDataset::BlindingStrategy::kRandomizePolarity);}
      entriesInCurrRun = data->N();

      cout << "AnitaDataset switched to run " << runToGet << endl;
      cout << "run " << runToGet << " has " << entriesInCurrRun << " entries, and I am starting at " << entryToStartAt;
      cout << ", so there are " << entriesInCurrRun - entryToStartAt << " left in the run" << endl;

      completedRunEvs = entry;
      runToGet++;

    }

    //get all the pointers set right
    data->getEntry(entryToGet);
    
    // the trig type needs bit masking for annoying reasons
    int trigType = data->header()->trigType&0x0F;

    //also I need the usefulAdu5Pat to determine if this is wais or not
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(data->gps());

    //Pick out things that you don't want
    //    So no Vpol triggered events, do want to keep things that aren't rf though
    //    also pay attention to the wais and minbias flags
    if ( 
	(entryFlag && (!data->header()->l3TrigPatternH && trigType == 1)) ||   // entry: if it is not Hpol and triggerred
	(minFlag && trigType==1) ||                                            // minbias: not RF triggered
	(waisFlag && !UCorrelator::isWAISHPol(usefulGPS,data->header(),config)) // wais: if it is wais
	 ) {
      skippedEvs++;
      skippedEvsInst++;
      continue;
    }

    processedEvs++;

    //update the gps data!
    gpsEvent = data->gps();


    //1) calibrate and then filter the event and get a FilteredAnitaEvent back
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    analyzer->analyze(filteredEvent, eventSummary); 

    //do the noise analysis (only update it if it is a min bias though
    if (trigType != 1) {
      noiseSummary->isMinBias = true;
      noiseMachine->updateMachine(analyzer,filteredEvent);
    }
    else {
      noiseSummary->isMinBias = false;
      noiseSummary->deleteHists(); //this might make it compress better if lots of things are zero;
    }
    delete filteredEvent;

    //fill stuff, but only if the noise machine has at least a single entry
    if (!noiseMachine->isJustInitialized()){
      noiseMachine->fillNoiseSummary(noiseSummary);
      noiseMachine->fillEventSummary(eventSummary);
    }

    /*==========
      Remember:
      (H=0, V=1)
    ============*/

    /*=============
    Coherent Waveform
    */

    templateSummary->zeroInternals();
    templateSummary_filtered->zeroInternals();

    int dir=0;//only going to do the first peak
    for (int poli=0; poli<2; poli++) {
      //not filtered
      const AnalysisWaveform *waveform = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,dir,false);
      templateMachine->doTemplateAnalysis(waveform,poli,dir,templateSummary);
      //filtered
      const AnalysisWaveform *waveform_filtered = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,dir,true);
      templateMachine->doTemplateAnalysis(waveform_filtered,poli,dir,templateSummary_filtered);
      
      const AnitaResponse::DeconvolutionMethod *deconv = analyzer->getResponseManager()->getDeconvolutionMethod();
      const AnalysisWaveform *waveform_deconv = analyzer->getDeconvolved((AnitaPol::AnitaPol_t)poli,dir,false);
      templateMachine->doDeconvolvedTemplateAnalysis(waveform_deconv,deconv,poli,dir,templateSummary);
      
      const AnalysisWaveform *waveform_deconv_filtered = analyzer->getDeconvolved((AnitaPol::AnitaPol_t)poli,dir,true);
      templateMachine->doDeconvolvedTemplateAnalysis(waveform_deconv_filtered,deconv,poli,dir,templateSummary_filtered);
      
      
    }
      

    //do the geomagnetic stuff
    double phiRad = eventSummary->peak[0][0].phi * TMath::DegToRad();
    double thetaRad = eventSummary->peak[0][0].theta * TMath::DegToRad();
    geomagExpect = GeoMagnetic::getExpectedPolarisation(*usefulGPS,phiRad,thetaRad);
    geomagExpectUp = GeoMagnetic::getExpectedPolarisationUpgoing(*usefulGPS,phiRad,thetaRad,2000);
	


    //okay all done, now I can write out and move to the next event

    outFile->cd();
    outTree->Fill();
    //  outTree->FlushBaskets(); //maybe will make writing at the end faster?


    delete usefulGPS;
    noiseSummary->deleteHists(); //I don't want to save histograms for every event I guess

    analyzer->clearInteractiveMemory();
  }

  cout << "Processed " << processedEvs << " events, skipped " << skippedEvs << endl;
  cout << endl << "Final Processing Rate: " << float(totalEntriesToDo)/totalTimeSec << "ev/sec" << endl;

  outFile->cd();
  cout << "Writing out to file..." << endl;
  outTree->Write();
  outFile->Close();


  delete eventSummary;

  cout << "Physics complete!  See ya later buddy :)" << endl;

  FFTtools::saveWisdom(wisdomDir.str().c_str());

  return 1;
  
}


  
