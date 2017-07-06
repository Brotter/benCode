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

//my stuff
#include "windowing.h"
#include "templates.h"
#include "AnitaTemplateResults.h"
#include "AnitaNoiseSummary.h"

//lets handle SIGINT correctly so I can close the files
#include "signal.h"


TFile* outFile = NULL; //the pointer has to be global otherwise the inturrupts can't call it
void emergencyClose(int sig) {
  
  /*
  if (paused) {
    unpause();
    return;
  }
  */
  if (sig==15) cout << endl << "emergencyClose(): SIGTERM Signal caught!" << endl;
  else if (sig==3) cout << endl << "emergencyClose(): SIGQUIT Signal caught!" << endl;
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


void fillStokes(int length,TGraph *window, TGraph *windowXpol, AnitaEventSummary::WaveformInfo wave) {

    //calculate the windowed stokes parameters
    double *I = new double;
    double *Q = new double;
    double *U = new double;
    double *V = new double;
    double *hilbertForStokes = FFTtools::getHilbertTransform(length,window->GetY());
    double *hilbertForStokesXpol = FFTtools::getHilbertTransform(length,windowXpol->GetY());
    FFTtools::stokesParameters(length,window->GetY(),hilbertForStokes,windowXpol->GetY(),hilbertForStokesXpol,
			       I,Q,U,V);

    delete[] hilbertForStokes;
    delete[] hilbertForStokesXpol;

    wave.I = *I;
    wave.Q = *Q;
    wave.U = *U;
    wave.V = *V;

    delete I;
    delete Q;
    delete U;
    delete V;

  return;
}




int main(int argc, char** argv) {

  //load wisdom crap that doesn't seem to do anything
  char* homeDir = getenv("HOME");
  stringstream wisdomDir;
  wisdomDir.str("");
  wisdomDir << homeDir << "/macros/fftWisdom.dat";
  FFTtools::loadWisdom(wisdomDir.str().c_str());

  //handle inturrupt signals for real
  signal(SIGTERM, emergencyClose); 
  signal(SIGQUIT, emergencyClose); 
  signal(SIGINT, pauseSwitch);
    
  string outFileName;
  int startEntry,endEntry,totalEntriesToDo;

  //figure out what you're doing
  if (argc==4) {
    outFileName = argv[1];
    startEntry = atoi(argv[2]);
    endEntry = atoi(argv[3]);
    cout << "Hello!  Let us do some physics mate" << endl;
    char serverName[64];
    gethostname(serverName,64);
    cout << "For reference, this process is running on " << serverName << endl;
    cout << "Doing entry " << startEntry << " to " << endEntry << endl;
    totalEntriesToDo = endEntry - startEntry;
    cout << "This will be a total of " << totalEntriesToDo << " events to process" << endl;
    }
  else {
    cout << "Usage: " << argv[0] << " [output base filename] [start entry] [end entry]" << endl;
    return -1;
  }


  stringstream name;
  //also make a Tree of event headers that pass cuts
  name.str("");
  name << outFileName << ".root";
  cout << "Using " << name.str() << " as output file" << endl;
  outFile = TFile::Open(name.str().c_str(),"recreate");
  outFile->SetCompressionLevel(9); //LZMA is the "fancy" compression, and 9 is the strongest
  outFile->SetCompressionAlgorithm(ROOT::kLZMA);


  outFile->cd();
  TTree *outTree = new TTree("summaryTree","summaryTree");

  //Lets make the summary object that I can shove into the output tree
  AnitaEventSummary *eventSummary = new AnitaEventSummary; 

  outTree->Branch("eventSummary",&eventSummary);



  cout << "Making the filters" << endl;
  //Make a filter strategy
  //  with a debug file
  //  name.str("");
  //  name << outFileName << "_filtOutFile.root";
  //  TFile *filterOutFile = TFile::Open(name.str().c_str(),"recreate"); 
  //  FilterStrategy strategy(filterOutFile);
  //  without a debug file
  FilterStrategy *strategy = new FilterStrategy();


  //the spectrum average is used for a couple of filters to make them sort of "adaptive"
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);


  //Sine subtract alghorithm (this is the complicated way to do it)
  UCorrelator::SineSubtractFilter *sineSub = new UCorrelator::SineSubtractFilter(0.05,2);
  sineSub->makeAdaptive(specAvgLoader);
  //  strategy->addOperation(sineSub);
  // This seems like it should work and is easier
  //add "adsinsub_2_5_13" (default in MagicDisplay)
  //  UCorrelator::fillStrategyWithKey(strategy,"sinsub_05_1_ad_1");
  

  //Brick wall filter, should be way faster
  // UCorrelator::AdaptiveBrickWallFilter(const UCorrelator::SpectrumAverageLoader * spec, double thresh=2, bool fillNotch = true);  
  // Don't fill in the noise because whats the point of that really
  UCorrelator::AdaptiveBrickWallFilter *brickWall = new UCorrelator::AdaptiveBrickWallFilter(specAvgLoader,2,false);
  strategy->addOperation(brickWall);


  //  with abby's list of filtering
  //  UCorrelator::applyAbbysFilterStrategy(&strategy);



  cout << "Making the Analyzer" << endl;
  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  //set the response to my "single" response
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;

  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); ;


  //The length every waveform I look at from now on will be length, until I window it I guess.
  const int length = 2048;



  cout << "Adding more stuff to output tree" << endl;
  //and add a thing to store template results
  AnitaTemplateResults *templateResults = new AnitaTemplateResults();
  outTree->Branch("templateResults",&templateResults);

  //and the noise summary too
  AnitaNoiseSummary *noiseSummary = new AnitaNoiseSummary();
  outTree->Branch("noiseSummary",&noiseSummary);

  //the thing that calculates the summary is persistant between events
  AnitaNoiseMachine *noiseMachine = new AnitaNoiseMachine();
  noiseMachine->fillMap = true;

  /*=================
  //get the templates
  */
  cout << "Getting the templates" << endl;
  FFTWComplex *theImpTemplateFFT;
  TGraph *theImpTemplate;
  getImpulseResponseTemplate(length,theImpTemplateFFT,theImpTemplate);
  cout << "Got Impulse Response Template" << endl;
  cout << theImpTemplateFFT << " " << theImpTemplate << endl;

  FFTWComplex *theWaisTemplateFFT;
  TGraph *theWaisTemplate;
  getWaisTemplate(length,theWaisTemplateFFT,theWaisTemplate);
  outFile->cd();
  cout << theWaisTemplate << " " << theWaisTemplateFFT << endl;
  theWaisTemplate->Write();
  cout << "Got Wais Template" << endl;

  const int numCRTemplates = templateResults->numCRTemplates;
  cout << "numCRTemplates = " << numCRTemplates << endl;
  FFTWComplex *theCRTemplateFFTs[numCRTemplates];
  TGraph *theCRTemplates[numCRTemplates];
  getCRTemplates(length,numCRTemplates,theCRTemplateFFTs,theCRTemplates);

  //also write them
  outFile->cd();
  theImpTemplate->Write();

  for (int i=0; i<numCRTemplates; i++) {
    theCRTemplates[i]->Write();
  }
  cout << "Wrote templates to output file" << endl;

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
 
    if (entry%refreshRate==0 || entry<10) {
      int timeElapsed = watch.RealTime(); //+1 to prevent divide by zero error
      totalTimeSec += timeElapsed;

      double totalTimeMin = float(totalTimeSec)/60.;
      double totalRateSec = float(entry-skippedEvs)/totalTimeSec;
      double totalRateMin = float(entry-skippedEvs)/totalTimeMin;
      double instRateSec = float(refreshRate-skippedEvsInst)/timeElapsed;
      double minutesLeft =float(totalEntriesToDo-entry)/totalRateMin;
      double secondsLeft =float(totalEntriesToDo-entry)/totalRateSec;

      cout << std::setprecision(3);
      cout << entry << "/" << totalEntriesToDo << " evs in " << totalTimeMin << "mins,";
      cout << " skipped " << skippedEvs << "(+" << skippedEvsInst << ") ";
      if (timeElapsed != 0) {
	cout << totalRateSec << "ev/sec (" << instRateSec << " ev/sec instant)";
      }
      cout << secondsLeft << " seconds (" << minutesLeft << " mins) remain";
      cout << " {" << entryToGet << "/" << entriesInCurrRun << "}";

      cout << endl;
      fflush(stdout);
      watch.Start();
      skippedEvsInst=0;
    }


    if (entry == 0 || entryToGet > entriesInCurrRun-1) {
      if (data != NULL) delete data;
      data = new AnitaDataset(runToGet,false);
      data->setStrategy(AnitaDataset::BlindingStrategy::kRandomizePolarity);

      cout << "AnitaDataset switched to run " << runToGet << endl;
      
      completedRunEvs = entry;
      runToGet++;
      entriesInCurrRun = data->N();
      if (entry != 0) entryToStartAt = 0; //startEntryInRun becomes zero after the first runswitch

      cout << "New run has " << entriesInCurrRun << " entries, and I am starting at " << entryToStartAt;
      cout << ", so there are " << entriesInCurrRun - entryToStartAt << " left in the run" << endl;


    }

    //get all the pointers set right
    data->getEntry(entryToGet);
    
    // the trig type needs bit masking for annoying reasons
    int trigType = data->header()->trigType&0x0F;

    //0) If I'm going to re-run this again, I want to at least cut it in half
    //    So no Vpol triggered events
    //    do want to keep things that aren't rf though
    if (!data->header()->l3TrigPatternH && trigType == 1) {
      skippedEvs++;
      skippedEvsInst++;
      continue;
    }


    //1) calibrate and then filter the event and get a FilteredAnitaEvent back
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    analyzer->analyze(filteredEvent, eventSummary); 

    //do the noise analysis (only update it if it is a min bias though
    if (trigType != 1) {
      noiseSummary->isMinBias = true;
      noiseMachine->fillAvgRMSNoise(filteredEvent);
      noiseMachine->fillAvgMapNoise(analyzer);
      noiseMachine->fillNoiseSummary(noiseSummary);
    }
    else {
      noiseSummary->isMinBias = false;
    }
    delete filteredEvent;

    /*==========
      Remember:
      (H=0, V=1)
    ============*/

    /*=============
    Coherent Waveform
    */
    for (int poli=0; poli<2; poli++) {
      //get coherently aligned waveform
      const AnalysisWaveform *coherentAnalysis = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,0,false);
      //I actually want to do SOME filtering though... so sine subtract a single one?
      //      sineSub->processOne(coherentAnalysis,data->header(),);
      const TGraphAligned *coherentAligned = coherentAnalysis->even();
      TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      TGraph *coherentPad = FFTtools::padWaveToLength(coherent,length);
      delete coherent;


      //-------stokes stuff
      //window it for stokes parameters
      int peakHilbert = -1;
      TGraph *windowed = windowDispersed(coherentPad,peakHilbert);
      //      delete coherentPad; //don't delete the padded one because its better for the template
      int windLength = windowed->GetN();


      //and xpol for stokes
      const AnalysisWaveform *coherentAnalysisXpol = analyzer->getCoherentXpol((AnitaPol::AnitaPol_t)poli,0,false);
      //I actually want to do SOME filtering though... so sine subtract a single one?
      //      sineSub->processOne(coherentAnalysisXpol);
      const TGraphAligned *coherentAlignedXpol = coherentAnalysisXpol->even();
      TGraph *coherentXpol = new TGraph(coherentAlignedXpol->GetN(),coherentAlignedXpol->GetX(),coherentAlignedXpol->GetY());
      //make sure it is the same length as the template
      TGraph *coherentPadXpol = FFTtools::padWaveToLength(coherentXpol,length);
      delete coherentXpol;
      TGraph *windowedXpol = windowDispersed(coherentPadXpol,peakHilbert);
      delete coherentPadXpol;
      if (windLength != windowedXpol->GetN()) cout << "xpol is a different length" << endl;

      //and do stokes for that, then you can throw out the windowed things because thats all I use them for
      fillStokes(windowed->GetN(),windowed,windowedXpol,eventSummary->coherent[poli][1]);
      delete windowedXpol;
      delete windowed;
      //--------done with stokes stuff


      //==== Template stuff
      //normalize coherently summed waveform
      TGraph *normCoherent = normalizeWaveform(coherentPad);
      delete coherentPad;
      FFTWComplex *coherentFFT=FFTtools::doFFT(length,normCoherent->GetY());
      delete normCoherent; 

      double *dCorr = getCorrelationFromFFT(length,theImpTemplateFFT,coherentFFT);
      double max = TMath::MaxElement(length,dCorr);
      double min = TMath::Abs(TMath::MinElement(length,dCorr));

      double *dCorrWais = getCorrelationFromFFT(length,theWaisTemplateFFT,coherentFFT);
      double maxWais = TMath::MaxElement(length,dCorrWais);
      double minWais = TMath::Abs(TMath::MinElement(length,dCorrWais));


      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal ) {
	templateResults->coherentH->templateImp  = TMath::Max(max,min);
	templateResults->coherentH->templateWais = TMath::Max(maxWais,minWais);
      }
      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
	templateResults->coherentV->templateImp  = TMath::Max(max,min);
	templateResults->coherentV->templateWais = TMath::Max(maxWais,minWais);
      }

      double maxCR[numCRTemplates];
      double minCR[numCRTemplates];
      for (int i=0; i<numCRTemplates; i++) {
	double *dCorrCR = getCorrelationFromFFT(length,theCRTemplateFFTs[i],coherentFFT);
	maxCR[i] = TMath::MaxElement(length,dCorrCR);
	minCR[i] = TMath::Abs(TMath::MinElement(length,dCorrCR));
	delete[] dCorrCR;
	if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
	  templateResults->coherentV->templateCRay[i] = TMath::Max(maxCR[i],minCR[i]);
	}
	else if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal )  {
	  templateResults->coherentH->templateCRay[i] = TMath::Max(maxCR[i],minCR[i]);
	}
	else {
	  cout << "Something is horrible wrong!  I don't know where to save the templates results! ";
	  cout << "poli=" << poli << endl;
	}

      }
      
      delete[] coherentFFT;
      delete[] dCorr;
      delete[] dCorrWais;


      //      p.inc(entry, totalEntriesToDo); //BenS's status bar (I like mine better)
    }
    /*
      ----------*/

    /*=============
      Deconvolved
      fucked up right now
    //do everything once more for the deconvolved coherent waveform
    for (int poli=0; poli<2; poli++) {

      //get coherently aligned waveform
      const TGraphAligned *coherentAligned = analyzer->getDeconvolved((AnitaPol::AnitaPol_t)poli,0)->even();
      TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      TGraph *coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      int peakHilbert = -1;
      TGraph *windowed = windowEField(coherent2,peakHilbert);
      delete coherent2;
      int windLength = windowed->GetN();

      //and xpol for stokes
      coherentAligned = analyzer->getDeconvolvedXpol((AnitaPol::AnitaPol_t)poli,0)->even();
      coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      TGraph *windowedXpol = windowEField(coherent2,peakHilbert);
      delete coherent2;
      if (windLength != windowedXpol->GetN()) cout << "xpol has a different length" << endl;

      //and do stokes for that
      fillStokes(windowed->GetN(),windowed,windowedXpol,eventSummary->deconvolved[poli][1]);
      delete windowedXpol;
       
      //normalize it
      TGraph *normCoherent = normalizeWaveform(windowed);
      delete windowed;
      FFTWComplex *coherentFFT=FFTtools::doFFT(windLength,normCoherent->GetY());
      delete normCoherent; 


      double maxCR[numCRTemplates];
      double minCR[numCRTemplates];
      for (int i=0; i<numCRTemplates; i++) {
	double *dCorrCR = getCorrelationFromFFT(windLength,theCRTemplates[i+10],coherentFFT);
	maxCR[i] = TMath::MaxElement(windLength,dCorrCR);
	minCR[i] = TMath::Abs(TMath::MinElement(windLength,dCorrCR));
	delete[] dCorrCR;
	if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
	  templateCRayV[i][1] = TMath::Max(maxCR[i],minCR[i]);
	}
	if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal )  {
	  templateCRayH[i][1] = TMath::Max(maxCR[i],minCR[i]);
	}
      }
      
      delete[] coherentFFT;


      //      p.inc(entry, totalEntriesToDo); //BenS's status bar (I like mine better)
    }
*/

    //okay all done, now I can write out and move to the next event

    outFile->cd();
    outTree->Fill();
    //  outTree->FlushBaskets(); //maybe will make writing at the end faster?

    noiseSummary->deleteHists(); //I don't want to save histograms for every event I guess

    analyzer->clearInteractiveMemory();
  }


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


  
