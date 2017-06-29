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

//my stuff
#include "windowing.h"


//lets handle SIGINT correctly so I can close the files
#include "signal.h"

TFile* outFile = NULL;
void emergencyClose(int sig) {
  
  cout << endl << "^^^ emergencyClose(): SIGINT Signal caught!" << endl;
  cout << endl << "Okay sounds good, Let me save things first though :) " << endl;
  if (outFile != NULL) {
    outFile->cd();
    outFile->Write();
    outFile->Close();
  }
  else {
    cout << "Warning!  I don't have a file open so I didn't actually do anthing right there" << endl;
  }

  cout << "Goodnight! :D" << endl;

  exit(EXIT_SUCCESS);
}
    


using namespace std;

/*

  Ben Rotter - March 2017 - University of Hawaii at Manoa

  FINAL STRETCH

  I have a solid impulse response for the system now, and a good template for CR events.  So lets search based on that

  Also switch over to the DataSet type format, which makes everything go real quickly.

 */


TGraph *normalizeWaveform(TGraph *inGraph) {
  
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}




double *getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2) 
{


    int newLength=(length/2)+1;
//     cout << "newLength " << newLength << endl;
    FFTWComplex *tempStep = new FFTWComplex [newLength];
    int no2=length>>1;
    for(int i=0;i<newLength;i++) {
	double reFFT1=theFFT1[i].re;
	double imFFT1=theFFT1[i].im;
	double reFFT2=theFFT2[i].re;
	double imFFT2=theFFT2[i].im;

	//Real part of output 
	tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2)/double(no2/2);
	//Imaginary part of output 
	tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2)/double(no2/2);
    }
//    cout << "finished messing around" << endl;
    double *theOutput=FFTtools::doInvFFT(length,tempStep);
//    cout << "got inverse" << endl;
    delete [] tempStep;
    return theOutput;

}



FFTWComplex* getImpulseResponseTemplate(int length) {
  
  //and get the "averaged" impulse response as the template"
  char* templateDir = getenv("ANITA_UTIL_INSTALL_DIR");
  stringstream name;
  name.str("");
  name << templateDir << "/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp";
  TGraph *grTemplateRaw = new TGraph(name.str().c_str());
  //waveforms are normally a little over 1024 so lets pad to 2048 (defined in length above)
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateRaw,length);
  delete grTemplateRaw;
  //then cut it back down with a window function
  int peakHilb = -1;
  TGraph *grTemplateCut = windowDispersed(grTemplatePadded,peakHilb);
  delete grTemplatePadded;
  //and finally normalize it (last step!)
  TGraph *grTemplate = normalizeWaveform(grTemplateCut);
  delete grTemplateCut;

  //save it quickly
  outFile->cd();
  grTemplate->SetName("templateImp");
  grTemplate->Write();

  //and get the FFT of it as well, since we don't want to do this every single event
  FFTWComplex *theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  delete grTemplate;
  
  return theTemplateFFT;
}

void getCRTemplates(int length, const int numTemplates,FFTWComplex** theTemplates) {

  stringstream name;


  //  TFile *inFile = TFile::Open("/Users/brotter/benCode/ZHAiresReader/convolveCRWithSigChain.root");
  TFile *inFile = TFile::Open("convolveCRWithSigChain.root");
  
  for (int disp=0; disp<2; disp++) {
    for (int i=0; i<numTemplates; i++) {
      //want to get graphs 13 through 24 (like in makeTemplate.C)
      int wave = i+13; //peak seems to be at around the 13th one, then by 23 it is basically zero
      name.str("");
      if (disp==0) name << "disp";
      if (disp==1) name << "efield";
      name << wave;
      TGraph *grTemplateRaw = (TGraph*)inFile->Get(name.str().c_str());

      //waveforms are super long so we can just cut it to the window dimentions
      TGraph *grTemplateCut;
      int peakHilb = -1;
      if (disp==0) grTemplateCut = windowDispersed(grTemplateRaw,peakHilb);
      if (disp==1) grTemplateCut = windowEField(grTemplateRaw,peakHilb);
      delete grTemplateRaw;
       

      //and finally normalize it (last step!)
      TGraph *grTemplate = normalizeWaveform(grTemplateCut);
      delete grTemplateCut;

      //save it quickly
      outFile->cd();
      grTemplate->SetName(name.str().c_str());
      grTemplate->Write();

      //and get the FFT of it as well, since we don't want to do this every single event
      FFTWComplex *theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
      delete grTemplate;
      
      theTemplates[i+numTemplates*disp] = theTemplateFFT;
    }
  }
  inFile->Close();

  return;
}

    
FFTWComplex* getWaisTemplate(int length) {
  
  //and get the "averaged" impulse response as the template"
  TFile *inFile = TFile::Open("waisTemplate.root");
  TGraph *grTemplateRaw = (TGraph*)inFile->Get("wais01TH");
  //the wais waveform is like N=2832, but most of it is dumb, so cut off the beginning
  //actually just window it!
  int peakHilb = -1;
  TGraph *grTemplateCut = windowDispersed(grTemplateRaw,peakHilb);
  delete grTemplateRaw;

  //and then normalize it
  TGraph *grTemplate = normalizeWaveform(grTemplateCut);
  delete grTemplateCut;

  //save it quickly  
  outFile->cd();
  grTemplate->SetName("templateWais");
  grTemplate->Write();

  //and get the FFT of it as well, since we don't want to do this every single event
  FFTWComplex *theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  delete grTemplate;
  
  return theTemplateFFT;
}


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

  char* homeDir = getenv("HOME");
  stringstream wisdomDir;
  wisdomDir.str("");
  wisdomDir << homeDir << "/macros/fftWisdom.dat";
  FFTtools::loadWisdom(wisdomDir.str().c_str());

  //handle inturrupt signals for real
  signal(SIGINT, emergencyClose); 
                                                                                                 
  string outFileName;
  int startEntry,endEntry,totalEntriesToDo;

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


  //the spectrum average is used for a couple of filters to make them sort of "adaptive"
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);


  //Sine subtract alghorithm (this is the complicated way to do it)
  //  UCorrelator::SineSubtractFilter *sineSub = new UCorrelator::SineSubtractFilter(0.05,2);
  //  sineSub->makeAdaptive(specAvgLoader);
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



  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  //set the response to my "single" response
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;

  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); ;


  //The length every waveform I look at from now on will be length, until I window it I guess.
  const int length = 2048;



  /*=================
  //get the templates
  */
  FFTWComplex *theTemplateFFT = getImpulseResponseTemplate(length);
  FFTWComplex *theWaisTemplateFFT = getWaisTemplate(length);
  const int numCRTemplates = 10;
  FFTWComplex *theCRTemplates[numCRTemplates*2]; //LSB=dispersed, MSB=efield
  getCRTemplates(length,numCRTemplates,theCRTemplates);
  /*
    --------*/



  //and add a thing to store whatever it finds
  Double_t templateImpV = 0;
  Double_t templateImpH = 0;
  outTree->Branch("templateImpV",&templateImpV);
  outTree->Branch("templateImpH",&templateImpH);

  //one for the WAIS template too
  Double_t templateWaisV = 0;
  Double_t templateWaisH = 0;
  outTree->Branch("templateWaisV",&templateWaisV);
  outTree->Branch("templateWaisH",&templateWaisH);

  //and for the bigger multi-coherence-angle one
  Double_t templateCRayV[numCRTemplates][2];
  Double_t templateCRayH[numCRTemplates][2];
  name.str("");
  name << "templateCRayV[" << numCRTemplates << "][2]/D";
  outTree->Branch("templateCRayV",&templateCRayV,name.str().c_str());
  name.str("");
  name << "templateCRayH[" << numCRTemplates << "][2]/D";
  outTree->Branch("templateCRayH",&templateCRayH,name.str().c_str());
        
  
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

  Acclaim::ProgressBar p(totalEntriesToDo);
  TStopwatch watch; //!< ROOT's stopwatch class, used to time the progress since object construction
  watch.Start(kTRUE);
  int timeElapsed;
  int totalTime;
  for (Long64_t entry=0; entry<totalEntriesToDo; entry++) {
    //entry - tracks how far you are in the requested range
    //entryToGet - which entry you're calling from the run
    int entryToGet = entry + entryToStartAt - completedRunEvs;

    //little bit longer progress bar that normal (Acclaim's is good too but I want my OWN)
    const int refreshRate = 100;
    if (entry%refreshRate==0 || entry<10) {
      timeElapsed = watch.RealTime(); //+1 to prevent divide by zero error
      totalTime += timeElapsed;
      cout << entry << "/" << totalEntriesToDo << " evs in " << float(totalTime)/60 << " mins ";
      if (timeElapsed != 0) {
	cout << float(entry)/totalTime << "ev/sec (" << float(refreshRate)/timeElapsed << " ev/sec instant)";
      }
      cout << " {" << entryToGet << "/" << entriesInCurrRun << "}";

      cout << endl;
      fflush(stdout);
      watch.Start();
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
    
    //0) If I'm going to re-run this again, I want to at least cut it in half
    //    So no Vpol triggered events
    if (!data->header()->l3TrigPatternH) continue;


    //1) calibrate and then filter the event and get a FilteredAnitaEvent back
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    analyzer->analyze(filteredEvent, eventSummary); 
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
      const TGraphAligned *coherentAligned = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,0,false)->even();
      TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      TGraph *coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      int peakHilbert = -1;
      TGraph *windowed = windowDispersed(coherent2,peakHilbert);
      delete coherent2;
      int newLength = windowed->GetN();

      //and xpol for stokes
      coherentAligned = analyzer->getCoherentXpol((AnitaPol::AnitaPol_t)poli,0,false)->even();
      coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      TGraph *windowedXpol = windowDispersed(coherent2,peakHilbert);
      delete coherent2;
      if (newLength != windowedXpol->GetN()) cout << "xpol is a different length" << endl;

      //and do stokes for that
      fillStokes(windowed->GetN(),windowed,windowedXpol,eventSummary->coherent[poli][1]);
      delete windowedXpol;


      //normalize it
      TGraph *normCoherent = normalizeWaveform(windowed);
      delete windowed;
      FFTWComplex *coherentFFT=FFTtools::doFFT(newLength,normCoherent->GetY());
      delete normCoherent; 

      double *dCorr = getCorrelationFromFFT(newLength,theTemplateFFT,coherentFFT);
      double max = TMath::MaxElement(newLength,dCorr);
      double min = TMath::Abs(TMath::MinElement(newLength,dCorr));

      double *dCorrWais = getCorrelationFromFFT(newLength,theWaisTemplateFFT,coherentFFT);
      double maxWais = TMath::MaxElement(newLength,dCorrWais);
      double minWais = TMath::Abs(TMath::MinElement(newLength,dCorrWais));


      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal ) {
	templateImpH = TMath::Max(max,min);
	templateWaisH = TMath::Max(maxWais,minWais);
      }
      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
	templateImpV = TMath::Max(max,min);
	templateWaisV = TMath::Max(maxWais,minWais);
      }

      double maxCR[numCRTemplates];
      double minCR[numCRTemplates];
      for (int i=0; i<numCRTemplates; i++) {
	double *dCorrCR = getCorrelationFromFFT(newLength,theCRTemplates[i],coherentFFT);
	maxCR[i] = TMath::MaxElement(newLength,dCorrCR);
	minCR[i] = TMath::Abs(TMath::MinElement(newLength,dCorrCR));
	delete[] dCorrCR;
	if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) {
	  templateCRayV[i][0] = TMath::Max(maxCR[i],minCR[i]);
	}
	else if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal )  {
	  templateCRayH[i][0] = TMath::Max(maxCR[i],minCR[i]);
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
    */
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
      int newLength = windowed->GetN();

      //and xpol for stokes
      coherentAligned = analyzer->getDeconvolvedXpol((AnitaPol::AnitaPol_t)poli,0)->even();
      coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      TGraph *windowedXpol = windowEField(coherent2,peakHilbert);
      delete coherent2;
      if (newLength != windowedXpol->GetN()) cout << "xpol has a different length" << endl;

      //and do stokes for that
      fillStokes(windowed->GetN(),windowed,windowedXpol,eventSummary->deconvolved[poli][1]);
      delete windowedXpol;
       
      //normalize it
      TGraph *normCoherent = normalizeWaveform(windowed);
      delete windowed;
      FFTWComplex *coherentFFT=FFTtools::doFFT(newLength,normCoherent->GetY());
      delete normCoherent; 

      double maxCR[numCRTemplates];
      double minCR[numCRTemplates];
      for (int i=0; i<numCRTemplates; i++) {
	double *dCorrCR = getCorrelationFromFFT(newLength,theCRTemplates[i+10],coherentFFT);
	maxCR[i] = TMath::MaxElement(newLength,dCorrCR);
	minCR[i] = TMath::Abs(TMath::MinElement(newLength,dCorrCR));
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


    //okay all done, now I can write out and move to the next event

    outFile->cd();
    outTree->Fill();
    //  outTree->FlushBaskets(); //maybe will make writing at the end faster?

    analyzer->clearInteractiveMemory();
  }


  cout << endl << "Final Processing Rate: " << float(totalEntriesToDo)/totalTime << "ev/sec" << endl;

  outFile->cd();
  cout << "Writing out to file..." << endl;
  outTree->Write();
  outFile->Close();

  delete eventSummary;

  cout << "Physics complete!  See ya later buddy :)" << endl;

  FFTtools::saveWisdom(wisdomDir.str().c_str());

  return 1;
  
}


  
