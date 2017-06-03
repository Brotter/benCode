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

#include "ProgressBar.h"

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


  


int main(int argc, char** argv) {


  //  FFTtools::loadWisdom("/home/brotter/macros/fftWisdom.dat");
  FFTtools::loadWisdom("/Users/brotter/macros/fftWisdom.dat");

                                                                                                 
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
  AnitaDataset *data = new AnitaDataset(runNum,true);

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
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;

  //and create an analyzer object
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); ;
  
  
  //and get the "averaged" impulse response as the template"
  int length = 2048;
  char* templateDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << templateDir << "/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp";
  TGraph *grTemplateRaw = new TGraph(name.str().c_str());
  //waveforms are normally a little over 1024 so lets pad to 2048 (defined in length above)
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateRaw,length);
  delete grTemplateRaw;
  //and normalize it
  TGraph *grTemplate = normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;

  //and get the FFT of it as well, since we don't want to do this every single event
  FFTWComplex *theTemplateFFT=FFTtools::doFFT(length,grTemplate->GetY());
  delete grTemplate;


  //and add a thing to store whatever it finds
  Double_t templateValueV = 0;
  Double_t templateValueH = 0;
  outTree->Branch("templateValueV",&templateValueV);
  outTree->Branch("templateValueH",&templateValueH);


  //**loop through entries
  //option to have less entries! (-1 is the default, so in case you don't specify)
  if (lenEntries == -1) lenEntries = numEntries;

  int entryIndex=0;

  cout << "templateSearch(): starting event loop" << endl;

  Acclaim::ProgressBar p(lenEntries);
  for (Long64_t entry=0; entry<lenEntries; entry++) {
    //    if (entry%1==0) {
    //    cout << entry << "/" << lenEntries << "(" << entryIndex << ")" << endl;
    fflush(stdout);
    //    }
    //get all the pointers set right
    data->getEntry(entry);
    
    //1) calibrate and then filter the event and get a FilteredAnitaEvent back
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();
    //2) then analyze the filtered event!
    analyzer->analyze(filteredEvent, eventSummary); 
    delete filteredEvent;

    //(H=0, V=1)
    //pull the peak coherence graph out of the analyzer so we can compare it vs the template
    for (int poli=0; poli<2; poli++) {
      const TGraphAligned *coherentAligned = analyzer->getCoherent((AnitaPol::AnitaPol_t)poli,0)->even();
      TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
      //make sure it is the same length as the template
      TGraph *coherent2 = FFTtools::padWaveToLength(coherent,length);
      delete coherent;
      //normalize it
      TGraph *normCoherent = normalizeWaveform(coherent2);
      delete coherent2;
      FFTWComplex *coherentFFT=FFTtools::doFFT(length,normCoherent->GetY());
      delete normCoherent;

      double *dCorr = getCorrelationFromFFT(length,theTemplateFFT,coherentFFT);
      double max = TMath::MaxElement(length,dCorr);
      double min = TMath::Abs(TMath::MinElement(length,dCorr));

      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kHorizontal ) templateValueH = TMath::Max(max,min);
      if ( (AnitaPol::AnitaPol_t)poli == AnitaPol::kVertical ) templateValueV = TMath::Max(max,min);
      delete[] coherentFFT;
      delete[] dCorr;

      p.inc(entry, lenEntries);
    }    

    outFile->cd();
    outTree->Fill();
    
    analyzer->clearInteractiveMemory();
  }


  outFile->cd();
  outTree->Write();
  outFile->Close();

  delete eventSummary;

  //  cout << "Physics complete!  See ya later buddy :)" << endl;

  FFTtools::saveWisdom("/Users/brotter/macros/fftWisdom.dat");

  return 1;
  
}


  
