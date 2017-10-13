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
#include "AnitaConventions.h"
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
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"

#include "AnitaTemplates.h"
#include "AnitaNoiseSummary.h"
#include "AnitaNoiseMachine.h"


using namespace std;

double calcRMS(AnitaGeomTool *geom, const UCorrelator::AntennaPositions *ap, double peakPhi, AnitaNoiseSummary *noiseSum) {
/*
  from makeSNRFile.C

  Should be re-calculating this because it is fast instead of loading a huge tree index into memory
 */
  const int antennasInSum = 15; //I'm using 15 antennas
  const int phiSectorsInSum = antennasInSum/3;

  int peakPhiSector;
  int closest[antennasInSum];
      

  int nant = ap->getClosestAntennas(peakPhi, antennasInSum, closest, 0, AnitaPol::kHorizontal);
  if (nant != antennasInSum) {
    cout << "Only using " << nant << " antennas to compute noise" << endl;
  }

  int phiSectors[phiSectorsInSum];
  int phiSector = 0;
  peakPhiSector = -1;
  for (int ant=0; ant<antennasInSum; ant++) {
    if (closest[ant] < 16) { //0->15 is top ring
      phiSectors[phiSector] = geom->getPhiFromAnt(closest[ant]);
      if (peakPhiSector == -1) peakPhiSector = phiSectors[phiSector];
      phiSector++;
    }
  }
  
  double noiseAvg = 0;
  for (phiSector=0; phiSector<phiSectorsInSum; phiSector++) {
    for (int ring=0; ring<3; ring++) {
      double currNoise = noiseSum->avgRMSNoise[phiSectors[phiSector]][ring][0];
      noiseAvg += currNoise;
    }
  }
  noiseAvg /= noiseSum->fifoLength;
  noiseAvg /= antennasInSum;
  
  return noiseAvg;
}







int main(int argc,char **argv) {

  /*

    The dedispersion has a pretty significant problem:
       the responses are tuned so that they look nice, but this puts their zero-delay at like 20ns and pushes the dedispersed waveforms off
       the front of the window, truncating them and making all their results useless.  Lets fix that

    Also it looks like the 50ns Stokes window is too short?  Or the ifirst parameter in Analyzer is wrong?  Either way I have some 
    garbage values for that.  I can re-do that here, no windowing, lets try that.

    Also the SNR value is dumb, replace that with the one I like from the pre-calculated RMS tables in the snrFile.root

    Also lets just unblind right here.  BOOM LETS DO IT NO BLINDING.

    So I'm overwriting:
    - stokes parameters
    - dedispersed waveforms info
    - snr stuff
    - template correlation stuff

   */
  stringstream name;

  int split,numSplits;
  if (argc==1) {
    split = 0;
    numSplits = 1;
    cout << "Defaulting to doing the entire dataset" << endl;
  }
  else if (argc==3) {
    numSplits = atoi(argv[1]);
    split = atoi(argv[2]);
    cout << "Doing " << numSplits << " splits, of which this is number " << split << endl;
  }
  else {
    cout << "Usage: ./reCalcKeyValues [numSplits] [split]" << endl;
}


  //grab only the events that I care about
  char* dataDir = getenv("ANITA3_RESULTSDIR");
  name << dataDir << "/templateSearch/09.27.17_19h/goodEvents.root";
  TFile *inFile = TFile::Open(name.str().c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  int totEntries = summaryTree->GetEntries();
  cout << "goodEvents.root has " << totEntries << " entries" << endl;


  //I'll need to split this up onto the servers, because it takes forever
  int startEntry,stopEntry,lenEntries;
  if (numSplits == 1) {
    startEntry=0;
    stopEntry=totEntries;
    lenEntries=totEntries;
  }
  else {
    lenEntries = totEntries/numSplits;
    startEntry = split*lenEntries;
    stopEntry = (split+1)*lenEntries;
    if (split==numSplits-1) {
      stopEntry = totEntries;
      lenEntries = totEntries-startEntry;
    }
    cout << "Splitting into " << numSplits << " sections, which means " << lenEntries << " events for this section" << endl;
    cout << "Doing section: " << split << ", starting at entry " << startEntry << " and stopping at " << stopEntry << endl;
  }

  //open an output file
  name.str("");
  name << dataDir << "/templateSearch/09.27.17_19h/reCalcKeyValues/goodEvents_corr";
  if (numSplits > 1) name << "_" << split;
  name << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  if (!outFile->IsOpen()) {
    cout << "Didn't open output file!  Quitting." << endl;
    return -1;
  }
  TTree *outTree = new TTree("summaryTree","summaryTree");  

  //link it all together
  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *templateSummary = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&templateSummary);
  outTree->Branch("template",&templateSummary);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);

  //need all the full root data too.  Start at 130, no decimation, and no blinding strat
  AnitaDataset *data = new AnitaDataset(130,false);

  //some config stuff
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig();
  //set the response to my "single" response, which I can shift really easily
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;
  //and to make the previous one obsolete, make it calculate the offsets compared to some central point
  config->delay_to_center = true;


  //filtering strat
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");


  //get the responses
  AnitaResponse::ResponseManager *responses = new AnitaResponse::ResponseManager(UCorrelator::AnalysisConfig::getResponseString(config->response_option),config->response_npad, config->deconvolution_method);


  //waveform combiner.  (15 ants, npad=3 is default, filtered, deconvolved,responses)
  UCorrelator::WaveformCombiner *wfcomb = new UCorrelator::WaveformCombiner(15,config->combine_npad,true,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_xpol = new UCorrelator::WaveformCombiner(15,config->combine_npad,true,true,responses);

  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_xpol_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);


  wfcomb->setDelayToCenter(config->delay_to_center);
  wfcomb_xpol->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setDelayToCenter(config->delay_to_center);
  wfcomb_xpol_filtered->setDelayToCenter(config->delay_to_center);

  wfcomb->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_xpol->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_filtered->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_xpol_filtered->setGroupDelayFlag(config->enable_group_delay); 

  //Template Stuff
  //The length every waveform I look at from now on will be length, until I window it I guess.
  const int length = 2048;

  AnitaTemplateMachine *templateMachine = new AnitaTemplateMachine(length);
  templateMachine->loadTemplates();
  templateMachine->deconvolveTemplates(apd);

  //need some stuff for the snr thing
  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 







  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */

  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */

  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */
  //Loop through all the data!
  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  int skipped = 0;
  for (int entry=startEntry; entry<stopEntry; entry++) {
    //printout
    int printEntry = entry-startEntry;
    if (printEntry%100==0 && printEntry!=0) {
      cout << printEntry << "/" << lenEntries << " | ";
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(printEntry)/totalTimeSec;
      double processedFrac = 1.-float(skipped)/printEntry;
      double remaining = ((float(lenEntries-printEntry)*processedFrac)/rate)/60.;
      cout << std::setprecision(4) << 1./timeElapsed << "Hz <" << rate << ">";
      cout << "(skipped: " << skipped << " - " << processedFrac*100 << "%) ";
      cout << remaining << " minutes left, " << totalTimeSec/60. << " minutes elapsed" << endl;
      watch.Start();
    }


    summaryTree->GetEntry(entry);

    int eventNumber = evSum->eventNumber;

    int entryCheck = data->getEvent(eventNumber);
    if (entryCheck < 0) {
      cout << "Whoops I couldn't find eventNumber " << eventNumber << ".  Gotta fix me! Saving and Quitting..." << endl;
      outTree->Write();
      outFile->Close();
      return -1;
    }

    //things with shitty hwTrigAngles are still in here, so if this is one just skip it
    if (TMath::Abs(evSum->peak[0][0].hwAngle) > 45) {
      skipped++;
      continue;
    }


    //otherwise filter it
    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    //combine the info
    AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;//whatever just define it
    wfcomb->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    wfcomb_xpol->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);

    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    wfcomb_xpol_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);
    
    //grab the waveforms
    AnalysisWaveform *coherent = new AnalysisWaveform(*wfcomb->getCoherent()); 
    AnalysisWaveform *coherent_xpol = new AnalysisWaveform(*wfcomb_xpol->getCoherent()); 
    AnalysisWaveform *deconvolved = new AnalysisWaveform(*wfcomb->getDeconvolved()); 
    AnalysisWaveform *deconvolved_xpol = new AnalysisWaveform(*wfcomb_xpol->getDeconvolved()); 

    AnalysisWaveform *coherent_filtered = new AnalysisWaveform(*wfcomb_filtered->getCoherent()); 
    AnalysisWaveform *coherent_filtered_xpol = new AnalysisWaveform(*wfcomb_xpol_filtered->getCoherent()); 
    AnalysisWaveform *deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    AnalysisWaveform *deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_xpol_filtered->getDeconvolved()); 



    //redo stokes

    //define stuff
    AnalysisWaveform *wf;
    AnalysisWaveform *wf_xpol;
    AnitaEventSummary::WaveformInfo *info;

    //coherent
    wf = coherent;
    wf_xpol = coherent_xpol;
    info = &evSum->coherent[0][0];
    FFTtools::stokesParameters(wf->even()->GetN(),
                               wf->even()->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               wf_xpol->even()->GetY(), 
                               wf_xpol->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 

    //deconvolved
    wf = deconvolved;
    wf_xpol = deconvolved_xpol;
    info = &evSum->deconvolved[0][0];
    FFTtools::stokesParameters(wf->even()->GetN(),
                               wf->even()->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               wf_xpol->even()->GetY(), 
                               wf_xpol->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 

    //coherent filtered
    wf = coherent_filtered;
    wf_xpol = coherent_filtered_xpol;
    info = &evSum->coherent_filtered[0][0];
    FFTtools::stokesParameters(wf->even()->GetN(),
                               wf->even()->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               wf_xpol->even()->GetY(), 
                               wf_xpol->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 

    //deconvolved filtered
    wf = deconvolved_filtered;
    wf_xpol = deconvolved_filtered_xpol;
    info = &evSum->deconvolved_filtered[0][0];
    FFTtools::stokesParameters(wf->even()->GetN(),
                               wf->even()->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               wf_xpol->even()->GetY(), 
                               wf_xpol->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 




    //write over the snr with the better snr value
    double max,min;

    double rms = calcRMS(geom,ap,evSum->peak[0][0].phi,noiseSum);

    if (rms != 0) {
      max = TMath::MaxElement(coherent->even()->GetN(),coherent->even()->GetY());
      min = TMath::MinElement(coherent->even()->GetN(),coherent->even()->GetY());
      evSum->coherent[0][0].snr = ((max-min)/2.)/rms;
      
      max = TMath::MaxElement(coherent_filtered->even()->GetN(),coherent_filtered->even()->GetY());
      min = TMath::MinElement(coherent_filtered->even()->GetN(),coherent_filtered->even()->GetY());
      evSum->coherent_filtered[0][0].snr = ((max-min)/2.)/rms;
      
      max = TMath::MaxElement(deconvolved->even()->GetN(),deconvolved->even()->GetY());
      min = TMath::MinElement(deconvolved->even()->GetN(),deconvolved->even()->GetY());
      evSum->deconvolved[0][0].snr = ((max-min)/2.)/rms;
      
      max = TMath::MaxElement(deconvolved_filtered->even()->GetN(),deconvolved_filtered->even()->GetY());
      min = TMath::MinElement(deconvolved_filtered->even()->GetN(),deconvolved_filtered->even()->GetY());
      evSum->deconvolved_filtered[0][0].snr = ((max-min)/2.)/rms;
    }
    else {
      cout << "rms = 0" << endl;
    }


    //template stuff
    templateSummary->zeroInternals();

    //waveform,poli,dir,summary
    templateMachine->doTemplateAnalysis(coherent,0,0,templateSummary);    
    templateMachine->doDeconvolvedTemplateAnalysis(deconvolved,apd,0,0,templateSummary);
    
      
     
    //cool done with that then, save it
    outFile->cd();
    outTree->Fill();

    //delete stuff
    delete filtered;
    delete coherent;
    delete coherent_filtered;
    delete deconvolved;
    delete deconvolved_filtered;
    delete coherent_xpol;
    delete coherent_filtered_xpol;
    delete deconvolved_xpol;
    delete deconvolved_filtered_xpol;


  }

  cout << "all done looping! saving..." << endl;

  outFile->cd();
  outTree->Write();
  outFile->Close();

  
  cout << "Saved!  Goodbye!" << endl;

  return 1;

}


