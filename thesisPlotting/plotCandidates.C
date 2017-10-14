#include "UsefulAdu5Pat.h"
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
#include "SystemResponse.h"
#include "UCUtil.h"
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"

/*


  Lets make nice plots of the candidates!


 */



void plotCandidates(string candidateFilename="/Users/brotter/anita16/benCode/templateSearch/macros/passHarshCandidates_oct11.root") {
  /*

    I want to plot the single channels that went into a waveform, the coherent sum, and the map

    So lets make that into three plots?
    1) Interferometric map
    2) Channels that go into it
    3) Coherent sum and spectral magnitude (with slope!)
    4) Values?

   */
  stringstream name;


  //need the candidates
  TFile *inFile = TFile::Open(candidateFilename.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  if (!summaryTree) {
    cout << "There wasn't a summaryTree in " << candidateFilename << "! Quitting." << endl;
    return;
  }
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  int lenCands = summaryTree->GetEntries();

  
  //also need the data for waveforms
  AnitaDataset *data = new AnitaDataset(130,false);
  

  //filtering strat
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");


  //some config stuff
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig();
  //set the response to my "single" response, which I can shift really easily
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;
  //and to make the previous one obsolete, make it calculate the offsets compared to some central point
  config->delay_to_center = true;

  
  //get the responses
  AnitaResponse::ResponseManager *responses = new AnitaResponse::ResponseManager(UCorrelator::AnalysisConfig::getResponseString(config->response_option),config->response_npad, config->deconvolution_method);
  
  
  //WaveformCombiners for coherent summed waveforms
  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_xpol_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);
  wfcomb_filtered->setDelayToCenter(config->delay_to_center);
  wfcomb_xpol_filtered->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_xpol_filtered->setGroupDelayFlag(config->enable_group_delay); 


  TH2D** corrMaps = (TH2D**)malloc(sizeof(TH2D*)*lenCands);
  TGraphAligned** coherent = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** coherent_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** deconvolved = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** deconvolved_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  
  UCorrelator::Correlator *corr = new UCorrelator::Correlator(360,0,360,100,-60,40);

  
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  //  lenCands = 2;
  
  for (int event=0; event<lenCands; event++) {
    summaryTree->GetEntry(event);

    int eventNumber = evSum->eventNumber;
    data->getEvent(eventNumber);
    cout << "Candidate " << event << " - eventNumber " << eventNumber << endl;



    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    corr->compute(filtered,AnitaPol::kHorizontal);
    corrMaps[event] = new TH2D(*corr->getHist());

    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    wfcomb_xpol_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);

    name.str("");
    name << "ev" << eventNumber;
    coherent[event] = new TGraphAligned(*wfcomb_filtered->getCoherent()->even());
    coherent[event]->SetTitle(name.str().c_str());
    coherent_xPol[event] = new TGraphAligned(*wfcomb_xpol_filtered->getCoherent()->even());
    coherent_xPol[event]->SetTitle(name.str().c_str());
    deconvolved[event] = new TGraphAligned(*wfcomb_filtered->getDeconvolved()->even());
    deconvolved[event]->SetTitle(name.str().c_str());
    deconvolved_xPol[event] = new TGraphAligned(*wfcomb_xpol_filtered->getDeconvolved()->even());
    deconvolved_xPol[event]->SetTitle(name.str().c_str());	

    delete filtered;
  }    

  //everything should start at zero, maybe that will make this work?
  for (int i=0; i<lenCands; i++) {
    double startX = coherent[i]->GetX()[0];
    for (int pt=0;pt<coherent[i]->GetN();pt++) coherent[i]->GetX()[pt] -= startX;
    startX = deconvolved[i]->GetX()[0];
    for (int pt=0;pt<deconvolved[i]->GetN();pt++) deconvolved[i]->GetX()[pt] -= startX;
  }

  double deltaT = coherent[0]->GetX()[1]-coherent[0]->GetX()[0];
  cout << "deltaT=" << deltaT << endl;

  //I need to shift them so they align (by peak I guess?)
  for (int i=1; i<lenCands; i++) {
    TGraph *gCorr = FFTtools::getCorrelationGraph(coherent[0],coherent[i]);
    //    int peakBin = FFTtools::getPeakBin(gCorr);
    //    Int_t offsetPt = peakBin-(gCorr->GetN()/2);
    int offsetPt = FFTtools::getPeakBin(coherent[0]) - FFTtools::getPeakBin(coherent[i]);
    double offset = offsetPt*deltaT;
    delete gCorr;
    cout << "coherent"<< i << " " << offsetPt << " " << offset << endl;
    for (int pt=0;pt<coherent[i]->GetN();pt++) coherent[i]->GetX()[pt] += offset;
  }


  //I need to shift them so they align
  for (int i=1; i<lenCands; i++) {
    TGraph *gCorr = FFTtools::getCorrelationGraph(deconvolved[0],deconvolved[i]);
    //    int peakBin = FFTtools::getPeakBin(gCorr);
    //    Int_t offsetPt = peakBin-(gCorr->GetN()/2);
    int offsetPt = FFTtools::getPeakBin(deconvolved[0]) - FFTtools::getPeakBin(deconvolved[i]);
    double offset = offsetPt*deltaT;
    delete gCorr;
    cout << "deconvolved" << i << " " << offsetPt << " " << offset << endl;
    for (int pt=0;pt<deconvolved[i]->GetN();pt++) deconvolved[i]->GetX()[pt] += offset;
  }



  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(1,2);

  c1->cd(1);
  for (int i=0; i<lenCands; i++) {
    if (i==0)coherent[i]->Draw("al");
    else coherent[i]->Draw("lSame");
  }

  c1->cd(2);
  for (int i=0; i<lenCands; i++) {
    if (i==0)deconvolved[i]->Draw("al");
    else deconvolved[i]->Draw("lSame");
  }


  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->Divide(5,3);
  for (int i=0; i<lenCands; i++) {
    c2->cd(i+1);
    deconvolved[i]->GetXaxis()->SetRangeUser(0,40);
    if (i==2) deconvolved[i]->SetLineColor(kRed);
    deconvolved[i]->Draw("al");
  }
  


}
