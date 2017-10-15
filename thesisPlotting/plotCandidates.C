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



void plotCandidates(string candidateFilename="/Users/brotter/anita16/benCode/templateSearch/macros/candidatesHarsh_oct13.root") {
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
  FilterStrategy *fNone = new FilterStrategy();

  //some config stuff
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig();
  //set the response to my "single" response, which I can shift really easily
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;
  //and to make the previous one obsolete, make it calculate the offsets compared to some central point
  //  config->delay_to_center = true;

  
  //get the responses
  AnitaResponse::ResponseManager *responses = new AnitaResponse::ResponseManager(UCorrelator::AnalysisConfig::getResponseString(config->response_option),config->response_npad, config->deconvolution_method);
  
  
  //WaveformCombiners for coherent summed waveforms
  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_xpol_filtered = new UCorrelator::WaveformCombiner(15,config->combine_npad,false,true,responses);
  //  wfcomb_filtered->setDelayToCenter(config->delay_to_center);
  //  wfcomb_xpol_filtered->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setGroupDelayFlag(false);
  wfcomb_xpol_filtered->setGroupDelayFlag(false); 


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



    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fNone, data->gps(), data->header());

    corr->compute(filtered,AnitaPol::kHorizontal);
    corrMaps[event] = new TH2D(*corr->getHist());

    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    if (eventNumber != 15717147) {
      wfcomb_xpol_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);}
    else {
      wfcomb_xpol_filtered->combine(evSum->peak[pol][0].phi-8, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);}

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

    coherent_xPol[event]->SetLineStyle(3);
    deconvolved_xPol[event]->SetLineStyle(3);

    //highlight ev15717147
    if (eventNumber == 15717147) {
      coherent[event]->SetLineColor(kRed);
      coherent_xPol[event]->SetLineColor(kBlue);

      deconvolved[event]->SetLineColor(kRed);
      deconvolved_xPol[event]->SetLineColor(kBlue);
    }

    delete filtered;
  }    



  double deltaT = coherent[0]->GetX()[1]-coherent[0]->GetX()[0];
  cout << "deltaT=" << deltaT << endl;
  
  //Plot all of the deconvolved waveforms separately
  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->Divide(4,4);
  for (int i=0; i<lenCands; i++) {
    c2->cd(i+1);
    deconvolved[i]->GetXaxis()->SetRangeUser(0,25);
    deconvolved[i]->Draw("al");
    deconvolved_xPol[i]->Draw("lSame");
  }
  
  //Plot all of the coherent waveforms separately too
  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->Divide(4,4);
  for (int i=0; i<lenCands; i++) {
    c3->cd(i+1);
    //    coherent[i]->GetXaxis()->SetRangeUser(0,25);
    coherent[i]->Draw("al");
    coherent_xPol[i]->Draw("lSame");
  }


  TCanvas *cTau = new TCanvas("cTau","cTau",1000,600);
  cTau->Divide(1,2);
  cTau->cd(1);
  coherent[2]->Draw("al");
  coherent_xPol[2]->Draw("lSame");
  cTau->cd(2);
  deconvolved[2]->Draw("al");
  deconvolved_xPol[2]->Draw("lSame");
    

  TCanvas *cX = new TCanvas("cX","cX",1000,600);
  cX->Divide(1,2);
  cX->cd(1);
  coherent[8]->Draw("al");
  coherent_xPol[8]->Draw("lSame");
  cX->cd(2);
  deconvolved[8]->Draw("al");
  deconvolved_xPol[8]->Draw("lSame");
  


  //I need to shift them so they align to overlay them (by correlation)
  TGraphAligned** aligned_coherent = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** aligned_coherent_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** aligned_deconvolved = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** aligned_deconvolved_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  for (int i=0; i<lenCands; i++) {
    aligned_coherent[i] = (TGraphAligned*)coherent[i]->Clone();
    aligned_coherent_xPol[i] = (TGraphAligned*)coherent_xPol[i]->Clone();
    aligned_deconvolved[i] = (TGraphAligned*)deconvolved[i]->Clone();
    aligned_deconvolved_xPol[i] = (TGraphAligned*)deconvolved_xPol[i]->Clone();
  }

  //everything should start at zero, maybe that will make this work?
  for (int i=0; i<lenCands; i++) {
    double startX = aligned_coherent[i]->GetX()[0];
    for (int pt=0;pt<aligned_coherent[i]->GetN();pt++) aligned_coherent[i]->GetX()[pt] -= startX;
    for (int pt=0;pt<aligned_coherent_xPol[i]->GetN();pt++) aligned_coherent_xPol[i]->GetX()[pt] -= startX;
    startX = aligned_deconvolved[i]->GetX()[0];
    for (int pt=0;pt<aligned_deconvolved[i]->GetN();pt++) aligned_deconvolved[i]->GetX()[pt] -= startX;
    for (int pt=0;pt<aligned_deconvolved_xPol[i]->GetN();pt++) aligned_deconvolved_xPol[i]->GetX()[pt] -= startX;
  }
  //correlate and align coherent
  for (int i=1; i<lenCands; i++) {
    TGraph *gCorr = FFTtools::getCorrelationGraph(aligned_coherent[0],aligned_coherent[i]);
    if (i==2) for (int pt=0; pt<gCorr->GetN(); pt++) gCorr->GetY()[pt] *= -1; //tau
    int peakBin = FFTtools::getPeakBin(gCorr);
    Int_t offsetPt = peakBin-(gCorr->GetN()/2);
    double offset = offsetPt*deltaT;
    delete gCorr;
    cout << "aligned_coherent"<< i << " " << offsetPt << " " << offset << endl;
    for (int pt=0;pt<aligned_coherent[i]->GetN();pt++) aligned_coherent[i]->GetX()[pt] += offset;
    for (int pt=0;pt<aligned_coherent_xPol[i]->GetN();pt++) aligned_coherent_xPol[i]->GetX()[pt] += offset;
  }
  //correlate and align deconvolved
  for (int i=1; i<lenCands; i++) {
    TGraph *gCorr = FFTtools::getCorrelationGraph(aligned_deconvolved[0],aligned_deconvolved[i]);
    if (i==2) for (int pt=0; pt<gCorr->GetN(); pt++) gCorr->GetY()[pt] *= -1; //tau
    int peakBin = FFTtools::getPeakBin(gCorr);
    Int_t offsetPt = peakBin-(gCorr->GetN()/2);
    double offset = offsetPt*deltaT;
    delete gCorr;
    cout << "aligned_deconvolved" << i << " " << offsetPt << " " << offset << endl;
    for (int pt=0;pt<aligned_deconvolved[i]->GetN();pt++) aligned_deconvolved[i]->GetX()[pt] += offset;
    for (int pt=0;pt<aligned_deconvolved_xPol[i]->GetN();pt++) aligned_deconvolved_xPol[i]->GetX()[pt] += offset;
  }



  //Overlay all of the waveforms
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(1,2);
  // aligned_coherent on top
  c1->cd(1);
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_coherent[i]->Draw("al");
    else aligned_coherent[i]->Draw("lSame");
  }
  // deconvolved on bottom
  c1->cd(2);
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_deconvolved[i]->Draw("al");
    else aligned_deconvolved[i]->Draw("lSame");
  }
  


}
