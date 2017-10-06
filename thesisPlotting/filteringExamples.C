#include "/Users/brotter/anita16/local/bin/runMagicDisplayEv.C" //for the runToEv thing
#include "AnitaConventions.h"
#include "UCFilters.h"


void sineSubFiteringPlots(int eventNumber) {

  /*

    I need a plot that shows what sine subtraction does to a waveform

    MagicDisplay is nice, but it looks like MagicDisplay and isn't really thesis worthy

   */


  int run = evToRun(eventNumber);

  AnitaDataset *data = new AnitaDataset(run,false);
  data->getEvent(eventNumber);
  

  //Make a filter strategy
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");
  FilterStrategy *fStrat_none = UCorrelator::getStrategyWithKey("");
  
  
  /* CONFIGURATION FOR THE ANALYSIS!!!!!! */
  //and a configuration for the analysis
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 

  //set the response to my "single" response
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  //  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;

  //set the stokes windowing to use the entire waveform
  config->windowStokes = true;
  config->stokesWindowLength = 500; //in "points" which are 100ps I think, so this is 50ns
  
  //lets try to do only the first peak (I never use the others)
  config->nmaxima = 2;

  //lets also try to always get the bottom ring to be the "first" waveform in the coherent sum
  config->set_bottom_first = true;

  //and to make the previous one obsolete, make it calculate the offsets compared to some central point
  config->delay_to_center = true;

  //and maybe 5 phi sectors instead of 4
  config->combine_nantennas = 15;

  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); //interactive needs to be true

  AnitaEventSummary *eventSummary = new AnitaEventSummary();

  //sine sub filtered
  FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());
  analyzer->analyze(filteredEvent, eventSummary); 

  UCorrelator::Correlator *corr = analyzer->getCorrelator();
  TH2D *hist = (TH2D*)corr->getHist()->Clone();
  hist->SetTitle("Sine Subtraction Filtering;Payload Azimuth (degrees); Payload Elevation (degrees)");

  TH2D *hist_zoomed = (TH2D*)corr->computeZoomed(eventSummary->wais.phi,-eventSummary->wais.theta,100,0.2,100,0.2)->Clone();
  hist_zoomed->SetTitle("Zoomed Filtered;Payload Azimuth (degrees); Payload Elevation (degrees)");


  //no filtering
  FilteredAnitaEvent *filteredEvent_none = new FilteredAnitaEvent(data->useful(), fStrat_none, data->gps(), data->header());
  analyzer->clearInteractiveMemory();
  analyzer->analyze(filteredEvent_none, eventSummary); 

  corr = analyzer->getCorrelator();
  TH2D *hist_none = (TH2D*)corr->getHist()->Clone();
  hist_none->SetTitle("No Filtering;Payload Azimuth (degrees); Payload Elevation (degrees)");

  TH2D *hist_none_zoomed = (TH2D*)corr->computeZoomed(eventSummary->wais.phi,-eventSummary->wais.theta,100,0.2,100,0.2)->Clone();
  hist_none_zoomed->SetTitle("Zoomed Unfiltered;Payload Azimuth (degrees); Payload Elevation (degrees)");

  TGraph *wais = new TGraph();
  wais->SetPoint(0,eventSummary->wais.phi,-eventSummary->wais.theta);
  wais->SetMarkerStyle(29);//star
  wais->SetMarkerSize(1);

  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->Divide(2,2);
  c1->cd(1);
  hist->Draw("colz");
  wais->Draw("pSame");
  c1->cd(2);
  hist_zoomed->Draw("colz");
  wais->Draw("pSame");

  c1->cd(3);
  hist_none->Draw("colz");
  wais->Draw("pSame");
  c1->cd(4);
  hist_none_zoomed->Draw("colz");
  wais->Draw("pSame");

}


const FFTtools::SineSubtract* sineSubStepByStep(int eventNumber) {

  /*

    I need a plot that shows what sine subtraction does to a waveform

    I want to plot each iteration I guess

   */


  int run = evToRun(eventNumber);

  AnitaDataset *data = new AnitaDataset(run,false);
  data->getEvent(eventNumber);
  

  //Make a filter strategy
  UCorrelator::SineSubtractFilter *sineSubtract = new UCorrelator::SineSubtractFilter(0.1,3);
  sineSubtract->setInteractive(true);
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::SpectrumAverageLoader *specAvgLoader = new UCorrelator::SpectrumAverageLoader(specAvgDir);
  sineSubtract->makeAdaptive(specAvgLoader,2); 
  sineSubtract->setVerbose(true);

  FilterStrategy *fStrat = new FilterStrategy();
  fStrat->addOperation(sineSubtract);
  


  FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

  const FFTtools::SineSubtract *ss12 = sineSubtract->sinsub(AnitaPol::kHorizontal,13);
  
  //  ss12->makePlots();
  ss12->makeSlides("event56765803","slides","./","png",false);


  return ss12;

}

    

void filteringExamples() {
  cout << "loaded filteringExamples.C" << endl;
  //  sineSubFiteringPlots(56765803);
}
