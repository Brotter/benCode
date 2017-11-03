#include "Analyzer.h"

/*

  I want to re-calculate the stokes parameters because they are surely broken


  This might do it :)

 */


void allocateStokes(int N, double* Iavg,double *Qavg,double* Uavg,double* Vavg, double* Iins,double *Qins,double* Uins,double* Vins) {
  /*

    Allocate a bunch of huge buckets of memory

  */

  Iavg = new double[N];
  Qavg = new double[N];
  Uavg = new double[N];
  Vavg = new double[N];

  Iins = new double[N];
  Qins = new double[N];
  Uins = new double[N];
  Vins = new double[N];

  return;
}

void deleteStokes(double* Iavg,double *Qavg,double* Uavg,double* Vavg, double* Iins,double *Qins,double* Uins,double* Vins) {
  /*
    
    deletes all these things

   */

  if (Iavg) delete [] Iavg;
  if (Qavg) delete [] Qavg;
  if (Uavg) delete [] Uavg;
  if (Vavg) delete [] Vavg;

  if (Iins) delete [] Iins;
  if (Qins) delete [] Qins;
  if (Uins) delete [] Uins;
  if (Vins) delete [] Vins;
}


void debugDrawStokes(string name, int N, 
		     double* Iavg,double *Qavg,double* Uavg,double* Vavg, double* Iins,double *Qins,double* Uins,double* Vins) {
  /*

    Creates a new canvas and draws it with all these parameters

   */
  string gName;

  TGraph *gIavg = new TGraph();
  gName = "g"+name+"_Iavg";
  gIavg->SetName(gName.c_str());
  gIavg->SetTitle("Average");
  TGraph *gQavg = new TGraph();
  gName = "g"+name+"_Qavg";
  gQavg->SetName(gName.c_str());
  TGraph *gUavg = new TGraph();
  gName = "g"+name+"_Uavg";
  gUavg->SetName(gName.c_str());
  TGraph *gVavg = new TGraph();
  gName = "g"+name+"_Vavg";
  gVavg->SetName(gName.c_str());

  TGraph *gIins = new TGraph();
  gName = "g"+name+"_Iins";
  gIins->SetName(gName.c_str());
  gIins->SetTitle("Instantaneous");
  TGraph *gQins = new TGraph();
  gName = "g"+name+"_Qins";
  gQins->SetName(gName.c_str());
  TGraph *gUins = new TGraph();
  gName = "g"+name+"_Uins";
  gUins->SetName(gName.c_str());
  TGraph *gVins = new TGraph();
  gName = "g"+name+"_Vins";
  gVins->SetName(gName.c_str());


  gQavg->SetMarkerColor(kBlue);
  gUavg->SetMarkerColor(kRed);
  gVavg->SetMarkerColor(kGreen);
  
  gQins->SetMarkerColor(kBlue);
  gUins->SetMarkerColor(kRed);
  gVins->SetMarkerColor(kGreen);


  
  for (int pt=0; pt<N; pt++) {
    gIavg->SetPoint(pt,pt,Iavg[pt]);
    gQavg->SetPoint(pt,pt,Qavg[pt]);
    gUavg->SetPoint(pt,pt,Uavg[pt]);
    gVavg->SetPoint(pt,pt,Vavg[pt]);

    gIins->SetPoint(pt,pt,Iins[pt]);
    gQins->SetPoint(pt,pt,Qins[pt]);
    gUins->SetPoint(pt,pt,Uins[pt]);
    gVins->SetPoint(pt,pt,Vins[pt]);
  }

  string cName = "c"+name;
  TCanvas *c1 = new TCanvas(cName.c_str(),cName.c_str(),1000,500);
  c1->Divide(1,2);

  c1->cd(1);
  gIavg->Draw("al");
  gQavg->Draw("lsame");
  gUavg->Draw("lsame");
  gVavg->Draw("lsame");
  TLegend *leg1 = new TLegend(0.2,0.7,0.4,0.8);
  leg1->AddEntry(gIavg,"I","l");
  leg1->AddEntry(gQavg,"Q","l");
  leg1->AddEntry(gUavg,"U","l");
  leg1->AddEntry(gVavg,"V","l");
  leg1->Draw();

  c1->cd(2);
  gIins->Draw("al");
  gQins->Draw("lsame");
  gUins->Draw("lsame");
  gVins->Draw("lsame");
  TLegend *leg2 = new TLegend(0.5,0.7,0.7,0.8);
  leg2->AddEntry(gIins,"I","l");
  leg2->AddEntry(gQins,"Q","l");
  leg2->AddEntry(gUins,"U","l");
  leg2->AddEntry(gVins,"V","l");
  leg2->Draw();

  return;
}


void redoStokesFromSummaryTree(string inFileName,bool debugDraw=false) {
  /*

    Input: Give it a summary tree file name
    Output: It will output a new one with a "_newStokes" at the end that has recalculated values

    opt: 
    debugDraw - will draw a stokes picture for each plot
   */


  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(inFileName.c_str());
  int lenEntries = summaryTree->GetEntries();
  if (!lenEntries) {
    cout << "No entries in that tree!" << endl;
    return;
  }
  cout << "found " << lenEntries << " events to redo" << endl;

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);

  
  //  need all the full root data too.  Start at 130, no decimation, and no blinding strat
  AnitaDataset *data = new AnitaDataset(130,false);

  //  some config stuff
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig();
  //set the response to my "individual" response, since that appropriately handles all the channels
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseIndividualBRotter;
  //all pass deconvolution
  AnitaResponse::AllPassDeconvolution *apd = new AnitaResponse::AllPassDeconvolution();
  config->deconvolution_method = apd;
  //and to make the previous one obsolete, make it calculate the offsets compared to some central point
  config->delay_to_center = true;
  //set it to 15 antennas
  const int combine_nantennas = 15;
  config->combine_nantennas = combine_nantennas;
  //get the responses
  AnitaResponse::ResponseManager *responses = new AnitaResponse::ResponseManager(UCorrelator::AnalysisConfig::getResponseString(config->response_option),config->response_npad, config->deconvolution_method);


  //  filtering strat
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");



  //waveform combiner.  (15 ants, npad=3 is default, filtered, deconvolved,responses)
  UCorrelator::WaveformCombiner *wfcomb = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
									    true,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
										     false,true,responses);


  wfcomb->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setDelayToCenter(config->delay_to_center);

  wfcomb->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_filtered->setGroupDelayFlag(config->enable_group_delay); 





  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */

  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */

  /*  LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START */
  //Loop through all the data!
  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%100)) cout << entry << " / " << lenEntries << endl;

    summaryTree->GetEntry(entry);

    int eventNumber = evSum->eventNumber;

    int entryCheck = data->getEvent(eventNumber);
    if (entryCheck < 0) {
      cout << "Whoops I couldn't find eventNumber " << eventNumber << ".  Gotta fix me! Saving and Quitting..." << endl;
      outTree->Write();
      outFile->Close();
      return -1;
    }

    //filter that waveform
    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    //combine the channels and grab the waveforms
    AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;//whatever just define it

    wfcomb->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    AnalysisWaveform *coherent = new AnalysisWaveform(*wfcomb->getCoherent()); 
    AnalysisWaveform *deconvolved = new AnalysisWaveform(*wfcomb->getDeconvolved()); 
    wfcomb->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);
    AnalysisWaveform *coherent_xpol = new AnalysisWaveform(*wfcomb_xpol->getCoherent()); 
    AnalysisWaveform *deconvolved_xpol = new AnalysisWaveform(*wfcomb_xpol->getDeconvolved());

    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    AnalysisWaveform *coherent_filtered = new AnalysisWaveform(*wfcomb_filtered->getCoherent()); 
    AnalysisWaveform *deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);
    AnalysisWaveform *coherent_filtered_xpol = new AnalysisWaveform(*wfcomb_xpol_filtered->getCoherent()); 
    AnalysisWaveform *deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_xpol_filtered->getDeconvolved());     



    //redo stokes
    //  This is where things change a lot.  The FFTtools::stokesParameters() function has a TON of functionality I want to test


    //define stuff
    AnalysisWaveform *wf;
    AnalysisWaveform *wf_xpol;
    double *Iavg, *Qavg, *Uavg, *Vavg;
    double *Iins, *Qins, *Uins, *Vins;
    int N;

    //coherent
    wf = coherent;
    wf_xpol = coherent_xpol;
    N = wf->even->GetN();
    allocateStokes(N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               Iavg,Qavg,Uavg,Vavg,
                               Iins,Qins,Uins,Vins);
			       
    if (debugDraw) debugDrawStokes("coherent", N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    deleteStokes(Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    //deconvolved
    wf = deconvolved;
    wf_xpol = deconvolved_xpol;
    N = wf->even->GetN();
    allocateStokes(N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               Iavg,Qavg,Uavg,Vavg,
                               Iins,Qins,Uins,Vins);
			       
    if (debugDraw) debugDrawStokes("deconvolved", N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    deleteStokes(Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);


    //coherent_filtered
    wf = coherent_filtered;
    wf_xpol = coherent_filtered_xpol;
    N = wf->even->GetN();
    allocateStokes(N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               Iavg,Qavg,Uavg,Vavg,
                               Iins,Qins,Uins,Vins);
			       
    if (debugDraw) debugDrawStokes("coherent_filtered", N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    deleteStokes(Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);


    //deconvolved_filtered
    wf = deconvolved_filtered;
    wf_xpol = deconvolved_filtered_xpol;
    N = wf->even->GetN();
    allocateStokes(N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               Iavg,Qavg,Uavg,Vavg,
                               Iins,Qins,Uins,Vins);
			       
    if (debugDraw) debugDrawStokes("deconvolved_filtered", N, Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);

    deleteStokes(Iavg,Qavg,Uavg,Vavg, Iins,Qins,Uins,Vins);



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


  cout << "cool whatever" << endl;
  return;

}

void newStokes() {

  cout << "loaded newStokes.C" << endl;
  return;

}
