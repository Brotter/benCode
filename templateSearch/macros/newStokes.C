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
#include "Polarimetry.h"

/*

  I want to re-calculate the stokes parameters because they are surely broken


  10GS/s
  smoothing (2ns, 500MHz LP) - peter found most power between 300 and 450MHz
  dedispersed

  


  This might do it :)

 */

#include "loadAll.C"


void reCalcStokes(polarimetry::StokesAnalysis *stokesAnalysis, double &Iout, double &Qout, double &Uout, double &Vout) {
  /*

    Get the integral of the peak of the instantaneous stokes parameters

   */


  Iout = 0;
  Qout = 0;
  Uout = 0;
  Vout = 0;

  //cool, lets take the integral of the useful region as well
  TGraph instI = stokesAnalysis->instI();
  TGraph instQ = stokesAnalysis->instQ();
  TGraph instU = stokesAnalysis->instU();
  TGraph instV = stokesAnalysis->instV();
  int N = instI.GetN();
  
  Long64_t locMaxI = TMath::LocMax(N,instI.GetY());
  double maxI = TMath::MaxElement(N,instI.GetY());
  
  int i = locMaxI;
  //  cout << "middle:" << i << " ";
  while ( instI.GetY()[i]/maxI > 0.25) {
    Iout += instI.GetY()[i];
    Qout += instQ.GetY()[i];
    Uout += instU.GetY()[i];
    Vout += instV.GetY()[i];
    i++;
  }
  //  cout << "upperEdge:" << i << " ";
  i = locMaxI-1;
  while ( instI.GetY()[i]/maxI > 0.25) {
    Iout += instI.GetY()[i];
    Qout += instQ.GetY()[i];
    Uout += instU.GetY()[i];
    Vout += instV.GetY()[i];
    i--;
    
  }
  //  cout << "upperEdge:" << i << endl;

  return;
}


void newStokes(string inFileName,string outBaseName="", int numSplits=1, int splitNum=0, bool copyAll=true) {
  /*
    
    Generate from instantaneous stokes 20%-height integral


    Input: Give it a summary tree file name (_ALL_ will do all ~3M quality events)
    Output: It will output a new one with a "_newStokes" at the end that has recalculated values

    opt:
    numSplits: number of splits to divide the data up into in case you want to cluster run it
    splitNum: number of the split this specific instance will run (starting at zero)
    copyAll: copy the whole summaryTree analysis file?  or just make a "AnitaEventSummary::WaveformSummary *newStokes" friendable branch
    outBaseName: The output file name base, which will have a "_#.root" if you are splitting
   */

  stringstream name;

  TChain *summaryTree;

  if (inFileName=="_ALL_"){
    //good events that have been fixed once (template)
    summaryTree = loadReKey(false);
  }
  
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(inFileName.c_str());
  }

  int totEntries = summaryTree->GetEntries();
  if (!totEntries) {
    cout << "No entries in that tree!" << endl;
    return;
  }
  cout << "found " << totEntries << " events to redo" << endl;



  //I'll need to split this up onto the servers, because it takes forever
  int startEntry,stopEntry,lenEntries;
  if (numSplits == 1) {
    startEntry=0;
    stopEntry=totEntries;
    lenEntries=totEntries;
  }
  else {
    lenEntries = totEntries/numSplits;
    startEntry = splitNum*lenEntries;
    stopEntry = (splitNum+1)*lenEntries;
    if (splitNum==numSplits-1) {
      stopEntry = totEntries;
      lenEntries = totEntries-startEntry;
    }
    cout << "Splitting into " << numSplits << " sections, which means " << lenEntries << " events for this section" << endl;
    cout << "Doing section: " << splitNum << ", starting at entry " << startEntry << " and stopping at " << stopEntry << endl;
  }


  
  //open an output file
  name.str("");
  size_t pos = inFileName.find(".root");
  string baseName = inFileName.substr(0,pos);
  if (outBaseName=="") { name << baseName << "_newStokes"; }
  name << outBaseName;
  if (numSplits > 1) { name << "_" << splitNum; }
  name << ".root";
  cout << "Using output file " << name.str() << endl;
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
  summaryTree->SetBranchAddress("template",&templateSummary);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  if (copyAll) {
    outTree->Branch("eventSummary",&evSum);
    outTree->Branch("template",&templateSummary);
    outTree->Branch("noiseSummary",&noiseSum);
    outTree->Branch("gpsEvent",&gps);
  }

  //new stuff
  AnitaEventSummary::WaveformInfo *newTop = new AnitaEventSummary::WaveformInfo();
  AnitaEventSummary::WaveformInfo *newMid = new AnitaEventSummary::WaveformInfo();
  AnitaEventSummary::WaveformInfo *newBot = new AnitaEventSummary::WaveformInfo();
  AnitaEventSummary::WaveformInfo *newAll = new AnitaEventSummary::WaveformInfo();
  //  outTree->Branch("newStokesTop",&newTop);
  //  outTree->Branch("newStokesMid",&newMid);
  //  outTree->Branch("newStokesBot",&newBot);
  outTree->Branch("newStokes",&newAll);
  

  
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


  //waveform combiner.  (15 ants, npad=3 is default, filtered, deconvolved,responses)
  UCorrelator::WaveformCombiner *wfcomb = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
									    true,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
										     false,true,responses);
  wfcomb->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setDelayToCenter(config->delay_to_center);

  wfcomb->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_filtered->setGroupDelayFlag(config->enable_group_delay); 



  //  filtering strat
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

  SimplePassBandFilter *lpFilt = new SimplePassBandFilter(0.2,0.5);
  fStrat->addOperation(lpFilt);

  //lets save each ring's polarimetry separately, so make some masks
  ULong64_t allowedTop = 0ul;
  for (int i=0; i<16; i++) {
    allowedTop |= (1ul << i);
  }
  ULong64_t disallowedTop = ~allowedTop;

  ULong64_t allowedMid = 0ul;
  for (int i=0; i<16; i++) {
    allowedMid |= (1ul << (i+16));
  }
  ULong64_t disallowedMid = ~allowedMid;

  ULong64_t allowedBot = 0ul;
  for (int i=0; i<16; i++) {
    allowedBot |= (1ul << (i+32));
  }
  ULong64_t disallowedBot = ~allowedBot;

  cout << "bitMasks:" << endl;
  cout << std::bitset<64>(disallowedTop) << endl;
  cout << std::bitset<64>(allowedTop) << endl;
  cout << std::bitset<64>(disallowedMid) << endl;
  cout << std::bitset<64>(allowedMid) << endl;
  cout << std::bitset<64>(disallowedBot) << endl;
  cout << std::bitset<64>(allowedBot) << endl;





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
    
    data->getEvent(evSum->eventNumber);

    //filter that waveform
    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    //just define this because I use it a lot
    AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;


    AnalysisWaveform *deconvolved_filtered, *deconvolved_filtered_xpol;
    polarimetry::StokesAnalysis *stokesAnalysis;

    //all antennas
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);
    deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved());     

    stokesAnalysis = new polarimetry::StokesAnalysis(deconvolved_filtered,deconvolved_filtered_xpol);
    reCalcStokes(stokesAnalysis,newAll->I,newAll->Q,newAll->U,newAll->V);

    delete deconvolved_filtered;
    delete deconvolved_filtered_xpol;
    delete stokesAnalysis;

    /*
    //top antennas
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal,disallowedTop);
    deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical,disallowedTop);
    deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved());     
    stokesAnalysis = new polarimetry::StokesAnalysis(deconvolved_filtered,deconvolved_filtered_xpol);
    reCalcStokes(stokesAnalysis,newTop->I,newTop->Q,newTop->U,newTop->V);

    delete deconvolved_filtered;
    delete deconvolved_filtered_xpol;
    delete stokesAnalysis;


    //mid antennas
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal,disallowedMid);
    deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical,disallowedMid);
    deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved());     
    stokesAnalysis = new polarimetry::StokesAnalysis(deconvolved_filtered,deconvolved_filtered_xpol);
    reCalcStokes(stokesAnalysis,newMid->I,newMid->Q,newMid->U,newMid->V);

    delete deconvolved_filtered;
    delete deconvolved_filtered_xpol;
    delete stokesAnalysis;

    
    //bot antennas
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal,disallowedBot);
    deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical,disallowedBot);
    deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved());     
    stokesAnalysis = new polarimetry::StokesAnalysis(deconvolved_filtered,deconvolved_filtered_xpol);
    reCalcStokes(stokesAnalysis,newBot->I,newBot->Q,newBot->U,newBot->V);
    
    delete deconvolved_filtered;
    delete deconvolved_filtered_xpol;
    delete stokesAnalysis;
    */


    outFile->cd();
    outTree->Fill();
    
    delete filtered;
  
  }

  outFile->cd();
  outTree->Write();
  outFile->Close();


  return;
}




/*testing stuff */

void allocateStokes(int N,double** Iins,double** Qins,double** Uins,double** Vins) {
  /*

    Allocate a bunch of huge buckets of memory

  */
  		  
  *Iins = (double*)calloc(N,sizeof(double));
  *Qins = (double*)calloc(N,sizeof(double));
  *Uins = (double*)calloc(N,sizeof(double));
  *Vins = (double*)calloc(N,sizeof(double));


  cout << "stokes allocated to N=" << N << endl;
  return;
}

void deleteStokes(double** Iins,double** Qins,double** Uins,double** Vins) {
		  
  /*
    
    deletes all these things

   */
            
  free(*Iins);
  free(*Qins);
  free(*Uins);
  free(*Vins);
}


void debugDrawStokes(string name, int N, double* Iins,double *Qins,double* Uins,double* Vins) {
		     
  /*

    Creates a new canvas and draws it with all these parameters

   */
  string gName;

  TGraph *gIins = new TGraph();
  gName = "g"+name+"_Iins";
  gIins->SetName(gName.c_str());
  gName = "Instantaneous "+name;
  gIins->SetTitle(gName.c_str());
  TGraph *gQins = new TGraph();
  gName = "g"+name+"_Qins";
  gQins->SetName(gName.c_str());
  TGraph *gUins = new TGraph();
  gName = "g"+name+"_Uins";
  gUins->SetName(gName.c_str());
  TGraph *gVins = new TGraph();
  gName = "g"+name+"_Vins";
  gVins->SetName(gName.c_str());
  
  gQins->SetLineColor(kBlue);
  gUins->SetLineColor(kRed);
  gVins->SetLineColor(kGreen);


  
  for (int pt=0; pt<N; pt++) {

    gIins->SetPoint(pt,pt,Iins[pt]);
    gQins->SetPoint(pt,pt,Qins[pt]);
    gUins->SetPoint(pt,pt,Uins[pt]);
    gVins->SetPoint(pt,pt,Vins[pt]);
  }

  string cName = "c"+name;
  TCanvas *c1 = new TCanvas(cName.c_str(),cName.c_str(),1000,500);
  gIins->Draw("al");
  gQins->Draw("lsame");
  gUins->Draw("lsame");
  gVins->Draw("lsame");
  TLegend *leg2 = new TLegend(0.1,0.7,0.3,0.9);
  leg2->AddEntry(gIins,"I","l");
  leg2->AddEntry(gQins,"Q","l");
  leg2->AddEntry(gUins,"U","l");
  leg2->AddEntry(gVins,"V","l");
  leg2->Draw();

  return;
}


void testStokesFromSummaryTree(string inFileName,bool debugDraw=false) {
  /*

    Messing around with calculating it by hand or using Cosmin's new polarimetry routine

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


  //maybe mask out the top since it is flipped and has weird stokes?
  ULong64_t allowed = 0ul;
  for (int i=0; i<16; i++) {
    allowed |= (1ul << i);
  }
  ULong64_t disallowed = ~allowed;
  //disallowed is the last argument to WaveformCombiner::combine();


  //waveform combiner.  (15 ants, npad=3 is default, filtered, deconvolved,responses)
  UCorrelator::WaveformCombiner *wfcomb = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
									    true,true,responses);
  UCorrelator::WaveformCombiner *wfcomb_filtered = new UCorrelator::WaveformCombiner(combine_nantennas,config->combine_npad,
										     false,true,responses);


  wfcomb->setDelayToCenter(config->delay_to_center);
  wfcomb_filtered->setDelayToCenter(config->delay_to_center);

  wfcomb->setGroupDelayFlag(config->enable_group_delay); 
  wfcomb_filtered->setGroupDelayFlag(config->enable_group_delay); 



  //for wais pulses specifically
  TH2D *h2Wais = new TH2D("h2Wais","wais plane of polarization",1000,0,120000,100,-20,20);



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
      return -1;
    }

    //filter that waveform
    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    //combine the channels and grab the waveforms
    AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;//whatever just define it

    wfcomb->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal,disallowed);
    AnalysisWaveform *coherent = new AnalysisWaveform(*wfcomb->getCoherent()); 
    AnalysisWaveform *deconvolved = new AnalysisWaveform(*wfcomb->getDeconvolved()); 
    wfcomb->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical,disallowed);
    AnalysisWaveform *coherent_xpol = new AnalysisWaveform(*wfcomb->getCoherent()); 
    AnalysisWaveform *deconvolved_xpol = new AnalysisWaveform(*wfcomb->getDeconvolved());

    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal,disallowed);
    AnalysisWaveform *coherent_filtered = new AnalysisWaveform(*wfcomb_filtered->getCoherent()); 
    AnalysisWaveform *deconvolved_filtered = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved()); 
    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical,disallowed);
    AnalysisWaveform *coherent_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getCoherent()); 
    AnalysisWaveform *deconvolved_filtered_xpol = new AnalysisWaveform(*wfcomb_filtered->getDeconvolved());     



    //redo stokes
    //  This is where things change a lot.  The FFTtools::stokesParameters() function has a TON of functionality I want to test
    //    cout << "redoing stokes" << endl;
    
    //define stuff
    AnalysisWaveform *wf;
    AnalysisWaveform *wf_xpol;
    double Iavg,Qavg,Uavg,Vavg;
    double *Iins, *Qins, *Uins, *Vins;
    int N;

    //coherent
    cout << "coherent" << endl;
    wf = coherent;
    wf_xpol = coherent_xpol;
    N = wf->even()->GetN();
    allocateStokes(N, &Iins,&Qins,&Uins,&Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               &Iavg,&Qavg,&Uavg,&Vavg,
                               Iins,Qins,Uins,Vins, false);
			       
    cout << "computed stokes" << endl;
    if (debugDraw) debugDrawStokes("coherent", N, Iins,Qins,Uins,Vins);

    deleteStokes(&Iins,&Qins,&Uins,&Vins);

    //deconvolved
    cout << "deconvolved" << endl;
    wf = deconvolved;
    wf_xpol = deconvolved_xpol;
    N = wf->even()->GetN();
    allocateStokes(N, &Iins,&Qins,&Uins,&Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               &Iavg,&Qavg,&Uavg,&Vavg,
                               Iins,Qins,Uins,Vins, false);
			       
    if (debugDraw) debugDrawStokes("deconvolved", N, Iins,Qins,Uins,Vins);

    deleteStokes(&Iins,&Qins,&Uins,&Vins);


    //coherent_filtered
    cout << "coherent_filtered" << endl;
    wf = coherent_filtered;
    wf_xpol = coherent_filtered_xpol;
    N = wf->even()->GetN();
    allocateStokes(N, &Iins,&Qins,&Uins,&Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               &Iavg,&Qavg,&Uavg,&Vavg,
                               Iins,Qins,Uins,Vins, false);
			       
    if (debugDraw) debugDrawStokes("coherent_filtered", N, Iins,Qins,Uins,Vins);

    deleteStokes(&Iins,&Qins,&Uins,&Vins);


    //deconvolved_filtered
    cout << "deconvolved_filtered" << endl;
    wf = deconvolved_filtered;
    wf_xpol = deconvolved_filtered_xpol;
    N = wf->even()->GetN();
    allocateStokes(N, &Iins,&Qins,&Uins,&Vins);

    FFTtools::stokesParameters(N,
			       wf->even()->GetY(), wf->hilbertTransform()->even()->GetY(),
			       wf_xpol->even()->GetY(),wf_xpol->hilbertTransform()->even()->GetY(), 
                               &Iavg,&Qavg,&Uavg,&Vavg,
                               Iins,Qins,Uins,Vins, false);
			       
    if (debugDraw) debugDrawStokes("deconvolved_filtered", N, Iins,Qins,Uins,Vins);

    deleteStokes(&Iins,&Qins,&Uins,&Vins);

    



    //lets try it with Cosmin's new Polarimetry.h thing too

    polarimetry::StokesAnalysis *stokesAnalysis = new polarimetry::StokesAnalysis(deconvolved_filtered,deconvolved_filtered_xpol);

    if (debugDraw) {
      TCanvas *cInstStokes = new TCanvas("cInstStokes","cInstStokes",1000,500);
      stokesAnalysis->instGraphs().Draw("a pmc plc") ;
      //      cInstStokes->BuildLegend(0.7,0.7,0.9,0.9,"","l");
      
      
      TCanvas *cCumuStokes = new TCanvas("cCumuStokes","cCumuStokes",1000,500);
      stokesAnalysis->cumuGraphs().Draw("a pmc plc") ;
      //      cCumuStokes->BuildLegend(0.7,0.7,0.9,0.9,"","l");
    }


    //cool, lets take the integral of the useful region as well
    TGraph *instI = new TGraph(stokesAnalysis->instI());
    TGraph *instQ = new TGraph(stokesAnalysis->instQ());
    TGraph *instU = new TGraph(stokesAnalysis->instU());
    TGraph *instV = new TGraph(stokesAnalysis->instV());
    N = instI->GetN();

    Long64_t locMaxI = TMath::LocMax(N,instI->GetY());
    double maxI = TMath::MaxElement(N,instI->GetY());

    double integralI=0;
    double integralQ=0;
    double integralU=0;
    double integralV=0;
    
    for (int i=0; i<N; i++) {
      if ( instI->GetY()[i]/maxI > 0.2) {
	integralI += instI->GetY()[i];
	integralQ += instQ->GetY()[i];
	integralU += instU->GetY()[i];
	integralV += instV->GetY()[i];
      }
    }

    double polAngle = (TMath::ATan(integralU/integralQ)/2.)*TMath::RadToDeg();
    double polFrac =  TMath::Sqrt(pow(integralQ/integralI,2) + pow(integralU/integralI,2));

    //    cout << integralI << " " << integralQ << " " << integralU << " " << integralV << endl;
    cout << polAngle << " " << polFrac << endl;;

    
    h2Wais->Fill(entry,polAngle);
		 
		 
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

  h2Wais->Draw("colz");
      
  cout << "cool whatever" << endl;
  return;

}




void newStokes() {
  cout << "loaded newStokes.C" << endl;
  return;

}
