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
#include "AnitaTemplates.h"


#include "loadAll.C"
/*


  Lets make nice plots of the candidates!


 */



void plotCandidates(string candidateFilename="/Users/brotter/anita16/benCode/templateSearch/macros/trueCandidates_oct14.root",
		    bool save = true, bool blind = false) {
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
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("template",&tempSum);

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


  TH2D** corrMapsH = (TH2D**)malloc(sizeof(TH2D*)*lenCands);
  TH2D** corrMapsV = (TH2D**)malloc(sizeof(TH2D*)*lenCands);
  TGraphAligned** coherent = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** coherent_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** deconvolved = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  TGraphAligned** deconvolved_xPol = (TGraphAligned**)malloc(sizeof(TGraphAligned*)*lenCands);
  
  UCorrelator::Correlator *corr = new UCorrelator::Correlator(360,0,360,100,-60,40);

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  vector<int> eventNumbers;
  vector<int> polarity;

  //  lenCands = 2;
  
  for (int event=0; event<lenCands; event++) {
    summaryTree->GetEntry(event);

    int eventNumber = evSum->eventNumber;
    eventNumbers.push_back(eventNumber);
    data->getEvent(eventNumber);
    cout << "Candidate " << event << " - eventNumber " << eventNumber << endl;
    
    polarity.push_back(tempSum->deconvolved[0][0].cRay_pol[4]);

    FilteredAnitaEvent *filtered = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

    corr->compute(filtered,AnitaPol::kHorizontal);
    corrMapsH[event] = new TH2D(*corr->getHist());

    corr->compute(filtered,AnitaPol::kVertical);
    corrMapsV[event] = new TH2D(*corr->getHist());



    wfcomb_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kHorizontal);
    wfcomb_xpol_filtered->combine(evSum->peak[pol][0].phi, evSum->peak[pol][0].theta, filtered, AnitaPol::kVertical);


    name.str("");
    name << "ev" << eventNumber << ";Time (ns); Measured Voltage (mV)";
    coherent[event] = new TGraphAligned(*wfcomb_filtered->getCoherent()->even());
    coherent[event]->SetTitle(name.str().c_str());
    if (tempSum->deconvolved[0][0].cRay_pol[4] && blind) {
      for (int pt=0; pt<coherent[event]->GetN(); pt++) { coherent[event]->GetY()[pt] *= -1.; } }
    coherent_xPol[event] = new TGraphAligned(*wfcomb_xpol_filtered->getCoherent()->even());
    coherent_xPol[event]->SetTitle(name.str().c_str());
    deconvolved[event] = new TGraphAligned(*wfcomb_filtered->getDeconvolved()->even());
    deconvolved[event]->SetTitle(name.str().c_str());
    if (tempSum->deconvolved[0][0].cRay_pol[4] && blind) {
      for (int pt=0; pt<deconvolved[event]->GetN(); pt++) { deconvolved[event]->GetY()[pt] *= -1.; } }
    deconvolved_xPol[event] = new TGraphAligned(*wfcomb_xpol_filtered->getDeconvolved()->even());
    deconvolved_xPol[event]->SetTitle(name.str().c_str());	

    coherent_xPol[event]->SetLineColor(kBlue);
    deconvolved_xPol[event]->SetLineColor(kBlue);

    //highlight ev15717147
    if (eventNumber == 15717147 && !blind) {
      coherent[event]->SetLineColor(kRed);
      coherent_xPol[event]->SetLineColor(kBlue);

      deconvolved[event]->SetLineColor(kRed);
      deconvolved_xPol[event]->SetLineColor(kBlue);
    }


    //highlight above horizon
    if (eventNumber == 39599205 || eventNumber == 27142546) {
      coherent[event]->SetLineColor(kGreen+3);
      coherent_xPol[event]->SetLineColor(kBlue);

      deconvolved[event]->SetLineColor(kGreen+3);
      deconvolved_xPol[event]->SetLineColor(kBlue);
    }

    delete filtered;
  }    



  double deltaT = coherent[0]->GetX()[1]-coherent[0]->GetX()[0];
  cout << "deltaT=" << deltaT << endl;
  
  //Plot all of the deconvolved waveforms separately
  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->Divide(4,5);
  for (int i=0; i<lenCands; i++) {
    c2->cd(i+1);
    deconvolved[i]->GetXaxis()->SetRangeUser(-5,20);
    deconvolved[i]->Draw("al");
    deconvolved_xPol[i]->Draw("lSame");
  }
  
  //Plot all of the coherent waveforms separately too
  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->Divide(4,5);
  for (int i=0; i<lenCands; i++) {
    c3->cd(i+1);
    //    coherent[i]->GetXaxis()->SetRangeUser(0,25);
    coherent[i]->Draw("al");
    coherent_xPol[i]->Draw("lSame");
  }


  
  for (int i=0; i<eventNumbers.size(); i++) {
    int eventNumber = eventNumbers[i];
    TCanvas *cX = new TCanvas("cX","cX",1000,600);
    cX->Divide(1,2);
    cX->cd(1);
    name.str(""); name << "Ev" << eventNumber << ";Time (ns); Voltage (mv)";
    coherent[i]->SetTitle(name.str().c_str());
    coherent[i]->Draw("al");
    coherent_xPol[i]->Draw("lSame");
    cX->cd(2);
    name.str(""); name << ";Time (ns); Voltage (mv)";
    deconvolved[i]->Draw("al");
    deconvolved_xPol[i]->Draw("lSame");
    
    name.str(""); name << "Ev" << eventNumber << "_waveform.png";
    cX->SaveAs(name.str().c_str());
    delete cX;
  }


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
    if (polarity[i]) for (int pt=0; pt<gCorr->GetN(); pt++) gCorr->GetY()[pt] *= -1; //inverted pol
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
    if (polarity[i]) for (int pt=0; pt<gCorr->GetN(); pt++) gCorr->GetY()[pt] *= -1; //inverted pol
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
  aligned_coherent[0]->SetTitle("Coherently Summed Waveform; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_coherent[i]->Draw("al");
    else aligned_coherent[i]->Draw("lSame");
  }
  // deconvolved on bottom
  c1->cd(2);
  aligned_deconvolved[0]->SetTitle("De-Dispersed Waveform; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_deconvolved[i]->Draw("al");
    else aligned_deconvolved[i]->Draw("lSame");
  }
  
  //Also the Vpol
  TCanvas *c1_xpol = new TCanvas("c1_xpol","c1_xpol",1000,600);
  c1_xpol->Divide(1,2);
  // aligned_coherent on top
  c1_xpol->cd(1);
  aligned_coherent_xPol[0]->SetTitle("Coherently Summed Waveform, V-Pol; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_coherent_xPol[i]->Draw("al");
    else aligned_coherent_xPol[i]->Draw("lSame");
  }
  // deconvolved on bottom
  c1_xpol->cd(2);
  aligned_deconvolved_xPol[0]->SetTitle("De-Dispersed Waveform, V-Pol; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (i==0)aligned_deconvolved_xPol[i]->Draw("al");
    else aligned_deconvolved_xPol[i]->Draw("lSame");
  }
  
  //overlay the direct events
  TCanvas *c1_direct = new TCanvas("c1_direct","c1_direct",1000,600);
  c1_direct->Divide(1,2);
  // aligned_coherent on top
  c1_direct->cd(1);
  aligned_coherent[0]->SetTitle("Coherently Summed Waveform, H-Pol; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (eventNumbers[i]==15717147) {aligned_coherent[i]->Draw("al"); aligned_coherent[i]->GetYaxis()->SetRangeUser(-150,150); }
    else if (eventNumbers[i]==27142546 || eventNumbers[i]==39599205) aligned_coherent[i]->Draw("lSame");
  }
  // deconvolved on bottom
  c1_direct->cd(2);
  aligned_deconvolved[0]->SetTitle("De-Dispersed Waveform, H-Pol; Time (ns); Measured Voltage (mV)");
  for (int i=0; i<lenCands; i++) {
    if (eventNumbers[i]==15717147) { aligned_deconvolved[i]->Draw("al"); aligned_deconvolved[i]->GetYaxis()->SetRangeUser(-150,150); }
    else if (eventNumbers[i]==27142546 || eventNumbers[i]==39599205) aligned_deconvolved[i]->Draw("lSame");
  }
  




  /*    Interferometric Maps     */
  for (int i=0; i<lenCands; i++) {
    name.str("");
    name << "mapCan" << i;
    TCanvas *mapCan = new TCanvas(name.str().c_str(),name.str().c_str(),1000,600);
    mapCan->Divide(2);
    TPaveText *title = new TPaveText(0.4,0.90, 0.6,0.95);
    name.str(""); name << "Event " << eventNumbers[i];
    title->AddText(name.str().c_str());
    title->Draw();
    
    TVirtualPad *p1 = mapCan->cd(1);
    p1->SetPad(0,0,0.5,0.9);
    corrMapsH[i]->SetTitle("Horizontal Polarization; Payload Phi (degrees); Elevation (degrees)");
    corrMapsH[i]->SetStats(0);
    corrMapsH[i]->Draw("colz");
    TVirtualPad *p2 = mapCan->cd(2);
    p2->SetPad(0.5,0,1,0.9);
    corrMapsV[i]->SetTitle("Vertical Polarization; Payload Phi (degrees); Elevation (degrees)");
    corrMapsV[i]->SetStats(0);
    corrMapsV[i]->Draw("colz");


    name.str(""); name << "intMap_ev" << eventNumbers[i] << ".png";
    if (save) mapCan->SaveAs(name.str().c_str());


    delete mapCan;
    
  }


  /*    Interferometric Maps  but just HPol   */
  for (int i=0; i<lenCands; i++) {
    name.str("");
    name << "mapCan" << i;
    TCanvas *mapCan = new TCanvas(name.str().c_str(),name.str().c_str(),1000,600);
     TPaveText *title = new TPaveText(0.4,0.90, 0.6,0.95);
    name.str(""); name << "Event " << eventNumbers[i];
    title->AddText(name.str().c_str());
    title->Draw();
 
    name.str("");
    name << "Event " << eventNumbers[i] << "; Payload Phi (degrees); Elevation (degrees)";
    corrMapsH[i]->SetTitle(name.str().c_str());
    corrMapsH[i]->SetStats(0);
    corrMapsH[i]->Draw("colz");

    name.str(""); name << "intMap_Hpol_ev" << eventNumbers[i] << ".png";
    if (save) mapCan->SaveAs(name.str().c_str());


    delete mapCan;
    
  }





}



void saveLocalDistributions(string date="10.16.17_01h") {
  /*

    Saves the 2d histograms for cut parameters of events near each candidate, with the candidate highlighted

   */
  stringstream name;

  gStyle->SetOptStat(0);

  char* dataDir = getenv("ANITA3_RESULTSDIR");
  name.str(""); name << dataDir << "/cluster/" << date << "/";
  TChain *summaryTree = loadWhatever(name.str(),"clusterBackground",64,false);

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  summaryTree->BuildIndex("eventNumber");

  //to get the event numbers
  TFile *candTreeFile = TFile::Open("/home/brotter/anita16/benCode/templateSearch/macros/trueCandidates_oct14.root");
  TTree *candTree = (TTree*)candTreeFile->Get("summaryTree");
  AnitaEventSummary *candSum = NULL;
  candTree->SetBranchAddress("eventSummary",&candSum);

  TGraph *temp = new TGraph();
  temp->SetMarkerStyle(4);
  temp->SetMarkerSize(4);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

  int entry;
  for (int cand=0; cand<candTree->GetEntries(); cand++) {
    candTree->GetEntry(cand);
    int eventNumber = candSum->eventNumber;

    c1->Clear();
    c1->cd();
    name.str("");
    name << "ev" << eventNumber << ";Deconvoled Template Correlation; Coherent Sum Template Correlation";
    TH2D* hTempVsTemp = new TH2D("hTempVsTemp",name.str().c_str(),100,0,1,100,0,1);
    name.str(""); name << "seedEventNumber == " << eventNumber;
    summaryTree->Draw("template.coherent[0][0].cRay[4]:template.deconvolved[0][0].cRay[4] >> hTempVsTemp",name.str().c_str(),"colz");

    entry = summaryTree->GetEntryNumberWithBestIndex(eventNumber);
    if (entry>0) summaryTree->GetEntry(entry);
    else cout << "whoops, couldn't find event number " << eventNumber << endl;
    temp->SetPoint(0,tempSum->deconvolved[0][0].cRay[4],tempSum->coherent[0][0].cRay[4]);
    temp->Draw("p same");

    name.str(""); name << "ev" << eventNumber << "_tempVsTemp.png";
    c1->SaveAs(name.str().c_str());

    delete hTempVsTemp;

    c1->Clear();
    c1->cd();
    name.str("");
    name << "ev" << eventNumber << ";Map Peak; Map SNR";
    TH2D* hMap = new TH2D("hMap",name.str().c_str(),100,0,0.3,100,0,25);
    name.str(""); name << "seedEventNumber == " << eventNumber;
    summaryTree->Draw("peak[0][0].snr:peak[0][0].value >> hMap",name.str().c_str(),"colz");

    entry = summaryTree->GetEntryNumberWithBestIndex(eventNumber);
    if (entry>0) summaryTree->GetEntry(entry);
    else cout << "whoops, couldn't find event number " << eventNumber << endl;
    temp->SetPoint(0,evSum->peak[0][0].value,evSum->peak[0][0].snr);
    temp->Draw("p same");

    name.str(""); name << "ev" << eventNumber << "_interfMap.png";
    c1->SaveAs(name.str().c_str());

    delete hMap;

    c1->Clear();
    c1->cd();
    name.str("");
    name << "ev" << eventNumber << ";Coherent Sum Hilbert Peak; Coherent Linear Polarization Fraction";
    TH2D* hHilbLinFrac = new TH2D("hHilbLinFrac",name.str().c_str(),100,0,130,100,0,1);
    name.str(""); name << "seedEventNumber == " << eventNumber;
    summaryTree->Draw("coherent_filtered[0][0].linearPolFrac():coherent_filtered[0][0].peakHilbert >> hHilbLinFrac",name.str().c_str(),"colz");

    entry = summaryTree->GetEntryNumberWithBestIndex(eventNumber);
    if (entry>0) summaryTree->GetEntry(entry);
    else cout << "whoops, couldn't find event number " << eventNumber << endl;
    temp->SetPoint(0,evSum->coherent_filtered[0][0].peakHilbert,evSum->coherent_filtered[0][0].linearPolFrac());
    temp->Draw("p same");

    name.str(""); name << "ev" << eventNumber << "_hilbLinFrac.png";
    c1->SaveAs(name.str().c_str());
    
    delete hHilbLinFrac;

  }

  return;
}

