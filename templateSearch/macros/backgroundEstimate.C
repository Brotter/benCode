/*

  I have some background estimate code that should live here

  I was living in cluster.C for whatever reason.

 */
#include "Normalize.C" //makeNormCumulative
#include "AnitaConventions.h"
#include "GeoMagnetic.h"

#include "loadAll.C"

void drawQuadrent(TH1D *hBlackIn, TH1D *hRedIn) {

  TH1 *hRed = makeNormCumulative(hRedIn,false);
  TH1 *hBlack   = makeNormCumulative(hBlackIn);

  hBlack->SetStats(false);
  hRed->SetStats(false);

  hBlack->GetYaxis()->SetRangeUser(1e-2,1);
  hRed->GetYaxis()->SetRangeUser(1e-2,1);

  hRed->GetYaxis()->SetTitleColor(kRed);
  hRed->SetLineColor(kRed);
  hBlack->SetLineColor(kBlack);

  TPad *p1 = new TPad("p1","",0,0,1,1);
  p1->SetLogy();
  TPad *p2 = new TPad("p2","",0,0,1,1);
  p2->SetLogy();
  p2->SetFillStyle(4000);
  p2->SetFillColor(0);
  p2->SetFrameFillStyle(0);
  p1->Draw();
  p1->cd();
  hBlack->Draw();
  p2->Draw();
  p2->cd();
  hRed->Draw("Y+");

}

void drawBaseClusterHists() {
  /*	
	I made all those histograms for base clusters, so lets plot them 
  */
  stringstream name;
  
  const int numBases = 104; //thats just how many there are                                                                                     
  TCanvas *c1 = new TCanvas("c1","",1000,600);
  

  TFile *inFile = TFile::Open("mergeBaseClusters.root");

  for (int base=0; base<numBases; base++) {
    name.str("");
    name << "hCluster_"<< base;
    TH1D* currHist = (TH1D*)inFile->Get(name.str().c_str());

    if (currHist == NULL) continue;

    if (currHist->Integral() == 0) continue;
    
    currHist->Draw();
    
    name.str("");
    name << "baseDists/base" << base << ".png";
    c1->SaveAs(name.str().c_str());
  }
}


void drawBaseDistributionsWithCandidates(bool cumulative=false,bool doCut=false) {
  /*

    I want to see what the distributions of known base pointed events is vs the distributions of candidates

   */


  //get base pointed event summaries
  //  TFile *baseFile = TFile::Open("mergeBaseClusters.root");
  //  TTree *baseTree = (TTree*)baseFile->Get("summaryTree");
  TChain *baseTree = loadWhatever("/home/brotter/anita16/benCode/templateSearch/macros/clusteringOutputs/10.16.17_01h/backgroundClusterABCD",64,false);
  //get candidate event summaries
  TFile *candidateFile = TFile::Open("trueCandidates_oct14.root");
  TTree *candidateTree = (TTree*)candidateFile->Get("summaryTree");
  candidateTree->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TVirtualPad *pad;


  TCut allCuts;
  if (doCut) allCuts = "peak[0][0].value > 0.0435 && peak[0][0].snr > 9.05 && deconvolved_filtered[0][0].peakHilbert > 47.5 && template.coherent[0][0].cRay[4] > 0.666 && flags.pulser==0";
  //if (doCut) allCuts = "peak[0][0].value > 0.0695 && peak[0][0].snr > 13.35 && template.coherent[0][0].cRay[4] > 0.485 && deconvolved_filtered[0][0].peakHilbert > 37.5 && flags.pulser==0";
  //  if (doCut) allCuts = "peak[0][0].value > 0.0055 && peak[0][0].snr > 1.35 && template.coherent[0][0].cRay[4] > 0.115 && deconvolved_filtered[0][0].peakHilbert > 11.5 && flags.pulser==0";
  //  else allCuts = "flags.pulser==0 && flags.maxBottomToTopRatio < 3";
  //point at bases (from cluster.C)                                                                                        
  else allCuts = "flags.pulser==0 && flags.maxBottomToTopRatio[0] < 3 && TMath::Abs(peak[0][0].hwAngle) < 45 && eventNumber != 84405480 && eventNumber != 84114142";
  cout << "1/4" << endl;
  TH1D *hMapPeak = new TH1D("hMapPeak","Interferometric Map Peak;Interferometric Map Peak;count",250,0,0.5);
  TH1D *hMapPeakCand = new TH1D("hMapPeakCand","Interferometric Map Peak;Interferometric Map Peak;count",250,0,0.5);
  baseTree->Draw("peak[0][0].value >> hMapPeak",allCuts);
  int numCands = candidateTree->Draw("peak[0][0].value >> hMapPeakCand",allCuts,"same");
  cout << "numCands:" << numCands << endl;
  cout << "2/4" << endl;
  TH1D *hLinPolFrac = new TH1D("hLinPolFrac","Linear Polarization;Linear Polarization Fraction;count",250,0,1);
  TH1D *hLinPolFracCand = new TH1D("hLinPolFracCand","Linear Polarization;Linear Polarization Fraction;count",250,0,1);
  baseTree->Draw("coherent_filtered[0][0].linearPolFrac() >> hLinPolFrac",allCuts);
  candidateTree->Draw("coherent_filtered[0][0].linearPolFrac() >> hLinPolFracCand",allCuts,"same");
  cout << "3/4" << endl;
  TH1D *hPeakHilbertDF = new TH1D("hPeakHilbertDF","Deconvolved Hilbert Peak;Deconvolved Hilbert Peak;count",250,0,500);
  TH1D *hPeakHilbertDFCand = new TH1D("hPeakHilbertDFCand","Deconvolved Hilbert Peak;Deconvolved Hilbert Peak;count",250,0,500);
  baseTree->Draw("deconvolved_filtered[0][0].peakHilbert >> hPeakHilbertDF",allCuts);
  candidateTree->Draw("deconvolved_filtered[0][0].peakHilbert >> hPeakHilbertDFCand",allCuts,"same");
  cout << "4/4" << endl;
  TH1D *hTemplate = new TH1D("hTemplate","cRay +4 Correlation;cRay +4 Correlation;count",250,0,1);
  TH1D *hTemplateCand = new TH1D("hTemplateCand","cRay +4 Correlation;cRay +4 Correlation;count",250,0,1);
  baseTree->Draw("template.coherent[0][0].cRay[4] >> hTemplate",allCuts);
  candidateTree->Draw("template.coherent[0][0].cRay[4] >> hTemplateCand",allCuts,"same");

  hMapPeakCand->SetLineColor(kRed);
  hLinPolFracCand->SetLineColor(kRed);
  hPeakHilbertDFCand->SetLineColor(kRed);
  hTemplateCand->SetLineColor(kRed);

  c1->Clear();
  c1->Divide(2,2);

  pad = c1->cd(1);
  pad->SetLogy();
  if (cumulative) {
    drawQuadrent(hMapPeak,hMapPeakCand);
  }
  else {
    hMapPeak->Draw();
    hMapPeakCand->Draw("same"); }

  pad = c1->cd(2);
  pad->SetLogy();
  if (cumulative) {
    drawQuadrent(hLinPolFrac,hLinPolFracCand);
  }
  else {
    hLinPolFrac->Draw();
    hLinPolFracCand->Draw("same"); }  

  pad = c1->cd(3);
  pad->SetLogy();
  if (cumulative) {
    drawQuadrent(hPeakHilbertDF,hPeakHilbertDFCand);
  }
  else {
    hPeakHilbertDF->Draw();
    hPeakHilbertDFCand->Draw("same"); }

  pad = c1->cd(4);
  pad->SetLogy();
  if (cumulative) {
    drawQuadrent(hTemplate,hTemplateCand);
  }
  else {
    hTemplate->Draw();
    hTemplateCand->Draw("same"); }
   

  return;
}





void saveCandidatePValues(bool doCut=false,string fileName="") {

  /* 

     Basically, calculate the fraction of events above any given event in each distribution

     If they were all uncorrelated, Ptot = TT(Px)

     They aren't though, so you have to do some chain rule probability stuff


   */

  //get base pointed event summaries
  //  TFile *baseFile = TFile::Open("mergeBaseClusters.root");
  //  TTree *baseTree = (TTree*)baseFile->Get("summaryTree");
  TChain *baseTree = loadWhatever("/home/brotter/anita16/benCode/templateSearch/macros/clusteringOutputs/10.16.17_01h/backgroundClusterABCD",64,false);

  //get candidate event summaries
  TFile *candidateFile = TFile::Open("trueCandidates_oct14.root");
  TTree *candidateTree = (TTree*)candidateFile->Get("summaryTree");
  AnitaEventSummary *candSum = NULL;
  candidateTree->SetBranchAddress("eventSummary",&candSum);
  AnitaTemplateSummary *candTemp = NULL;
  candidateTree->SetBranchAddress("template",&candTemp);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TVirtualPad *pad;

  TCut allCuts;
  if (doCut) allCuts = "peak[0][0].value > 0.0435 && peak[0][0].snr > 9.05 && deconvolved_filtered[0][0].peakHilbert > 47.5 && template.coherent[0][0].cRay[4] > 0.666 && flags.pulser==0";
  //  if (doCut) allCuts = "peak[0][0].value > 0.0695 && peak[0][0].snr > 13.35 && template.coherent[0][0].cRay[4] > 0.485 && deconvolved_filtered[0][0].peakHilbert > 37.5 && flags.pulser==0";
  // if (doCut) allCuts = "peak[0][0].value > 0.0055 && peak[0][0].snr > 1.35 && template.coherent[0][0].cRay[4] > 0.115 && deconvolved_filtered[0][0].peakHilbert > 11.5 && flags.pulser==0";
  //  else allCuts = "flags.pulser==0 && flags.maxBottomToTopRatio < 3 && eventNumber != 11116669 && eventNumber != 11989349 && eventNumber != 16952229 && eventNumber != 33484995 && eventNumber != 58592863 && eventNumber != 62273732 && eventNumber != 63210848 && eventNumber != 80561103 && eventNumber != 83877990 && eventNumber != 84114142";
  else allCuts = "flags.pulser==0 && flags.maxBottomToTopRatio[0] < 3 && TMath::Abs(peak[0][0].hwAngle) < 45";
  
  
  cout << "1" << endl;
  TH1D *hMapPeak = new TH1D("hMapPeak","Interferometric Map Peak;Interferometric Map Peak;count",250,0,0.5);
  baseTree->Draw("peak[0][0].value >> hMapPeak",allCuts);
  cout << "2" << endl;
  TH1D *hLinPolFrac = new TH1D("hLinPolFrac","Linear Polarization;Linear Polarization Fraction;count",250,0,1);
  baseTree->Draw("coherent_filtered[0][0].linearPolFrac() >> hLinPolFrac",allCuts);
  cout << "3" << endl;
  TH1D *hPeakHilbert = new TH1D("hPeakHilbert","Coherent Hilbert Peak;Coherent Hilbert Peak;count",250,0,500);
  baseTree->Draw("coherent_filtered[0][0].peakHilbert >> hPeakHilbert",allCuts);
  cout << "4" << endl;
  TH1D *hTemplate = new TH1D("hTemplate","cRay +4 Correlation;cRay +4 Correlation;count",250,0,1);
  baseTree->Draw("template.coherent[0][0].cRay[4] >> hTemplate",allCuts);

  TH1 *hMapPeakC   = makeNormCumulative(hMapPeak);
  TH1 *hLinPolFracC   = makeNormCumulative(hLinPolFrac);
  TH1 *hPeakHilbertC   = makeNormCumulative(hPeakHilbert);
  TH1 *hTemplateC   = makeNormCumulative(hTemplate);

  TFile *outRoot = TFile::Open("pValues.root","recreate");

  hMapPeak->Write();
  hLinPolFrac->Write();
  hPeakHilbert->Write();
  hTemplate->Write();

  hMapPeakC->Write();
  hLinPolFracC->Write();
  hPeakHilbertC->Write();
  hTemplateC->Write();

  outRoot->Close();

  if (fileName=="") fileName="pValues.csv";
  ofstream outFile(fileName);
  outFile << "#eventNumber,mapPeakVal,PmapPeakVal,mapPeakSNR,PmapPeakSNR,peakHilbert,PpeakHilbert,cRayTemplate,PcRayTemplate" << endl;
  outFile << "#calculated using " << baseTree->GetEntries() << " events" << endl;

  for (int entry=0; entry<candidateTree->GetEntries(); entry++) {
    candidateTree->GetEntry(entry);
    
    bool skip=false;

    int bin;
    double mapPeakVal = candSum->peak[0][0].value;
    bin = hMapPeakC->GetXaxis()->FindBin(mapPeakVal);
    double PmapPeakVal = hMapPeakC->GetBinContent(bin);

    double linearPolFrac = candSum->coherent_filtered[0][0].linearPolFrac();
    bin = hLinPolFracC->GetXaxis()->FindBin(linearPolFrac);
    double PlinearPolFrac = hLinPolFracC->GetBinContent(bin);

    double peakHilbert = candSum->coherent_filtered[0][0].peakHilbert;
    bin = hPeakHilbertC->GetXaxis()->FindBin(peakHilbert);
    double PpeakHilbert = hPeakHilbertC->GetBinContent(bin);

    double cRayTemp = candTemp->coherent[0][0].cRay[4];
    bin = hTemplateC->GetXaxis()->FindBin(cRayTemp);
    double PcRayTemp = hTemplateC->GetBinContent(bin);


    //write event out to file
    outFile << candSum->eventNumber << ",";
    outFile << mapPeakVal << ",";
    outFile << PmapPeakVal << ",";
    outFile << linearPolFrac << ",";
    outFile << PlinearPolFrac << ",";
    outFile << peakHilbert << ",";
    outFile << PpeakHilbert << ",";
    outFile << cRayTemp << ",";
    outFile << PcRayTemp << ",";
    outFile <<  PmapPeakVal*PpeakHilbert*PcRayTemp << endl;
  }

  outFile.close();

}



void make2DCovarience(bool savePlots=false) {
  stringstream name,toPlot;

  //get base pointed event summaries
  TFile *baseFile = TFile::Open("mergeBaseClusters.root");
  TTree *baseTree = (TTree*)baseFile->Get("summaryTree");

  //things you tell to TTree::Draw() so that it draws
  string plotStrings[4] = {"peak[0][0].value","coherent_filtered.linearPolFrac()","deconvolved_filtered[0][0].peakHilbert","template.coherent[0][0].cRay[4]"};

  //if you want to save them, they should have descriptive names
  string plotNames[6] = {"mapPeak_linPol.png","mapPeak_hilb.png","mapPeak_temp.png","linPol_hilb.png","linPol_temp","hilb_temp.png"};

  string cutString = "flags.pulser==0 && flags.maxBottomToTopRatio < 3 && TMath::Abs(peak[0][0].hwAngle) < 45";

  TH2D *hCovariance = new TH2D("hCovariance","hCovariance",4,-0.5,3.5, 4,-0.5,3.5);

  int count=0;
  for (int i=0; i<4; i++) {
    for (int j=i+1; j<4; j++) {
      toPlot.str("");
      toPlot << plotStrings[i] << ":" << plotStrings[j];
      name.str("");
      name << "c" << i << j;
      TCanvas *c1 = new TCanvas(name.str().c_str(),"",1000,600);
      c1->SetLogz();
      TVirtualPad *pad1 = c1->cd(1);
      baseTree->Draw(toPlot.str().c_str(),cutString.c_str(),"colz");

      if (savePlots) c1->SaveAs(plotNames[count].c_str());

      TH2D* currHist = (TH2D*)pad1->GetPrimitive("htemp");
      double covar = currHist->GetCorrelationFactor(1,2);
      cout << toPlot.str() << " " << covar << endl;

      hCovariance->Fill(i,j,covar);

      count++;
    }
  }

  TCanvas *c1 = new TCanvas(name.str().c_str(),"",1000,600);
  hCovariance->Draw("colz");


}





void makeMinbiasBackgroundHist() {
  /*
    drawCanidateClusters needs background plots, gotta run this to make 'em
  */

  TFile *inFile = TFile::Open("07.28.17_17h_decimated.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't open file" << endl;
    return;
  }
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");

  int lenEntries = summaryTree->GetEntries();
  cout << "got " << lenEntries << " entries" << endl;

  TFile *outFile = TFile::Open("minbiasBackgrounds.root","recreate");

  cout << "1" << endl;
  TH2D* histMapPeakVsTemplate = new TH2D("histMapPeakVsTemplate","minbias",100,0,1,350,0,0.35);
  summaryTree->Draw("peak[0][0].value:template.coherent[0][0].cRay[4] >> histMapPeakVsTemplate","!flags.isRF","colz");
  histMapPeakVsTemplate->SetStats(0);
  histMapPeakVsTemplate->Write();

  cout << "2" << endl;
  TH2D* histMapSNRVsHilbert = new TH2D("histMapSNRVsHilbert","minbias",800,0,800,450,0,45);
  summaryTree->Draw("peak[0][0].snr:deconvolved_filtered[0][0].peakHilbert >> histMapSNRVsHilbert","!flags.isRF","colz");
  histMapSNRVsHilbert->SetStats(0);
  histMapSNRVsHilbert->Write();

  cout << "3" << endl;
  TH1D* histMapPeak = new TH1D("histMapPeak","minbias;Map Peak; Count",350,0,0.35);
  summaryTree->Draw("peak[0][0].value >> histMapPeak","!flags.isRF");
  histMapPeak->SetStats(0);
  histMapPeak->Write();

  cout << "4" << endl;
  TH1D* histMapSNR = new TH1D("histMapSNR","minbias;Map SNR; Count",450,0,45);
  summaryTree->Draw("peak[0][0].snr >> histMapSNR","!flags.isRF");
  histMapSNR->SetStats(0);
  histMapSNR->Write();

  cout << "5" << endl;
  TH1D* histTemplate = new TH1D("histTemplate","minbias; cRay +4 Template Correlation; Count",100,0,1);
  summaryTree->Draw("template.coherent[0][0].cRay[4] >> histTemplate","!flags.isRF");
  histTemplate->SetStats(0);
  histTemplate->Write();

  cout << "6" << endl;
  TH1D* histHilbert = new TH1D("histHilbert","minbias; Deconvolved Hilbert Peak; Count",800,0,800);
  summaryTree->Draw("deconvolved_filtered[0][0].peakHilbert >> histHilbert","!flags.isRF");
  histHilbert->SetStats(0);
  histHilbert->Write();

  outFile->Close();

  return;
}



TGraph* numberOfEventsExceedingCandidates() {
  /*

    Reads in the pseudoBaseCluster.root list, which includes all events from near events that pass cuts, and 
    returns the number of events that have reduced quantities higher than all the events in candidates.root

    Doesn't use TTree::Draw because you can only do that once per thing, so scanning over everything is faster

    Returns a TGraph with x=eventNumber and y=numberOfExceedingBackgroundEvs
   */

  TFile *candFile = TFile::Open("trueCandidates_oct14.root");
  TTree *candTree = (TTree*)candFile->Get("summaryTree");
  AnitaEventSummary *candEvSum = NULL;
  candTree->SetBranchAddress("eventSummary",&candEvSum);
  AnitaTemplateSummary *candTempSum = NULL;
  candTree->SetBranchAddress("template",&candTempSum);
  int numCandidates = candTree->GetEntries();
  cout << "found " << numCandidates << " candidates" << endl;


  //Cut values:
  vector<double> vTemplateCorr;
  vector<double> vMapPeak;
  vector<double> vHilbPeak;
  vector<double> vLinPolFrac;


  for (int i=0; i<numCandidates; i++) {
    candTree->GetEntry(i);
    vTemplateCorr.push_back(candTempSum->coherent[0][0].cRay[4]);
    vMapPeak.push_back(candEvSum->peak[0][0].value);
    vHilbPeak.push_back(candEvSum->coherent_filtered[0][0].peakHilbert);
    vLinPolFrac.push_back(candEvSum->coherent_filtered[0][0].linearPolFrac());
  }

  
  //  TFile *backFile = TFile::Open("pseudoBaseCluster.root");
  //  TTree *backTree = (TTree*)backFile->Get("summaryTree");
    TChain *backTree = loadWhatever("/home/brotter/anita16/benCode/templateSearch/macros/clusteringOutputs/10.16.17_01h/backgroundClusterABCD",64,false);
  AnitaEventSummary *backEvSum = NULL;
  backTree->SetBranchAddress("eventSummary",&backEvSum);
  AnitaTemplateSummary *backTempSum = NULL;
  backTree->SetBranchAddress("template",&backTempSum);
  
  
  int lenEntries = backTree->GetEntries();
  cout << "Found " << lenEntries << " background events from impulsive sources" << endl;
  
  int* numPassing = new int [numCandidates];
  for (int i=0; i<numCandidates; i++) {
    numPassing[i] = 0;
  }
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << endl;
    backTree->GetEntry(entry);
    
    double templateCorr = backTempSum->coherent[0][0].cRay[4];
    double mapPeak = backEvSum->peak[0][0].value;
    double hilbPeak = backEvSum->deconvolved_filtered[0][0].peakHilbert;
    double linPolFrac = backEvSum->coherent[0][0].linearPolFrac();
    
    for (int cand=0; cand<numCandidates; cand++) {
      if (vTemplateCorr[cand] < templateCorr &&
	  vMapPeak[cand] < mapPeak && 
	  vHilbPeak[cand] < hilbPeak &&
	  vLinPolFrac[cand] < linPolFrac) {
	numPassing[cand]++;
      }
    }
  }
  
  TGraph *outGraph = new TGraph();
  outGraph->SetName("gExceedingBkgdEvs");
  outGraph->SetTitle("Events Exceeding Candidates;eventNumber;# of events exceeding candidate values");
  for (int i=0; i<numCandidates; i++) {
    candTree->GetEntry(i);
    outGraph->SetPoint(i,candEvSum->eventNumber,numPassing[i]);
    cout << candEvSum->eventNumber << " " << numPassing[i] << endl;
  }

  

  delete [] numPassing;

  return outGraph;
}
  

void poissonConfidenceInterval() {
  /*

    One of the candidates has zero events that exceed it, so I need to use the Poisson Confidence Limit equation

    So this determines what the means of the distributions are from the pseudoBaseList I guess

   */

  
  TFile *backFile = TFile::Open("pseudoBaseCluster.root");
  TTree *backTree = (TTree*)backFile->Get("summaryTree");

  TH1D *htemp;

  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->Divide(2,2);
  TVirtualPad *p1 = c1->cd(1);
  backTree->Draw("template.coherent[0][0].cRay[4]");
  htemp = (TH1D*)p1->GetPrimitive("htemp");
  double templateCorr = htemp->GetMean();
  cout << "templateCorr=" << templateCorr << endl;
  
  TVirtualPad *p2 = c1->cd(2);  backTree->Draw("peak[0][0].value");
  htemp = (TH1D*)p2->GetPrimitive("htemp");
  double mapPeak = htemp->GetMean();
  cout << "mapPeak=" << mapPeak << endl;  

  TVirtualPad *p3 = c1->cd(3);  backTree->Draw("deconvolved_filtered[0][0].peakHilbert");
  htemp = (TH1D*)p3->GetPrimitive("htemp");
  double peakHilbert = htemp->GetMean();
  cout << "peakHilbert=" << peakHilbert << endl;

  TVirtualPad *p4 = c1->cd(4);  backTree->Draw("coherent[0][0].linearPolFrac()");
  htemp = (TH1D*)p4->GetPrimitive("htemp");
  double linPolFrac = htemp->GetMean();
  cout << "linPolFrac=" << linPolFrac << endl;

}



void eventsClusteredWithCandidates(string date = "10.05.17_14h") {

  return;
}


void ABCDMethod() {
  
  /*

    Okay so I need to do the ABCD method:
    A) Events that both pass cuts and cluster together into pseudobases
    B) Events that DO NOT pass cuts, but DO cluster
    C) Events that DO NOT pass cuts, and DO NOT cluster
    D) Events that pass cuts, but DO NOT cluster (candidates)

    I end up just doing this by hand with separateNotable

   */

  //load all the pseudobase events
  /*  TChain *pseudoTree = new TChain("summaryTree","summaryTree");

  stringstream name;
  string baseDir = "/Volumes/ANITA3Data/bigAnalysisFiles/cluster/10.03.17_22h/";
  for (int i=0; i<32; i++) {
    name.str("");
    name << baseDir << "pseudoBaseCluster_" << i << ".root";
    pseudoTree->Add(name.str().c_str());
  }
  int pseudoN = spTree->GetEntries();
  cout << "Found " << lenPseudo << " entries" << endl;

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);

  //how many of those cluster with the pseudobases (-1==WAIS and -2==LDB)
  //DO NOT pass, but DO cluster (B)
  int B = summaryTree->Draw("eventNumber","baseClusteredWith > 0","goff");
  
  */
  return;
  // 
}




void localUnnormalizedLikelihood(bool draw=false,string date="10.26.17_14h28m") {
  /*

    The ABCD method might be stupid.

    So, lets look at the distributions of events that fall near each of our candidates.

    To determine the "unnormalized likelihood", we can take the ratio of product of the mode of the noise values, to
    the product of the candidate values.  This will tell us how far "outside" the noise distributions the candidate is


   */
  stringstream name,name2;

  char* dataDir = getenv("ANITA3_RESULTSDIR");
  name.str(""); name << dataDir << "/cluster/" << date << "/clusterBackground";
  TChain *summaryTree = loadWhatever(name.str().c_str(),64);
  int lenEntries = summaryTree->GetEntries();
  if (!lenEntries) { cout << "No events found in that file, qutting." << endl; return; }
  else cout << "Found " << lenEntries << " events" << endl;

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("template",&tempSum);
  int seedEventNumber;
  summaryTree->SetBranchAddress("seedEventNumber",&seedEventNumber);


  //some vectors to store all the candidates you find
  //  This is just because I want this code to be as configurable as possible and driven by the backgroundCluster output
  //these should all be in the same order, because I don't want to learn maps
  vector<int> seedCandidates;
  std::vector<int>::iterator seedIt;  //need at least one iterator for the event number
  vector<TH1D*> productHists;
  vector<TH1D*> candHists;
  vector<TF1*> fitFuncs;
  //lets just make the antarctica maps here too
  vector<TH2DAntarctica*> backCountMaps;
  vector<TProfile2DAntarctica*> productMaps;
  vector<TGraphAntarctica*> candMaps;


  //also a TGraph because they are easy to grab single values from, easier than a TH1F
  TGraph *candGraph = new TGraph();
  candGraph->SetName("candGraph");

  //pointer for "current" histograms so it isn't always reallocated
  TH1D* currHist;
  TH2DAntarctica* currCountMap;
  TProfile2DAntarctica* currProductMap;
  TGraphAntarctica *currCandMap;

  //index of position in vectors
  int vectIndex;

  //Loop!
  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%10000)) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);


    //things that shouldn't have passed...
    if (TMath::Abs(evSum->peak[0][0].hwAngle) > 45) continue;
    if (evSum->flags.maxBottomToTopRatio[0] > 3) continue;
    if (!evSum->flags.isRF) continue;
    if (evSum->flags.isPayloadBlast) continue;

    //find the location in the vectors for that seed candidate, or make a new place for it
    seedIt = find(seedCandidates.begin(), seedCandidates.end(), seedEventNumber);
    //not found
    if (seedIt == seedCandidates.end()) {
      cout << "New eventNumber:" << seedEventNumber << endl;
      seedCandidates.push_back(seedEventNumber);

      name.str("");
      if (seedEventNumber > 0)        name << "ev" << seedEventNumber;
      else if (seedEventNumber == -1) name << "WAIS";
      else if (seedEventNumber == -2) name << "LDB";
      else                            name << "ev_negative" << seedEventNumber;

      currHist = new TH1D(name.str().c_str(),name.str().c_str(),1000,-5,5);
      productHists.push_back(currHist);
      
      candGraph->SetPoint(candGraph->GetN(),seedEventNumber,-999);

      name2.str("");
      name2 << name.str() << "_candidate";
      currHist = new TH1D(name2.str().c_str(),name2.str().c_str(),1000,-5,5);
      currHist->SetLineColor(kRed); //red :)
      candHists.push_back(currHist);

      name2.str("");
      name2 << name.str() << "_countMap";
      currCountMap = new TH2DAntarctica(name2.str().c_str(),name2.str().c_str(),1000,1000);
      backCountMaps.push_back(currCountMap);


      name2.str("");
      name2 << name.str() << "_productMap";
      currProductMap = new TProfile2DAntarctica(name2.str().c_str(),name2.str().c_str(),1000,1000);
      productMaps.push_back(currProductMap);


      name2.str("");
      name2 << name.str() << "_candMap";
      currCandMap = new TGraphAntarctica();
      currCandMap->SetName(name.str().c_str());
      candMaps.push_back(currCandMap);

      vectIndex = seedCandidates.size() - 1;

    }
    //already exists in vectors
    else {
      vectIndex = seedIt - seedCandidates.begin();
    }
    currHist = productHists[vectIndex];


    
    double product = 1.0;
    //    product *= evSum->peak[0][0].value;
    //    product *= evSum->peak[0][0].snr;
    product *= tempSum->coherent[0][0].cRay[4];
    product *= tempSum->deconvolved[0][0].cRay[4];
    //    product *= evSum->coherent_filtered[0][0].peakHilbert;
    //    product *= evSum->coherent_filtered[0][0].peakVal; //huh? why are you using this?
    product *= evSum->coherent_filtered[0][0].linearPolFrac();

    currHist->Fill(product);


    backCountMaps[vectIndex]->Fill(evSum->peak[0][0].longitude,evSum->peak[0][0].latitude);
    productMaps[vectIndex]->Fill(evSum->peak[0][0].longitude,evSum->peak[0][0].latitude,product);


    //if this is actually the seed event.  AKA the candidate! Assumes candidate will always seed itself
    if (seedEventNumber == evSum->eventNumber) {
      cout << "Found ev" << evSum->eventNumber << " with product: " << product << endl;
      name.str("");
      if (seedEventNumber > 0) name << "cand" << seedEventNumber;
      candHists[vectIndex]->Fill(TMath::Log10(product));
      candGraph->SetPoint(vectIndex,evSum->eventNumber,TMath::Log10(product));
      candMaps[vectIndex]->SetPoint(0,evSum->peak[0][0].longitude,evSum->peak[0][0].latitude);
      candMaps[vectIndex]->SetMarkerStyle(41);

    }

  }


      
  //write out all that stuff to a root file!
  //need to do this before the fit for whatever reason, otherwise it segfaults
  TFile *outFile = TFile::Open("localUnnormalizedLikelihood.root","recreate");
  for (int i=0;i<productHists.size();i++) {
    cout << "eventNumber:" << seedCandidates[i] << endl;
    productHists[i]->Write();
    candHists[i]->Write();
  }
  candGraph->Write();
  outFile->Close();

  /* Fitting and finding the final background estimate number (number of sigma outside of log normal fit */
  ofstream outTxtFile("localUnnormalizedLikelihood.txt");
  outTxtFile << "eventNumber candProduct distAmp distMean distSigma background" << endl;
  for (int i=0; i<productHists.size(); i++) {
    int eventNumber = candGraph->GetX()[i];

    name.str(""); name << "gausFit_" << eventNumber;
    TF1 *gausFit = new TF1(name.str().c_str(),"gaus(0)",-5,5);
    productHists[i]->Fit(gausFit);
    fitFuncs.push_back(gausFit);

    //gaussian has three parameters: [0]*exp(-0.5*((x-[1])/[2])**2)
    double amplitude = gausFit->GetParameter(0);//amplitude
    double mean = gausFit->GetParameter(1);//mean
    double sigma = gausFit->GetParameter(2);//sigma
    
    double candProduct = candGraph->GetY()[i];
    double background = amplitude*exp(-0.5*pow((candProduct-mean)/sigma,2));

    outTxtFile << eventNumber << " " << candProduct << " " << amplitude << " " << mean << " " << sigma << " " << background << endl;
  }
  outTxtFile.close();


  /* Drawing the plots so you can see them */
  if (draw) {
    for (int i=0;i<productHists.size();i++) {
      name.str(""); name << "localUnnormalizedLikelihood_" << i;
      TCanvas *c1 = new TCanvas(name.str().c_str(),"",1000,500);
      c1->SetLogy();
      cout << "eventNumber:" << seedCandidates[i] << endl;
      productHists[i]->Draw();
      fitFuncs[i]->Draw("same");
      candHists[i]->Draw("same");
      name << ".png";
      c1->SaveAs(name.str().c_str());


      /*
      TCanvas *cMap = new TCanvas("antarcticaMap","",1000,1000);
      backCountMaps[i]->Draw("colz");
      backCountMaps[i]->GetXaxis()->SetRangeUser(candMaps[i]->GetX()[0]-8e5,candMaps[i]->GetX()[0]+8e5);
      backCountMaps[i]->GetYaxis()->SetRangeUser(candMaps[i]->GetY()[0]-4e5,candMaps[i]->GetY()[0]+4e5);
      candMaps[i]->Draw("same");
      name.str(""); name << "localUnnormalizedLikelihood_antMap_" << i << ".png";
      c1->SaveAs(name.str().c_str());
      */
    }
  }

  
  

  
  cout << "Okay bye that was fun :D" << endl;
  return;
}



void saveLocalDistributions(string date="10.26.17_14h28m") {
  /*

    Saves the 2d histograms for cut parameters of events near each candidate, with the candidate highlighted

   */
  stringstream name;

  gStyle->SetOptStat(0);

  char* dataDir = getenv("ANITA3_RESULTSDIR");
  name.str(""); name << dataDir << "/cluster/" << date << "/clusterBackground";
  TChain *summaryTree = loadWhatever(name.str(),64,false);

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  summaryTree->BuildIndex("eventNumber");

  //to get the event numbers
  TFile *candTreeFile = TFile::Open("/home/brotter/anita16/benCode/templateSearch/macros/weakIsolated_oct14.root");
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





double getGeomagJP(double measPol, double expPol) {
  /*
    Determines the "JP" value for the geomagnetic
    Plateaus at 1 for some sigma value (2.5 default) and then falls off as 1/sigma past that

    opt:
    measPol & expPol: 

   */

  //constants
  const double expSigma = 2.0; //degrees (should be variable...)
  const double measSigma = 5.0; //degrees (from WAIS calibration)
  const double sigmaPlateau = 2.5; //how far away before you stop just being 1

  double totSigma = TMath::Sqrt(pow(expSigma,2)+pow(measSigma,2));

  double polDiff = TMath::Abs(measPol-expPol);
  double outJP;
  if ((polDiff/totSigma) < 2.5) {
    outJP = 1.0;
  }
  else {
      outJP = sigmaPlateau / (polDiff/totSigma); //the plateau value keeps it from being discontinuous at that point
  }


  return outJP;
  
}



void JPStatisticBackground(string inFileName) {
  /*

    reads in a signal file and a background file, then calculates and plots the JP statistic

    inSigName: "signal" file name, which is just plotted in red on top of the background
    inBackName: "background" file name, which is fit and plotted in black
   */


  TChain *summaryTree;
  if (inFileName == "all") {
  TChain *summaryTree = loadGeoAssociated();
  }
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(inFileName.c_str());
  }

  int lenEntries = summaryTree->GetEntries();
  cout << "found " << lenEntries << " signal entries" << endl;  


  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("template",&tempSum);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  bool fSeedEvent;
  summaryTree->SetBranchAddress("fSeedEvent",&fSeedEvent);

  //the new stokes has to be copied too
  AnitaEventSummary::WaveformInfo *newStokes = NULL;
  summaryTree->SetBranchAddress("newStokes",&newStokes);

  //and the geomag stuff
  double geoMagExp;
  summaryTree->SetBranchAddress("geoMagExp",&geoMagExp);





  TH1D *hTempCoher_back = new TH1D("hTempCoher_back","Template Correlation (coher) background;correlation;count",1000,0,1);
  TH1D *hTempDecon_back = new TH1D("hTempDecon_back","Template Correlation (decon) background;correlation;count",1000,0,1);

  TH2D *h2Geomag_back = new TH2D("h2Geomag_back","Geomagnetic Expected vs Measured;expected plane (degree);measured plane (degrees)",91,-45,45,91,-45,45);
  TH1D *hGeomag_back = new TH1D("hGeomag_back","Geomagnetic JP value (1/sigma);Geomag JP value; count",100,0,1);

  TH1D *jp_back = new TH1D("jp_back","JP Statistic (Background)",1000,0,1.1);
  TH1D *jpLog_back = new TH1D("jpLog_back","JP Statistic (Background)",1000,-5,0);


  TH1D *hTempCoher_sig = new TH1D("hTempCoher_sig","Template Correlation (coher) candidate;correlation;count",1000,0,1);
  TH1D *hTempDecon_sig = new TH1D("hTempDecon_sig","Template Correlation (decon) candidate;correlation;count",1000,0,1);

  TH2D *h2Geomag_sig = new TH2D("h2Geomag_sig","Geomagnetic Expected vs Measured;expected plane (degree);measured plane (degrees)",91,-45,45,91,-45,45);
  TH1D *hGeomag_sig = new TH1D("hGeomag_sig","Geomagnetic JP value (1/sigma) Candidate;Geomag JP value; count",100,0,1);

  TH1D *jp_sig = new TH1D("jp_sig","JP Statistic (Candidate)",1000,0,1.1);
  TH1D *jpLog_sig = new TH1D("jpLog_sig","JP Statistic (Candidate)",1000,-5,0);


  cout << "Starting Loop: "<< endl;
  for (int entry=0; entry<lenEntries; entry++) {
    //    if (!(entry%10)) {cout << entry << "/" << lenEntries << endl; fflush(stdout); }
    summaryTree->GetEntry(entry);

    //    cout << evSum->eventNumber << ":" << endl;

    double geoMeas = newStokes->linearPolAngle();
    double geoExp = geoMagExp;
    double geomagJP = getGeomagJP(geoMeas,geoExp);

    double product = 1.0;
    product *= tempSum->coherent[0][0].cRay[4];
    product *= tempSum->deconvolved[0][0].cRay[4];
    product *= pow(geomagJP,2);

    //    cout << tempSum->coherent[0][0].cRay[4] << " " << tempSum->deconvolved[0][0].cRay[4] << " " << geomagJP << " " << product << endl;

    if (!fSeedEvent) {
      h2Geomag_back->Fill(geoMeas,geoExp);
      hTempCoher_back->Fill(tempSum->coherent[0][0].cRay[4]);
      hTempDecon_back->Fill(tempSum->deconvolved[0][0].cRay[4]);
      hGeomag_back->Fill(geomagJP);
      jp_back->Fill(product);
      jpLog_back->Fill(TMath::Log10(product));
    }
    else {
      h2Geomag_sig->Fill(geoMeas,geoExp);
      hTempCoher_sig->Fill(tempSum->coherent[0][0].cRay[4]);
      hTempDecon_sig->Fill(tempSum->deconvolved[0][0].cRay[4]);
      hGeomag_sig->Fill(geomagJP);
      jp_sig->Fill(product);
      jpLog_sig->Fill(TMath::Log10(product));
    }

  }
    
  string outName;
  if (inFileName == "all") {
    outName = "allCandiates_JP.root";
  }
  else {
    size_t pos = inFileName.find(".root");
    string baseName = inFileName.substr(0,pos);
    outName = baseName+"_JP.root";
  }
  TFile *outFile = TFile::Open(outName.c_str(),"recreate");
  hTempCoher_back->Write();
  hTempDecon_back->Write();
  h2Geomag_back->Write();
  hGeomag_back->Write();
  jp_back->Write();
  jpLog_back->Write();

  hTempCoher_sig->Write();
  hTempDecon_sig->Write();
  h2Geomag_sig->Write();
  hGeomag_sig->Write();
  jp_sig->Write();
  jpLog_sig->Write();

  outFile->Close();



  return;
}   







void multipleCandidateJP(string filename = "trueCandidates_oct14_reMasked.root") {

  
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(filename.c_str());
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " events to save geoAssociated events for" << endl;
  
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  
  string name;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    cout << "multipleCanidateJP(): Doing " << evSum->eventNumber << "..." << endl;
    string filename = "geoAssociated/geoAssociated_ev" + to_string(evSum->eventNumber) + "_geomag.root";
    JPStatisticBackground(filename);
  }
  
  return;
}




void saveJPPlots(string filename = "trueCandidates_oct14_reMasked.root") {

  
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(filename.c_str());
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " events to save geoAssociated events for" << endl;
  
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->SetLogy();
  string name;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    cout << "multipleCanidateJP(): Doing " << evSum->eventNumber << "..." << endl;
    string filename = "geoAssociated/geoAssociated_ev" + to_string(evSum->eventNumber) + "_geomag_JP.root";
    TFile *inFile = TFile::Open(filename.c_str());
    TH1D *hBack = (TH1D*)inFile->Get("jp_back");
    TH1D *hSig = (TH1D*)inFile->Get("jp_sig");
    c1->cd();
    string title = "JP Statistic ev" + to_string(evSum->eventNumber) + ";JP Value; Count";
    hBack->SetTitle(title.c_str());
    hBack->Draw();
    hSig->SetLineColor(kRed);
    hSig->Draw("same");
    title = "JP_ev" + to_string(evSum->eventNumber) + ".png";
    c1->SaveAs(title.c_str());
    
  }
  
  return;
}



void fitAndReturnLikelihood(string inFileName) {
  
  size_t pos = inFileName.find(".root");
  string baseName = inFileName.substr(0,pos);
  string outName = baseName+"_JPlikelihood.txt";
  ofstream outFile(outName.c_str(),ofstream::out);

  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->SetLogy();
  string name;

  TLine *cutLine = new TLine(-0.35305,1e-3,-0.35305,3e3);
  cutLine->SetLineStyle(2);
  cutLine->SetLineColor(kRed);


  
  TFile *inFile = TFile::Open(inFileName.c_str());
  TH1D *hBack = (TH1D*)inFile->Get("jp_back");
  hBack->RebinX(4);//down to 250
  hBack->GetXaxis()->SetRangeUser(-3,0);
  hBack->GetYaxis()->SetRangeUser(1e-3,3e3);
  hBack->SetStats(0);
  TH1D *hSig = (TH1D*)inFile->Get("jp_sig");
  c1->cd();
  string title = "JP Statistic;JP Value; Count";
  hBack->SetTitle(title.c_str());
  double hBackMode = hBack->GetBinCenter(hBack->GetMaximumBin());
  hBack->Draw();
  hBack->Fit("expo","","",0.15,1);
  TF1 *gausFit = hBack->GetFunction("expo");
  gausFit->SetLineColor(kGreen);
  double constant = gausFit->GetParameter(0);
  double mean = gausFit->GetParameter(1);
  
  double candVal = hSig->GetBinCenter(hSig->GetMaximumBin());
  


  //  TF1 *unLog = new TF1("unLog","([0]/(x*[2]))*exp(-1*pow(log(x-[1]),2)/(2*pow([2],2)))",0.2,1);
  //  TF1 f1 = new TF1("logNormal","[2]*ROOT::Math::lognormal_pdf(x,[0],[1])",0.04,1);


  //  TCanvas *c2 = new TCanvas("c2","c2",1000,500);
  //  TH1D *hBackUnLog = (TH1D*)inFile->Get("jp_back");
  //  hBackUnLog->Draw("");
  //  unLog->Draw("same");


  //  ROOT::Math::WrappedTF1 wf1(*gausFit);
  //  ROOT::Math::GaussIntegrator ig;
  //  ig.SetFunction(wf1);
  //  ig.SetRelTolerance(0.001);
  
  //  double cutLineValue = TMath::Log10(0.666*0.666);
  
  //  cout << ig.Integral(cutLineValue,10) << " " << ig.Integral(candVal,10) << endl;
  
  //  double pCandVal = -constant*TMath::Sqrt(TMath::PiOver2()) * sigma * (1 + TMath::Erf((mean-candVal)/(TMath::Sqrt(2)*sigma)));
  //  double pBackground = -constant*TMath::Sqrt(TMath::PiOver2()) * sigma * (1 + TMath::Erf((mean-cutLineValue)/(TMath::Sqrt(2)*sigma)));
  
  //  outFile << constant << " " << mean << " " << sigma << " " << cutLineValue << " " << pBackground << " " << candVal << " " << pCandVal << endl;
  //  cout << constant << " " << mean << " " << sigma << " " << cutLineValue << " " << pBackground << " " << candVal << " " << pCandVal << endl;
  
  c1->cd();
  gausFit->SetRange(0.1,1);
  cutLine->Draw("same");
  hSig->SetLineColor(kRed);
  hSig->Draw("same");
  title = "JPLog.png";
  c1->SaveAs(title.c_str());
  
  
  outFile.close();

  return;
}


void fitAndReturnLikelihoodMultiple(string filename ="trueCandidates_oct14_reMasked.root") {
  
  return;
}


void mergeAllCandidates(string filename ="trueCandidates_oct14_reMasked.root") {

  
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(filename.c_str());
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " events to save geoAssociated events for" << endl;
  
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  
  ofstream outFile("JPlikelihood.txt",ofstream::out);

  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->SetLogy();
  string name;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    if (evSum->eventNumber == 39599205) continue;
    cout << "multipleCanidateJP(): Doing " << evSum->eventNumber << "..." << endl;
    string filename = "geoAssociated/geoAssociated_ev" + to_string(evSum->eventNumber) + "_geomag_JP.root";
    TFile *inFile = TFile::Open(filename.c_str());
    TH1D *hBack = (TH1D*)inFile->Get("jp_back");
    hBack->SetStats(0);
    TH1D *hSig = (TH1D*)inFile->Get("jp_sig");
    c1->cd();
    string title = "JP Statistic ev" + to_string(evSum->eventNumber) + ";Log10(JP Value); Count";
    hBack->SetTitle(title.c_str());
    hBack->GetXaxis()->SetRangeUser(0.1,0.40);
    hBack->Draw();
    hBack->Fit("expo");
    hBack->GetXaxis()->SetRangeUser(0.0,1.0);
    hBack->Draw();

    hBack->GetFunction("expo")->SetLineColor(kGreen);
    double constant = hBack->GetFunction("expo")->GetParameter(0);
    double mean = hBack->GetFunction("expo")->GetParameter(1);

    outFile << evSum->eventNumber << " " << constant << " " << mean << endl;

    hSig->SetLineColor(kRed);
    hSig->Draw("same");
    title = "JPLog_ev" + to_string(evSum->eventNumber) + ".png";
    c1->SaveAs(title.c_str());
    
  }
  
  outFile.close();

  return;
}

















