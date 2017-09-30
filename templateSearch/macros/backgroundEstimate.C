/*

  I have some background estimate code that should live here

  I was living in cluster.C for whatever reason.

 */
#include "Normalize.C" //makeNormCumulative
#include "AnitaConventions.h"

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
  TFile *baseFile = TFile::Open("mergeBaseClusters.root");
  TTree *baseTree = (TTree*)baseFile->Get("summaryTree");

  //get candidate event summaries
  TFile *candidateFile = TFile::Open("candidates.root");
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
  TFile *baseFile = TFile::Open("mergeBaseClusters.root");
  TTree *baseTree = (TTree*)baseFile->Get("summaryTree");

  //get candidate event summaries
  TFile *candidateFile = TFile::Open("candidates.root");
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
  TH1D *hPeakHilbertDF = new TH1D("hPeakHilbertDF","Deconvolved Hilbert Peak;Deconvolved Hilbert Peak;count",250,0,500);
  baseTree->Draw("deconvolved_filtered[0][0].peakHilbert >> hPeakHilbertDF",allCuts);
  cout << "4" << endl;
  TH1D *hTemplate = new TH1D("hTemplate","cRay +4 Correlation;cRay +4 Correlation;count",250,0,1);
  baseTree->Draw("template.coherent[0][0].cRay[4] >> hTemplate",allCuts);

  TH1 *hMapPeakC   = makeNormCumulative(hMapPeak);
  TH1 *hLinPolFracC   = makeNormCumulative(hLinPolFrac);
  TH1 *hPeakHilbertDFC   = makeNormCumulative(hPeakHilbertDF);
  TH1 *hTemplateC   = makeNormCumulative(hTemplate);

  if (fileName=="") fileName="pValues.csv";
  ofstream outFile(fileName);
  outFile << "#eventNumber,mapPeakVal,PmapPeakVal,mapPeakSNR,PmapPeakSNR,peakHilbertDF,PpeakHilbertDF,cRayTemplate,PcRayTemplate" << endl;

  for (int entry=0; entry<candidateTree->GetEntries(); entry++) {
    candidateTree->GetEntry(entry);
    
    bool skip=false;

    if (TMath::Abs(candSum->peak[0][0].hwAngle) > 45) {
      cout << candSum->eventNumber << " hw angle doesn't match peak angle" << endl;
      continue;
    }
    if (candSum->flags.maxBottomToTopRatio[0] > 3 ) {
      cout << candSum->eventNumber << " blast" << endl;
      continue;
    }
   
    if (candSum->eventNumber == 84114142 || candSum->eventNumber == 84405480) {
      cout << candSum->eventNumber << " doesn't have right geomagnetic" << endl;
      continue;
    }


    int bin;
    double mapPeakVal = candSum->peak[0][0].value;
    bin = hMapPeakC->GetXaxis()->FindBin(mapPeakVal);
    double PmapPeakVal = hMapPeakC->GetBinContent(bin);

    double linearPolFrac = candSum->coherent_filtered[0][0].linearPolFrac();
    bin = hLinPolFracC->GetXaxis()->FindBin(linearPolFrac);
    double PlinearPolFrac = hLinPolFracC->GetBinContent(bin);

    double peakHilbert = candSum->deconvolved_filtered[0][0].peakHilbert;
    bin = hPeakHilbertDFC->GetXaxis()->FindBin(peakHilbert);
    double PpeakHilbert = hPeakHilbertDFC->GetBinContent(bin);

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


int numberOfEventsExceedingCandidate(int candNum) {
  /*

    Reads in the pseudoBaseCluster.root list, which includes all events from near events that pass cuts, and 
    returns the number of events that have reduced quantities higher than the candidate selected

    The candNum is just the index it occurs in, 
    if you pick an index higher than the total number of candidates it will yell at you

   */

  TFile *candFile = TFile::Open("candidates.root");
  TTree *candTree = (TTree*)candFile->Get("summaryTree");
  AnitaEventSummary *candEvSum = NULL;
  candTree->SetBranchAddress("eventSummary",&candEvSum);
  AnitaTemplateSummary *candTempSum = NULL;
  candTree->SetBranchAddress("template",&candTempSum);

  if (candNum >= candTree->GetEntries()) {
    cout << "selected candidate number " << candNum << " is greater than number of candidates (" << candTree->GetEntries() << ")";
    cout << endl << "Quitting" << endl;
    return -1;
  }

  candTree->GetEntry(candNum);


  //Cut values:
  double templateCorr = candTempSum->coherent[0][0].cRay[4];
  double mapPeak = candEvSum->peak[0][0].value;
  double hilbPeak = candEvSum->deconvolved_filtered[0][0].peakHilbert;
  double linPolFrac = candEvSum->coherent[0][0].linearPolFrac();
  

  
  TFile *backFile = TFile::Open("pseudoBaseCluster.root");
  TTree *backTree = (TTree*)backFile->Get("summaryTree");
  AnitaEventSummary *backEvSum = NULL;
  backTree->SetBranchAddress("eventSummary",&backEvSum);
  AnitaTemplateSummary *backTempSum = NULL;
  backTree->SetBranchAddress("template",&backTempSum);
  
  
  int lenEntries = backTree->GetEntries();
  cout << "Found " << lenEntries << " background events from impulsive sources" << endl;
  
  stringstream name;
  name.str("");
  name << "template.coherent[0][0].cRay[4] > " << templateCorr;
  name << " && ";
  name << "peak[0][0].value > " << mapPeak;
  name << " && ";
  name << "deconvolved_filtered[0][0].peakHilbert > " << hilbPeak;
  name << " && ";
  name << "coherent_filtered[0][0].linearPolFrac() > " << linPolFrac;
  
  int numPassing = backTree->Draw("eventNumber",name.str().c_str(),"goff");

  return numPassing;
}
  
