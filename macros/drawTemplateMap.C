#include "AnitaEventSummary.h"

void drawTemplateMap() {

  stringstream name;

  //  string basedir = "/Volumes/ANITA3Data/bigAnalysisFiles/templateSearch/06.11.17_19h/";
  string basedir = "~/nfsShared/results/templateSearch/06.11.17_19h/";


  TChain *summaryTree = new TChain("summaryTree");
  for (int run=130; run<440; run++) {
    name.str("");
    name << basedir << run << ".root";
    summaryTree->Add(name.str().c_str());
  }

  gROOT->ProcessLine(".x setupProof.C");
  summaryTree->SetProof();

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  
  double crTemplate[10];
  double waisTemplate,irTemplate;
  summaryTree->SetBranchAddress("templateImpH",&irTemplate);
  summaryTree->SetBranchAddress("templateWaisH",&waisTemplate);
  summaryTree->SetBranchAddress("templateCRayH",&crTemplate);

  int lenEntries = summaryTree->GetEntries();

  cout << "Found " << lenEntries << " entries" << endl;
  
  TH2D *hWaterfall = new TH2D("hWaterfall","template correlation value waterfall",13,-0.5,12.5, lenEntries,0,lenEntries);

  TH1D *hMaxes = new TH1D("hMaxes","maximum correlation values",100,0,1);
  TH1D *hMaxLocs = new TH1D("hMaxLocs","maximum correlation value templates",13,-0.5,12.5);

  TH2D *hMaxVsLoc = new TH2D("hMaxVsLoc","maximum value vs location",13,-0.5,12.5, 100,0,1);

  TH2D *h2Maxes = new TH2D("h2Maxes","maximum value per location",13,-0.5,12.5,100,0,1);

  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << endl; fflush(stdout);
    summaryTree->GetEntry(entry);

    if (summary->flags.pulser != 0) continue;

    for (int i=0; i<10; i++) {
      hWaterfall->Fill(i,entry,crTemplate[i]);
      h2Maxes->Fill(i,crTemplate[i]);
    }
    hWaterfall->Fill(11,entry,irTemplate);
    hWaterfall->Fill(12,entry,waisTemplate);

    h2Maxes->Fill(11,irTemplate);
    h2Maxes->Fill(12,waisTemplate);


    double allTempValues[12];
    for (int i=0; i<10; i++) {
      allTempValues[i] = crTemplate[i];
    }
    allTempValues[10] = irTemplate;
    allTempValues[11] = waisTemplate;

    double max = TMath::MaxElement(12,allTempValues);
    int maxLoc = TMath::LocMax(12,allTempValues);
    
    if (maxLoc > 9) maxLoc++;

    hMaxLocs->Fill(maxLoc);

    hMaxVsLoc->Fill(maxLoc,max);

    


  }

  h2Maxes->Draw("colz");


  TFile *outFile = TFile::Open("drawTemplateMap.root","recreate");
  hWaterfall->Write();
  hMaxes->Write();
  hMaxLocs->Write();
  hMaxVsLoc->Write();
  h2Maxes->Write();
  outFile->Close();

}
      
