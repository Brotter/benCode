#include "AnitaEventSummary.h"

void globInteresting(string date = "06.07.17_17h"){

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  stringstream name;
  for (int run=130; run<433; run++) {
    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    name << resultsDir << date << "/" << run << ".root";

    summaryTree->Add(name.str().c_str());

  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  
  TFile *outFile = TFile::Open("interesting.root","recreate");
  TTree *outTree = new TTree("interestingTree","interestingTree");

  double templateCRayH[10];
  double templateCRayV[10];
  summaryTree->SetBranchAddress("templateCRayH",&templateCRayH);
  summaryTree->SetBranchAddress("templateCRayV",&templateCRayV);
  outTree->Branch("templateCRayH",&templateCRayH,"templateCRayH[10]/D");
  outTree->Branch("templateCRayV",&templateCRayV,"templateCRayV[10]/D");

  AnitaEventSummary::SourceHypothesis wais,ldb,sun;
  summaryTree->SetBranchAddress("sun",&sun);
  summaryTree->SetBranchAddress("wais",&wais);
  summaryTree->SetBranchAddress("ldb",&ldb);
  outTree->Branch("sun",&sun);
  outTree->Branch("wais",&wais);
  outTree->Branch("ldb",&ldb);

  AnitaEventSummary::PointingHypothesis peakH,peakV;
  outTree->Branch("peakH",&peakH);
  outTree->Branch("peakV",&peakV);

  AnitaEventSummary::WaveformInfo coherentH,coherentV;
  outTree->Branch("coherentH",&coherentH);
  outTree->Branch("coherentV",&coherentV);

  int eventNumber,realTime;
  outTree->Branch("realTime",&realTime);
  outTree->Branch("eventNumber",&eventNumber);


  int lenEntries = summaryTree->GetEntries();
  
  cout << "Found " << lenEntries << " entries" << endl;

  for (int entry=0; entry<lenEntries; entry++) {
    if (entry % 10000 == 0) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);

    realTime = summary->realTime;
    eventNumber = summary->eventNumber;
    
    peakH = summary->peak[0][0];
    peakV = summary->peak[1][0];

    coherentH = summary->coherent[0][0];
    coherentV = summary->coherent[1][0];

    outFile->cd();
    outTree->Fill();
  }

  outTree->Write();
  outFile->Close();


  return;

}


