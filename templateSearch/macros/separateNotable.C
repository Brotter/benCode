/*

Separates out the notable events (wais, ldb, things that have good simple cut values) and saves them to a different file

*/



void separateNotable() {

  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  int lenEntries = summaryTree->GetEntries();
  cout << lenEntries << " total entries found" << endl;

  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  AnitaNoiseSummary *noiseSummary = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);


  TFile *outFile = TFile::Open("notableEvents.root","recreate");
  
  TTree *waisTree = new TTree("waisSummary","waisSummary");
  waisTree->Branch("eventSummary",&eventSummary);
  waisTree->Branch("template",&templateSummary);
  waisTree->Branch("noiseSummary",&noiseSummary);
  

  TTree *ldbTree = new TTree("ldbSummary","ldbSummary");
  ldbTree->Branch("eventSummary",&eventSummary);
  ldbTree->Branch("template",&templateSummary);
  ldbTree->Branch("noiseSummary",&noiseSummary);
    

  TTree *cutTree = new TTree("cutSummary","cutSummary");
  cutTree->Branch("eventSummary",&eventSummary);
  cutTree->Branch("template",&templateSummary);
  cutTree->Branch("noiseSummary",&noiseSummary);

  TTree *minbiasTree = new TTree("decMinBiasSummary","decMinBiasSummary");
  minbiasTree->Branch("eventSummary",&eventSummary);
  minbiasTree->Branch("template",&templateSummary);
  minbiasTree->Branch("noiseSummary",&noiseSummary);
  
  

  int waisCnt = 0;
  int ldbCnt = 0;
  int cutCnt = 0;
  int minbiasCnt = 0;
  
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%1000 == 0) {
      cout << entry << "/" << lenEntries;
      cout << "(" << waisCnt << " | " << ldbCnt << " | " << cutCnt << ")" << endl;
    }

    summaryTree->GetEntry(entry);

    if (eventSummary->flags.pulser == 1) {
      outFile->cd();
      waisTree->Fill();
      waisCnt++;
    }
    if (eventSummary->flags.pulser == 2) {
      outFile->cd();
      ldbTree->Fill();
      ldbCnt++;
    }

    if ( (eventSummary->flags.pulser == 0) &&
	 (eventSummary->peak[0][0].value > 0.04) &&
	 (templateSummary->coherent[0][0].cRay[5] > 0.7) &&
	 ((FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->wais.phi)) > 2.5) ||
	  (eventSummary->wais.distance > 1000e3) ) &&
	 ((FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->ldb.phi)) > 2.5) ||
	  (eventSummary->ldb.distance > 1000e3) ) ) {
      outFile->cd();
      cutTree->Fill();
      cutCnt++;
    }

    if ( !eventSummary->flags.isRF) {
      minbiasCut++;
      if (minbiasCut % 10 == 0) {
	outFile->cd();
	minbiasTree->Fill();
      }
    }

  }
  
  outFile->cd();
  waisTree->Write();
  ldbTree->Write();
  cutTree->Write();
  outFile->Close();


  cout << "Done!" << endl;

  cout << "Totals:" << endl << "Wais: " << waisCnt << endl << "LDB: " << ldbCnt << endl << "Cuts: " << cutCnt << endl;


  return;

}
