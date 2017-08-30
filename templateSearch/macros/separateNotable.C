/*

Separates out the notable events (wais, ldb, things that have good simple cut values) and saves them to a different file

*/


void separateNotable_fromScratch() {
/*
  This one goes through the full dataset and grabs the ones that are important
*/

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
      minbiasCnt++;
      if (minbiasCnt % 10 == 0) {
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




void separateNotable_fromFile(string fileName="makeCuts_weak.csv") {
  /*
    Take the output of a Scan file (that has had the *'s replaced with nothing and the first few lines cut off)
    command:
    sed -e 's|*||g' makeCuts_weak.txt | sed 1,4d >> makeCuts_weak.csv

    And saves those files to a new root file
   */


  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");
  summaryTree->BuildIndex("eventNumber");

  TGraph *evNums = new TGraph(fileName.c_str());

  TFile *outFile = TFile::Open("cuts.root","recreate");
  TTree *cutTree = new TTree("summaryTree","summaryTree");
  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  cutTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  cutTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  cutTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  cutTree->Branch("gpsEvent",&gps);



  for (int ev=0; ev<evNums->GetN(); ev++) {
    cout << ev << "/" << evNums->GetN() << endl;
    int evNum = evNums->GetY()[ev];
    int entry = summaryTree->GetEntryNumberWithIndex(evNum);
    summaryTree->GetEntry(entry);
    outFile->cd();
    cutTree->Fill();
  }
  
  cutTree->Write();
  outFile->Close();
}


void separateNotable_hardCoded() {


  TFile *inFile = TFile::Open("cuts.root");  
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  summaryTree->BuildIndex("eventNumber");

  TFile *outFile = TFile::Open("unclusteredEvs.root","recreate");
  TTree *cutTree = new TTree("summaryTree","summaryTree");
  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  cutTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  cutTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  cutTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  cutTree->Branch("gpsEvent",&gps);


  vector<int> goodEvs;
  goodEvs.push_back(8135326);
  goodEvs.push_back(9097075);
  goodEvs.push_back(11116669);
  goodEvs.push_back(11989349);
  goodEvs.push_back(15717147);
  goodEvs.push_back(16952229);
  goodEvs.push_back(19459851);
  goodEvs.push_back(23695286);
  goodEvs.push_back(32907848); //BenS
  goodEvs.push_back(33484995); //BenS
  goodEvs.push_back(41529195); //BenS
  goodEvs.push_back(58592863); //BenS
  goodEvs.push_back(62273732);
  goodEvs.push_back(62365441);
  goodEvs.push_back(63210848);
  goodEvs.push_back(64201621);
  goodEvs.push_back(66313844);
  goodEvs.push_back(68298837);
  goodEvs.push_back(70013898);
  goodEvs.push_back(73726742);
  goodEvs.push_back(75277769);
  goodEvs.push_back(80561103);
  goodEvs.push_back(80973610);
  goodEvs.push_back(83877990);
  goodEvs.push_back(84114142);
  goodEvs.push_back(84405480);


  for (int ev=0; ev<goodEvs.size(); ev++) {
    int entry = summaryTree->GetEntryNumberWithIndex(goodEvs[ev]);
    summaryTree->GetEntry(entry);
    cutTree->Fill();
  }

  cutTree->Write();
  outFile->Close();
}




void separateNotable_singleEv(int evNum) {
			      
  /*
    Finds the selected event number and saves the AnitaEventSummary and Adu5Pat from it
   */


  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);


  summaryTree->BuildIndex("eventNumber");
  int entry = summaryTree->GetEntryNumberWithIndex(evNum);
  cout << "event number " << evNum << " is entry " << entry << endl;
  summaryTree->GetEntry(entry);
 

  stringstream name;
  name.str("");
  name << "ev" << evNum << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  cout << "saving as " << name.str() << endl;
  outFile->WriteObject(evSum,"eventSummary");
  outFile->WriteObject(gps,"gpsEvent");
    
  outFile->Close();

  return;
}
