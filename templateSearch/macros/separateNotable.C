/*

Separates out the notable events (wais, ldb, things that have good simple cut values) and saves them to a different file

*/

#include "loadAll.C"



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


void separateTrueCandidates() {
  /*
    Okay I've been dealing with these dumb lists from forever ago, so lets make a new one with direct and reflected and only the non-blast things.
   */

  

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add("candidates.root");
  summaryTree->Add("aboveHorizon.root");

  TFile *outFile = TFile::Open("trueCandidates.root","recreate");
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

  int lenEntries = summaryTree->GetEntries();
  cout << "number of events in tree: " << lenEntries << endl;

  int count=0;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    cout << entry << "/" << lenEntries << " : " << evSum->eventNumber << endl;

    if (evSum->flags.maxBottomToTopRatio[0] > 6) {
      cout << "nope, blast event" << endl;
      continue;
    }
    if (TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0)) < 6 && evSum->ldb.distance < 650e3) {
      cout << "nope, points at ldb" << endl;
      continue;
    }
    if (TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)) < 6 && evSum->wais.distance < 650e3) {
      cout << "nope, points at wais" << endl;
      continue;
    }
    
    if (evSum->eventNumber >= 80299371 && evSum->eventNumber <= 80700153) {
      cout << "these are the above horizon clustered ones" << endl;
      continue;
    }

    if (evSum->eventNumber == 9855625) {
      cout << "points at LDB" << endl;
    }

    if (evSum->eventNumber == 69050312 || evSum->eventNumber == 78946778) {
      cout << "too close to the horizon" << endl;
      continue;
    }

    count++;

    outFile->cd();
    cutTree->Fill();

    cout << "good" << endl;
  }

  cout << "total saved: " << count << endl;

  outFile->cd();
  cutTree->Write();
  outFile->Close();


  return;

}



void saveOnlyGoodEvents(string date="09.27.17_19h") {
  /*

    lots of events seem to point either above zero degrees, or are like blasts or minbias

    lets save all the ones that aren't

   */

  TChain *summaryTree = loadAll(date,false);
  int lenEntries = summaryTree->GetEntries();

  char* dataDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  name.str("");
  name << dataDir << "templateSearch/" << date << "/goodEvents.root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
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
  
  
  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  int savedCount = 0;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0 && entry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entry)/totalTimeSec;
      double remaining = (float(lenEntries-entry)/rate)/60.;
      watch.Start();
      cout << entry << "/" << lenEntries << " (" << savedCount << ") ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;
    }
    summaryTree->GetEntry(entry);

    if (evSum->peak[0][0].theta < 0 || !evSum->flags.isRF || evSum->flags.maxBottomToTopRatio[0] > 3) {
      continue;
    }
    savedCount++;
    outFile->cd();
    cutTree->Fill();
  }
  cout << "done looping, saving..." << endl;


  outFile->cd();
  cutTree->Write();
  outFile->Close();
    

  cout << "Done!" << endl;

  return;

}

void saveNonClusteredEvents(string clusterFileName="cluster.root", string outFileName="candidates.root",
			    double threshold = 10.,string rootFileName=""){
  /*

    This takes the output of cluster.root, specifically the gClosest TGraph, and saves events with
    clustering values greater than the input threshold

    I was doing 40 before, but 10 is probably more correct

  */


  TFile *inFile = TFile::Open(clusterFileName.c_str());
  TGraph *gClosest = (TGraph*)inFile->Get("gClosest");
  int lenImp = gClosest->GetN();
  cout << "Found " << lenImp << " events in the cut list" << endl;

  TChain *summaryTree;
  if (rootFileName=="") summaryTree = loadAllDefault_noproof();
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(rootFileName.c_str());
  }
  cout << "building index... this might take awhile... "; fflush(stdout);
  summaryTree->BuildIndex("eventNumber");
  cout << "okay done!" << endl;

  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
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

  int savedEvs = 0;
  for (int imp=0; imp<lenImp; imp++) {
    if (gClosest->GetY()[imp] < threshold) continue;
    int eventNumber = gClosest->GetX()[imp];
    int entry = summaryTree->GetEntryNumberWithIndex(eventNumber);
    if (entry==-1) {
      cout << "I couldn't find event number " << eventNumber << " in the dataset oh no!!!" << endl;
      continue;
    }
    summaryTree->GetEntry(entry);
    outFile->cd();
    cutTree->Fill();
    savedEvs++;
  }
  
  outFile->cd();
  cutTree->Write();
  outFile->Close();
      
  cout << "done :) saved " << savedEvs << "events" << endl;
}



void savePassingEvents(string outFileName, int strength=4, bool save=true) {
  /*

    the makeCuts -> separateNotable way of doing this is stupid

    This might be better.  Its like separateNotable_fromScratch but newer

    I've decided on doing two distinct cut strengths. Right now 4 (historically numbered) is weak cuts

   */

  cout << "Starting savePassingEvents: using strength " << strength;
  cout << " and output file " << outFileName << endl;

  //  TChain *summaryTree = loadAll("09.27.17_19h",false);
  TFile *inFile = TFile::Open("/home/brotter/nfsShared/results/templateSearch/09.27.17_19h/goodEvents.root");
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  int lenEntries = summaryTree->GetEntries();

  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
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

  
  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  int savedCount = 0;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0 && entry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entry)/totalTimeSec;
      double remaining = (float(lenEntries-entry)/rate)/60.;
      watch.Start();
      cout << entry << "/" << lenEntries << " (" << savedCount << ") ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;
    }
    summaryTree->GetEntry(entry);
    
    /* \/\/\/\/\/   CUTS!    \/\/\/\/\/\ */

    /* Varying strength cuts go here */
    if (strength == 4 ) {
      if (tempSum->coherent[0][0].cRay[4] < 0.5) continue;
      //      if (tempSum->coherent[0][0].cRay[4] > 0.67) continue;
      if (evSum->coherent_filtered[0][0].peakHilbert < 25) continue;
      if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.6) continue;
      if (evSum->peak[0][0].value < 0.0435) continue;
      if (evSum->peak[0][0].snr < 8.95) continue;
    }
    /* --------------- */
    
    /* Cuts that should _always_ be made for candidates*/
    // not flagged as a pulser
    if (evSum->flags.pulser != 0) continue;
    // needs an rf trigger
    if (!evSum->flags.isRF) continue;
    //not pointed at ldb when it is nearby (~700km)
    if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0)),2) + pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->ldb.theta,360,0)),2))  < 6 
	 && evSum->ldb.distance  < 700e3)) continue;
    //not pointed at wais when it is nearby (~700km)
    if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)),2) + pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->wais.theta,360,0)),2))  < 6 
	 && evSum->wais.distance  < 700e3)) continue;
    //not a blast event
    if (evSum->flags.maxBottomToTopRatio[0] > 3) continue;
    //not pointing "above" zero (- == up), since not even direct CRs will be above zero (no atmosphere!)
    if (evSum->peak[0][0].theta < 0) continue;

    /* ------------ */
    /* \/\/\/\/\/\/\/\/\/\/\ */

    outFile->cd();
    cutTree->Fill();
    savedCount++;
  }
  cout << "done scanning, found " << savedCount << " events" << endl;

  outFile->cd();
  cutTree->Write();
  outFile->Close();

  cout << "Bye!" << endl;

  return;

}



void makeHarshCuts(string inFileName, string outFileName) {
  /*

    Takes a summaryTree and makes two really hard cuts on the events

   */

  //  TChain *summaryTree = loadAll("09.27.17_19h",false);
  TFile *inFile = TFile::Open(inFileName.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");
  int lenEntries = summaryTree->GetEntries();

  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
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

  
  for (int i=0; i<lenEntries; i++) {
    summaryTree->GetEntry(i);
    if (tempSum->coherent[0][0].cRay[4] < 0.78) continue;
    if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.75) continue; //0.8 is too harsh, 0.75 looks good though (only lose one)

    outFile->cd();
    cutTree->Fill();

  }

  outFile->cd();
  cutTree->Write();
  outFile->Close();


  cout << "Done!" << endl;
  return;
}





void mergeTwoSummaries(string filename1, string filename2, string filenameOut) {
  /*

    In case you want to turn two files with summaryTrees into one

   */
  TFile *outFile = TFile::Open(filenameOut.c_str(),"recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(filename1.c_str());
  summaryTree->Add(filename2.c_str());


  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);

  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " to put into output file " << filenameOut << endl;

  for (int i=0; i<lenEntries; i++) {
    if (i%10000) cout << i << "/" << lenEntries << endl;
    summaryTree->GetEntry(i);
    outFile->cd();
    outTree->Fill();
  }

  outFile->cd();
  outTree->Write();
  outFile->Close();

  cout << "done!" << endl;

  return;

}




void combineWAISTrees() {
  /*

    takes in the waisEvents.root file with all the wais events, and shoves the snr and efficiency branches
    into them too.  This is only required for the snr branch really, since it isn't the same dimensions

   */

  TFile *waisFile = TFile::Open("waisEvents.root");
  TTree *summaryTree = (TTree*)waisFile->Get("summaryTree");



  TFile *outFile = TFile::Open("waisEvents_comb.root","recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);

  //the snr tree syncs up with the full data set
  TFile *snrFile = TFile::Open("/home/brotter/nfsShared/results/templateSearch/09.27.17_19h/SNRs.root");
  TTree *snrTree = (TTree*)snrFile->Get("snrTree");
  cout << "snrTree has " << snrTree->GetEntries() << " entries" << endl;
  cout << "Building index..."; fflush(stdout);
  snrTree->BuildIndex("eventNumber");
  cout << "done!" << endl;
  double snr,snr_filtered;
  snrTree->SetBranchAddress("snr",&snr);
  outTree->Branch("snr",&snr);
  snrTree->SetBranchAddress("snr_filtered",&snr_filtered);
  outTree->Branch("snr_filtered",&snr_filtered);

  //the efficiency tree is per wais event which is sort of annoying
  TFile *effFile = TFile::Open("/home/brotter/anita16/benCode/templateSearch/macros/waisEfficiency.root");
  TTree *effTree = (TTree*)effFile->Get("waisEfficiency");
  summaryTree->AddFriend(effTree);
  int waisCount; //waisCount is over 1000 seconds
  effTree->SetBranchAddress("count",&waisCount); 
  outTree->Branch("waisCount",&waisCount);


  int notFound=0;
  int found = 0;
  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    if (!entry%10000) cout << entry << "/" << effTree->GetEntries() << endl;

    summaryTree->GetEntry(entry);
    
    int eventNumber = evSum->eventNumber;
    int snrEntry = snrTree->GetEntryNumberWithIndex(eventNumber);
    if (snrEntry == -1) {
      cout << "couldn't find evNum:" << eventNumber << " (" << notFound << "/" << found << ")" << endl;
      notFound++;
    }
    else found++;
    snrTree->GetEntry(entry);


    //okay they're all together, save em`
    outFile->cd();
    outTree->Fill();
  }

  cout << "Saving..." << endl;

  outFile->cd();
  outTree->Write();
  outFile->Close();

  cout << "Done!" << endl;

}


void separateWAIS() {
  /*

    splits out the WAIS events

   */


  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  char* dataDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  for (int i=160; i<185; i++) { //wais runs go from ~160 to ~185
    name.str("");
    name << dataDir << "templateSearch/09.27.17_19h/" << i << ".root";
    summaryTree->Add(name.str().c_str());
  }
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " in full set" << endl;

  TFile *outFile = TFile::Open("waisEvents.root","recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);

  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  int savedCount = 0;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0 && entry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entry)/totalTimeSec;
      double remaining = (float(lenEntries-entry)/rate)/60.;
      watch.Start();
      cout << entry << "/" << lenEntries << " (" << savedCount << ") ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;
    }
    summaryTree->GetEntry(entry);

    int eventNumber = evSum->eventNumber;

    //if it isn't a wais pulser then skip it
    if (evSum->flags.pulser != 1) continue;

    //like 100 events don't point right and should be excluded from now on
    if (TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)) > 6) continue;
  
    //otherwise you're good!  save it.
    outFile->cd();
    outTree->Fill();
    savedCount++;
  }

  cout << "Found " << savedCount << " entries, Saving..." << endl;

  outFile->cd();
  outTree->Write();
  outFile->Close();

  cout << "Done!" << endl;

}


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




void separateNotable_fromCSV(string fileName="makeCuts_weak.csv",string outFileName="") {
  /*
    Take the output of a Scan file (that has had the *'s replaced with nothing and the first few lines cut off)
    command:
    sed -e 's|*||g' makeCuts_weak.txt | sed 1,4d >> makeCuts_weak.csv

    And saves those files to a new root file
   */

  TGraph *evNums = new TGraph(fileName.c_str());
  cout << "Opened " << fileName << " as input file, found " << evNums->GetN() << " events to save" << endl;

  TChain *summaryTree = loadAll("09.27.17_19h");
  summaryTree->BuildIndex("eventNumber");

  if (outFileName=="") outFileName = "cuts.root";
  cout << "Using: " << outFileName << " as output file" << endl;
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
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


void separateNotable_hardCodedEvNums() {


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



void separateNotable() {
  cout << "loaded separateNotable.C" << endl;
  return;
}
