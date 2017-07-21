void getDateFromRealTime(int realTime,char* formattedDate,const int length) {
  
  time_t t = realTime;
  struct tm *tm = localtime(&t);
  strftime(formattedDate, sizeof(char)*length, "%m/%d %H:%M:%S", tm);
  
  return;
}

void showImagesFromTChain(TChain *summaryTree) {

  AnitaNoiseSummary *noiseSum = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);

  int lenEntries = summaryTree->GetEntry();

  TCanvas *c1 = new TCanvas("c1","c1",1920,1080);
  c1->Draw();
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    c1->cd();
    if (noiseSum->avgMapProf[0] != NULL) { 
      noiseSum->avgMapProf[0]->Draw("colz");
      c1->Modified();
      c1->Update();
      cout << "drawn, waiting for enter" << endl;
      cin.get(); 
    }
    else { cout << "empty" << endl; }
  }

  return;
}


void saveImagesFromTChain(TChain *summaryTree,string prefix="",bool saveImages=true) {
  
  stringstream name;
  
  AnitaNoiseSummary *noiseSum = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  
  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TGraph *gSun = new TGraph();
  gSun->SetMarkerStyle(kStar);


  //counter since I only want to save every 60
  int cnt = 0;
 
  const int numZeros = 5;
  
  TChain *gpsTree = new TChain("adu5PatTree","adu5PatTree");
  for (int run=130; run<440; run++) {
    name.str("");
    char* anitaData = getenv("ANITA3_DATA");
    name << anitaData << "/run" << run << "/gpsEvent" << run << ".root";

    gpsTree->Add(name.str().c_str());
  }

  cout << "Building gps index" << endl;
  gpsTree->BuildIndex("eventNumber");

  Adu5Pat *pat = NULL;
  gpsTree->SetBranchAddress("pat",&pat);

  TProfile2D *mapProfile = NULL;

  int lenEntries = summaryTree->GetEntries();
  cout << "lenEntries:" << lenEntries << endl;



  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    
    if (cnt%60 == 0) {
      mapProfile->Draw("colz");
      
      double sunTheta = eventSummary->sun.theta;
      double sunPhi = eventSummary->sun.phi - heading;
      while (sunPhi < 0) sunPhi += 360;
      gSun->SetPoint(0,sunPhi,-sunTheta);
      gSun->Draw("pSame");
      
      
      name.str("");
      name << "avgMaps/" << prefix << "_";
      name << setfill('0') << setw(numZeros) << cnt/60 << ".png";
      if (saveImages) c1->SaveAs(name.str().c_str());
    }
    
    cnt++;

  }
}

  return;
}



TChain *loadWholeCluster(string date="07.05.17_22h/") {

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  stringstream name;

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  for (int core=0; core<256; core++) {
    name.str("");
    name << resultsDir << date << core << ".root";
    summaryTree->Add(name.str().c_str());
  }
  
  return summaryTree;
}


TChain *loadSingle(string name) {
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(name.c_str());

  return summaryTree;
}


void drawAvgMaps(int core=-1) {

  if (core==-1) {
    cout << "No core selected, doing nothing" << endl;
    return;
  }

  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  name.str("");
  name << resultsDir << "07.05.17_22h/" << core << ".root";
  TChain *summaryTree = loadSingle(name.str());

  name.str("");
  name << "core" << setfill('0') << setw(3) << core;
  saveImagesFromTChain(summaryTree,name.str());

  return;
}
