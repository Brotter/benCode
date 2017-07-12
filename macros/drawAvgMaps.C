
void getDateFromRealTime(int realTime,char* formattedDate,const int length) {
  
  time_t t = realTime;
  struct tm *tm = localtime(&t);
  strftime(formattedDate, sizeof(char)*length, "%m/%d %H:%M:%S", tm);
  
  return;
}


void saveImagesFromTChain(TChain *summaryTree,string prefix="") {
  
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
  

  int lenEntries = summaryTree->GetEntries();
  cout << "lenEntries:" << lenEntries << endl;

  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    
    
    if (noiseSum->isMinBias) {
      cnt++;
      if (cnt%60 != 0) continue; //only do every 60

      cout << "entry:" << entry << "/" << lenEntries << " (" << cnt << ")" <<  endl;

      c1->cd();
      c1->Clear();
      
      

      noiseSum->avgMapProf[0]->SetStats(0);
      noiseSum->avgMapProf[0]->GetZaxis()->SetRangeUser(-0.02,0.05);
      name.str("");
      char currTime[64];
      getDateFromRealTime(eventSummary->realTime,currTime,64);
      cout << currTime << endl;
      name << "Average Interferometric Map - " << currTime;
      noiseSum->avgMapProf[0]->SetTitle(name.str().c_str());

      noiseSum->avgMapProf[0]->Draw("colz");

      double sunTheta = eventSummary->sun.theta;
      double sunPhi = eventSummary->sun.phi;
      if (sunPhi < 0) sunPhi += 360;
      //      cout << sunPhi << " " << sunTheta << endl;
      gSun->SetPoint(0,sunPhi,-sunTheta);
      gSun->Draw("pSame");
      

      name.str("");
      name << "avgMaps/" << prefix << "_";
      name << setfill('0') << setw(numZeros) << cnt/60 << ".png";
      c1->SaveAs(name.str().c_str());

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
  name << "core" << core;
  saveImagesFromTChain(summaryTree,name.str());

  return;
}
