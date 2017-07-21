void getDateFromRealTime(int realTime,char* formattedDate,const int length) {
  
  time_t t = realTime;
  struct tm *tm = localtime(&t);
  strftime(formattedDate, sizeof(char)*length, "%m/%d %H:%M:%S", tm);
  
  return;
}

TProfile2D *TH2DtoTProfile2D(TH2D* inTH2D,double realTime) {

  int nBinX = inTH2D->GetNbinsX();
  int nBinY = inTH2D->GetNbinsY();
  int xMin  = inTH2D->GetXaxis()->GetBinLowEdge(1);
  int yMin  = inTH2D->GetYaxis()->GetBinLowEdge(1);
  int xMax  = inTH2D->GetXaxis()->GetBinUpEdge(nBinX);
  int yMax  = inTH2D->GetYaxis()->GetBinUpEdge(nBinY);
  TProfile2D *mapProfile = new TProfile2D("temp","temp",nBinX,xMin,xMax,nBinY,yMin,yMax);
  mapProfile->SetStats(0);
  stringstream name;
  name.str("");
  char currTime[64];
  getDateFromRealTime(realTime,currTime,64);
  cout << currTime << endl;
  name << "Average Interferometric Map - " << currTime;
  mapProfile->SetTitle(name.str().c_str());

  return mapProfile;

}

void fillTProfile2DWithTH2D(TProfile2D *prof, TH2D* hist) {

  int nBinX = hist->GetNbinsX();
  int nBinY = hist->GetNbinsY();

  for (int binX=0; binX<nBinX+1; binX++) {
    double x = hist->GetXaxis()->GetBinCenter(binX+1);
    for (int binY=0; binY<nBinY+1; binY++) {
      double y = hist->GetYaxis()->GetBinCenter(binY+1);
          prof->Fill(x,y,hist->GetBinContent(binX,binY));
    }
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

  TProfile2D *mapProfile = NULL;

  int lenEntries = summaryTree->GetEntries();
  cout << "lenEntries:" << lenEntries << endl;


  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    
    
    if (cnt%60 == 0) {
      if (mapProfile) {
	delete mapProfile;
      }
      mapProfile = TH2DtoTProfile2D(noiseSum->avgMapProf[0],eventSummary->realTime);
    }

    fillTProfile2DWithTH2D(mapProfile,noiseSum->avgMapProf[0]);

    if (cnt%60 == 59) {
      mapProfile->Draw("colz");
            
      name.str("");
      name << "avgMaps/" << prefix << "_";
      name << setfill('0') << setw(numZeros) << cnt/60 << ".png";
      if (saveImages) c1->SaveAs(name.str().c_str());
    }
    
    cnt++;

  }

  return;
}



TChain *loadWholeCluster(string date="07.19.17_16h/") {

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  stringstream name;

  char* resultsDir = getenv("ANITA3_RESULTSDIR");

  for (int run=130; run<440; run++) {
    name.str("");
    name << resultsDir << "noiseSummary/" << date << run << ".root";
    summaryTree->Add(name.str().c_str());
  }
  
  return summaryTree;
}


TChain *loadSingle(string name) {
  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(name.c_str());

  return summaryTree;
}


void drawAvgMaps(int run=-1) {

  if (run==-1) {
    cout << "No run selected, doing nothing" << endl;
    return;
  }

  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  name.str("");
  name << resultsDir << "07.19.17_16h/" << run << ".root";
  TChain *summaryTree = loadSingle(name.str());

  name.str("");
  name << "run" << setfill('0') << setw(3) << run;
  saveImagesFromTChain(summaryTree,name.str());

  return;
}
