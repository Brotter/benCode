TH2D* rotateMap(TH2D *inHist, double heading) {

  /* Returns a rotated map where zero is north */

  TH2D *outHist = (TH2D*)inHist->Clone();
  outHist->Reset();

  for (int binX=1; binX < inHist->GetNbinsX()+1; binX++) {
    double binCenterX = inHist->GetXaxis()->GetBinCenter(binX);
    binCenterX -= heading;
    if (binCenterX <= 0) binCenterX += 360;
    
    for (int binY=1; binY < inHist->GetNbinsY()+1; binY++) {
      double binCenterY = inHist->GetYaxis()->GetBinCenter(binY);
      double value = inHist->GetBinContent(binX,binY);
      
      outHist->Fill(binCenterX,binCenterY,value);
    }
  }

  return outHist;
}

  


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
  
  summaryTree->GetEntry(0);
  AnitaDataset *data = new AnitaDataset(eventSummary->run);


  int lenEntries = summaryTree->GetEntries();
  cout << "lenEntries:" << lenEntries << endl;

  TProfile2D *mapProfile = NULL;
  
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);

    if (noiseSum->isMinBias) {
      cout << "entry:" << entry << "/" << lenEntries << " (" << cnt << ")" <<  endl;

      c1->cd();
      c1->Clear();
    
      int run = eventSummary->run;
      if (data->currRun != run) {
	data->loadRun(run);
      }
      data->getEvent(eventSummary->eventNumber);

      double heading = data->gps()->heading;
      //TH2D *rotatedMap = rotateMap(noiseSum->avgMapProf[0],heading);
      TH2D *rotatedMap = (TH2D*)noiseSum->avgMapProf[0]->Clone();

      if (cnt%60 == 0) {
	if (mapProfile != NULL) {
	  delete mapProfile;
	}
	int nBinX = rotatedMap->GetNbinsX();
	int nBinY = rotatedMap->GetNbinsY();
	int xMin  = rotatedMap->GetXaxis()->GetBinLowEdge(1);
	int yMin  = rotatedMap->GetYaxis()->GetBinLowEdge(1);
	int xMax  = rotatedMap->GetXaxis()->GetBinUpEdge(nBinX);
	int yMax  = rotatedMap->GetYaxis()->GetBinUpEdge(nBinY);
	mapProfile = new TProfile2D("temp","temp",nBinX,xMin,xMax,nBinY,yMin,yMax);
	mapProfile->SetStats(0);
	name.str("");
	char currTime[64];
	getDateFromRealTime(eventSummary->realTime,currTime,64);
	cout << currTime << endl;
	name << "Average Interferometric Map - " << currTime;
	mapProfile->SetTitle(name.str().c_str());
      }

      /*
      rotatedMap->SetStats(0);
      rotatedMap->GetZaxis()->SetRangeUser(-0.02,0.04);
      name.str("");
      char currTime[64];
      getDateFromRealTime(eventSummary->realTime,currTime,64);
      cout << currTime << endl;
      name << "Average Interferometric Map - " << currTime;
      rotatedMap->SetTitle(name.str().c_str());
      */
      int nBinX = rotatedMap->GetNbinsX();
      int nBinY = rotatedMap->GetNbinsY();
      for (int binX=0; binX<nBinX+1; binX++) {
	double x = rotatedMap->GetXaxis()->GetBinCenter(binX+1);
	for (int binY=0; binY<nBinY+1; binY++) {
	  double y = rotatedMap->GetYaxis()->GetBinCenter(binY+1);
	  mapProfile->Fill(x,y,rotatedMap->GetBinContent(binX,binY));
	}
      }
      

      if (cnt%60 == 59) {
	mapProfile->Draw("colz");

	double sunTheta = eventSummary->sun.theta;
	double sunPhi = heading - eventSummary->sun.phi;
	while (sunPhi < 0) sunPhi += 360;
	gSun->SetPoint(0,sunPhi,-sunTheta);
	gSun->Draw("pSame");
	

	name.str("");
	name << "avgMaps/" << prefix << "_";
	name << setfill('0') << setw(numZeros) << cnt/60 << ".png";
	c1->SaveAs(name.str().c_str());
      }

      delete rotatedMap;

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
