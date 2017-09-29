/*

  cluster.C is super long, and I'm not doing base clustering anymore really


  This contains all the stuff that relates events to minor bases!


 */


//probably the only useful thing in here, since I'll need it for WAIS and LDB
double calcBaseDistance(AnitaEventSummary *event, UsefulAdu5Pat *gps,double lat, double lon,double alt,
			bool debug=false) {
  /*

    For calculating the clustering to bases which only have one definate known location

   */
  const double Re = 6371e3; //km
  double distToHorizon = TMath::Sqrt(event->anitaLocation.altitude*(2*Re+event->anitaLocation.altitude));

  if (gps->getDistanceFromSource(lat,lon,alt) > distToHorizon) {
    //    if (debug) cout << "calcBaseDistance(): too far, " << gps->getDistanceFromSource(lat,lon,alt) << endl;
    return -9999;
  }


  //where the base is
  double thetaBase,phiBase;
  gps->getThetaAndPhiWave(lon,lat,alt,thetaBase,phiBase);
  thetaBase *= TMath::RadToDeg();
  phiBase   *= TMath::RadToDeg();

  //where the event is
  double phi =  event->peak[0][0].phi;
  double theta = event->peak[0][0].theta;
  
  //difference between A->a and A->b
  if (debug) cout << "theta:" << theta << " - " << thetaBase;
  if (debug) cout << " | phi:" << phi << " - " << phiBase;
  double diffTheta = TMath::Abs(theta - thetaBase);
  double diffPhi = TMath::Abs(FFTtools::wrap(phi - phiBase,360,0));
  if (diffPhi > 180) diffPhi = 360 - diffPhi;
  

  const double sigmaTheta = 0.2193;
  const double sigmaPhi = 0.5429;
      
  double diff = TMath::Sqrt(pow(diffTheta/sigmaTheta,2) + pow(diffPhi/sigmaPhi,2));
      
  if (debug) cout << " | diff: " << diff << endl;
  
  return diff;

}




/*
  Clustering with bases

*/

TTree *getBaseTree() {

  stringstream name;

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");

  return baseCampTree;
}

TTree *getBaseTree(double *lat, double *lon, double *alt, string** baseName) {
  TTree* baseTree = getBaseTree();
  baseTree->SetBranchAddress("fullLat",lat);
  baseTree->SetBranchAddress("fullLong",lon);
  baseTree->SetBranchAddress("alt",alt);
  baseTree->SetBranchAddress("name",baseName);
  return baseTree;
}
  

TGraph* getBaseGraph() {
  /*
    First get the bases!

    Returns them as a TGraph (lat,lon) that is supposed to be easy to go through

    This is dumber than just doing the trees
    
  */


  stringstream name;

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  TTree *awsCampTree = (TTree*)fBaseList->Get("awsTree");
  double fullLat,fullLon;
  double fullLataws,fullLongaws;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLon);

  awsCampTree->SetBranchAddress("fullLat",&fullLataws);
  awsCampTree->SetBranchAddress("fullLong",&fullLongaws);

  TGraph *gBaseList = new TGraph();
  gBaseList->SetName("gBaseList_latlon");

  cout << "baseCampTree->GetEntries() = " << baseCampTree->GetEntries() << endl;
  for (int entry=0; entry<baseCampTree->GetEntries(); entry++) {
    baseCampTree->GetEntry(entry);
    gBaseList->SetPoint(entry,fullLat,fullLon);
  }


  cout << "awsCampTree->GetEntries() = " << awsCampTree->GetEntries() << endl;
  for (int entry=0; entry<awsCampTree->GetEntries(); entry++) {
    awsCampTree->GetEntry(entry);
    gBaseList->SetPoint(gBaseList->GetN(),fullLat,fullLon);
  }

  return gBaseList;
}




void saveEventsNearBases(double threshold=40.,int numSplits=1,int split=0, string outFileName="") {
/*
  This will take up some space, but otherwise histograms take forever to generate

  Save every event that passes a clustering threshold, default was 40 but probably should be different, near a recorded
  base from jruss's list, or on one of the events that passed all the other cuts (which I also should probably tune!!)

  Basically identical to saveEventsNearCandidates


 */
  stringstream name,title;


  //get ALL the events
  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C(\"07.28.17_17h/\",false)");

  int lenEntries = summaryTree->GetEntries();
  cout << lenEntries << " total entries found" << endl;

  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  AnitaNoiseSummary *noiseSummary = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSummary);

  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);


  //get the bases
  TTree* baseTree = getBaseTree();
  double lat,lon,alt;
  baseTree->SetBranchAddress("fullLat",&lat);
  baseTree->SetBranchAddress("fullLong",&lon);
  baseTree->SetBranchAddress("alt",&alt);
  string *baseName = NULL;
  baseTree->SetBranchAddress("name",&baseName);

  int numBases = baseTree->GetEntries();

  cout << "Found " << numBases << " Bases~" << endl;

  //make an output file and trees, and save all the candidate info
  if (outFileName == "") outFileName = "baseClustering.root";
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");
  outTree->Branch("eventSummary",&eventSummary);
  outTree->Branch("template",&templateSummary);
  outTree->Branch("noiseSummary",&noiseSummary);
  outTree->Branch("gpsEvent",&gps);

  //and some new stuff that is important to remember
  int baseNum;
  double sigmaDist;
  outTree->Branch("baseNum",&baseNum);
  outTree->Branch("sigmaDist",&sigmaDist);

  //histograms to record distance distributions past threshold
  TH1D *hCluster[numBases];
  for (int base=0; base<numBases; base++) {
    baseTree->GetEntry(base);
    name.str("");
    name << "hCluster_" << base;
    title.str("");
    title << "Event Cluster Distance to " << *baseName;
    hCluster[base] = new TH1D(name.str().c_str(),title.str().c_str(),1000,0,threshold*4);
  }


  int close[numBases];
  for (int base=0; base<numBases; base++) {
    close[base] = 0;
  }


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;



  int startEntry,stopEntry;
  if (numSplits == 1) {
    startEntry=0;
    stopEntry=lenEntries;
  }
  else {
    lenEntries /= numSplits;
    startEntry = split*lenEntries;
    stopEntry = (split+1)*lenEntries;
    cout << "Splitting into " << numSplits << " sections, which means " << lenEntries << " events per section" << endl;
    cout << "Doing section: " << split << ", starting at entry " << startEntry << " and stopping at " << stopEntry << endl;
  }

  
  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10000 == 0) {
      int printEntry = entry-startEntry;
      cout << printEntry << "/" << lenEntries << " ( ";	
      for (int base=0; base<numBases; base++) {
	cout << close[base] << " ";
      }
      cout << ") ";
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(printEntry)/totalTimeSec;
      double remaining = (float(lenEntries-printEntry)/rate)/60.;
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left" << endl;
      watch.Start();
    }
    
    summaryTree->GetEntry(entry);
    
    if (notNotable(eventSummary)) continue;

    //where event b was captured from (B)
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);


    bool fill=false;
    
    int closestBase = -1;
    double closestBaseDist  = -1;

    //loop through bases
    for (int base=0; base<numBases; base++) {
      baseTree->GetEntry(base);
      double tempAlt = alt;
      if (tempAlt<0) tempAlt=0;

      double dist = calcBaseDistance(eventSummary,usefulGPS,lat,lon,tempAlt);

      hCluster[base]->Fill(dist);

      if (dist <= threshold && dist != -9999) {
	close[base]++;
	fill=true;
	if (closestBase==-1 || closestBaseDist > dist) {
	  closestBase = base;
	  closestBaseDist = dist;
	}
      }
    }
    if (fill) {
      baseNum = closestBase;
      sigmaDist = closestBaseDist;
      outTree->Fill();
    }

    delete usefulGPS;
  }      

  cout << "Finished!" << endl;

  for (int base=0; base<numBases; base++) {
    hCluster[base]->Write();  
  }
  outTree->Write();
  outFile->Close();

  return;

}




void mergeBaseClusterBoth(int numCores=32,int numBases=104,string date="08.28.17_22h") {
  /*
    I do clustering with the cluster servers so I gotta merge the results by hand
   */


  stringstream name;

  TChain *eventSummary = new TChain("summaryTree","summaryTree");


  TH1D *hCluster[numBases];
  TList *histList[numBases];
  for (int base=0; base<numBases; base++) {
    histList[base] = new TList;
  }

  char* basedir = getenv("ANITA3_RESULTSDIR");

  for (int i=0; i<numCores; i++) {
    name.str("");
    name << basedir << "cluster/" << date << "/baseCluster_" << i << ".root";
    cout << "loading: " << name.str() << endl;

    eventSummary->Add(name.str().c_str());

    TFile *inFile = TFile::Open(name.str().c_str());

    for (int base=0; base<numBases; base++) {
      name.str("");
      name << "hCluster_" << base;
      TH1D *currHist = (TH1D*)inFile->Get(name.str().c_str());
      if (i==0) hCluster[base] = (TH1D*)currHist->Clone();
      cout << name.str() << " " << currHist->GetEntries() << endl;
      histList[base]->Add(currHist);
    }

  }

  cout << "merging" << endl;
  for (int base=0; base<numBases; base++) {
    cout << "base:" << base << endl;
    hCluster[base]->Reset();
    hCluster[base]->Merge(histList[base]);

  }
  
    
  TFile *outFile = TFile::Open("mergeBaseClusters.root","recreate");
  for (int base=0; base<numBases; base++) {
    hCluster[base]->Write();
  }
  
  cout << "Copying summary" << endl;
  eventSummary->CloneTree(-1,"fast");
  outFile->Write();
  

  outFile->Close();

  return;
}




void printCandidateVsBases(string inFileName="candidates.root",bool debug=false) {
  /*
    print out whether the candidates cluster with any bases
   */

  TFile *inFile = TFile::Open(inFileName.c_str());
  if (!inFile->IsOpen()) {
    cout << "Coudln't open " << inFileName << endl;
    return;
  }

  TTree* summaryTree = (TTree*)inFile->Get("summaryTree");
  if (summaryTree == NULL) {
    cout << "File didn't have a summaryTree! quitting" << endl;
    return;
  }


  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  
  int lenEntries = summaryTree->GetEntries();
  cout << "Found " << lenEntries << " candidates" << endl;

  /* Old way, which skips a bunch I think
  double lat,lon,alt;
  string *baseName = NULL;
  TTree* baseTree = getBaseTree(&lat,&lon,&alt,&baseName);
  */

  //new way to get bases from Cosmin!
  BaseList::makeBaseList();


  int count = 0;

  stringstream name;
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    double closest = 9999;
    double closestNum = -1;
    UsefulAdu5Pat *useful = new UsefulAdu5Pat(gps);

    for (int baseNum=0; baseNum < BaseList::getNumBases(); baseNum++) {
      BaseList::base base = BaseList::getBase(baseNum);

      if (!base.isValid(evSum->realTime)) continue;

      AntarcticCoord coord = base.getPosition(evSum->realTime);
      coord.to(AntarcticCoord::CoordType::WGS84);
      double dist = calcBaseDistance(evSum,useful,coord.x,coord.y,coord.z,debug);
      if (closest > dist && dist != -9999) {
	if (debug) {
	  cout << "evCoords:" << evSum->anitaLocation.latitude << " " << evSum->anitaLocation.longitude << endl;
	  cout << "baseCoords:" << coord.x << " " << coord.y << " " <<  coord.z << endl;}
	closest = dist;
	closestNum = baseNum;
	name.str("");
	name << base.getName();
      }
    }

    for (int baseNum=0; baseNum < BaseList::getNumPaths(); baseNum++) {
      BaseList::path base = BaseList::getPath(baseNum);

      if (!base.isValid(evSum->realTime)) continue;
      
      AntarcticCoord coord = base.getPosition(evSum->realTime);
      coord.to(AntarcticCoord::CoordType::WGS84);
      double dist = calcBaseDistance(evSum,useful,coord.x,coord.y,coord.z,debug);
      if (closest > dist && dist != -9999) {
	if (debug) {
	  cout << "evCoords:" << evSum->anitaLocation.latitude << " " << evSum->anitaLocation.longitude << endl;
	  cout << "baseCoords:" << coord.x << " " << coord.y << " " <<  coord.z << endl; }
	closest = dist;
	closestNum = baseNum;
	name.str("");
	name << base.getName();
      }
    }

    if (closest > 40) count++;

    cout << "ev" << evSum->eventNumber << " : " << closest << " " << name.str() << " " << closestNum << endl;
    delete useful;
  }

	
  cout << "Only " << count << " were more than 40 away from a base" << endl;
  
  return;
}




