#include "AntarcticaMapPlotter.h"
#include "AnitaEventSummary.h"
#include "AnitaConventions.h"
#include "BaseList.h"
#include "AntarcticaGeometry.h"
#include <Riostream.h>

/*

  I want to plot the events onto the continent, probably with the base list, so that I can see how likely they are an event


 */


int evToRun(int ev) {

  ifstream runToEv("/Users/brotter/anita16/benMacros/runToEv.txt");
  int run,evLow,evHigh;
  int runOut = -1;
  while (runToEv >> run >> evLow >> evHigh)
    if (ev >= evLow && ev <= evHigh) {
      runOut = run;
    }

  runToEv.close();
  
  return runOut;

}

void placeHiCalOnMap(Acclaim::AntarcticaMapPlotter *aMap,TGraph *gHiCal) {
  
  TFile *hiCalFile = TFile::Open("/Users/brotter/anita16/benCode/HiCal/anita_hical_bearing_all.root");
  TTree *hiCalTree = (TTree*)hiCalFile->Get("HCbearingTree");
  double lat,lon;
  hiCalTree->SetBranchAddress("hicalLat",&lat);
  hiCalTree->SetBranchAddress("hicalLon",&lon);
  
  
  double x,y;
  for (int entry=0; entry<hiCalTree->GetEntries(); entry++) {
    if (entry%10==0) {
      hiCalTree->GetEntry(entry);
      aMap->getRelXYFromLatLong(lat,lon,x,y);
      gHiCal->SetPoint(gHiCal->GetN(),x,y);
    }
  }

  gHiCal->SetMarkerColor(kBlue);

  return;
}


void placeBasesOnMap(Acclaim::AntarcticaMapPlotter *aMap,TGraph *gBaseList) {

  BaseList::makeBaseList();


  double x,y;
  for (int baseNum=0; baseNum < BaseList::getNumBases(); baseNum++) {
    BaseList::base base = BaseList::getBase(baseNum);
    AntarcticCoord coord = base.getPosition(0);
    aMap->getRelXYFromLatLong(coord.x,coord.y, x,y);
    gBaseList->SetPoint(gBaseList->GetN(),x,y);
  }

  return;
}


void placeBasesOnMap_old(Acclaim::AntarcticaMapPlotter *aMap,TGraph *gBaseList) {

  stringstream name;

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  TTree *awsCampTree = (TTree*)fBaseList->Get("awsTree");
  double fullLat,fullLong;
  double fullLataws,fullLongaws;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLong);

  awsCampTree->SetBranchAddress("fullLat",&fullLataws);
  awsCampTree->SetBranchAddress("fullLong",&fullLongaws);

  cout << "baseCampTree->GetEntries() = " << baseCampTree->GetEntries() << endl;
  for (int entry=0; entry<baseCampTree->GetEntries(); entry++) {
    baseCampTree->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(fullLat,fullLong,x,y);
    gBaseList->SetPoint(entry,x,y);
  }


  cout << "awsCampTree->GetEntries() = " << awsCampTree->GetEntries() << endl;
  for (int entry=0; entry<awsCampTree->GetEntries(); entry++) {
    awsCampTree->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(fullLataws,fullLongaws,x,y);
    gBaseList->SetPoint(gBaseList->GetN(),x,y);
  }

  return;
}

void drawCandidatesOnAntarctica_hardcoded(string fileName="cuts.root") {

  /*

    Draws all events from a bunch of AnitaEventSummaries in a file

   */

  stringstream name;

  TFile *inFile = TFile::Open(fileName.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("eventSummary");

  if (summaryTree == NULL) {
    cout << "Couldn't find cutSummary in file " << fileName << endl;
  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  cout << "found " << summaryTree->GetEntries() << " events in " << fileName << endl;

  
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  aMap->addTGraph("anitaPosition","anitaPosition");
  TGraph *anitaPosition = aMap->getCurrentTGraph();
  anitaPosition->SetMarkerColor(kWhite);

  aMap->addTGraph("eventLocations","eventLocations");
  TGraph *allEvs = aMap->getCurrentTGraph();
  allEvs->SetMarkerStyle(29); //star
  allEvs->SetMarkerSize(1);
  allEvs->SetMarkerColor(kOrange);

  vector<TArrow*> arrows;
  

  //from clustering results
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
    
  vector<TArrow*> highlightArrows;
  vector<TGraph*> highlights;


  int eventCount=0;
  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);

    //cut on dense event numbers
    //    if (summary->eventNumber < 10e6 || summary->eventNumber > 50e6) continue;

    //cut on things that don't point at the continent
    if (summary->peak[0][0].latitude < -999) continue;

    //    if (TMath::Abs(summary->deconvolved[0][0].linearPolAngle()) > 20) continue;
    //    if (TMath::Abs(summary->deconvolved[0][0].linearPolFrac()) < 0.6) continue;

    eventCount++;

    double xA,yA,latEv,lonEv,xEv,yEv;

    //fill the position graph with ANITA's location
    latEv = summary->anitaLocation.latitude;
    lonEv = summary->anitaLocation.longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xA,yA);
    anitaPosition->SetPoint(entry,xA,yA);

    //plot the event source location
    latEv = summary->peak[0][0].latitude;
    lonEv = summary->peak[0][0].longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.003,"|>");

    //highlight graphs
    if (std::find(goodEvs.begin(),goodEvs.end(),summary->eventNumber)!=goodEvs.end()) {
      name.str("");
      name << "event" << summary->eventNumber;
      cout << name.str() << " " << latEv << "|" << lonEv << endl;
      TGraph *gHighlight = new TGraph();
      gHighlight->SetName(name.str().c_str());
      gHighlight->SetTitle(name.str().c_str());
      gHighlight->SetMarkerStyle(29); //29=star
      gHighlight->SetMarkerColor(kGreen);
      gHighlight->SetPoint(0,xEv,yEv);
      highlights.push_back(gHighlight);

      currArrow->SetLineWidth(2);
      currArrow->SetLineColor(kGreen);
      currArrow->SetFillColor(kGreen);
      highlightArrows.push_back(currArrow);
    }
    //clustered events
    else {
      allEvs->SetPoint(allEvs->GetN(),xEv,yEv);
      arrows.push_back(currArrow);
    }



  }
  
  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(1);
  gBaseList->SetMarkerColor(kRed);
  placeBasesOnMap(aMap,gBaseList);


  cout << "Found " << eventCount << " entries " << endl;



  aMap->setCurrentTGraph("anitaPosition");
  aMap->DrawTGraph("p");

  aMap->setCurrentTGraph("baseList");
  aMap->DrawTGraph("psame");

  for (int i=0; i<arrows.size(); i++) {
  arrows[i]->Draw("");
  }

  aMap->setCurrentTGraph("eventLocations");
  aMap->DrawTGraph("psame");

  
  for (int i=0; i<highlights.size(); i++) {
    highlightArrows[i]->Draw();
    highlights[i]->Draw("psame");
  }
  




  return;
}




void drawCandidatesOnAntarctica(string fileName="candidates.root",bool moreCuts = true) {

  /*

    Draws all events from a bunch of AnitaEventSummaries in a file

   */

  stringstream name;

  TFile *inFile = TFile::Open(fileName.c_str());
  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");

  if (summaryTree == NULL) {
    cout << "Couldn't find summaryTree in file " << fileName << endl;
  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);

  cout << "found " << summaryTree->GetEntries() << " events in " << fileName << endl;

  
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  aMap->addTGraph("anitaPosition","anitaPosition");
  TGraph *anitaPosition = aMap->getCurrentTGraph();
  anitaPosition->SetMarkerColor(kWhite);


  cout << "Getting HiCal" << endl;
  aMap->addTGraph("hiCal","hiCal");
  TGraph *gHiCal = aMap->getCurrentTGraph();
  placeHiCalOnMap(aMap,gHiCal);


  vector<TArrow*> arrows;

  vector<TArrow*> highlightArrows;
  vector<TGraph*> highlights;


  int count=0;
  for (int entry=0; entry<summaryTree->GetEntries(); entry++) {
    summaryTree->GetEntry(entry);

    if (moreCuts) {
      //blasts
      if (summary->flags.maxBottomToTopRatio[0] > 3) {
	cout << summary->eventNumber << " is a blast" << endl;
	continue;
      }
    
      //point at bases (from cluster.C) NOT DOING THIS ANYMORE!
      /*      if (summary->eventNumber == 11116669 || 
	  summary->eventNumber == 11989349 || 
	  summary->eventNumber == 16952229 || 
	  summary->eventNumber == 33484995 || 
	  summary->eventNumber == 58592863 ||
	  summary->eventNumber == 62273732 || 
	  summary->eventNumber == 63210848 ||
	  summary->eventNumber == 80561103 || 
	  summary->eventNumber == 83877990 || 
	  summary->eventNumber == 84114142) continue;
      */

      //hardware trigger angle
      if (TMath::Abs(summary->peak[0][0].hwAngle) > 45) continue;

      //geomagnetic > 10 (from geomagnetic.C)
      if (summary->eventNumber == 84114142 || summary->eventNumber == 84405480) continue;
      
     
    }



    count++;
    double xA,yA,latEv,lonEv,xEv,yEv;

    //fill the position graph with ANITA's location
    latEv = summary->anitaLocation.latitude;
    lonEv = summary->anitaLocation.longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xA,yA);
    anitaPosition->SetPoint(entry,xA,yA);

    //plot the event source location
    latEv = summary->peak[0][0].latitude;
    lonEv = summary->peak[0][0].longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.003,"|>");

    name.str("");
    name << "event" << summary->eventNumber;
    cout << name.str() << " " << latEv << "|" << lonEv << endl;
    TGraph *gHighlight = new TGraph();
    gHighlight->SetName(name.str().c_str());
    gHighlight->SetTitle(name.str().c_str());
    gHighlight->SetMarkerStyle(29); //29=star
    gHighlight->SetMarkerSize(1);
    gHighlight->SetMarkerColor(kGreen);
    gHighlight->SetPoint(0,xEv,yEv);
    highlights.push_back(gHighlight);
    
    currArrow->SetLineWidth(2);
    //    currArrow->SetLineColor(kGreen);
    //    currArrow->SetFillColor(kGreen);
    highlightArrows.push_back(currArrow);

  }
  
  cout << count << " were not blasts and didn't point at bases" << endl;

  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(1);
  gBaseList->SetMarkerColor(kRed);
  placeBasesOnMap(aMap,gBaseList);


  aMap->setCurrentTGraph("anitaPosition");
  aMap->DrawTGraph("p");

  aMap->setCurrentTGraph("baseList");
  //  aMap->DrawTGraph("psame");

  aMap->setCurrentTGraph("hiCal");
  aMap->DrawTGraph("psame");
  
  for (int i=0; i<highlights.size(); i++) {
    highlights[i]->Draw("psame");
    highlightArrows[i]->Draw();

  }
  




  return;
}



void drawFlightPathOnAntarctica() {

  /*
    I need a plot of the flight path in chapter two
   */
  stringstream name;

  char* dataDir = getenv("ANITA_ROOT_DATA");

  
  TChain* gpsTree = new TChain("adu5PatTree","adu5PatTree");
  for (int i=120; i<440; i++) {
    name.str("");
    name << dataDir << "/run" << i << "/gpsFile" << i << ".root";
    gpsTree->Add(name.str().c_str());
  }
  Adu5Pat *gps = NULL;
  gpsTree->SetBranchAddress("pat",&gps);

  int lenEntries = gpsTree->GetEntries();


  double x,y;

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  aMap->addTGraph("anitaPosition","anitaPosition");
  TGraph *anitaPosition = aMap->getCurrentTGraph();
  //  anitaPosition->SetMarkerColor(kWhite);
  
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << "(" << (100.*entry)/lenEntries << "%)" << endl;
    gpsTree->GetEntry(entry);
    //    if (eventSummary->flags.isRF) continue;
    //    aMap->getRelXYFromLatLong(eventSummary->anitaLocation.latitude,eventSummary->anitaLocation.longitude,x,y);
    aMap->getRelXYFromLatLong(gps->latitude,gps->longitude,x,y);    
    
    anitaPosition->SetPoint(entry,x,y);
  }



  aMap->DrawTGraph("p");

}


void drawBaseListOnAntarctica() {

  /*
    I need a base map that tells me the name, and number, for each base so I can figure out what is what

    unfinished

   */
  stringstream name;

  
  //  vector <TGraph*> gBaseList;

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  aMap->addTGraph("baseList","baseList");
  TGraph *currGraph = aMap->getCurrentTGraph();
  TGraph *gBaseList = currGraph;
  currGraph->SetMarkerStyle(20); //20=filled circle
  currGraph->SetMarkerSize(1);
  currGraph->SetMarkerColor(kRed);
  placeBasesOnMap(aMap,gBaseList);
  aMap->DrawTGraph("p");

}




void drawClusteredBases(bool thermalCut=false) {

  /*

    I have the code that clusters things to known bases, so lets plot those onto the continent too!

  */


  
  TFile *inFile = TFile::Open("mergeBaseClusters.root");
  TTree* summaryTree = (TTree*)inFile->Get("summaryTree");
  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);
  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  int lenEntries = summaryTree->GetEntries();


  double x,y;

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  aMap->addHistogram("baseCluster","baseCluster",1000,1000);
  TH2D *baseCluster = aMap->getCurrentHistogram();
  aMap->addTGraph("anitaPosition","anitaPosition");
  TGraph *anitaPosition = aMap->getCurrentTGraph();
  anitaPosition->SetMarkerColor(kWhite);
  
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << "(" << (100.*entry)/lenEntries << "%)" << endl;
    summaryTree->GetEntry(entry);

    if (thermalCut) {
      if (eventSummary->peak[0][0].value < 0.0055 ||
	  eventSummary->peak[0][0].snr < 1.35 ||
	  eventSummary->deconvolved_filtered[0][0].peakHilbert < 11.5 ||
	  templateSummary->coherent[0][0].cRay[4] < 0.115) { 
	continue;
      } }

    aMap->getRelXYFromLatLong(eventSummary->peak[0][0].latitude,eventSummary->peak[0][0].longitude,x,y);    
    baseCluster->Fill(x,y);

    aMap->getRelXYFromLatLong(eventSummary->anitaLocation.latitude,eventSummary->anitaLocation.longitude,x,y);
    anitaPosition->SetPoint(entry,x,y);
  }


  aMap->setCurrentHistogram("baseCluster");
  aMap->DrawHist("colz");

  aMap->setCurrentTGraph("anitaPosition");
  aMap->DrawTGraph("psame");


  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(0.2);
  gBaseList->SetMarkerColor(kRed);
  placeBasesOnMap(aMap,gBaseList);


  aMap->setCurrentTGraph("baseList");
  aMap->DrawTGraph("psame");


}


void plotBasesWithClusteredEvents() {

  /*

    I generated a huge list of all events that cluster with bases, so lets plot those onto Antarctica!

   */
  stringstream name;


  //file with all the events
  TFile *summaryFile = TFile::Open("mergeBaseClusters.root");
  TTree* summaryTree = (TTree*)summaryFile->Get("summaryTree");
  int lenEntries = summaryTree->GetEntries();
  cout << "Got summaryTree with " << lenEntries << " entries" << endl;
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  int baseNum;
  summaryTree->SetBranchAddress("baseNum",&baseNum);

  //get bases
  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  double fullLat,fullLong;
  string *baseName = NULL;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLong);
  baseCampTree->SetBranchAddress("name",&baseName);

  //make an Antarctica map
  Acclaim::AntarcticaMapPlotter *aMap;
  double x,y;

  //loop through bases
  int lenBases = baseCampTree->GetEntries();
  for (int base=0; base<lenBases; base++) {
    cout << "Base " << base << endl;

    baseCampTree->GetEntry(base);

    name.str("");
    name << "baseNum == " << base;
    if (summaryTree->Draw("baseNum>>hist",name.str().c_str(),"goff") == 0) continue;
    

    aMap = new Acclaim::AntarcticaMapPlotter();

    aMap->addTGraph("gBase",*baseName);
    TGraph *gBase = aMap->getCurrentTGraph();
    gBase->SetMarkerStyle(29);//star
    gBase->SetMarkerStyle(kRed);
    gBase->SetMarkerSize(1);
    aMap->getRelXYFromLatLong(fullLat,fullLong,x,y);
    gBase->SetPoint(0,x,y);

    aMap->addTGraph("anitaPosition","anitaPosition");
    TGraph *anitaPosition = aMap->getCurrentTGraph();
    anitaPosition->SetMarkerColor(kWhite);

    aMap->addHistogram("hBase",*baseName,1000,1000);
    TH2D* hBase = aMap->getCurrentHistogram();
    for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << "(" << (100.*entry)/lenEntries << "%)" << endl;
      summaryTree->GetEntry(entry);
      if (baseNum != base) continue;
      aMap->getRelXYFromLatLong(evSum->peak[0][0].latitude,evSum->peak[0][0].longitude,x,y);
      hBase->Fill(x,y);
      aMap->getRelXYFromLatLong(evSum->anitaLocation.latitude,evSum->anitaLocation.longitude,x,y);
      anitaPosition->SetPoint(anitaPosition->GetN(),x,y);

    }


    TCanvas *canvas =aMap->DrawHist("colz");
    aMap->setCurrentTGraph("anitaPosition");
    aMap->DrawTGraph("psame");
    aMap->setCurrentTGraph("hBase");
    aMap->DrawTGraph("psame");
    name.str("");
    name << "baseDists/baseMap_" << base << ".png";
    canvas->SaveAs(name.str().c_str());

    delete aMap;
  }

  return;
}

  

       
void drawDirectEvents(){
  /*
    Direct events are events too!  Lets just point them in some generalized direction I guess?
   */
  stringstream name;
  
  //make an Antarctica map
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  cout << "Getting base list" << endl;
  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  placeBasesOnMap(aMap,gBaseList);
  //  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerColor(kGreen);
  gBaseList->SetMarkerSize(0.05);
  cout << "got " << gBaseList->GetN() << " bases" << endl;


  cout << "Getting HiCal" << endl;
  aMap->addTGraph("hiCal","hiCal");
  TGraph *gHiCal = aMap->getCurrentTGraph();
  placeHiCalOnMap(aMap,gHiCal);

  vector<TGraph*> candidates;
  vector<TArrow*> directArrows;

  TFile *directFile = TFile::Open("aboveHorizon.root");
  if (!directFile->IsOpen()) {
    cout << "Couldn't open aboveHorizon.root" << endl;
    return;
  }

  TTree* directTree = (TTree*)directFile->Get("summaryTree");
  int lenEntries = directTree->GetEntries();
  cout << "Found " << lenEntries << " entries" << endl;
  AnitaEventSummary *directSum = NULL;
  directTree->SetBranchAddress("eventSummary",&directSum);
  Adu5Pat *directGPS = NULL;
  directTree->SetBranchAddress("gpsEvent",&directGPS);

  for (int entry=0; entry<lenEntries; entry++) {
    directTree->GetEntry(entry);

    cout << "entry:" << entry << endl;

    //2 sigma above horizon cut
    double Re = 6371e3;
    if (TMath::ACos(Re/(directSum->anitaLocation.altitude+Re))*TMath::RadToDeg() - directSum->peak[0][0].theta < 0.5) continue;

    double latEv,lonEv,altEv,theta_adj;
    double xEv,yEv,xA,yA;
    
    UsefulAdu5Pat *directUseful = new UsefulAdu5Pat(directGPS);
    int status = directUseful->traceBackToContinent(TMath::DegToRad()*directSum->peak[0][0].phi,
					      TMath::DegToRad()*6.5,
					      &lonEv,&latEv,&altEv,&theta_adj);

    double latA = directSum->anitaLocation.latitude;
    double lonA = directSum->anitaLocation.longitude;

    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    aMap->getRelXYFromLatLong(latA, lonA, xA, yA);

    cout << latEv << " " << lonEv << endl;
    cout << latA  << " " << lonA  << endl;
    cout << "x,y: " << xEv << "," << yEv << endl;

    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.01,"|>");
    //    currArrow->SetLineColor(entry);
    //    currArrow->SetFillColor(entry);
    currArrow->SetLineWidth(2);
    directArrows.push_back(currArrow);

    name.str("");
    name << "ev" << directSum->eventNumber;
    TGraph *currEv = new TGraph();
    currEv->SetName(name.str().c_str());
    currEv->SetPoint(0,xEv,yEv);
    candidates.push_back(currEv);

  }

  aMap->DrawTGraph("p");
  
  for (int i=0; i<directArrows.size(); i++) {
    directArrows[i]->Draw();
    candidates[i]->Draw("same");
  }
  
  return;
}
  

void drawDirectAndReflected() {
  /*
    
    do em both to see if there are weird clusters like that

   */
  stringstream name;


  //make an Antarctica map
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  cout << "Getting base list" << endl;
  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  placeBasesOnMap(aMap,gBaseList);
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerColor(kGreen);
  gBaseList->SetMarkerSize(1);
  cout << "got " << gBaseList->GetN() << " bases" << endl;


  vector<TGraph*> evNumGraphs;

  vector<TArrow*> directArrows;

  TFile *directFile = TFile::Open("aboveHorizon.root");
  if (!directFile->IsOpen()) {
    cout << "Couldn't open aboveHorizon.root" << endl;
    return;
  }

  TTree* directTree = (TTree*)directFile->Get("summaryTree");
  int lenEntries = directTree->GetEntries();
  cout << "Found " << lenEntries << " entries" << endl;
  AnitaEventSummary *directSum = NULL;
  directTree->SetBranchAddress("eventSummary",&directSum);
  Adu5Pat *directGPS = NULL;
  directTree->SetBranchAddress("gpsEvent",&directGPS);

  for (int entry=0; entry<lenEntries; entry++) {
    directTree->GetEntry(entry);

    cout << "entry:" << entry << endl;

    double latEv,lonEv,altEv,theta_adj;
    double xEv,yEv,xA,yA;
    
    UsefulAdu5Pat *directUseful = new UsefulAdu5Pat(directGPS);
    int status = directUseful->traceBackToContinent(TMath::DegToRad()*directSum->peak[0][0].phi,
					      TMath::DegToRad()*6.5,
					      &lonEv,&latEv,&altEv,&theta_adj);

    double latA = directSum->anitaLocation.latitude;
    double lonA = directSum->anitaLocation.longitude;

    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    aMap->getRelXYFromLatLong(latA, lonA, xA, yA);

    cout << latEv << " " << lonEv << endl;
    cout << latA  << " " << lonA  << endl;
    cout << "x,y: " << xEv << "," << yEv << endl;

    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.01,"|>");
    directArrows.push_back(currArrow);

    name.str("");
    name << "ev" << directSum->eventNumber;
    aMap->addTGraph(name.str().c_str(),name.str().c_str());
    TGraph *currGraph = aMap->getCurrentTGraph();
    currGraph->SetPoint(0,xA,yA);
    evNumGraphs.push_back(currGraph);

  }
   
  vector<TArrow*> reflectArrows;

  TFile *reflectFile = TFile::Open("candidates.root");
  if (!reflectFile->IsOpen()) {
    cout << "Couldn't open candidates.root" << endl;
    return;
  }

  TTree* reflectTree = (TTree*)reflectFile->Get("summaryTree");
  int reflectEntries = reflectTree->GetEntries();
  cout << "Found " << lenEntries << " entries" << endl;
  AnitaEventSummary *reflectSum = NULL;
  reflectTree->SetBranchAddress("eventSummary",&reflectSum);
  Adu5Pat *reflectGPS = NULL;
  reflectTree->SetBranchAddress("gpsEvent",&reflectGPS);

  for (int entry=0; entry<reflectEntries; entry++) {
    reflectTree->GetEntry(entry);

    double xA,yA,latEv,lonEv,xEv,yEv;

    //fill the position graph with ANITA's location
    latEv = reflectSum->anitaLocation.latitude;
    lonEv = reflectSum->anitaLocation.longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xA,yA);

    //plot the event source location
    latEv = reflectSum->peak[0][0].latitude;
    lonEv = reflectSum->peak[0][0].longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.003,"|>");
    currArrow->SetLineWidth(2);
    currArrow->SetLineColor(kRed);
    currArrow->SetFillColor(kRed);
    reflectArrows.push_back(currArrow);

    name.str("");
    name << "ev" << reflectSum->eventNumber;
    aMap->addTGraph(name.str().c_str(),name.str().c_str());
    TGraph *currGraph = aMap->getCurrentTGraph();
    currGraph->SetPoint(0,xA,yA);
    evNumGraphs.push_back(currGraph);


  }

  aMap->setCurrentTGraph("baseList");
  aMap->DrawTGraph("p");
  
  for (int i=0; i<directArrows.size(); i++) {
    directArrows[i]->Draw();
  }
  
  for (int i=0; i<reflectArrows.size(); i++) {
    reflectArrows[i]->Draw();
  }

  for (int i=0; i<evNumGraphs.size(); i++) {
    evNumGraphs[i]->SetMarkerStyle(kDiamond);
    evNumGraphs[i]->SetMarkerSize(1);
    evNumGraphs[i]->Draw("psame");
  }
  
  return;
}


void drawTemplateMap() {
  /*
    Draws a 2d histogram on antarctica with the template correlation value
   */


  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();

  aMap->addProfile("templateProf","tempalteProf",1000,1000);
  
  TFile *inFile = TFile::Open("mergeBaseClusters.root");
  TTree* summaryTree = (TTree*)inFile->Get("summaryTree");
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  AnitaTemplateSummary *tempSum = NULL;
  summaryTree->SetBranchAddress("template",&tempSum);


  int lenEntries = summaryTree->GetEntries();

  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0) cout << "Hi CRAB!!! " << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);
    aMap->Fill(evSum->peak[0][0].latitude,evSum->peak[0][0].longitude,tempSum->coherent[0][0].cRay[4]);
  }

  aMap->DrawHist("colz");

  return;
}

  
void drawHiCalOnAntarctica() {
  /*

    Jess sent me a bunch of HiCal stuff, so I should plot it for my thesis!

  */

  

  return;
}

    

void drawSingleEvOnAntarctica(int evNum, string inFileName = "") {
  /*
    Draw a single event on to the continent if you can
    if inFileName is filled, it loads that file looking for a summaryTree to get the event from.
    Otherwise it loads the whole thing.
  */
  
  TChain *summaryTree;
  if (inFileName == "") {
    summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");
  }
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(inFileName.c_str());
  }

  summaryTree->BuildIndex("eventNumber");
  int entry=summaryTree->GetEntryNumberWithIndex(evNum);
  if (entry==-1) {
    cout << "Couldn't find event number " << evNum << ". Quitting" << endl;
    return;
  }

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  
  summaryTree->GetEntry(entry);

  AntarcticaBackground *aBkgd = new AntarcticaBackground();
  
  TGraphAntarctica *gEv = new TGraphAntarctica();
  gEv->SetPoint(0,evSum->peak[0][0].longitude,evSum->peak[0][0].latitude);
  
  
  
  aBkgd->Draw();
  gEv->Draw("psame");

  return;
}



void checkGPS(int evNum, string inFileName = "") {
  /*

    Cosmin fucked up getSourceLatLon or something stupid

    if inFileName is filled, it loads that file looking for a summaryTree to get the event from.
    Otherwise it loads the whole thing.
  */
  
  TChain *summaryTree;
  if (inFileName == "") {
    summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");
  }
  else {
    summaryTree = new TChain("summaryTree","summaryTree");
    summaryTree->Add(inFileName.c_str());
  }

  summaryTree->BuildIndex("eventNumber");
  int entry=summaryTree->GetEntryNumberWithIndex(evNum);
  if (entry==-1) {
    cout << "Couldn't find event number " << evNum << ". Quitting" << endl;
    return;
  }

  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gps);

  summaryTree->GetEntry(entry);

  UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gps);
  double lat,lon,alt,theta_adj;
  int hitsCont = usefulGPS->traceBackToContinent(TMath::DegToRad()*evSum->peak[0][0].phi,
						 TMath::DegToRad()*evSum->peak[0][0].theta,
						 &lon,&lat,&alt,&theta_adj);
  cout << "hitsCont: " << hitsCont << endl;
  cout << "eventSum: " << evSum->peak[0][0].latitude << " : " << evSum->peak[0][0].longitude << " : " << evSum->peak[0][0].altitude << endl;
  cout << "now: " << lat << " : " << lon << " : " << alt << endl;
  
  
  return;
}

  


void drawEvOnAntarctica() {
  cout << "loaded drawEvOnAntarctica.C" << endl;
}
