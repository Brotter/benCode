#include "AntarcticaMapPlotter.h"
#include "AnitaEventSummary.h"
#include "AnitaConventions.h"
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


void fillBaseList(Acclaim::AntarcticaMapPlotter *aMap,TGraph *gBaseList) {

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

void drawEvOnAntarctica(string fileName="cuts.root") {

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
  
  vector<int> goodEvs;
  goodEvs.push_back(11116669);
  goodEvs.push_back(11989349);
  goodEvs.push_back(15717147);
  goodEvs.push_back(16952229);
  goodEvs.push_back(19459851);
  goodEvs.push_back(23695286);
  goodEvs.push_back(32907848);
  goodEvs.push_back(33484995);
  goodEvs.push_back(41529195);

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


    //plot the event source location
    double xEv,yEv;
    double latEv = summary->peak[0][0].latitude;
    double lonEv = summary->peak[0][0].longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    allEvs->SetPoint(allEvs->GetN(),xEv,yEv);

    //fill the position graph with ANITA's location
    double xA,yA;
    latEv = summary->anitaLocation.latitude;
    lonEv = summary->anitaLocation.longitude;
    aMap->getRelXYFromLatLong(latEv,lonEv,xA,yA);
    anitaPosition->SetPoint(entry,xA,yA);

    TArrow *currArrow = new TArrow(xA,yA,xEv,yEv,0.003,"|>");
    arrows.push_back(currArrow);

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

      TArrow *highlightArrow = (TArrow*)currArrow->Clone();
      highlightArrow->SetLineWidth(2);
      highlightArrow->SetLineColor(kGreen);
      highlightArrow->SetFillColor(kGreen);
      highlightArrows.push_back(highlightArrow);
    }



  }

  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(1);
  gBaseList->SetMarkerColor(kRed);
  fillBaseList(aMap,gBaseList);


  cout << "Found " << eventCount << " entries " << endl;

  //  gBaseList->Draw("psame");


  aMap->setCurrentTGraph("anitaPosition");
  aMap->DrawTGraph("p");

  for (int i=0; i<arrows.size(); i++) {
    arrows[i]->Draw("");
  }

  aMap->setCurrentTGraph("eventLocations");
  aMap->DrawTGraph("psame");

  for (int i=0; i<highlights.size(); i++) {
    highlights[i]->Draw("psame");
    highlightArrows[i]->Draw();
  }





  return;
}



