#include "AntarcticaMapPlotter.h"
#include "AnitaEventSummary.h"
#include "AnitaConventions.h"

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



void drawOnAntarctica(string fileName="passingAwk.txt") {

  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  string date="06.11.17_19h/";

  TChain *summaryTree = new TChain("summaryTree","summaryTree");

  stringstream name;
  for (int run=130;run<440;run++) {
    name.str("");
    name << resultsDir << date << run << ".root";
    summaryTree->Add(name.str().c_str());
  }
  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);

  cout << "There are " << summaryTree->GetEntries() << " entries in the summary tree" << endl;

  cout << "Building index..."; 
  fflush(stdout);
  summaryTree->BuildIndex("eventNumber");
  cout << " done!" << endl;

  
  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter();
  aMap->addTGraph("event","event");
  TGraph *gEv = aMap->getCurrentTGraph();
  gEv->SetMarkerStyle(29); //29=star
  gEv->SetMarkerColor(kGreen);

  TH2D *myHist = aMap->addHistogram("hist","hist",1000,1000);
  

  TGraph *passingEvs = new TGraph(fileName.c_str(),"%lg %lg*");
  cout << "found " << passingEvs->GetN() << " passing events in " << fileName << endl;

  for (int ev=0; ev<passingEvs->GetN(); ev++) {
    int summaryEntry = summaryTree->GetEntryNumberWithBestIndex(passingEvs->GetY()[ev]);
    cout << ev << " summaryEntry=" << summaryEntry << endl;
    summaryTree->GetEntry(summaryEntry);

    double xEv,yEv;
    double latEv = summary->peak[0][0].latitude;
    double lonEv = summary->peak[0][0].longitude;
    cout << "event position: " << latEv << " , " << lonEv << endl;
    aMap->getRelXYFromLatLong(latEv,lonEv,xEv,yEv);
    gEv->SetPoint(ev,xEv,yEv);
    aMap->Fill(latEv,lonEv);
  }

  aMap->addTGraph("baseList","baseList");
  TGraph *gBaseList = aMap->getCurrentTGraph();
  gBaseList->SetMarkerStyle(20); //20=filled circle
  gBaseList->SetMarkerSize(1);
  gBaseList->SetMarkerColor(kRed);

  char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  name.str("");
  name << anitaInstallDir << "/share/anitaCalib/baseListA3.root";
  cout << "looking for base list in " << name.str() << endl;
  TFile *fBaseList = TFile::Open(name.str().c_str());
  TTree *baseCampTree = (TTree*)fBaseList->Get("baseCampTree");
  double fullLat,fullLong;
  baseCampTree->SetBranchAddress("fullLat",&fullLat);
  baseCampTree->SetBranchAddress("fullLong",&fullLong);

  cout << "baseCampTree->GetEntries() = " << baseCampTree->GetEntries() << endl;
  for (int entry=0; entry<baseCampTree->GetEntries(); entry++) {
    baseCampTree->GetEntry(entry);
    double x,y;
    aMap->getRelXYFromLatLong(fullLat,fullLong,x,y);
    gBaseList->SetPoint(entry,x,y);
  }

  aMap->setCurrentTGraph("event");
  aMap->DrawTGraph("pSame");
  gBaseList->Draw("pSame");




  return;
}
