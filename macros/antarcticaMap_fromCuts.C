#include "AnitaEventSummary.h"
#include "AntarcticaMapPlotter.h"

void antarcticaMap_fromCuts(string date="06.11.17_19h/") {

  
  TChain *summaryTree = new TChain("summaryTree");


  char* resultsDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  for (int run=130; run<440; run++) {
    name.str("");
    name << resultsDir << date << run << ".root";
    summaryTree->Add(name.str().c_str());
  }
  
  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  double templateCRayH[10];
  double templateCRayV[10];

  summaryTree->SetBranchAddress("templateCRayH",&templateCRayH);
  summaryTree->SetBranchAddress("templateCRayV",&templateCRayV);

  summaryTree->BuildIndex("eventNumber");

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter("mapPlotter","mapPlotter",1000,1000);
  aMap->addHistogram("pointMap","pointMap",250,250);


  int lenEntries = summaryTree->GetEntries();
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry % 10000 == 0) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);

    if ( templateCRayH[5] > 0.5 
	 && summary->peak[0][0].value > 0.06
	 && summary->flags.isRF ) {

      double lat = summary->peak[0][0].latitude;
      double lon = summary->peak[0][0].longitude;
      if (lat < -999 || lon < -999) continue;

      aMap->Fill(lat,lon);
    }    


  }

  aMap->DrawHist("colz");

  return;
}
