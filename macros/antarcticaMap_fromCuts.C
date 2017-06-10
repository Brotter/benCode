#include "AnitaEventSummary.h"
#include "AntarcticaMapPlotter.h"

void antarcticaMap_fromCuts() {

  
  TChain *summaryTree = new TChain("summaryTree");

  summaryTree->Add("/Volumes/ANITA3Data/bigAnalysisFiles/templateSearch/06.05.17_14h/all.root");
  
  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  double templateValueH,templateValueV;
  summaryTree->SetBranchAddress("templateValueH",&templateValueH);
  summaryTree->SetBranchAddress("templateValueV",&templateValueV);

  summaryTree->BuildIndex("eventNumber");

  Acclaim::AntarcticaMapPlotter *aMap = new Acclaim::AntarcticaMapPlotter("mapPlotter","mapPlotter",1000,1000);
  aMap->addHistogram("pointMap","pointMap",250,250);


  int lenEntries = summaryTree->GetEntries();
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry % 1000 == 0) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);

    if ( templateValueH > 0.7 && !summary->flags.isVPolTrigger && summary->flags.isRF ) {

      double lat = summary->peak[0][0].latitude;
      double lon = summary->peak[0][0].longitude;
      if (lat < -999 || lon < -999) continue;

      aMap->Fill(lat,lon);
    }    


  }

  aMap->DrawHist("colz");

  return;
}
