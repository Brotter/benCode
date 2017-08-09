#include "AnitaVersion.h"

/*


  Do you want to split out coherently sumemd wais pulses into a text file?  Look no further!


*/



void waisTxtSaver() {

  AnitaVersion::set(3);


  const int run=342; //a good wais run


  AnitaDataset *data = new AnitaDataset(342);

  
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 

  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); //interactive needs to be true

  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

  AnitaEventSummary *eventSummary = new AnitaEventSummary();

  const int numWaveforms = 100;

  TGraph *waveforms[numWaveforms];

  int numSavedWaveforms = 0;
  int entry=0;

  while (numSavedWaveforms < numWaveforms) {
    data->getEntry(entry);
    entry++;

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(data->gps());

    if (UCorrelator::isWAISHPol(usefulGPS,data->header(),config)) {
      
      FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), fStrat, data->gps(), data->header());

      analyzer->analyze(filteredEvent, eventSummary); 

      waveforms[numSavedWaveforms] = new TGraph(*analyzer->getCoherent(AnitaPol::kHorizontal,0,true)->even());

      numSavedWaveforms++;
      cout << numSavedWaveforms << " " << entry << endl;
      analyzer->clearInteractiveMemory();
      delete filteredEvent;
    }
    
    delete usefulGPS;

  }

  ofstream outFile("waisEvents.txt");
  for (int pt=0; pt<waveforms[0]->GetN(); pt++) {
    outFile << waveforms[0]->GetX()[pt] << " ";
    for (int wave=0; wave<numWaveforms; wave++) {
      if (waveforms[wave]->GetN() > pt) outFile << waveforms[wave]->GetY()[pt] << " ";
      else                              outFile <<              0              << " ";
    }
    outFile << endl;
  }
  outFile.close();

}
