#include "AnitaConventions.h"


void runMagicDisplayRun(int run);

void runMagicDisplay(int run=1002) {
  cout << "runMagicDisplay(" << run << ")\n";



  TChain *fred=0; //Will this work?
  runMagicDisplayRun(run);
}


void runMagicDisplayRun(int run) {

  MagicDisplay *magicPtr = new MagicDisplay("/Volumes/ANITA3Data/root",run,WaveCalType::kDefault);

  //magicPtr->startSurfDisplay();
  //  magicPtr->startAvgSurfDisplay();
  //   magicPtr->startTurfDisplay();
  //  magicPtr->startSumTurfDisplay();
  magicPtr->startEventDisplay();
  // magicPtr->applyCut("triggerTimeNs>999.9e6 && triggerTimeNs<1000e6 && eventNumber%2==1");
  //  magicPtr->applyCut("eventNumber%2==1");

  //  magicPtr->startGpsDisplay();
  //  magicPtr->startControlPanel();

}
