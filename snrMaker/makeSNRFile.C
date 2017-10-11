/*

  Whoops!  I need to pre-compile this list before running over things because the wais events aren't in order

 */

#include "AnitaConventions.h"
#include "AntennaPositions.h"
#include "UsefulAdu5Pat.h"

#include "loadAll.C"

string date="09.27.17_19h";


void makeSNRFileSplit(int numSplits,int split) {
  bool debug=false;
  /*

    This should go through the AnitaNoiseSummaries, and in conjunction with the peak pointing direction and
    peak of the waveform (I guess?), determine the SNR of that waveform.

    Problem: I don't save the minimum!  I can't do peak to peak DAMNIT.  I'll have to fix that

   */


  TChain *summaryTree = loadAll(date,false);
  int lenEntries = summaryTree->GetEntries();

  //I'll need to split this up onto the servers, because it takes forever
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


  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gpsEvent = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gpsEvent);
  AnitaNoiseSummary *noiseSum = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);


  char* dataDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  name.str("");
  name << "data/SNRs";
  if (numSplits>1) name << "_" << split;
  name << ".root";
  cout << "using " << name.str() << " as output file" << endl;
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  TTree *snrTree = new TTree("snrTree","snrTree");
  int eventNumber;
  snrTree->Branch("eventNumber",&eventNumber,"eventNumber/I");
  double snr,snr_filtered,rms;
  int peakPhiSector;
  snrTree->Branch("snr",&snr,"snr/D");
  snrTree->Branch("snr_filtered",&snr_filtered,"snr_filtered/D");
  snrTree->Branch("rms",&rms,"rms/D");
  snrTree->Branch("peakPhiSector",&peakPhiSector);

  const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 
  const int antennasInSum = 15; //I'm using 15 antennas
  const int phiSectorsInSum = antennasInSum/3;
  int closest[antennasInSum];

  AnitaGeomTool *geom = AnitaGeomTool::Instance();


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  for (int entry=startEntry; entry<stopEntry; entry++) {
    int printEntry = entry-startEntry;
    if (printEntry%10000 == 0 && printEntry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(printEntry)/totalTimeSec;
      double remaining = (float(stopEntry-printEntry)/rate)/60.;
      watch.Start();
      cout << printEntry << "/" << lenEntries << " | ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;
    }
    summaryTree->GetEntry(entry);

    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(gpsEvent);

    double peakPhi = evSum->peak[0][0].phi;
    
    int nant = ap->getClosestAntennas(peakPhi, antennasInSum, closest, 0, AnitaPol::kHorizontal);
    if (nant != antennasInSum) {
      cout << "Only using " << nant << " antennas to compute noise" << endl;
    }

    int phiSectors[phiSectorsInSum];
    int phiSector = 0;
    peakPhiSector = -1;
    if (debug) cout << "phi=" << peakPhi << " | phi sectors in use: ";
    for (int ant=0; ant<antennasInSum; ant++) {
      if (closest[ant] < 16) { //0->15 is top ring
	phiSectors[phiSector] = geom->getPhiFromAnt(closest[ant]);
	if (debug) cout << phiSectors[phiSector] << " ";
	if (peakPhiSector == -1) peakPhiSector = phiSectors[phiSector];
	phiSector++;
      }
    }
    if (debug) cout << endl;

    double noiseAvg = 0;
    for (phiSector=0; phiSector<phiSectorsInSum; phiSector++) {
      for (int ring=0; ring<3; ring++) {
	double currNoise = noiseSum->avgRMSNoise[phiSectors[phiSector]][ring][0];
	if (debug) cout << phiSectors[phiSector] << " " << ring << " " << currNoise << endl;
	noiseAvg += currNoise;
      }
    }
    if (debug) cout << "rms:" << rms << " " << noiseSum->fifoLength << endl;
    noiseAvg /= noiseSum->fifoLength;
    noiseAvg /= antennasInSum;

    rms = noiseAvg;
    eventNumber = evSum->eventNumber;
    snr_filtered = evSum->coherent_filtered[0][0].peakVal/noiseAvg;
    snr = evSum->coherent[0][0].peakVal/noiseAvg;
    if (debug) cout << "evNum:" << eventNumber << " noiseAvg:" << noiseAvg << endl;
    if (debug) cout << "peakVal:" << evSum->coherent[0][0].peakVal << " snr:" << snr << endl;
    if (debug) cout << "peakVal_f:" << evSum->coherent_filtered[0][0].peakVal << " snr_f:" << snr << endl;

    outFile->cd();
    snrTree->Fill();

    delete usefulGPS;

  }

  cout << "Done looping!  Just saving now" << endl;
  
  outFile->cd();
  snrTree->Write();
  outFile->Close();
  
  cout << "Done completely :)  Bye!" << endl;


  return;
}


void combineSNRFiles(int numSplits) {
  
  TChain *inTree = new TChain("snrTree","snrTree");

  stringstream name;
  for (int i=0; i<numSplits; i++) {
    name.str("");
    name << "data/SNRs_" << i << ".root";
    inTree->Add(name.str().c_str());
  }

  TFile *outFile = TFile::Open("snrFile.root","recreate");
  inTree->CloneTree(-1,"fast");
  outFile->Write();
  outFile->Close();

  return;
}
  
  


void makeSNRFile(int numSplits, int split) {
  makeSNRFileSplit(numSplits,split);
  return;
}


void makeSNRFile() {
  cout << "loaded makeSNRFile.C" << endl;
  return;
}
