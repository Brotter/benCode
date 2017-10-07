/*

  Whoops!  I need to pre-compile this list before running over things because the wais events aren't in order

 */

#include "AnitaConventions.h"
#include "AntennaPositions.h"
#include "UsefulAdu5Pat.h"

#include "loadAll.C"

void makeSNRFile(string date="09.27.17_19h") {
  bool debug=false;
  /*

    This should go through the AnitaNoiseSummaries, and in conjunction with the peak pointing direction and
    peak of the waveform (I guess?), determine the SNR of that waveform.

    Problem: I don't save the minimum!  I can't do peak to peak DAMNIT.  I'll have to fix that

   */


  TChain *summaryTree = loadAll(date,false);
  int lenEntries = summaryTree->GetEntries();
  AnitaEventSummary *evSum = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  Adu5Pat *gpsEvent = NULL;
  summaryTree->SetBranchAddress("gpsEvent",&gpsEvent);
  AnitaNoiseSummary *noiseSum = NULL;
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);


  char* dataDir = getenv("ANITA3_RESULTSDIR");
  stringstream name;
  name.str("");
  name << dataDir << "/templateSearch/" << date << "/SNRs.root";
  cout << "using " << name.str() << " as output file" << endl;
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  TTree *snrTree = new TTree("snrTree","snrTree");
  int eventNumber;
  snrTree->Branch("eventNumber",&eventNumber,"eventNumber/I");
  double snr,snr_filtered;
  snrTree->Branch("snr",&snr,"snr/D");
  snrTree->Branch("snr_filtered",&snr_filtered,"snr_filtered/D");

  const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 
  const int antennasInSum = 15; //I'm using 15 antennas
  const int phiSectorsInSum = antennasInSum/3;
  int closest[antennasInSum];

  AnitaGeomTool *geom = AnitaGeomTool::Instance();


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0 && entry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entry)/totalTimeSec;
      double remaining = (float(lenEntries-entry)/rate)/60.;
      watch.Start();
      cout << entry << "/" << lenEntries << " | ";
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
    if (debug) cout << "phi=" << peakPhi << " | phi sectors in use: ";
    for (int ant=0; ant<antennasInSum; ant++) {
      if (closest[ant] < 16) { //0->16 is top
	phiSectors[phiSector] = geom->getPhiFromAnt(closest[ant]);
	if (debug) cout << phiSectors[phiSector] << " ";
	phiSector++;
      }
    }
    if (debug) cout << endl;

    double noiseAvg = 0;
    for (phiSector=0; phiSector<phiSectorsInSum; phiSector++) {
      for (int ring=0; ring<3; ring++) {
	noiseAvg += (noiseSum->avgRMSNoise[phiSectors[phiSector]][ring][0]/noiseSum->fifoLength)/antennasInSum;
      }
    }

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
