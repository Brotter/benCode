/*


  Just make the cuts here I guess so I can save it

  draw:
       true = draw the eventNumber distribution and return the number of events passing
       false = _scan_ the eventNumbers into a file makeCuts.txt

  strength:
       -1 = Very Weak, 
       0 = Weak, 0.1% of WAIS pulses are cut
       1 = Nominal, 1% of WAIS pulses are cut
       2 = Strong, 10% of WAIS pulses are cut

*/


void makeCuts(bool draw=true,int strength=0) {




  //wais stuff
  //  TFile *inFile = TFile::Open("waisEvents.root");
  //  TTree *summaryTree = (TTree*)inFile->Get("waisSummary");

  //decimated stuff
  //  TFile *inFile = TFile::Open("07.28.17_17h_decimated.root");
  //  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");


  //full data set
  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  if (summaryTree == NULL) {
    cout << "Couldn't find tree, quitting" << endl;
    return -1;
  }
  else {
    cout << "found tree with " << summaryTree->GetEntries() << " entries" << endl;
  }


  vector<string> cuts;

  string waveformString = "deconvolved_filtered[0][0].";


  /* by hand picked out values when wais pulses have 99% efficiency  */
  /*
  cuts.push_back(waveformString+"impulsivityMeasure > 0.62");
  cuts.push_back(waveformString+"peakHilbert > 62");
  cuts.push_back(waveformString+"peakVal > 39");
  cuts.push_back(waveformString+"snr > 8.25");
  cuts.push_back(waveformString+"linearPolFrac() > 0.606");
  cuts.push_back("peak[0][0].value > 0.06");
  cuts.push_back("peak[0][0].snr > 16");
  cuts.push_back("(peak[0][1].value/peak[0][0].value) > 0.1");
  cuts.push_back("template.coherent[0][0].wais > 0.72");
  cuts.push_back("TMath::Abs("+waveformString+"linearPolAngle()) < 15"); 
  */

  /*  Generated from getCutsFromValue() @ 10^-3 cut of WAIS(default) (weakest)*/
  /*     Gets you 6121 total passing events */
  if (strength==0) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.666");
    cuts.push_back(waveformString+"peakHilbert > 48.5");
    cuts.push_back("peak[0][0].value > 0.0435");
    cuts.push_back("peak[0][0].snr > 9.45");
  }

  /*  Generated from getCutsFromValue() @ 10^-2 cut of WAIS ("best" value?)*/
  /*     Gets you 3570 total passing events */
  if (strength==1) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.726");
    cuts.push_back(waveformString+"peakHilbert > 63.5");
    cuts.push_back("peak[0][0].value > 0.0605");
    cuts.push_back("peak[0][0].snr > 11.75");
  }

  /* Generated from getCutsFromValue() @ 10^-1 cut of WAIS (strongest cut) */
  /*    Gets you 313 total passing events */
  if (strength==2) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.794");
    cuts.push_back(waveformString+"peakHilbert > 81.5");
    cuts.push_back("peak[0][0].value > 0.0945");
    cuts.push_back("peak[0][0].snr > 16.15");
  }


  /* Cuts that should _always_ be made */
  cuts.push_back("flags.pulser == 0");
  cuts.push_back("flags.isRF");
  cuts.push_back("(TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(peak[0][0].phi - ldb.phi,360,0)),2) + pow(TMath::Abs(FFTtools::wrap(peak[0][0].theta - ldb.theta,360,0)),2))  > 4 || ldb.distance  > 1000e3)");
  cuts.push_back("(TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(peak[0][0].phi - wais.phi,360,0)),2)+pow(TMath::Abs(FFTtools::wrap(peak[0][0].theta - wais.theta,360,0)),2))  > 4 || wais.distance  > 1000e3)");



  string allCuts = "";
  for (int i=0; i<cuts.size(); i++) {
    if (i>0) allCuts += " && ";
    allCuts += cuts[i];
  }
  cout << allCuts << endl;

  /*  Event number start and end */
  const int firstEv = 5428201;
  const int lastEv = 84415499;

  TH1D *hEv;
  if (draw) {
    cout << "Drawing" << endl;
    hEv = new TH1D("hEv","Event Numbers passing cuts",10000,firstEv,lastEv);
    summaryTree->Draw("eventNumber >> hEv",allCuts.c_str());
    double numSurviving = hEv->GetEntries();
  }
  else {
    cout << "Scanning" << endl;
    ((TTreePlayer*)(summaryTree->GetPlayer()))->SetScanRedirect(true);
    string outFileName;
    if (strength==0) outFileName = "makeCuts_weak.txt";
    if (strength==1) outFileName = "makeCuts_nominal.txt";
    if (strength==2) outFileName = "makeCuts_strong.txt";
    ((TTreePlayer*)(summaryTree->GetPlayer()))->SetScanFileName(outFileName.c_str());
    summaryTree->Scan("eventNumber",allCuts.c_str());
  } 

  return;
}
