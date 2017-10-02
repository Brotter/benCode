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


//I need to make this a library or something
TH1* makeNormCumulative(TH1* inHist,Bool_t forward=true) {
  /*
    Makes a "normalized cumulative cut fraction" graph I think

    Should let me set cuts on a specific amount of "reduction"

    forward specifies the direction of the cumulative sum
   */


  TH1* copyHist = (TH1*)inHist->Clone();
  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);
  TH1* outHist = (TH1*)copyHist->GetCumulative(forward);

  for (int i=0; i<copyHist->GetNbinsX()+1; i++) {
    double value = outHist->GetBinContent(i);
    outHist->SetBinContent(i,1.-value);

  }
  delete copyHist;

  if (forward) outHist->GetYaxis()->SetTitle("Surviving Fraction");
  else         outHist->GetYaxis()->SetTitle("Fraction of Events Cut");
  return outHist;
}


double findCrossing(TH1* hist,double value) {
  double outValue;
  for (int bin=0; bin<hist->GetNbinsX(); bin++) {
    if (hist->GetBinContent(bin+1) < value) {
      outValue = hist->GetBinCenter(bin+1);
      break;
    }
  }
  
  return outValue;
}
  

string getThermalCuts(double reduction=1e-3) {
  /*
    I need a way to get the cut values from thermal events.  This should do it!

    uses the histograms from makeMinbiasBackgroundHist() in cluster.C

   */


  TFile *inFile = TFile::Open("minbiasBackgrounds.root");
  if (!inFile->IsOpen()) {
    cout << "Couldn't find minbiasBackgrounds.root" << endl;
    return "";
  }
  
  TH1* hMapPeak = makeNormCumulative((TH1D*)inFile->Get("histMapPeak"));
  double mapPeakVal = findCrossing(hMapPeak,reduction);

  TH1* hMapSNR = makeNormCumulative((TH1D*)inFile->Get("histMapSNR"));
  double mapPeakSNR = findCrossing(hMapSNR,reduction);

  TH1* hTemplate = makeNormCumulative((TH1D*)inFile->Get("histTemplate"));
  double templateVal = findCrossing(hTemplate,reduction);

  TH1* hHilbert = makeNormCumulative((TH1D*)inFile->Get("histHilbert"));
  double hilbertPeak = findCrossing(hHilbert,reduction);  
  
  stringstream name;
  name << "peak[0][0].value > " << mapPeakVal << " && ";
  name << "peak[0][0].snr > " << mapPeakSNR << " && ";
  name << "template.coherent[0][0].cRay[4] > " << templateVal << " && ";
  name << "deconvolved_filtered[0][0].peakHilbert > " << hilbertPeak;

  return name.str();

   
  
}

void printPassingEvents(bool draw=true,int strength=0) {




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

  /*  Generated from drawThingsScan.C::getCutsFromValue() @ 10^-3 cut of WAIS(default) (weakest)*/

  /*    From the newest run with some fixed things, gets you 6031 */
  if (strength==3) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.67");
    cuts.push_back(waveformString+"peakHilbert > 43.5");
    cuts.push_back(waveformString+"linearPolFrac() > 0.3055");
    cuts.push_back("peak[0][0].value > 0.0435");
    cuts.push_back("peak[0][0].snr > 8.95");
  }
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
  else {
    cout << "I don't know what strength setting that is so I refuse to do anything till you get your shit together" << endl;
    return;
  }

  /* Cuts that should _always_ be made for candidates*/
  // not flagged as a pulser
  cuts.push_back("flags.pulser == 0");
  // needs an rf trigger
  cuts.push_back("flags.isRF");
  //not pointed at ldb when it is nearby (~700km)
  cuts.push_back("(TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(peak[0][0].phi - ldb.phi,360,0)),2) + pow(TMath::Abs(FFTtools::wrap(peak[0][0].theta - ldb.theta,360,0)),2))  > 4 || ldb.distance  > 700e3)");
  //not pointed at wais when it is nearby (~700km)
  cuts.push_back("(TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(peak[0][0].phi - wais.phi,360,0)),2)+pow(TMath::Abs(FFTtools::wrap(peak[0][0].theta - wais.theta,360,0)),2))  > 4 || wais.distance  > 700e3)");
  //not a blast event
  cuts.push_back("flags.maxBottomToTopRatio > 3");
  //not pointing "above" zero (- == up), since not even direct CRs will be above zero (no atmosphere!)
  cuts.push_back("peak[0][0].theta > 0");



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
    if (strength==3) outFileName = "makeCuts_final.txt");
    ((TTreePlayer*)(summaryTree->GetPlayer()))->SetScanFileName(outFileName.c_str());
    summaryTree->Scan("eventNumber",allCuts.c_str());
  } 

  return;
}






void makeCuts() {
  cout << "loaded makeCuts.C" << endl;
  return;
}
