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

#include "loadAll.C"
#include "AnitaEventSummary.h"
#include "AnitaTemplates.h"

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
  /*

    

   */



  //wais stuff
  //  TFile *inFile = TFile::Open("waisEvents.root");
  //  TTree *summaryTree = (TTree*)inFile->Get("waisSummary");

  //decimated stuff
  //  TFile *inFile = TFile::Open("07.28.17_17h_decimated.root");
  //  TTree *summaryTree = (TTree*)inFile->Get("summaryTree");


  //full data set
  TChain *summaryTree = loadAllDefault();

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

  /* for the ABCD background estimate, grab the events that fall between some lower level and (3) */
  if (strength==4) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.5");
    cuts.push_back("template.coherent[0][0].cRay[4] < 0.67");
    cuts.push_back(waveformString+"peakHilbert > 43.5");
    cuts.push_back(waveformString+"linearPolFrac() > 0.3055");
    cuts.push_back("peak[0][0].value > 0.0435");
    cuts.push_back("peak[0][0].snr > 8.95");
  }
  /*    From the newest run with some fixed things, gets you 12037 */
  else if (strength==3) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.67");
    cuts.push_back(waveformString+"peakHilbert > 43.5");
    cuts.push_back(waveformString+"linearPolFrac() > 0.3055");
    cuts.push_back("peak[0][0].value > 0.0435");
    cuts.push_back("peak[0][0].snr > 8.95");
  }
  /*     Gets you 6121 total passing events */
  else if (strength==0) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.666");
    cuts.push_back(waveformString+"peakHilbert > 48.5");
    cuts.push_back("peak[0][0].value > 0.0435");
    cuts.push_back("peak[0][0].snr > 9.45");
  }

  /*  Generated from getCutsFromValue() @ 10^-2 cut of WAIS ("best" value?)*/
  /*     Gets you 3570 total passing events */
  else if (strength==1) {
    cuts.push_back("template.coherent[0][0].cRay[4] > 0.726");
    cuts.push_back(waveformString+"peakHilbert > 63.5");
    cuts.push_back("peak[0][0].value > 0.0605");
    cuts.push_back("peak[0][0].snr > 11.75");
  }

  /* Generated from getCutsFromValue() @ 10^-1 cut of WAIS (strongest cut) */
  /*    Gets you 313 total passing events */
  else if (strength==2) {
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
  cuts.push_back("flags.maxBottomToTopRatio[0] < 3");
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
    if (strength==3) outFileName = "makeCuts_final.txt";
    if (strength==4) outFileName = "makeCuts_background.txt";
    ((TTreePlayer*)(summaryTree->GetPlayer()))->SetScanFileName(outFileName.c_str());
    summaryTree->Scan("eventNumber",allCuts.c_str());
  } 

  return;
}




void printCutStrengths(int strength=1) {
  /*

    For the table that shows how many events each cut parameter eliminates.
    strength:
    0 = quality
    1 = impulsive
    3 = signal


   */

  TChain *summaryTree;
  if (strength==0) summaryTree = loadAllDefault_noproof();
  if (strength==1) summaryTree = loadReKey(false);
  if (strength==3) summaryTree = new TChain("summaryTree"); summaryTree->Add("cutsClust_oct14.root");
  int lenEntries = summaryTree->GetEntries();
  if (!lenEntries) {
    cout << "didn't find any entries" << endl;
    return;
  }
  else cout << "found " << lenEntries << " entries" << endl;
  AnitaEventSummary *evSum=NULL;
  AnitaTemplateSummary *tempSum=NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);


  int evsLeft = lenEntries;

  //strength cuts
  int CmapValue = 0;
  int CmapSNR = 0;
  int Chilbert = 0;
  int CtemplateCorr = 0;
  int CtemplateDCorr = 0;
  int ClinearPolFrac = 0;

  int mapValue = 0;
  int mapSNR = 0;
  int hilbert = 0;
  int templateCorr = 0;
  int templateDCorr = 0;
  int linearPolFrac = 0;

  //quality cuts
  int pulserLDB = 0;
  int pulserWAIS = 0;
  int baseLDB = 0;
  int baseWAIS = 0;
  int aboveHorizontal = 0;
  int payloadBlast = 0;
  int notRF = 0;

  int CpulserLDB = 0;
  int CpulserWAIS = 0;
  int CbaseLDB = 0;
  int CbaseWAIS = 0;
  int CaboveHorizontal = 0;
  int CpayloadBlast = 0;
  int CnotRF = 0;


  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    if (!(entry%10000)) {
      cout << entry << "/" << lenEntries << " " << evsLeft << " " << baseWAIS << " " << baseLDB << " |";
      if (strength != 0) { 
	cout << mapValue << " " << mapSNR << " " << hilbert << " ";
	cout << linearPolFrac << " " << templateCorr << " " << templateDCorr << endl;
      }
      else { 
	cout << pulserLDB << " " << pulserWAIS << " "  << baseLDB << " " << baseWAIS << " ";
	cout << aboveHorizontal << " " << payloadBlast << " " << notRF << endl;
      }
    }
	
    //start by removing a remaining event.  If it makes it to the end without failing, it will be added back
    evsLeft--;


    // the "good" events still have WAIS and LDB in them it seems, and maybe hw trigger angle
    // so this needs to break the loop if you aren't doing the strength(0) quality cuts, so you don't count bad things
    if (strength != 0) {
      if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0)),2) + 
		       pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->ldb.theta,360,0)),2))  < 6 
	   && evSum->ldb.distance  < 700e3)) {baseLDB++; continue;}
      //not pointed at wais when it is nearby (~700km)
      if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)),2) + 
		       pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->wais.theta,360,0)),2))  < 6 
	   && evSum->wais.distance  < 700e3)) {baseWAIS++; continue;}
      if (evSum->peak[0][0].theta < 0) continue;
      if (evSum->flags.maxBottomToTopRatio[0] > 3) continue;
      if (!evSum->flags.isRF) continue;
      if (evSum->flags.pulser != 0) continue;
    }



    //Quality cuts
    if (strength == 0) {

      if (evSum->flags.pulser == 1) pulserWAIS++;
      if (evSum->flags.pulser == 2) pulserLDB++;
      if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0)),2) + 
		       pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->ldb.theta,360,0)),2))  < 6 
	   && evSum->ldb.distance  < 700e3)) {
	baseLDB++;}
      if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)),2) + 
			    pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->wais.theta,360,0)),2))  < 6 
		&& evSum->wais.distance  < 700e3)) {
	baseWAIS++; }
      if (evSum->peak[0][0].theta < 0) aboveHorizontal++;
      if (evSum->flags.maxBottomToTopRatio[0] > 3) payloadBlast++;
      if (!evSum->flags.isRF) notRF++;
    

      if (evSum->flags.pulser == 1) pulserWAIS++;
      else if (evSum->flags.pulser == 2) pulserLDB++;
      else if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0)),2) + 
			    pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->ldb.theta,360,0)),2))  < 6 
		&& evSum->ldb.distance  < 700e3)) {
	CbaseLDB++; }
      else if ((TMath::Sqrt(pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0)),2) + 
			  pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].theta - evSum->wais.theta,360,0)),2))  < 6 
	      && evSum->wais.distance  < 700e3)) {
	CbaseWAIS++;}
      else if (evSum->peak[0][0].theta < 0) aboveHorizontal++;
      else if (evSum->flags.maxBottomToTopRatio[0] > 3) payloadBlast++;
      else if (!evSum->flags.isRF) notRF++;
      else evsLeft++;
    }


    //weak cuts. These are for determining the anthropogenic base distribution
    else if (strength == 1 ) {
      if (evSum->peak[0][0].value < 0.0435) mapValue++;
      if (evSum->peak[0][0].snr < 8.95) mapSNR++;
      if (evSum->coherent_filtered[0][0].peakHilbert < 25) hilbert++;
      if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.5) linearPolFrac++;
      if (tempSum->coherent[0][0].cRay[4] < 0.5) templateCorr++;
      if (tempSum->deconvolved[0][0].cRay[4] < 0.5) templateDCorr++;

      if (evSum->peak[0][0].value < 0.0435) CmapValue++;
      else if (evSum->peak[0][0].snr < 8.95) CmapSNR++;
      else if (evSum->coherent_filtered[0][0].peakHilbert < 25) Chilbert++;
      else if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.5) ClinearPolFrac++;
      else if (tempSum->coherent[0][0].cRay[4] < 0.5) CtemplateCorr++;
      else if (tempSum->deconvolved[0][0].cRay[4] < 0.5) CtemplateDCorr++;
      else evsLeft++;
    }


    //signal cuts
    else if (strength==3) {
      if (evSum->peak[0][0].value < 0.0435) mapValue++;
      if (evSum->peak[0][0].snr < 9.05) mapSNR++;
      if (evSum->coherent_filtered[0][0].peakHilbert < 31.1) hilbert++;
      if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.6) linearPolFrac++;
      if (tempSum->coherent[0][0].cRay[4] < 0.666) templateCorr++;
      if (tempSum->deconvolved[0][0].cRay[4] < 0.666) templateDCorr++;

      if (evSum->peak[0][0].value < 0.0435) CmapValue++;
      else if (evSum->peak[0][0].snr < 9.05) CmapSNR++;
      else if (evSum->coherent_filtered[0][0].peakHilbert < 31.1) Chilbert++;
      else if (evSum->coherent_filtered[0][0].linearPolFrac() < 0.6) ClinearPolFrac++;
      else if (tempSum->coherent[0][0].cRay[4] < 0.666) CtemplateCorr++;
      else if (tempSum->deconvolved[0][0].cRay[4] < 0.666) CtemplateDCorr++;
      else evsLeft++;
    }
  }


  if (strength == 0) {
    cout << pulserLDB << " " << CpulserLDB << endl;
    cout << pulserWAIS << " " <<  CpulserWAIS << endl;
    cout << baseLDB << " " << CbaseLDB << endl;
    cout << baseWAIS << " " << CbaseWAIS << endl; 
    cout << aboveHorizontal << " " << CaboveHorizontal << endl;
    cout << payloadBlast << " " << CpayloadBlast << endl;
    cout << notRF << " " << CnotRF << endl;
  }

    else{
    cout << mapValue << " " << CmapValue << endl;
    cout << mapSNR << " " << CmapSNR << endl;
    cout << hilbert << " " << Chilbert << endl;
    cout << linearPolFrac << " " << ClinearPolFrac << endl;
    cout << templateCorr << " " << CtemplateCorr << endl;
    cout << templateDCorr << " " << CtemplateDCorr << endl;
    cout << evsLeft << evsLeft << endl;
  }


  return;

}


void makeCuts() {
  cout << "loaded makeCuts.C" << endl;
  return;
}




