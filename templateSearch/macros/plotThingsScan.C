void makeReducedHistograms(TH1D **histograms,string subName) {

  //PointingHypothesis
  //1
  histograms[0] = new TH1D("peakPhi"+TString(subName),"Map Peak Phi;Map Peak Phi;Count",360,0,360);

  //2
  histograms[1] = new TH1D("peakTheta"+TString(subName),"Map Peak Theta;Map Peak Theta;Count",181,-90,90);

  //3
  histograms[2] = new TH1D("mapPeak"+TString(subName),"Map Peak Value;Map Peak Value;Count",1000,0,1);

  //4
  histograms[3] = new TH1D("mapSNR"+TString(subName),"Map SNR; Map SNR;; Count",180,-90,90);
  
  //5
  histograms[4] = new TH1D("peakRatio"+TString(subName),"Peak Ratio - p2/p1; Peak Ratio (p2/p1); Count",1000,0,100);

  //WaveformInfo (Unfiltered)
  //6
  histograms[5] = new TH1D("linPolFrac"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //7
  histograms[6] = new TH1D("linPolAng"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",180,-90,90);
  
  //8
  histograms[7] = new TH1D("waveSNR"+TString(subName),"Coherently Summed Waveform SNR;SNR;Count",1000,0,100);
  
  //9
  histograms[8] = new TH1D("wavePeakVal"+TString(subName),"Coherently Summed Waveform Peak Value; Peak (mV);Count",1000,0,100);
  
  //10
  histograms[9] = new TH1D("wavePeakHilb"+TString(subName),"Coherently Summed Waveform Peak Hilbert;Peak Hilbert;Count",1000,0,1000);
  
  //11
  histograms[10] = new TH1D("impulsivity"+TString(subName),"Impulsivity Measurement;Impulsivity;Count",500,0,5);

  //WaveformInfo (Filtered)
  //12
  histograms[11] = new TH1D("linPolFrac_F"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //13
  histograms[12] = new TH1D("linPolAng_F"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",180,-90,90);
  
  //14
  histograms[13] = new TH1D("waveSNR_F"+TString(subName),"Coherently Summed Waveform SNR;SNR;Count",1000,0,100);
  
  //15
  histograms[14] = new TH1D("wavePeakVal_F"+TString(subName),"Coherent_Filteredly Summed Waveform Peak Value; Peak (mV);Count",1000,0,100);
  
  //16
  histograms[15] = new TH1D("wavePeakHilb_F"+TString(subName),"Coherently Summed Filtered Waveform Peak Hilbert;Peak Hilbert;Count",1000,0,1000);
  
  //17
  histograms[16] = new TH1D("impulsivity_F"+TString(subName),"Impulsivity Measurement (Filtered);Impulsivity;Count",500,0,5);

  //WaveformInfo (Deconvolved)
  //18
  histograms[17] = new TH1D("linPolFrac_D"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //19
  histograms[18] = new TH1D("linPolAng_D"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",180,-90,90);

  //20
  histograms[19] = new TH1D("waveSNR_D"+TString(subName),"Deconvolved Waveform SNR;SNR;Count",1000,0,100);
  
  //21
  histograms[20] = new TH1D("wavePeakVal_D"+TString(subName),"Deconvolved Waveform Peak Value; Peak (mV);Count",1000,0,100);

  //22
  histograms[21] = new TH1D("wavePeakHilb_D"+TString(subName),"Deconvolved Waveform Peak Hilbert;Peak Hilbert;Count",1000,0,1000);

  //23
  histograms[22] = new TH1D("impulsivity_D"+TString(subName),"Impulsivity Measurement;Impulsivity;Count",500,0,5);

  //TemplateInfo
  //24
  histograms[23] = new TH1D("template_Wais"+TString(subName),"Wais Template Correlation Value;Correlation Value; Count",250,0,1);
  //25
  histograms[24] = new TH1D("template_cRay"+TString(subName),"CRay +4 Template Correlation Value;Correlation Value; Count",250,0,1);
  


return;
}


void fillHistograms(AnitaEventSummary *eventSummary, AnitaTemplateSummary *templateSummary,
		    TH1D** histograms) {

    histograms[0]->Fill(eventSummary->peak[0][0].phi);
    histograms[1]->Fill(eventSummary->peak[0][0].theta);
    histograms[2]->Fill(eventSummary->peak[0][0].value);
    histograms[3]->Fill(eventSummary->peak[0][0].snr);
    histograms[4]->Fill(eventSummary->peak[0][1].value/eventSummary->peak[0][0].value);
    
    histograms[5]->Fill(eventSummary->coherent[0][0].linearPolFrac());
    histograms[6]->Fill(eventSummary->coherent[0][0].linearPolAngle());
    histograms[7]->Fill(eventSummary->coherent[0][0].snr);
    histograms[8]->Fill(eventSummary->coherent[0][0].peakVal);
    histograms[9]->Fill(eventSummary->coherent[0][0].peakHilbert);
    histograms[10]->Fill(eventSummary->coherent[0][0].impulsivityMeasure);

    histograms[11]->Fill(eventSummary->coherent_filtered[0][0].linearPolFrac());
    histograms[12]->Fill(eventSummary->coherent_filtered[0][0].linearPolAngle());
    histograms[13]->Fill(eventSummary->coherent_filtered[0][0].snr);
    histograms[14]->Fill(eventSummary->coherent_filtered[0][0].peakVal);
    histograms[15]->Fill(eventSummary->coherent_filtered[0][0].peakHilbert);
    histograms[16]->Fill(eventSummary->coherent_filtered[0][0].impulsivityMeasure);

    histograms[17]->Fill(eventSummary->deconvolved[0][0].linearPolFrac());
    histograms[18]->Fill(eventSummary->deconvolved[0][0].linearPolAngle());
    histograms[19]->Fill(eventSummary->deconvolved[0][0].snr);
    histograms[20]->Fill(eventSummary->deconvolved[0][0].peakVal);
    histograms[21]->Fill(eventSummary->deconvolved[0][0].peakHilbert);
    histograms[22]->Fill(eventSummary->deconvolved[0][0].impulsivityMeasure);

    histograms[23]->Fill(templateSummary->coherent[0][0].wais);
    histograms[24]->Fill(templateSummary->coherent[0][0].cRay[4]);

}

void plotThingsScan() {

  /*
    Proof is stupid and keeps dropping huge amounts of events for some reason.

    This is another way to make a bunch of histograms, though it does take awhile.

   */

  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  const int numHists = 25;

  //four different sets of data
  // 1) things that aren't pulsers and don't point at bases
  TH1D *thermal[numHists];
  makeReducedHistograms(thermal,"_thermal");
  // 2) things that are wais pulsers
  TH1D *waisPulses[numHists];
  makeReducedHistograms(waisPulses,"_WAIS");
  // 3) things that are ldb pulsers
  TH1D *ldbPulses[numHists];
  makeReducedHistograms(ldbPulses,"_LDB");
  // 4) things that point at bases but aren't bases
  TH1D *basePointed[numHists];
  makeReducedHistograms(basePointed,"_basePointed");


  AnitaEventSummary *eventSummary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&eventSummary);

  AnitaTemplateSummary *templateSummary = NULL;
  summaryTree->SetBranchAddress("template",&templateSummary);

  int lenEntries = summaryTree->GetEntries();

  cout << "found " << lenEntries << " entries" << endl;


  TStopwatch watch;
  watch.Start(true);
  int totalTime = 0;

  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);

    int entriesRemaining = lenEntries-entry;
    if (entry%10000 == 0 && entry!=0) {
      totalTime += watch.RealTime();
      double rate = entry/totalTime;
      watch.Start();
      cout << entry << "/" << lenEntries << " (" << rate << " ev/sec) ";
      cout << (entriesRemaining/rate)/(3600.) << " hours remain" << endl;

      if (eventSummary->flags.pulser == 0) {
	double waisDiff = FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->wais.phi));
	double ldbDiff  = FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->ldb.phi));
	if ( (eventSummary->wais.distance < 1000e3 && waisDiff < 5) || 
	     (eventSummary->ldb.distance < 1000e3 && ldbDiff < 5) ) {
	  fillHistograms(eventSummary,templateSummary,basePointed);
	}
	else {
	  fillHistograms(eventSummary,templateSummary,thermal);
	}
	
      }
      else if (eventSummary->flags.pulser == 1) {
	fillHistograms(eventSummary,templateSummary,waisPulses);
      }
      else if (eventSummary->flags.pulser == 2) {
	fillHistograms(eventSummary,templateSummary,ldbPulses);
      }


    }

  }
  cout << "Done! Writing to file" << endl;

  TFile *outFile = TFile::Open("pointThingsScan.root","recreate");
  for (int i=0; i<numHists; i++) {
    outFile->cd();
    basePointed[i]->Write();
    thermal[i]->Write();
    waisPulses[i]->Write();
    ldbPulses[i]->Write();
  }

  outFile->Close();

  return;
  }
