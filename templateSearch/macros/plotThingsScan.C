void makeReducedHistograms(TH1D **histograms,string subName) {

  //PointingHypothesis
  //1
  histograms[0] = new TH1D("peakPhi"+TString(subName),"Map Peak Phi;Map Peak Phi;Count",360,0,360);

  //2
  histograms[1] = new TH1D("peakTheta"+TString(subName),"Map Peak Theta;Map Peak Theta;Count",141,-60,80);

  //3
  histograms[2] = new TH1D("mapPeak"+TString(subName),"Map Peak Value;Map Peak Value;Count",600,0,0.6);

  //4
  histograms[3] = new TH1D("mapSNR"+TString(subName),"Map SNR; Map SNR;; Count",1000,0,100);
  
  //5
  histograms[4] = new TH1D("peakRatio"+TString(subName),"Peak Ratio - p2/p1; Peak Ratio (p2/p1); Count",300,0,30);

  //WaveformInfo (Unfiltered)
  //6
  histograms[5] = new TH1D("linPolFrac"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //7
  histograms[6] = new TH1D("linPolAng"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",181,-45,45);
  
  //8
  histograms[7] = new TH1D("waveSNR"+TString(subName),"Coherently Summed Waveform SNR;SNR;Count",1000,0,300);
  
  //9
  histograms[8] = new TH1D("wavePeakVal"+TString(subName),"Coherently Summed Waveform Peak Value; Peak (mV);Count",1000,0,1000);
  
  //10
  histograms[9] = new TH1D("wavePeakHilb"+TString(subName),"Coherently Summed Waveform Peak Hilbert;Peak Hilbert;Count",2000,0,2000);
  
  //11
  histograms[10] = new TH1D("impulsivity"+TString(subName),"Impulsivity Measurement;Impulsivity;Count",500,0,1);

  //WaveformInfo (Filtered)
  //12
  histograms[11] = new TH1D("linPolFrac_F"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //13
  histograms[12] = new TH1D("linPolAng_F"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",181,-45,45);
  
  //14
  histograms[13] = new TH1D("waveSNR_F"+TString(subName),"Coherently Summed Waveform SNR;SNR;Count",1000,0,300);
  
  //15
  histograms[14] = new TH1D("wavePeakVal_F"+TString(subName),"Coherent_Filteredly Summed Waveform Peak Value; Peak (mV);Count",1000,0,1000);
  
  //16
  histograms[15] = new TH1D("wavePeakHilb_F"+TString(subName),"Coherently Summed Filtered Waveform Peak Hilbert;Peak Hilbert;Count",2000,0,2000);
  
  //17
  histograms[16] = new TH1D("impulsivity_F"+TString(subName),"Impulsivity Measurement (Filtered);Impulsivity;Count",500,0,1);

  //WaveformInfo (Deconvolved)
  //18
  histograms[17] = new TH1D("linPolFrac_D"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //19
  histograms[18] = new TH1D("linPolAng_D"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",181,-45,45);

  //20
  histograms[19] = new TH1D("waveSNR_D"+TString(subName),"Deconvolved Waveform SNR;SNR;Count",1000,0,300);
  
  //21
  histograms[20] = new TH1D("wavePeakVal_D"+TString(subName),"Deconvolved Waveform Peak Value; Peak (mV);Count",1000,0,1000);

  //22
  histograms[21] = new TH1D("wavePeakHilb_D"+TString(subName),"Deconvolved Waveform Peak Hilbert;Peak Hilbert;Count",2000,0,2000);

  //23
  histograms[22] = new TH1D("impulsivity_D"+TString(subName),"Impulsivity Measurement;Impulsivity;Count",500,0,1);


  //WaveformInfo (Deconvolved Filtered)
  //18
  histograms[23] = new TH1D("linPolFrac_DF"+TString(subName),"Linear Polarization Fraction; Linear Polarization Fraction; Count",1000,0,1);
  
  //19
  histograms[24] = new TH1D("linPolAng_DF"+TString(subName),"Linear Polarization Angle; Linear Polarization Angle; Count",181,-45,45);

  //20
  histograms[25] = new TH1D("waveSNR_DF"+TString(subName),"Deconvolved Waveform SNR;SNR;Count",1000,0,300);
  
  //21
  histograms[26] = new TH1D("wavePeakVal_DF"+TString(subName),"Deconvolved Waveform Peak Value; Peak (mV);Count",1000,0,1000);

  //22
  histograms[27] = new TH1D("wavePeakHilb_DF"+TString(subName),"Deconvolved Waveform Peak Hilbert;Peak Hilbert;Count",2000,0,2000);

  //23
  histograms[28] = new TH1D("impulsivity_DF"+TString(subName),"Impulsivity Measurement;Impulsivity;Count",500,0,1);


  //TemplateInfo
  //24
  histograms[29] = new TH1D("template_Wais"+TString(subName),"Wais Template Correlation Value;Correlation Value; Count",250,0,1);
  //25
  histograms[30] = new TH1D("template_cRay"+TString(subName),"CRay +4 Template Correlation Value;Correlation Value; Count",250,0,1);

  


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

    histograms[23]->Fill(eventSummary->deconvolved_filtered[0][0].linearPolFrac());
    histograms[24]->Fill(eventSummary->deconvolved_filtered[0][0].linearPolAngle());
    histograms[25]->Fill(eventSummary->deconvolved_filtered[0][0].snr);
    histograms[26]->Fill(eventSummary->deconvolved_filtered[0][0].peakVal);
    histograms[27]->Fill(eventSummary->deconvolved_filtered[0][0].peakHilbert);
    histograms[28]->Fill(eventSummary->deconvolved_filtered[0][0].impulsivityMeasure);

    histograms[29]->Fill(templateSummary->coherent[0][0].wais);
    histograms[30]->Fill(templateSummary->coherent[0][0].cRay[4]);

}

void plotThingsScan() {

  /*
    Proof is stupid and keeps dropping huge amounts of events for some reason.

    This is another way to make a bunch of histograms, though it does take awhile.

   */

  //  TChain *summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add("07.28.17_17h_decimated.root");

  const int numHists = 31;

  //four different sets of data
  // 1) things that aren't pulsers and don't point at bases
  TH1D *events[numHists];
  makeReducedHistograms(events,"_events");
  // 2) things that are wais pulsers
  TH1D *waisPulses[numHists];
  makeReducedHistograms(waisPulses,"_WAIS");
  // 3) things that are ldb pulsers
  TH1D *ldbPulses[numHists];
  makeReducedHistograms(ldbPulses,"_LDB");
  // 4) things that point at bases but aren't bases
  TH1D *basePointed[numHists];
  makeReducedHistograms(basePointed,"_basePointed");
  // 5) minbias
  TH1D *minbias[numHists];
  makeReducedHistograms(minbias,"_minbias");

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
    }
  
    //if no pulser
    if (eventSummary->flags.pulser == 0) {
      double waisDiff = FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->wais.phi));
      double ldbDiff  = FFTtools::wrap(TMath::Abs(eventSummary->peak[0][0].phi - eventSummary->ldb.phi));
      //if pointed at (and near) camps
      if ( (eventSummary->wais.distance < 1000e3 && waisDiff < 5) || 
	   (eventSummary->ldb.distance < 1000e3 && ldbDiff < 5) ) {
	fillHistograms(eventSummary,templateSummary,basePointed);
      }
      //if minbias
      else if (!eventSummary->flags.isRF) {
	fillHistograms(eventSummary,templateSummary,minbias);
      }	
      //otherwise it is thermal
      else {
	fillHistograms(eventSummary,templateSummary,events);
      }	
    }

    //if wais
    else if (eventSummary->flags.pulser == 1) {
      fillHistograms(eventSummary,templateSummary,waisPulses);
    }

    //if ldb
    else if (eventSummary->flags.pulser == 2) {
      fillHistograms(eventSummary,templateSummary,ldbPulses);
    }

  }

  cout << "Done! Writing to file" << endl;

  TFile *outFile = TFile::Open("plotThingsScan.root","recreate");
  for (int i=0; i<numHists; i++) {
    outFile->cd();
    basePointed[i]->Write();
    events[i]->Write();
    waisPulses[i]->Write();
    ldbPulses[i]->Write();
    minbias[i]->Write();
  }

  outFile->Close();

  return;
}


TH1* makeNormCumulative(TH1* inHist) {
  /*
    Makes a "normalized cumulative cut fraction" graph I think

    Should let me set cuts on a specific amount of "reduction"
   */


  TH1* copyHist = (TH1*)inHist->Clone();
  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);
  TH1* outHist = (TH1*)copyHist->GetCumulative();

  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = outHist->GetBinContent(i);
    outHist->SetBinContent(i,1.-value);
  }


  delete copyHist;

  outHist->GetYaxis()->SetTitle("Surviving Fraction");
  return outHist;
}


TH1* makeNormHist(TH1* inHist) {
  /*

    Make a normalized histogram

  */

  TH1* copyHist = (TH1*)inHist->Clone();
  copyHist->GetYaxis()->SetTitle("Occupancy");

  
  double integral = 0;
  for (int i=0; i<copyHist->GetNbinsX(); i++) {
    double value = copyHist->GetBinContent(i);
    integral += value;
  }

  copyHist->Scale(1./integral);

  return copyHist;
}

void drawOneDHistos(string inFileName="plotThingsScan.root") {
  
  stringstream name;

  TFile *inFile = TFile::Open(inFileName.c_str());

  const int numLines = 5;
  string subStrings[numLines] = {"thermal","WAIS","LDB","basePointed","minbias"};

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);


  const int numGraphs = 31;
  string graphNames[numGraphs] = {"peakPhi","peakTheta","mapPeak","mapSNR","peakRatio",
			   "linPolFrac","linPolAng","waveSNR","wavePeakVal","wavePeakHilb","impulsivity",
			   "linPolFrac_F","linPolAng_F","waveSNR_F","wavePeakVal_F","wavePeakHilb_F","impulsivity_F",
			   "linPolFrac_D","linPolAng_D","waveSNR_D","wavePeakVal_D","wavePeakHilb_D","impulsivity_D",
			   "linPolFrac_DF","linPolAng_DF","waveSNR_DF","wavePeakVal_DF","wavePeakHilb_DF","impulsivity_DF",
			   "template_Wais","template_cRay"};

  for (int graph=0; graph<numGraphs; graph++) {
    c1->Clear();
    c1->SetGrid(1);
    c1->SetLogy();
    TLegend *leg = new TLegend(0.8,0.7,0.95,0.95);
    for (int i=0; i<numLines; i++) {
      name.str("");
      name << graphNames[graph] << "_" << subStrings[i];
      TH1 *hist = makeNormCumulative((TH1D*)inFile->Get(name.str().c_str()));
      hist->SetLineColor(i+1);
      hist->GetYaxis()->SetRangeUser(1e-7,1);
      hist->SetStats(0);
      if (i==0) hist->Draw();
      else      hist->Draw("same");
      if (strcmp(subStrings[i].c_str(),"thermal")==0) leg->AddEntry(hist,"events","l");
      else leg->AddEntry(hist,subStrings[i].c_str(),"l");
    }

    leg->Draw();
    c1->SaveAs(TString("plots/")+graphNames[graph]+TString("_cumulative.png"));
  }


  for (int graph=0; graph<numGraphs; graph++) {
    c1->Clear();
    c1->SetGrid(1);
    c1->SetLogy();
    TLegend *leg = new TLegend(0.8,0.7,0.95,0.95);
    for (int i=0; i<numLines; i++) {
      name.str("");
      name << graphNames[graph] << "_" << subStrings[i];
      TH1 *hist = makeNormHist((TH1D*)inFile->Get(name.str().c_str()));
      hist->SetLineColor(i+1);
      hist->GetYaxis()->SetRangeUser(1e-8,1);
      hist->SetStats(0);
      if (i==0) hist->Draw();
      else      hist->Draw("same");
      if (strcmp(subStrings[i].c_str(),"thermal")==0) leg->AddEntry(hist,"events","l");
      else leg->AddEntry(hist,subStrings[i].c_str(),"l");
    }

    leg->Draw();
    c1->SaveAs(TString("plots/")+graphNames[graph]+TString("_hist.png"));
  }

  return;
}
  

TH1 *dividedCumulative(TH1D *h1, TH1D *h2) {

  TH1* h1Cum = makeNormCumulative(h1);
  TH1* h2Cum = makeNormCumulative(h2);

  h1Cum->Divide(h2Cum);
  delete h2Cum;
  
  return h1Cum;
}
