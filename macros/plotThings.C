
TH1* makeNormCumulative(TH1* inHist) {

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

  return outHist;
}

/*==========
  I like gifs*/
void makeMovies(TChain *summaryTree) {

  //remember that gifs are stupid so make a bunch of pngs to stick together with ffmpeg
  
  cout << "making movies!" << endl;

  const int numFrames = 1000;
  const int evsPerFrame = summaryTree->GetEntries()/numFrames;

  const int startEvNum = 45671358;

  cout << "there are going to be " << evsPerFrame << "events per frame" << endl;

  stringstream name,title;
    
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  for (int frame=0; frame<numFrames; frame++ ) {
    int startEntry = frame*evsPerFrame + startEvNum;
    int endEntry = (1+frame)*evsPerFrame + startEvNum;
    name.str("");
    name << "flags.pulser == 0 && eventNumber > " << startEntry << " && eventNumber <= " <<  endEntry;
    title.str("");
    title << "Cosmic Ray Template +4 Correlation - No Pulsers - ev " << startEntry << "; Interferometric Peak; Template Corr";
    cout << startEntry << endl;

    TH2D *currHist = new TH2D("currHist",title.str().c_str(),500,0,0.2,500,0,1);
      
    c1->cd();
    summaryTree->Draw("templateResults.coherentH[0].templateCRay[5]:peak[0][0].value >> currHist",name.str().c_str(),"colz");


    cout << "There are " << currHist->GetEntries() << " entries in that histogram" << endl;


    name.str("");
    name << "movies/template_map" << frame << ".png";
    cout << name.str() << endl;
    c1->SaveAs(name.str().c_str());
    c1->Clear();

    delete currHist;
  }



  return;
}



void plotCorr(TChain *summaryTree,TFile *outFile) {
  /*================
    Template Correlation Stuff
  */
  cout << "Doing template correlation stuff" << endl;

  TH2D *noiseTmplt = new TH2D("noiseTmplt","Cosmic Ray Template +4 Correlation - No Pulsers; Interferometric Peak; Template Corr",
			      500,0,0.5,500,0,1);
  TH2D *waisTmplt = new TH2D("waisTmplt","Cosmic Ray Template +4 Correlation - WAIS; Interferometric Peak; Template Corr",
			     500,0,0.5,500,0,1);
  TH2D *ldbTmplt = new TH2D("ldbTmplt","Cosmic Ray Template +4 Correlation - LDB; Interferometric Peak; Template Corr",
			    500,0,0.5,500,0,1);
  TH2D *cutTmplt = new TH2D("cutTmplt","Cosmic Ray Template +4 Correlation - Simple Cut; Interferometric Peak; Template Corr",
			    500,0,0.5,500,0,1);

  summaryTree->Draw("templateResults.coherentH[0].templateCRay[5]:peak[0][0].value >> noiseTmplt",
		    "flags.pulser == 0","colz");
  summaryTree->Draw("templateResults.coherentH[0].templateCRay[5]:peak[0][0].value >> waisTmplt",
		    "flags.pulser == 1","colzSame");
  summaryTree->Draw("templateResults.coherentH[0].templateCRay[5]:peak[0][0].value >> ldbTmplt",
		    "flags.pulser == 2","colzSame");
  summaryTree->Draw("templateResults.coherentH[0].templateCRay[5]:peak[0][0].value >> cutTmplt",
		    "flags.pulser == 0 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05","colz");

  outFile->cd();
  noiseTmplt->Write();
  waisTmplt->Write();
  ldbTmplt->Write();
  cutTmplt->Write();

 

  TH1D *h1Tmplt = noiseTmplt->ProjectionY("h1Tmplt",0,500);
  h1Tmplt->SetTitle("Cosmic Ray Template +4 Correlation - No Pulsers;Template Corr; Occupancy");
  TH1D *h1TmpltWais = waisTmplt->ProjectionY("h1TmpltWais",0,500);
  h1TmpltWais->SetTitle("Cosmic Ray Template +4 Correlation - Wais;Template Corr; Occupancy");
  TH1D *h1TmpltLDB = ldbTmplt->ProjectionY("h1TmpltLDB",0,500);
  h1TmpltLDB->SetTitle("Cosmic Ray Template +4 Correlation - LDB;Template Corr; Occupancy");

  TH1D *h1IPeak = noiseTmplt->ProjectionX("h1IPeak",0,500);
  h1IPeak->SetTitle("Interferometric Map Peak - No Pulsers; Map Peak; Occupancy");
  TH1D *h1IPeakWais = waisTmplt->ProjectionX("h1IPeakWais",0,500);
  h1IPeakWais->SetTitle("Interferometric Map Peak - WAIS; Map Peak; Occupancy");
  TH1D *h1IPeakLDB = ldbTmplt->ProjectionX("h1IPeakLDB",0,500);
  h1IPeakLDB->SetTitle("Interferometric Map Peak - LDB; Map Peak; Occupancy");



  outFile->cd();
  h1Tmplt->Write();
  h1TmpltWais->Write();
  h1TmpltLDB->Write();
  h1IPeak->Write();
  h1IPeakWais->Write();
  h1IPeakLDB->Write();

  h1Tmplt->Scale(1./h1Tmplt->GetIntegral()[500]);
  TH1 *h1TmpltCum = h1Tmplt->GetCumulative();
  h1TmpltCum->Write();
  TH1 *h1TmpltCumWais = makeNormCumulative(h1TmpltWais);
  h1TmpltCumWais->Write();
  TH1 *h1TmpltCumLDB = makeNormCumulative(h1TmpltLDB);
  h1TmpltCumLDB->Write();


  //LDB pulses as time are actually kinda neat since we used different pulser configs
  TH2D *ldbVsTime = new TH2D("ldbVsTime","LDB Pulser Template Correlations;eventnumber; Template Correlation",
			     200,5400e3,5550e3,200,0,1);

  summaryTree->Draw("templateCRayH[5]:eventNumber >> ldbVsTime","flags.pulser==2","colz");
  outFile->cd();
  ldbVsTime->Write();

  return;
  /*
    --------------*/
}



  /* =========================
     SNR Cuts
  */
void plotSNR(TChain *summaryTree, TFile *outFile) {

  cout << "doing SNR cut stuff" << endl;


  TH1D *coherentSNR = new TH1D("coherentSNR","Coherent SNR;Coherent SNR; Occupancy",100,0,100);
  TH1D *coherentSNRWais = new TH1D("coherentSNRWais","Coherent SNR Wais;Coherent SNR; Occupancy",100,0,100);
  TH1D *coherentSNRLDB = new TH1D("coherentSNRLDB","Coherent SNR LDB;Coherent SNR; Occupancy",100,0,100);
  
  summaryTree->Draw("coherent_filtered[0][0].snr >> coherentSNR","flags.pulser == 0");
  summaryTree->Draw("coherent_filtered[0][0].snr >> coherentSNRWais","flags.pulser == 1");
  summaryTree->Draw("coherent_filtered[0][0].snr >> coherentSNRLDB","flags.pulser == 2");

  outFile->cd();
  coherentSNR->Write();
  coherentSNRWais->Write();
  coherentSNRLDB->Write();

  TH1 *coherentSNRCum = makeNormCumulative(coherentSNR);
  coherentSNRCum->Write();
  TH1 *coherentSNRCumWais = makeNormCumulative(coherentSNRWais);
  coherentSNRCum->Write();
  TH1 *coherentSNRCumLDB = makeNormCumulative(coherentSNRLDB);
  coherentSNRCumLDB->Write();



  TH1D *deconvolvedSNR = new TH1D("deconvolvedSNR","Deconvolved SNR;Deconvolved SNR; Occupancy",100,0,100);
  TH1D *deconvolvedSNRWais = new TH1D("deconvolvedSNRWais","Deconvolved SNR Wais;Deconvolved SNR; Occupancy",100,0,100);
  TH1D *deconvolvedSNRLDB = new TH1D("deconvolvedSNRLDB","Deconvolved SNR LDB;Deconvolved SNR; Occupancy",100,0,100);

  TH1 *deconvolvedSNRCum = makeNormCumulative(deconvolvedSNR);
  deconvolvedSNRCum->Write();
  TH1 *deconvolvedSNRCumWais = makeNormCumulative(deconvolvedSNRWais);
  deconvolvedSNRCum->Write();
  TH1 *deconvolvedSNRCumLDB = makeNormCumulative(deconvolvedSNRLDB);
  deconvolvedSNRCumLDB->Write();
  
  outFile->cd();
  deconvolvedSNR->Write();
  deconvolvedSNRWais->Write();
  deconvolvedSNRLDB->Write();


  TH1D *deconvFiltSNR = new TH1D("deconvFiltSNR","Deconvolved SNR;Deconvolved SNR; Occupancy",100,0,100);
  TH1D *deconvFiltSNRWais = new TH1D("deconvFiltSNRWais","Deconvolved SNR Wais;Deconvolved SNR; Occupancy",100,0,100);
  TH1D *deconvFiltSNRLDB = new TH1D("deconvFiltSNRLDB","Deconvolved SNR LDB;Deconvolved SNR; Occupancy",100,0,100);

  summaryTree->Draw("deconvolved_filtered[0][0].snr >> deconvFiltSNR","flags.pulser == 0");
  summaryTree->Draw("deconvolved_filtered[0][0].snr >> deconvFiltSNRWais","flags.pulser == 1");
  summaryTree->Draw("deconvolved_filtered[0][0].snr >> deconvFiltSNRLDB","flags.pulser == 2");

  outFile->cd();
  deconvFiltSNR->Write();
  deconvFiltSNRWais->Write();
  deconvFiltSNRLDB->Write();

  TH1 *deconvFiltSNRCum = makeNormCumulative(deconvFiltSNR);
  deconvFiltSNRCum->Write();
  TH1 *deconvFiltSNRCumWais = makeNormCumulative(deconvFiltSNRWais);
  deconvFiltSNRCum->Write();
  TH1 *deconvFiltSNRCumLDB = makeNormCumulative(deconvFiltSNRLDB);
  deconvFiltSNRCumLDB->Write();


  TH2D *bothSNR = new TH2D("bothSNR","Both SNR;Coherent SNR; Deconvolved SNR",100,0,100,100,0,100);
  TH2D *bothSNRWais = new TH2D("bothSNRWais","Coherent SNR Wais;Both SNR; Deconvolved SNR",100,0,100,100,0,100);
  TH2D *bothSNRLDB = new TH2D("bothSNRLDB","Both SNR LDB;Coherent SNR; Deconvolved SNR",100,0,100,100,0,100);
    
  summaryTree->Draw("coherent_filtered[0][0].snr:deconvolved_filtered[0][0].snr >> bothSNR","flags.pulser == 0","colz");
  summaryTree->Draw("coherent_filtered[0][0].snr:deconvolved_filtered[0][0].snr >> bothSNRWais","flags.pulser == 1","colz");
  summaryTree->Draw("coherent_filtered[0][0].snr:deconvolved_filtered[0][0].snr >> bothSNRLDB","flags.pulser == 2","colz");

  outFile->cd();
  bothSNR->Write();
  bothSNRWais->Write();
  bothSNRLDB->Write();

  /*
    ----------------*/

  return;

}


  /*=================
    Polarization Plots
  */

void plotPol(TChain* summaryTree,TFile* outFile) {

  cout << "Doing polarization stuff" << endl;


  TH2D *polariz = new TH2D("polariz","Polarization Info - Noise; Linear Polarization Fraction; Linear Polarization Angle",
			   500,0,1,360,-45,45);
  TH2D *polarizWais = new TH2D("polarizWais","Polarization Info - Wais; Linear Polarization Fraction; Linear Polarization Angle",
			       500,0,1,360,-45,45);
  TH2D *polarizLDB = new TH2D("polarizLDB","Polarization Info - LDB; Linear Polarization Fraction; Linear Polarization Angle",
			      500,0,1,360,-45,45);



  summaryTree->Draw("deconvolved_filtered[0][1].linearPolAngle():deconvolved_filtered[0][1].linearPolFrac() >> polariz","flags.pulser == 0","colz");

  summaryTree->Draw("deconvolved_filtered[0][1].linearPolAngle():deconvolved_filtered[0][1].linearPolFrac() >> polarizWais","flags.pulser == 1","colzSame");

  summaryTree->Draw("deconvolved_filtered[0][1].linearPolAngle():deconvolved_filtered[0][1].linearPolFrac() >> polarizLDB","flags.pulser == 2","colzSame");



  outFile->cd();
  polariz->Write();
  polarizWais->Write();
  polarizLDB->Write();

  TH1D *linPolFrac = polariz->ProjectionX("linPolFrac",0,500);
  linPolFrac->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");
  TH1D *linPolFracWais = polarizWais->ProjectionX("linPolFracWais",0,500);
  linPolFracWais->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");
  TH1D *linPolFracLDB = polarizLDB->ProjectionX("linPolFracLDB",0,500);
  linPolFracLDB->SetTitle("Linear Polarization Fraction - No Pulsers; Linear Polarization Fraction; Occupancy");

  TH1D *linPolAng = polariz->ProjectionY("linPolAng",0,360);
  linPolAng->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");
  TH1D *linPolAngWais = polarizWais->ProjectionY("linPolAngWais",0,360);
  linPolAngWais->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");
  TH1D *linPolAngLDB = polarizLDB->ProjectionY("linPolAngLDB",0,360);
  linPolAngLDB->SetTitle("Linear Polarization Angtion - No Pulsers; Linear Polarization Angle; Occupancy");


  outFile->cd();
  linPolFrac->Write();
  linPolFracWais->Write();
  linPolFracLDB->Write();

  TH1 *linPolFracCum = makeNormCumulative(linPolFrac);
  linPolFracCum->Write();
  TH1 *linPolFracCumWais = makeNormCumulative(linPolFracWais);
  linPolFracCum->Write();
  TH1 *linPolFracCumLDB = makeNormCumulative(linPolFracLDB);
  linPolFracCumLDB->Write();


  linPolAng->Write();
  linPolAngWais->Write();
  linPolAngLDB->Write();

  TH1 *linPolAngCum = makeNormCumulative(linPolAng);
  linPolAngCum->Write();
  TH1 *linPolAngCumWais = makeNormCumulative(linPolAngWais);
  linPolAngCum->Write();
  TH1 *linPolAngCumLDB = makeNormCumulative(linPolAngLDB);
  linPolAngCumLDB->Write();


  return;
  /*
    ---------*/
}


void PlotCutsAndBases(TChain* summaryTree, TFile *outFile) {  
  TH1D *bothCutsWais = new TH1D("bothCutsWais","bothCutsWais",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCutsWais","flags.pulser == 1 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 ","colz");


  TH1D *bothCuts = new TH1D("bothCuts","bothCuts",1000,0,9e7);

  summaryTree->Draw("eventNumber >> bothCuts","flags.pulser == 0 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5","colz");


  TH2D *cutLatLon = new TH2D("cutLatLon","cutLatLon",250,-90,-65,360,-180,180);

  summaryTree->Draw("peak[0][0].longitude:peak[0][0].latitude >> cutLatLon","flags.pulser == 0 && TMath::Abs(TMath::RadToDeg()*TMath::ATan(coherent[0][0].U/coherent[0][0].Q)/2) < 15 && TMath::Sqrt(pow(coherent[0][0].Q,2)+pow(coherent[0][0].U,2))/coherent[0][0].I > 0.4 && templateCRayH[5] > 0.5 && peak[0][0].value > 0.05 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && TMath::Abs(peak[0][0].phi - ldb.phi) > 5 && peak[0][0].latitude > -999 && !( (peak[0][0].latitude < -78 && peak[0][0].latitude > -82) && (peak[0][0].longitude < -105 && peak[0][0].longitude > -115) )","colz");

  TFile *baseList = TFile::Open("~/anita16/local/share/anitaCalib/baseListA3.root");
  TTree *baseTree = (TTree*)baseList->Get("baseCampTree");
  TTree *awsTree  = (TTree*)baseList->Get("awsTree");

  baseTree->SetMarkerStyle(3);
  baseTree->SetMarkerColor(kRed);
  awsTree->SetMarkerStyle(3);
  awsTree->SetMarkerColor(kBlue);

  baseTree->Draw("fullLong:fullLat","","same");
  awsTree->Draw("fullLong:fullLat","","same");
  
    


  outFile->Close();

  return;
}


void waisVsWais(TChain* summaryTree, TFile* outFile) {

  TH2D *waisVsWais = new TH2D("waisVsWais","wais pulser events vs wais template",500,0,0.5, 500,0,1);

  summaryTree->Draw("templateWaisH:peak[0][0].value >> waisVsWais","flags.pulser==1","colz");

  waisVsWais->Write();


}


void plotThings(int movie=false) {

  TFile *outFile = TFile::Open("plotThings.root","recreate");
  TChain* summaryTree = (TChain*)gROOT->ProcessLine(".x loadAll.C");

  if (movie) makeMovies(summaryTree);
  plotPol(summaryTree,outFile);
  plotCorr(summaryTree,outFile);
  plotSNR(summaryTree,outFile);

  outFile->Close();

  return;
}

