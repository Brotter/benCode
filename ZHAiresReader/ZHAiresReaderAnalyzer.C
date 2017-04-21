#include "FFTtools.h"

/*

Okay so that other code (ZHIresReader) just turns the .dat files into a bunch of stupid histograms which 
I now have to deal with.  As a note, it also seems to do it wrong.

So lets:
1) Make everything into TGraphs (ZHAiresParser does this now)
2) find the "peak cone" antenna location (whatever has the peak amplitude)
3) Convolve all the waveofmrs in the histogram with the system impulse response
4) Correlate them all and see how bad a template trigger will be at correlation due to shower geometry


 */


void ZHAiresReaderAnalyzer() {

  TFile *inFile = TFile::Open("ZHAiresReader_zen60.root");



  TH2D *hWaves = (TH2D*)inFile->Get("hWavesZ");

  const int antsPerLeg = 9*4;
  
  const int numAntennas = 16*9;
  TGraph* gWavesX[numAntennas];


  //import all the simulated pulses and generate their time bases
  const int initNumPts = 2000;
  double dT = 0.3; //from sim .inp file

  for (int ant=0; ant<numAntennas; ant++) {
    gWavesX[ant] = new TGraph();
    int maxBin = hWaves->ProjectionX("temp",ant,ant)->GetMaximumBin();
    for (int pt=0; pt<initNumPts; pt++) {
      int adjustedPt = pt+maxBin-100;
      gWavesX[ant]->SetPoint(pt,dT*pt,hWaves->GetBinContent(adjustedPt,ant)*100);
    }
  }
  //interpolate them
  TGraph* gWavesInterp[numAntennas];
  for (int ant=0; ant<numAntennas; ant++) {
    gWavesInterp[ant] = FFTtools::getInterpolatedGraph(gWavesX[ant],0.1);
  }
  dT = 0.1;  //dT is different once I interpolate!
  //zero pad them a ton
  const int numPts = 10000;
  for (int ant=0; ant<numAntennas; ant++) {
    while (gWavesInterp[ant]->GetN() < numPts) {
      gWavesInterp[ant]->SetPoint(gWavesInterp[ant]->GetN(),gWavesInterp[ant]->GetN()*dT,0);
    }
  }

  //find the maximum pulse height of all those
  int x,maxAnt,z;
  hWaves->GetBinXYZ(hWaves->GetMaximumBin(),x,maxAnt,z);


  //get the impulse response
  TGraph *impulse = new TGraph("/Users/brotter/benCode/impulseResponse/integratedTF/transferFunctions/normTF_01BV.txt");

  /*
  cout << "impulse dT:" << impulseRaw->GetX()[1]-impulseRaw->GetX()[0] << endl;
  //shift it because it is wrapped stupidly
  TGraph *impulse = new TGraph();
  for (int pt=700; pt<impulseRaw->GetN(); pt++) {
    impulse->SetPoint(impulse->GetN(),impulse->GetN()*dT,impulseRaw->GetY()[pt]);
  }
  for (int pt=0; pt<700; pt++) {
    impulse->SetPoint(impulse->GetN(),impulse->GetN()*dT,impulseRaw->GetY()[pt]);
  }
  */
  //zero pad it so it is the same length as the sim
  while (impulse->GetN() < gWavesInterp[1]->GetN()) {
    impulse->SetPoint(impulse->GetN(),impulse->GetN()*dT,0);
  }
  


  //convolve the impulse response with the simulated pulse
  TGraph* gWavesConv[numAntennas];
  for (int ant=0; ant<numAntennas; ant++) {
    gWavesConv[ant] = FFTtools::getConvolution(impulse,gWavesInterp[ant]);
    for (int pt=0; pt<4000; pt++) {
      gWavesConv[ant]->GetY()[pt] = gWavesConv[ant]->GetY()[pt+5000];
      gWavesConv[ant]->GetY()[pt+5000] = 0;
    }
  }


  //Draw all that stuff!
  TCanvas *cConv = new TCanvas("cConv","cConv",1000,800);
  stringstream filename;
  for (int drawAnt = 1; drawAnt<36; drawAnt++) {
    //    int drawAnt = 14;
  cConv->Clear();
  cConv->Divide(1,3);
  TLegend *cConvLeg = new TLegend(0.6,0.7,0.9,0.9);
  cConv->cd(1);

  impulse->SetLineColor(kBlue);
  impulse->SetTitle("System Impulse Response; ; Amplitude (Volts-ish)");
    impulse->GetXaxis()->SetRangeUser(0,150);
  cConvLeg->AddEntry(impulse,"System response","l");
  impulse->Draw("al");

  cConv->cd(2);
  cConvLeg->AddEntry(gWavesConv[drawAnt],"Convolved Signal","l");
  gWavesConv[drawAnt]->GetXaxis()->SetRangeUser(0,150);
  gWavesConv[drawAnt]->SetTitle("Convolution; ;Amplitude (Volts-ish)");
  gWavesConv[drawAnt]->Draw("al");
  //gWavesConv[drawAnt]->GetYaxis()->SetRangeUser(-0.3,0.3);

  cConv->cd(3);
  gWavesInterp[drawAnt]->SetLineColor(kRed);
  gWavesInterp[drawAnt]->SetTitle("Simulated E_y Field; Time (ns); Electric Field (V/m)");
  gWavesInterp[drawAnt]->GetXaxis()->SetRangeUser(0,150);
  cConvLeg->AddEntry(gWavesInterp[drawAnt],"Simulated Ey field","l");
  gWavesInterp[drawAnt]->Draw("al");
  gWavesInterp[drawAnt]->GetYaxis()->SetRangeUser(-0.3,0.3);

  cConv->cd(1);
  cConvLeg->Draw("Same");

  if (drawAnt==35) cConv->SaveAs("cone.gif++");
  else cConv->SaveAs("cone.gif+");
  }


  //correlate the maximum waveform with the rest of them
  TGraph* gWavesCorr[numAntennas];
  for (int ant=0; ant<numAntennas; ant++) {
    gWavesCorr[ant] = new TGraph(numPts,gWavesConv[ant]->GetX(),
				 FFTtools::getCorrelation(gWavesConv[ant],gWavesConv[maxAnt],0,numPts));
  }

  //make a TGraph of the maximum correlations
  TGraph* maxCorr = new TGraph();
  TGraph* minCorr = new TGraph();
  for (int ant=0; ant<numAntennas; ant++) {
    maxCorr->SetPoint(ant,ant,TMath::MaxElement(gWavesCorr[ant]->GetN(),gWavesCorr[ant]->GetY()));
    minCorr->SetPoint(ant,ant,TMath::MinElement(gWavesCorr[ant]->GetN(),gWavesCorr[ant]->GetY()));
  }

  //normalize that whole thing
  double maxCorrValue = TMath::MaxElement(numAntennas,maxCorr->GetY());
  double minCorrValue = TMath::MinElement(numAntennas,minCorr->GetY());
  for (int ant=0; ant<numAntennas; ant++) {
    maxCorr->GetY()[ant] /= maxCorrValue;
    minCorr->GetY()[ant] /= maxCorrValue;
  }

  //separate the correlations into "legs" so you can tell what is happening
  for (int leg=0; leg<4; leg++) {
    for (int ant=0; ant<antsPerLeg; ant++) {
      int index = leg*antsPerLeg+ant;
      if (leg==0) gWavesCorr[index]->SetLineColor(kRed);
      if (leg==1) gWavesCorr[index]->SetLineColor(kBlue);
      if (leg==2) gWavesCorr[index]->SetLineColor(kGreen);
      if (leg==3) gWavesCorr[index]->SetLineColor(kBlack);
      for (int pt=0; pt<gWavesCorr[index]->GetN(); pt++) {
	gWavesCorr[index]->GetY()[pt] += 0.00005*index;
      }
      int lenPts = gWavesCorr[index]->GetN();
      double endT = gWavesCorr[index]->GetX()[lenPts-1];
      for (int pt=0; pt<lenPts/2; pt++) {
	gWavesCorr[index]->SetPoint(gWavesCorr[index]->GetN(),
				    gWavesCorr[index]->GetX()[pt]+endT,
				    gWavesCorr[index]->GetY()[pt]);
      }
      
    }//end ant
  }//end leg


  //separate the peaks into legs too
  TGraph *legMaxCorr[4];
  TGraph *legMinCorr[4];
  for (int leg=0; leg<4; leg++) {
    legMaxCorr[leg] = new TGraph();
    legMinCorr[leg] = new TGraph();
    legMinCorr[leg]->SetMarkerStyle(kStar);
    legMinCorr[leg]->SetMarkerSize(1);
    for (int ant=0; ant<antsPerLeg; ant++) {
      int index = leg*antsPerLeg+ant;
      legMaxCorr[leg]->SetPoint(ant,maxCorr->GetX()[index],maxCorr->GetY()[index]);
      legMinCorr[leg]->SetPoint(ant,minCorr->GetX()[index],-minCorr->GetY()[index]);
    }
  }
  legMaxCorr[0]->SetMarkerColor(kRed);
  legMaxCorr[1]->SetMarkerColor(kBlue);
  legMaxCorr[2]->SetMarkerColor(kGreen);
  legMaxCorr[3]->SetMarkerColor(kBlack);
  legMaxCorr[0]->SetLineColor(kRed);
  legMaxCorr[1]->SetLineColor(kBlue);
  legMaxCorr[2]->SetLineColor(kGreen);
  legMaxCorr[3]->SetLineColor(kBlack);
  legMinCorr[0]->SetMarkerColor(kRed);
  legMinCorr[1]->SetMarkerColor(kBlue);
  legMinCorr[2]->SetMarkerColor(kGreen);
  legMinCorr[3]->SetMarkerColor(kBlack);
  legMinCorr[0]->SetLineColor(kRed);
  legMinCorr[1]->SetLineColor(kBlue);
  legMinCorr[2]->SetLineColor(kGreen);
  legMinCorr[3]->SetLineColor(kBlack);


  


  //Draw the "correlation fraction" graph
  TCanvas *c2 = new TCanvas("c2","c2",800,400);
  c2->Divide(2);
  c2->cd(1);
  gWavesCorr[1]->Draw("al");
  gWavesCorr[1]->GetYaxis()->SetRangeUser(0,0.01);
  //  gWavesCorr[0]->GetXaxis()->SetRangeUser(200,350);
  for (int ant=2; ant<144; ant++) {
    gWavesCorr[ant]->Draw("lSame");
  }
  c2->cd(2);
  legMaxCorr[0]->GetXaxis()->SetLimits(0,numAntennas);
  legMaxCorr[0]->Draw("alp");
  legMaxCorr[0]->GetYaxis()->SetRangeUser(-1.1,1.1);
  legMaxCorr[1]->Draw("lpSame");
  legMaxCorr[2]->Draw("lpSame");
  legMaxCorr[3]->Draw("lpSame");
  legMinCorr[0]->Draw("lpSame");
  legMinCorr[1]->Draw("lpSame");
  legMinCorr[2]->Draw("lpSame");
  legMinCorr[3]->Draw("lpSame");


  //make a color coded antenna map
  // **Note, the z's returned are from Insertion Alt, which is 100km for these sims, so invert them
  //       x and y are from shower core on ground.
  TGraph2D *fullAntMap = (TGraph2D*)inFile->Get("antPos");
  TGraph2D *legAntMap[4];
  for (int leg=0; leg<4; leg++) {
    legAntMap[leg] = new TGraph2D();
    for (int ant=0; ant<antsPerLeg; ant++) {
      int index = leg*antsPerLeg+ant;
      fullAntMap->GetZ()[index] -= 100000;
      fullAntMap->GetZ()[index] *= -1;
      legAntMap[leg]->SetPoint(ant,fullAntMap->GetX()[index],fullAntMap->GetY()[index],
			       fullAntMap->GetZ()[index] );
    }
  }

  legAntMap[0]->SetMarkerColor(kRed);
  legAntMap[1]->SetMarkerColor(kBlue);
  legAntMap[2]->SetMarkerColor(kGreen);
  legAntMap[3]->SetMarkerColor(kBlack);



  //Draw that map
  TCanvas *cAntMap = new TCanvas("cAntMap","cAntMap",800,800);
  fullAntMap->Draw("P0");
  legAntMap[0]->Draw("PSame");
  legAntMap[1]->Draw("PSame");
  legAntMap[2]->Draw("PSame");
  legAntMap[3]->Draw("PSame");



  //Also draw a crude "heatmap" of correlation
  TH2D *heatMap = new TH2D("heatMap","Correlation heatMap",61,-1500,1500, 70,-6500,-3000);
  for (int ant=0; ant<numAntennas; ant++) {
    double x = fullAntMap->GetX()[ant];
    double y = fullAntMap->GetY()[ant];
    double peak = maxCorr->GetY()[ant];
    heatMap->Fill(x,y,peak);
  }
  TCanvas *cHeatMap = new TCanvas("cHeatMap","cHeatMap",800,800);
  cHeatMap->SetLogz();
  heatMap->Draw("colz");
  

  //Also draw a crude "heatmap" of peak amplitude?
  TH2D *heatMap2 = new TH2D("heatMap2","Amplitude heatMap",61,-1500,1500, 70,-6500,-3000);
  for (int ant=0; ant<numAntennas; ant++) {
    double x = fullAntMap->GetX()[ant];
    double y = fullAntMap->GetY()[ant];
    double peak = TMath::MaxElement(gWavesInterp[ant]->GetN(),gWavesInterp[ant]->GetY());
    heatMap2->Fill(x,y,peak);
  }
  TCanvas *cHeatMap2 = new TCanvas("cHeatMap2","cHeatMap2",800,800);
  cHeatMap2->SetLogz();
  heatMap2->Draw("colz");



  //slant depth from :http://inspirehep.net/record/1125935/files/arXiv%3A1208.1171.pdf
  // X = X0*exp(-h/h0) , X0 = 1030g/cm^2 , h0 = 8400m
  //for zen60, X = 380.844, 
  // h = -h0*ln(X/X0)

  double xMaxAlt = TMath::Log(380.844/1030.);
  cout << "xMaxAlt = " << xMaxAlt << endl; 
  

  //get a wais pulse and compare against that
  TFile *waisFile = TFile::Open("/Users/brotter/anita16/benPrograms/waisPulses/waisCorrelationTest.root");
  TGraph *waisPulse = (TGraph*)waisFile->Get("Phi15Ring2");
  cout << waisPulse->GetN() << endl;
  while (waisPulse->GetN() < numPts) {
    waisPulse->SetPoint(waisPulse->GetN(),waisPulse->GetX()[waisPulse->GetN()-1]+0.1,0);
  }


  

  TGraph *waisCorr = new TGraph();

  double* waisAutoCorr = FFTtools::getCorrelation(waisPulse,waisPulse,0,numPts);
  double waisAutoCorrMax = TMath::MaxElement(numPts,waisAutoCorr);

  for (int ant=0; ant<numAntennas; ant++) {
    double* convAutoCorr = FFTtools::getCorrelation(gWavesConv[ant],gWavesConv[ant],0,numPts);
    double convAutoCorrMax = TMath::MaxElement(numPts,convAutoCorr);
    double *temp = FFTtools::getCorrelation(waisPulse,gWavesConv[ant],0,numPts);
    double scale = (convAutoCorrMax*waisAutoCorrMax);
    if (scale == 0) continue;
    waisCorr->SetPoint(ant,ant,TMath::MaxElement(numPts,temp)*scale);
  }

    

  TCanvas* cWais = new TCanvas("cWais","cWais",800,600);
  cWais->Divide(1,3);
  cWais->cd(1);
  waisPulse->Draw("alp");
  cWais->cd(2);
  gWavesConv[maxAnt]->Draw("alp");
  cWais->cd(3);
  waisCorr->Draw("alp");




  return;
}
