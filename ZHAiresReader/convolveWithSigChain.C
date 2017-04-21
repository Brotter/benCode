#include "FFTtools.h"

/*

  does the same thing as ZHAiresReaderAnalyzer, except with the ZHAiresParser TGraph outputs which are actually right

  Ben Rotter - Apr 2017 - University of Hawaii at Manoa
 */


void convolveWithSigChain(int energy=195) {


  stringstream name;
  name.str("");
  name << "ZHAiresParser_En" << energy << ".root";

  cout << "opening: " << name.str() << endl;
  TFile *inFile = TFile::Open(name.str().c_str());

  const int numSubs = 4;
  const int numAntennas = 8;
  const int numWaveforms = numSubs*numAntennas;

  TGraph* gWavesX[numAntennas];


  //import all the simulated pulses and generate their time bases
  const int initNumPts = 2000;
  double dT = 0.3; //from sim .inp file

  for (int sub=0; sub<numSubs; sub++) {
    for (int ant=0; ant<numAntennas; ant++) {
      int index = sub*numAntennas + ant;

      name.str("");
      name << "off0sub" << sub << "ant" << ant+1 << "_X";

      gWavesX[index] = (TGraph*)inFile->Get(name.str().c_str());

      cout << index << " " << name.str() << " " << gWavesX[index] << endl;

    }
  }


  //interpolate them
  TGraph* gWavesInterp[numWaveforms];
  for (int wave=0; wave<numWaveforms; wave++) {
    //    cout << wave << " " << gWavesX[wave]->GetN() << endl;

    gWavesInterp[wave] = FFTtools::getInterpolatedGraph(gWavesX[wave],0.1);
  }
  dT = 0.1;  //dT is different once I interpolate!

  
  //find the maximum one (for plotting purposes)

  int maxWave;
  double maxVal = 0;
  for (int wave=0; wave<numWaveforms; wave++) {
    double val = TMath::Abs(TMath::MinElement(gWavesInterp[wave]->GetN(),gWavesInterp[wave]->GetY()));
    //    cout << wave << " " << val << " " << maxVal << endl;
    if (val > maxVal) {
      maxWave = wave;
      maxVal = val;
    }
  }

  cout << "maxWave=" << maxWave << " maxVal=" << maxVal << endl;
      


  //zero pad them a ton
  const int numPts = 10000;
  for (int wave=0; wave<numWaveforms; wave++) {
    while (gWavesInterp[wave]->GetN() < numPts) {
      gWavesInterp[wave]->SetPoint(gWavesInterp[wave]->GetN(),gWavesInterp[wave]->GetN()*dT,0);
    }
  }


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
  TGraph* gWavesConv[numWaveforms];
  for (int wave=0; wave<numWaveforms; wave++) {
    gWavesConv[wave] = FFTtools::getConvolution(impulse,gWavesInterp[wave]);
    for (int pt=0; pt<4000; pt++) {
      gWavesConv[wave]->GetY()[pt] = gWavesConv[wave]->GetY()[pt+5000];
      gWavesConv[wave]->GetY()[pt+5000] = 0;
    }
  }


  //Draw all that stuff!
  TCanvas *cConv = new TCanvas("cConv","cConv",1000,800);
  stringstream filename;
  for (int drawWave = 1; drawWave<numWaveforms; drawWave++) {
    //    int drawWave = 14;
  cConv->Clear();
  cConv->Divide(1,3);
  TLegend *cConvLeg = new TLegend(0.6,0.7,0.9,0.9);
  cConv->cd(1);

  impulse->SetLineColor(kBlue);
  impulse->SetTitle("System Impulse Response; ; Amplitude (Volts-ish)");
  impulse->GetXaxis()->SetRangeUser(0,100);
  cConvLeg->AddEntry(impulse,"System response","l");
  impulse->Draw("al");

  cConv->cd(2);
  cConvLeg->AddEntry(gWavesConv[drawWave],"Convolved Signal","l");
  gWavesConv[drawWave]->GetXaxis()->SetRangeUser(200,400);
  gWavesConv[drawWave]->SetTitle("Convolution; ;Amplitude (Volts-ish)");
  gWavesConv[drawWave]->Draw("al");
  gWavesConv[drawWave]->GetYaxis()->SetRangeUser(-.03,.03);

  cConv->cd(3);
  gWavesInterp[drawWave]->SetLineColor(kRed);
  gWavesInterp[drawWave]->SetTitle("Simulated E_y Field; Time (ns); Electric Field (V/m)");
  gWavesInterp[drawWave]->GetXaxis()->SetRangeUser(200,400);
  cConvLeg->AddEntry(gWavesInterp[drawWave],"Simulated Ey field","l");
  gWavesInterp[drawWave]->Draw("al");
  gWavesInterp[drawWave]->GetYaxis()->SetRangeUser(-0.001,0.0002);

  cConv->cd(1);
  cConvLeg->Draw("Same");

  if (drawWave==35) cConv->SaveAs("cone.gif++");
  else cConv->SaveAs("cone.gif+");
  }


  //correlate the maximum waveform with the rest of them
  TGraph* gWavesCorr[numWaveforms];
  for (int wave=0; wave<numWaveforms; wave++) {
    gWavesCorr[wave] = new TGraph(numPts,gWavesConv[wave]->GetX(),
				 FFTtools::getCorrelation(gWavesConv[wave],gWavesConv[maxWave],0,numPts));
  }

  //make a TGraph of the maximum correlations
  TGraph* maxCorr = new TGraph();
  TGraph* minCorr = new TGraph();
  for (int wave=0; wave<numWaveforms; wave++) {
    maxCorr->SetPoint(wave,wave,TMath::MaxElement(gWavesCorr[wave]->GetN(),gWavesCorr[wave]->GetY()));
    minCorr->SetPoint(wave,wave,TMath::MinElement(gWavesCorr[wave]->GetN(),gWavesCorr[wave]->GetY()));
  }

  //normalize that whole thing
  double maxCorrValue = TMath::MaxElement(numWaveforms,maxCorr->GetY());
  double minCorrValue = TMath::MinElement(numWaveforms,minCorr->GetY());
  for (int wave=0; wave<numWaveforms; wave++) {
    maxCorr->GetY()[wave] /= maxCorrValue;
    minCorr->GetY()[wave] /= maxCorrValue;
  }



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

  for (int wave=0; wave<numWaveforms; wave++) {
    double* convAutoCorr = FFTtools::getCorrelation(gWavesConv[wave],gWavesConv[wave],0,numPts);
    double convAutoCorrMax = TMath::MaxElement(numPts,convAutoCorr);
    double *temp = FFTtools::getCorrelation(waisPulse,gWavesConv[wave],0,numPts);
    double scale = (convAutoCorrMax*waisAutoCorrMax);
    if (scale == 0) continue;
    waisCorr->SetPoint(wave,wave,TMath::MaxElement(numPts,temp)*scale);
  }

    

  TCanvas* cWais = new TCanvas("cWais","cWais",800,600);
  cWais->Divide(1,3);
  cWais->cd(1);
  waisPulse->Draw("alp");
  cWais->cd(2);
  gWavesConv[maxWave]->Draw("alp");
  cWais->cd(3);
  waisCorr->Draw("alp");

  name.str("");
  name << "peakTemplate_En" << energy << ".txt";
  ofstream out(name.str().c_str());
  for (int pt=0; pt<gWavesConv[maxWave]->GetN(); pt++) {
    out << gWavesConv[maxWave]->GetX()[pt] << " " << gWavesConv[maxWave]->GetY()[pt] << endl;
  }


  return;
}
