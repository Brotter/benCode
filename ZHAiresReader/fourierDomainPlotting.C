#include "FFTtools.h"


void fourierDomainPlotting(int energy=195){

  /*

    Lets plot the fourier transforms of the different showers?


   */




  
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

    gWavesInterp[wave] = FFTtools::getInterpolatedGraphFreqDom(gWavesX[wave],0.1);
  }
  dT = 0.1;  //dT is different once I interpolate!



  TFile *outFile = TFile::Open("fourierDomainPlotting.root","recreate");


  TF1 *f1 = new TF1("f1","[0]*TMath::Power(1,x)*(TMath::Exp(-1)/TMath::Factorial(x))", 0, 1000);

  TCanvas *cConv = new TCanvas("cConv","cConv",1000,800);
  TGraph *gWavesSpecMag[numWaveforms];
  for (int wave=0; wave<numWaveforms; wave++) {
    gWavesSpecMag[wave] = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(gWavesInterp[wave]);
    gWavesSpecMag[wave]->Draw("alp");
    gWavesSpecMag[wave]->Fit("f1");
    gWavesSpecMag[wave]->GetXaxis()->SetRangeUser(0,2000);
    
    name.str("");
    name << "wave" << wave;
    gWavesSpecMag[wave]->SetName(name.str().c_str());
    outFile->cd();
    gWavesSpecMag[wave]->Write();
    //    if (wave==numWaveforms-1) cConv->SaveAs("spectrum.gif++");
    //    else cConv->SaveAs("spectrum.gif+");
    cConv->Update();
    cConv->Clear();

  }

  outFile->Close();



  return;
}


  












