/*

  The calibration data has runs with terminated inputs, 5dB ENR, and 15dB ENR.  

  From those and the average event RMS I should be able to do do some simple y-factor analysis

  Also... what if I do some NOT simple y-factor analysis?  

  Like, RMS is basically energy in the waveform, so you should be able to that as a function of frequency from the fft magnitude

  So these will be 2d histograms, x axis frequency, y axis magnitude at that frequency per waveform
  
  need to evenly sample them first though and zero padding is wrong (no noise) so lets _decrease_ the number of samples
  (shorter record length = larger f bins ; any binning is better than the rms)

 */

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <random>
//root
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TEllipse.h" 
#include "TMarker.h" 
#include "TStyle.h" 
#include "TCanvas.h"
//anita
#include "AnitaDataset.h"
#include "FFTtools.h"
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"

using namespace std;


//this can be a global since it should never change, just the "shortest" a wave can be
int recoLen = 230;

TGraph* makePowerSpectrumMilliVoltsNanoSecondsdBFreqDomain(int length, FFTWComplex *theFFT, double deltaT=1./2.6) {

  int newLength=(length/2)+1;

  double *newY = new double [newLength];
  double *newX = new double [newLength];
  
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  //    std::cout  << deltaF << "\t" << deltaT << "\t" << length << std::endl;

  double tempF=0;
  for(int i=0;i<newLength;i++) {
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    power /= pow(length,2); //normalize by N^2
    if (length%2 == 0) power *=2; // if it is even all the bins get multiplied by two
    else if (i>0 && i<newLength-1) power*=2; //Changed form 2 by RJN 29/01/10 //account for symmetry
    //  power*=deltaT/(length); //For time-integral squared amplitude
    //	power/=(1e3*1e3*1e9); //This bit converts from mv*mv*ns to v*v*s
    //    power/=deltaF;//Just to normalise bin-widths
    //Ends up the same as dt^2, need to integrate the power (multiply by df)
    //to get a meaningful number out.	
    
    if (power>0 ) power=10*TMath::Log10(power);
    else power=-1000; //no reason
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;
  }
  
  
  TGraph *grPower = new TGraph(newLength,newX,newY);
  delete [] newY;
  delete [] newX;
  return grPower;
  
}


int yFactorFreqDomain(int run, int surf, int chan, int maxEntries,
		      TH2D *currHist,TH2D *currHistLin, TH1D *hRMS) {

  int chanIndex =  AnitaGeomTool::getChanIndex(surf,chan);

  AnitaDataset *data = new AnitaDataset(run,false,WaveCalType::kOnlyDTs);
  if (!data->fRunLoaded) {
    cout << "WARNING! yFactorFreqDomain(): Couldn't load run " << run << endl;
    return -1;
  }

  int lenEntries;
  if (maxEntries<0) {
    lenEntries = data->N();}
  else {
    lenEntries = maxEntries; }


  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10 == 0) cout << entry << "/" << lenEntries << "\r"; fflush(stdout);
    data->getEntry(entry);

    UsefulAnitaEvent *useful = new UsefulAnitaEvent(data->calibrated());

    hRMS->Fill(useful->rms[chanIndex]);

    TGraph *waveformRaw = useful->getGraphFromSurfAndChan(surf,chan);
    delete useful;
    TGraph *waveform = FFTtools::getInterpolatedGraph(waveformRaw,1./2.6);
    delete waveformRaw;
    //      if (waveform->GetN() < recoLen) cout << "entry " << entry << " record length too low!" << "(" << waveform->GetN() << ")" << endl;
    //      while (waveform->GetN() > recoLen) waveform->RemovePoint(waveform->GetN() - 1); //this is super inefficient
    
    FFTWComplex *fft = FFTtools::doFFT(recoLen,waveform->GetY());
    delete waveform;
    
    TGraph *magnitude = makePowerSpectrumMilliVoltsNanoSecondsdBFreqDomain(recoLen,fft);
    delete [] fft;
    //      if ( magnitude->GetN() != (recoLen/2 + 1) ) cout << " magnitude->GetN() != (recoLen/2 + 1)!! " <<  magnitude->GetN() << " " <<  (recoLen/2 + 1) << endl;
    //      cout << "magnitude->GetN() = " << magnitude->GetN() << endl;
    
    
    double* y = magnitude->GetY();
    double* x = magnitude->GetX();
    for (int pt=0; pt<magnitude->GetN(); pt++) {
      currHist->Fill(x[pt],y[pt]);
      double yLin = TMath::Sqrt(pow(10.,y[pt]/10.));
      currHistLin->Fill(x[pt],yLin);
    }
    
    delete magnitude;
    
  }
  delete data;

  cout << endl;

  return 1;
  
}
      


int main(int argc, char** argv) {
  
  int maxEntries = -1;
  string fileNameBase = "yFactorFreqDomain";

  if (argc == 2) {
    maxEntries = atoi(argv[1]);
  }
  
  if (argc == 3) {
    maxEntries = atoi(argv[1]);
    fileNameBase = argv[2];
  }
  else {
    cout << "Usage: " << argv[0] << " [opt: maxEntries] [opt: fileNameBase" << endl;
  }
  
  cout << "Starting!  Using " << maxEntries << " events per channel and output file name: " << fileNameBase << ".root" << endl;




  FFTtools::loadWisdom("/Users/brotter/macros/fftWisdom.dat");  

  //this is calibration data, and AnitaDataset looks to ${ANITA3_DATA} for the location, not where it actually is
  char dataDir[512] = {"ANITA_ROOT_DATA=/Volumes/ANITA3Data/ANITA3_calibrationFiles/antarctica14/root/"};
  putenv(dataDir);


  //get the run list of noise source measurements
  ifstream ifsLog("usedRunList_single.txt");
  string line;
  
  //set some variables to be filled
  int surf=0;
  int chan=0;
  int run=0;
  string antName;
  string calType;
  string temp;

  //create the file
  stringstream name;
  name.str("");
  name << fileNameBase << ".root";
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");

  while ( ifsLog.good()) {
    getline(ifsLog,line);
    stringstream lineStream(line);
    lineStream >> antName;
    lineStream >> temp;
    surf = atoi(temp.c_str())-1;
    lineStream >> temp;
    chan = atoi(temp.c_str())-1;
    lineStream >> calType;
    lineStream >> temp;
    run = atoi(temp.c_str());
    
    cout << "got: " << antName << " surf:" << surf << " chan:" << chan << " run:" << run << " calType:" << calType << endl;
    name.str("");
    name << "surf" << surf << "Chan" << chan << "_" << calType;
    TH2D *currHist = new TH2D(name.str().c_str(),name.str().c_str(),recoLen/2. + 1,0,1301, 201,-100,100);
    name << "_lin";
    TH2D *currHistLin = new TH2D(name.str().c_str(),name.str().c_str(),recoLen/2. + 1, 0,1301,250,0,100);

    name.str("");
    name << "surf" << surf << "Chan" << chan << "_" << calType;
    name << "_rms";
    TH1D *hRMS = new TH1D(name.str().c_str(),name.str().c_str(),200,0,100);

    int status = yFactorFreqDomain(run,surf,chan,maxEntries,
				   currHist,currHistLin,hRMS);
    if (status > 0) {
      cout << "Writing to file..." << endl;
      outFile->cd();
      currHist->Write();
      currHistLin->Write();
      hRMS->Write();
    }
    else {
      ofstream brokenRuns;
      brokenRuns.open("brokenRuns.txt",ofstream::out | ofstream::app);
      brokenRuns << "run" << run << endl;
      brokenRuns.close();
    }
    cout << "Done with that channel!" << endl;
    
    delete currHist;
    delete currHistLin;
    delete hRMS;

  }
  cout << "Done with loop" << endl;

  
  outFile->Close();
  cout << "Closed output file" << endl;

  cout << "Writing wisdom" << endl;
  FFTtools::saveWisdom("/Users/brotter/macros/fftWisdom.dat");

  cout << "done" << endl;


  return 1;

}
