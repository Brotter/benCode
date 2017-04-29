/*


  Plotting and analysis tools to use on the output file from yFactorFreqDomain


 */


string fileNameGlob = "yFactorFreqDomain.root";

TF1 *makeRiceFit() {

  TF1 *riceFit = new TF1("riceFit","[3]^2 * (x) * TMath::Exp( -1*((x)^2 + [2]^2) / (2*[0]^2)) * TMath::BesselI0((x)*[2]/[0]^2)",0,100);
  riceFit->SetParameters(6.8,10,3);
  riceFit->SetParLimits(0,0,20);
  riceFit->SetParLimits(1,0,20);
  riceFit->SetParLimits(1,0,20);
  return riceFit;
}
 


TGraph *histToMeanGraph(TH2D* inGraphRaw) {
  
  TH2D *inGraph = (TH2D*)inGraphRaw->Clone();
  TGraph *outGraph = new TGraph();
  int numBins = inGraph->GetNbinsX();


  //rician fit function

  for (int bin=0; bin<numBins; bin++) {
    TH1D *temp = inGraph->ProjectionY("temp",bin+1,bin+1);
    TF1 *riceFit = makeRiceFit();
    temp->Fit(riceFit,"QN");
    double peak = riceFit->GetMaximumX();
    delete riceFit;
    //    delete temp;
    //    cout << bin << " " << peak << endl;
    outGraph->SetPoint(bin,inGraph->GetXaxis()->GetBinCenter(bin+1),peak);
  }

  delete inGraph;
  
  return outGraph;

}


void plotSurfChan(int surf, int chan, string fileName=fileNameGlob) {

  TFile *inFile = TFile::Open(fileName.c_str());

  stringstream name;

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Big_lin";
  TH2D *hBig = (TH2D*)inFile->Get(name.str().c_str());
  hBig->SetTitle("+15dBENR;Frequency (MHz); Linear Power (mW)");
  hBig->SetMarkerColor(kRed);
  hBig->GetYaxis()->SetRangeUser(0,80);
  hBig->SetStats(0);
  

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Little_lin";
  TH2D *hLittle = (TH2D*)inFile->Get(name.str().c_str());
  hLittle->SetTitle("5dBENR;Frequency (MHz); Linear Power (mW)");
  hLittle->SetMarkerColor(kBlue);
  hLittle->GetYaxis()->SetRangeUser(0,80);
  hLittle->SetStats(0);

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_None_lin";
  TH2D *hNone = (TH2D*)inFile->Get(name.str().c_str());
  hNone->SetTitle("Terminated;Frequency (MHz); Linear Power (mW)");
  hNone->GetYaxis()->SetRangeUser(0,80);
  hNone->SetStats(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(3);
  for (int i=0; i<3; i++) {
    c1->GetPad(i+1)->SetGridx(1);
    c1->GetPad(i+1)->SetGridy(1);
    c1->GetPad(i+1)->SetLogz();

  }

  cout << "poopy Pants" << endl;

  c1->cd(1);
  hNone->GetYaxis()->SetRangeUser(0,70);
  hNone->Draw("colz");
  TGraph *gNone = histToMeanGraph(hNone);
  gNone->Draw("lpSame");
  c1->cd(2);
  hLittle->GetYaxis()->SetRangeUser(0,70);
  hLittle->Draw("colz");
  TGraph *gLittle = histToMeanGraph(hLittle);
  //  for (int pt=0; pt<gLittle->GetN(); pt++) cout << gLittle->GetY()[pt] << endl;
  gLittle->Draw("lpSame");
  c1->cd(3);
  hBig->GetYaxis()->SetRangeUser(0,70);
  hBig->Draw("colz");
  TGraph *gBig = histToMeanGraph(hBig);
  gBig->Draw("lpSame");


  return;

}



void yFactor(int surf, int chan, 
	     string fileName=fileNameGlob, bool savePlots=false, bool showPlots=false,  //options
	     TGraph* gGain=NULL, TGraph* gNoiseT=NULL) { //options


  TFile *inFile = TFile::Open(fileName.c_str());

  stringstream name;

  //input graphs and temps
  name.str("");
  name << "surf" << surf << "Chan" << chan << "_None_lin";
  TH2D *hNone = (TH2D*)inFile->Get(name.str().c_str());
  TGraph *gNone = histToMeanGraph(hNone);
  double noneK = 290.0;

  
  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Little_lin";
  TH2D *hLittle = (TH2D*)inFile->Get(name.str().c_str());
  TGraph *gLittle = histToMeanGraph(hLittle);
  double littleK = noneK*(1+pow(10.0 , 5.0/10.0));

  //  cout << "noneK=" << noneK << " littleK=" << littleK << endl;

  //output graphs
  TGraph *gSlope = new TGraph();
  gSlope->SetName("gSlope");
  gSlope->SetTitle("Slope; Frequency (MHz); Slope (k*B*G)");

  TGraph *gYint = new TGraph();
  gYint->SetName("gYint");
  gYint->SetTitle("Y Intercept; Frequency (MHz); Yint (Watts)");

  if (gGain==NULL) gGain = new TGraph();
  gGain->SetName("gGain");
  gGain->SetTitle("Gain; Frequency (MHz); Gain (dB)");
  if (gGain==NULL) gNoiseT = new TGraph();
  gNoiseT->SetName("gNoiseT");
  gNoiseT->SetTitle("Noise Temperature; Frequency (MHz); Noise Temperature (K)");


  double binWidth = gLittle->GetXaxis()->GetBinWidth(0) * 1e6; //bin width in Hz
  double kBoltz = 1.38e-23; //boltzmann's constant (W*s)

  //  cout << "binWidth=" << binWidth << endl;

  if ( gNone->GetN() != gLittle->GetN() ) cout << "WARNING: uh they don't have the same length" << endl;
  for (int pt=1; pt<gLittle->GetN(); pt++) {

    double binX = gLittle->GetX()[pt];
    double yLow = pow(gNone->GetY()[pt]/(1e3),2)/230.;
    double yHigh = pow(gLittle->GetY()[pt]/(1e3),2)/230.;

    double Y = yHigh/yLow;
    double T = (littleK - Y*noneK) / (Y-1);

    double slope = (yHigh - yLow) / (littleK - noneK);
    double yInt = yLow - slope*noneK;

    //    cout << pt << " " << binX << " " << yLow << " " << yHigh << " " << slope << " " << yInt << " " << Y << " " << T << endl;
    
    gSlope->SetPoint(pt,binX,slope);
    gYint->SetPoint(pt,binX,yInt);

    double gain = slope / (binWidth*kBoltz);
    gGain->SetPoint(pt,binX,10*TMath::Log10(gain));
    double noise = (yInt / (binWidth*kBoltz*gain)) + 20;
    //    double noise = T;
    gNoiseT->SetPoint(pt,binX,noise);
  
  }


  //lets do RMS too...
  name.str("");
  name << "surf" << surf << "Chan" << chan << "_None_rms";
  TH1D *rmsNone = (TH1D*)inFile->Get(name.str().c_str());
  rmsNone->SetTitle("Surf1 Chan1 ADC RMS; RMS (ADC counts); Count");
  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Little_rms";
  TH1D *rmsLittle = (TH1D*)inFile->Get(name.str().c_str());
  rmsLittle->SetLineColor(kRed);


  gSlope->GetYaxis()->SetRangeUser(0,100e-12);
  gYint->GetYaxis()->SetRangeUser(0,10e-9);
  gGain->GetYaxis()->SetRangeUser(40,60);
  gNoiseT->GetYaxis()->SetRangeUser(0,300);

  if (showPlots && savePlots) {
    TCanvas *c1 = new TCanvas("c1","c1",1000,600);
    c1->Divide(3);
    TVirtualPad *pad1 = c1->cd(1);
    pad1->Divide(1,2);
    pad1->cd(1);
    gSlope->Draw("alp");
    pad1->cd(2);
    gGain->Draw("alp");
    TVirtualPad *pad2 = c1->cd(2);
    pad2->Divide(1,2);
    pad2->cd(1);
    gYint->Draw("alp");
    pad2->cd(2);
    gNoiseT->Draw("alp");
    TVirtualPad *pad3 = c1->cd(3);
    //  pad3->Divide(1,2);
    //  pad3->cd(1);
    TLegend *leg = new TLegend(0.4,0.8,0.9,0.9);
    rmsNone->SetStats(0);
    rmsNone->GetXaxis()->SetRangeUser(0,50);
    rmsNone->Draw();
    rmsLittle->Draw("same");
    leg->AddEntry(rmsNone,"Terminated","l");
    leg->AddEntry(rmsLittle,"5dBENR Noise","l");
    leg->SetTextSize(0.04);
    leg->Draw();
    
    
    if (savePlots) {
      name.str("");
      name << "surf" << surf << "chan" << chan << ".png";
      c1->SaveAs(name.str().c_str());
    }
  }
  inFile->Close();

  return;

}
  






/****************************************
  Test stuff
  ******************************************/

void plotMean(int surf, int chan, int bin, int calInt=0, string fileName=fileNameGlob) {

  string calType;
  if (calInt==0) {
    calType = "None"; }
  else if (calInt==1) {
    calType = "Little"; }
  else if (calInt==2) {
    calType = "Big"; }
  else {
    cout << "ERROR in plotMean(): " << calInt << " not a valid cal type!" << endl; }

  TFile *inFile = TFile::Open(fileName.c_str());

  stringstream name;

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_" << calType << "_lin";
  TH2D *hBig = (TH2D*)inFile->Get(name.str().c_str());

  TH1D *projection = hBig->ProjectionY("projection",bin+1,bin+1);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);

  
    //  TF1 *f1 = new TF1("fitFunction","TMath::Log10( (1./(TMath::Sqrt(2*TMath::Pi()*[0]^2))) * TMath::Exp(-1*(x-[1])^2/(2*[0]^2)))",0,100);
  //  f1->Draw("same");
  TF1 *riceFit = makeRiceFit();
  projection->Fit(riceFit);
  //  cout << "peak: " << riceFit->GetMaximumX() << endl;
  projection->Draw();
  riceFit->Draw("same");
  
  


  return;

}


 
void plotRice() {
 
  const int numLines=10;
  TF1 *f1[numLines];
  //  f1[0] = new TF1("fitFunction","[0]^2 * ([1]-x) * TMath::Exp( -1*(([1]-x)^2 + [2]^2) / (2*[0]^2)) * TMath::BesselI0(([1]-x)*[2]/[0]^2)",-10,75);
  //  f1[0]->SetParameters(5,75,0);

  //  f1[0] = new TF1("fitFunction","TMath::Log10( (1./(TMath::Sqrt(2*TMath::Pi()*[0]^2))) * TMath::Exp(-1*(x-[1])^2/(2*[0]^2)))",0,100);
  //  f1[0]->SetParameters(10,79);

  f1[0]->Draw("");
  
  for (int line=1; line<numLines; line++) {
    f1[line] = (TF1*)f1[0]->Clone();
    f1[line]->SetParameters(10,79+line);
    f1[line]->SetLineColor(line);
    f1[line]->Draw("same");
  }

  return;


}


void saveAll() {

  

  for (int surf=0; surf<12; surf++) {
    for (int chan=0; chan<8; chan++) {
      yFactor(surf,chan,fileNameGlob,true,true);
    }
  }


  return;

}



void compAll() {


  TGraph *gGain[96];
  TGraph *gNoiseT[96];

  TH2D *hGain;
  TH2D *hNoiseT;

  TH2D *hGain2;
  TH2D *hNoiseT2;
   
  for (int surf=0; surf<12; surf++) {
    for (int chan=0; chan<8; chan++) {
      int chanIndex = surf*8 + chan;
      cout << "surf" << surf << " chan" << chan << " chanIndex:" << chanIndex << endl;
      gGain[chanIndex] = new TGraph();
      gNoiseT[chanIndex] = new TGraph();
      TGraph *gGainCurr = gGain[chanIndex];
      TGraph *gNoiseTCurr = gNoiseT[chanIndex];
      yFactor(surf,chan,fileNameGlob,false,false,gGainCurr,gNoiseTCurr);

      int lenGGain = gGainCurr->GetN();
      if (surf==0 && chan==0) {
	hGain = new TH2D("hGain","ANITA3 gain from y-factor; Frequency (MHz); chanIndex (surf*8 + chan); Gain (dB)",
			 lenGGain,0,TMath::MaxElement(lenGGain,gGainCurr->GetX()), 96,-0.5,95.5);
	hNoiseT = new TH2D("hNoiseT","ANITA3 noise figure from y-factor; Frequency (MHz); chanIndex (surf*8 + chan); Noise Temperature (K)",
			   lenGGain,0,TMath::MaxElement(lenGGain,gNoiseTCurr->GetX()), 96,-0.5,95.5);

	hGain2 = new TH2D("hGain2","ANITA3 gain from y-factor; Frequency (MHz); gain (dB)",
			 lenGGain,0,TMath::MaxElement(lenGGain,gGainCurr->GetX()), 30,40,70);

	hNoiseT2 = new TH2D("hNoiseT2","ANITA3 noise figure from y-factor; Frequency (MHz); Noise Temperature (K)",
			   lenGGain,0,TMath::MaxElement(lenGGain,gNoiseTCurr->GetX()),50,0,500);

      }

      for (int pt=0; pt<lenGGain; pt++) {
	hGain->Fill(gGainCurr->GetX()[pt],surf*8+chan,gGainCurr->GetY()[pt]);
	hNoiseT->Fill(gNoiseTCurr->GetX()[pt],surf*8+chan,gNoiseTCurr->GetY()[pt]);
	hGain2->Fill(gGainCurr->GetX()[pt],gGainCurr->GetY()[pt]);
	hNoiseT2->Fill(gNoiseTCurr->GetX()[pt],gNoiseTCurr->GetY()[pt]);
      }
      
    }
  }


  TCanvas *c1 = new TCanvas("cCompAll","cCompAll",1000,600);
  c1->Divide(2);
  c1->cd(1);
  hGain->Draw("colz");
  c1->cd(2);
  hNoiseT->Draw("colz");


  TCanvas *c2 = new TCanvas("cCompAll2","cCompAll2",1000,600);
  c2->Divide(2);
  for (int chanIndex=0; chanIndex<96; chanIndex++) {
    if (chanIndex==0) {
      c2->cd(1);      
      gGain[chanIndex]->Draw("alp");
      cout << "farts" << endl;
      c2->cd(2);
      gNoiseT[chanIndex]->Draw("alp");
      cout << "poops" << endl;
    }
    else {
      c2->cd(1);
      gGain[chanIndex]->Draw("lpSame");
      c2->cd(2);
      gNoiseT[chanIndex]->Draw("lpSame");
    }

  }

  TCanvas *c3 = new TCanvas("cCompAll3","cCompAll3",1000,600);
  c3->Divide(2);
  c3->cd(1);
  hGain2->Draw("colz");
  c3->cd(2);
  hNoiseT2->Draw("colz");


  return;

}
