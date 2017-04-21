/*


  Plotting and analysis tools to use on the output file from yFactorFreqDomain


 */


TGraph *histToMeanGraph(TH2D* inGraph) {
  
  TGraph *outGraph = new TGraph();
  int numBins = inGraph->GetNbinsX();

  for (int bin=0; bin<numBins; bin++) {
    TH1D *temp = inGraph->ProjectionY("temp",bin+1,bin+1);
    double mean = temp->GetMean();
    delete temp;
    outGraph->SetPoint(bin,inGraph->GetXaxis()->GetBinCenter(bin+1),mean);
  }

  return outGraph;

}


void plotSurfChan(int surf, int chan, string fileName="yf_rms.root") {

  TFile *inFile = TFile::Open(fileName.c_str());

  stringstream name;

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Big";
  TH2D *hBig = (TH2D*)inFile->Get(name.str().c_str());
  hBig->SetMarkerColor(kRed);
  hBig->GetYaxis()->SetRangeUser(0,80);
  hBig->SetStats(0);
  

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Little";
  TH2D *hLittle = (TH2D*)inFile->Get(name.str().c_str());
  hLittle->SetMarkerColor(kBlue);
  hLittle->GetYaxis()->SetRangeUser(0,80);
  hLittle->SetStats(0);

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_None";
  TH2D *hNone = (TH2D*)inFile->Get(name.str().c_str());
  hNone->GetYaxis()->SetRangeUser(0,80);
  hNone->SetStats(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(3);
  for (int i=0; i<3; i++) {
    c1->GetPad(i+1)->SetGridx(1);
    c1->GetPad(i+1)->SetGridy(1);
  }

  c1->cd(1);
  hNone->Draw("colz");
  TGraph *gNone = histToMeanGraph(hNone);
  gNone->Draw("lpSame");
  c1->cd(2);
  hLittle->Draw("colz");
  TGraph *gLittle = histToMeanGraph(hLittle);
  gLittle->Draw("lpSame");
  c1->cd(3);
  hBig->Draw("colz");
  TGraph *gBig = histToMeanGraph(hBig);
  gBig->Draw("lpSame");


  return;

}







/*
  Test stuff
 */

void plotMean(int surf, int chan, string fileName="yf_all.root") {


  TFile *inFile = TFile::Open(fileName.c_str());

  stringstream name;

  name.str("");
  name << "surf" << surf << "Chan" << chan << "_Big_lin";
  TH2D *hBig = (TH2D*)inFile->Get(name.str().c_str());

  TH1D *projection = hBig->ProjectionY("projection",30,30);

  TGraph *projUnlog = new TGraph();
  projUnlog->SetName("projUnlog");

  for (int bin=0; bin<projection->GetNbinsX(); bin++) {
    double xValue = TMath::Sqrt(pow(10.,projection->GetXaxis()->GetBinCenter(bin+1)/10.))/230;
    double yValue = projection->GetBinContent(bin+1);
    cout << bin << " " << xValue << " " << yValue << endl;
    projUnlog->SetPoint(bin,xValue,yValue);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2);
  c1->cd(1);

  TF1 *f1 = new TF1("fitFunction","[3]^2 * ([1]-x) * TMath::Exp( -1*(([1]-x)^2 + [2]^2) / (2*[0]^2)) * TMath::BesselI0(([1]-x)*[2]/[0]^2)",0,100);
    //  TF1 *f1 = new TF1("fitFunction","TMath::Log10( (1./(TMath::Sqrt(2*TMath::Pi()*[0]^2))) * TMath::Exp(-1*(x-[1])^2/(2*[0]^2)))",0,100);
  f1->SetParameters(6.8,79,2.3,4.6);
  f1->Draw("same");
  projection->Fit(f1);
  f1->Draw("same");
  projection->Draw();

  c1->cd(2);
  TF1 *f2 = new TF1("fitFunction2","[3]^2 * (x*[1]) * TMath::Exp( -1*((x*[1])^2 + [2]^2) / (2*[0]^2)) * TMath::BesselI0((x*[10])*[2]/[0]^2)",0,100);
  f2->SetParameters(20,1,0,2.67);
  projUnlog->Fit(f2);
  projUnlog->Draw("alp");
  f2->Draw("same");
  


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
