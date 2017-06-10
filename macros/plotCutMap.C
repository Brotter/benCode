

void plotCutMap(string date = "04.09.17_20h") {


  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  
  gROOT->ProcessLine(".x setupProof.C");

  stringstream name;
  for (int run=320; run<354; run++) {
    //    if ( (run==130) || (run==144) || (run==150) || (run==186) || (run==198) )continue;
    name.str("");
    //    name << "/home/brotter/nfsShared/results/templateSearch/" << date << "/" << run << ".root";
    name << "/Volumes/ANITA3data/bigAnalysisFiles/templateSearch/" << date << "/" << run << ".root";
    summaryTree->Add(name.str().c_str());

  }

  AnitaEventSummary *summary = NULL;
  summaryTree->SetBranchAddress("eventSummary",&summary);
  Double_t templateValue;
  summaryTree->SetBranchAddress("templateValue",&templateValue);

  TH2D *noiseH_HpolTrig    = new TH2D("noiseH_HpolTrig","noiseH_HpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *noiseV_HpolTrig    = new TH2D("noiseV_HpolTrig","noiseV_HpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *noiseH_VpolTrig    = new TH2D("noiseH_VpolTrig","noiseH_VpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *noiseV_VpolTrig    = new TH2D("noiseV_VpolTrig","noiseV_VpolTrig",500,0,0.5,1000,0,10e-6);
	                   
  TH2D *sigH_HpolTrig      = new TH2D("sigH_HpolTrig","sigH_HpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *sigV_HpolTrig      = new TH2D("sigV_HpolTrig","sigV_HpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *sigH_VpolTrig      = new TH2D("sigH_VpolTrig","sigH_VpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *sigV_VpolTrig      = new TH2D("sigV_VpolTrig","sigV_VpolTrig",500,0,0.5,1000,0,10e-6);
	                   
  TH2D *noiseH_BothpolTrig = new TH2D("noiseH_BothpolTrig","noiseH_BothpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *noiseV_BothpolTrig = new TH2D("noiseV_BothpolTrig","noiseV_BothpolTrig",500,0,0.5,1000,0,10e-6);
  TH2D *sigH_BothpolTrig   = new TH2D("sigH_BothpolTrig","sigH_BothpolTrig"    ,500,0,0.5,1000,0,10e-6);
  TH2D *sigV_BothpolTrig   = new TH2D("sigV_BothpolTrig","sigV_BothpolTrig"    ,500,0,0.5,1000,0,10e-6);

  int lenEntries = summaryTree->GetEntries();
  for (int entry=0; entry<lenEntries; entry++) {
    summaryTree->GetEntry(entry);
    if (entry%10000 == 0) cout << entry << "/" << lenEntries << "\r"; fflush(stdout);

    double hTemplate = templateValue / summary->coherent[0][0].totalPower;
    double vTemplate = templateValue / summary->coherent[1][0].totalPower;

    double hPeak = summary->peak[0][0].value;
    double vPeak = summary->peak[1][0].value;

    int pulser = summary->flags.pulser;

    int isHPolTrigger = summary->flags.isHPolTrigger;
    int isVPolTrigger = summary->flags.isVPolTrigger;

    //H only
    if ( (isHPolTrigger == 1) && (isVPolTrigger == 0) ) {
      if (pulser ==0) {
	noiseH_HpolTrig->Fill(hPeak,hTemplate);
	noiseV_HpolTrig->Fill(vPeak,vTemplate);
	  }
      else {
	sigH_HpolTrig->Fill(hPeak,hTemplate);
	sigV_HpolTrig->Fill(vPeak,vTemplate);
	  }
    }

    //V only
    if ( (isHPolTrigger == 0) && (isVPolTrigger == 1) ) {
      if (pulser ==0) {
	noiseH_VpolTrig->Fill(hPeak,hTemplate);
	noiseV_VpolTrig->Fill(vPeak,vTemplate);
	  }
      else {
	sigH_VpolTrig->Fill(hPeak,hTemplate);
	sigV_VpolTrig->Fill(vPeak,vTemplate);
	  }
    }
    

    //Both
    if ( (isHPolTrigger == 1) && (isVPolTrigger == 1) ) {
      if (pulser ==0) {
	noiseH_BothpolTrig->Fill(hPeak,hTemplate);
	noiseV_BothpolTrig->Fill(vPeak,vTemplate);
	  }
      else {
	sigH_BothpolTrig->Fill(hPeak,hTemplate);
	sigV_BothpolTrig->Fill(vPeak,vTemplate);
	  }
    }
    
  }


  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,2);
  c1->cd(1);
  noiseH_HpolTrig->Draw("colz");
  c1->cd(2);
  sigH_HpolTrig->Draw("colz");
  c1->cd(3);
  noiseH_VpolTrig->Draw("colz");
  c1->cd(4);
  sigH_VpolTrig->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->Divide(2,2);
  c2->cd(1);
  noiseV_HpolTrig->Draw("colz");
  c2->cd(2);
  sigV_HpolTrig->Draw("colz");
  c2->cd(3);
  noiseV_VpolTrig->Draw("colz");
  c2->cd(4);
  sigV_VpolTrig->Draw("colz");


  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->Divide(2,2);
  c3->cd(1);
  noiseH_BothpolTrig->Draw("colz");
  c3->cd(2);
  noiseV_BothpolTrig->Draw("colz");
  c3->cd(3);
  sigH_BothpolTrig->Draw("colz");
  c3->cd(4);
  sigV_BothpolTrig->Draw("colz");
    

  TFile *outFile = TFile::Open("plotCutMap.root","recreate");
  noiseH_HpolTrig->Write();    
  noiseV_HpolTrig->Write();    
  noiseH_VpolTrig->Write();    
  noiseV_VpolTrig->Write();    
                   
  sigH_HpolTrig->Write();      
  sigV_HpolTrig->Write();      
  sigH_VpolTrig->Write();      
  sigV_VpolTrig->Write();      
                   
  noiseH_BothpolTrig->Write(); 
  noiseV_BothpolTrig->Write(); 
  sigH_BothpolTrig->Write();   
  sigV_BothpolTrig->Write();   
  outFile->Close();



return;
}


void drawFromTree(TChain* summaryTree){

  summaryTree->SetProof();

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  TCanvas *c2_V = new TCanvas("c2_V","c2_V",1000,600);
  TCanvas *c2_H = new TCanvas("c2_H","c2_H",1000,600);
  TCanvas *c4 = new TCanvas("c4","c4",1000,600);
  TCanvas *c3 = new TCanvas("c3","c3",1000,600);


  c1->cd();
  summaryTree->Draw("templateValue/coherent[1][0]->totalPower:peak[1][0]->value",
		    "flags->isHPolTrigger == 0 && flags->isVPolTrigger == 1 && flags->pulser == 0","colz");

  c2->cd();
  summaryTree->Draw("templateValue/coherent[0][0]->totalPower:peak[0][0]->value",
		    "flags->isHPolTrigger == 1 && flags->isVPolTrigger == 0 && flags->pulser == 0","colz");


  c2_V->cd();
  summaryTree->Draw("templateValue/coherent[0][0]->totalPower:peak[0][0]->value",
		    "flags->isHPolTrigger == 1 && flags->isVPolTrigger == 1 && flags->pulser == 0","colz");

  c2_H->cd();
  summaryTree->Draw("templateValue/coherent[1][0]->totalPower:peak[1][0]->value",
		    "flags->isHPolTrigger == 1 && flags->isVPolTrigger == 1 && flags->pulser == 0","colz");


  c3->cd();
  summaryTree->Draw("templateValue/coherent[0][0]->totalPower:peak[0][0]->value",
		    "flags->isHPolTrigger == 1 && flags->isVPolTrigger == 0 && flags->pulser != 0","colz");


  c4->cd();
  summaryTree->Draw("templateValue/coherent[1][0]->totalPower:peak[1][0]->value",
		    "flags->isHPolTrigger == 0 && flags->isVPolTrigger == 1 && flags->pulser != 0","colz");
  
 

  return;

}
