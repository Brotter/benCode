/*


  How about instead of making a billion copies, I include a string entry with each event that describes it


 */









void labelEvents(string inFileName = "cutsClust_oct14.root") {
  /*

    Reads in a file, outputs a file with "_labeled" tacked on before the .root

    gotta give it the file with clusterValue included though (from cluster.C)

    Current labels added:
    - Below Horizon CR Candidate
    - Above Horizon CR Candidate
    - Inverted Candidate

    - Passing Above Horizon
    - Failing Above Horizon

    - Dirty Dozen (failing isolated)

    - Failing Clustered
    - Passing Clustered


   */


  TChain *summaryTree = new TChain("summaryTree","summaryTree");
  summaryTree->Add(inFileName.c_str());

  int lenEntries = summaryTree->GetEntries();
  cout << " Found " << lenEntries << " events" << endl;


  size_t pos = inFileName.find(".root");
  string baseName = inFileName.substr(0,pos);
  stringstream outFileName;
  outFileName << baseName << "_labeled.root";

  TFile *outFile = TFile::Open(outFileName.str().c_str(),"recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  outTree->Branch("gpsEvent",&gps);
  double clusterValue;
  summaryTree->SetBranchAddress("clusterValue",&clusterValue);
  outTree->Branch("clusterValue",&clusterValue);
  TString *labelString = NULL;
  outTree->Branch("label",&labelString);
  labelString = new TString();

  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%10)) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);

    labelString->Clear();

    // Some basic info for splitting them up
    bool fClustered = false; //is it clustered? (L < 40)
    if (clusterValue < 40) {fClustered = true;}
    
    bool fStrongCuts = false; //does it pass strong cuts?
    //nominal cuts that I initially developed.  Based on WAIS pulser signals.
    //see presentation on August 16th
    //these are the final signal cuts, combined with the clustering results from the weak cuts
    if ( (evSum->peak[0][0].value > 0.0435)
	 && (evSum->peak[0][0].snr > 9.05)
	 && (evSum->coherent_filtered[0][0].peakHilbert > 31.1)
	 && (evSum->coherent_filtered[0][0].linearPolFrac() > 0.60)
	 && (tempSum->coherent[0][0].cRay[4] > 0.666) 
	 && (tempSum->deconvolved[0][0].cRay[4] > 0.666) ) { fStrongCuts = true; }
    

    bool fAboveHorizon = false; //is it above the horizon?
    if (evSum->peak[0][0].altitude <= -999 || evSum->peak[0][0].theta_adjustment_needed > 0) { fAboveHorizon = true; }
    
    
    //if it is above the horizon clustering is stupid
    if (fAboveHorizon) {

      //if it is strong, it is either a candidate or clustered
      if (fStrongCuts) {
	//"Above Horizon Candidate" : these two
	if (evSum->eventNumber == 39599205 || evSum->eventNumber == 27142546) {
	  labelString->Append("Above Horizon Candidiate");
	}
	//"Above Horizon Passing" : the rest that passed
	else {
	  labelString->Append("Above Horizon Passing");
	}
      }
      //"Above Horizon Failing" : if it is weak then it is just a failing above horizon
      else {
	labelString->Append("Above Horizon Failing");
      }
    }

    //otherwise it is below horizon!
    else {


      // "Passing Clustered" : clustered, but pass signal cuts (~800)
      if (fClustered && fStrongCuts) {
	labelString->Append("Passing Clustered");
      }

      // "Failing Clustered" : clustered and weak
      if (fClustered && !fStrongCuts) {
	labelString->Append("Failing Clustered");
      }


      // "Dirty Dozen" : isolated, but wak
      else if (!fClustered && !fStrongCuts) {
	labelString->Append("Dirty Dozen");
      }


      // Candidates!  Not clustered, passes strong cuts
      else if (!fClustered && fStrongCuts) {
	
	// "Inverted Candidate" : try to ignore this one
	if (evSum->eventNumber == 15717147) {
	  labelString->Append("Inverted Candidate");
	}
	
	// "Below Horizon Candidate" : whatevers left I guess
	else {
	  labelString->Append("Below Horizon Candidate");
	}
      }
    }
    
    outTree->Fill();

  } 

  
  cout << "Saving..." << endl;
  outFile->cd();
  outTree->Write();
  outFile->Close();



  cout << "Done!" << endl;


  return labelEvents;

}
