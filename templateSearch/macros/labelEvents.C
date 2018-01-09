/*


  How about instead of making a billion copies, I include a string entry with each event that describes it


 */

#include "AnitaEventSummary.h"
#include "loadAll.C"





void labelEvents() {
  /*

    Reads in a file, outputs a file with "_labeled" tacked on before the .root

    gotta give it the file with clusterValue included though (from cluster.C)

    This is the complex version.  It wants to read _all_ the quality data.
     It will label the "backgroundClustered" stuff too
     So this means every possible events gets a category



    Current labels added:
    - "Below Horizon Candidate" (isolated passing)
    - "Inverted Candidate"
    - "Dirty Dozen" (aka isolated failing)
    - "Above Horizon Candidate"
    - "Above Horizon Passing"
    - "Above Horizon Failing"
    - "Above Horizon Not Impulsive" - not sure if there will be any of these
    - "Clustered Passing" : clustered with another "impulsive" event
    - "Clustered Failing"

    - "Clustered Not Impulsive" : failing impulsivity cuts, but close to something impulsive
    - int seedEventNumber : the event number that it clusters with (-1=wais, -2=ldb, -999= not clustered)
    - double clusterValue : distance from nearest impulsive event ( -999 if not clustered )
    - "The Wasteland" : not clustered, not impulsive.  Nothing.

    - bool fMajorBase : whether it comes from within 3 sigma of a "major" base.  All pulsers should have this flag
          (sigma= 0.75:0.3 | phi:theta)
   */

  //input tree, reKey is the last set of "quality" data
  TChain *summaryTree = loadReKey(false);
  int lenEntries = summaryTree->GetEntries();
  cout << " Found " << lenEntries << " events" << endl;

  //load the impulsive vs impulsive cluster so you know what is "isolated"
  // 5997 entries, all should be impulsive...
  TChain *clusterTree = new TChain("summaryTree");
  clusterTree->Add("cutsClust_oct14.root");

  clusterTree->BuildIndex("eventNumber");

  //load the last backgroundCluster list of events that cluster with Impulsives
  TChain *backgroundTree = loadWhatever("/Volumes/ANITA3Data/bigAnalysisFiles/cluster/10.31.17_11h59m/clusterBackground",64,false);
  backgroundTree->BuildIndex("eventNumber");
  

  //make a save file
  //  size_t pos = inFileName.find(".root");
  //  string baseName = inFileName.substr(0,pos);
  //  stringstream outFileName;
  //  outFileName << baseName << "_labeled.root";

  TFile *outFile = TFile::Open("labelEvents_All.root","recreate");
  TTree *outTree = new TTree("summaryTree","summaryTree");

  AnitaEventSummary *evSum = NULL;
  AnitaTemplateSummary *tempSum = NULL;
  AnitaNoiseSummary *noiseSum = NULL;
  Adu5Pat *gps = NULL;
  summaryTree->SetBranchAddress("eventSummary",&evSum);
  //  outTree->Branch("eventSummary",&evSum);
  summaryTree->SetBranchAddress("template",&tempSum);
  //  outTree->Branch("template",&tempSum);
  summaryTree->SetBranchAddress("noiseSummary",&noiseSum);
  //  outTree->Branch("noiseSummary",&noiseSum);
  summaryTree->SetBranchAddress("gpsEvent",&gps);
  //  outTree->Branch("gpsEvent",&gps);
  double backClusterValue;
  backgroundTree->SetBranchAddress("clusterValue",&backClusterValue);
  outTree->Branch("backClusterValue",&backClusterValue);
  double impClusterValue;
  clusterTree->SetBranchAddress("clusterValue",&impClusterValue);
  outTree->Branch("impClusterValue",&impClusterValue);

  //  double exGeoPol;
  //  summaryTree->SetBranchAddress("exGeoPol",&exGeoPol);
  //  outTree->Branch("exGeoPol",&exGeoPol);
  int seedEventNumber;
  backgroundTree->SetBranchAddress("seedEventNumber",&seedEventNumber);
  outTree->Branch("seedEventNumber",&seedEventNumber);

  bool fInBackgroundFile,fInClusterFile;
  outTree->Branch("fInBackgroundFile",&fInBackgroundFile);
  outTree->Branch("fInClusterFile",&fInClusterFile);
  
  bool fMajorBase;
  outTree->Branch("fMajorBase",&fMajorBase);
  
  //this file is disconnected, so add the event number
  int eventNumber;
  outTree->Branch("eventNumber",&eventNumber);

  // label!
  TString *labelString = NULL;
  outTree->Branch("label",&labelString);
  labelString = new TString();


  TStopwatch watch;
  watch.Start(kTRUE);
  int totalTimeSec = 0;
  int savedCount = 0;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10000 == 0 && entry>0) {
      int timeElapsed = watch.RealTime();
      totalTimeSec += timeElapsed;
      double rate = float(entry)/totalTimeSec;
      double remaining = (float(lenEntries-entry)/rate)/60.;
      watch.Start();
      cout << entry << "/" << lenEntries << " (" << savedCount << ") ";
      cout << 10000./timeElapsed << "Hz <" << rate << "> " << remaining << " minutes left, ";
      cout << totalTimeSec/60. << " minutes elapsed" << endl;
    }

    //clear the string
    labelString->Clear();
    
    //get the next entry
    summaryTree->GetEntry(entry);
    eventNumber = evSum->eventNumber;

    //find if it is in the background file
    int backIndex = backgroundTree->GetEntryNumberWithIndex(eventNumber);
    if (backIndex > 0) {
      fInBackgroundFile = true;
      backgroundTree->GetEntry(backIndex);
    }
    else {
      fInBackgroundFile = false;
      backClusterValue = -999;
      seedEventNumber = -999;
    }

    //find if it is in the cluster file
    int clustIndex = clusterTree->GetEntryNumberWithIndex(eventNumber);
    if (clustIndex > 0) {
      fInClusterFile = true;
      clusterTree->GetEntry(clustIndex);
    }
    else {
      impClusterValue = -999;
      fInClusterFile = false;
    }


    // Some basic info for splitting them up: 

    //is it from a major base? (sigma= 0.75:0.3 | phi:theta)
    fMajorBase = false;
    if ( (TMath::Sqrt( pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.phi,360,0))/0.75,2) + 
		       pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->wais.theta,360,0))/0.3,2) ) < 3.0 
	  && evSum->wais.distance < 1000e3 ) ||
	 (TMath::Sqrt( pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.phi,360,0))/0.75,2) + 
		       pow(TMath::Abs(FFTtools::wrap(evSum->peak[0][0].phi - evSum->ldb.theta,360,0))/0.3,2) ) < 3.0 
	  && evSum->ldb.distance < 1000e3 ) ) {
      fMajorBase = true;
    }

    //is it geoexcluded? (L < 40 from nearest impulsive event)
    // if it is impulsive or more, then check the cluster file
    // if it is not impulsive, then check the background file
    bool fGeoExcluded = false; //default no
    if (fInClusterFile) { //impulsive
      if (impClusterValue < 40) {fGeoExcluded = true; }
    }
    else if (fInBackgroundFile) { fGeoExcluded = true; }

    bool fImpulsive = false;
    // impulsive cuts, from separateNotable::savePassingEvents :
    //"weak cuts. These are for determining the anthropogenic base distribution"
    if ( (evSum->peak[0][0].value > 0.0435)
	 && (evSum->peak[0][0].snr > 8.95)
	 && (evSum->coherent_filtered[0][0].peakHilbert > 25)
	 && (evSum->coherent_filtered[0][0].linearPolFrac() > 0.5)
	 && (tempSum->coherent[0][0].cRay[4] > 0.5) )
	 //	 && (tempSum->deconvolved[0][0].cRay[4] > 0.5) ) //I don't think I used this one
      { fImpulsive = true; }
	 
	  
    bool fSignal = false; 
    //nominal cuts that I initially developed.  Based on WAIS pulser signals.
    //see presentation on August 16th
    //these are the final signal cuts, combined with the clustering results from the weak cuts
    if ( (evSum->peak[0][0].value > 0.0435)
	 && (evSum->peak[0][0].snr > 9.05)
	 && (evSum->coherent_filtered[0][0].peakHilbert > 31.1)
	 && (evSum->coherent_filtered[0][0].linearPolFrac() > 0.60)
	 && (tempSum->coherent[0][0].cRay[4] > 0.666) 
	 && (tempSum->deconvolved[0][0].cRay[4] > 0.666) ) 
      { fSignal = true; }


    bool fAboveHorizon = false; //is it above the horizon?
    if (evSum->peak[0][0].altitude <= -999 || evSum->peak[0][0].theta_adjustment_needed > 0) 
      { fAboveHorizon = true; }


    //checks
    if (!fInBackgroundFile && fInClusterFile) {
      cout << "ev" << eventNumber << " is in the cluster file but NOT the background file!?!?!" << endl;
      cout << "fAboveHorizon:" << fAboveHorizon << endl;
    }

    if (fSignal && !fImpulsive) {
      cout << "ev" << eventNumber << " is signal but NOT impulsive!?!?!" << endl;
    }
    
    
    //if it is above the horizon
    // clustering is stupid here, so it is only a strength based thing
    //  two events were determined by hand to not cluster and are the aboveHorizonCandidates
    if (fAboveHorizon) {

      //if it is strong, it is either a candidate or clustered
      if (fSignal) {
	//"Above Horizon Candidate" : these two
	if (evSum->eventNumber == 39599205 || evSum->eventNumber == 27142546) {
	  labelString->Append("Above Horizon Candidate");
	}
	//"Above Horizon Passing" : the rest that passed
	else {
	  labelString->Append("Above Horizon Passing");
	}
      }

      //"Above Horizon Failing" : if it is weak then it is just a failing above horizon
      else if (fImpulsive) {
	labelString->Append("Above Horizon Failing");
      }
      else {
	labelString->Append("Above Horizon Not Impulsive");
      }
    }



    //otherwise it is below horizon!
    // These things can cluster
    else {
      

      //Things that cluster with Impulsives
      if (fGeoExcluded) {

	// "Clustered Passing" : clustered, but pass signal cuts (~800)
	if (fSignal) {
	  labelString->Append("Clustered Passing");
	}
	// "Clustered Failing" : clustered and weak
	else if (fImpulsive) {
	  labelString->Append("Clustered Failing");
	}	
	// "Clustered Not Impulsive" 
	else {
	  labelString->Append("Clustered Not Impulsive");
	}
      }


      //Things that do _not_ cluster with impulsives
      else {

	// Candidates!  Not clustered, passes strong cuts
	if (fSignal) {
	
	  // "Inverted Candidate" : try to ignore this one, but just label it for now
	  if (evSum->eventNumber == 15717147) {
	    labelString->Append("Inverted Candidate");
	  }
	  
	  // "Below Horizon Candidate" : whatevers left I guess (aboveHorz are above)
	  else {
	    labelString->Append("Below Horizon Candidate");
	  }
	}

	// "Dirty Dozen" : isolated, but not signal cuts, just impulsive
	else if (fImpulsive) {
	  labelString->Append("Dirty Dozen");
	}
	//Isolated and weak, the wasteland
	else { 
	  labelString->Append("The Wasteland");
	  savedCount++;
	}
      }

    }//endelse below horizon  
      

    outTree->Fill();

  } 

  
  cout << "Saving..." << endl;
  outFile->cd();
  outTree->Write();
  outFile->Close();



  cout << "Done!" << endl;


  return labelEvents;

}








/********************

 *******************/




void labelEventsSimple(string inFileName = "cutsClust_oct14_geo.root") {
  /*

    Reads in a file, outputs a file with "_labeled" tacked on before the .root

    gotta give it the file with clusterValue included though (from cluster.C)

    Current labels added:
    - "Below Horizon Candidate"
    - "Above Horizon Candidate"
    - "Inverted Candidate"
    - "Dirty Dozen" (aka isolated failing)
    - "Above Horizon Passing"
    - "Above Horizon Failing"
    - "Clustered Failing"
    - "Clustered Passing"

    
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
  double exGeoPol;
  summaryTree->SetBranchAddress("exGeoPol",&exGeoPol);
  outTree->Branch("exGeoPol",&exGeoPol);
 
  TString *labelString = NULL;
  outTree->Branch("label",&labelString);
  labelString = new TString();

  for (int entry=0; entry<lenEntries; entry++) {
    if (!(entry%10)) cout << entry << "/" << lenEntries << endl;
    summaryTree->GetEntry(entry);

    labelString->Clear();

    // Some basic info for splitting them up
    bool fGeoExcluded = false; //is it clustered? (L < 40)
    if (clusterValue < 40) {fGeoExcluded = true;}
    
    bool fSignal = false; //does it pass strong cuts?
    //nominal cuts that I initially developed.  Based on WAIS pulser signals.
    //see presentation on August 16th
    //these are the final signal cuts, combined with the clustering results from the weak cuts
    if ( (evSum->peak[0][0].value > 0.0435)
	 && (evSum->peak[0][0].snr > 9.05)
	 && (evSum->coherent_filtered[0][0].peakHilbert > 31.1)
	 && (evSum->coherent_filtered[0][0].linearPolFrac() > 0.60)
	 && (tempSum->coherent[0][0].cRay[4] > 0.666) 
	 && (tempSum->deconvolved[0][0].cRay[4] > 0.666) ) { fSignal = true; }
    

    bool fAboveHorizon = false; //is it above the horizon?
    if (evSum->peak[0][0].altitude <= -999 || evSum->peak[0][0].theta_adjustment_needed > 0) { fAboveHorizon = true; }
    
    
    //if it is above the horizon clustering is stupid
    if (fAboveHorizon) {

      //if it is strong, it is either a candidate or clustered
      if (fSignal) {
	//"Above Horizon Candidate" : these two
	if (evSum->eventNumber == 39599205 || evSum->eventNumber == 27142546) {
	  labelString->Append("Above Horizon Candidate");
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


      // "Clustered Passing" : clustered, but pass signal cuts (~800)
      if (fGeoExcluded && fSignal) {
	labelString->Append("Clustered Passing");
      }

      // "Clustered Failing" : clustered and weak
      if (fGeoExcluded && !fSignal) {
	labelString->Append("Clustered Failing");
      }


      // "Dirty Dozen" : isolated, but wak
      else if (!fGeoExcluded && !fSignal) {
	labelString->Append("Dirty Dozen");
      }


      // Candidates!  Not clustered, passes strong cuts
      else if (!fGeoExcluded && fSignal) {
	
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
