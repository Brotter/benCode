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
#include "TH2D.h"
//anita
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "Adu5Pat.h"
#include "CalibratedAnitaEvent.h"

using namespace std;

/*
  Ben Rotter September 2016

  This finds the WAIS pulses from the ANITA3 flight and creates a root file with all their headers.
  With that I can easily parse through them and plot things and determine impulse response characteristics,
  as well as compare the thermal noise environment to real signal events.

  I did this before, but now I am doing it again with the standardized functions (so I'm not at fault for
  calculating the distance wrong (though I did it right).

  It also uses the gpsEvent.root files that Cosmin made and are super useful!
  http://www.phys.hawaii.edu/~gorham/ANITA/GPSEventFiles011916.pdf

 */


int main(int argc, char** argv) {

  cout << "Hello!  Let us do some physics mate" << endl;

  stringstream name;
  //okay lets start by grabbing the header files
  TChain *headTree = new TChain("headTree","headTree");  
  TChain *eventTree = new TChain("eventTree","eventTree");
  //I also need the gps files to know where the heck ANITA is!
  TChain *patTree = new TChain("adu5PatTree","adu5PatTree");  
  for (int runNum=330; runNum<331; runNum++) {
    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << runNum << "/headFile" << runNum << ".root";
    headTree->Add(name.str().c_str());

    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << runNum << "/calEventFile" << runNum << ".root";
    eventTree->Add(name.str().c_str());

    name.str("");
    name << "/Volumes/ANITA3Data/root/run" << runNum << "/gpsEvent" << runNum << ".root";
    patTree->Add(name.str().c_str());
  }

  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);
  
  CalibratedAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);

  Adu5Pat *pat = NULL;
  patTree->SetBranchAddress("pat",&pat);

  int numEntries = headTree->GetEntries();

  cout << "I found " << numEntries << " header entries" << endl;
  cout << "I found " << eventTree->GetEntries() << " event entries" << endl;
  cout << "I found " << patTree->GetEntries() << " gps entries" << endl;
  //this is annoying cuz I gotta sort it (maybe not...)
  //  patTree->BuildIndex("eventNumber");
  

  //create the storage graphs
  TH2D *trigTimeNs = new TH2D("trigTimeNs","trigTimeNs",1000,0,numEntries,1000,1e6,2e6);
  TH2D *nsOffset = new TH2D("nsOffset","nsOffset",1000,0,numEntries,1000,-1500,600);
  TGraph *gNsExpect = new TGraph();
  gNsExpect->SetName("gNsExpect");

  //also make a Tree of event headers that pass cuts
  TFile *outHeadFile = TFile::Open("headFileWais.root","recreate");
  TTree *outHeadTree = new TTree("headTree","headTree");
  outHeadTree->Branch("header","RawAnitaHeader",&head);

  TFile *outEventFile = TFile::Open("calEventFileWais.root","recreate");
  TTree *outEventTree = new TTree("eventTree","eventTree");
  outEventTree->Branch("event","CalibratedAnitaEvent",&event);


 
  //loop through entries
  for (int entry=0; entry<numEntries; entry++) {
    if (entry%1000==0) {
      cout << entry << "/" << numEntries << "\r";
      fflush(stdout);
    }
    //get all the pointers set right
    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);
    patTree->GetEntry(entry);
    //Then find the info we can cut on
    UsefulAdu5Pat *usefulPat = new UsefulAdu5Pat(pat);
    long waisTriggerTimeNs = usefulPat->getWaisDivideTriggerTimeNs();
    long expectDiff = head->triggerTimeNs - waisTriggerTimeNs;
    delete usefulPat;

    ///////////////////////////////////////////////////////////
    //           CUTS(this gets them all, lol so easy!)      //
    if ((head->trigType&0x0F) != 1) continue;
    //The ground GPS loses sync a lot I think, though they are all in a 2.5uS band (.00025% of flight)
    if ((expectDiff > 1000) || (expectDiff < -1500)) continue; 
    ///////////////////////////////////////////////////////////

    //fill up the graphs and histograms
    trigTimeNs->Fill(entry,head->triggerTimeNs);
    nsOffset->Fill(entry,expectDiff);
    if (entry%1000==0) gNsExpect->SetPoint(gNsExpect->GetN(),entry,waisTriggerTimeNs);
    
    //write that header to the tree
    outHeadFile->cd();
    outHeadTree->Fill();
    outEventFile->cd();
    outEventTree->Fill();
    
  }
  //the header is full so close it up
  outHeadFile->cd();
  outHeadTree->Write();
  outHeadFile->Close();
  outEventFile->cd();
  outEventTree->Write();
  outEventFile->Close();

  //open up a file to write the graphs, write them, then close it
  TFile *output = TFile::Open("findWaisPulses.root","recreate");
  trigTimeNs->Write();
  nsOffset->Write();
  gNsExpect->Write();
  output->Close();


  cout << "Physics complete!  See ya later buddy :)" << endl;

  return 1;
  
}


  
