#include "BedmapReader.h"

void bedMapPlot(){

  BedmapReader *bedmap = new BedmapReader(0);

  const double latLow = -82;
  const double latHigh = -73;
  const double lonLow = -140;
  const double lonHigh = -90;

  const int numLatBins = 101;
  double latPerBin = (latHigh-latLow)/numLatBins;
  const int numLonBins = 101;
  double lonPerBin = (lonHigh-lonLow)/numLonBins;

  cout <<  latPerBin << " " << lonPerBin << endl;

  TProfile2D *hist = new TProfile2D("bedmap","bedmap;lon;lat",numLonBins-1,lonLow,lonHigh,numLatBins-1,latLow,latHigh);

  for (int latBin=0; latBin<numLatBins; latBin++) {
    double lat = latLow+latBin*latPerBin;

    for (int lonBin=0; lonBin<numLatBins; lonBin++) {
      double lon = lonLow+lonBin*lonPerBin;

      double value = bedmap->Surface(lon,lat);
      cout << value << endl;
      hist->Fill(lon,lat,value);
    }
  }

  hist->GetZaxis()->SetRangeUser(6374e3,6380e3);
  hist->Draw("colz");


}
    
  
