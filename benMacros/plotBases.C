{


  cout<<"initializing base list..."<<endl;
  std::string nameString;
  double input_lat;
  double input_lon;
  char lat_sign; //A character: this one should always be 'S'.
  char lon_sign; //A character: either 'E' or 'W'
  int misc_field;
  int ctr=0;
  ifstream base_file;
  char filename[150];
  char *pEnv;

  TGraph *bases = new TGraph();
  bases->SetTitle("antarcticBases");
  bases->SetName("antarcticBases");

  pEnv = getenv("ANITA_ANALYSIS_AG_BASE");
  sprintf(filename,"%s/baseLocations/all_base_locations_new.txt",pEnv);
  base_file.open(filename);
  if (!base_file.is_open())
    {
      std::cerr<<"Error!  Could not open "
	       <<("all_base_locations_new.txt")<<" to read in base locations!\n";
    } //end if
  while (base_file >> nameString >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {

      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      
      cout << "lat:" << input_lat << " lon:" << input_lon << endl;
      bases->SetPoint(bases->GetN(),input_lon,input_lat);
      ctr++;
    }
  
  //now do pseudo bases

  TGraph *pbases = new TGraph();
  pbases->SetTitle("pseudoBases");
  pbases->SetName("pseudoBases");


  char filenamePseudo[150];
    sprintf(filenamePseudo,"%s/baseLocations/pseudoBases.txt",pEnv);
  //std::string filenamePseudo="/rh5stuff/64bit/src/anita/analysis/agoodhue/baseLocations/pseudoBases.txt";
  ifstream base_filePseudo;
  base_filePseudo.open(filenamePseudo);
  if (!base_filePseudo.is_open())
    {
      std::cerr<<"Error!  Could not open "
	       <<("pseudoBases.txt")<<" to read in base locations!\n";
    } //end if
  while (base_filePseudo >> nameString >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {
      
      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      pbases->SetPoint(pbases->GetN(),input_lon,input_lat);
      ctr++;
    }

  //now there is room here to add my own hard-coded bases or traverses like stephen did
  cout<<"done initializing base list"<<endl;

}
