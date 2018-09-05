double CalculateX0(double theta){
	double X0min = 0.00001;
	double X0max = 10;  
	double tolerance = 1e-8;
	double X0 = (X0min+X0max)/2.0;	
	double error = X0max;

	double Ze;		// atomic number of target material
	Ze=par[1];

	double E = par[0];		// Energy in GeV
	double mass = par[2];	// mass of the particle in GeV
	double p = TMath::Sqrt(E*E-mass*mass);
	double beta = p/E;
	//cout<<"beta = "<<beta<<endl;
	//cout<<"p = "<<p<<endl;
	for (int i = 1; TMath::Abs(error) > tolerance; i++){
		X0 = (X0min+X0max)/2.0;
		error = theta/1000 - (13.6*0.001/(beta*p)*Ze*TMath::Sqrt(X0)*(1+0.038*TMath::Log(X0))); 	
		//cout<<"error = "<<error<<endl; 
		if (error <= 0){
			X0max = X0;
		}else{
			X0min = X0;
		}		
		if(i > 20){		//maximum  estimation accuracy for 30 iterations is 10^-8    
			break;
			cout<<"Solution for X/X0 did not converge"<<endl; 
		}
	}
	return X0;
}

TString FormFilename(TString path, int run){	
		TString fileStart = "run000";
		TString fileEnd = "-gblkinkestimator.root";
		TString out;
		if(run < 10){
			out = path + fileStart + Form("00%d",run) + fileEnd;
		}else if (run >= 10 && run < 100){
			out = path + fileStart + Form("0%d",run) + fileEnd;
		}else{
			out = path + fileStart + Form("%d",run) + fileEnd;
		}		
		return out;
	}
	
void SetColorPalette(){
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  // Viridis
  Double_t red[9]   = { 246./255., 144./255.,  74./255.,  35./255.,  28./255.,  33./255.,  43./255., 51./255., 26./255.};
  Double_t green[9] = {  222./255., 200./255.,  180./255.,  150./255., 118./255., 87./255., 55./255., 24./255., 9./255.};
  Double_t blue[9]  = { 0./255., 35./255., 72./255., 101./255., 112./255., 114./255.,  112./255.,  96./255.,   30./255.};

  Int_t Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255);
  gStyle->SetNumberContours(100);

  return;
}
