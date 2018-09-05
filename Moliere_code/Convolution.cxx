{
	TCanvas *c1 = new TCanvas("c1","c1",1000,500);
	c1->Divide(2,1);
	c1->cd(1);
	TH1F * h1 = new TH1F("h1","h1",100,-3,3);
	TH1F * h_conv = new TH1F("conv","conv",100,-3,3);
	double val[100];
	double h1_val = 0;
	double h2_val = 0;
	
	double FWHM1 = 2;
	double sigma1 = FWHM1/2.355;

	double FWHM2 = 1;
	double sigma2 = FWHM2/2.355;
	
	for(int i = 0; i <100; i++){
		val[i] = -3 + 6*i/100.; 	
	}
	for(int i = 0; i <100; i++){
		h1_val = TMath::Gaus(val[i],0,sigma1,1);
		h1->SetBinContent(i+1,h1_val);
	}
	h1->Scale(1/h1->Integral("width"));
	h1->Draw();
	double conv_integral = 0;

	for(int i = 0; i <100; i++){
		for(int j = 0; j <100; j++){	
			conv_integral += h1->GetBinContent(j+1)*ROOT::Math::tdistribution_pdf(val[i],0.2);
			cout<<i<<"\t"<<conv_integral<<endl;
		}
		h_conv->SetBinContent(i+1,conv_integral);	
		conv_integral = 0;
	}
	c1->cd(2);
	h_conv->Scale(h_conv->Integral("width"));
	h_conv->SetLineColor(2);
	h_conv->Draw();
    h1->Scale(h_conv->GetMaximum()/(h1->GetMaximum()));
	h1->Draw("SAME");
}
