#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include <cstdlib>
#include "TChain.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "ParListAl.h"


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

void ImgProcess_s(void){
    double xl = -9;
	double xr = 9;
	double yl = -4.5;
	double yr = 4.5;
	double step = 0.1;
	int nbinx = (int)((xr-xl)/step);
	int nbiny = (int)((yr-yl)/step);
	cout <<"Frame size: "<< nbinx <<"x"<< nbiny <<"\t"<<"Pixel size: "<<step<<" mm"<<endl;
	
	double xr_real = xl+step*(double)nbinx;
	double yr_real = yl+step*(double)nbiny;
	//cout << xr_real<<"\t" << yr_real <<endl;
	TH2D* pix_map = new TH2D("pix_map","pix_map", nbinx, xl, xr_real, nbiny, yl, yr_real);
	TString filename = "/afs/desy.de/user/z/zakharos/TelescopeExp/Si_det_analysis/output/histograms/run000035-gblkinkestimator.root";
	//TString filename = "/nfs/dust/atlas/user/arling/forStepan/20180418_analysis_atlasX0_Feb2017_TB21_alu3mm_1p0GeV.root";
	TFile* file = TFile::Open(filename,"READ");
	TTree * tree = (TTree*)file->Get("KinkEstimator/KinkAngles"); // for Al add to the GBL before KinkEstimator 
	
	std::vector<double>* x = 0;
	std::vector<double>* y = 0;
	std::vector<double>* kinkx = 0;
	std::vector<double>* kinky = 0;
	
	tree->SetBranchAddress("x",&x);
	tree->SetBranchAddress("y",&y);
	tree->SetBranchAddress("kinkx",&kinkx);
	tree->SetBranchAddress("kinky",&kinky);
	
	TH1F **xkink_h = new TH1F*[nbinx*nbiny];
	TH1F **ykink_h = new TH1F*[nbinx*nbiny];

	for(int i = 0; i < nbinx*nbiny; i++){
			xkink_h[i] = new TH1F(Form("xkink x=%d,y=%d",i%nbinx,i/nbinx), "xkink_h",1000, -0.05, 0.05); //LIMITS!!!!
			ykink_h[i] = new TH1F(Form("ykink x=%d,y=%d",i%nbinx,i/nbinx),"ykink_h",1000, -0.05, 0.05);

	}
	double xr_cur = xl;
	double yr_cur = yl;
	int binx = 0;
	int biny = 0;
	int tsize = tree->GetEntries();
	double count = 0;
	//cout<<"Tree size is "<<tsize<<endl;
	for(int i = 0; i < tsize; ++i){
		tree->GetEvent(i);
		//cout<<"New entry "<<i<<endl;
		for ( int j = 0; j < x->size(); ++j){
			//cout<<x->at(j)<<"\t"<<y->at(j)<<"\t"<<endl;
			while( xr_cur<= x->at(j) && xr_cur < xr_real-step){ 
				binx++;
				xr_cur = xl + step * binx;	
				//cout<<"binx =  "<< binx <<endl;		
			}

			while(y->at(j) >= yr_cur && yr_cur < yr_real-step){
				biny++;
				yr_cur = yl + step * biny;	
				//cout<<"biny =  "<< biny <<endl;		
			}
			if (xr_cur >= xl && xr_cur <= xr_real && yr_cur >= yl && yr_cur <= yr_real){
				//cout<<"Adding event "<<count<<"\t"<<"Px= "<<xr_cur<<"\t"<<"Py= "<<yr_cur<<"\t"<<"X= "<<x->at(j)<<"\t"<<"Y="<<y->at(j)<<endl;				 
				xkink_h[ biny * nbinx + binx ]->Fill(kinkx->at(j));
			 	ykink_h[ biny * nbinx + binx ]->Fill(kinky->at(j));		
			}
			binx = 0;
			biny = 0;
			xr_cur = xl;
	        yr_cur = yl;
		}
		//printf("Data processing: %2.1lf %\n",count/tsize*100);		
		count++;		
	} 	
	double xtheta = 0;
	double ytheta = 0;
	double theta = 0;
	double relX = 0;
	double avHitX = 0;
	double avHitY = 0;
	double xq[2]={0.1,0.9};
	double yq[2] ={0,0};
	

	for(int i = 0; i < nbinx*nbiny; i++){
		xkink_h[i]->GetQuantiles(2,yq,xq);
		xkink_h[i]->GetXaxis()->SetRangeUser(yq[0],yq[1]);
		
		ykink_h[i]->GetQuantiles(2,yq,xq);
		ykink_h[i]->GetXaxis()->SetRangeUser(yq[0],yq[1]);

		xtheta = xkink_h[i]->GetRMS();
		ytheta = ykink_h[i]->GetRMS();
		avHitX += xkink_h[i]->GetEntries();
		avHitY += ykink_h[i]->GetEntries();
		theta = TMath::Sqrt(1/2.*(xtheta*xtheta + ytheta*ytheta));
		relX = CalculateX0(theta*1000);	//theta in milirads
		pix_map->SetBinContent(i%nbinx+1,i/nbinx+1,relX);
	}
	TFile *out = new TFile("output.root","recreate");
	gStyle->SetPalette(kBlueGreenYellow);	
	pix_map->Draw("colz");
	pix_map->Write("X0map");
	out->Close();
	cout<<"Average hits/px number = "<<avHitX/(nbinx*nbiny)<<endl;
	/*int xent = 0;
	int yent = 0;
	for(int i = 0; i < nbinx*nbiny; i++){
		if(i%10 = 0){
			x_ent = xkink_h[i]->GetEntries();	
			y_ent = ykink_h[i]->GetEntries();
		}
		//cout<<"Nent(hx"<<i<<")="<<x_ent<<"\n"<<"Nent(hy"<<i<<")="<<y_ent<<"\n"<<endl;	
	}*/

}
