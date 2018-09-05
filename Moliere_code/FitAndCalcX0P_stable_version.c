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
#include <cstdlib>
#include "TChain.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"


double par[7] = {1, 1, 0.000512, 2.7, 13, 26.98, 0.025};

double GetThetaOld(TString filename, double data_percentage){
	TFile* file = TFile::Open(filename);
	TDirectory * dir;
	dir = file->GetDirectory("Fitter06/GBL");	
	dir->cd();
	TH1F *hist = (TH1F*)gDirectory->Get("gblsumkxandsumky");
	hist->SetAxisRange(-4., 4.,"X");
	hist->Draw();		//Draw experimental kink angle histo
	double integral = hist->Integral();
	double frac_integral = 0;
	int center_bin = 0; 
	center_bin = hist->FindBin(0); 
	
	int i = 0; 	
	for(i = 1; frac_integral <= data_percentage * integral; i++){
		frac_integral = hist->Integral(center_bin-i,center_bin+i,"");
	}

	double l_lim = hist->GetBinCenter(center_bin-i);
	double r_lim = hist->GetBinCenter(center_bin+i-1);
	//cout<<"% = "<<frac_integral/integral<<endl;	
	
	hist->Fit("gaus","","",l_lim,r_lim);
	TF1 *gausf = hist->GetFunction("gaus");
	double sigma = gausf->GetParameter(2);
	return sigma;
}
double GetTheta(TString filename, double data_percentage, double xl, double xr, double yb, double yt){
	TFile* file = TFile::Open(filename);
	if(file->IsOpen()){
		TTree * tree = (TTree*)file->Get("Fitter06/KinkAngles"); // for Al add to the GBL before KinkEstimator 
	
		std::vector<double>* x = 0;
		std::vector<double>* y = 0;
		std::vector<double>* kinkx = 0;
		std::vector<double>* kinky = 0;

		tree->SetBranchAddress("x",&x);	
		tree->SetBranchAddress("y",&y);
		tree->SetBranchAddress("kinkx",&kinkx);
		tree->SetBranchAddress("kinky",&kinky);
	
		TH1F * xkink_h = new TH1F("xkink_h", "xkink_h",1000, -0.05, 0.05);
		TH1F * ykink_h = new TH1F("ykink_h", "ykink_h",1000, -0.05, 0.05);
	
		int tsize = tree->GetEntries();
		for(int i = 0; i <tsize; i++){
			tree->GetEvent(i);
			for (int j = 0; j < kinkx->size(); j++)
			{
				if (x->at(j) >= xl && x->at(j) < xr && y->at(j) >= yb && y->at(j) < yt){			
					xkink_h->Fill(kinkx->at(j)-0.000434);
					ykink_h->Fill(kinky->at(j));
				}
			}
		}
		TCanvas * kink_plotX = new TCanvas("kink_agle_x","kink_angle_x",800,600);
		TCanvas * kink_plotY = new TCanvas("kink_agle_y","kink_angle_y",800,600);
		//kink_plotX->SetLogy();
		
		xkink_h->GetXaxis()->SetRangeUser(-0.005,0.005);
		xkink_h->GetXaxis()->SetTitle("Kink angle, rad");
		xkink_h->GetYaxis()->SetTitle("Number of events");
		ykink_h->GetXaxis()->SetRangeUser(-0.02,0.02);
		ykink_h->GetXaxis()->SetTitle("Kink angle, rad");
		ykink_h->GetYaxis()->SetTitle("Number of events");
		//gStyle->SetOptStat(0);
		
		kink_plotX->cd();
		xkink_h->Draw();
		kink_plotY->cd();
		xkink_h->Draw();

		
		double xq[2]={0.3,0.7};
		double yq[2] ={0.,0.};
		double meanx = 0;
		double meany = 0;

		xkink_h->GetQuantiles(2,yq,xq);
		
		TLine * llx = new TLine(yq[0],0,yq[0],xkink_h->GetMaximum());		
		TLine * lrx = new TLine(yq[1],0,yq[1],xkink_h->GetMaximum());
		llx->SetLineColor(kRed);
		lrx->SetLineColor(kRed);
		llx->SetLineWidth(2);
		lrx->SetLineWidth(2);
		kink_plotX->cd();		
		llx->Draw();
		lrx->Draw();
		kink_plotX->Update();
		
		kink_plotY->cd();
		ykink_h->GetQuantiles(2,yq,xq);
		ykink_h->Fit("gaus","","",yq[0],yq[1]);
		cout<<"q1 = "<<yq[0]<<"\t"<<"q2 = "<<yq[1]<<endl;
		
		TF1 *gausf_x = xkink_h->GetFunction("gaus");
		TF1 *gausf_y = ykink_h->GetFunction("gaus");
		double sigmax = gausf_x->GetParameter(2);
		double sigmay = gausf_y->GetParameter(2);

		return TMath::Sqrt((sigmax*sigmax + sigmay*sigmay)/2.);
	}else{ return -1;}
}

double CalculateX0(double theta){
	double X0min = 0.0001;
	double X0max = 1;  
	double tolerance = 1e-8;
	double X0 = (X0min+X0max)/2.0;	
	double error = X0max;

	double Ze;		// atomic number of target material
	Ze=par[1];

	double E = par[0];		// Energy in GeV
	double mass = par[2];	// mass of the particle in GeV
	double p = TMath::Sqrt(E*E-mass*mass);
	double beta = p/E;
	cout<<"beta = "<<beta<<endl;
	cout<<"p = "<<p<<endl;
	for (int i = 1; TMath::Abs(error) > tolerance; i++){
		X0 = (X0min+X0max)/2.0;
		error = theta/1000 - (13.6*0.001/(beta*p)*Ze*TMath::Sqrt(X0)*(1+0.038*TMath::Log(X0))); //Highland formula
		//error = theta/1000 - 13.6*0.001/p*TMath::Power(X0,0.555); //Approximation 20% error
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

void FitAndCalcX0(void){ 
	
	double theta = 0;
	double theta_bkg = 0;
	double theta_eff = 0;
	double X0;
	double absX0 = 88.9; //mm
	const int n_graph = 5;
	const int n_dots = 6;
	par[0] = 2.4;	

	double xl = -9; // frame set in mm
	double xr = 9;
	double yb = -4;
	double yt = 4;		
	TString filename = "/nfs/dust/atlas/user/michaela/ITk/X0-analysis/analysis/output/histograms/run000110-GBLKinkEstimator_kappa100_testtwoPointScat_decorr.root";
	theta = GetTheta(filename.Data(),0.95,xl,xr,yb,yt);
	theta_bkg = 0;//GetTheta(filename_bkg.Data(),0.95);
	theta_eff = TMath::Sqrt(theta * theta - theta_bkg * theta_bkg);
	X0 = CalculateX0(theta_eff);
	cout<<"X0_calc = "<< X0 <<endl;	
	
}
	
