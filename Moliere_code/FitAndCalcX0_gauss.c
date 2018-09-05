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
#include "ParList.h"

double GetTheta(TString filename, double data_percentage){
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
	TString FormFilename(int run){	
		TString fileStart = "/afs/desy.de/user/z/zakharos/kappa075_2kink/run0000";
		TString fileEnd = "-GBLKinkEstimator_kappa075_2kink.root";
		TString out;
		if(run < 10){
			out = fileStart + Form("0%d",run) + fileEnd;
		}else{
			out = fileStart + Form("%d",run) + fileEnd;
		}	
		return out;
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
void DrawMultiEGraph(double *X0, const int n_graph){
		//Init and Fill array of graphs for different energy
	TGraph *gr[n_graph];
	int gr_num = 0; 
	int graph_it[n_graph];
	for (int i = 0; i < n_graph; i++){ 
		graph_it[i] = 0;
		gr[i] = new TGraph;	
	}
	for (int i = 0; i < 30; i++){
		double x_real = al[i].width;
		double x_calc = X0[i];
		gr_num = al[i].E-1;	
		gr[gr_num]->SetPoint(graph_it[gr_num],x_real,x_calc);
		graph_it[gr_num]++;
	}
	
	TCanvas *c2 = new TCanvas("c2","corr",800,600);
	c2->cd();
	TMultiGraph *mg = new TMultiGraph();
	
	TLegend* legend = new TLegend(0.1,0.6,0.3,0.9);
	legend->SetHeader("X corellation"); // option "C" allows to center the header
	legend->Draw();
	for (int i = 0; i < n_graph; i++){
		gr[i]->SetMarkerStyle(22);
		gr[i]->SetMarkerSize(1.5);
		gr[i]->SetMarkerColor(kGreen+i);
		gr[i]->Fit("pol1","NP");
		TString entery_name =Form("E= %d GeV",i+1);
		legend->AddEntry(gr[i],entery_name,"p");
		mg->Add(gr[i]);
	}
	mg->Draw("AP");
	legend->Draw();
	c2->SetLogx();
	c2->SetLogy();
	c2->Update();
} 

void FitAndCalcX0_gauss(void){ 
	dataSet* al = initSampleAl();
	TString filename;
	TString filename_bkg;
	filename = FormFilename(48);
	filename_bkg = FormFilename(3);
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	//cout<<filename.Data()<<endl;
	double theta = 0;
	double theta_bkg = 0;
	double theta_eff = 0;
	double X0[30];
	double absX0 = 88.9; //mm
	const int n_graph = 5;
	const int n_dots = 6;
	for (int i = 0; i < 30; i++){ 	//Fit experimental data and get array of X0 values;
			par[0] = al[i].E;			
			filename = FormFilename(al[i].run);
			filename_bkg = FormFilename(al[i].run_bkg);
	
			theta = GetTheta(filename.Data(),0.95);
			theta_bkg = GetTheta(filename_bkg.Data(),0.95);
			theta_eff = TMath::Sqrt(theta * theta - theta_bkg * theta_bkg);
			X0[i]= absX0 * CalculateX0(theta_eff);
			cout<<"x_real = "<<al[i].width<<"\tx_calc = "<<X0[i]<<endl;	
	
	}
	c1->Close();
	DrawMultiEGraph(X0, n_graph);
	

	double mean[n_dots];
	double sigma[n_dots];
   	for (int i = 0; i < n_dots; i++){
		mean[i] = 0;
		sigma[i] = 0;
	} 
	for (int i = 0; i < 30; i++){ 
		switch(al[i].width) {
			case 25 : mean[0] += X0[i];
				     break;     
			case 50 : mean[1] += X0[i];
				     break;
			case 100 : mean[2] += X0[i];
				     break;
			case 200 : mean[3] += X0[i];
				     break;
			case 1000 : mean[4] += X0[i];
				     break;
			case 10000 : mean[5] += X0[i];
				     break;
		}
	}
	for (int i = 0; i < n_dots; i++){
		mean[i] = mean[i]/((double)n_graph);
	} 

	for (int i = 0; i < 30; i++){ 
		switch(al[i].width) {
			case 25 : sigma[0] += (mean[0] - X0[i])*(mean[0] - X0[i]);
				     break;     
			case 50 : sigma[1] += (mean[1] - X0[i])*(mean[1] - X0[i]);
				     break;
			case 100 : sigma[2] += (mean[2] - X0[i])*(mean[2] - X0[i]);
				     break;
			case 200 : sigma[3] += (mean[3] - X0[i])*(mean[3] - X0[i]);
				     break;
			case 1000 : sigma[4] += (mean[4] - X0[i])*(mean[4] - X0[i]);
				     break;
			case 10000 : sigma[5] += (mean[5] - X0[i])*(mean[5] - X0[i]);
				     break;
		}
	}
	for (int i = 0; i < n_dots; i++){
		sigma[i] = TMath::Sqrt(sigma[i]/((double)n_graph-1.));
	} 
	double width[6] = { 0.025, 0.05, 0.1, 0.2, 1.,10.};
	double zeros[6] = { 0,0,0,0,0,0 };
	TCanvas *c3 = new TCanvas("c3","corr. graph",1600,1200);
	c3->cd();
	c3->SetLogx();
	c3->SetLogy();
	TGraphErrors* gr_err = new TGraphErrors(n_dots, width,mean,zeros,sigma);
	gr_err->SetTitle("");
	gr_err->SetMarkerColor(4);
    gr_err->SetMarkerStyle(21);
	gr_err->SetMarkerSize(1);
	gr_err->SetFillColor(kGreen-9);
	gr_err->GetXaxis()->SetTitle("Actual Width [mm]");
	gr_err->GetYaxis()->SetTitle("Measured Width [mm]");
	gr_err->GetXaxis()->SetTitleOffset(1.2);
	gr_err->Draw("A3");
	gr_err->Draw("PE");
	gr_err->Fit("pol1");
	
	gStyle->SetOptFit(111);
	TPaveStats* st =(TPaveStats*)(gr_err->GetListOfFunctions()->FindObject("stats"));
	st->SetX1NDC(.1);
	st->SetX2NDC(.4);
	st->SetY1NDC(.7);
	st->SetY2NDC(.9);
	c3->SetGridx();
	c3->SetGridy();
	/*TLine *line = new TLine(0,0,10,10);
	line->SetLineWidth(2);
	line->SetLineColor(kViolet);
	line->Draw("SAME");*/
	gPad->Modified();
	
}
	
