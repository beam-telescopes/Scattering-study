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
#include "X0Lib.h"

void ImageX0(void){
    double xl = -1; // frame set in mm
	double xr = 0;
	double yb = -1;
	double yt = 0;
	double step = 0.005;
	
	int nbinx = (int)((xr-xl)/step);
	int nbiny = (int)((yt-yb)/step);
	cout <<"Frame size: "<< nbinx <<"x"<< nbiny <<"\t"<<"Pixel size: "<<step<<" mm"<<endl;
	
	double xr_real = xl+step*(double)nbinx;
	double yt_real = yb+step*(double)nbiny;
	//cout << xr_real<<"\t" << yt_real <<endl;
	TH2D* pix_map = new TH2D("pix_map","pix_map", nbinx, xl, xr_real, nbiny, yb, yt_real);
	
	int startrun = 529;
	int endrun = 729;
	//int run_arr[9] = {194,195,196,198,199,200,201,202,203}; //frame 2 
	//int run_arr[11] = {35,36,37,38,39,40,41,42,43,44,45};

	TH1F **xkink_h = new TH1F*[nbinx*nbiny]; // kink angle histograms unique for all files;
	TH1F **ykink_h = new TH1F*[nbinx*nbiny];
	
	for(int i = 0; i < nbinx*nbiny; i++){
		xkink_h[i] = new TH1F(Form("xkink x=%d,y=%d",i%nbinx,i/nbinx), "xkink_h",1000, -0.05, 0.05); //LIMITS!!!!
		ykink_h[i] = new TH1F(Form("ykink x=%d,y=%d",i%nbinx,i/nbinx),"ykink_h",1000, -0.05, 0.05);
	}

	TString path = "/afs/desy.de/user/z/zakharos/dust/GEM_2_run/output/histograms/";	
	double n_events = 0;
 	for(int k = startrun; k <= endrun; k++){ // loop for working with multiple files
		TString filename  = FormFilename(path,k);
		TFile* file = TFile::Open(filename,"READ");
		if(file != NULL){
			TTree * tree = (TTree*)file->Get("KinkEstimator/KinkAngles"); 
			if(tree != NULL){
	
				std::vector<double>* x = 0;
				std::vector<double>* y = 0;
				std::vector<double>* kinkx = 0;
				std::vector<double>* kinky = 0;
	
				tree->SetBranchAddress("x",&x);
				tree->SetBranchAddress("y",&y);
				tree->SetBranchAddress("kinkx",&kinkx);
				tree->SetBranchAddress("kinky",&kinky);

				double xr_cur = xl;
				double yt_cur = yb;
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

						while(y->at(j) >= yt_cur && yt_cur < (yt_real-step)){
							biny++;
							yt_cur = yb + step * biny;	
							//cout<<"biny =  "<< biny <<endl;		
						}
						if (xr_cur >= xl && xr_cur < xr_real && yt_cur >= yb && yt_cur < yt_real){
							//cout<<"Adding event "<<count<<"\t"<<"Px= "<<xr_cur<<"\t"<<"Py= "<<yt_cur<<"\t"<<"X= "<<x->at(j)<<"\t"<<"Y="<<y->at(j)<<endl;
							//cout<<"binx = " << binx <<"\t"<<"biny = " << biny <<endl;				 
							xkink_h[ biny * nbinx + (nbinx-binx-1) ]->Fill(kinkx->at(j));
						 	ykink_h[ biny * nbinx + (nbinx-binx-1) ]->Fill(kinky->at(j));
							n_events++;			
						}
						binx = 0;
						biny = 0;
						xr_cur = xl;
						yt_cur = yb;
					}
					if ((int)(count) % (int)(100000) == 1) printf("run: %d\t Data processing: %2.0lf %% \n",k,count/tsize*100);		
					count++;	
					
				}
				delete tree;
				file->Close();
			}
		}
	}
 		
	double xtheta = 0;
	double ytheta = 0;
	double theta = 0;
	double relX = 0;
	double avHitX = 0;
	double avHitY = 0;
	double xq[2]={0.15,0.85};
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
		//xkink_h[i]->Delete();
		//ykink_h[i]->Delete();
	}

	TFile *out = new TFile("output.root","recreate");
	gStyle->SetOptStat(0);
	SetColorPalette();
	TLatex text;
	text.SetTextSize(0.04);
	text.SetTextAlign(13);
	text.DrawLatex(1,1.07,"X/X_{0}");
	//pix_map->GetZaxis()->SetRangeUser(0.0, 0.1);
	pix_map->GetZaxis()->SetTitle("X/Xo,rel");
	pix_map->GetZaxis()->SetTitle("X/Xo,rel");
	pix_map->GetXaxis()->SetTitle("x, mm");
	pix_map->GetYaxis()->SetTitle("y, mm");
	pix_map->GetYaxis()->SetTitleOffset(1);
	pix_map->Draw("colz");
	pix_map->Write("X0map");
	cout<< "hits/px = " <<n_events/(nbinx*nbiny)<<endl;
	cout<<"Total event number: "<<n_events<<endl;
	out->Close();
}
