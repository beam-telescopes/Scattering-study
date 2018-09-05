#include <cmath>
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

void GaussianFit(void){ 
	TFile* file = TFile::Open("/afs/desy.de/user/z/zakharos/kappa075_2kink/run000033-GBLKinkEstimator_kappa075_2kink.root");
	TDirectory * dir;
	dir = file->GetDirectory("Fitter06/GBL");	
	dir->cd();
	TH1F *hist = (TH1F*)gDirectory->Get("gblsumkxandsumky");
	hist->SetAxisRange(-4., 4.,"X");
	hist->Draw();		//Draw experimental kink angle histo
	TCanvas *c1 = new TCanvas("c1","c1",800,600);	
	c1->SetLogy();
	hist->Fit("gaus","","",-1,1);
	TF1 *gausf = hist->GetFunction("gaus");
	double sigma = gausf->GetParameter(2);
	double theta = sigma*1.1775;		//HWHM
	cout<<"<Theta> = "<<theta<<endl;
}


