{
//Draw experimental kink angle histo
TFile* file = TFile::Open("/afs/desy.de/user/z/zakharos/kappa075_2kink/run000033-GBLKinkEstimator_kappa075_2kink.root");
Fitter06->cd();
GBL->cd();
gblsumkxandsumky->SetAxisRange(-4., 4.,"X");
gblsumkxandsumky->Draw();
c1->SetLogy();
// Plot Moliere function
TFile * fitHist = TFile::Open("/afs/desy.de/user/z/zakharos/Moliere_code/Moliere_hist.root");
moliere->SetLineColor(kRed);
moliere->Scale(0.55e6);
moliere->Draw("SAME");
moliere->Fit("gaus");
}

