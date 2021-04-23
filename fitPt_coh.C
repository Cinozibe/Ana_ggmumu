#include <iostream>
#include <string>
#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLatex.h"
#include "TF1.h"
#include "TStyle.h"
#include "TF1Convolution.h"
//#include "FitFunctions.C"
#include <cmath>
#include <vector>
#include <stdio.h>
#include "TMath.h"
#include "TObjArray.h"

Double_t PowerLaw(Double_t *x, Double_t *par)
{
    return x[0] * par[0] * TMath::Power( 1 / TMath::Power(1+x[0]/par[1],par[2]) ,par[3]);
}


Double_t LogNormal(Double_t *x, Double_t *par)
{
    //par[1]=mean
    //par[2]=width
    //par[0]=normalization

    Double_t numExp =  (TMath::Log(x[0]) - par[1]) * (TMath::Log(x[0]) - par[1])  ;
    Double_t denExp = 2 * par[2] * par[2];
    
    return par[0]/(x[0] * par[2] * TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(- numExp / denExp) ;
    
}

Double_t Gaussian(Double_t *x, Double_t *par)
{
    //par[0]=normalization
    //par[1]=mean
    //par[2]=width
    
    return x[0] * par[0] * TMath::Exp( -0.5 * (x[0]-par[1]) * (x[0]-par[1]) / (par[2] * par[2]) );
}

Double_t TailExpo(Double_t *x, Double_t *par)
{
    return x[0] * par[0] * TMath::Exp( - x[0] / par[1]);
}


Double_t sum(Double_t *x, Double_t *par)
{
    return PowerLaw(x,par) + Gaussian(x, &par[4]) + Gaussian(x, &par[7]);
}

/*
Double_t sum(Double_t *x, Double_t *par)
{
    return PowerLaw(x,par) + Gaussian(x, &par[4]) + TailExpo(x, &par[7]);
}
*/
void fitPt_coh()
{
 
 
 TFile* fin = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/AnalysisResults_cohJpsi_weighted.root");
 //Dimuon Histo
 TList* lCAnyDiMuonHist = (TList*) fin->Get("ReconstructedHistos_CAny");
 lCAnyDiMuonHist->ls();

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCAnyDiMuonHist->FindObject("fHistoDiMuonOS");
 
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(2.6));
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 
 Double_t par[10];
 
 TF1* PowerLawfct = new TF1("PowerLawfct", PowerLaw, 0.,1.,4);
 PowerLawfct->SetParameters(1.,1.,1.,1.);
 PowerLawfct->SetLineStyle(4);
 PowerLawfct->SetLineColor(kRed);
 
 TF1* LogNormalfct = new TF1("LogNormalfct", LogNormal, 0.,1.,3);
 LogNormalfct->SetParameters(1.,1.,1.);
 LogNormalfct->SetLineStyle(4);
 LogNormalfct->SetLineColor(kGreen);

 TF1* Gaussianfct = new TF1("Gaussianfct", Gaussian, 0.,1.,3);
 Gaussianfct->SetParameters(1.,1.,1.);
 Gaussianfct->SetLineStyle(4);
 Gaussianfct->SetLineColor(kBlack);

 TF1* Gaussianfct2 = new TF1("Gaussianfct2", Gaussian, 0.,1.,3);
 Gaussianfct2->SetParameters(1.,1.,1.);
 Gaussianfct2->SetLineStyle(4);
 Gaussianfct2->SetLineColor(kGreen);

 TF1* Tail = new TF1("Tail", TailExpo, 0., 1., 2);
 Tail->SetParameters(1.,1.);
 Tail->SetLineColor(kMagenta);
 Tail->SetLineStyle(4);


 TF1* GlobalFit = new TF1("GlobalFit", sum, 0.,1.,10);
 GlobalFit->SetParameters(1.,1.,1.,1. ,1.,1.,1., 1.,1.,1.);
 //GlobalFit->SetParLimits(0,1.,500.);
 GlobalFit->SetParLimits(2,0.,10.);
 //GlobalFit->SetParLimits(7,0.,1.e7);
 //GlobalFit->SetParLimits(5,0.,2.);
 //GlobalFit->SetParLimits(8,0.,2.);
 //GlobalFit->SetParLimits(9,0.,1000.);
 GlobalFit->SetLineColor(kBlue);


 TH1D* h1 = fHistoDiMuonOS->Projection(1);
 h1->SetTitle("Fit PowerLaw + Gaus + Gaus - Pt Distribution - mass range [1.9 - 2.6] - coh; p_{T} GeV/c; # DiMuons");
 //h1->Fit("GlobalFit","R","");

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(2.6),fHistoDiMuonOS->GetAxis(0)->FindBin(3.4));
 TH1D* h2 = fHistoDiMuonOS->Projection(1);
 h2->SetTitle("Fit PowerLaw + Gaus + Gaus - Pt Distribution - mass range [2.6 - 3.4] - coh; p_{T} GeV/c; # DiMuons");
 //GlobalFit->SetParameters(1.,1.,1.,1.,1.,1.,1.);
 //h2->Fit("GlobalFit","R","");

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(3.4),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 TH1D* h3 = fHistoDiMuonOS->Projection(1);
 h3->SetTitle("Fit PowerLaw + Gaus + Gaus - Pt Distribution - mass range [3.4 - 5] - coh; p_{T} GeV/c; # DiMuons");
 //GlobalFit->SetParameters(1.,1.,1.,1.,1.,1.,1.);
 //h3->Fit("GlobalFit","R","");

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 TH1D* h4 = fHistoDiMuonOS->Projection(1);
 h4->SetTitle("Fit PowerLaw + Gaus + Gaus - Pt Distribution - mass range [1.9 - 5] - coh; p_{T} GeV/c; # DiMuons");
 //h4->Fit("GlobalFit","R","");
 //GlobalFit->SetParameters(1.,1.,1.,1.,1.,1.,1.);
 

 TCanvas* c1 = new TCanvas("FitPowerLawTwoGaus_Pt Distribution_1.9_5._coh","mass fct pt coh");
 c1->cd();
 gPad->SetLogy();
 gStyle->SetOptFit(1111);
 h3->Draw("e");
 h3->Fit("GlobalFit","EM","",0.,1.);
 
 const char *formula = "f( p_{T} ) = p_{T} p0 #left( #frac{ 1 }{ #left(1 + #frac{ p_{T} }{ p1 } #right)^{p2} } #right)^{p3}+ p_{T} p4 exp#left(- #frac{#left(p_{T} - p5#right)^{2}}{2 p6^{2}}#right) + p7p_{T} exp#left(- #frac{#left(p_{T} - p8#right)^{2}}{2 p9^{2}}#right)";
 TLatex latex;
 latex.SetTextSize(0.025);
 latex.DrawLatex(.1,100.,formula);
 
 GlobalFit->GetParameters(par);
 PowerLawfct->SetParameters(par);
 Gaussianfct->SetParameters(&par[4]);
 Gaussianfct2->SetParameters(&par[7]);
 PowerLawfct->Draw("same");
 Gaussianfct->Draw("same");
 Gaussianfct2->Draw("same");


/*
 
 c1->cd(2);gPad->SetLogy();
 h2->Draw("e");
 c1->cd(3);gPad->SetLogy();
 h3->Draw("e");
 c1->cd(4);gPad->SetLogy();
 h4->Draw("e");

 TFile *fout = TFile::Open("FittedLogN_coh.root", "RECREATE");
 fout->cd(); 
 h1->Write("PtDistribution Mass 19 26 coh");
 h2->Write("PtDistribution Mass 26 34 coh");
 h3->Write("PtDistribution Mass 34 50 coh");
 h4->Write("PtDistribution Mass 19 50 coh");
 c1->Write("PtDistribution coh");
 
 delete fout;
 */
}