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
#include <cmath>
#include <vector>
#include <stdio.h>
#include "TMath.h"
#include "TObjArray.h"

Double_t PowerLaw(Double_t *x, Double_t *par)
{
    return x[0] * par[0] *TMath::Power(x[0],par[1]);
}

Double_t LogNormal(Double_t *x, Double_t *par)
{
    //par[1]=mean
    //par[2]=width
    //par[0]=normalization

    Double_t numExp =  ( (TMath::Log(x[0])) - par[1] ) * ( (TMath::Log(x[0])) - par[1] )  ;
    Double_t denExp = 2 * par[2] * par[2];
    
    return par[0]/(x[0] * par[2] * TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(- numExp / denExp) ;
    
}

Double_t TailExpo(Double_t *x, Double_t *par)
{
    return x[0] * par[0] * TMath::Exp( - x[0] / par[1]);
}

Double_t Gaussian(Double_t *x, Double_t *par)
{
    return x[0] * par[0] * TMath::Exp( -0.5 * (x[0]-par[1]) * (x[0]-par[1]) / (par[2] * par[2]) );
}
Double_t GammaContribution(Double_t *x, Double_t *par)
{
    return LogNormal(x,par) + TailExpo(x, &par[3]); 
}

Double_t CohJPsiContribution(Double_t *x, Double_t *par)
{
    return LogNormal(x,par) + TailExpo(x, &par[3]);
}
 
Double_t IncohJPsiContribution(Double_t *x, Double_t *par)
{
    return x[0] * par[0] * TMath::Exp( -0.5 * (x[0]-par[1]) * (x[0]-par[1]) / (par[2] * par[2]) );
}

Double_t Linear(Double_t *x, Double_t *par)
{
    return par[0] * x[0];
}

Double_t FitFunction(Double_t *x, Double_t *par)
{
    return GammaContribution(x, par) + CohJPsiContribution(x, &par[5]) 
    + IncohJPsiContribution(x, &par[10]) + Linear(x, &par[13]);
}


void fitPt_data()
{
 
 TF1* fgamma     = new TF1("fgamma", GammaContribution, 0., 3., 5);
 TF1* fcohjpsi   = new TF1("fcohjpsi", CohJPsiContribution, 0.,3.,5);
 TF1* fincohjpsi = new TF1("fincohjpsi", IncohJPsiContribution, 0., 3., 3);
 TF1* fhadro     = new TF1("fhadro", Linear, 0., 3., 1); //hadronic contribution

 //global fit function
 TF1* Fitfct     = new TF1("Fitfct",FitFunction,0.,3.,14);
 Double_t par[14];
 for(int i=0; i<14;i++){Fitfct->SetParameters(i,1.);}
 Fitfct->SetLineColor(kRed);
 Fitfct->SetLineWidth(4);

 //gamma gamma contribution
 fgamma->SetParameters(1.,1.,1.,1.,1.);
 fgamma->SetLineStyle(4);
 fgamma->SetLineWidth(4);
 fgamma->SetLineColor(kBlue);
 
 //jpsi coh contribution
 fcohjpsi->SetParameters(1.,1.,1., 1.,1.);
 fcohjpsi->SetLineStyle(4);
 fcohjpsi->SetLineWidth(4);
 fcohjpsi->SetLineColor(kMagenta);

 //jpsi incoh contribution
 fincohjpsi->SetParameters(1.,1.,1.);
 fincohjpsi->SetLineStyle(4);
 fincohjpsi->SetLineWidth(4);
 fincohjpsi->SetLineColor(kGreen+1);

 //hadro
 fhadro->SetParameters(0,1.);
 fhadro->SetLineStyle(4);
 fhadro->SetLineWidth(4);
 fhadro->SetLineColor(kBlack);

//Parameters for mass range [3.4-5.]
 
 Fitfct->SetParLimits(0,1.,100.);
 Fitfct->SetParLimits(1,-4.,-1.);  //fgamma p1
 Fitfct->SetParLimits(2,0.4,1.);   //fgamma p2
 Fitfct->SetParLimits(3,0.,300.);  //fgamma p3
 Fitfct->FixParameter(4,0.1922);   //fgamma p4

 Fitfct->FixParameter(5,0.0475);    //fcoh p0 norm
 Fitfct->FixParameter(6,-0.7683);    //fcoh p1
 Fitfct->FixParameter(7,0.468);      //fcoh p2
 Fitfct->FixParameter(8,0.21825);   //fcoh p3 norm
 Fitfct->FixParameter(9,0.3596);     //fcoh p4
 
 Fitfct->FixParameter(10,0.2524);   //fincoh p0 norm
 Fitfct->FixParameter(11,-0.1312);   //fincoh p1
 Fitfct->FixParameter(12,0.5448);    //fincoh p2

 Fitfct->SetParLimits(13,30.71,32.65); //fhadro p0

 // file data
 TFile* fin = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/AnalysisResults_DataMerged.root");
 
 TList* lCMUL7DiMuonHist = (TList*) fin->Get("DiMuonHistos_CMUL7");
 lCMUL7DiMuonHist->ls();

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist->FindObject("fHistoDiMuonOS");
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 fHistoDiMuonOS->GetAxis(3)->SetRange(fHistoDiMuonOS->GetAxis(3)->FindBin(70.), fHistoDiMuonOS->GetAxis(3)->FindBin(90.));


 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(2.6));
 TH1D* h1 = fHistoDiMuonOS->Projection(1);
 h1->SetTitle("Fit - Pt Distribution - mass range [1.9 - 2.6] - data; p_{T} GeV/c; # DiMuons");
 h1->SetLineColor(kBlack);
 h1->SetMarkerStyle(20);
 h1->SetMarkerSize(0.4);

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(2.6),fHistoDiMuonOS->GetAxis(0)->FindBin(3.4));
 TH1D* h2 = fHistoDiMuonOS->Projection(1);
 h2->SetTitle("Fit - Pt Distribution - mass range [2.6 - 3.4] - data; p_{T} GeV/c; # DiMuons");
 h2->SetLineColor(kBlack);
 h2->SetMarkerStyle(20);
 h2->SetMarkerSize(0.4);

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(3.4),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 TH1D* h3 = fHistoDiMuonOS->Projection(1);
 h3->SetTitle("Fit - Pt Distribution - mass range [3.4 - 5] - data; p_{T} GeV/c; # DiMuons");
 h3->SetLineColor(kBlack);
 h3->SetMarkerStyle(20);
 h3->SetMarkerSize(0.4);

 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 TH1D* h4 = fHistoDiMuonOS->Projection(1);
 h4->SetTitle("Fit - Pt Distribution - mass range [1.9 - 5] - data; p_{T} GeV/c; # DiMuons");
 h4->SetLineColor(kBlack);
 h4->SetMarkerStyle(20);
 h4->SetMarkerSize(0.4);



 TCanvas* c1 = new TCanvas("PtDistribution data","mass fct pt data");

 const char *formula = "f(p_{T}) = p0 p_{T}";
 TLatex latex;
 latex.SetTextSize(0.025); gStyle->SetOptFit(1111);
 latex.DrawLatex(.2,10.,formula);

 c1->cd(); gPad->SetLogy();
 gStyle->SetOptFit(1111);

 h3->Fit("Fitfct","EM","",0.,1.);
 h3->Draw("esame");
 
 Fitfct->GetParameters(par);
 fgamma->SetParameters(par);
 fcohjpsi->SetParameters(&par[5]);
 fincohjpsi->SetParameters(&par[10]);
 fhadro->SetParameters(&par[13]);
 
 fgamma->Draw("same");
 fcohjpsi->Draw("same");
 fincohjpsi->Draw("same");
 fhadro->Draw("same");
 c1->BuildLegend(0.23,0.68,0.6,0.88);
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(3.4),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));

 Double_t N_JPsiCoh = fcohjpsi->Integral(0,0.3)/0.05;
 cout << "N_JPsiCoh= " << N_JPsiCoh << endl;
 Double_t N_JPsiIncoh = fincohjpsi->Integral(0,0.3)/0.05;
 cout << "N_JPsiIncoh= " << N_JPsiIncoh << endl;
 Double_t N_GG = fgamma->Integral(0,0.3)/0.05;
 cout << "N_GG= " << N_GG << endl;
 Double_t N_hadro = fhadro->Integral(0,0.3)/0.05;
 cout << "N_hadro= " << N_hadro << endl;

}