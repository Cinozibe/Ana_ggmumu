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


Double_t CrystalBallExtended(Double_t *x, Double_t *par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=nL for the left tail
	//par[5]=alphaR  Alpha for the right tail
	//par[6]=nR for the right tail


	Double_t t = (x[0]-par[1])/par[2];
	if (par[2] < 0) t = -t;
	
	Double_t absAlphaL = TMath::Abs(par[3]);
  	Double_t absAlphaR = TMath::Abs(par[5]);

  	if (t > -absAlphaL && t < absAlphaR) // gaussian core
  	{
   		return par[0]*(TMath::Exp(-0.5*t*t));
 	}
		
	
	if (t <= -absAlphaL) //left tail
  	{
    		Double_t A =  TMath::Power(par[4]/absAlphaL,par[4])*TMath::Exp(-0.5*absAlphaL*absAlphaL);
    		Double_t B = par[4]/absAlphaL - absAlphaL;
	
   		return par[0]*(A/TMath::Power(B - t, par[4]));
 	}
  
  	
  	
  	if (t >= absAlphaR) //right tail
 	{
 		Double_t C =  TMath::Power(par[6]/absAlphaR,par[6])*TMath::Exp(-0.5*absAlphaR*absAlphaR);
    		Double_t D = par[6]/absAlphaR - absAlphaR;
    		
 		return par[0]*(C/TMath::Power(D + t, par[6]));
  	}
  
  	return 0. ; 
} 

Double_t poly6(Double_t *x, Double_t *par)
{
	return par[0] + par[1] *x[0]+ par[2] * TMath::Power(x[0],2) + par[3] * TMath::Power(x[0],3)
	+ par[4] * TMath::Power(x[0],4)+ par[5] * TMath::Power(x[0],5)+ par[6] * TMath::Power(x[0],6);
}

Double_t sum(Double_t *x, Double_t *par)
{
	return  poly6(x, par) + CrystalBallExtended(x,&par[7]) + CrystalBallExtended(x,&par[14]);
}



void fitMass_data()
{

 // file data
 TFile* fin = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/AnalysisResults_DataMerged.root");
 
 TList* lCMUL7DiMuonHist = (TList*) fin->Get("DiMuonHistos_CMUL7");
 lCMUL7DiMuonHist->ls();

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist->FindObject("fHistoDiMuonOS");
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(.3));
 fHistoDiMuonOS->GetAxis(3)->SetRange(fHistoDiMuonOS->GetAxis(3)->FindBin(70.), fHistoDiMuonOS->GetAxis(3)->FindBin(90.));


 //contribution JPsi
 //CBExtended
 Double_t AlphaL_coh = 0.8096;
 Double_t nL_coh = 135.747;
 Double_t AlphaR_coh= 2.3248;
 Double_t nR_coh = 7.3228;

 Double_t AlphaL_incoh = 0.8597;
 Double_t nL_incoh = 12.76 ; 
 Double_t AlphaR_incoh= 2.3387;
 Double_t nR_incoh = 4.2698;

 TF1* Fitfct = new TF1("Fitfct", sum, 1.9,5.,21);
 Double_t par[21];
 for(int i=0; i<21;i++){Fitfct->SetParameters(i,1.);}
 Fitfct->SetLineColor(kRed);

 TF1* CBE_coh = new TF1("CBE_coh",CrystalBallExtended,1.9,5.,7);
 CBE_coh->SetParameters(1.,1.,1.,1.,1.,1.,1.);
 CBE_coh->SetLineStyle(4);
 CBE_coh->SetLineColor(kMagenta);
 /*
 Fitfct->FixParameter(1,7053);
 Fitfct->FixParameter(2,3379);
 Fitfct->FixParameter(3,-1672);
 Fitfct->FixParameter(4,-75.32);
 Fitfct->FixParameter(5,95.63);
 Fitfct->FixParameter(6,-9.487);
*/

 Fitfct->FixParameter(7,171.4);
 Fitfct->FixParameter(8,3.108);
 Fitfct->FixParameter(9,0.05666);
 Fitfct->FixParameter(10,AlphaL_coh);
 Fitfct->FixParameter(11,nL_coh);
 Fitfct->FixParameter(12,AlphaR_coh);
 Fitfct->FixParameter(13,nR_coh);


 TF1* CBE_incoh = new TF1("CBE_incoh",CrystalBallExtended,1.9,5.,7);
 CBE_incoh->SetParameters(10.,1.,1.,1.,1.,1.,1.);
 CBE_incoh->SetLineStyle(4);
 CBE_incoh->SetLineColor(kGreen);

 Fitfct->FixParameter(15,3.108);
 Fitfct->FixParameter(16,-0.05647);
 Fitfct->FixParameter(17,AlphaL_incoh);
 Fitfct->FixParameter(18,nL_incoh);
 Fitfct->FixParameter(19,AlphaR_incoh);
 Fitfct->FixParameter(20,nR_incoh);

 //contribution gamma gamma
 TF1* fgamma = new TF1("fgamma", poly6, 1.9,5.,7);
 fgamma->SetParameters(1.,1.,1.,1.,1.,1.,1.);
 fgamma->SetLineStyle(4);
 fgamma->SetLineColor(kBlue);

 TH1D* h1 = fHistoDiMuonOS->Projection(0);
 h1->SetTitle("Fit - Mass Distribution - pT [0. - 0.3] GeV/c - data; M_{#mu #mu} GeV/c^{2}; # DiMuons");
 
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.3),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 TH1D* h2 = fHistoDiMuonOS->Projection(0);
 h2->SetTitle("Fit - Mass Distribution - pT [0.3 - 1.] GeV/c - data; M_{#mu #mu} GeV/c^{2}; # DiMuons");

 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 TH1D* h3 = fHistoDiMuonOS->Projection(0);
 h3->SetTitle("Fit - Mass Distribution - pT [0. - 1.] GeV/c - data; M_{#mu #mu} GeV/c^{2}; # DiMuons");

 TCanvas* c1 = new TCanvas("MassDistribution data","mass fct pt data");
 c1->cd();gPad->SetLogy();
 gStyle->SetOptFit(1111);
 h1->Fit("Fitfct","REM","");
 h1->Draw("e");
 Fitfct->GetParameters(par);
 fgamma->SetParameters(par);
 CBE_coh->SetParameters(&par[7]);
 CBE_incoh->SetParameters(&par[14]);
 fgamma->Draw("same");
 CBE_coh->Draw("same");
 CBE_incoh->Draw("same");
 c1->BuildLegend(0.67,0.14,0.90,0.28);
 

 Double_t bin_width = h1->GetBinWidth(1);
 cout << "Bin width: " << bin_width << endl;

 Double_t N_jpsi_coh = CBE_coh->Integral(3.4,5.)/bin_width;
 cout << "N_Jpsi_coh= " << N_jpsi_coh << endl;

 Double_t N_jpsi_incoh = CBE_incoh->Integral(3.4,5.)/bin_width;
 cout << "N_Jpsi_incoh= " << N_jpsi_incoh << endl;
 
}