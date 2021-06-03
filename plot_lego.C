#include <iostream>
#include <string>
#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLatex.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include <cmath>
#include <vector>
#include <stdio.h>
#include "TMath.h"
#include "TObjArray.h"
#include "TRatioPlot.h"

TH1* gHisto_coh = nullptr; // Define a global object
TH1* gHisto_incoh = nullptr; // Define a global object

Double_t hist_coh(Double_t* x, Double_t* par) {
  // The fit function has only 1 free parameter, which is the normalisation.
  Double_t xx = x[0];
  // We use the global histo in the function, which has to be correctly initialised in the main.
  int ibin = gHisto_coh->GetXaxis()->FindBin(xx);
  return par[0] * gHisto_coh->GetBinContent(ibin);
}

Double_t hist_incoh(Double_t* x, Double_t* par) {
  // The fit function has only 1 free parameter, which is the normalisation.
  Double_t xx = x[0];
  // We use the global histo in the function, which has to be correctly initialised in the main.
  int ibin = gHisto_incoh->GetXaxis()->FindBin(xx);
  return par[0] * gHisto_incoh->GetBinContent(ibin);
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

Double_t GammaContribution(Double_t *x, Double_t *par)
{
    return LogNormal(x,par) + TailExpo(x, &par[3]); 
}


Double_t Linear(Double_t *x, Double_t *par)
{
    return par[0] * x[0];
}

Double_t FitFunctionPt(Double_t *x, Double_t *par)
{
    return GammaContribution(x, par) + hist_coh(x, &par[5]) 
    + hist_incoh(x, &par[6]) + Linear(x, &par[7]);
}

Double_t CrystalBallExtended(Double_t *y, Double_t *par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=nL for the left tail
	//par[5]=alphaR  Alpha for the right tail
	//par[6]=nR for the right tail


	Double_t t = (y[0]-par[1])/par[2];
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

Double_t poly6(Double_t *y, Double_t *par)
{
	return par[0] + par[1] *y[0]+ par[2] * TMath::Power(y[0],2) + par[3] * TMath::Power(y[0],3)
	+ par[4] * TMath::Power(y[0],4)+ par[5] * TMath::Power(y[0],5)+ par[6] * TMath::Power(y[0],6);
}

Double_t FitFunctionMass(Double_t *y, Double_t *par)
{
	return  poly6(y, par) + CrystalBallExtended(y,&par[7]) + CrystalBallExtended(y,&par[14]) + CrystalBallExtended(y,&par[21]);// HADRO JPSI
}

Double_t FitFunctionTotal(Double_t *x,  Double_t *par)
{
    return  FitFunctionPt(x,par) + FitFunctionMass(&x[1],&par[8]);
}



void plot_lego(){

 TFile* file_coh = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResults_cohJpsi_weighted.root");
 TList* lCMUL7DiMuonHist_coh = (TList*) file_coh->Get("ReconstructedHistos_CAny");

 THnSparseT<TArrayD>* fHistoDiMuonOS_coh = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_coh->FindObject("fHistoDiMuonOS");
 fHistoDiMuonOS_coh->GetAxis(1)->SetRange(fHistoDiMuonOS_coh->GetAxis(1)->FindBin(0.),fHistoDiMuonOS_coh->GetAxis(1)->FindBin(1.));
 fHistoDiMuonOS_coh->GetAxis(0)->SetRange(fHistoDiMuonOS_coh->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS_coh->GetAxis(0)->FindBin(5.));

 gHisto_coh = static_cast<TH1*>(fHistoDiMuonOS_coh->Projection(1));
 gHisto_coh->SetDirectory(0); // This is needed because otherwise gHisto belongs to TFile and it is deleted when TFIle is closed
 delete file_coh;


 TFile* file_incoh = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResults_incohJpsi_weighted.root");
 TList* lCMUL7DiMuonHist_incoh = (TList*) file_incoh->Get("ReconstructedHistos_CAny");

 THnSparseT<TArrayD>* fHistoDiMuonOS_incoh = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_incoh->FindObject("fHistoDiMuonOS");
 fHistoDiMuonOS_incoh->GetAxis(1)->SetRange(fHistoDiMuonOS_incoh->GetAxis(1)->FindBin(0.),fHistoDiMuonOS_incoh->GetAxis(1)->FindBin(1.));
 fHistoDiMuonOS_incoh->GetAxis(0)->SetRange(fHistoDiMuonOS_incoh->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS_incoh->GetAxis(0)->FindBin(5.));

 gHisto_incoh = static_cast<TH1*>(fHistoDiMuonOS_incoh->Projection(1));
 gHisto_incoh->SetDirectory(0); // This is needed because otherwise gHisto belongs to TFile and it is deleted when TFIle is closed
 delete file_incoh;

//fit function created with 36 parameters
 TF2* Fitfct = new TF2("Fitfct",FitFunctionTotal,0.,1.,1.9,5.,36);
 Double_t par[36];
 for(int i=0;i<36;i++){Fitfct->SetParameters(i,1.);}


//CB Extended parameters for the 3 different jpsi contributions (coh, incoh, hadro)
 Double_t AlphaL_coh = 0.8096;
 Double_t nL_coh = 135.747;
 Double_t AlphaR_coh= 2.3248;
 Double_t nR_coh = 7.3228;

 Double_t AlphaL_incoh = 0.8597;
 Double_t nL_incoh = 12.76 ; 
 Double_t AlphaR_incoh= 2.3387;
 Double_t nR_incoh = 4.2698;

 Double_t AlphaL_hadro = 0.890;
 Double_t nL_hadro = 8.735;
 Double_t AlphaR_hadro = 3.043;
 Double_t nR_hadro = 15.323;

//Parameters fixed for mass range [1.9-5] GeV/c2 and pt [0-1] GeV/c

 Fitfct->FixParameter(0,47.53);       //fgamma p0
 Fitfct->FixParameter(1,-2.167);      //fgamma p1
 Fitfct->FixParameter(2,0.7125);      //fgamma p2
 Fitfct->FixParameter(3,2.365e-6);    //fgamma p3
 Fitfct->FixParameter(4,0.1846);      //fgamma p4

 Fitfct->FixParameter(5,0.00059688);  //fcoh p0 norm template
 
 Fitfct->FixParameter(6,0.000319);    //fincoh p0 norm template

 Fitfct->FixParameter(7,569.3);       //fhadro p0

 //gg contribution
 Fitfct->FixParameter(8,-191.5); 
 Fitfct->FixParameter(9,262.3);
 Fitfct->FixParameter(10,-48.64);
 Fitfct->FixParameter(11,-12.86);
 Fitfct->FixParameter(12,1.643);
 Fitfct->FixParameter(13,0.7639);
 Fitfct->FixParameter(14,-0.1094);

 //coherent jpsi contribution
 Fitfct->FixParameter(15,174.595);
 Fitfct->FixParameter(16,3.108);
 Fitfct->FixParameter(17,0.05666);
 Fitfct->FixParameter(18,AlphaL_coh);
 Fitfct->FixParameter(19,nL_coh);
 Fitfct->FixParameter(20,AlphaR_coh);
 Fitfct->FixParameter(21,nR_coh);

 //incoherent jpsi contribution
 Fitfct->FixParameter(22,85.08);
 Fitfct->FixParameter(23,3.108);
 Fitfct->FixParameter(24,-0.05647);
 Fitfct->FixParameter(25,AlphaL_incoh);
 Fitfct->FixParameter(26,nL_incoh);
 Fitfct->FixParameter(27,AlphaR_incoh);
 Fitfct->FixParameter(28,nR_incoh);

 //hadronic jpsi contribution
 Fitfct->FixParameter(29,114.7);
 Fitfct->FixParameter(30,3.108);
 Fitfct->FixParameter(31,-0.05647);
 Fitfct->FixParameter(32,AlphaL_hadro);
 Fitfct->FixParameter(33,nL_hadro);
 Fitfct->FixParameter(34,AlphaR_hadro);
 Fitfct->FixParameter(35,nR_hadro);


 // file data
 TFile* fin = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResults_DataMerged.root");
 
 TList* lCMUL7DiMuonHist = (TList*) fin->Get("DiMuonHistos_CMUL7");

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist->FindObject("fHistoDiMuonOS");
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 fHistoDiMuonOS->GetAxis(3)->SetRange(fHistoDiMuonOS->GetAxis(3)->FindBin(70.), fHistoDiMuonOS->GetAxis(3)->FindBin(90.));

// make the 3D fit
 TH2D* h1 = fHistoDiMuonOS->Projection(0,1);
 h1->SetTitle("p_{T} M_{#mu#mu}; p_{T} GeV/c; M_{#mu#mu} GeV/c^{2}; entries");

 TH1D* hData = (TH1D*) fHistoDiMuonOS->Projection(1); //histo 1D for pt distri
 hData->SetLineColor(kBlue+1);
 hData->SetMarkerStyle(20);
 hData->SetMarkerSize(0.4);

 TH1D* hData_mass = (TH1D*) fHistoDiMuonOS->Projection(0); //histo 1D for mass distri
 hData_mass->SetLineColor(kBlue+1);
 hData_mass->SetMarkerStyle(20);
 hData_mass->SetMarkerSize(0.4);

 TCanvas* c1 = new TCanvas("c1","plot lego");
 c1->cd();
 gStyle->SetOptFit(1111);
 h1->Fit("Fitfct","REM","");
 h1->Draw("lego2");

//draw the pt and mass distribution
 TCanvas* c2 = new TCanvas("c2","proj pt");
 c2->Divide(1,2);
 c2->cd(1);

 TH1D* hProjPt = h1->ProjectionX("Pt_Distri",0,-1,"");
 hProjPt->SetTitle("p_{T} Distribution");
 hProjPt->Draw("e");

 c2->cd(2);

 TH1D* hProjMass = h1->ProjectionY("Mass_Distri",0,-1,"");
 hProjMass->SetTitle("M_{#mu#mu} Distribution");
 hProjMass->Draw("e");


// comparison between data and pt distribution 3D fit
 TH2D* hRandFit_Pt = new TH2D("hRandFit_Pt", "histo filled with fitfct; p_{T} GeV/c ;entries", 21, 0., 1., 156, 1.9, 5.);
 hRandFit_Pt->FillRandom("Fitfct",8777000);
 hRandFit_Pt->Scale(0.001);

 TH1D* h2 = hRandFit_Pt->ProjectionX();
 h2->SetLineColor(kRed);
 h2->SetMarkerStyle(4);
 h2->SetMarkerSize(0.4);

 hData->SetTitle("Data");

 TCanvas* c3 = new TCanvas("c3", "test rand fit");
 c3->cd();
 h2->Draw("e");
 h2->GetYaxis()->SetRangeUser(150.,1000);
 hData->Draw("esame");
 c3->BuildLegend(0.30,0.60,0.65,0.75);
 h2->SetTitle("Comparaison data V 2D fit function");


//comparison between data and mass distribution 3D fit
 TH2D* hRandFit_Mass = new TH2D("hRandFit_Mass", "histo filled with fitfct; p_{T} GeV/c ;M_{#mu#mu} GeV/c^{2}", 21, 0., 1., 156, 1.9, 5.);
 hRandFit_Mass->FillRandom("Fitfct",8777000);

 //scale to have a good comparison between data and fit
 hRandFit_Mass->Scale(0.001);

 TH1D* h3 = hRandFit_Mass->ProjectionY();
 h3->SetLineColor(kRed);
 h3->SetMarkerStyle(4);
 h3->SetMarkerSize(0.4);

 TCanvas* c4 = new TCanvas("c4", "test rand fit mass");
 c4->cd();
 h3->Draw("e");
 h3->GetYaxis()->SetRangeUser(0.001,500.);
 hData_mass->SetTitle("Data mass distribution");
 hData_mass->Draw("esame");
 
 c4->BuildLegend(0.30,0.60,0.65,0.75);
 h3->SetTitle("Comparaison data V 2D fit function");
 
 TCanvas* c5 = new TCanvas("c5", "test 2 windows",1300.,500.);
 c5->Divide(2,1);
 c5->cd(1);
 //h3->SetTitle("histo filled with fit function (rescaled)");
 h3->Draw("e");

 c5->cd(2);
 hData_mass->Draw("e");
}