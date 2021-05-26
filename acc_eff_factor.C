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

void acc_eff_factor(){

//get reconstructed file
 TFile* file_reco = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResultsGammaGammaMedium.root");
 TList* lCMUL7DiMuonHist_reco = (TList*) file_reco->Get("ReconstructedHistos_CAny");

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_reco->FindObject("fHistoDimuonReconstructed");
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(2.));
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(1.9),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));

//get generated file
 TFile* file_gen = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResultsGammaGammaMedium.root");
 TList* lCMUL7DiMuonHist_gen = (TList*) file_gen->Get("GeneratedHistos_CAny");

 THnSparseT<TArrayD>* fHistoJPsiGenerated = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_gen->FindObject("fHistoDimuonGenerated");
 fHistoJPsiGenerated->GetAxis(1)->SetRange(fHistoJPsiGenerated->GetAxis(1)->FindBin(0.),fHistoJPsiGenerated->GetAxis(1)->FindBin(2.));
 fHistoJPsiGenerated->GetAxis(0)->SetRange(fHistoJPsiGenerated->GetAxis(0)->FindBin(1.9),fHistoJPsiGenerated->GetAxis(0)->FindBin(5.));

 TH1D* hGen  = (TH1D*) fHistoJPsiGenerated->Projection(1);
 hGen->SetTitle("Generated");
 hGen->SetLineColor(kBlack);
 hGen->SetMarkerStyle(20);
 hGen->SetMarkerSize(0.4);

 TH1D* hReco = (TH1D*) fHistoDiMuonOS->Projection(1);
 hReco->SetTitle("Reconstructed");
 hReco->SetLineColor(kRed+1);
 hReco->SetMarkerStyle(4);
 hReco->SetMarkerSize(0.4);

 TH1D* hComp = (TH1D*) hGen->Clone("comp");
 hComp->SetTitle("Comparaison gen v rec - ggm");
 //hComp->Clear();

 //hReco->Scale(10);
 TH1D* hAccEff = (TH1D*) hReco->Clone("hAccEff");
 hAccEff->SetTitle("Acceptance efficacity - ggm; p_{T} GeV/c ; Acc x #epsilon");
 hAccEff->Divide(hReco,hGen);

 TCanvas* c1 = new TCanvas("c1", "acceff");
 c1->cd();
 hAccEff->SetStats(0);
 hAccEff->Draw("e");



 TCanvas* c2 = new TCanvas("c2", "generated v reconstructed");
 c2->cd();
 hComp->Draw("");
 hGen->Draw("esame");
 hReco->Draw("esame");

 c2->BuildLegend(0.4,0.65,0.7,0.85);

}