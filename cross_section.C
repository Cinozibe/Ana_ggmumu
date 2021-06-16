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



void cross_section()
{

 const char *rapidity = "-4 < y < -2.5";
 const char *mass = "3.4 < m_{#mu#mu} < 5";
 TLatex latex;
 latex.SetTextSize(0.025);
 
//get reconstructed file
 TFile* file_reco = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResultsGammaGammaMedium.root");
 TList* lCMUL7DiMuonHist_reco = (TList*) file_reco->Get("ReconstructedHistos_CAny");

 THnSparseT<TArrayD>* fHistoDiMuonOS = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_reco->FindObject("fHistoDimuonReconstructed");
 fHistoDiMuonOS->GetAxis(1)->SetRange(fHistoDiMuonOS->GetAxis(1)->FindBin(0.),fHistoDiMuonOS->GetAxis(1)->FindBin(1.));
 fHistoDiMuonOS->GetAxis(0)->SetRange(fHistoDiMuonOS->GetAxis(0)->FindBin(3.4),fHistoDiMuonOS->GetAxis(0)->FindBin(5.));

//get generated file
 TFile* file_gen = TFile::Open("/Users/nicolasbize/Desktop/Stage_ALICE/analysisfiles/AnalysisResultsGammaGammaMedium.root");
 TList* lCMUL7DiMuonHist_gen = (TList*) file_gen->Get("GeneratedHistos_CAny");

 THnSparseT<TArrayD>* fHistoJPsiGenerated = (THnSparseT<TArrayD>*) lCMUL7DiMuonHist_gen->FindObject("fHistoDimuonGenerated");
 fHistoJPsiGenerated->GetAxis(1)->SetRange(fHistoJPsiGenerated->GetAxis(1)->FindBin(0.),fHistoJPsiGenerated->GetAxis(1)->FindBin(1.));
 fHistoJPsiGenerated->GetAxis(0)->SetRange(fHistoJPsiGenerated->GetAxis(0)->FindBin(3.4),fHistoJPsiGenerated->GetAxis(0)->FindBin(5.));

 TH1D* hGen  = (TH1D*) fHistoJPsiGenerated->Projection(1);
 hGen->SetTitle("Generated");
 hGen->Rebin(2); //2 bins are merged into one
 hGen->SetLineColor(kRed+1);
 hGen->SetMarkerStyle(20);
 hGen->SetMarkerSize(0.4);

 TH1D* hReco = (TH1D*) fHistoDiMuonOS->Projection(1);
 hReco->SetTitle("Reconstructed");
 hReco->Rebin(2);
 hReco->SetLineColor(kRed+1);
 hReco->SetMarkerStyle(20);
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
 hAccEff->GetXaxis()->SetRangeUser(0.,.3);
 hAccEff->SetStats(0);
 hAccEff->Draw("e");
 hAccEff->GetYaxis()->SetRangeUser(0.001,.65);
 latex.DrawLatex(.25,.4,mass);

 TH1D* hCS = (TH1D*) hAccEff->Clone();
 //TH1D* hCS = new TH1D("hCS","Cross Section", 3, 0.,0.3);

 TH1D* hN_gg = (TH1D*) hAccEff->Clone("hN");
 hN_gg->SetTitle("N_{#gamma#gamma->#mu#mu}");
 hN_gg->SetStats(1);
 hN_gg->Reset();

/*
//[1.9-2.6] GeV/c^2
 for (int i_0 = 0; i_0<147 ; i_0 ++){hN_gg->Fill(0);}
 for (int i_1 = 0; i_1<103  ; i_1 ++){hN_gg->Fill(0.1);}
 for (int i_2 = 0; i_2<52  ; i_2 ++){hN_gg->Fill(0.2);}

 //[2.6-3.4] GeV/c^2
 for (int i_0 = 0; i_0<92 ; i_0 ++){hN_gg->Fill(0);}
 for (int i_1 = 0; i_1<158  ; i_1 ++){hN_gg->Fill(0.1);}
 for (int i_2 = 0; i_2<85  ; i_2 ++){hN_gg->Fill(0.2);}
*/
 //[3.4-5] GeV/c^2
 for (int i_0 = 0; i_0<179 ; i_0 ++){hN_gg->Fill(0);}
 for (int i_1 = 0; i_1<92  ; i_1 ++){hN_gg->Fill(0.1);}
 for (int i_2 = 0; i_2<18  ; i_2 ++){hN_gg->Fill(0.2);}
 
 //for (int i=0;i<10;i++){hUnit->AddBinContent(i);}

 TCanvas* c2 = new TCanvas("c2", "test N");
 c2->cd();
 hN_gg->GetXaxis()->SetRangeUser(0.,0.3);
 hN_gg->GetYaxis()->SetRangeUser(0.,200.);
 hN_gg->Draw("e");
 latex.DrawLatex(.25,130.,mass);

 hCS->Divide(hN_gg,hAccEff);
 hCS->SetTitle("Cross section ggtomu - [3.4-5.] GeV/c^{2};p_{T} GeV/c; #frac{d#sigma_{#gamma#gamma to #mu#mu}}{dy} #mub");
 hCS->SetLineColor(kRed+1);
 hCS->SetMarkerStyle(20);
 hCS->SetMarkerSize(0.4);

 Double_t BR = 1.;
 Double_t L_int = 758.7;
 Double_t Delta_y = 1.5;
 Double_t a = 1./ BR / L_int / Delta_y;

 TCanvas* c3 = new TCanvas("c3", "cross section");
 c3->cd();
 gPad->SetLogy();
 hCS->SetStats(0);
 hCS->GetXaxis()->SetRangeUser(0.,0.3);
 hCS->Scale(a);
 hCS->Draw("e");

 latex.DrawLatex(.25,1.,rapidity);
 latex.DrawLatex(.25,0.5,mass);

}