#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "StopTrees.h"

#include <fstream>
#include <iostream>


void JetPt30_Study() {

  // Declare your minitrees
  TString FileName[2] = {"../../minitreesPT/nominal/Stop/TTTo2L2Nu_ext1.root",
			 "../../minitreesPT/nominal/Stop/T2tt_mStop-150to1200.root"};
  int const nvar = 18;
  TString svar[nvar]= {"jet1phi", "jet1pt", "jet1eta", "jet2pt", "jet2phi", "jet2eta", "lep1eta", "lep1phi", "lep1pt", "lep2eta", "lep2phi", "lep2pt", "njet", "Met_ht", "Met_sqrt_ht", "Met_meff", "Met_sqrt_meff", "mt2ll"};
  TString xTitle[nvar] = {"jet1phi [Gev]", "jet1pt [Gev]", "jet1eta [Gev]", "jet2pt [Gev]", "jet2phi [Gev]", "jet2eta [Gev]", "lep1eta [Gev]", "lep1phi [Gev]", "lep1pt [Gev]", "lep2eta [Gev]", "lep2phi [Gev]", "lep2pt [Gev]", "njet", "Met/ht", "Met/sqrt(ht)", "Met/meff", "Met/sqrt(meff)","MT2ll [GeV]"};  

  TString Region[6] = {"R1", "R2", "R3", "R4", "R5", "R6"};

  // Book histograms
  TH1F *DISCR[7][nvar];
   // R1 = MStop < 240 && abs(MStop - Mlsp -Mtop) < 10 (GeV) == GAP 
   // R2 = MStop > 240 && abs(MStop - Mlsp -Mtop) < 10 (GeV)
   // R3 = MStop < 240 && MStop - Mlsp -Mtop < -10 (GeV)
   // R4 = MStop > 240 && MStop - Mlsp -Mtop < -10 (GeV)
   // R5 = MStop < 350 && MStop - Mlsp -Mtop >  10 (GeV)
   // R03 = MStop > 350 && MStop - Mlsp -Mtop >  10 (GeV) 
  
 
  for (int R= 0; R < 6; R++){
    for (int dt = 0; dt<2; dt++) {
     if (dt == 0 && R>0) continue;
     TString HistoName = (dt==0) ? "_Top" : ("_Stop_"+ Region[R]);
     cout<< HistoName <<endl; 
 
     DISCR[dt+dt*R][0]  = new TH1F("jet1phi"          + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][1]  = new TH1F("jet1pt"         + HistoName, "", 60, 0., 300.);
     DISCR[dt+dt*R][2]  = new TH1F("jet1eta"         + HistoName, "", 60, -3., 3.);
     DISCR[dt+dt*R][3]  = new TH1F("jet2pt"          + HistoName, "", 60, 0., 300.);
     DISCR[dt+dt*R][4]  = new TH1F("jet2phi"         + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][5]  = new TH1F("jet2eta"         + HistoName, "", 60, -3., 3.);
     DISCR[dt+dt*R][6]  = new TH1F("lep1eta"         + HistoName, "", 60, -3., 3.);
     DISCR[dt+dt*R][7]  = new TH1F("lep1phi"         + HistoName, "", 80, -4., 4.);
     DISCR[dt+dt*R][8]  = new TH1F("lep1pt"          + HistoName, "", 60, 0., 300.);
     DISCR[dt+dt*R][9]  = new TH1F("lep2eta"         + HistoName, "", 80, -4., 4.);
     DISCR[dt+dt*R][10]  = new TH1F("lep2phi"        + HistoName, "", 40, 0., 4.);
     DISCR[dt+dt*R][11]  = new TH1F("lep2pt"         + HistoName, "", 60, 0., 300.);
     DISCR[dt+dt*R][12]  = new TH1F("njet"           + HistoName, "", 60, 0., 10.);
     DISCR[dt+dt*R][13]  = new TH1F("met/ht"         + HistoName, "", 60, 0., 0.7);
     DISCR[dt+dt*R][14]  = new TH1F("met/sqrt(ht)"   + HistoName, "", 60, 0., 15);
     DISCR[dt+dt*R][15]  = new TH1F("met/meff"       + HistoName, "", 60, 0., 0.7);
     DISCR[dt+dt*R][16]  = new TH1F("met/sqrt(meff)" + HistoName, "", 60, 0., 15);
     DISCR[dt+dt*R][17]  = new TH1F("MT2ll"         + HistoName, "", 100, 20., 120.);
    // DISCR[dt+dt*R][18]  = new TH1F("MT2bb"         + HistoName, "", 80, 0., 800.);
    // DISCR[dt+dt*R][19]  = new TH1F("MT2lblb"       + HistoName, "", 80, 0., 800.);
    

    }
 }
 
 for (int dt = 0; dt<2; dt++) {
    
    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    
    for (Int_t i = 0; i<nentries; i++) {
      
      MiniTree->GetEntry(i);
      
      // Apply selection
      if (nbjet20csvv2m<1) continue;
      //if (nbjet20csvv2m<1) continue;
   

      // Fill histograms
      
      V1 = metPfType1/ht;  // ht = MET.Et() + Lepton1.v.Pt() + Lepton2.v.Pt() + for (int ijet=0; ijet<_njet; ijet++) {_ht += AnalysisJets[ijet].v.Pt();} 
      V2 = metPfType1/sqrt(ht);  
      V3 = metPfType1/meff; // meff = MET.Pt() + Lepton1.v.Pt() + Lepton2.v.Pt() + AnalysisJets[0].v.Pt() + AnalysisJets[1].v.Pt(); 
      V4 = metPfType1/sqrt(meff);  
           
      float var [nvar] = { jet1phi, jet1pt, jet1eta, jet2pt, jet2phi, jet2eta, lep1eta, lep1phi, lep1pt, lep2eta, lep2phi, lep2pt, njet, V1, V2, V3, V4, mt2ll};
    
      if (njet==1 && jet1pt < 30.) continue; 
  
      if (njet > 1 && jet1pt < 30. && jet2pt < 30.) continue; 
 
      bool R1, R2, R3, R4, R5, R6;
      if (dt>0) {
        R1 = Mstop <= 240 && abs(Mstop - Mlsp - Mtop) <= 10; 
        R2 = Mstop > 240 && abs(Mstop - Mlsp - Mtop) <= 10; 
        R3 = Mstop <= 240 && Mstop - Mlsp -Mtop < -10;
        R4 = Mstop > 240 && Mstop - Mlsp -Mtop < -10;
        R5 = Mstop <= 350 && Mstop - Mlsp -Mtop >  10; 
        R6 = Mstop > 350 && Mstop - Mlsp -Mtop >  10;
       }
       
      for (int vv = 0;vv<nvar; vv++) {

        if (njet<2 && svar[vv].Contains("jet2")) continue; 
        if (njet<2 && svar[vv].Contains("meff")) continue;
  
        
        if (dt==0){
        	DISCR[0][vv]->Fill(var[vv], eventW); // Top
        } else {
  

          if (R1) DISCR[1][vv]->Fill(var[vv], eventW);
          if (R2) DISCR[2][vv]->Fill(var[vv], eventW);
          if (R3) DISCR[3][vv]->Fill(var[vv], eventW);
          if (R4) DISCR[4][vv]->Fill(var[vv], eventW);
          if (R5) DISCR[5][vv]->Fill(var[vv], eventW);
          if (R6) DISCR[6][vv]->Fill(var[vv], eventW);

       }
      }
    }
  }

  TCanvas *CC = new TCanvas("CC", "", 1200, 1000);
  CC->Divide(2, 3);
  gStyle->SetOptStat("");

   
  for (int vv = 0;vv<nvar; vv++) { 
 
  DISCR[0][vv]  ->SetLineColor(1); DISCR[0][vv]->SetLineWidth(2);
  float integral = DISCR[0][vv]->Integral(0, DISCR[0][vv]->GetNbinsX()+1);
  DISCR[0][vv]->Scale(1/integral);
  
  TString hTitle[6] = {"Mstop <= 240 && |Mstop-Mlsp-Mtop| <= 10" , " Mstop > 240 && |Mstop-Mlsp-Mtop| <= 10", "Mstop <= 240 && Mstop-Mlsp-Mtop < -10", "Mstop > 240 && Mstop-Mlsp-Mtop < -10", "Mstop <= 350 && Mstop-Mlsp-Mtop > 10", "Mstop > 350 && Mstop-Mlsp-Mtop > 10" };
  TString lHeader[6] = { "R1", "R2", "R3", "R4", "R5", "R6"};
 
  for (int R = 0; R <6 ; R++){

    TPad *PD1 = (TPad*)CC->GetPad(R+1); PD1->SetGridx(); PD1->SetGridy();
    PD1->cd();

    TLegend  *leg = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
    leg->SetTextColor(1); leg->SetTextSize(0.045);
    leg->SetTextFont(62); 
    leg->SetHeader(svar[vv] + "  " + lHeader[R]); //leg_fit->SetMargin(0.2); 
    leg->AddEntry((TObject*)0, " ", "");
    leg->AddEntry(DISCR[0][vv], " TTbar", "f");
    leg->AddEntry(DISCR[R+1][vv], " T2tt", "f");
    

    DISCR[R+1][vv]->SetLineColor(2); DISCR[R+1][vv]->SetLineWidth(2);
    integral = DISCR[R+1][vv]->Integral(0, DISCR[R+1][vv]->GetNbinsX()+1);
    DISCR[R+1][vv]->Scale(1/integral);

    //DISCR[0][vv]->SetTitle (hTitle[R]);
    //DISCR[0][vv]->SetXTitle (xTitle[vv]);
    int fp1=0, fp2=R+1;
    if (DISCR[0][vv]->GetMaximum() < DISCR[R+1][vv]->GetMaximum()) { fp1 =R+1; fp2 =0; }
    DISCR[fp1][vv]->SetTitle (hTitle[R]);
    DISCR[fp1][vv]->SetXTitle (xTitle[vv]);
    DISCR[fp1][vv]->DrawCopy("histo");
    DISCR[fp2][vv]->DrawCopy("histsame");
    leg->Draw();
  }

  CC->Print("./Pt_Study/pt30/" + svar[vv] + ".png");
  CC->Print("./Pt_Study/pt30/" + svar[vv] + ".pdf");
 } 
}
