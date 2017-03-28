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

#include "../Razor/Razor.C"
#include "../Razor/SuperRazor.C"

void StopTrees() {

  // Declare your minitrees
  TString FileName[2] = {"../../minitrees/nominal/Stop/TTTo2L2Nu_ext1.root",
			 "../../minitrees/nominal/Stop/T2tt_mStop-150to1200.root"};
  int const nvar = 22;
  TString Region[6] = {"R1", "R_2", "R3", "R4", "R5", "R6"};
  TString svar[nvar]= {"MT2ll", "MT2bb","MT2lblb", "dyll", "ptbll", "dphimetptbll", "m2l", "dphimetjet", "mllbb", "meff", "dphillmet", "dphill", "dphilmet1", "dphijet1met", "ht", "htjets", "metPfType1", "metPfType1Phi","V1","V2", "V3", "V4"};
/*,"MR", "R2", "Rpt", "invGamma", "Mdr", "DeltaPhiRll"*/
  TString xTitle[nvar] = {"MT2ll [GeV]", "MT2bb [GeV]","MT2lblb [GeV]", "dyll", "ptbll [GeV]", "dphimetptbll [rad]", "m2l [GeV]", "dphimetjet [rad]", "mllbb [GeV]", "meff [GeV]", "dphillmet [rad]", "dphill [rad]", "dphilmet1 [rad]", "dphijet1met [rad]", "ht [GeV]", "htjets [GeV]", "metPfType1 [GeV]", "metPfType1Phi [GeV]","V1 []", "V2 []", "V3 []", "V4 []"};
//, "MR [GeV]", "R2", "Rpt", "invGamma ", "Mdr [GeV]", "DeltaPhiRll [rad]"};  
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
 
     DISCR[dt+dt*R][0]  = new TH1F("MT2ll"         + HistoName, "", 60, 0., 600.);
     DISCR[dt+dt*R][1]  = new TH1F("MT2bb"         + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][2]  = new TH1F("MT2lblb"       + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][3]  = new TH1F("dyll"          + HistoName, "", 80, 0., 5.);
     DISCR[dt+dt*R][4]  = new TH1F("ptbll"         + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][5]  = new TH1F("dphimetptbll"  + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][6]  = new TH1F("m2l"           + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][7]  = new TH1F("dphimetjet"    + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][8]  = new TH1F("mllbb"         + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][9]  = new TH1F("meff"          + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][10] = new TH1F("dphillmet"     + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][11] = new TH1F("dphill"        + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][12] = new TH1F("dphilmet1"     + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][13] = new TH1F("dphijet1met"   + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][14] = new TH1F("ht"            + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][15] = new TH1F("htjets"        + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][16] = new TH1F("metPfType1"    + HistoName, "", 80, 0., 800.);
     DISCR[dt+dt*R][17] = new TH1F("metPfType1Phi" + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][17] = new TH1F("metPfType1Phi" + HistoName, "", 80, 0., 4.);
     DISCR[dt+dt*R][18] = new TH1F("V1=Met/Meff"   + HistoName, "", 80, 0., 1.);
     DISCR[dt+dt*R][19] = new TH1F("V1=Met/sqrt(Meff)"+ HistoName, "", 80, 0., 1.);
     DISCR[dt+dt*R][20] = new TH1F("V1=Met/ht"      + HistoName, "", 80, 0., 1.);
     DISCR[dt+dt*R][21] = new TH1F("V4=Met/sqrt(ht)"+ HistoName, "", 80, 0., 1.);
    
    }
 }
 
 for (int dt = 0; dt<2; dt++) {
    
    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    
    for (Int_t i = 0; i<nentries; i++) {
      
      MiniTree->GetEntry(i);
      
      // Apply selection
      if (njet<2) continue;
      if (nbjet30csvv2m<1) continue;
    
      
      // Fill histograms
      
      V1 = metPfType1/meff;  
      V2 = metPfType1/sqrt(meff);  
      V3 = metPfType1/ht;  
      V4 = metPfType1/sqrt(ht);  
           
      float var [nvar] = { mt2ll, mt2bb, mt2lblb, dyll, ptbll, dphimetptbll, m2l, dphimetjet, mllbb, meff, dphillmet, dphill, dphilmet1, dphijet1met, ht, htjets, metPfType1, metPfType1Phi,V1, V2, V3, V4};
//, MR, R2, Rpt, invGamma, Mdr, DeltaPhiRll};
     
      bool R1, R_2, R3, R4, R5, R6;
      if (dt>0) {
        R1 = Mstop <= 240 && abs(Mstop - Mlsp - Mtop) <= 10; 
        R_2 = Mstop > 240 && abs(Mstop - Mlsp - Mtop) <= 10; 
        R3 = Mstop <= 240 && Mstop - Mlsp -Mtop < -10;
        R4 = Mstop > 240 && Mstop - Mlsp -Mtop < -10;
        R5 = Mstop <= 350 && Mstop - Mlsp -Mtop >  10; 
        R6 = Mstop > 350 && Mstop - Mlsp -Mtop >  10;
       }
       
      for (int vv = 0;vv<nvar; vv++) {

        if (dt==0){
        	DISCR[0][vv]->Fill(var[vv], eventW); // Top
        } else {
  

          if (R1) DISCR[1][vv]->Fill(var[vv], eventW);
          if (R_2) DISCR[2][vv]->Fill(var[vv], eventW);
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
  float Integral = DISCR[0][vv]->Integral(0, DISCR[0][vv]->GetNbinsX()+1);
  DISCR[0][vv]->Scale(1/Integral);
  
  TString hTitle[6] = {"Mstop <= 240 && |Mstop-Mlsp-Mtop| <= 10" , " Mstop > 240 && |Mstop-Mlsp-Mtop| <= 10", "Mstop <= 240 && Mstop-Mlsp-Mtop < -10", "Mstop > 240 && Mstop-Mlsp-Mtop < -10", "Mstop <= 350 && Mstop-Mlsp-Mtop > 10", "Mstop > 350 && Mstop-Mlsp-Mtop > 10" };
  TString lHeader[6] = { "R1", "R_2", "R3", "R4", "R5", "R6"};
 
  for (int R = 0; R <6 ; R++){

    TPad *PD1 = (TPad*)CC->GetPad(R+1); PD1->SetLogy(); PD1->SetGridx(); PD1->SetGridy();
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
    Integral = DISCR[R+1][vv]->Integral(0, DISCR[R+1][vv]->GetNbinsX()+1);
    DISCR[R+1][vv]->Scale(1/Integral);

    DISCR[0][vv]->SetTitle (hTitle[R]);
    DISCR[0][vv]->SetXTitle (xTitle[vv]);
    DISCR[0][vv]->DrawCopy("histo");
    DISCR[R+1][vv]->DrawCopy("histsame");
    leg->Draw();
  }
  CC->Print("./sminitrees/" + svar[vv] + ".png");
  CC->Print("./sminitrees/" + svar[vv] + ".pdf");
 } 

}
