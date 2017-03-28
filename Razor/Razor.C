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

void Razor() {

  // Declare your minitrees
  TString FileName[2] = {"../../minitrees/nominal/Stop/TTTo2L2Nu_ext1.root",
			 "../../minitrees/nominal/Stop/T2tt_mStop-150to1200.root"};
  int const nvar = 6;
  TString Region[6] = {"R1", "R_2", "R3", "R4", "R5", "R6"};
  TString svar[nvar]= {"MR", "R2", "Rpt", "invGamma", "Mdr", "DeltaPhiRll"};
  TString xTitle[nvar] = { "MR [GeV]", "R2", "Rpt", "invGamma ", "Mdr [GeV]", "DeltaPhiRll [rad]"};  
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
 
     DISCR[dt+dt*R][0] = new TH1F("MR"             + HistoName, "", 80, 0., 1000.);
     DISCR[dt+dt*R][1] = new TH1F("R2"             + HistoName, "", 80, 0., 1.4);
     DISCR[dt+dt*R][2] = new TH1F("Rpt"            + HistoName, "", 80, 0., 1.1);
     DISCR[dt+dt*R][3] = new TH1F("invGamma"       + HistoName, "", 80, 0., 1.5);
     DISCR[dt+dt*R][4] = new TH1F("Mdr"            + HistoName, "", 80, 0., 450.);
     DISCR[dt+dt*R][5] = new TH1F("DeltaPhiRll"    + HistoName, "", 80, 0, 3.5);
    
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
    
            TLorentzVector MET; MET.SetPtEtaPhiM(metPfType1, 0., metPfType1Phi, 0.);

      vector<TLorentzVector> Jets;
      TLorentzVector Jet1; Jet1.SetPtEtaPhiM(tjet1pt, tjet1eta, tjet1phi, tjet1mass);
      TLorentzVector Jet2; Jet2.SetPtEtaPhiM(tjet2pt, tjet2eta, tjet2phi, tjet2mass);
      Jets.push_back(Jet1);
      Jets.push_back(Jet2);

      vector<TLorentzVector> Leps;
      TLorentzVector Lep1; Lep1.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.);
      TLorentzVector Lep2; Lep2.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.);
      Leps.push_back(Lep1);
      Leps.push_back(Lep2);

      vector<TLorentzVector> Hemispheres = getHemispheres(Jets, Leps);

      float MR = computeMR(Hemispheres[0], Hemispheres[1]);
      float R2 = computeR2(Hemispheres[0], Hemispheres[1], MET);

      //cout << dt << " " << MR << " " << R2 << endl;
      
      TVector3 vBETA_z, pT_CM, vBETA_T_CMtoR, vBETA_R;
      double SHATR,  dphi_LL_vBETA_T,  dphi_L1_L2; 
      double gamma_R,  dphi_vBETA_R_vBETA_T;
      double MDELTAR,  costhetaRp1;
      
      SuperRazor(Leps, MET, vBETA_z, pT_CM,
      vBETA_T_CMtoR, vBETA_R,
      SHATR, dphi_LL_vBETA_T, dphi_L1_L2,
      gamma_R, dphi_vBETA_R_vBETA_T,
      MDELTAR, costhetaRp1);

      float Rpt = pT_CM.Mag()/(pT_CM.Mag() + SHATR/4.);
      float invGamma = 1./gamma_R;
      float Mdr = SHATR/gamma_R;
      float DeltaPhiRll = dphi_LL_vBETA_T;

      
      // Fill histograms
      
      V1 = metPfType1/meff;  
      V2 = metPfType1/sqrt(meff);  
      V3 = metPfType1/ht;  
      V4 = metPfType1/sqrt(ht);  
           
      float var [nvar] = { MR, R2, Rpt, invGamma, Mdr, DeltaPhiRll};
     
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

  TCanvas *CC = new TCanvas("CC", "", 1200,1000);
  CC->Divide(2, 3);
  gStyle->SetOptStat("");
 
  for (int vv = 0;vv<nvar; vv++) { 

  DISCR[0][vv]  ->SetLineColor(1); DISCR[0][vv]->SetLineWidth(2);
  float Integral = DISCR[0][vv]->Integral(0, DISCR[0][vv]->GetNbinsX()+1);
  DISCR[0][vv]->Scale(1/Integral);
  
  TString hTitle[6] = {"Mstop <= 240 && |Mstop-Mlsp-Mtop| <= 10" , " Mstop > 240 && |Mstop-Mlsp-Mtop| <= 10", "Mstop <= 240 && Mstop-Mlsp-Mtop < -10", "Mstop > 240 && Mstop-Mlsp-Mtop < -10", "Mstop <= 350 && Mstop-Mlsp-Mtop > 10", "Mstop > 350 && Mstop-Mlsp-Mtop > 10" };
  TString lHeader[6] = { "R1", "R_2", "R3", "R4", "R5", "R6"};
 
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
    Integral = DISCR[R+1][vv]->Integral(0, DISCR[R+1][vv]->GetNbinsX()+1);
    DISCR[R+1][vv]->Scale(1/Integral);

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
  CC->Print("./RazorPlots/" + svar[vv] + ".png");
  CC->Print("./RazorPlots/" + svar[vv] + ".pdf");
 } 
}
