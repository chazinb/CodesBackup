#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
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


void JetsStudy() {

  // Declare your minitrees
     TString FileName[2] = {"../../minitrees/nominal/Stop/TTTo2L2Nu_ext1.root",
                            "../../minitrees/nominal/Stop/T2tt_mStop-350to400.root"};
  // Declare and initialize your counters
    
     TH1F* S_hcounter = new TH1F ("Signal_hcounter", "", 9, 0, 8);
     TH1F* B_hcounter = new TH1F ("Background_hcounter", "", 9, 0, 8);
     
     TH1F* significance;
     
        

//     float S_counters [12]; = { S_bjet20, S_bjet25, S_bjet30, S_b20Jet20, S_b20Jet25, S_b20Jet30, S_b25Jet20, S_b25Jet25, S_b25Jet30, S_b30Jet20, S_b30Jet25, S_b30Jet30 };
//      
//     float B_counters [12] = { B_bjet20, B_bjet25, B_bjet30, B_b20Jet20, B_b20Jet25, B_b20Jet30, B_b25Jet20, B_b25Jet25, B_b25Jet30, B_b30Jet20, B_b30Jet25, B_b30Jet30 };
//
//     for (int i = 0; i < 12 ; i++ ) {
//
//       S_counters[i] = 0; 
//       B_counters[i] = 0; 
//     }

   //  float S_bjet20 = 0;
   //  float S_bjet25 = 0; 
   //  float S_bjet30 = 0; 

   //  float S_b20Jet20 = 0;                           
   //  float S_b20Jet25 = 0;                           
   //  float S_b20Jet30 = 0;
   //    
   //  float S_b25Jet20 = 0;                           
   //  float S_b25Jet25 = 0;                           
   //  float S_b25Jet30 = 0;

   //  float S_b30Jet20 = 0;
   //  float S_b30Jet25 = 0;
   //  float S_b30Jet30 = 0;

     //float S_counters [12] = { S_bjet20, S_bjet25, S_bjet30, S_b20Jet20, S_b20Jet25, S_b20Jet30, S_b25Jet20, S_b25Jet25, S_b25Jet30, S_b30Jet20, S_b30Jet25, S_b30Jet30 };

  // Open the minitrees for your signal

     TFile *MiniTreeFile = TFile::Open("../../minitrees/nominal/Stop/T2tt_mStop-350to400.root"); 
     TTree *MiniTree     = GetMiniTree(MiniTreeFile);
     Int_t nentries      = (Int_t) MiniTree -> GetEntries();
    
     for (Int_t i = 0; i<nentries; i++){

       MiniTree -> GetEntry(i);

       if (Mstop != 350 || Mlsp != 175) continue; 

       if (nbjet20csvv2m >= 1 ) S_hcounter-> Fill (0, eventW); 
       if (nbjet20csvv2m > 1 && jet2pt > 20) { S_b20Jet20 += eventW ; S_hcounter-> Fill (1, eventW);} 
       //if (nbjet20csvv2m > 1 && jet2pt > 25) S_b20Jet25 += 1 && S_hcounter-> Fill (3); 
       //if (nbjet20csvv2m > 1 && jet2pt > 30) S_b20Jet30 += 1 && S_hcounter-> Fill (4); 
       
       if (nbjet25csvv2m > 1)  S_hcounter-> Fill (2, eventW);     
       if (nbjet25csvv2m > 1 && jet2pt > 20) S_hcounter-> Fill (3, eventW);
       if (nbjet25csvv2m > 1 && jet2pt > 25) S_hcounter-> Fill (4, eventW);
       //if (nbjet25csvv2m > 1 && jet2pt > 30) S_b25Jet30 += 1 && S_hcounter-> Fill (8);

       if (nbjet30csvv2m > 1 ) S_hcounter-> Fill (5, eventW);     
       if (nbjet30csvv2m > 1 && jet2pt > 20) S_hcounter-> Fill (6, eventW);
       if (nbjet30csvv2m > 1 && jet2pt > 25) S_hcounter-> Fill (7, eventW);
       if (nbjet30csvv2m > 1 && jet2pt > 30) S_hcounter-> Fill (8, eventW);

     }


    // float B_bjet20 = 0;
    // float B_bjet25 = 0;
    // float B_bjet30 = 0;

    // float B_b20Jet20 = 0;
    // float B_b20Jet25 = 0;
    // float B_b20Jet30 = 0;

    // float B_b25Jet20 = 0;
    // float B_b25Jet25 = 0;
    // float B_b25Jet30 = 0;

    // float B_b30Jet20 = 0;
    // float B_b30Jet25 = 0;
    // float B_b30Jet30 = 0;

    // float B_counters [12] = { B_bjet20, B_bjet25, B_bjet30, B_b20Jet20, B_b20Jet25, B_b20Jet30, B_b25Jet20, B_b25Jet25, B_b25Jet30, B_b30Jet20, B_b30Jet25, B_b30Jet30 };

    // Open the minitrees for your background

     TFile *B_MiniTreeFile = TFile::Open("../../minitrees/nominal/Stop/TTTo2L2Nu_ext1.root");
     TTree *B_MiniTree     = GetMiniTree(B_MiniTreeFile);
     Int_t B_nentries      = (Int_t) B_MiniTree -> GetEntries();

     for (Int_t i = 0; i<nentries; i++){

       MiniTree -> GetEntry(i);

       if (nbjet20csvv2m > 1 ) B_hcounter-> Fill (0, eventW);;
       if (nbjet20csvv2m > 1 && jet2pt > 20) B_hcounter-> Fill (1, eventW);
       //if (nbjet20csvv2m > 1 && jet2pt > 25) B_hcounter-> Fill (3);
       //if (nbjet20csvv2m > 1 && jet2pt > 30) B_hcounter-> Fill (4);
       
       if (nbjet25csvv2m > 1)  B_hcounter-> Fill (2, eventW);
       if (nbjet25csvv2m > 1 && jet2pt > 20) B_hcounter-> Fill (3, eventW);
       if (nbjet25csvv2m > 1 && jet2pt > 25) B_hcounter-> Fill (4, eventW);
       //if (nbjet25csvv2m > 1 && jet2pt > 20) B_b25Jet30 += 1 && B_hcounter-> Fill (8);

       if (nbjet30csvv2m > 1 ) B_hcounter-> Fill (5, eventW);
       if (nbjet30csvv2m > 1 && jet2pt > 20) B_hcounter-> Fill (6, eventW);
       if (nbjet30csvv2m > 1 && jet2pt > 25) B_hcounter-> Fill (7, eventW);
       if (nbjet30csvv2m > 1 && jet2pt > 30) B_hcounter-> Fill (8, eventW);

      }

    // Draw and save  a counter histo

    TCanvas *CCounter = new TCanvas("CCounter", "", 1200, 1000);
    gStyle->SetOptStat("");
    //CCounter -> Draw();
    S_hcounter-> Draw ();
   // hcounter[]->SetBarWidth(0.2);
   // hcounter->SetBarOffset(0.1);
    B_hcounter-> Draw("same");  
    S_hcounter->SetFillColor(49);
    //B_hcounter->SetFillColor(34);
    CCounter->Print("kkkk.png");
     
    // Fill and Draw a significance graphs
    // S/sqrt(S+B)

     gr1 = new TGraph();
     gr2 = new TGraph();
     for (int i = 0; i < 12; i ++){
  
       if (B_counters[i] == 0) i == -1 && B_counters[i] == 1; 
 
       cout<< "B_counters" << i << B_counters[i] << endl;

       float x = i;
       float y1 = (S_counters[i] / sqrt(S_counters[i] + B_counters[i])); 
       float y2 = (S_counters[i] / B_counters[i]); 
       gr1->SetPoint(i, x, y1);
       gr2->SetPoint(i, x, y2);

    }
    
    TCanvas *mycanvas = new TCanvas();

    gr1 -> Draw("APE");
    //gr -> SetTitle ( variable + "\t channel \t " + schannel[i]);
    //gr -> GetXaxis() -> SetTitle(variable + "(GeV)");
    gr1 -> GetYaxis() -> SetTitle("S/sqrt(S+B)");
    gr1 -> SetMarkerStyle (kFullDotSmall);
    gr1 -> SetMarkerSize (12);
    gr1 -> SetMarkerColor(kBlue);
                            
    //leg = new TLegend(0.7,0.8,0.9,0.9);
    //leg->SetHeader(icut); 	     
    //leg->Draw();
    mycanvas->SaveAs("significance.png");

}                          
