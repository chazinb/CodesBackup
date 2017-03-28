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
#include "TAxis.h"
#include "StopTrees.h"



#include <fstream>
#include <iostream>

//----------------------------------------------
// Global Variables
// ---------------------------------------------

     TString FileAddress = "../AnalysisCMS/minitreesDY/nominal/Stop/";

     const int nsample = 4;
     TString SampleName [nsample] = { "01_Data.root", "07_ZJets.root", "02_WZTo3LNu.root", "03_VZ.root"}; 
     TString  _sample   [nsample] = {"Data", "ZJets", "WZ", "VZ"};
     // Data  = runG
     // ZJets = DYJetsToLL_M-10to50.root DYJetsToLL_M-50_*.root
     // WZ    = WZTo3LNu.root
     // VZ    = ZZTo2L2Q*.root WZTo2L2Q__part*.root

     // ######################################################################################
     enum { ee, mm, em, nchannel}; const TString schannel[nchannel] = {"ee", "mm","em"};
 
     TH2D* h_metPfType1_m2l [nchannel];




                  
void FillTH2() {


  // Define Histograms //
  // --------------------------------------------------------------------------------

 
     for (int ichannel = 0; ichannel < nchannel; ichannel++) {

       h_metPfType1_m2l[ichannel] = new TH2D("h_metPfType1_m2l_" + schannel[ichannel], "", 400, 0, 400, 120, 20, 140); 
     }


  // Loop //
  // ---------------------------------------------------------------------------------

  for (int s = 0; s < 1 ; s++){
//  for (int s = 0; s < nsample ; s++){

 
     TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]); 
     TTree *MiniTree      = GetMiniTree(MiniTreeFile);
     Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

    for (int ichannel = 0; ichannel < nchannel; ichannel++) {
       h_metPfType1_m2l[ichannel]-> Reset(); 
    }
    for (Int_t en = 0; en<nentries; en++){
       MiniTree -> GetEntry(en);

       // Block the overflow //
       // ---------------------------------------------------------------------------
         if (metPfType1 >= 400) metPfType1 = 400 - 0.01;
         if (m2l >= 140) m2l = 140 - 0.01; 
       // --------------------------------------------------------------------------

        int i=-1;

        if (channel == 3) { i = ee; }
        else if (channel == 4) { i = mm; }
        else if (channel == 5) { i = em; }
            
	if (i==-1) { cout << "Warning! # l < 2 " <<  channel <<  endl; continue;}	

        
        if (leadingPtCSVv2M <  20.) h_metPfType1_m2l[i]->Fill(metPfType1, m2l, eventW);
     }
     for (int ch = 0; ch < nchannel; ch ++){
      TFile* OutFileName = new TFile("TH2metPfType1_m2l/h_metPfType1_m2l_" + _sample[s] + ".root" , "update");
      h_metPfType1_m2l[ch] -> Write();
      OutFileName -> Close();
     }
   }

 }         
         

  


