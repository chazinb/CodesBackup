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

     const int nbck  = 3;
     TString SampleName [1] = {"T2tt_mStop-350to400.root"};// Sample mass point [1][2] = {(350,175)}
     TString BckName    [3] = {"TTTo2L2Nu_ext1.root", "DYJetsToLL.root", "WWTo2L2Nu.root"};
     TString LegName    [3] = {"tt", "DY", "WW"};
 
     const int ncut  = 4;
     float varCut [ncut] = { nbjet20csvv2m , metPfType1, V4, dphimetjet};
     float valCut [ncut] = { 1, 80,  5, 0.25};
     TString scut[ncut]  = { "nbjet20csvv2m" , "80met", "5met_sqrtHt", "025dphimetjet"};

     TH1F *S_jet2Pt [ncut],    * S_met [ncut],     *S_dphimetjet [ncut],     *S_met_sqrtHt [ncut],    *S_met_meff [ncut],     *S_Ht [ncut];
     TH1F *B_jet2Pt [ncut][3], * B_met [ncut][3],  *B_dphimetjet [ncut][3],  *B_met_sqrtHt [ncut][3], *B_met_meff [ncut][3],  *B_Ht [ncut][3];
    
     TCanvas *CC     [ncut];
     TString massPoint = "Mstop, Mlsp = (350, 225)";
     TString FileAddress = "../../minitreesPT/nominal/Stop/";   


//-----------------------------------------------------
// Set Legend 
//-----------------------------------------------------

TLegend *myLegend ( TString title ){

    TLegend *leg1;
    leg1 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg1->SetFillColor(kWhite); leg1->SetBorderSize(0.);
    leg1->SetTextColor(1); leg1->SetTextSize(0.020);
    leg1->SetTextFont(50);
    //leg1->SetHeader(svar[vv] + "  " + lHeader[R]); //leg1_fit->SetMargin(0.2); 
    leg1->AddEntry((TObject*)0, " ", "");
    
    return leg1; 
     }


//-------------------------------------------------------
// Save Histograms as .root File
// ------------------------------------------------------

void saveHistoRoot ( TString rootfileName, TH1F* histoName){
     
     TFile* OutFileName = new TFile( rootfileName, "recreate");
     histoName   -> Write();
     OutFileName -> Close();
}


//-------------------------------------------------------
// Normalize the Histogram to Unity 
// ------------------------------------------------------

void drawShape (TH1F *histoName){
    float integral;
    integral = histoName -> Integral(0, histoName -> GetNbinsX()+1);
    histoName -> Scale(1/integral);
}


//---------------------------------------------------------
// Plot Signal and Background in a same Canvas
// --------------------------------------------------------

void Plot ( TH1F* S_hName[ncut], TH1F* B_hName[ncut][nbck], TString varName, TString Xtitle){

     TCanvas *CC     [ncut]; 

     for (int c=0; c<ncut; c++){
    
       CC[c]  = new TCanvas(varName + scut[c], "", 1200, 1000);  
       TLegend*  _leg = myLegend ("leg_bck");
       
       for (int i = 0; i < 3; i++){
         int x = 2+i;
         _leg -> AddEntry(B_hName[c][i], LegName [i] , "f");
         B_hName[c][i]->SetLineColor(x); B_hName[c][i]->SetLineWidth(2);
         drawShape(B_hName[c][i]);
        }
      
        _leg -> AddEntry(S_hName[c], massPoint , "f");
        S_hName[c] -> SetLineColor(7); S_hName[c]->SetLineWidth(2);
        drawShape ( S_hName[c]);      
 
//       B_hName[c][0] -> SetMaximum(0.15);
        B_hName[c][0] -> SetTitle (varName + scut[c]);
        B_hName[c][0] -> SetXTitle (Xtitle);
        B_hName[c][0] -> DrawCopy("histo");
        B_hName[c][1] -> DrawCopy("histsame");
        B_hName[c][2] -> DrawCopy("histsame");
        S_hName[c]    -> DrawCopy("histsame");
        _leg ->Draw();

        CC[c] -> Print("shape_study/test/"+ varName + "2" + "_" + scut[c] + ".png");
  
     }
    }


//----------------------------------------------------------------------------------
// Read minitrees, Select the events by Cut and Fill, Save and Plot your Histograms 
// ---------------------------------------------------------------------------------

                  
void shape_study() {

         

     for (int i = 0; i < ncut; i++) {
        S_jet2Pt     [i] = new TH1F ("jet2Pt_Stop_"     + scut[i],  "", 60 , 0. , 300.);
        S_met        [i] = new TH1F ("Met_Stop_"        + scut[i],  "", 80 , 0. , 800.);
        S_dphimetjet [i] = new TH1F ("dphimetjet_Stop_" + scut[i],  "", 80 , 0. , 4.  );
        S_met_meff   [i] = new TH1F ("met_meff_Stop"    + scut[i],  "", 100, 0., 1.  );
        S_met_sqrtHt [i] = new TH1F ("met_sqrtHt_Stop_" + scut[i],  "", 50 , 0. , 15. );

        for (int j = 0; j < 3; j++){
           B_jet2Pt     [i][j] = new TH1F ("jet2Pt_"     + scut[i] + "_"+ BckName [j], "", 60 , 0., 300.);
           B_met        [i][j] = new TH1F ("Met_"        + scut[i] + "_"+ BckName [j], "", 80 , 0., 800.);
           B_dphimetjet [i][j] = new TH1F ("dphimetet_"  + scut[i] + "_"+ BckName [j], "", 80 , 0., 4.  );
           B_met_meff   [i][j] = new TH1F ("met_meff_"   + scut[i] + "_"+ BckName [j], "", 100, 0., 1.  );
           B_met_sqrtHt [i][j] = new TH1F ("met_sqrtHt_" + scut[i] + "_"+ BckName [j], "", 50 , 0. , 15.);
         }
      }




  // Open the minitrees for your signal/ Fill the histo counter/ Save it as .root histogram

 
     TFile *S_MiniTreeFile  = TFile::Open( FileAddress +"/T2tt_mStop-350to400.root"); 
     TTree *S_MiniTree      = GetMiniTree(S_MiniTreeFile);
     Int_t  S_nentries      = (Int_t) S_MiniTree -> GetEntries();

     //for (Int_t i = 0; i<500; i++){
     for (Int_t i = 0; i<S_nentries; i++){
       S_MiniTree -> GetEntry(i);

       if (Mstop != 350 || Mlsp != 225) continue; 
   
       float varCut [ncut] = { nbjet20csvv2m , metPfType1, V4, dphimetjet};   
       htjets = ht - metPfType1 - lep1pt - lep1pt; 
       V3 = metPfType1/meff; 
       V4 = metPfType1/sqrt(htjets);

     /*  bool pass;
       int j; 
       for ( j = 0; j < ncut; j++){ 
         
         if (j==0) { pass =  nbjet20csvv2m >= valCut [0];}
         if (j==1) { pass &= metPfType1    >= valCut [1];}
         if (j==2) { pass &= V4            >= valCut [2];}
         if (j==3) { pass &= dphimetjet    >= valCut [3];}
         if (pass) {S_jet2Pt [j] -> Fill (jet2pt, eventW); S_met [j] -> Fill ( metPfType1, eventW); S_dphimetjet [j] -> Fill (dphimetjet,  eventW); S_met_sqrtHt [j] -> Fill ( V4, eventW);                                     S_met_meff[j] -> Fill (V3, eventW);}
       } 
*/     
       bool pass = true;  
       for (int j = 0; j < ncut; j++){
           pass &= varCut[j] >= valCut [j];
           if (pass) { S_jet2Pt [j] -> Fill (jet2pt, eventW); S_met [j] -> Fill ( metPfType1, eventW); S_dphimetjet [j] -> Fill (dphimetjet,  eventW); S_met_sqrtHt [j] -> Fill ( V4, eventW);                                      S_met_meff[j] -> Fill (V3, eventW);}
       }
     }

   
  // Open the minitrees for your background/ Fill the histo counters/ Save them as .root histogram

    for (int j = 0; j < nbck ; j++){

     TFile *B_MiniTreeFile  = TFile::Open(FileAddress + BckName[j]);
     TTree *B_MiniTree      = GetMiniTree(B_MiniTreeFile);
     Int_t  B_nentries      = (Int_t) B_MiniTree -> GetEntries();

     for (Int_t i = 0; i<B_nentries; i++){

       B_MiniTree -> GetEntry(i);
       
       float varCut [ncut] = { nbjet20csvv2m , metPfType1, V4, dphimetjet};   
       htjets = ht - metPfType1 - lep1pt - lep1pt;
       V3 = metPfType1/meff;
       V4 = metPfType1/sqrt(htjets);

   //    bool pass; 
   //    for (int k = 0; k < ncut; k++ ){
   //       
   //       if (k==0) { pass =  nbjet20csvv2m >= valCut [0];}
   //       if (k==1) { pass &= metPfType1    >= valCut [1];}
   //       if (k==2) { pass &= V4            >= valCut [2];}
   //       if (k==3) { pass &= dphimetjet    >= valCut [3];}


       bool pass = true;
       for (int k = 0; k < ncut; k++){
         pass &= varCut[k] >= valCut [k];
         if (pass) { B_jet2Pt[k][j] -> Fill (jet2pt, eventW); B_met[k][j] -> Fill ( metPfType1, eventW); B_dphimetjet[k][j] -> Fill (dphimetjet,  eventW); B_met_sqrtHt[k][j] -> Fill ( V4, eventW);
                     B_met_meff[k][j] -> Fill (V3, eventW); }
        }
       }
      
    } 

    // TString OutFileName;


     for (int i = 0; i < ncut; i++){ 
         saveHistoRoot ("S_jet2Pt"     + scut[i] + ".root", S_jet2Pt     [i]);  
         saveHistoRoot ("S_met"        + scut[i] + ".root", S_met        [i]);  
         saveHistoRoot ("S_dphimetjet" + scut[i] + ".root", S_dphimetjet [i]);  
         saveHistoRoot ("S_met_meff"   + scut[i] + ".root", S_met_meff   [i]);  
         saveHistoRoot ("S_met_sqrtHt" + scut[i] + ".root", S_met_sqrtHt [i]);  

       for (int j = 0; j < nbck; j++){
          saveHistoRoot ("B_jet2Pt_"     + scut[i] + "_Cut_" +  BckName [j] +".root", B_jet2Pt     [i][j]);
          saveHistoRoot ("B_met_"        + scut[i] + "_Cut_" +  BckName [j] +".root", B_met        [i][j]);
          saveHistoRoot ("B_dphimetjet_" + scut[i] + "_Cut_" +  BckName [j] +".root", B_dphimetjet [i][j]);
          saveHistoRoot ("B_met_meff_"   + scut[i] + "_Cut_" +  BckName [j] +".root", B_met_meff   [i][j]);
          saveHistoRoot ("B_met_sqrtHt"  + scut[i] + "_Cut_" +  BckName [j] +".root", B_met_sqrtHt [i][j]);
        }
      }
    
     Plot ( S_jet2Pt,     B_jet2Pt,     "jet2Pt",     "trailing jet pt [GeV]");
     Plot ( S_met,        B_met,        "met",        "met [GeV]");
     Plot ( S_met_sqrtHt, B_met_sqrtHt, "met_sqrtHt", "met_sqrtHt [GeV]");
     Plot ( S_met_meff,   B_met_meff,   "met_meff",   "met_meff [GeV]");
     Plot ( S_dphimetjet, B_dphimetjet, "dphijetmet", "dphijetmet [rad]");
     
} 
