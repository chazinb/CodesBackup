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

#include "../Razor/Razor.C"
#include "../Razor/SuperRazor.C"


#include <fstream>
#include <iostream>

//----------------------------------------------
// Global Variables
// ---------------------------------------------

     TString FileAddress = "../../minitrees/nominal/Stop/";
     //TString massPoint = "Ms350_Mx225";
     TString  Reg = "R4";

     const int nsample = 7;
     TString SampleName [nsample] = {"16_T2tt_mStop-150to1200.root", "04_TTTo2L2Nu.root", "05_ST.root",  "07_ZJets.root", "06_WW.root", "09_TTV.root", "03_VZ.root"}; 
     TString  _sample   [nsample] = {"T2tt", "TTbar", "ST", "DYJets", "WW", "TTV", "VZ"};
     //----------------------------------------------------------------------------------------------------------------
     // T2tt  = T2tt_mStop-150to1200.root
     // TTbar = TTTo2L2Nu_ext1__part*.root 
     // ST    = ST_tW_antitop.root ST_tW_top.root
     // ZJets = DYJetsToLL_M-10to50.root DYJetsToLL_M-50_000*.root;
     // WW    = WWTo2L2Nu.root GluGluWWTo2L2Nu_MCFM.root
     // TTV   = TTWJetsToLNu.root TTZjets.root
     // VZ    = ZZTo2L2Q*.root WZTo2L2Q__part*.root
     // -----------------------------------------------------------------------------------------------------------------
     
     // const int nsample = 6;
     // TString SampleName [nsample] = {"T2tt_mStop-150to1200.root", "TTTo2L2Nu_ext1.root", "ST_tW.root",  "DYJetsToLL.root", "WWTo2L2Nu.root", "TTZjets.root"}; 
     //----------------------------------------------------------------------------------------------------------------
     // T2tt  = T2tt_mStop-150to1200.root
     // TTbar = TTTo2L2Nu_ext1__part*.root 
     // ST    = ST_tW_antitop.root ST_tW_top.root
     // ZJets = DYJetsToLL_M-10to50.root DYJetsToLL_M-50_000*.root DYJetsToTT_MuEle_M-50.root
     // WW    = WWTo2L2Nu.root 
     // TTZ   = TTZjets.root
     // -----------------------------------------------------------------------------------------------------------------
     // TString  _sample   [nsample] = {"T2tt", "TTbar", "ST", "DYJets", "WW", "TTZ"}; 
        
 
     //enum { nbjet20csvv2m , metPfType1, V4, dphimetjet, ncut};

     // ######################################################################################
     const int ncut = 2;
     //TString scut [ncut] = { "-1nbjetptcsvv2m" , "20nbjetptcsvv2m", "30nbjetptcsvv2m"};
     //TString scut [ncut] = { "-999jet2ptChecktag", "20jet2ptChecktag", "50metChecktag", "20nbjetptcsvv2mChecktag"};
     TString scut [ncut] = { "bVETO", "bSEL"};
     //TString scut [ncut] = { "-999jet2ptCheckMet", "50metCheckMet", "20jet2ptCheckMet","20nbjetptcsvv2mCheckMet"};
     //float valCut [ncut] = { -1, 20, 30};
     float valCut [ncut] = { 20, 20};
     TString svarCut [ncut] = { "bVETO", "bSEL"};
     //float valCut [ncut] = {-999, 50, 20, 20};
     //TString svarCut [ncut] = { "-999jet2ptCheckMet", "50metCheckMet", "20jet2ptCheckMet", "20nbjetptcsvv2mCheckMet"};
     //TString svarCut [ncut] = { "nbjetptcsvv2m" , "nbjetptcsvv2m", "nbjetptcsvv2m"};

     // ################  REmenber Change varCut in your function!!!! ########################

     //float varCut [ncut] = { jet2pt, metPfType1, jet2pt, LeadingPtCSVv2M};
     //float varCut [ncut] = { jet2pt, jet2pt,  metPfType1, LeadingPtCSVv2M};
     //float varCut [ncut] = { LeadingPtCSVv2M , LeadingPtCSVv2M, LeadingPtCSVv2M};
     //float varCut [ncut] = { LeadingPtCSVv2M , jet2pt, metPfType1};

     // ######################################################################################
     enum { ee, mm, em, nchannel}; const TString schannel[nchannel] = {"ee", "mm","em"};
 
     float max_jet1pt, max_jet2pt, max_mlb1, max_mlb2, max_min_mlb, max_max_mlb, max_lep1pt, max_lep2pt, max_lep1phi, max_lep2phi, max_metPfType1, max_dphimetjet, max_met_sqrtHt, max_met_meff, max_Ht, max_Ht_visible, max_MT2ll, max_MT2bb, max_MT2lblb, max_dyll, max_ptbll, max_dphimetptbll, max_m2l, max_mllbb, max_meff, max_dphillmet, max_dphill, max_dphilmet1, max_dphijet1met, max_htjets, max_metPfType1Phi, max_MR, max_R2, max_Rpt, max_invGamma, max_Mdr, max_DeltaPhiRll;    



     TH1F* h_jet1pt         [nchannel][ncut];
     TH1F* h_jet2pt         [nchannel][ncut];
     TH1F* h_mlb1           [nchannel][ncut];    
     TH1F* h_mlb2           [nchannel][ncut];    
     TH1F* h_min_mlb        [nchannel][ncut];    
     TH1F* h_max_mlb        [nchannel][ncut];    
 
     TH1F* h_lep1pt         [nchannel][ncut];
     TH1F* h_lep2pt         [nchannel][ncut];
     TH1F* h_lep1phi        [nchannel][ncut];
     TH1F* h_lep2phi        [nchannel][ncut];

     TH1F* h_metPfType1     [nchannel][ncut]; 
     TH1F* h_dphimetjet     [nchannel][ncut];
     TH1F* h_met_sqrtHt     [nchannel][ncut];
     TH1F* h_met_meff       [nchannel][ncut];
     TH1F* h_Ht             [nchannel][ncut];
     TH1F* h_Ht_visible     [nchannel][ncut];   

     TH1F* h_MT2ll          [nchannel][ncut];
     TH1F* h_MT2bb          [nchannel][ncut];
     TH1F* h_MT2lblb        [nchannel][ncut];
     TH1F* h_dyll           [nchannel][ncut];
     TH1F* h_ptbll          [nchannel][ncut];
     TH1F* h_dphimetptbll   [nchannel][ncut];
     TH1F* h_m2l            [nchannel][ncut];
     TH1F* h_mllbb          [nchannel][ncut];
     TH1F* h_meff           [nchannel][ncut];
     TH1F* h_dphillmet      [nchannel][ncut];
     TH1F* h_dphill         [nchannel][ncut];
     TH1F* h_dphilmet1      [nchannel][ncut];
     TH1F* h_dphijet1met    [nchannel][ncut];
     TH1F* h_htjets         [nchannel][ncut];
     TH1F* h_metPfType1Phi  [nchannel][ncut];

     TH1F* h_MR             [nchannel][ncut];
     TH1F* h_R2             [nchannel][ncut];
     TH1F* h_Rpt            [nchannel][ncut];
     TH1F* h_invGamma       [nchannel][ncut];
     TH1F* h_Mdr            [nchannel][ncut];
     TH1F* h_DeltaPhiRll    [nchannel][ncut];


TString massPoint(float Smass, float Xmass){

   TString sSmass = ""; sSmass += int(Smass); TString sXmass = ""; sXmass += int(Xmass);
   TString masspoint = "_Sm" + sSmass + "_" + "Xm" + sXmass;
   return masspoint; 

}


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

void saveHistoRoot (TString irootfileName, TH1F* histoName){
     TFile* OutFileName = new TFile( irootfileName , "update");    
     histoName   -> Write();
     OutFileName -> Close();
}

void save2dHistoRoot (TString irootfileName, TH2F* histoName){
     TFile* OutFileName = new TFile( irootfileName , "update");
     histoName   -> Write();
     OutFileName -> Close();
}

//-------------------------------------------------------
// Normalize the Histogram to Unity 
// ------------------------------------------------------

void drawShape (TH1F *histoName){
    cout << "n  1" << endl;
    float integral;
    cout << "n  2" << endl;
    integral = histoName -> Integral(0, histoName -> GetNbinsX()+1);
    cout << "n  3" << endl;
    cout  << integral << endl; 
    histoName -> Scale(1/integral);
    cout << "n  4" << endl;
}


void resetHistograms(int ichannel, int icut, TString suffix){

         h_jet1pt         [ichannel][icut]->Reset();
         h_jet2pt         [ichannel][icut]->Reset();
         h_mlb1           [ichannel][icut]->Reset();    
         h_mlb2           [ichannel][icut]->Reset();    
         h_min_mlb        [ichannel][icut]->Reset();
         h_max_mlb        [ichannel][icut]->Reset();
 
         h_lep1pt         [ichannel][icut]->Reset();
         h_lep2pt         [ichannel][icut]->Reset();
         h_lep1phi        [ichannel][icut]->Reset();
         h_lep2phi        [ichannel][icut]->Reset();
     
         h_metPfType1     [ichannel][icut]->Reset();
         h_dphimetjet     [ichannel][icut]->Reset();
         h_met_sqrtHt     [ichannel][icut]->Reset();
         h_met_meff       [ichannel][icut]->Reset();
         h_Ht             [ichannel][icut]->Reset();
         h_Ht_visible     [ichannel][icut]->Reset();

         h_MT2ll          [ichannel][icut]->Reset();           
         h_MT2bb          [ichannel][icut]->Reset();
         h_MT2lblb        [ichannel][icut]->Reset();
         h_dyll           [ichannel][icut]->Reset();
         h_ptbll          [ichannel][icut]->Reset();
         h_dphimetptbll   [ichannel][icut]->Reset();
         h_m2l            [ichannel][icut]->Reset();
         h_mllbb          [ichannel][icut]->Reset();
         h_meff           [ichannel][icut]->Reset();
         h_dphillmet      [ichannel][icut]->Reset();
         h_dphill         [ichannel][icut]->Reset();
         h_dphilmet1      [ichannel][icut]->Reset();
         h_dphijet1met    [ichannel][icut]->Reset();
         h_htjets         [ichannel][icut]->Reset();
         h_metPfType1Phi  [ichannel][icut]->Reset();


         h_MR             [ichannel][icut]->Reset();
         h_R2             [ichannel][icut]->Reset();
         h_Rpt            [ichannel][icut]->Reset();
         h_invGamma       [ichannel][icut]->Reset();
         h_Mdr            [ichannel][icut]->Reset();
         h_DeltaPhiRll    [ichannel][icut]->Reset();


         

}

// --------------------------------
// Define Histograms
// -----------------------------------
void defineHistograms(int ichannel, int icut, TString suffix){

         h_jet1pt         [ichannel][icut] = new TH1F ("h_jet1pt"       + suffix,  "", 45 , 0. , 450.);
         max_jet1pt = h_jet1pt [ichannel][icut] -> GetBinLowEdge(h_jet1pt [ichannel][icut] -> GetNbinsX()+1);         
         
         h_jet2pt         [ichannel][icut] = new TH1F ("h_jet2pt"       + suffix,  "", 45 , 0. , 450.);
         max_jet2pt = h_jet2pt [ichannel][icut] -> GetBinLowEdge(h_jet2pt [ichannel][icut] -> GetNbinsX()+1);         
        
         h_mlb1           [ichannel][icut] = new TH1F ("h_mlb1"           + suffix,  "", 40 , 0. , 400.);
         max_mlb1 = h_mlb1 [ichannel][icut]-> GetBinLowEdge(h_mlb1 [ichannel][icut]-> GetNbinsX()+1);         
         
         h_mlb2           [ichannel][icut] = new TH1F ("h_mlb2"           + suffix,  "", 40 , 0. , 400.);
         max_mlb2 = h_mlb2 [ichannel][icut] -> GetBinLowEdge(h_mlb2 [ichannel][icut] -> GetNbinsX()+1);         
         
         h_min_mlb        [ichannel][icut] = new TH1F ("h_min_mlb"        + suffix,  "", 40 , 0. , 400.);
         max_min_mlb = h_min_mlb [ichannel][icut] -> GetBinLowEdge(h_min_mlb [ichannel][icut] -> GetNbinsX()+1);         
         
         h_max_mlb        [ichannel][icut] = new TH1F ("h_max_mlb"        + suffix,  "", 40 , 0. , 400.);
         max_max_mlb = h_max_mlb [ichannel][icut] -> GetBinLowEdge(h_max_mlb [ichannel][icut] -> GetNbinsX()+1);         



         h_lep1pt         [ichannel][icut] = new TH1F ("h_lep1pt"       + suffix,  "", 45 , 0. , 450.);
         max_lep1pt = h_lep1pt [ichannel][icut] -> GetBinLowEdge(h_lep1pt [ichannel][icut] -> GetNbinsX()+1);         
         
         h_lep2pt         [ichannel][icut] = new TH1F ("h_lep2pt"       + suffix,  "", 45 , 0. , 450.);
         max_lep2pt = h_lep2pt [ichannel][icut] -> GetBinLowEdge(h_lep2pt [ichannel][icut] -> GetNbinsX()+1);         
         
         h_lep1phi        [ichannel][icut] = new TH1F ("h_lep1phi"      + suffix,  "", 200, -3.2, 3.2);       
         max_lep1phi = h_lep1phi [ichannel][icut] -> GetBinLowEdge(h_lep1phi [ichannel][icut] -> GetNbinsX()+1);         
         
         h_lep2phi        [ichannel][icut] = new TH1F ("h_lep2phi"      + suffix,  "", 200, -3.2, 3.2);       
         max_lep2phi = h_lep2phi [ichannel][icut] -> GetBinLowEdge(h_lep2phi [ichannel][icut]-> GetNbinsX()+1);         



         h_metPfType1     [ichannel][icut] = new TH1F ("h_metPfType1"   + suffix,  "", 80 , 0. , 800.);
         max_metPfType1 = h_metPfType1 [ichannel][icut] -> GetBinLowEdge(h_metPfType1 [ichannel][icut] -> GetNbinsX()+1);         

         h_dphimetjet     [ichannel][icut] = new TH1F ("h_dphimetjet"   + suffix,  "", 40 , 0. , 4.  );        
         max_dphimetjet = h_dphimetjet [ichannel][icut] -> GetBinLowEdge(h_dphimetjet [ichannel][icut] -> GetNbinsX()+1);         

         h_met_sqrtHt     [ichannel][icut] = new TH1F ("h_met_sqrtHt"   + suffix,  "", 60 , 0. , 300. );
         max_met_sqrtHt = h_met_sqrtHt [ichannel][icut]-> GetBinLowEdge(h_met_sqrtHt [ichannel][icut] -> GetNbinsX()+1);         

         h_met_meff       [ichannel][icut] = new TH1F ("h_met_meff"     + suffix,  "", 50, 0., 1.  );
         max_met_meff = h_met_meff [ichannel][icut] -> GetBinLowEdge(h_met_meff [ichannel][icut] -> GetNbinsX()+1);         

         h_Ht             [ichannel][icut] = new TH1F ("h_Ht"           + suffix,  "", 100 , 0. , 1000. );
         max_Ht = h_Ht [ichannel][icut] -> GetBinLowEdge(h_Ht [ichannel][icut] -> GetNbinsX()+1);         

         h_Ht_visible     [ichannel][icut] = new TH1F ("h_Ht_visible"   + suffix,  "", 100 , 0. , 1000. );
         max_Ht_visible = h_Ht_visible [ichannel][icut] -> GetBinLowEdge(h_Ht_visible [ichannel][icut] -> GetNbinsX()+1);         



         h_MT2ll          [ichannel][icut] = new TH1F("h_MT2ll"           + suffix, "", 60, 0., 300.);
         max_MT2ll = h_MT2ll [ichannel][icut] -> GetBinLowEdge(h_MT2ll [ichannel][icut] -> GetNbinsX()+1);         
         
         h_MT2bb          [ichannel][icut] = new TH1F("h_MT2bb"           + suffix, "", 80, 0., 800.);
         max_MT2bb = h_MT2bb [ichannel][icut] -> GetBinLowEdge(h_MT2bb [ichannel][icut] -> GetNbinsX()+1);         

         h_MT2lblb        [ichannel][icut] = new TH1F("h_MT2lblb"         + suffix, "", 80, 0., 800.);
         max_MT2lblb = h_MT2lblb [ichannel][icut] -> GetBinLowEdge(h_MT2lblb [ichannel][icut] -> GetNbinsX()+1);         

         h_dyll           [ichannel][icut] = new TH1F("h_dyll"            + suffix, "", 50, 0., 5.);
         max_dyll = h_dyll [ichannel][icut] -> GetBinLowEdge(h_dyll [ichannel][icut] -> GetNbinsX()+1);         

         h_ptbll          [ichannel][icut] = new TH1F("h_ptbll"           + suffix, "", 80, 0., 800.);
         max_ptbll = h_ptbll [ichannel][icut]  -> GetBinLowEdge(h_ptbll [ichannel][icut] -> GetNbinsX()+1);         

         h_dphimetptbll   [ichannel][icut] = new TH1F("h_dphimetptbll"    + suffix, "", 40, 0., 4.);
         max_dphimetptbll = h_dphimetptbll [ichannel][icut] -> GetBinLowEdge(h_dphimetptbll [ichannel][icut] -> GetNbinsX()+1);         

         h_m2l            [ichannel][icut] = new TH1F("h_m2l"             + suffix, "", 80, 0., 800.);
         max_m2l = h_m2l [ichannel][icut] -> GetBinLowEdge(h_m2l [ichannel][icut] -> GetNbinsX()+1);         

         h_mllbb          [ichannel][icut] = new TH1F("h_mllbb"           + suffix, "", 80, 0., 800.);
         max_mllbb = h_mllbb [ichannel][icut] -> GetBinLowEdge(h_mllbb [ichannel][icut] -> GetNbinsX()+1);         

         h_meff           [ichannel][icut] = new TH1F("h_meff"            + suffix, "", 75, 0., 1500.);
         max_meff = h_meff [ichannel][icut] -> GetBinLowEdge(h_meff [ichannel][icut] -> GetNbinsX()+1);         

         h_dphillmet      [ichannel][icut] = new TH1F("h_dphillmet"       + suffix, "", 40, 0., 4.);
         max_dphillmet = h_dphillmet [ichannel][icut] -> GetBinLowEdge(h_dphillmet [ichannel][icut] -> GetNbinsX()+1);         

         h_dphill         [ichannel][icut] = new TH1F("h_dphill"          + suffix, "", 40, 0., 4.);
         max_dphill = h_dphill [ichannel][icut] -> GetBinLowEdge(h_dphill [ichannel][icut] -> GetNbinsX()+1);         

         h_dphilmet1      [ichannel][icut] = new TH1F("h_dphilmet1"       + suffix, "", 40, 0., 4.);
         max_dphilmet1 = h_dphilmet1 [ichannel][icut] -> GetBinLowEdge(h_dphilmet1 [ichannel][icut] -> GetNbinsX()+1);         

         h_dphijet1met    [ichannel][icut] = new TH1F("h_dphijet1met"     + suffix, "", 40, 0., 4.);
         max_dphijet1met = h_dphijet1met [ichannel][icut] -> GetBinLowEdge(h_dphijet1met [ichannel][icut] -> GetNbinsX()+1);         

         h_htjets         [ichannel][icut] = new TH1F("h_htjets"          + suffix, "", 55, 0., 1100.);
         max_htjets = h_htjets [ichannel][icut] -> GetBinLowEdge(h_htjets [ichannel][icut] -> GetNbinsX()+1);         

         h_metPfType1Phi  [ichannel][icut] = new TH1F("h_metPfType1Phi"   + suffix, "", 40, 0., 4.);
         max_metPfType1Phi = h_metPfType1Phi [ichannel][icut] -> GetBinLowEdge(h_metPfType1Phi [ichannel][icut] -> GetNbinsX()+1);         



         h_MR             [ichannel][icut] = new TH1F("h_MR"              + suffix, "", 100, 0., 1000.);
         max_MR = h_MR [ichannel][icut] -> GetBinLowEdge(h_MR [ichannel][icut] -> GetNbinsX()+1);         

         h_R2             [ichannel][icut] = new TH1F("h_R2"              + suffix, "", 28, 0., 1.4);
         max_R2 = h_R2 [ichannel][icut] -> GetBinLowEdge(h_R2 [ichannel][icut] -> GetNbinsX()+1);         

         h_Rpt            [ichannel][icut] = new TH1F("h_Rpt"             + suffix, "", 22, 0., 1.1);
         max_Rpt = h_Rpt [ichannel][icut] -> GetBinLowEdge(h_Rpt [ichannel][icut] -> GetNbinsX()+1);         

         h_invGamma       [ichannel][icut] = new TH1F("h_invGamma"        + suffix, "", 15, 0., 1.5);
         max_invGamma = h_invGamma [ichannel][icut] -> GetBinLowEdge(h_invGamma [ichannel][icut] -> GetNbinsX()+1);         

         h_Mdr            [ichannel][icut] = new TH1F("h_Mdr"             + suffix, "", 45, 0., 450.);
         max_Mdr = h_Mdr [ichannel][icut] -> GetBinLowEdge(h_Mdr [ichannel][icut] -> GetNbinsX()+1);         

         h_DeltaPhiRll    [ichannel][icut] = new TH1F("h_DeltaPhiRll"     + suffix, "", 35, 0, 3.5);
         max_DeltaPhiRll = h_DeltaPhiRll [ichannel][icut] -> GetBinLowEdge(h_DeltaPhiRll [ichannel][icut] -> GetNbinsX()+1);         

            
}





//----------------------------------------------------------------------------------
// Cut Study
// ---------------------------------------------------------------------------------



                  
void cut_study(float Smass, float Xmass) {

  TString massP;

  if (Smass == 1. && Xmass == 1.){
       massP = Reg;
      }else{
       massP = massPoint(Smass, Xmass);
      }

  cout << " rm -rf irootfiles    ? " << endl; 

  //                   //
  // Define Histograms //
  //                   //

 
  for (int icut = 0; icut < ncut; icut++){
     for (int ichannel = 0; ichannel < nchannel; ichannel++) {

       TString suffix = "_" + schannel[ichannel] + "_" + scut[icut];
       defineHistograms(ichannel, icut, suffix);
     }
  }


  //      //
  // Loop //
  //      //

  for (int s = 0; s < 1 ; s++){
  //for (int s = 0; s < nsample ; s++){

     if (_sample[s] == "T2tt") _sample[s] += massP;
     gSystem->mkdir("irootfile/" + _sample[s] + "/", kTRUE);

 
     TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]); 
     TTree *MiniTree      = GetMiniTree(MiniTreeFile);
     Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

     for (int icut = 0; icut < ncut; icut++){
        for (int ichannel = 0; ichannel < nchannel; ichannel++) {

          TString suffix = "_" + schannel[ichannel] + "_" + scut[icut];
          resetHistograms(ichannel, icut, suffix);
       }
     }
     
     //for (Int_t en = 0; en<1000; en++){
     for (Int_t en = 0; en<nentries; en++){
       MiniTree -> GetEntry(en);

       eventW *= 36.;
       bool stopsample = true;

       if (_sample[s].Contains("T2tt")) {
        if (Smass==1. && Xmass==1.) {
          if ( Mstop < 240 || Mstop - Mlsp -Mtop > -10) stopsample = false;
        } else {
          if ( Mstop != Smass || Mlsp != Xmass) stopsample = false;
        }
       }


       if (stopsample == false) continue;

       //                           //
       // Met Variable Calculations //
       //                           //

       float HtJets = ht - metPfType1 - lep1pt - lep1pt;
       V3 = metPfType1/meff;
       V4 = metPfType1/sqrt(HtJets);
       V5 = ht - metPfType1;
       V6 = min(mlb1,mlb2);
       V7 = max(mlb1,mlb2); 
      
       //                    //
       // Razor Calculations //
       //                    //
       

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


       // Block the overflow //
       // ---------------------------------------------------------------------------
       
       if (jet1pt       >= max_jet1pt)         jet1pt        = max_jet1pt       -0.01;
       if (jet2pt       >= max_jet2pt)         jet2pt        = max_jet2pt       -0.01; 
       if (mlb1         >= max_mlb1)           mlb1          = max_mlb1         -0.01;
       if (mlb2         >= max_mlb2)           mlb2          = max_mlb2         -0.01;
       if (V6           >= max_min_mlb)        V6            = max_min_mlb      -0.01;
       if (V7           >= max_max_mlb)        V7            = max_max_mlb      -0.01;
       if (lep1pt       >= max_lep1pt)         lep1pt        = max_lep1pt       -0.01;
       if (lep2pt       >= max_lep2pt)         lep2pt        = max_lep2pt       -0.01;
       if (lep1phi      >= max_lep1phi)        lep1phi       = max_lep1phi      -0.01; 
       if (lep2phi      >= max_lep2phi)        lep2phi       = max_lep2phi      -0.01;
       if (metPfType1   >= max_metPfType1)     metPfType1    = max_metPfType1   -0.01;
       if (dphimetjet   >= max_dphimetjet)     dphimetjet    = max_dphimetjet   -0.01;
       if (V4           >= max_met_sqrtHt)     V4            = max_met_sqrtHt   -0.01;
       if (V3           >= max_met_meff)       V3            = max_met_meff     -0.01;
       if (ht           >= max_Ht)             ht            = max_Ht           -0.01;
       if (V5           >= max_Ht_visible)     V5            = max_Ht_visible   -0.01; 
       if (mt2ll        >= max_MT2ll)          mt2ll         = max_MT2ll        -0.01;
       if (mt2bb        >= max_MT2bb)          mt2bb         = max_MT2bb        -0.01;
       if (mt2lblb      >= max_MT2lblb)        mt2lblb       = max_MT2lblb      -0.01; 
       if (dyll         >= max_dyll)           dyll          = max_dyll         -0.01;
       if (ptbll        >= max_ptbll)          ptbll         = max_ptbll        -0.01;
       if (dphimetptbll >= max_dphimetptbll)  dphimetptbll  = max_dphimetptbll -0.01;
       if (m2l          >= max_m2l)            m2l           = max_m2l          -0.01;
       if (mllbb        >= max_mllbb)          mllbb         = max_mllbb        -0.01;
       if (meff         >= max_meff)           meff          = max_meff         -0.01;
       if (dphillmet    >= max_dphillmet)      dphillmet     = max_dphillmet    -0.01;
       if (dphillmet    >= max_dphill)         dphillmet     = max_dphill       -0.01;
       if (dphilmet1    >= max_dphilmet1)      dphilmet1     = max_dphilmet1    -0.01;
       if (dphijet1met  >= max_dphijet1met)    dphijet1met   = max_dphijet1met  -0.01;
       if (htjets       >= max_htjets)         htjets        = max_htjets       -0.01;
       if (metPfType1Phi>= max_metPfType1Phi)  metPfType1Phi = max_metPfType1Phi-0.01;
       if (MR           >= max_MR)             MR            = max_MR           -0.01;
       if (R2           >= max_R2)             R2            = max_R2           -0.01;
       if (Rpt          >= max_Rpt)            Rpt           = max_Rpt          -0.01;
       if (invGamma     >= max_invGamma)       invGamma      = max_invGamma     -0.01;
       if (Mdr          >= max_Mdr)            Mdr           = max_Mdr          -0.01;
       if (DeltaPhiRll  >= max_DeltaPhiRll)    DeltaPhiRll   = max_DeltaPhiRll  -0.01;
         
      //-----------------------------------------------------------------------------






       //               //
       // varCut [ncut] //
       //               //

       // ####################  Change varCut !!!!!! ################################## 

       float varCut [ncut] = {LeadingPtCSVv2M, LeadingPtCSVv2M};
       //float varCut [ncut] = { jet2pt, metPfType1, jet2pt, LeadingPtCSVv2M};
       //float varCut [ncut] = { jet2pt, jet2pt,  metPfType1, LeadingPtCSVv2M};
       //float varCut [ncut] = { LeadingPtCSVv2M , LeadingPtCSVv2M, LeadingPtCSVv2M};
       //float varCut [ncut] = { LeadingPtCSVv2M , jet2pt, metPfType1};

       // #############################################################################

       //         //
       // Channel //
       //         //


        int i=-1;

        if (channel == 3) { i = ee; }
        else if (channel == 4) { i = mm; }
        else if (channel == 5) { i = em; }
            
	if (i==-1) { cout << "Warning! # l < 2 " <<  channel <<  endl; continue;}	

       //              //
       // Cut and Fill //
       //              //
 
        bool pass = true;  
        for (int j = 0; j < ncut; j++){

         //  if (channel!=5 || scut[j]!="nbjetptcsvv2m")//////
                //pass = varCut[j] >= valCut [j];
           if (j == 0) pass = varCut[j] < valCut [j];
           if (j == 1) pass = varCut[j] > valCut [j];
           //if (channel == 5 && varCut[j] == nbjet20csvv2m) continue;
           if (pass) { h_jet1pt         [i][j] -> Fill (jet1pt,       eventW);
                       h_jet2pt         [i][j] -> Fill (jet2pt,       eventW);

                       h_mlb1           [i][j] -> Fill (mlb1,       eventW);
                       h_mlb2           [i][j] -> Fill (mlb2,       eventW);
                       h_min_mlb        [i][j] -> Fill (V6,       eventW);
                       h_max_mlb        [i][j] -> Fill (V7,       eventW);

                       h_lep1pt         [i][j] -> Fill (lep1pt,       eventW);                       
                       h_lep2pt         [i][j] -> Fill (lep2pt,       eventW);                       
                       h_lep1phi        [i][j] -> Fill (lep1phi,      eventW);
                       h_lep2phi        [i][j] -> Fill (lep2phi,      eventW);
       
                                      
                       h_metPfType1     [i][j] -> Fill ( metPfType1,  eventW); 
                       h_dphimetjet     [i][j] -> Fill (dphimetjet,   eventW); 
                       h_met_sqrtHt     [i][j] -> Fill (V4,           eventW); 
                       h_met_meff       [i][j] -> Fill (V3,           eventW); 
                       h_Ht             [i][j] -> Fill (ht,           eventW);
                       h_Ht_visible     [i][j] -> Fill (V5,           eventW);

                       h_MT2ll          [i][j] -> Fill(mt2ll ,        eventW);
                       h_MT2bb          [i][j] -> Fill(mt2bb ,        eventW);
                       h_MT2lblb        [i][j] -> Fill(mt2lblb ,      eventW);
                       h_dyll           [i][j] -> Fill(dyll ,         eventW);
                       h_ptbll          [i][j] -> Fill(ptbll ,        eventW);
                       h_dphimetptbll   [i][j] -> Fill(dphimetptbll,   eventW);
                       h_m2l            [i][j] -> Fill(m2l ,          eventW);
                       h_mllbb          [i][j] -> Fill(mllbb ,        eventW);
                       h_meff           [i][j] -> Fill(meff ,         eventW);
                       h_dphillmet      [i][j] -> Fill(dphillmet ,    eventW);
                       h_dphill         [i][j] -> Fill(dphill,        eventW);
                       h_dphilmet1      [i][j] -> Fill(dphilmet1,     eventW);
                       h_dphijet1met    [i][j] -> Fill(dphijet1met,   eventW);
                       h_htjets         [i][j] -> Fill(htjets,        eventW);
                       h_metPfType1Phi  [i][j] -> Fill(metPfType1Phi, eventW);

                       h_MR             [i][j] -> Fill(MR,            eventW);
                       h_R2             [i][j] -> Fill(R2,            eventW);
                       h_Rpt            [i][j] -> Fill(Rpt,           eventW);
                       h_invGamma       [i][j] -> Fill(invGamma,      eventW);
                       h_Mdr            [i][j] -> Fill(Mdr,           eventW);
                       h_DeltaPhiRll    [i][j] -> Fill(DeltaPhiRll,   eventW);



           }
        }

      }         
         
          for (int ch = ee; ch < nchannel; ch ++){
            for (int cut = 0; cut < ncut; cut++){

  
              gSystem->mkdir("irootfile/" + _sample[s] + "/" + scut[cut] + "/", kTRUE);

              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_jet1pt"         + ".root", h_jet1pt     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_jet2pt"         + ".root", h_jet2pt     [ch][cut]);

              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_lep1pt"         + ".root", h_lep1pt      [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_lep2pt"         + ".root", h_lep2pt      [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_lep1phi"        + ".root", h_lep1phi     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_lep2phi"        + ".root", h_lep2phi     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_mlb1"           + ".root", h_mlb1        [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_mlb2"           + ".root", h_mlb2        [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_min_mlb"        + ".root", h_min_mlb     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_max_mlb"        + ".root", h_max_mlb     [ch][cut]);

              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_metPfType1"     + ".root", h_metPfType1 [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_dphimetjet"     + ".root", h_dphimetjet [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_met_meff"       + ".root", h_met_meff   [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_met_sqrtHt"     + ".root", h_met_sqrtHt [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_Ht"             + ".root", h_Ht         [ch][cut]); 
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_Ht_visible"     + ".root", h_Ht_visible [ch][cut]); 
    
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_MT2ll"         + ".root", h_MT2ll         [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_MT2bb"         + ".root", h_MT2bb         [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_MT2lblb"       + ".root", h_MT2lblb       [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dyll"          + ".root", h_dyll          [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_ptbll"         + ".root", h_ptbll         [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dphimetptbll"  + ".root", h_dphimetptbll  [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_m2l"           + ".root", h_m2l           [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_mllbb"         + ".root", h_mllbb         [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_meff"          + ".root", h_meff          [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dphillmet"     + ".root", h_dphillmet     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dphill"        + ".root", h_dphill        [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dphilmet1"     + ".root", h_dphilmet1     [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_dphijet1met"   + ".root", h_dphijet1met   [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_htjets"        + ".root", h_htjets        [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_metPfType1Phi" + ".root", h_metPfType1Phi [ch][cut]);
     
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_MR"            + ".root", h_MR            [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_R2"            + ".root", h_R2            [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_Rpt"           + ".root", h_Rpt           [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_invGamma"      + ".root", h_invGamma      [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_Mdr"           + ".root", h_Mdr           [ch][cut]);
              saveHistoRoot ("irootfile/" + _sample[s] + "/" + scut[cut] + "/" +  "h_DeltaPhiRll"   + ".root", h_DeltaPhiRll   [ch][cut]);


            }
         }
   

    }
}





//------------------------------------------------------------------------------------------------------------------------
// Shape Study
//-----------------------------------------------------------------------------------------------------------------------



void shape_study(float Smass, float Xmass){

      TString massP;
      if (Smass == 1. && Xmass == 1.){
         massP = Reg;
      }else{
         massP = massPoint(Smass, Xmass);
      }
 
     cout << " rm -rf irootfiles    ? " << endl;

     //                   //
     // Define Histograms //
     //                   //

      for (int ichannel = 0; ichannel < nchannel; ichannel++) {
      
        cout <<1<<endl;  
        const int icut = 0;       
        TString suffix = "_" + schannel[ichannel] ;
        defineHistograms(ichannel, icut, suffix);
     
       }

     
      //      //
      // Loop // 
      //      //

      for (int s = 0; s < nsample ; s++){
      cout << 2<< endl;
      if (_sample[s] == "T2tt") _sample[s] += massP; // if (_sample[s].Contains("T2tt") { _sample[s] = "T2tt" + massP;}
      gSystem->mkdir("irootfile/" + _sample[s] + "/shape/", kTRUE);

      TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]);
      TTree *MiniTree      = GetMiniTree(MiniTreeFile);
      Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

      for (int ichannel = 0; ichannel < nchannel; ichannel++) {
      cout<<3<<endl;
        const int icut = 0;
        TString suffix = "_" + schannel[ichannel] + "_";
      cout << 3.1 << endl;  
      resetHistograms(ichannel, icut, suffix);
      cout << 3.2 << endl;
      }

      //for (Int_t en = 0; en<500; en++){
      for (Int_t en = 0; en<nentries; en++){
        MiniTree -> GetEntry(en);
       cout << 4<< endl; 
        eventW *= 36.; 
        bool stopsample = true;

        if (_sample[s].Contains("T2tt")) {
        if (Smass==1. && Xmass==1.) {
          if ( Mstop < 240 || Mstop - Mlsp -Mtop > -10) stopsample = false;
        } else {
          if ( Mstop != Smass || Mlsp != Xmass) stopsample = false;
        }
       }

        if (stopsample == false) continue;

        //                           //
        // Met variable Calculations //
        //                           //

         float HtJets = ht - metPfType1 - lep1pt - lep1pt;
         V3 = metPfType1/meff;
         V4 = metPfType1/sqrt(HtJets);
         V5 = ht - metPfType1;
         V6 = min(mlb1,mlb2);
         V7 = max(mlb1,mlb2);

 
        //                    //
        // Razor Calculations //
        //                    //

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

        //         //
        // Channel //
        //         //

        //const int icut = 0;
        const int j = 0;
        int i=-1;
        if (channel == 3) { i = ee; }
        else if (channel == 4) { i = mm; }
        else if (channel == 5) { i = em; }


        //      //
        // Fill //
        //      //


        h_jet1pt         [i][j] -> Fill (jet1pt,           eventW);
        h_jet2pt         [i][j] -> Fill (jet2pt,           eventW);
 
        h_mlb1           [i][j] -> Fill (mlb1,       eventW);
        h_mlb2           [i][j] -> Fill (mlb2,       eventW);
        h_min_mlb        [i][j] -> Fill (V6,         eventW);
        h_max_mlb        [i][j] -> Fill (V7,         eventW);

        h_lep1pt         [i][j] -> Fill (lep1pt,       eventW);   
        h_lep2pt         [i][j] -> Fill (lep2pt,       eventW);
        h_lep1phi        [i][j] -> Fill (lep1phi,      eventW);  
        h_lep2phi        [i][j] -> Fill (lep2phi,      eventW);  

        h_metPfType1     [i][j] -> Fill (metPfType1,       eventW);
        h_dphimetjet     [i][j] -> Fill (dphimetjet,       eventW);
        h_met_sqrtHt     [i][j] -> Fill (V4,               eventW);
        h_met_meff       [i][j] -> Fill (V3,               eventW);
        h_Ht             [i][j] -> Fill (ht,               eventW);
        h_Ht_visible     [i][j] -> Fill (V5,           eventW);

        h_MT2ll          [i][j] -> Fill(mt2ll ,        eventW);
        h_MT2bb          [i][j] -> Fill(mt2bb ,        eventW);
        h_MT2lblb        [i][j] -> Fill(mt2lblb ,      eventW);
        h_dyll           [i][j] -> Fill(dyll ,         eventW);
        h_ptbll          [i][j] -> Fill(ptbll ,        eventW);
        h_dphimetptbll   [i][j] -> Fill(dphimetptbll,  eventW);
        h_m2l            [i][j] -> Fill(m2l ,          eventW);
        h_mllbb          [i][j] -> Fill(mllbb ,        eventW);
        h_meff           [i][j] -> Fill(meff ,         eventW);
        h_dphillmet      [i][j] -> Fill(dphillmet ,    eventW);
        h_dphill         [i][j] -> Fill(dphill,        eventW);
        h_dphilmet1      [i][j] -> Fill(dphilmet1,     eventW);
        h_dphijet1met    [i][j] -> Fill(dphijet1met,   eventW);
        h_htjets         [i][j] -> Fill(htjets,        eventW);
        h_metPfType1Phi  [i][j] -> Fill(metPfType1Phi, eventW);

        h_MR             [i][j] -> Fill(MR,               eventW);
        h_R2             [i][j] -> Fill(R2,               eventW);
        h_Rpt            [i][j] -> Fill(Rpt,              eventW);
        h_invGamma       [i][j] -> Fill(invGamma,         eventW);
        h_Mdr            [i][j] -> Fill(Mdr,              eventW);
        h_DeltaPhiRll    [i][j] -> Fill(DeltaPhiRll,      eventW);
 

      
       }
 
       //      //
       // Save //
       //      //
        
       for (int ch = ee; ch < nchannel; ch ++){
  
         const int icut = 0;
         
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_jet1pt"        + ".root", h_jet1pt        [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_jet2pt"        + ".root", h_jet2pt        [ch][icut]);       
     
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_lep1pt"         + ".root", h_lep1pt        [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_lep2pt"         + ".root", h_lep2pt        [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_lep1phi"        + ".root", h_lep1phi       [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_lep2phi"        + ".root", h_lep2phi       [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_mlb1"           + ".root", h_mlb1          [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_mlb2"           + ".root", h_mlb2          [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_min_mlb"        + ".root", h_min_mlb       [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" + "h_max_mlb"        + ".root", h_max_mlb       [ch][icut]);
     

         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_metPfType1"    + ".root", h_metPfType1    [ch][icut]);      
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphimetjet"    + ".root", h_dphimetjet    [ch][icut]);      
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_met_sqrtHt"    + ".root", h_met_sqrtHt    [ch][icut]);      
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_met_meff"      + ".root", h_met_meff      [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_Ht"            + ".root", h_Ht            [ch][icut]);
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_Ht_visible"    + ".root", h_Ht_visible    [ch][icut]);
       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_MT2ll"         + ".root", h_MT2ll         [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_MT2bb"         + ".root", h_MT2bb         [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_MT2lblb"       + ".root", h_MT2lblb       [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dyll"          + ".root", h_dyll          [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_ptbll"         + ".root", h_ptbll         [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphimetptbll"  + ".root", h_dphimetptbll  [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_m2l"           + ".root", h_m2l           [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_mllbb"         + ".root", h_mllbb         [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_meff"          + ".root", h_meff          [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphillmet"     + ".root", h_dphillmet     [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphill"        + ".root", h_dphill        [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphilmet1"     + ".root", h_dphilmet1     [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_dphijet1met"   + ".root", h_dphijet1met   [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_htjets"        + ".root", h_htjets        [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_metPfType1Phi" + ".root", h_metPfType1Phi [ch][icut]);
       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_MR"            + ".root", h_MR            [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_R2"            + ".root", h_R2            [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_Rpt"           + ".root", h_Rpt           [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_invGamma"      + ".root", h_invGamma      [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_Mdr"           + ".root", h_Mdr           [ch][icut]);       
         saveHistoRoot ("irootfile/" + _sample[s] + "/shape/" +  "h_DeltaPhiRll"   + ".root", h_DeltaPhiRll   [ch][icut]);       
         
       }

}
}





//---------------------------------------------------------------------------------------
// Ingredients for Combined Higgs Tool -> Made for 1 variable, for 1 final cut selection
// --------------------------------------------------------------------------------------




void Ingredients(float Smass, float Xmass, TString varName, TString cutName){

	TString massP = massPoint(Smass, Xmass);

        gSystem->mkdir("datacard/T2tt" + massP + "/", kTRUE);
   
        // T2tt_massP_shape_varName_THF1.root creation 

 
         TH1F *_histo[nchannel][nsample];

        for (int s = 0; s < nsample; s++){ 
        
           if (_sample[s].Contains("T2tt")) { _sample[s] = "T2tt" + massP;} 
           TFile* _H = new TFile( "irootfile/" + _sample[s] + "/" + cutName + "/" + "h_" + varName  +".root","read");
           for (int ch = 0; ch<nchannel; ch++) {  
             _histo[ch][s] = (TH1F*)_H-> Get ("h_" + varName + "_" + schannel[ch] + "_" + cutName);
             _histo[ch][s] -> SetName(_sample[s] + "_" + schannel[ch]);
            cout << _histo[ch][s] -> GetNbinsX() << endl; 
             //_histo[ch][s] -> Rebin(4);
            cout << _histo[ch][s] -> GetNbinsX() << endl;
//             _histo[ch][s] -> GetXaxis() -> SetRangeUser(0, 140)
           }
        }
    
       TFile* OutFileName = new TFile( "datacard/T2tt" + massP + "/" + "T2tt_" + massP + "_shapes_Met50_" + varName + "_TH1F.root" , "recreate");
       for (int ch = 0; ch < nchannel; ch ++){ 
          for (int s = 0; s < nsample; s++){ _histo[ch][s] -> Write();}
           _histo[ch][1]->SetName("data_obs_" + schannel[ch]);
           _histo[ch][1]-> Write();  
        }
        OutFileName-> Close();

       
        // Datacard .txt writting
      
       std::ofstream inFile("datacard/T2tt" + massP + "/" + "T2tt_" + massP + "_shapes_Met50_" + varName + ".txt",std::ios::out);

       inFile << "\n";
       inFile << "   imax  "<<"     " << nchannel   << endl;
       inFile << "   jmax "<<"     " << nsample -1 << endl;
       inFile << "   kmax "<<"     " <<  "*"       << endl;
       inFile << "-----------------------------------------------------------------------------------" << endl;
       inFile << "  shapes * * T2tt_" +  massP + "_shapes_Met50_" + varName + "_TH1F.root  $PROCESS_$CHANNEL" << endl;
       inFile << "-----------------------------------------------------------------------------------" << endl;
       inFile << "   bin  "<<"      "<< "     ee          "<<"     mm             "<<"      em             "<<endl;
       inFile << "observation           "; 
       for (int ch = 0; ch < nchannel; ch++){ float NBck = _histo[ch][1] -> Integral(0, _histo[ch][1] -> GetNbinsX()+1); 
       inFile <<"" <<  NBck << "         ";
       }
       inFile << ""<< endl;
       inFile << "-----------------------------------------------------------------------------------" << endl;
       inFile << "   bin  "<<"        "; 
       for (int p = 0; p < nsample; p++){ 
                                inFile << "   ee        "<<"     mm         "<<"      em      ";
       }
       inFile << " " <<endl;
       inFile <<" process "<<"           "; 
       for (int s = 0; s < nsample; s++){ for (int ch = 0; ch < nchannel; ch++){ 
                                 inFile << " " << s  << "             ";
         }
       }
       inFile << " " <<endl;
       inFile <<" process              ";
       for (int s = 0; s < nsample; s++){
         for(int ch = 0; ch < nchannel; ch++){
                                 inFile << "   " << _sample[s] << "   ";  
         }
       }
       inFile << "" <<endl;
       inFile << " rate        ";
       for (int s = 0; s < nsample; s++){ for (int ch = 0; ch < nchannel; ch++){float integral = _histo[ch][s] -> Integral(0, _histo[ch][s] -> GetNbinsX()+1);
       inFile << integral<<"      ";}}      
       inFile.close();
}






//---------------------------------------------------------
// Plot Cut in Signal and Background in a same Canvas
// --------------------------------------------------------




void Plot_cut(float Smass, float Xmass, TString varName, TString xtitle){

      TString massP;
      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }
     
 
      TH1F* histo[nsample][nchannel][ncut];
      TFile* Drawhisto[nsample][ncut];     

      for (int s = 0; s < nsample; s++){
        for (int cut = 0; cut < ncut;  cut++){
          if (_sample[s] == "T2tt"){
            if (Smass==1. && Xmass==1.) { _sample[s] += Reg;
             } else {
              _sample[s] += massP;
             }
           }
          
           Drawhisto[s][cut] = new TFile( "irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_" + varName  +".root","read");
        }
      }
  
      for (int cut = 0; cut < ncut; cut++){
          gSystem->mkdir("iplots/" + massP + "/" + scut[cut] + "/", kTRUE);
        
        for (int ch = 0; ch < nchannel; ch++){
           
          TCanvas* CC; 
          CC = new TCanvas (varName +  "_" + schannel[ch] + "_" + scut[cut], "", 1200,1000);
          gStyle->SetOptStat(""); 
          TLegend*  _leg = myLegend ("");  

          float MaxPlot = 0.001; 
          for (int s = 0; s< nsample; s++){
             cout << _sample[s] << endl;
             histo[s][ch][cut] = (TH1F*)Drawhisto[s][cut]-> Get ("h_" + varName + "_" + schannel[ch] + "_" + scut[cut]);
             cout << "h_" + varName + "_" + schannel[ch] + "_" + scut[cut] << endl;  
             _leg -> AddEntry(histo[s][ch][cut], _sample[s]  , "f");
             drawShape(histo[s][ch][cut]);
             if (histo[s][ch][cut]->GetMaximum()>MaxPlot) {
               MaxPlot = histo[s][ch][cut]->GetMaximum();
             }

          }
     
 
          histo[0][ch][cut] -> SetLineColor(2); histo[0][ch][cut]->SetLineWidth(2);         
          histo[0][ch][cut] -> SetMaximum(MaxPlot*1.1);
          histo[0][ch][cut] -> SetTitle (varName);
          histo[0][ch][cut] -> SetXTitle (xtitle); 
          for (int sma = 1; sma < nsample; sma++){histo[sma][ch][cut] -> SetLineColor(2+sma); histo[sma][ch][cut]->SetLineWidth(2);}
          
          histo[0][ch][cut] -> Draw("histo");  
          for (int sma = 1; sma < nsample; sma++){histo [sma][ch][cut]-> DrawCopy("histsame");}
          _leg ->Draw();
          
          CC -> Print("iplots/" + massP + "/" + scut[cut] + "/" +  varName + "_" + schannel[ch] + ".png");
         }
       }

       gSystem->Exec("cp  ./index.php ./iplots/; done");
       gSystem->Exec("for dir in $(find ./iplots/ -type d); do cp -n ./index.php $dir/; done");
}
 





//--------------------------------------------------------------
// Plot Shape and Give the yield Table ( Signal and Background)
// -------------------------------------------------------------




void Plot_shape(float Smass, float Xmass, TString varName, TString xtitle){

      TString massP;
      
      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }

      TH1F* histo[nsample][nchannel];
      TFile* Drawhisto[nsample];
      const int icut = 0;

      for (int s = 0; s < nsample; s++){
          
       if (_sample[s] == "T2tt"){
            if (Smass==1. && Xmass==1.) { _sample[s] += Reg;
             } else {
              _sample[s] += massP;
             }
           }

        Drawhisto[s] = new TFile( "irootfile/" + _sample[s] + "/shape/" + "h_" + varName  +".root","read");
      }

        gSystem->mkdir("iplots/shape/"+ massP +"/" , kTRUE);
        for (int ch = 0; ch < nchannel; ch++){

          TCanvas* CC;
          CC = new TCanvas (varName +  "_" + schannel[ch], "", 1200,1000);
          gStyle->SetOptStat("");
          TLegend*  _leg = myLegend ("");

          
          for (int s = 0; s< nsample; s++){ histo[s][ch] = (TH1F*)Drawhisto[s]-> Get ("h_" + varName + "_" + schannel[ch]);}
             
 
             // Table of Yields //
            
        
             cout << "                                                                         " << endl;
             cout << "------------------------------------------------------------------------ " << endl;
             cout << "                                                                         " << endl;
             cout << " channel  " << schannel[ch] << endl;
             cout << "------------------------------------------------------------------------ " << endl;
             cout << "                                                                         " << endl;
             cout << " Yields  "  << "    #  events   " << "    S/sqrt(S+B) " << "      S/B    " <<endl;
             cout << "                                                                         " << endl;
             cout << "------------------------------------------------------------------------ " << endl;

             float SignalYields = histo[0][ch] -> Integral(0, histo[0][ch] -> GetNbinsX()+1);

             cout << _sample[0] << "       " <<  SignalYields << endl;

             float BackgroundYield = 0.;
             
             for (int b = 1; b < nsample; b++) {
             float ThisBckYield = histo[b][ch] -> Integral(0, histo[b][ch]-> GetNbinsX()+1);
             float Signb = SignalYields / sqrt(SignalYields + ThisBckYield);
             float Signa = SignalYields / ThisBckYield;
             BackgroundYield += ThisBckYield;
             cout <<  _sample[b] << "            " << ThisBckYield << "      " << Signb << "          "  << Signa << endl;
             cout << "                                                                                          " << endl;

             }

             float Signtb = SignalYields / sqrt(SignalYields + BackgroundYield); 
             float Signta = SignalYields / sqrt(SignalYields + BackgroundYield);
             cout << "ToTal MC              " << BackgroundYield << "          " <<Signtb << "         "  << Signta <<endl;
 
             // End Table of Yields //

          float MaxPlot = 0.001;
          for (int s = 0; s< nsample; s++){ 
             _leg -> AddEntry(histo[s][ch], _sample[s]  , "f"); drawShape(histo[s][ch]);
             if (histo[s][ch]->GetMaximum()>MaxPlot) {
               MaxPlot = histo[s][ch]->GetMaximum();
             } 
          } 
          

          histo[0][ch] -> SetLineColor(2); histo[0][ch] -> SetLineWidth(2);
          histo[0][ch] -> SetMaximum(MaxPlot*1.1);
          histo[0][ch] -> SetTitle (varName);
          histo[0][ch] -> SetXTitle (xtitle);
          for (int sma = 1; sma < nsample; sma++){histo[sma][ch] -> SetLineColor(2+sma); histo[sma][ch]->SetLineWidth(2);}

          histo[0][ch] -> Draw("histo");
          for (int sma = 1; sma < nsample; sma++){histo [sma][ch]-> DrawCopy("histsame");}
          _leg ->Draw();

       
          CC -> Print("iplots/shape/" + massP +"/" +  varName + "_" + schannel[ch] + ".png");
         }

       gSystem->Exec("cp  ./index.php ./iplots/; done");
       gSystem->Exec("for dir in $(find ./iplots/ -type d); do cp -n ./index.php $dir/; done");
       //gSystem->Exec("rm -f index.php");
     
}



//-------------------------------------------------------------------------
// Significance Plot and Calculations
//-------------------------------------------------------------------------


void significance (float Smass, float Xmass, int ncut, TString varName){

      TString massP;
      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }


     TFile* Drawhisto   [nsample]            [ncut];
     TH1F*  histo2      [nsample]  [nchannel][ncut];
   

     const int  nbckg = nsample-1;
 
     TH1F*  signf_sqrt  [nbckg][nchannel];             
     TH1F*  signf_s     [nbckg][nchannel];        
     TH1F*  signft_sqrt [nchannel];             
     TH1F*  signft_s    [nchannel];        
     
    
     for (int ch = 0; ch < nchannel; ch++){

       signft_sqrt[ch] = new TH1F ("S_sqrt(S+B)_Significance"       + schannel[ch], "", ncut , 0., ncut); // Total S/srt(S +B)  event 
       signft_s   [ch] = new TH1F ("S_B_Significance"               + schannel[ch], "", ncut , 0., ncut); // Total S/B  event counter 
       for (int s = 0; s < nbckg; s++){
     	 signf_sqrt [s][ch]  = new TH1F ("S_sqrt(S+B)_Significance" + schannel[ch], "", ncut , 0., ncut); // Each Bckg S/srt(S +B)  event counter
     	 signf_s    [s][ch]  = new TH1F ("S_B_Significance"         + schannel[ch], "", ncut , 0., ncut); // Each Bckg S/B  event counter
       }
     }

    const int label_ncut = ncut+1;  
    for (int i = 1; i <= ncut; i++){
      for (int ch = 0; ch < nchannel; ch++){
    //     
       signft_sqrt [ch] -> GetXaxis()-> SetBinLabel (i, scut[i-1]);
       signft_s    [ch]-> GetXaxis()-> SetBinLabel (i, scut[i-1]);
       for (int s = 0; s < nbckg; s++){
         signf_sqrt [s][ch]-> GetXaxis()-> SetBinLabel (i, scut[i-1]);
         signf_s    [s][ch]-> GetXaxis()-> SetBinLabel (i, scut[i-1]);
       }
      }
     }
      for (int cut = 0; cut < ncut; cut++)    {
        for (int s = 0; s< nsample; s++)    { 
       
	   if (_sample[s] == "T2tt"){
            if (Smass==1. && Xmass==1.) { _sample[s] += Reg;
             } else {
              _sample[s] += massP;
             }
           }

           Drawhisto[s][cut] = new TFile( "irootfile/" + _sample[s] + "/" + scut[cut] + "/" + "h_" + varName  +".root","read");  
       
	 }        
      }

      for (int ch = 0; ch < nchannel; ch++) {

     //   signft_sqrt[ch] -> Reset();
     //   signft_s   [ch] -> Reset();
     //   for (int s = 0; s < nbckg; s++){
     //     signf_sqrt [s][ch] -> Reset();
     //     signf_s    [s][ch] -> Reset();
     //   }     
 
      	for (int cut = 0; cut < ncut; cut++)    {
        
          //if (ch == 3 && svarCut [cut] == "nbjet20csvv2m" ) continue;          

          for (int s = 0; s< nsample; s++)    {
       
	     histo2[s][ch][cut] = (TH1F*)Drawhisto[s][cut]-> Get ("h_" + varName + "_" + schannel[ch] + "_" + scut[cut]);
          
	  }             
         
        cout << "                                                                            " << endl;
        cout << "--------------------------------------------------------------------------- " << endl;
        cout << " Cut:     " << svarCut [cut] << "  >=  " << valCut[cut] << endl;
        cout << "                                                                            " << endl;
	cout << " channel  " << schannel[ch] << endl;  
        cout << "--------------------------------------------------------------------------- " << endl;
        cout << "                                                                            " << endl;  
        cout << " Yields  "  << "    #  events   " << "    S/sqrt(S+B) " << "       S/B      " <<endl; 
        cout << "                                                                            " << endl;
        cout << "--------------------------------------------------------------------------- " << endl;
    
        float SignalYields = histo2[0][ch][cut] -> Integral(0, histo2[0][ch][cut] -> GetNbinsX()+1);

        cout << _sample[0] << "       " <<  SignalYields << endl;
        cout << "                                                      " << endl;  
        
        
        float BackgroundYield = 0.;     
        for (int b = 1; b < nsample; b++) {

          float ThisBckYield = histo2[b][ch][cut] -> Integral(0, histo2[b][ch][cut] -> GetNbinsX()+1);
          BackgroundYield += ThisBckYield;
          float Sign_sqrt = SignalYields / sqrt(SignalYields + ThisBckYield);
          float Sign_s   = SignalYields / ThisBckYield;
        
          signf_sqrt [b-1][ch] -> SetBinContent (cut+1, Sign_sqrt);       
          signf_s    [b-1][ch] -> SetBinContent (cut+1, Sign_s);       
  
          cout <<  _sample[b] << "            " << ThisBckYield << "      " << Sign_sqrt << "          "        << Sign_s << endl; 
          cout << "                                                                                             " << endl;
        }

        float Signt_sqrt = SignalYields / sqrt(SignalYields + BackgroundYield);
        float Signt_s    = SignalYields / sqrt(SignalYields + BackgroundYield);
 
        signft_sqrt [ch] -> SetBinContent (cut+1, Signt_sqrt);
        signft_s    [ch] -> SetBinContent (cut+1, Signt_s);
 
       cout << "ToTal MC              " << BackgroundYield << "          " <<Signt_sqrt << "         "  << Signt_s <<endl;
     }


     

     // Set the maximun of the distribution
     // ----------------------------------


	float MaxPlot_sqrt = 0.001;
	float MaxPlot_s = 0.001;
          for (int s = 1; s< nsample; s++){
           //  _leg -> AddEntry(histo[s][ch], _sample[s]  , "f"); //drawShape(histo[s][ch]);
             if (signf_sqrt[s-1][ch]->GetMaximum()>MaxPlot_sqrt) {
               MaxPlot_sqrt = signf_sqrt[s-1][ch]->GetMaximum();
             }
             if (signf_s[s-1][ch]->GetMaximum()>MaxPlot_s) {
               MaxPlot_s = signf_s[s-1][ch]->GetMaximum();
             }

          }



  
     gSystem->mkdir("iplots/significance" + massP +"/", kTRUE);

     // CSignificance sqrt Total background -> Draw()
     // ----------------------------------------------
    
     TCanvas *CSignif2_total = new TCanvas("CSignif2_total" + schannel[ch], "", 1200, 1000);
     gStyle -> SetOptStat("");
     signft_sqrt[ch] -> SetFillColor(1); signft_sqrt[ch] -> SetTitle  ("S/sqrt(S + B) total    " + schannel[ch]); signft_sqrt[ch] -> SetYTitle ("S/sqrt(S + B)"); 
     signft_sqrt[ch] -> Draw();     
     CSignif2_total  -> Print ("iplots/significance" + massP +"/" + "CSignif2_total_" + schannel[ch]+ ".png");
     
     // CSignificance sqrt Each background -> Draw()
     // -----------------------------------------------    

     TCanvas *CSignif2_each = new TCanvas("CSignif2_each" + schannel[ch], "", 1200, 1000);
     gStyle -> SetOptStat("");
     TLegend*  _leg1 = myLegend ("");
     for (int b = 1; b<nsample; b++){signf_sqrt[b-1][ch]-> SetLineColor(1+b); signf_sqrt [b-1][ch] ->SetLineWidth(5); _leg1 -> AddEntry(signf_sqrt[b-1][ch],_sample[b], "f"); }
     signf_sqrt[0][ch] -> SetTitle  ("S/sqrt(S + B) each bck_" + schannel[ch]);
     signf_sqrt[0][ch] -> SetMaximum (MaxPlot_sqrt*1.1);
     signf_sqrt[0][ch] -> Draw("histo");
     for (int b = 1; b<nsample; b++){signf_sqrt[b-1][ch]->Draw("histsame");}
     _leg1->Draw();
    
     CSignif2_each     -> Print ("iplots/significance" + massP +"/" +  "CSignif2_eachBck_" + schannel[ch]+ ".png");   


     // CSignificance S/B Total background -> Draw()
     // ----------------------------------------------

     TCanvas *CSignif1_total = new TCanvas("CSignif1_total" + schannel[ch], "", 1200, 1000);
     gStyle -> SetOptStat("");
     signft_s[ch] -> SetFillColor(1); signft_s[ch] -> SetTitle  ("S/B total_"+ schannel[ch]); signft_s[ch] -> SetYTitle ("S/B");
     signft_s[ch] -> Draw();
     CSignif1_total    -> Print ("iplots/significance" + massP +"/" + "CSignif1_total_" + schannel[ch]+ ".png");
 
     // CSignificance S/B Each background -> Draw()
     // ---------------------------------------------
     
     TCanvas *CSignif1_each = new TCanvas("CSignif1_each" + schannel[ch], "", 1200, 1000);
     gStyle -> SetOptStat("");
     TLegend*  _leg2 = myLegend ("");
     for (int b = 1; b<nsample; b++){signf_s[b-1][ch]-> SetLineColor(1+b); signf_s [b-1][ch] ->SetLineWidth(5); _leg2 -> AddEntry(signf_s[b-1][ch],_sample[b], "f"); }
     signf_s[0][ch] -> SetTitle  ("S/B each bck     " + schannel[ch]);
     signf_s[0][ch] -> SetMaximum (MaxPlot_s*1.1);
     signf_s[0][ch] -> Draw("histo");  
     for (int b = 1; b<nsample; b++){signf_s[b-1][ch]->Draw("histsame");} 
     _leg2 -> Draw();

     CSignif1_each-> Print ("iplots/significance" + massP +"/" +  "CSignif1_eachBck_" + schannel[ch]+ ".png");
   }
  
  gSystem->Exec("cp  ./index.php ./iplots/; done");
  gSystem->Exec("for dir in $(find ./iplots/ -type d); do cp -n ./index.php $dir/; done");

}



//--------------------------------------------------------
// Plot everything
//--------------------------------------------------------



void All_shapePlots( float Smass, float Xmass){

     Plot_shape(Smass, Xmass, "jet2pt",       "");
     Plot_shape(Smass, Xmass, "jet1pt",       "");
     Plot_shape(Smass, Xmass, "lep1pt",       "");
     Plot_shape(Smass, Xmass, "lep2pt",       "");
     Plot_shape(Smass, Xmass, "lep1phi",      "");
     Plot_shape(Smass, Xmass, "lep2phi",      "");
     Plot_shape(Smass, Xmass, "mlb1",         "");
     Plot_shape(Smass, Xmass, "mlb2",         "");
     Plot_shape(Smass, Xmass, "min_mlb",      "");
     Plot_shape(Smass, Xmass, "max_mlb",      "");
     Plot_shape(Smass, Xmass, "metPfType1",   "");
     Plot_shape(Smass, Xmass, "dphimetjet",   "");
     Plot_shape(Smass, Xmass, "met_sqrtHt",   "");
     Plot_shape(Smass, Xmass, "met_meff",     "");
     Plot_shape(Smass, Xmass, "Ht",           "");
     Plot_shape(Smass, Xmass, "MT2ll",        "");
     Plot_shape(Smass, Xmass, "MT2bb",        "");
     Plot_shape(Smass, Xmass, "MT2lblb",      "");
     Plot_shape(Smass, Xmass, "dyll",         "");
     Plot_shape(Smass, Xmass, "ptbll",        "");
     Plot_shape(Smass, Xmass, "dphimetptbll", "");
     Plot_shape(Smass, Xmass, "m2l",          "");
     Plot_shape(Smass, Xmass, "mllbb",        "");
     Plot_shape(Smass, Xmass, "meff",         "");
     Plot_shape(Smass, Xmass, "dphillmet",    "");
     Plot_shape(Smass, Xmass, "dphill",       "");
     Plot_shape(Smass, Xmass, "dphilmet1",    "");
     Plot_shape(Smass, Xmass, "dphijet1met",  "");
     Plot_shape(Smass, Xmass, "htjets",       "");
     Plot_shape(Smass, Xmass, "metPfType1Phi","");
     Plot_shape(Smass, Xmass, "MR",           "");
     Plot_shape(Smass, Xmass, "R2",           "");
     Plot_shape(Smass, Xmass, "Rpt",          "");
     Plot_shape(Smass, Xmass, "invGamma",     "");
     Plot_shape(Smass, Xmass, "Mdr",          "");
     Plot_shape(Smass, Xmass, "DeltaPhiRll",  "");


}

void All_cutPlots(float Smass, float Xmass){

     Plot_cut(Smass, Xmass, "jet2pt",       "");
     Plot_cut(Smass, Xmass, "jet1pt",       "");
     Plot_cut(Smass, Xmass, "lep1pt",       "");
     Plot_cut(Smass, Xmass, "lep2pt",       "");
     Plot_cut(Smass, Xmass, "lep1phi",      "");
     Plot_cut(Smass, Xmass, "lep2phi",      "");
     Plot_cut(Smass, Xmass, "mlb1",         "");
     Plot_cut(Smass, Xmass, "mlb2",         "");
     Plot_cut(Smass, Xmass, "min_mlb",      "");
     Plot_cut(Smass, Xmass, "max_mlb",      "");

     Plot_cut(Smass, Xmass, "metPfType1",   "");
     Plot_cut(Smass, Xmass, "dphimetjet",   "");
     Plot_cut(Smass, Xmass, "met_sqrtHt",   "");
     Plot_cut(Smass, Xmass, "met_meff",     "");
     Plot_cut(Smass, Xmass, "Ht",           "");
     Plot_cut(Smass, Xmass, "MT2ll",        "");
     Plot_cut(Smass, Xmass, "MT2bb",        "");
     Plot_cut(Smass, Xmass, "MT2lblb",      "");
     Plot_cut(Smass, Xmass, "dyll",         "");
     Plot_cut(Smass, Xmass, "ptbll",        "");
     Plot_cut(Smass, Xmass, "dphimetptbll", "");
     Plot_cut(Smass, Xmass, "m2l",          "");
     Plot_cut(Smass, Xmass, "mllbb",        "");
     Plot_cut(Smass, Xmass, "mllbb",        "");
     Plot_cut(Smass, Xmass, "meff",         "");
     Plot_cut(Smass, Xmass, "dphillmet",    "");
     Plot_cut(Smass, Xmass, "dphill",       "");
     Plot_cut(Smass, Xmass, "dphilmet1",    "");
     Plot_cut(Smass, Xmass, "dphijet1met",  "");
     Plot_cut(Smass, Xmass, "htjets",       "");
     Plot_cut(Smass, Xmass, "metPfType1Phi","");
     Plot_cut(Smass, Xmass, "MR",           "");
     Plot_cut(Smass, Xmass, "R2",           "");
     Plot_cut(Smass, Xmass, "Rpt",          "");
     Plot_cut(Smass, Xmass, "invGamma",     "");
     Plot_cut(Smass, Xmass, "Mdr",          "");
//     Plot_cut(Smass, Xmass, "DeltaPhiRll",  "");
//     Plot_cut(Smass, Xmass, "mlb1",         "");
//     Plot_cut(Smass, Xmass, "mlb2",         "");
//     Plot_cut(Smass, Xmass, "min_mlb",      "");
//     Plot_cut(Smass, Xmass, "max_mlb",      "");
//     Plot_cut(Smass, Xmass, "Ht_visible",   "");
}

void All_datacard (float Smass, float Xmass, TString CutName){

     Ingredients( Smass,  Xmass, "MT2ll",        CutName);
     Ingredients( Smass,  Xmass, "MT2bb",        CutName);
     Ingredients( Smass,  Xmass, "MT2lblb",      CutName);
     Ingredients( Smass,  Xmass, "mllbb",        CutName);
     Ingredients( Smass,  Xmass, "Mdr"  ,        CutName);
     Ingredients( Smass,  Xmass, "DeltaPhiRll" , CutName);
     Ingredients( Smass,  Xmass, "invGamma",     CutName);
     Ingredients( Smass,  Xmass, "Rpt",          CutName);
     Ingredients( Smass,  Xmass, "R2",           CutName);
     Ingredients( Smass,  Xmass, "MR",           CutName);
     Ingredients( Smass,  Xmass, "dphillmet",     CutName);
     Ingredients( Smass,  Xmass, "dphijet1met",  CutName);
     Ingredients( Smass,  Xmass, "dphimetptbll", CutName);
     Ingredients( Smass,  Xmass, "dyll",         CutName);
     Ingredients( Smass,  Xmass, "metPfType1",   CutName);
     Ingredients( Smass,  Xmass, "met_meff",     CutName);
     Ingredients( Smass,  Xmass, "jet2pt",       CutName);
     Ingredients( Smass,  Xmass, "jet1pt",       CutName);
     Ingredients( Smass,  Xmass, "htjets",       CutName);
     Ingredients( Smass,  Xmass, "ptbll",        CutName);
     Ingredients( Smass,  Xmass, "mlb1",         CutName);
     Ingredients( Smass,  Xmass, "mlb2",         CutName);
     Ingredients( Smass,  Xmass, "min_mlb",      CutName);
     Ingredients( Smass,  Xmass, "max_mlb",      CutName);
     Ingredients( Smass,  Xmass, "Ht_visible",   CutName);

}
       
//int main() // CORRECTO: la forma adecuada de hacerlo
//{





//     return 0;
//}     
