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
#include <typeinfo>
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

     TString FileAddress = "../../minitreesG/nominal/Stop/";
     
     TString  Reg = "R4";

     const int nsample = 7;
     TString SampleName [nsample] = { "16_T2tt_mStop-150to1200.root", "04_TTTo2L2Nu.root", "05_ST.root",  "07_ZJets.root", "06_WW.root", "09_TTV.root", "03_VZ.root"}; 
     TString  _sample   [nsample] = {"T2tt", "TTbar", "ST", "DYJets", "WW", "TTV", "VZ"};

     //-----------------------------------------------------------------------------------------
     // 28 Nov 2016
     //-----------------------------------------------------------------------------------------
     // T2tt  = T2tt_mStop-150to1200.root
     // TTbar = TTTo2L2Nu_ext1__part*.root 
     // ST    = ST_tW_antitop.root ST_tW_top.root
     // ZJets = DYJetsToLL_M-10to50.root DYJetsToLL_M-50_000*.root;
     // WW    = WWTo2L2Nu.root GluGluWWTo2L2Nu_MCFM.root
     // TTV   = TTWJetsToLNu.root TTZjets.root
     // VZ    = ZZTo2L2Q*.root WZTo2L2Q__part*.root
     // ----------------------------------------------------------------------------------------
     // 21 Nov 2016
     // ----------------------------------------------------------------------------------------     
     // const int nsample = 6;
     // TString SampleName [nsample] = {"T2tt_mStop-150to1200.root", "TTTo2L2Nu_ext1.root", "ST_tW.root",  "DYJetsToLL.root", "WWTo2L2Nu.root", "TTZjets.root"}; 
     // TString  _sample   [nsample] = {"T2tt", "TTbar", "ST", "DYJets", "WW", "TTZ"};    
     // -----------------------------------------------------------------------------------------
     // T2tt  = T2tt_mStop-150to1200.root
     // TTbar = TTTo2L2Nu_ext1__part*.root 
     // ST    = ST_tW_antitop.root ST_tW_top.root
     // ZJets = DYJetsToLL_M-10to50.root DYJetsToLL_M-50_000*.root DYJetsToTT_MuEle_M-50.root
     // WW    = WWTo2L2Nu.root 
     // TTZ   = TTZjets.root
     // -----------------------------------------------------------------------------------------                             
     
     // -----------------------------------------------------------------------------------------
     // 29 Nov 2016
     // -----------------------------------------------------------------------------------------

       const int ncut = 3;
       TString scut [ncut] = { "MET140", "MT2ll40", "tagL"}; // tried: "MET50", "MET80", "MET200"  
       float valCut [ncut] = {140,40,20}; //In this case the valCut is an upper limit for the variable varCut ( the cut is: if LeadingPtCSVv2L < 20 -> pass = true -> the code will fill the histograms)  

     // -----------------------------------------------------------------------------------------
     // 28 Nov 2016
     // -----------------------------------------------------------------------------------------
     
     // const int ncut = 4;
     // TString scut [ncut] = { "-999jet2ptCheckMet", "50metCheckMet", "20jet2ptCheckMet", "20nbjetptcsvv2mCheckMet"};
     // float valCut [ncut] = {-999, 50, 20, 20}; //In this case the valCut in a low limit for the variable varCut (the cut is: if varCut[i] > valCut[i] -> pass = true -> the code will fill the histograms)

     //  ------------------------------------------------------------------------------------------
     //  21 Nov 2016
     //  -----------------------------------------------------------------------------------------
     
     // const int ncut = 3;
     // TString scut [ncut] = { "nbjetptcsvv2mXCheck0Met" , "20jet2ptXCheck0Met", "50metXCheck0Met"};
     // float valCut [ncut] = { -1, 20, 0}; //In this case the valCut in a low limit for the variable varCut (the cut is: if varCut[i] > valCut[i] -> pass = true -> the code will fill the histograms)

     // ------------------------------------------------------------------------------------------
     

     enum { ee, mm, em, nchannel}; const TString schannel[nchannel] = {"ee", "mm","em"};
 
     enum discriminant { _jet2pt,  _metPfType1, _dphijetmet, _met_sqrtHt, _met_meff, _Ht, _MT2ll, _MT2bb, _MT2lblb, _dyll, _ptbll, _dphimetptbll, _m2l, _mllbb,  _meff, _dphillmet, _dphill, _dphilmet1, _dphijet1met, _htjets, _metPfType1Phi, _MR, _R2, _Rpt, _invGamma, _Mdr, _DeltaPhiRll, _Ht_visible, nD};
     const TString sD [nD] = {"_jet2pt", "_metPfType1", "_dphimetjet", "_met_sqrtHt", "_met_meff", "_Ht", "_MT2ll", "_MT2bb", "_MT2lblb", "_dyll", "_ptbll", "_dphimetptbll", "_m2l", "_mllbb", "_meff", "_dphillmet", "_dphill", "_dphilmet1", "_dphijet1met", "_htjets", "_metPfType1Phi", "_MR", "_R2", "_Rpt", "_invGamma", "_Mdr", "_DeltaPhiRll", "_Ht_visible"};

     float binning [nD] = {60, 40, 20, 50, 5, 100, 30, 80, 80, 80, 80, 20, 80, 80, 80, 20, 20, 20, 20, 80, 20, 100, 80, 80, 80, 18, 80, 100};
     float xmin         = 0.;
     float xmax    [nD] = {450., 800., 4., 300., 1., 1000., 600., 800., 800., 5., 800., 4., 800., 800., 1500., 4., 4., 4., 4., 1100., 4., 1000., 1.4, 1.1, 1.5, 450., 3.5, 1000};
     
     TString histoName [nD];   
     TH2F*   histo2d   [nD][nchannel][ncut]; // nD NOT excludes the h_i_i, i = _jet2pt, _metPfType1,...., _DeltaPhiRll -> You need nD-1;    

     bool ReopenRoot = false;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Histo Definitions
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



void Histo2dName( char Discriminant ){

    
     TString variable = sD[Discriminant];

     for (int i = 0; i < nD; i++){

       if ( i == Discriminant) continue; 
       histoName[i] = "h" + variable +  sD[i];          
       cout << histoName[i] << endl;
     }
    
}     


void define2dHistograms (char Discriminant, int ichannel, int icut, TString suffix){
  
     Histo2dName( Discriminant);

     for (int i = 0; i < nD; i++){
        //cout << histoName [i]<< "                         " << sD[i] << "                   "  <<  binning[i] << "                      " << xmax[i] <<endl;
        histo2d [i][ichannel][icut] = new  TH2F( histoName [i]   + suffix, "", binning[Discriminant],  0.,  xmax[Discriminant], binning[i], 0., xmax[i]);
     }
}

void reset2dHistograms(char Discriminant, int ichannel, int icut, TString suffix){

     for (int i = 0; i < nD; i++){
       histo2d [i][ichannel][icut] -> Reset();
     }
}


//------------------------------------------
// Mass Point String
// ----------------------------------------


TString massPoint(float Smass, float Xmass){

   TString sSmass = ""; sSmass += int(Smass); TString sXmass = ""; sXmass += int(Xmass);
   TString masspoint = "_" + sSmass + "_" + sXmass;
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
    leg1->SetHeader(title); //leg1_fit->SetMargin(0.2); 
    leg1->AddEntry((TObject*)0, " ", "");
    
    return leg1; 
     }


//-------------------------------------------------------
// Save Histograms as .root File
// ------------------------------------------------------


void save2dHistoRoot (TString irootfileName, TH2F* histoName){
     TFile* OutFileName = new TFile( irootfileName , "update");
     histoName   -> Write();
     OutFileName -> Close();
}

void saveHistoRoot (TString irootfileName, TH1F* histoName){
     TString Operation = "update";
     if (ReopenRoot) Operation = "recreate";
     TFile* OutFileName = new TFile( irootfileName , Operation);
     histoName   -> Write();
     OutFileName -> Close();
     ReopenRoot = false;
}

//-------------------------------------------------------
// Normalize the Histogram to Unity 
// ------------------------------------------------------

void drawShape (TH1F *histoName){
    float integral;
    integral = histoName -> Integral(0, histoName -> GetNbinsX()+1);
    cout  << integral << endl; 
    histoName -> Scale(1/integral);
}


//----------------------------------------------------------------------------------
// Cut Study
// ---------------------------------------------------------------------------------



                  
void cut_study(char Discriminant, float Smass, float Xmass) {

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
       define2dHistograms( Discriminant, ichannel, icut, suffix);
     }
  }


  //      //
  // Loop //
  //      //

 // for (int s = 3; s < 4 ; s++){
  for (int s = 0; s < nsample ; s++){

     if (_sample[s] == "T2tt") _sample[s] += massP;
     gSystem->mkdir("irootfile/" + _sample[s] + "/2D/", kTRUE);

 
     TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]); 
     TTree *MiniTree      = GetMiniTree(MiniTreeFile);
     Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

     for (int icut = 0; icut < ncut; icut++){
        for (int ichannel = 0; ichannel < nchannel; ichannel++) {

          TString suffix = "_" + schannel[ichannel] + "_" + scut[icut];
          reset2dHistograms(Discriminant, ichannel, icut, suffix);
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


       //                           //
       // Variables                 //
       //                           //


        float D [nD] = { jet2pt,  metPfType1, dphimetjet, V4, V3, ht, mt2ll, mt2bb, mt2lblb, dyll, ptbll, dphimetptbll, m2l, mllbb,  meff, dphillmet, dphill, dphilmet1, dphijet1met, htjets, metPfType1Phi, MR, R2, Rpt, invGamma, Mdr, DeltaPhiRll, V5};


       //               //
       // varCut [ncut] //
       //               //


       //float varCut [ncut] = { LeadingPtCSVv2M , jet2pt, metPfType1};
         float varCut [ncut] = { LeadingPtCSVv2L};
       //float varCut [ncut] = { LeadingPtCSVv2L, LeadingPtCSVv2M, LeadingPtCSVv2T};


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

           if (channel!=5 || scut[j]!="nbjetptcsvv2m")//////
             pass &= varCut[j] >= valCut [j];
           
           if (pass) { for (int d = 0; d < nD; d++){
                       if (d == Discriminant) continue;
                       //if (mt2ll > 200) cout << sD [Discriminant] << "     " <<  D[Discriminant] << "\n" << mt2ll<< endl; //sD [d] << "       " << D[d] << endl; 
                       //if (eventW<0) cout <<eventW<<endl;
                       histo2d[d][i][j] -> Fill( D[Discriminant], D[d], eventW); 
                     }

           }
        }

      }         
         
          for (int ch = ee; ch < nchannel; ch ++){
            for (int cut = 0; cut < ncut; cut++){
  
              gSystem->mkdir("irootfile/" + _sample[s] + "/2D/" + scut[cut] + "/", kTRUE);

              for (int nhist = 0; nhist < nD; nhist++){
                 save2dHistoRoot ("irootfile/" + _sample[s] + "/2D/" + scut[cut] + "/" + histoName[nhist]         + ".root", histo2d                 [nhist][ch][cut]);       
              }
            }
         }
       }
}





//------------------------------------------------------------------------------------------------------------------------
// Shape Study
//-----------------------------------------------------------------------------------------------------------------------



void shape_study(char Discriminant, float Smass, float Xmass){

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
        
        const int icut = 0;       
        TString suffix = "_" + schannel[ichannel] ;
        define2dHistograms( Discriminant, ichannel, icut, suffix);
     
       }

     
      //      //
      // Loop // 
      //      //

      for (int s = 0; s < nsample ; s++){

      if (_sample[s] == "T2tt") _sample[s] += massP;
      gSystem->mkdir("irootfile/" + _sample[s] + "/shape/2D/", kTRUE);

      TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]);
      TTree *MiniTree      = GetMiniTree(MiniTreeFile);
      Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

      for (int ichannel = 0; ichannel < nchannel; ichannel++) {
      
        const int icut = 0;
        TString suffix = "_" + schannel[ichannel] + "_";
        reset2dHistograms(Discriminant, ichannel, icut, suffix);
      }

      //for (Int_t en = 0; en<500; en++){
      for (Int_t en = 0; en<nentries; en++){
        MiniTree -> GetEntry(en);
        cout << eventW << endl;  
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



        //                            //
        // Variables                  //
        //                            //

    
        float D [nD] = { jet2pt,  metPfType1, dphimetjet, V4, V3, ht, mt2ll, mt2bb, mt2lblb, dyll, ptbll, dphimetptbll, m2l, mllbb,  meff, dphillmet, dphill, dphilmet1, dphijet1met, htjets, metPfType1Phi, MR, R2, Rpt, invGamma, Mdr, DeltaPhiRll,V5};


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

         for (int d = 0; d < nD; d++){
           if (d == Discriminant) continue;
              histo2d[d][i][j] -> Fill(D[Discriminant], D[d], eventW);
          }


      
       }
 
       //      //
       // Save //
       //      //
        
       for (int ch = ee; ch < nchannel; ch ++){
  
         const int icut = 0;
       
         for (int nhist = 0; nhist < nD; nhist++){
                 save2dHistoRoot ("irootfile/" + _sample[s] + "/shape/2D/" + histoName[nhist]         + ".root", histo2d                 [nhist][ch][icut]);
              }

       }

}
}



//-----------------------------------------------------------------------------------
// Shape Projection 
//-----------------------------------------------------------------------------------


void shape_projection(char Discriminant1, char Discriminant2, int icut, float Smass, float Xmass, int nbinV1, float xmaxV1, int ncutsV2, float xmaxV2 = -1.){

     ReopenRoot = true;

     if (xmaxV2==-1) xmaxV2 = xmax[Discriminant2];

     TString massP;

      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }

 
     //                   //
     // Define Histograms //
     //                   //

     
     TH1F*   shape_pro [nchannel][nsample];
                            
     for (int ichannel = 0; ichannel < nchannel; ichannel++) {
       for (int s = 0; s < nsample ; s++){
      //  TString suffix = "_shape_pro_" + sD[Discriminant2] + "_" + schannel[ichannel] + "_" + scut[icut];
       if (_sample[s].Contains("T2tt")) _sample[s] = "T2tt" +  massP; 
       shape_pro [ichannel][s] = new TH1F(_sample[s] + "_" + schannel[ichannel], "", nbinV1*ncutsV2, 0., xmaxV1*ncutsV2);
       cout << _sample[s] <<endl;
       shape_pro [ichannel][s] -> Sumw2();
       }
     }

     //                    //
     // Loop               //
     //                    //

     for (int s = 0; s < nsample ; s++){

     gSystem->mkdir("datacard/T2tt" + massP + "/", kTRUE); 

     TFile *MiniTreeFile  = TFile::Open(FileAddress + SampleName [s]);
     TTree *MiniTree      = GetMiniTree(MiniTreeFile);
     Int_t  nentries      = (Int_t) MiniTree -> GetEntries();

     for (int ichannel = 0; ichannel < nchannel; ichannel++) {
     //  TString suffix = "_shape_pro_" + schannel[ichannel] + "_" + scut[icut];
       shape_pro [ichannel][s] -> Reset();
     }


     for (Int_t en = 0; en<nentries; en++){
       MiniTree -> GetEntry(en);

       //eventW = 1; 
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


       //                           //
       // Variables                 //
       //                           //

       float D [nD] = { jet2pt,  metPfType1, dphimetjet, V4, V3, ht, mt2ll, mt2bb, mt2lblb, dyll, ptbll, dphimetptbll, m2l, mllbb,  meff, dphillmet, dphill, dphilmet1, dphijet1met, htjets, metPfType1Phi, MR, R2, Rpt, invGamma, Mdr, DeltaPhiRll, V5};

       if (D[Discriminant1]>=xmaxV1) D[Discriminant1] = xmaxV1 - 0.01;
       if (D[Discriminant2]>=xmaxV2) D[Discriminant2] = xmaxV2 - 0.01;

      //               //
      // varCut [ncut] //
      //               //

//TString scut [ncut] = { "-999jet2ptCheckTag", "50metCheckTag", "20jet2ptCheckTag", "20nbjetptcsvv2mCheckTag"};

    //  float varCut [ncut] = { jet2pt, metPfType1, jet2pt, LeadingPtCSVv2M}; // 28 Nov 2016
        float varCut [ncut] = { metPfType1, mt2ll, LeadingPtCSVv2L}; // 30 Nov 2016 Looking for two CR one with b veto loose one with 1 b tag loose 
        //float varCut [ncut] = { metPfType1, LeadingPtCSVv2L}; // 29 Nov 2016 Looking for two CR one with b veto loose one with 1 b tag loose 
        // float varCut[ncut] = { LeadingPtCSVv2L, LeadingPtCSVv2M, LeadingPtCSVv2T};// 29 Nov 2016 Looking for two CR one with b veto loose one with 1 b tag loose and check the better WP for this region 
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
      float cortes = xmaxV2/ncutsV2;


      // Concatenated cuts     

      for (int j = 0; j <= icut; j++){

          // if (channel!=5 || scut[j]!="nbjetptcsvv2m") 
              pass &= varCut[j] >= valCut [j];
      }
  
      //  + tagLVETO
     // pass &= LeadingPtCSVv2L < 20;
      // 

       if (pass) {
                 for (int N = 1; N <= ncutsV2; N++){ 
                   if ( (N-1)*cortes <= D[Discriminant2] && D[Discriminant2]< cortes*N) {
                      if (D[Discriminant1]+xmaxV1*(N-1)>=xmaxV1*ncutsV2) {       
                         cout << D[Discriminant1]<< " "<< D[Discriminant2]<<" " <<mt2ll << " " << Mdr << endl;
                      }
  				shape_pro [i][s] -> Fill (D[Discriminant1]+xmaxV1*(N-1), eventW);
                   }
                 }    
        }

    }
   
   for (int ch = 0; ch < nchannel; ch ++){
    saveHistoRoot ("datacard/T2tt" + massP + "/" + "T2tt_" + massP + "_shapes_proj_G_Met140MT2ll40tagL_" + sD[Discriminant1] + "_" + sD[Discriminant2] + "_TH1F.root", shape_pro [ch][s] );
    if ( s== 1) { shape_pro[ch][1]->SetName("data_obs_" + schannel[ch]); 
    saveHistoRoot ("datacard/T2tt" + massP + "/" + "T2tt_" + massP + "_shapes_proj_G_Met140MT2ll40tagL_" + sD[Discriminant1] + "_" + sD[Discriminant2] + "_TH1F.root", shape_pro [ch][1] ); 
    }
   }
  }
   std::ofstream inFile("datacard/T2tt" + massP + "/" + "T2tt_" + massP + "_shapes_proj_G_Met140MT2ll40tagL_" + sD[Discriminant1] + "_" + sD[Discriminant2] + ".txt",std::ios::out);

   inFile << "\n";
   inFile << "   imax  "<<"     " << nchannel   << endl;
   inFile << "   jmax "<<"     " << nsample -1 << endl;
   inFile << "   kmax "<<"     " <<  "*"       << endl;
   inFile << "-----------------------------------------------------------------------------------" << endl;
   inFile << "  shapes * * T2tt_" +  massP + "_shapes_proj_G_Met140MT2ll40tagL_" + sD[Discriminant1] + "_" + sD[Discriminant2]  + "_TH1F.root  $PROCESS_$CHANNEL" << endl;
   inFile << "-----------------------------------------------------------------------------------" << endl;
   inFile << "   bin  "<<"      "<< "     ee          "<<"     mm             "<<"      em             "<<endl;
   inFile << "observation           ";
   for (int ch = 0; ch < nchannel; ch++){ float NBck = shape_pro [ch][1] -> Integral(0, shape_pro[ch][1] -> GetNbinsX()+1);
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
   for (int s = 0; s < nsample; s++){ for (int ch = 0; ch < nchannel; ch++){float integral = shape_pro[ch][s] -> Integral(0, shape_pro[ch][s] -> GetNbinsX()+1);
   inFile << integral<<"      ";}}
   inFile.close();

 
 
}

//------------------------------------------------------------------------------------
// Plot the TH2F histograms
//------------------------------------------------------------------------------------




void d2Plot ( float Smass, float Xmass, int icut, TString sampleName, TString d2VarName, TString xtitle){

      TString massP;

      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }

      if (sampleName == "T2tt"){
         if (Smass==1. && Xmass==1.) { sampleName += Reg;
         } else {
          sampleName += massP;
         }
      }


      TH1F* histo[nchannel];
      TFile* Drawhisto;

      if (icut == 0){
         Drawhisto = new TFile( "irootfile/" + sampleName + "/shape/2D/" + "h_" + d2VarName  +".root","read");
      }else{
         Drawhisto = new TFile( "irootfile/" + sampleName + "/2D/" + scut[icut-1] + "/" + "h_" + d2VarName  +".root","read");         
      }

      gSystem->mkdir("iplots/2Dplots/"+ massP +"/" , kTRUE);
      
      for (int ch = 0; ch < nchannel; ch++){

          TCanvas* CC;
          CC = new TCanvas (d2VarName +  "_" + schannel[ch], "", 1200,1000);
          gStyle->SetOptStat("");
          TLegend*  _leg = myLegend ("");

          if (icut == 0){
             histo[ch] = (TH1F*)Drawhisto-> Get ("h_" + d2VarName + "_" + schannel[ch]);
          }else {
             histo[ch] = (TH1F*)Drawhisto-> Get ("h_" + d2VarName + "_" + schannel[ch] + "_" + scut[icut-1]);
          }

          int g = d2VarName.First("_");
          cout << g << endl;
          TString xtitle = d2VarName;
          xtitle.Remove(g);
          TString ytitle = d2VarName;
          ytitle.ReplaceAll( xtitle + "_", "");

          histo[ch] -> SetXTitle (xtitle);  
          histo[ch] -> SetYTitle (ytitle);  
          histo[ch] -> SetTitle (d2VarName);
          histo[ch] -> Draw("box" );
          _leg      -> Draw();
          CC -> Print("iplots/2Dplots/" + massP +"/" + sampleName + "_" + d2VarName + "_" + schannel[ch] + ".png");

          gSystem->Exec("cp  ./index.php ./iplots/; done");
          gSystem->Exec("for dir in $(find ./iplots/ -type d); do cp -n ./index.php $dir/; done");

      }

 
}


void d1_datacard_Plot(float Smass, float Xmass, TString V1, TString V2, TString xtitle){

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
        Drawhisto[s] = new TFile( "datacard/T2tt" + massP + "/T2tt_" + massP + "_shapes_proj_Met50__" + V1 + "__" + V2 + "_TH1F.root",  "read");
      }

        gSystem->mkdir("iplots/shape_proj/"+ massP +"/" , kTRUE);
        for (int ch = 0; ch < nchannel; ch++){

          TCanvas* CC;
          CC = new TCanvas (V1 + "_" + V2  +  "_" + schannel[ch], "", 1200,1000);
          gStyle->SetOptStat("");
          TLegend*  _leg = myLegend ("");


          for (int s = 0; s< nsample; s++){ histo[s][ch] = (TH1F*)Drawhisto[s]-> Get ( _sample[s] + "_" + schannel[ch]);}

          float MaxPlot = 0.001;
          for (int s = 0; s< nsample; s++){
             _leg -> AddEntry(histo[s][ch], _sample[s]  , "f"); drawShape(histo[s][ch]);
             if (histo[s][ch]->GetMaximum()>MaxPlot) {
               MaxPlot = histo[s][ch]->GetMaximum();
             }
          }


          histo[0][ch] -> SetLineColor(2); histo[0][ch] -> SetLineWidth(2);
          histo[0][ch] -> SetMaximum(MaxPlot*1.1);
          histo[0][ch] -> SetTitle (V1 + "_" + V2);
          histo[0][ch] -> SetXTitle (xtitle);
          for (int sma = 1; sma < nsample; sma++){histo[sma][ch] -> SetLineColor(2+sma); histo[sma][ch]->SetLineWidth(2);}

          histo[0][ch] -> Draw("histo");
          for (int sma = 1; sma < nsample; sma++){histo [sma][ch]-> DrawCopy("histsame");}
          _leg ->Draw();


          CC -> Print("iplots/shape_proj/" + massP +"/" +  V1 + "_" + V2 + "_" + schannel[ch] + ".png");
         }

       gSystem->Exec("cp  ./index.php ./iplots/; done");
       gSystem->Exec("for dir in $(find ./iplots/ -type d); do cp -n ./index.php $dir/; done");

}
//###########################################################################################
// Working in Progress
//#########################################################################################



void Projections ( float Smass, float Xmass, int icut, TString samplename1, TString d2VarName, float lowLimit, float upLimit, bool XonY ){

      TString massP;

      if (Smass == 1. && Xmass == 1.){
        massP = Reg;
      }else{
        massP = massPoint(Smass, Xmass);
      }

      if (samplename1 == "T2tt"){
         if (Smass==1. && Xmass==1.) { samplename1 += Reg;
         } else {
          samplename1 += massP;
         }
      }

      gSystem->mkdir("iplots/2Dplots/Proyections"+ massP +"/" , kTRUE);
      TFile* file;
      TH2F* h_ab;
      TH1F*  pI;
      TH1F*  low_pII;
      TH1F*  up_pII;      
      TString proyName;          
      TString legName;          
  
      if (icut == 0){
        file = new TFile( "irootfile/" + samplename1 + "/shape/2D/" + "h_" + d2VarName  +".root","read");
      }else{
        file = new TFile( "irootfile/" + samplename1 + "/2D/" + scut[icut-1] + "/" + "h_" + d2VarName  +".root","read");
      }

      for (int ch = 0; ch < nchannel; ch++){
 
       
          if (icut == 0){
             h_ab = (TH2F*) file -> Get ("h_" + d2VarName + "_" + schannel[ch]);
          }else {
             h_ab = (TH2F*) file -> Get ("h_" + d2VarName + "_" + schannel[ch] + "_" + scut[icut-1]);
          }

    
        if (XonY){pI = (TH1F*) h_ab -> ProjectionX("x on y");
        }else{pI = (TH1F*) h_ab -> ProjectionY("y on x");}

        float bin_low = pI -> FindBin(lowLimit);
        float bin_up  = pI -> FindBin(upLimit);
        float binLast = pI -> FindLastBinAbove(0,1);
      
        if (XonY){
        low_pII = (TH1F*) h_ab  -> ProjectionX("low x on y", bin_low, bin_up-1);
        up_pII  = (TH1F*) h_ab  -> ProjectionX("up x on y", bin_up,binLast);
        }else{      
        low_pII = (TH1F*) h_ab  -> ProjectionY("low y on x", bin_low, bin_up-1);
        up_pII  = (TH1F*) h_ab  -> ProjectionY("up  y on x", bin_up,   binLast);
        }

        int g = d2VarName.First("_");
        cout << g << endl;
        TString xtitle = d2VarName;
        xtitle.Remove(g);
        TString ytitle = d2VarName;
        ytitle.ReplaceAll( xtitle + "_", "");

 
        TCanvas* cf = new TCanvas ("T2tt_ttbar 0_100","",1200,1000);
        cf -> cd();
        gStyle->SetOptStat("");
        if(XonY){proyName = xtitle; legName = xtitle + "On" + ytitle;}else{ proyName = ytitle; legName = ytitle + "On" + xtitle;}
        TLegend* _leg = myLegend(samplename1 + "\n" + legName );
        
        _leg -> AddEntry(low_pII, "bin0", "f");// "bins: [" , bin_low , "," , bin_up  , ")", "f");  
        _leg -> AddEntry(up_pII,  "bin1", "f");//" bins: [" , bin_up  , "," , binLast , ")", "f"); 
        low_pII -> SetLineColor(1);
        up_pII  -> SetLineColor(6);
        low_pII -> SetXTitle(proyName);
        low_pII -> Draw();
        up_pII  -> Draw("same");
        _leg     -> Draw();
      }
}












       
void All_2DPlots() // CORRECTO: la forma adecuada de hacerlo
{

   for (int s = 0; s < nsample; s++){
 
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_DeltaPhiRll", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_Ht", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_MR", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_MT2bb", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_MT2lblb", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_Mdr", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_R2", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_Rpt", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphijet1met", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphill", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphillmet", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphilmet1", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphimetjet", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dphimetptbll", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_dyll", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_htjets", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_m2l", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_meff", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_metPfType1", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_metPfType1Phi", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_met_meff", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_met_sqrtHt", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_mllbb", "");
   d2Plot ( 350, 225, 2, _sample[s], "MT2ll_ptbll", "");
  
}





}     
