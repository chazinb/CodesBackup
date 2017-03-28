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

void saveHistoRoot ( TString rootfileName, TH1F* histoName){
     
     TFile* OutFileName = new TFile( rootfileName, "recreate");
     histoName   -> Write();
     OutFileName -> Close();
}

void drawShape (TH1F *histoName){
    float integral;
    integral = histoName -> Integral(0, histoName -> GetNbinsX()+1);
    histoName -> Scale(1/integral);
}
 
                  
void shape_study(float csvv2mcut = 1, float jet2ptcut = -9999., float metcut = 80., float V4cut = 5.) {


  // Minitrees address / Cuts parameters

     TString FileAddress = "../../minitrees/nominal/Stop/"; 
     int nbck  = 3; 
     const int ncut        = 4; 
     TString scut[ncut]    = {"20csvv2m", "jet2pt", "80met", "5met_sqrtHt"}; 
//"20csvv2m" , "jet2pt", "0.28 met_meff", "5met_sqrtHt", "025dphimetjet"};  
     float valueCut [ncut] = { csvv2mcut, jet2ptcut, metcut, V4cut};


  // Declare your minitrees

     TString SampleName [1] = {"T2tt_mStop-350to400.root"};// Sample mass point [1][2] = {(350,175)}
     TString BckName    [3] = {"TTTo2L2Nu_ext1.root", "DYJetsToLL.root", "WWTo2L2Nu.root"};
     TString LegName    [3] = {"tt", "DY", "WW"};

  // Declare and inintialize your variable 's histograms

     TH1F *S_jet2Pt [ncut],    * S_met [ncut],     *S_dphimetjet [ncut],     *S_met_sqrtHt [ncut],    *S_met_meff [ncut];  
     TH1F *B_jet2Pt [ncut][3], * B_met [ncut][3],  *B_dphimetjet [ncut][3],  *B_met_sqrtHt [ncut][3], *B_met_meff [ncut][3];
     

     for ( int i = 0; i < ncut; i++) { 
 	S_jet2Pt     [i] = new TH1F ("jet2Pt_Stop_"     + scut[i],  "", 40 , 0. , 200.);// [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
        S_met        [i] = new TH1F ("Met_Stop_"        + scut[i],  "", 20 , 0. , 400.);      
        S_dphimetjet [i] = new TH1F ("dphimetjet_Stop_" + scut[i],  "", 80 , 0. , 4.  );
        S_met_meff   [i] = new TH1F ("met_meff_Stop"    + scut[i],  "", 40, 0., 0.8  );
        S_met_sqrtHt [i] = new TH1F ("met_sqrtHt_Stop_" + scut[i],  "", 25 , 0. , 25. );
        
        for (int j = 0; j < 3; j++){
	   B_jet2Pt     [i][j] = new TH1F ("jet2Pt_"     + scut[i] + "_"+ BckName [j], "", 40 , 0., 200.);
           B_met        [i][j] = new TH1F ("Met_"        + scut[i] + "_"+ BckName [j], "", 20 , 0., 400.);
           B_dphimetjet [i][j] = new TH1F ("dphimetet_"  + scut[i] + "_"+ BckName [j], "", 80 , 0., 4.  );
           B_met_meff   [i][j] = new TH1F ("met_meff_"   + scut[i] + "_"+ BckName [j], "", 40, 0., 0.8  );
           B_met_sqrtHt [i][j] = new TH1F ("met_sqrtHt_" + scut[i] + "_"+ BckName [j], "", 25 , 0. , 25.);
       }
     }    
 

  // Open the minitrees for your signal/ Fill the histo counter/ Save it as .root histogram

 
     TFile *S_MiniTreeFile  = TFile::Open( FileAddress +"/T2tt_mStop-350to400.root"); 
     TTree *S_MiniTree      = GetMiniTree(S_MiniTreeFile);
     Int_t  S_nentries      = (Int_t) S_MiniTree -> GetEntries();

     for (Int_t i = 0; i<S_nentries; i++){

       S_MiniTree -> GetEntry(i);

       if (Mstop != 350 || Mlsp != 225) continue; 
      
       //htjets = ht - metPfType1 - lep1pt - lep1pt; 
       float V3 = metPfType1/meff; 
       float V4 = metPfType1/sqrt(htjets);

       bool pass;
       int j; 
       for ( j = 0; j < ncut; j++){ 
         
         if (j==0) { pass =  nbjet20csvv2m >= valueCut [0];}
         if (j==1) { pass &= jet2pt        >= valueCut [1];}
         if (j==2) { pass &= metPfType1    >= valueCut [2];}
         //pass &= V3            >= valueCut [2];}
         if (j==3) { pass &= V4            >= valueCut [3];}
         if (pass) {
	   S_jet2Pt [j] -> Fill (jet2pt, eventW); 
	   S_met [j] -> Fill ( metPfType1, eventW); 
	   S_dphimetjet [j] -> Fill (dphimetjet,  eventW); 
	   S_met_sqrtHt [j] -> Fill ( V4, eventW);                              
	   S_met_meff[j] -> Fill (V3, eventW);
	 }

       } 

     }

    // Open the minitrees for your background/ Fill the histo counters/ Save them as .root histogram

    for (int j = 0; j < nbck ; j++){

     TFile *B_MiniTreeFile  = TFile::Open(FileAddress + BckName[j]);
     TTree *B_MiniTree      = GetMiniTree(B_MiniTreeFile);
     Int_t  B_nentries      = (Int_t) B_MiniTree -> GetEntries();

     for (Int_t i = 0; i<B_nentries; i++){

       B_MiniTree -> GetEntry(i);
       
       //htjets = ht - metPfType1 - lep1pt - lep1pt;
       float V3 = metPfType1/meff;
       float V4 = metPfType1/sqrt(htjets);

       bool pass; 
       for (int k = 0; k < ncut; k++ ){
          
          if (k==0) { pass =  nbjet20csvv2m >= valueCut [0];}
          if (k==1) { pass &= jet2pt        >= valueCut [1];}
          if (k==2) { pass &= metPfType1    >= valueCut [2];}
          //if (k==2) { pass &= V3            >= valueCut [2];}
          if (k==3) { pass &= V4            >= valueCut [3];}

	  //if (eventW<0) cout << "ll " << eventW << " " << BckName[j] << endl;
	  if (pass) { 
	    B_jet2Pt[k][j] -> Fill (jet2pt, eventW); 
	    B_met[k][j] -> Fill ( metPfType1, eventW); 
	    B_dphimetjet[k][j] -> Fill (dphimetjet,  eventW); 
	    B_met_sqrtHt[k][j] -> Fill ( V4, eventW); 
	    B_met_meff[k][j]  -> Fill (V3, eventW); 
	  }
       }
       
     }
    }
    

    // TString OutFileName;

    /*
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
    */
    //}

    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001");

    for (int c = 0; c<ncut; c++) {

      float SignalYields = S_met[c]->Integral(0, S_met[c]->GetNbinsX()+1), BackgroundYield = 0.;
      cout << "Cut " << scut[c] << " " << valueCut[c] << " Signal yields " << SignalYields << endl;

      for (int b = 0; b<3; b++) {

	//float ThisBackYield = B_jet2Pt[c][b]->GetEntries();//Integral(0, B_jet2Pt[c][b]->GetNbinsX()+1);
	float ThisBackYield = B_met[c][b]->Integral(0, B_met[c][b]->GetNbinsX()+1);
	float Signb = SignalYields / sqrt(SignalYields + ThisBackYield);
	cout << BckName[b] << " " << ThisBackYield << " " << Signb << endl;
	BackgroundYield += ThisBackYield;
    
      }

      float Signt = SignalYields / sqrt(SignalYields + BackgroundYield);
      cout << "Total B" << " " << BackgroundYield << " " << Signt << endl <<  endl;

    }

   // Draw and save .png the  histos

    float integral;  

    for (int c=0; c<ncut; c++){ 
      // S_jet2Pt / B_jet2Pt

      
      TCanvas *CC1 = new TCanvas("CC", "", 1200, 1000);
      gStyle->SetOptStat(0000);
      TLegend*  legi = myLegend (" Hola"); 
      
      for (int i = 0; i < 3; i++){
        int x = 2+i; 
        legi->AddEntry(B_jet2Pt[c][i], LegName [i] , "f");
        B_jet2Pt[c][i]->SetLineColor(x); B_jet2Pt[c][i]->SetLineWidth(2);           
        drawShape(B_jet2Pt[c][i]);
      }
      
      
      legi->AddEntry(S_jet2Pt[c], scut[c]+"stop", "f");
      S_jet2Pt[c]->SetLineColor(7); S_jet2Pt[c]->SetLineWidth(2);
      drawShape ( S_jet2Pt[c]);
      
      B_jet2Pt[c][0] -> SetMaximum(0.15);
      B_jet2Pt[c][0] -> SetTitle ("jet2pt" + scut[c]);
      B_jet2Pt[c][0] -> SetXTitle ("trailing jet pt [GeV]");
      B_jet2Pt[c][0] -> DrawCopy("histo");
      B_jet2Pt[c][1] -> DrawCopy("histsame");
      B_jet2Pt[c][2] -> DrawCopy("histsame");
      S_jet2Pt[c]    -> DrawCopy("histsame");
      legi->Draw();
      CC1->Print("shape_study/jet2Pt_" + scut[c]+ "_.png");
      
      
      //
      // S_met / B_met
      
      TCanvas *CC2  = new TCanvas("CC2", "", 1200, 1000);
      gStyle->SetOptStat(0000);
      TLegend *leg2 = myLegend ("");
      
      for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg2->AddEntry(B_met[c][i], LegName [i] , "f");
        B_met[c][i]->SetLineColor(x); B_met[c][i]->SetLineWidth(2);
        drawShape(B_met[c][i]);
      }

      
      leg2->AddEntry(S_met[c], "(Mstop = 350, Mlsp = 225)", "f");
      S_met[c]->SetLineColor(7); S_met[c]->SetLineWidth(2);
      drawShape ( S_met[c]);
      
      B_met[c][0] -> SetMaximum(0.6);
      B_met[c][0] -> SetTitle ("met");
      B_met[c][0] -> SetXTitle ("metPfType1 [GeV]");
      B_met[c][0] -> DrawCopy("histo");
      B_met[c][1] -> DrawCopy("histsame");
      B_met[c][2] -> DrawCopy("histsame");
      S_met[c]    -> DrawCopy("histsame");
      leg2->Draw();
      CC2->Print("shape_study/metPfType1_" + scut[c]+ "_.png");
      
      // S_dphimetjet / B_dphimetjet
      
      TCanvas *CC3  = new TCanvas("CC3", "", 1200, 1000);
      TLegend *leg3 = myLegend("");
      
      for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg3->AddEntry(B_dphimetjet[c][i], LegName [i] , "f");
        B_dphimetjet[c][i]->SetLineColor(x); B_dphimetjet[c][i]->SetLineWidth(2);
        drawShape ( B_dphimetjet[c][i]);
      }
      
      
      leg3->AddEntry(S_met[c], "(Mstop = 350, Mlsp = 225)", "f");
      S_dphimetjet[c]->SetLineColor(7); S_dphimetjet[c]->SetLineWidth(2);
      drawShape (S_dphimetjet[c]);
      
      
      B_dphimetjet[c][0] -> SetMaximum(0.12);
      B_dphimetjet[c][0] -> SetTitle ("dphimetjet");
      B_dphimetjet[c][0] -> SetXTitle ("dphimetjet [rad]");
      B_dphimetjet[c][0] -> DrawCopy("histo");
      B_dphimetjet[c][1] -> DrawCopy("histsame");
      B_dphimetjet[c][2] -> DrawCopy("histsame");
      S_dphimetjet[c]    -> DrawCopy("histsame");
      leg3->Draw();
      CC3->Print("shape_study/dphimetjet_" + scut[c]+ "_.png");
      
      // S_met_sqrtHt / B_met_sqrtHt
      
      TCanvas *CC4  = new TCanvas("CC4", "", 1200, 1000);
      gStyle->SetOptStat(0000);
      TLegend *leg4 = myLegend("");
      
      for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg4->AddEntry(B_met_sqrtHt[c][i], LegName [i] , "f");
        B_met_sqrtHt[c][i]->SetLineColor(x); B_met_sqrtHt[c][i]->SetLineWidth(2);
        drawShape( B_met_sqrtHt[c][i]);
      }
      
      leg4->AddEntry(S_met_sqrtHt[c], "(Mstop = 350, Mlsp = 225)", "f");
      S_met_sqrtHt[c]->SetLineColor(7); S_met_sqrtHt[c]->SetLineWidth(2);
      drawShape(S_met_sqrtHt[c]);
      
      B_met_sqrtHt[c][0] -> SetMaximum(0.30);
      B_met_sqrtHt[c][0] -> SetTitle ("met_sqrtHt");
      B_met_sqrtHt[c][0] -> SetXTitle ("met_sqrtHt");
      B_met_sqrtHt[c][0] -> DrawCopy("histo");
      B_met_sqrtHt[c][1] -> DrawCopy("histsame");
      B_met_sqrtHt[c][2] -> DrawCopy("histsame");
      S_met_sqrtHt[c]    -> DrawCopy("histsame");
      leg4->Draw();
      CC4->Print("shape_study/met_sqrtHt_" + scut[c]+ "_.png");
      
      // S_met_meff / B_met_meff
      
      TCanvas *CC5  = new TCanvas("CC5", "", 1200, 1000);
      gStyle->SetOptStat(0000);
      TLegend *leg5 = myLegend("");
      
      for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg5->AddEntry(B_met_meff[c][i], LegName [i] , "f");
        B_met_meff[c][i]->SetLineColor(x); B_met_meff[c][i]->SetLineWidth(2);
        drawShape( B_met_meff[c][i]);
      }
      
      leg5->AddEntry(S_met_meff[c], "(Mstop = 350, Mlsp = 225)", "f");
      S_met_meff[c]->SetLineColor(7); S_met_meff[c]->SetLineWidth(2);
      drawShape(S_met_meff[c]);
      
      B_met_meff[c][0] -> SetMaximum(0.24);
      B_met_meff[c][0] -> SetTitle ("met_meff");
      B_met_meff[c][0] -> SetXTitle ("met_meff");
      B_met_meff[c][0] -> DrawCopy("histo");
      B_met_meff[c][1] -> DrawCopy("histsame");
      B_met_meff[c][2] -> DrawCopy("histsame");
      S_met_meff[c]    -> DrawCopy("histsame");
      leg5->Draw();
      CC5->Print("shape_study/met_meff_" + scut[c]+ "_.png");

    }
    
}



                       
