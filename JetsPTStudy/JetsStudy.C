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

void MyLegend (TLegend *leg1, TString title ){

    leg1 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg1->SetFillColor(kWhite); leg1->SetBorderSize(0.);
    leg1->SetTextColor(1); leg1->SetTextSize(0.020);
    leg1->SetTextFont(50);
    //leg1->SetHeader(svar[vv] + "  " + lHeader[R]); //leg1_fit->SetMargin(0.2); 
    leg1->AddEntry((TObject*)0, " ", "");
    
     }
                   
void JetsStudy() {


  // Minitrees address / Counters parameters

     TString FileAddress = "../../minitreesPT/nominal/Stop/"; 
     int nbck  = 3; 
     int nCuts = 12; 
     TString sCuts[12] = { "20csvv2m", "20csvv2m_jet2pt > 20", "20csvv2m_jet2pt > 25", "20csvv2m_jet2pt > 30","25csvv2m", "25csvv2m_jet2pt > 20", "25csvv2m_jet2pt > 25", "25csvv2m_jet2pt > 30", "30csvv2m", "30csvv2m_jet2pt > 20", "30csvv2m_jet2pt > 25", "30csvv2m_jet2pt > 30"};


  // Declare your minitrees

     TString SampleName [1] = {"T2tt_mStop-350to400.root"};
     // Sample mass point [1][2] = {(350,175)}
     TString BckName    [3] = {"TTTo2L2Nu_ext1.root", "DYJetsToLL.root", "WWTo2L2Nu.root"};
     TString LegName    [3] = {"tt", "DY", "WW"};

  // Declare and inintialize your variable 's histograms

     TH1F* S_jet2Pt; S_jet2Pt = new TH1F ("jet2Pt_Stop","", 60, 0., 300.);// [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
     TH1F* B_jet2Pt [3]; for (int i = 0; i < 3; i++){ B_jet2Pt [i] = new TH1F ("jet2Pt" + BckName [i],"",60,  0., 300.);}
     
     TH1F* S_met; S_met = new TH1F ("Met_Stop","", 80, 0., 800.); // [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
     TH1F* B_met [3]; for (int i = 0; i < 3; i++){ B_met [i] = new TH1F ("Met" + BckName [i],"", 80, 0., 800.);}
     
     TH1F* S_dphimetjet; S_dphimetjet = new TH1F ("dphimetjet_Stop","", 80, 0., 4.);// [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
     TH1F* B_dphimetjet [3]; for (int i = 0; i < 3; i++){ B_dphimetjet [i] = new TH1F ("dphimetjet" + BckName [i],"", 80, 0., 4.);}
      
     TH1F* S_met_sqrtHt; S_met_sqrtHt = new TH1F ("met_sqrtHt_Stop","",60, 0., 15.);// [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
     TH1F* B_met_sqrtHt [3]; for (int i = 0; i < 3; i++){ B_met_sqrtHt [i] = new TH1F ("met_sqrtHt" + BckName [i],"",60, 0., 15.);}

     TH1F* S_met_Ht; S_met_Ht = new TH1F ("met_Ht_Stop","",70, 0., 0.7);// [3]; for (int i = 0; i < 4; i++){ jet2Pt [i] = new THF1 ("jet2Pt" + SampleName ,"", 0., 300.);}
     TH1F* B_met_Ht [3]; for (int i = 0; i < 3; i++){ B_met_Ht [i] = new TH1F ("met_Ht" + BckName [i],"",70, 0., 0.7);}

     TH1F* S_met_meff; S_met_meff = new TH1F ("met_meff_Stop","",70, 0., 0.7);
     
     TH1F* B_met_meff [3]; for (int i = 0; i < 3; i++){ B_met_meff [i] = new TH1F ("met_meff" + BckName [i],"",70, 0., 0.7);}
 
  // Declare and initialize your event counters
    
    // 1 cut == 1 bin -> 17/10/16 -> 9 cuts == 9 bins

     TH1F* S_hcounter;     
     TH1F* B_hcounter    [4]; 
     TH1F* significance1 [4];
     TH1F* significance2 [4]; 

     S_hcounter        = new TH1F ("Signal_hcounter",          "",12 , 0., 12.); // Signal event counter
     B_hcounter    [3] = new TH1F ("Background_hcounter",      "",12 , 0., 12.); // Total Background event counter  
     significance1 [3] = new TH1F ("S_B_Significance",         "",12 , 0., 12.); // Total S/B  event counter
     significance2 [3] = new TH1F ("S_sqrt(B+S)_Significance", "",12 , 0., 12.); // Total S/ sqrt(S+B) event counter

     for (int i = 0; i < 3; i++){

        B_hcounter [i]    = new TH1F ("Background_hcounter"      + BckName [i],  "",12, 0., 12.); // Each background event counter
	significance1 [i] = new TH1F ("S_B_Significance"         + BckName [i],  "",12, 0., 12.); //       ""           ""
	significance2 [i] = new TH1F ("S_sqrt(B+S)_Significance" + BckName [i],  "",12, 0., 12.); //       ""           ""
     }
     // Set the Xaxis label
    
     for (int i = 1; i < 13; i++){ 
	
       S_hcounter -> GetXaxis()-> SetBinLabel (i, sCuts[i-1]); 
       
       for (int j = 0; j < 4; j++){

         B_hcounter    [j]  -> GetXaxis()-> SetBinLabel (i, sCuts[i-1]); 
	 significance1 [j]  -> GetXaxis()-> SetBinLabel (i, sCuts[i-1]); 
	 significance2 [j]  -> GetXaxis()-> SetBinLabel (i, sCuts[i-1]);
	}
      }

  // Open the minitrees for your signal/ Fill the histo counter/ Save it as .root histogram

     TFile *S_MiniTreeFile  = TFile::Open("../../minitreesPT/nominal/Stop/T2tt_mStop-350to400.root"); 
     TTree *S_MiniTree      = GetMiniTree(S_MiniTreeFile);
     Int_t  S_nentries      = (Int_t) S_MiniTree -> GetEntries();
    
     for (Int_t i = 0; i<S_nentries; i++){

       S_MiniTree -> GetEntry(i);

       if (Mstop != 350 || Mlsp != 225) continue; 

       V3 = metPfType1/meff;
       V4 = metPfType1/sqrt(ht);

       if (nbjet20csvv2m >= 1) { S_hcounter-> Fill(0.,  eventW); S_jet2Pt -> Fill (jet2pt, eventW); S_met -> Fill ( metPfType1, eventW); S_dphimetjet -> Fill (dphimetjet,  eventW); S_met_sqrtHt -> Fill ( V4, eventW); S_met_meff -> Fill (V3, eventW); S_met_Ht -> Fill (ht, eventW); } 
       
       if (nbjet20csvv2m >= 1 && jet2pt > 20)   S_hcounter-> Fill(1,   eventW); 
       if (nbjet20csvv2m >= 1 && jet2pt > 25)   S_hcounter-> Fill(2,   eventW);
       if (nbjet20csvv2m >= 1 && jet2pt > 30)   S_hcounter-> Fill(3,   eventW);      
 
       if (nbjet25csvv2m >= 1)                  S_hcounter-> Fill (4,  eventW);     
       if (nbjet25csvv2m >= 1 && jet2pt > 20)   S_hcounter-> Fill (5,  eventW);
       if (nbjet25csvv2m >= 1 && jet2pt > 25)   S_hcounter-> Fill (6,  eventW);
       if (nbjet25csvv2m >= 1 && jet2pt > 30)   S_hcounter-> Fill (7,  eventW);

       if (nbjet30csvv2m >= 1)                  S_hcounter-> Fill (8,  eventW);     
       if (nbjet30csvv2m >= 1 && jet2pt > 20)   S_hcounter-> Fill (9,  eventW);
       if (nbjet30csvv2m >= 1 && jet2pt > 25)   S_hcounter-> Fill (10, eventW);
       if (nbjet30csvv2m >= 1 && jet2pt > 30)   S_hcounter-> Fill (11, eventW);

     }

     TFile *OutFile1 = new TFile("S_hcounter_PtCuts.root", "recreate");
     S_hcounter -> Write();
     OutFile1 -> Close();


    // Open the minitrees for your background/ Fill the histo counters/ Save them as .root histogram

    for (int j = 0; j < nbck ; j++){

     TFile *B_MiniTreeFile  = TFile::Open(FileAddress + BckName[j]);
     TTree *B_MiniTree      = GetMiniTree(B_MiniTreeFile);
     Int_t  B_nentries      = (Int_t) B_MiniTree -> GetEntries();

     for (Int_t i = 0; i<B_nentries; i++){

       B_MiniTree -> GetEntry(i);

       V3 = metPfType1/meff;
       V4 = metPfType1/sqrt(ht);
 
       if (nbjet20csvv2m >= 1) {B_hcounter[j]-> Fill (0., eventW); B_hcounter[3]-> Fill (0., eventW); B_jet2Pt[j] -> Fill (jet2pt, eventW); B_met[j] -> Fill ( metPfType1, eventW); B_dphimetjet[j] -> Fill (dphimetjet,  eventW); B_met_sqrtHt[j] -> Fill ( V4, eventW); B_met_meff [j]-> Fill (V3, eventW); B_met_Ht [j] -> Fill (ht, eventW);} 
   //    if (nbjet20csvv2m >= 1) {B_hcounter[j]-> Fill (0., eventW); B_hcounter[3]-> Fill (0., eventW);}
       if (nbjet20csvv2m >= 1 && jet2pt > 20) { B_hcounter[j]-> Fill (1,  eventW); B_hcounter[3]-> Fill (1,  eventW);}
       if (nbjet20csvv2m >= 1 && jet2pt > 25) { B_hcounter[j]-> Fill (2,  eventW); B_hcounter[3]-> Fill (2,  eventW);}
       if (nbjet20csvv2m >= 1 && jet2pt > 30) { B_hcounter[j]-> Fill (3,  eventW); B_hcounter[3]-> Fill (3,  eventW);}      
 
       if (nbjet25csvv2m >= 1)                { B_hcounter[j]-> Fill (4,  eventW); B_hcounter[3]-> Fill (4,  eventW);}
       if (nbjet25csvv2m >= 1 && jet2pt > 20) { B_hcounter[j]-> Fill (5,  eventW); B_hcounter[3]-> Fill (5,  eventW);}
       if (nbjet25csvv2m >= 1 && jet2pt > 25) { B_hcounter[j]-> Fill (6,  eventW); B_hcounter[3]-> Fill (6,  eventW);}
       if (nbjet25csvv2m >= 1 && jet2pt > 30) { B_hcounter[j]-> Fill (7,  eventW); B_hcounter[3]-> Fill (7,  eventW);}

       if (nbjet30csvv2m >= 1)                { B_hcounter[j]-> Fill (8,  eventW); B_hcounter[3]-> Fill (8,  eventW);}
       if (nbjet30csvv2m >= 1 && jet2pt > 20) { B_hcounter[j]-> Fill (9,  eventW); B_hcounter[3]-> Fill (9,  eventW);}
       if (nbjet30csvv2m >= 1 && jet2pt > 25) { B_hcounter[j]-> Fill (10, eventW); B_hcounter[3]-> Fill (10, eventW);}
       if (nbjet30csvv2m >= 1 && jet2pt > 30) { B_hcounter[j]-> Fill (11, eventW); B_hcounter[3]-> Fill (11, eventW);}// B_jet2Pt[j] -> Fill (jet2pt, eventW); B_met[j] -> Fill ( metPfType1, eventW); B_dphimetjet[j] -> Fill (dphimetjet,  eventW); B_met_sqrtHt[j] -> Fill ( V4, eventW); }


      }
    }


     TFile *OutFile2 = new TFile("Background_hcounter_PtCuts.root", "recreate");
     B_hcounter [3] -> Write();
     OutFile2 -> Close();
     
     TString OutFileName;
 
     for (int i = 0; i < nbck; i++){

        OutFileName = " OutFile " + BckName [i];
        TFile* OutFileName = new TFile(BckName [i] +"_hcounter_PtCuts.root", "recreate");
        B_hcounter[i] -> Write();
        OutFileName -> Close();
     }

    // Fill significance* histograms / Save them as .root histogram
 
      for (int j = 0; j < 4; j++){
     	for (int i = 1; i < 13; i++){

       	  float signal  = S_hcounter    -> GetBinContent(i);
          float bck     = B_hcounter[j] -> GetBinContent(i);
       
          float signif1 = signal/bck;
          significance1[j] -> SetBinContent (i, signif1);
       
          float signif2 = signal/sqrt(signal + bck);
          significance2[j] -> SetBinContent (i, signif2);
       }
     }
 
     TFile *OutFile3 = new TFile("S_B_Significance.root", "recreate");
     significance1 [3] -> Write();
     OutFile3 -> Close();

     TFile *OutFile4 = new TFile("S_sqrt(B+S)_Significance.root", "recreate");
     significance2 [3] -> Write();
     OutFile4 -> Close();

     TString OutFileName1;
     TString OutFileName2;   

     for (int i = 0; i < nbck; i++){
    
       OutFileName1 = " OutFile1 " + BckName [i];
       OutFileName2 = " OutFile2 " + BckName [i];

       TFile* OutFileName1 = new TFile(BckName [i] +"_S_B_Significance.root", "recreate");
       significance1 [i] -> Write();
       OutFileName1 -> Close();

       TFile* OutFileName2 = new TFile(BckName [i] +"_S_sqrt(B+S)_Significance.root", "recreate");
       significance2 [i] -> Write();
       OutFileName2 -> Close();
     }
   
   // Draw and save .png the  histos

    //CCounter Total background -> Draw()
    TCanvas *CCounter_total = new TCanvas("CCounter_total", "", 1200, 1000);
    gStyle -> SetOptStat("");
 
    TLegend *zzz;
    zzz = new TLegend(0.75, 0.9, 0.9, 0.7);
    zzz -> SetFillColor(kWhite); zzz->SetBorderSize(0.);
    zzz -> SetTextColor(1); zzz ->SetTextSize(0.020);
    zzz -> SetTextFont(50);
    zzz -> AddEntry((TObject*)0, " ", "");
    zzz -> AddEntry(B_hcounter[3], "All Backgrounds", "f" );
    zzz -> AddEntry(S_hcounter,    "(mS=350,mX=225)*100", "f");
    
    B_hcounter [3] -> SetMinimum(0.);
    B_hcounter [3] -> SetTitle ( "BckgTotal_CutsEvents");
    B_hcounter [3] -> SetFillColor(39); 
    S_hcounter     -> SetFillColor(5);
    S_hcounter     -> Scale(100);
    B_hcounter [3] -> Draw ("hist");
    S_hcounter     -> Draw("histsame");
    zzz            -> Draw();
   
    CCounter_total -> Print("Pt_Study/R4/BckgTotal_CutsEvents.png");
/*  
    //CCounter Each background -> Draw() 
    TCanvas *CCounter_each = new TCanvas("CCounter_each", "", 1200, 1000);
    gStyle -> SetOptStat("");
    
    TLegend *hhh;
    hhh = new TLegend(0.75, 0.9, 0.9, 0.7);
    hhh->SetFillColor(kWhite); hhh->SetBorderSize(0.);
    hhh->SetTextColor(1); hhh->SetTextSize(0.020);
    hhh->SetTextFont(50);
    hhh->AddEntry((TObject*)0, " ", "");
    for (int i = 0; i < 2; i++) { hhh -> AddEntry( B_hcounter [i] , LegName [i] , "f");}
    hhh -> AddEntry( B_hcounter [2] , LegName [2] + "*100" , "f");
    hhh -> AddEntry(S_hcounter, "(mS=350,mX=175)*100", "f");
    
    B_hcounter [0] -> SetLineColor(39); B_hcounter [0] -> SetLineWidth(5);
    B_hcounter [1] -> SetLineColor(8);  B_hcounter [1] -> SetLineWidth(5);
    B_hcounter [2] -> SetLineColor(46); B_hcounter [2] -> SetLineWidth(5);
    S_hcounter     -> SetLineColor(5);  S_hcounter -> SetLineWidth(5);
 
    B_hcounter [0] -> SetMinimum(0.);
    B_hcounter [0] -> SetTitle ("EachBckg_CutsEvents");
    
    B_hcounter [0] -> Draw ("histo");           
    B_hcounter [1] -> DrawCopy ("histsame");
    B_hcounter [2] -> Scale (100);   
    B_hcounter [2] -> DrawCopy ("histsame");  
    S_hcounter     -> DrawCopy ("histsame");       
    hhh            -> Draw(); 
    
    CCounter_each -> Print("Pt_Study/R2/EachBckg_CutsEvents.png");

   //significance2 -> Draw();
 
    //CSignificance2 Total background -> Draw()
    TCanvas *CSignif2_total = new TCanvas("CSignif2_total", "", 1200, 1000);
    gStyle -> SetOptStat("");
    significance2 [3] -> Draw(); 
    significance2 [3] -> SetFillColor(1);
    significance2 [3] -> SetTitle  ("S/sqrt(S + B) total in Jets Pt Study");
    significance2 [3] -> SetYTitle ("S/sqrt(S + B)");
    CSignif2_total    -> Print ("Pt_Study/R2/CSignif2_total.png");

    //CSignificance2 Each background -> Draw()
    TCanvas *CSignif2_each = new TCanvas("CSignif2_each", "", 1200, 1000);
    gStyle -> SetOptStat("");

    TLegend *sss;
    sss = new TLegend(0.75, 0.9, 0.9, 0.7);
    sss->SetFillColor(kWhite); sss->SetBorderSize(0.);
    sss->SetTextColor(1); sss->SetTextSize(0.020);
    sss->SetTextFont(50);
    sss->AddEntry((TObject*)0, " ", "");
    for (int i = 0; i < 3; i++) { sss -> AddEntry( significance2 [i] , LegName [i] , "f");}
    
    significance2 [0] -> SetLineColor(39); significance2 [0] ->SetLineWidth(5);
    significance2 [1] -> SetLineColor(8);  significance2 [1] ->SetLineWidth(5);
    significance2 [2] -> SetLineColor(46); significance2 [2] ->SetLineWidth(5);
    
    significance2 [0] -> SetTitle  ("S/sqrt(S + B) each in Jets Pt Study");
    significance2 [0] -> SetYTitle ("S/sqrt(S + B)");
    significance2 [0] -> SetMinimum (0);
    significance2 [0] -> SetMaximum (5);
    
    significance2 [0] -> Draw("histo");  
    significance2 [1] -> DrawCopy("histsame");  
    significance2 [2] -> DrawCopy("histsame");  
    sss               -> Draw();

    CSignif2_each     -> Print ("Pt_Study/R2/CSignif2_each.png");



   //significance1 -> Draw();
    
    //CSignificance1 Total background -> Draw()

    TCanvas *CSignif1_total = new TCanvas("CSignif1_total", "", 1200, 1000);
    gStyle -> SetOptStat("");
    significance1 [3]-> Draw();
    significance1 [3]-> SetFillColor(1);
    significance1 [3]-> SetTitle  ("S/B total in Jets Pt Study");
    significance1 [3]-> SetYTitle ("S/B");
    CSignif1_total  -> Print ("Pt_Study/R2/CSignif1_total.png");

    //CSignificance1 Each background -> Draw()
    TCanvas *CSignif1_each = new TCanvas("CSignif1_each", "", 1200, 1000);
    gStyle -> SetOptStat("");

    TLegend *xxx;
    xxx = new TLegend(0.75, 0.9, 0.9, 0.7);
    xxx->SetFillColor(kWhite); xxx->SetBorderSize(0.);
    xxx->SetTextColor(1); xxx->SetTextSize(0.020);
    xxx->SetTextFont(50);
    xxx->AddEntry((TObject*)0, " ", "");
    for (int i = 0; i < 3; i++) { xxx -> AddEntry( significance1 [i] , LegName [i] , "f");}

    significance1 [0] -> SetLineColor(39); significance1 [0] ->SetLineWidth(5);
    significance1 [1] -> SetLineColor(8);  significance1 [1] ->SetLineWidth(5);
    significance1 [2] -> SetLineColor(46); significance1 [2] ->SetLineWidth(5);
    significance1 [0] -> SetTitle  ("S/B each in Jets Pt Study");
    significance1 [0] -> SetYTitle ("S/B");
    significance1 [0] -> SetMinimum (0);
    significance1 [0] -> SetMaximum (0.86);
    
    significance1 [0] -> Draw("histo");
    significance1 [1] -> DrawCopy("histsame");
    significance1 [2] -> DrawCopy("histsame");
    xxx               -> Draw();

    CSignif1_each     -> Print ("Pt_Study/R2/CSignif1_each.png");

  
    float integral;  

    // S_jet2Pt / B_jet2Pt

    TCanvas *CC1 = new TCanvas("CC", "", 1200, 1000);
    TLegend  *leg = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
    leg->SetTextColor(1); leg->SetTextSize(0.025);
    leg->SetTextFont(50);
    leg->AddEntry((TObject*)0, " ", "");
    
    for (int i = 0; i < 3; i++){
        int x = 2+i; 
        leg->AddEntry(B_jet2Pt[i], LegName [i] , "f");
        B_jet2Pt[i]->SetLineColor(x); B_jet2Pt[i]->SetLineWidth(2);           
        integral = B_jet2Pt[i]->Integral(0, B_jet2Pt[i]->GetNbinsX()+1);
        B_jet2Pt[i]->Scale(1/integral);
    }


    leg->AddEntry(S_jet2Pt, "(Mstop = 350, Mlsp = 175)", "f");
    S_jet2Pt->SetLineColor(7); S_jet2Pt->SetLineWidth(2);
    integral = S_jet2Pt->Integral(0, S_jet2Pt->GetNbinsX()+1);
    S_jet2Pt->Scale(1/integral);
    
    B_jet2Pt[0] -> SetMaximum(0.12);
    B_jet2Pt[0] -> SetTitle ("jet2pt");
    B_jet2Pt[0] -> SetXTitle ("trailing jet pt [GeV]");
    B_jet2Pt[0] -> DrawCopy("histo");
    B_jet2Pt[1] -> DrawCopy("histsame");
    B_jet2Pt[2] -> DrawCopy("histsame");
    S_jet2Pt    -> DrawCopy("histsame");
    leg->Draw();
    CC1->Print("Pt_Study/R2/jet2Pt.png");
//
    // S_met / B_met

    TCanvas *CC2  = new TCanvas("CC2", "", 1200, 1000);
    TLegend *leg2 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg2->SetFillColor(kWhite); leg2->SetBorderSize(0.);
    leg2->SetTextColor(1); leg2->SetTextSize(0.025);
    leg2->SetTextFont(50);
    leg2->AddEntry((TObject*)0, " ", "");

    for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg2->AddEntry(B_met[i], LegName [i] , "f");
        B_met[i]->SetLineColor(x); B_met[i]->SetLineWidth(2);
        integral = B_met[i]->Integral(0, B_met[i]->GetNbinsX()+1);
        B_met[i]->Scale(1/integral);
    }


    leg2->AddEntry(S_met, "(Mstop = 350, Mlsp = 175)", "f");
    S_met->SetLineColor(7); S_met->SetLineWidth(2);
    integral = S_met->Integral(0, S_met->GetNbinsX()+1);
    S_met->Scale(1/integral);

    B_met[0] -> SetMaximum(0.12);
    B_met[0] -> SetTitle ("met");
    B_met[0] -> SetXTitle ("metPfType1 [GeV]");
    B_met[0] -> DrawCopy("histo");
    B_met[1] -> DrawCopy("histsame");
    B_met[2] -> DrawCopy("histsame");
    S_met    -> DrawCopy("histsame");
    leg2->Draw();
    CC2->Print("Pt_Study/R2/metPfType1.png");

    // S_dphimetjet / B_dphimetjet

    TCanvas *CC3  = new TCanvas("CC3", "", 1200, 1000);
    TLegend *leg3 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg3->SetFillColor(kWhite); leg3->SetBorderSize(0.);
    leg3->SetTextColor(1); leg3->SetTextSize(0.025);
    leg3->SetTextFont(50);
    leg3->AddEntry((TObject*)0, " ", "");

    for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg3->AddEntry(B_dphimetjet[i], LegName [i] , "f");
        B_dphimetjet[i]->SetLineColor(x); B_dphimetjet[i]->SetLineWidth(2);
        integral = B_dphimetjet[i]->Integral(0, B_dphimetjet[i]->GetNbinsX()+1);
        B_dphimetjet[i]->Scale(1/integral);
    }


    leg3->AddEntry(S_met, "(Mstop = 350, Mlsp = 175)", "f");
    S_dphimetjet->SetLineColor(7); S_dphimetjet->SetLineWidth(2);
    integral = S_dphimetjet->Integral(0, S_dphimetjet->GetNbinsX()+1);
    S_dphimetjet->Scale(1/integral);

    B_dphimetjet[0] -> SetMaximum(0.12);
    B_dphimetjet[0] -> SetTitle ("dphimetjet");
    B_dphimetjet[0] -> SetXTitle ("dphimetjet [rad]");
    B_dphimetjet[0] -> DrawCopy("histo");
    B_dphimetjet[1] -> DrawCopy("histsame");
    B_dphimetjet[2] -> DrawCopy("histsame");
    S_dphimetjet    -> DrawCopy("histsame");
    leg3->Draw();
    CC3->Print("Pt_Study/R2/dphimetjet.png");
  
    // S_met_sqrtHt / B_met_sqrtHt

    TCanvas *CC4  = new TCanvas("CC4", "", 1200, 1000);
    TLegend *leg4 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg4->SetFillColor(kWhite); leg4->SetBorderSize(0.);
    leg4->SetTextColor(1);      leg4->SetTextSize(0.025);
    leg4->SetTextFont(50);
    leg4->AddEntry((TObject*)0, " ", "");

    for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg4->AddEntry(B_met_sqrtHt[i], LegName [i] , "f");
        B_met_sqrtHt[i]->SetLineColor(x); B_met_sqrtHt[i]->SetLineWidth(2);
        integral = B_met_sqrtHt[i]->Integral(0, B_met_sqrtHt[i]->GetNbinsX()+1);
        B_met_sqrtHt[i]->Scale(1/integral);
    }


    leg4->AddEntry(S_met_sqrtHt, "(Mstop = 350, Mlsp = 175)", "f");
    S_met_sqrtHt->SetLineColor(7); S_met_sqrtHt->SetLineWidth(2);
    integral = S_met_sqrtHt->Integral(0, S_met_sqrtHt->GetNbinsX()+1);
    S_met_sqrtHt->Scale(1/integral);

    B_met_sqrtHt[0] -> SetMaximum(0.12);
    B_met_sqrtHt[0] -> SetTitle ("met_sqrtHt");
    B_met_sqrtHt[0] -> SetXTitle ("met_sqrtHt");
    B_met_sqrtHt[0] -> DrawCopy("histo");
    B_met_sqrtHt[1] -> DrawCopy("histsame");
    B_met_sqrtHt[2] -> DrawCopy("histsame");
    S_met_sqrtHt    -> DrawCopy("histsame");
    leg4->Draw();
    CC4->Print("Pt_Study/R2/met_sqrtHt.png");


   // S_met_Ht / B_met_Ht
   
    TCanvas *CC5  = new TCanvas("CC5", "", 1200, 1000);
    TLegend *leg5 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg5->SetFillColor(kWhite); leg5->SetBorderSize(0.);
    leg5->SetTextColor(1);      leg5->SetTextSize(0.025);
    leg5->SetTextFont(50);
    leg5->AddEntry((TObject*)0, " ", "");

    for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg5->AddEntry(B_met_Ht[i], LegName [i] , "f");
        B_met_Ht[i]->SetLineColor(x); B_met_Ht[i]->SetLineWidth(2);
        integral = B_met_Ht[i]->Integral(0, B_met_Ht[i]->GetNbinsX()+1);
        B_met_Ht[i]->Scale(1/integral);
    }


    leg5->AddEntry(S_met_Ht, "(Mstop = 350, Mlsp =175)", "f");
    S_met_Ht->SetLineColor(7); S_met_Ht->SetLineWidth(2);
    integral = S_met_Ht->Integral(0, S_met_Ht->GetNbinsX()+1);
    S_met_Ht->Scale(1/integral);

    B_met_Ht[0] -> SetMaximum(0.012);
    B_met_Ht[0] -> SetTitle ("met_Ht");
    B_met_Ht[0] -> SetXTitle ("met_Ht");
    B_met_Ht[0] -> DrawCopy("histo");
    B_met_Ht[1] -> DrawCopy("histsame");
    B_met_Ht[2] -> DrawCopy("histsame");
    S_met_Ht    -> DrawCopy("histsame");
    leg5->Draw();
    CC5->Print("Pt_Study/R2/met_Ht.png");

   // Met/meff 

   TCanvas *CC6  = new TCanvas("CC6", "", 1200, 1000);
    TLegend *leg6 = new TLegend(0.75, 0.9, 0.9, 0.7);
    leg6->SetFillColor(kWhite); leg6->SetBorderSize(0.);
    leg6->SetTextColor(1);      leg6->SetTextSize(0.025);
    leg6->SetTextFont(50);
    leg6->AddEntry((TObject*)0, " ", "");

    for (int i = 0; i < 3; i++){
        int x = 2+i;
        leg6->AddEntry(B_met_meff[i], LegName [i] , "f");
        B_met_meff[i]->SetLineColor(x); B_met_meff[i]->SetLineWidth(2);
        integral = B_met_meff[i]->Integral(0, B_met_meff[i]->GetNbinsX()+1);
        B_met_meff[i]->Scale(1/integral);
    }


    leg6->AddEntry(S_met_meff, "(Mstop = 350, Mlsp = 175)", "f");
    S_met_meff->SetLineColor(7); S_met_meff->SetLineWidth(2);
    integral = S_met_meff->Integral(0, S_met_meff->GetNbinsX()+1);
    S_met_meff->Scale(1/integral);

    B_met_meff[0] -> SetMaximum(0.12);
    B_met_meff[0] -> SetTitle ("met_meff");
    B_met_meff[0] -> SetXTitle ("met_meff");
    B_met_meff[0] -> DrawCopy("histo");
    B_met_meff[1] -> DrawCopy("histsame");
    B_met_meff[2] -> DrawCopy("histsame");
    S_met_meff  -> DrawCopy("histsame");
    leg6->Draw();
    CC6->Print("Pt_Study/R2/met_meff.png");

*/

 }



                       
