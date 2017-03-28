#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include "TStopwatch.h"

void EventCounter() {

//eosmount  eoscms;
// This function creates a 2D Histogram of float numbers to see the number of event for each Mstop-Mlsp model
TString A = "eoscms/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/LatinoTrees/";
//TString FileName [48] = { 
TString FileName [48] = { 
A + "latino_T2tt_mStop-150to250_0000__part0.root",
A + "latino_T2tt_mStop-150to250_0000__part1.root",
A + "latino_T2tt_mStop-150to250_0000__part2.root",
A + "latino_T2tt_mStop-150to250_0000__part3.root",
A + "latino_T2tt_mStop-150to250_0000__part4.root",
A + "latino_T2tt_mStop-150to250_0000__part5.root",
A + "latino_T2tt_mStop-150to250_0000__part6.root",
A + "latino_T2tt_mStop-150to250_0000__part7.root",
A + "latino_T2tt_mStop-150to250_0000__part8.root",
A + "latino_T2tt_mStop-150to250_0000__part9.root",
A + "latino_T2tt_mStop-150to250_0001__part0.root",
A + "latino_T2tt_mStop-150to250_0001__part1.root",
A + "latino_T2tt_mStop-150to250_0001__part2.root",
 
A + "latino_T2tt_mStop-250to350_0000__part0.root",
A + "latino_T2tt_mStop-250to350_0000__part1.root",
A + "latino_T2tt_mStop-250to350_0000__part2.root",
A + "latino_T2tt_mStop-250to350_0000__part3.root",
A + "latino_T2tt_mStop-250to350_0000__part4.root",
A + "latino_T2tt_mStop-250to350_0000__part5.root",
A + "latino_T2tt_mStop-250to350_0000__part6.root",
A + "latino_T2tt_mStop-250to350_0000__part7.root",
A + "latino_T2tt_mStop-250to350_0000__part8.root",
A + "latino_T2tt_mStop-250to350_0000__part9.root",
A + "latino_T2tt_mStop-250to350_0001__part0.root",
A + "latino_T2tt_mStop-250to350_0001__part1.root",
A + "latino_T2tt_mStop-250to350_0001__part2.root",
 
A + "latino_T2tt_mStop-350to400_0000__part0.root",
A + "latino_T2tt_mStop-350to400_0000__part1.root",
A + "latino_T2tt_mStop-350to400_0000__part2.root",
A + "latino_T2tt_mStop-350to400_0000__part3.root",
A + "latino_T2tt_mStop-350to400_0000__part4.root",
A + "latino_T2tt_mStop-350to400_0000__part5.root",
A + "latino_T2tt_mStop-350to400_0000__part6.root",
A + "latino_T2tt_mStop-350to400_0000__part7.root",
A + "latino_T2tt_mStop-350to400_0000__part8.root",
A + "latino_T2tt_mStop-350to400_0000__part9.root",
A + "latino_T2tt_mStop-350to400_0001__part0.root",
A + "latino_T2tt_mStop-350to400_0001__part1.root",
 
A + "latino_T2tt_mStop-400to1200__part0.root", 
A + "latino_T2tt_mStop-400to1200__part1.root", 
A + "latino_T2tt_mStop-400to1200__part2.root", 
A + "latino_T2tt_mStop-400to1200__part3.root", 
A + "latino_T2tt_mStop-400to1200__part4.root", 
A + "latino_T2tt_mStop-400to1200__part5.root", 
A + "latino_T2tt_mStop-400to1200__part6.root", 
A + "latino_T2tt_mStop-400to1200__part7.root", 
A + "latino_T2tt_mStop-400to1200__part8.root", 
A + "latino_T2tt_mStop-400to1200__part9.root"};

TH2F* MassPlane = new TH2F ("S_LSP_MassPlane", "", 1100, 150,1250,1250,0,1250);
float susyMstop, susyMLSP;

//std::ofstream inFile("Mass_table.txt",std::ios::out);

for ( int i = 0; i<48; i++){
  TFile *latino_File = TFile::Open(FileName[i]);
  TTree *latino_tree = (TTree*) latino_File->Get("latino");
//  inFile << "FileName = " << FileName[i] <<endl;
  cout << "FileName = " << FileName[i] <<endl;
//  TStopwatch t;
//  t.Start();
  
  Int_t nentries = (Int_t) latino_tree->GetEntries();
  latino_tree->SetBranchAddress("susyMstop",          &susyMstop);
  latino_tree->SetBranchAddress("susyMLSP",           &susyMLSP);
  
 // inFile << "nentries = " << nentries << "\n" << endl; 
  cout << "nentries = " << nentries << "\n" << endl; 

 // for (int j = 0; j<100; j++){
  for (int j = 0; j<nentries; j++){
      	
      latino_tree -> GetEntry(i);
      MassPlane -> Fill(susyMstop,susyMLSP);
   //   inFile << "susyMstop =" << susyMstop << "\n" <<  "susyMLSP =" << susyMLSP << endl; 
   }
 // t.Stop();
  //t.Print();
 // printf("RealTime=%f seconds, CpuTime=%f seconds\n",t.RealTime(),t.CpuTime()); 
}
//inFile.close();
TFile *OutFile = new TFile("histoMassPlane.root", "recreate");
MassPlane->Write();
OutFile->Close();
TCanvas *CC = new TCanvas();
//gStyle->SetOptStat("");
//TPad *PD1 = (TPad*)CC->GetPad(1); PD1->SetLogy(); PD1->SetGridx(); PD1->SetGridy();
//PD1->cd();
MassPlane->Draw(); 
CC->Print("./Plots/MassPlane.png");
}
