#include "TCanvas.h"
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

#include <fstream>
#include <iostream>

float eventW, metPfType1, metPfType1Phi, pfType1Met, pfType1Metphi; 
float njet, channel, nbjet30csvv2l, nbjet30csvv2m, nbjet30csvv2t;
float mt2ll, mt2bb, mt2bbtrue, mt2lblb, mt2lblbcomb, mt2lblbtrue;
float mlb1, mlb1true, mlb1comb, mlb1truecomb;
float mlb2, mlb2true, mlb2comb, mlb2truecomb;
float tjet1pt, tjet1phi, tjet1eta, tjet1mass, tjet1csvv2ivf, tjet1assignment;
float tjet2pt, tjet2phi, tjet2eta, tjet2mass, tjet2csvv2ivf, tjet2assignment;
float bjet1pt, bjet1phi, bjet1eta, bjet1mass, bjet1csvv2ivf;
float bjet2pt, bjet2phi, bjet2eta, bjet2mass, bjet2csvv2ivf;
float lep1pt, lep1phi, lep1eta, lep2pt, lep2phi, lep2eta;
float neutrino1px, neutrino1py, neutrino1pz, neutrino2px, neutrino2py, neutrino2pz;
float Mstop, Mlsp;

TTree *GetMiniTree(TFile *MiniTreeFile) {

  TTree *MiniTree = (TTree*) MiniTreeFile->Get("latino"); 

  MiniTree->SetBranchAddress("eventW",          &eventW);
  MiniTree->SetBranchAddress("metPfType1",      &metPfType1);
  MiniTree->SetBranchAddress("metPfType1Phi",   &metPfType1Phi);
  //MiniTree->SetBranchAddress("pfType1Met",      &pfType1Met);
  //MiniTree->SetBranchAddress("pfType1Metphi",   &pfType1Metphi);
  MiniTree->SetBranchAddress("njet",            &njet);
  MiniTree->SetBranchAddress("channel",         &channel);

  MiniTree->SetBranchAddress("Mstop",           &Mstop);
  MiniTree->SetBranchAddress("Mlsp",            &Mlsp);

  MiniTree->SetBranchAddress("nbjet30csvv2l",   &nbjet30csvv2l);
  MiniTree->SetBranchAddress("nbjet30csvv2m",   &nbjet30csvv2m);
  MiniTree->SetBranchAddress("nbjet30csvv2t",   &nbjet30csvv2t);

  MiniTree->SetBranchAddress("mt2ll",           &mt2ll);
  MiniTree->SetBranchAddress("mt2bb",           &mt2bb);
  MiniTree->SetBranchAddress("mt2lblb",         &mt2lblb);
  MiniTree->SetBranchAddress("mt2bbtrue",       &mt2bbtrue);
  MiniTree->SetBranchAddress("mt2lblbcomb",     &mt2lblbcomb);
  MiniTree->SetBranchAddress("mt2lblbtrue",     &mt2lblbtrue);

  MiniTree->SetBranchAddress("mlb1",            &mlb1);
  MiniTree->SetBranchAddress("mlb2",            &mlb2);
  MiniTree->SetBranchAddress("mlb1comb",        &mlb1comb);
  MiniTree->SetBranchAddress("mlb2comb",        &mlb2comb);
  MiniTree->SetBranchAddress("mlb1true",        &mlb1true);
  MiniTree->SetBranchAddress("mlb2true",        &mlb2true);
//   MiniTree->SetBranchAddress("mlb1truecomb",    &mlb1truecomb);
//   MiniTree->SetBranchAddress("mlb2truecomb",    &mlb2truecomb);

  MiniTree->SetBranchAddress("bjet1pt",         &bjet1pt);        
  MiniTree->SetBranchAddress("bjet1eta",        &bjet1eta);       
  MiniTree->SetBranchAddress("bjet1phi",        &bjet1phi);       
  MiniTree->SetBranchAddress("bjet1mass",       &bjet1mass);      
  MiniTree->SetBranchAddress("bjet1csvv2ivf",   &bjet1csvv2ivf);  
  MiniTree->SetBranchAddress("bjet2pt",         &bjet2pt);        
  MiniTree->SetBranchAddress("bjet2eta",        &bjet2eta);       
  MiniTree->SetBranchAddress("bjet2phi",        &bjet2phi);       
  MiniTree->SetBranchAddress("bjet2mass",       &bjet2mass);      
  MiniTree->SetBranchAddress("bjet2csvv2ivf",   &bjet2csvv2ivf);  
  MiniTree->SetBranchAddress("tjet1pt",         &tjet1pt);        
  MiniTree->SetBranchAddress("tjet1eta",        &tjet1eta);       
  MiniTree->SetBranchAddress("tjet1phi",        &tjet1phi);       
  MiniTree->SetBranchAddress("tjet1mass",       &tjet1mass);      
  MiniTree->SetBranchAddress("tjet1csvv2ivf",   &tjet1csvv2ivf);  
  MiniTree->SetBranchAddress("tjet1assignment", &tjet1assignment);
  MiniTree->SetBranchAddress("tjet2pt",         &tjet2pt);        
  MiniTree->SetBranchAddress("tjet2eta",        &tjet2eta);       
  MiniTree->SetBranchAddress("tjet2phi",        &tjet2phi);
  MiniTree->SetBranchAddress("tjet2mass",       &tjet2mass);      
  MiniTree->SetBranchAddress("tjet2csvv2ivf",   &tjet2csvv2ivf);  
  MiniTree->SetBranchAddress("tjet2assignment", &tjet2assignment);

  MiniTree->SetBranchAddress("lep1pt",          &lep1pt);
  MiniTree->SetBranchAddress("lep1phi",         &lep1phi);
  MiniTree->SetBranchAddress("lep1eta",         &lep1eta);
  MiniTree->SetBranchAddress("lep2pt",          &lep2pt);
  MiniTree->SetBranchAddress("lep2phi",         &lep2phi);
  MiniTree->SetBranchAddress("lep2eta",         &lep2eta);

//   MiniTree->SetBranchAddress("neutrino1px",     &neutrino1px);
//   MiniTree->SetBranchAddress("neutrino1py",     &neutrino1py);
//   MiniTree->SetBranchAddress("neutrino1pz",     &neutrino1pz);
//   MiniTree->SetBranchAddress("neutrino2px",     &neutrino2px);
//   MiniTree->SetBranchAddress("neutrino2py",     &neutrino2py);
//   MiniTree->SetBranchAddress("neutrino2pz",     &neutrino2pz);

  return MiniTree;

}

