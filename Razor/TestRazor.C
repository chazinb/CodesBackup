#include "GetMiniTree.C"
#include "Razor.C"
#include "SuperRazor.C"

void TestRazor() {

  TString FileName[2] = {"/afs/cern.ch/user/b/bchazinq/public/minitrees/nominal/Stop/TTTo2L2Nu_ext1.root",
			 "/afs/cern.ch/user/b/bchazinq/public/minitrees/nominal/Stop/T2tt_mStop-150to1200.root"};

  for (int dt = 0; dt<2; dt++) {

    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    
    int MaxEntries = -1; int GoodEntries = 0;

    for (Int_t i = 0; i<nentries && (MaxEntries==-1 || GoodEntries<MaxEntries) ; i++) {
      
      MiniTree->GetEntry(i);
      
      // Apply ttbar selection
      if (njet<2) continue;
      if (nbjet30csvv2m<1) continue;

      if (dt==1 && fabs(Mstop-Mlsp-172.44)>10.) continue;
      //if (dt==1 && Mstop-Mlsp-172.44>-10.) continue;

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
	
      double Rpt = pT_CM.Mag()/(pT_CM.Mag() + SHATR/4.);
      double invGamma = 1./gamma_R;
      double Mdr = SHATR/gamma_R;
      double DeltaPhiRll = dphi_LL_vBETA_T;

      //cout << dt << " " << MDELTAR << " " << Mdr << " " << DeltaPhiRll << " " << invGamma << endl;

      GoodEntries++;

    }

  }

}
