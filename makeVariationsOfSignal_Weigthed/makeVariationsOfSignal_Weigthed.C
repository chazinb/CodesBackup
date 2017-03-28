#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include <fstream>
#include <iostream>



void update(TRandom3 *b, TH1F *h) {

  cout << "----Starting iterations for histogram: " << h->GetName() << endl;
  for( int i = 0; i < h->GetNbinsX()+1; i++) {
    
    float lambda = h->GetBinContent(i);
    float lambdaerror = h->GetBinError(i);
    cout << lambda << " " << lambdaerror << endl;
    float N, factor;
    if(lambda > 0) {
      N = (lambda/lambdaerror) * (lambda/lambdaerror);
      factor = lambda/N;
    } else { 
      N = 0;
      factor = 0;
    }
    float newlambda = b->Poisson(N);
    h->SetBinContent(i, factor*newlambda);
    h->SetBinError(i, factor*sqrt(newlambda));
   cout << "Bin " << i << " had " << N << " events, with lambda " << lambda << " updating to " << factor*newlambda << endl;
  }
 cout << "--------------------------------------" << endl;
}


void run(std::string name, unsigned int number) {
 
  TString IterationName = "MT2ll_metPfType1_";
  TString massP = "350_225";
  int nchannel = 3; 
  int nsample = 6; 

  TRandom3 *b = new TRandom3(0);

  TFile *f = new TFile(name.c_str());
  TH1F	*T2tt_350_225_ee = (TH1F *)f->Get("T2tt_350_225_ee");
  TH1F	*T2tt_350_225_mm = (TH1F *)f->Get("T2tt_350_225_mm");
  TH1F	*T2tt_350_225_em = (TH1F *)f->Get("T2tt_350_225_em");
  TH1F	*TTbar_ee	 = (TH1F *)f->Get("TTbar_ee");
  TH1F	*data_obs_ee	 = (TH1F *)f->Get("data_obs_ee");
  TH1F	*TTbar_mm	 = (TH1F *)f->Get("TTbar_mm");
  TH1F	*data_obs_mm	 = (TH1F *)f->Get("data_obs_mm");
  TH1F	*TTbar_em	 = (TH1F *)f->Get("TTbar_em");
  TH1F	*data_obs_em	 = (TH1F *)f->Get("data_obs_em");
  TH1F	*ST_ee	         = (TH1F *)f->Get("ST_ee");
  TH1F	*ST_mm	         = (TH1F *)f->Get("ST_mm");
  TH1F	*ST_em	         = (TH1F *)f->Get("ST_em");
  TH1F	*DYJets_ee	         = (TH1F *)f->Get("DYJets_ee");
  TH1F	*DYJets_mm	         = (TH1F *)f->Get("DYJets_mm");
  TH1F	*DYJets_em	         = (TH1F *)f->Get("DYJets_em");
  TH1F	*WW_ee	         = (TH1F *)f->Get("WW_ee");
  TH1F	*WW_mm	         = (TH1F *)f->Get("WW_mm");
  TH1F	*WW_em	         = (TH1F *)f->Get("WW_em");
  TH1F	*TTZ_ee	         = (TH1F *)f->Get("TTZ_ee");
  TH1F	*TTZ_mm	         = (TH1F *)f->Get("TTZ_mm");
  TH1F	*TTZ_em	         = (TH1F *)f->Get("TTZ_em");


  for(unsigned int i = 0; i < number; i++) {
     TString  newName; 
     newName = IterationName; newName += i; 
     
    //char newName[50];
    //sprintf(newName, "MT2ll_metPfType1_%d", i);
    cout << "Generating file " << newName << endl;
    TFile *output = new TFile( newName + ".root", "RECREATE");
    output->cd();
    TH1F *S_ee = (TH1F *)T2tt_350_225_ee->Clone();
    TH1F *S_mm = (TH1F *)T2tt_350_225_mm->Clone();
    TH1F *S_em = (TH1F *)T2tt_350_225_em->Clone();
    update(b, S_ee);
    update(b, S_mm);
    update(b, S_em);
    S_ee->Write();
    S_mm->Write();
    S_em->Write();
    TTbar_ee ->Write();
    data_obs_ee->Write();
    TTbar_mm->Write();
    data_obs_mm->Write();
    TTbar_em->Write();
    data_obs_em->Write();
    ST_ee->Write();
    ST_mm->Write();
    ST_em->Write();
    DYJets_ee->Write();
    DYJets_mm->Write();
    DYJets_em->Write();
    WW_ee->Write();
    WW_mm->Write();
    WW_em->Write();
    TTZ_ee->Write();
    TTZ_mm->Write();
    TTZ_em->Write();

   // Write the datacard
   // --------------------------------------------------------------------------------------------------------

   std::ofstream inFile("T2tt_" + massP + "_" + newName +"_makeVariationsOfSignal_Weigthed.txt",std::ios::out);   
 
   inFile << "\n";
   inFile << "   max  "<<"     " << nchannel   << endl;
   inFile << "   jmax "<<"     " << nsample -1 << endl;
   inFile << "   kmax "<<"     " <<  "*"       << endl;
   inFile << "-----------------------------------------------------------------------------------" << endl;
   inFile << "  shapes * * " + newName  + ".root  $PROCESS_$CHANNEL" << endl;
   inFile << "-----------------------------------------------------------------------------------" << endl;
   inFile << "   bin  "<<"      "<< "     ee          "<<"     mm             "<<"      em             "<<endl;
   inFile << "observation           ";
   float data_ee = data_obs_ee -> Integral(0, data_obs_ee -> GetNbinsX()+1);
      inFile <<"" <<  data_ee << "         ";
   float data_mm = data_obs_mm ->  Integral(0, data_obs_mm -> GetNbinsX()+1);
      inFile <<"" <<  data_mm << "         ";
   float data_em = data_obs_em -> Integral(0, data_obs_em -> GetNbinsX()+1);
      inFile <<"" <<  data_em << "         ";
  
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
   inFile << "   " << "T2tt_350_225           T2tt_350_225              T2tt_350_225"  << "   " << "TTbar      TTbar     TTbar " << "   " <<  "   ST      ST     ST    " << "    " << " DYJets      DYJets      DYJets " << "    "  <<  "WW         WW       WW   "<< "     " << " TTZ    TTZ    TTZ";
   inFile << "" <<endl;
   inFile << " rate        ";

   float integral = S_ee -> Integral(0, S_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral = S_mm -> Integral(0, S_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral = S_em -> Integral(0, S_em -> GetNbinsX()+1);
   inFile << integral<<"      ";

    integral =TTbar_ee -> Integral(0, TTbar_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =TTbar_mm -> Integral(0, TTbar_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =TTbar_em -> Integral(0, TTbar_em -> GetNbinsX()+1);
   inFile << integral<<"      ";
   
   integral =ST_ee -> Integral(0, ST_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =ST_mm -> Integral(0, ST_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =ST_em -> Integral(0, ST_em -> GetNbinsX()+1);
   inFile << integral<<"      ";

    integral =DYJets_ee -> Integral(0, DYJets_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =DYJets_mm -> Integral(0, DYJets_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =DYJets_em -> Integral(0, DYJets_em -> GetNbinsX()+1);
   inFile << integral<<"      ";
  

   integral = WW_ee -> Integral(0, WW_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =WW_mm -> Integral(0, WW_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =WW_em -> Integral(0, WW_em -> GetNbinsX()+1);
   inFile << integral<<"      ";

   integral =TTZ_ee -> Integral(0, TTZ_ee -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =TTZ_mm -> Integral(0, TTZ_mm -> GetNbinsX()+1);
   inFile << integral<<"      ";
    integral =TTZ_em -> Integral(0, TTZ_em -> GetNbinsX()+1);
   inFile << integral<<"      ";
   inFile << "" <<endl;
  
   inFile.close();
//---------------------------------------------------------------------------------------------------------------------------------
    output->Close();

 }

  f->Close();

}
