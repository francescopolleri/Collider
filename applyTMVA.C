#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

// USAGE
// .L applyMVA.C
// applyMVA("BDT")
// applyMVA("Fisher")

void applyTMVA(TString method){

   //---------------------------------------------------------------
   TMVA::Tools::Instance();
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   Float_t official_mv2;
   Float_t ip2d_pb;
   Float_t ip2d_pu;
   Float_t ip3d_pu;
   Float_t ip3d_pb;
   Float_t sv1_efc;
   Float_t sv1_m;
   Float_t sv1_l2d;
   Float_t sv1_l3d;
   Float_t jf_sig3d;
   Float_t jf_deltaR;
   Float_t jf_nvtx;
   Float_t jf_nvtx1t;
   Int_t label;

   Int_t num=0;
   Int_t den=0;

   // --> Configure reader (same variables, in the same order as in the training) <--
  
   reader->AddVariable("ip2d_pb", &ip2d_pb);
   reader->AddVariable("ip2d_pu", &ip2d_pu);
   
   reader->AddVariable("sv1_efc", &sv1_efc);
   reader->AddVariable("sv1_m", &sv1_m);
   
   reader->AddVariable("jf_sig3d", &jf_sig3d);
   reader->AddVariable("jf_deltaR",&jf_deltaR);
   reader->AddVariable("jf_nvtx",&jf_nvtx);

   /*   reader->AddVariable("ip3d_pb", &ip3d_pb);
   reader->AddVariable("ip3d_pu", &ip3d_pu);
   reader->AddVariable("sv1_l2d", &sv1_l2d);
   reader->AddVariable("sv1_l3d", &sv1_l3d);
   reader->AddVariable("jf_nvtx1t",&jf_nvtx1t);*/
   
   // Loading config. files
   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassificationCategory";
   TString methodName = method + TString(" method");
   TString weightfile = dir + prefix + TString("_") + method + TString(".weights.xml");
   reader->BookMVA( methodName, weightfile ); 

   // Reading data
   TFile*  input   = new TFile("FlavTag.root");
   TTree*  theTree = (TTree*) input->Get("FlavTag");
   
   // --> Access to Tree variables (same name used in the reader) <--

   theTree->SetBranchAddress("mv2", &official_mv2);

   theTree->SetBranchAddress("ip2d_pb", &ip2d_pb);
   theTree->SetBranchAddress("ip2d_pu", &ip2d_pu);

   theTree->SetBranchAddress ("sv1_efc", &sv1_efc);
   theTree->SetBranchAddress ("sv1_m", &sv1_m);


   theTree->SetBranchAddress  ("jf_sig3d", &jf_sig3d);
   theTree->SetBranchAddress  ("jf_deltaR",&jf_deltaR);
   theTree->SetBranchAddress  ("jf_nvtx",&jf_nvtx);

   /* theTree->SetBranchAddress("ip3d_pb", &ip3d_pb);
   theTree->SetBranchAddress("ip3d_pu", &ip3d_pu);
   theTree->SetBranchAddress  ("sv1_l2d", &sv1_l2d);
   theTree->SetBranchAddress  ("sv1_l2d", &sv1_l2d);
   theTree->SetBranchAddress  ("jf_nvtx1t",&jf_nvtx1t);*/

   theTree->SetBranchAddress("label", &label);

   
   

   // --> Define histograms to check mva performance <--
   // s == signal ; b == background
   TH1F *hMVA = new TH1F( "MVA_"+method, "MVA_"+method, 100, -1, 1);
   TH1F *hmv2 = new TH1F( "official_mv2", "official_mv2", 100, -1, 1);      
   TH1F *hs   = new TH1F( "hs", "hs", 100, -1.5, 2.5);
   TH1F *hb   = new TH1F( "hb", "hb", 100, -1.5, 2.5);

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
     
      theTree->GetEntry(ievt);

      // mva_val -> value of the discriminant variable for this event      
      Float_t mva_val = reader->EvaluateMVA(methodName);

      // distribution of the MVA score
      hMVA->Fill(mva_val);
      hmv2->Fill(official_mv2);
      
      // How does mva_val compare with the official MV2 classifier? 
      // The mv2 score could be saved in the root file below for further comparisons
      
      // same distribution with truth label
      if (label == 5)
	hs->Fill(mva_val);
      else if (label == 0)
	hb->Fill(mva_val);
   }

   //hMVA->Draw();   
   // --> Compare results
   hb->Draw(); hs->Draw("SAME"); hs->SetLineColor(2);
   
   TFile *output = TFile::Open("MVA_plots.root", "RECREATE");
   output->cd();
   
   hs->Write();
   hb->Write();
   hmv2->Write();
   
   output->Close();

} 
