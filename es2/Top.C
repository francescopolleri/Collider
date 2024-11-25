#define Top_cxx
#include "Top.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;

TLorentzVector bjetl;
TLorentzVector bjeth;
TLorentzVector ljet0;
TLorentzVector ljet1;
TLorentzVector mu;
TLorentzVector nu;

const float WBoson_width = 2.085;
const float TopQuark_width = 1.35;
const float jetReso_constraint = 0.2;  //risoluzione su p_T jet
const double mW = 80.379; // PDG
const int verb = 0; //debug variable

//function to be minimized
void chi2_f(int& npar, double *deriv, double& f, double *par, int flag){

    //Setting the particle 4-vectors entering in the chi2: the b-jet from hadronic top,
    // the other from the leptonic top, the light jets, the neutrino

    TLorentzVector bjeth_f(bjeth); // _f vuol dire "fit"
    TLorentzVector bjetl_f(bjetl);
    TLorentzVector ljet0_f(ljet0);
    TLorentzVector ljet1_f(ljet1);
    TLorentzVector nu_f(nu); //starting guess of the neutrino from the external TlorentzVector
    //Initial parameters, some are unknown some are known
    bjetl_f.SetPtEtaPhiM(par[0],bjetl.Eta(),bjetl.Phi(), 5);  //5 fa riferimento al fatto che è jet b
    bjeth_f.SetPtEtaPhiM(par[1],bjeth.Eta(),bjeth.Phi(), 5);
    ljet0_f.SetPtEtaPhiM(par[2],ljet0.Eta(),ljet0.Phi(), 0);  //0 fa riferimento a light quark
    ljet1_f.SetPtEtaPhiM(par[3],ljet1.Eta(),ljet1.Phi(), 0);

    
    //  MET during the fit
    double nuPx = -(bjetl_f.Px() + bjeth_f.Px() +  ljet0_f.Px() + ljet1_f.Px() + mu.Px());
    double nuPy = -(bjetl_f.Py() + bjeth_f.Py() +  ljet0_f.Py() + ljet1_f.Py() + mu.Py());
    nu_f.SetPxPyPzE( nuPx, nuPy, par[4], sqrt(nuPx*nuPx+nuPy*nuPy+par[4]*par[4]) );

    //4 vector of top quark with leptonic decay and of corresponding W->lv
    TLorentzVector wl  = nu_f + mu;
    TLorentzVector tl  = wl + bjetl_f; 
    //need to add the rest: wh and th
    TLorentzVector wh = ljet0_f + ljet1_f;
    TLorentzVector th = wh +bjeth_f;
    
    f  = 0;
    //Example of a constraint on the resolution of a b-jets and a jet from the W decay
    f += pow((bjetl_f.Pt() - bjetl.Pt())/(jetReso_constraint*bjetl_f.Pt()),2);
    f += pow((ljet0_f.Pt() - ljet0.Pt())/(jetReso_constraint*ljet0_f.Pt()),2);
    f += pow((ljet1_f.Pt() - ljet1.Pt())/(jetReso_constraint*ljet1_f.Pt()),2);
    f += pow((bjeth_f.Pt() - bjeth.Pt())/(jetReso_constraint*bjeth_f.Pt()),2);

    //Example of constraint on the top and W masses using their width and the Gaussian approximation
    f += pow((tl.M() - par[5])/TopQuark_width,2);
    f += pow((wl.M() - mW)/ WBoson_width,2);
    //Add corresponding implementation for hadronic top decay...
    f += pow((th.M() - par[5])/TopQuark_width,2);
    f += pow((wh.M() - mW)/ WBoson_width,2);


}


void Top::Loop(string outFileStr="", bool isMC = true)
{
    //   In a ROOT session, you can do:
    //      root> .L Top.C+g // Compile the macro (with debug symbols)
    //      root> TFile f("myTreeFile.root")
    //      root> Top t((TTree*)f.Get("Top"));
    //      root> t.Loop("outFile.root");  // Loop on all entries and save the mTop output hist
    //
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   //Plots that can be used for testing the assignation
   TH1D* h_tophAssign0 = new TH1D("h_tophAssign0","",150, 0,300);
   TH1D* h_tophAssign1 = new TH1D("h_tophAssign1","",150, 0,300);
   TH1D* h_diff = new TH1D("h_diff","",100, -150,150);


   //Example of plots that can be used for testing the jet pre-fit and post fit against the truth
   TH1D* h_blPtReso = new TH1D("h_blPtReso","",100,-1,1);
   TH1D* h_blPtPosFitReso = new TH1D("h_blPtPosFitReso","",100,-1,1);
   TH1D* h_bhPtReso = new TH1D("h_bhPtReso","",100,-1,1);
   TH1D* h_bhPtPosFitReso = new TH1D("h_bhPtPosFitReso","",100,-1,1);
   TH1D* h_l0PtReso = new TH1D("h_l0PtReso","",100,-1,1);
   TH1D* h_l0PtPosFitReso = new TH1D("h_l0PtPosFitReso","",100,-1,1);
   TH1D* h_l1PtReso = new TH1D("h_l1PtReso","",100,-1,1);
   TH1D* h_l1PtPosFitReso = new TH1D("h_l1PtPosFitReso","",100,-1,1);
   

   //top masses, prefit and then post-fit:
   // NB: The x-range of mTop is not too large to facilitate the comparison of several templates
   // using a chi2 test-statistics, which is appropriate if bins are well populated 
   TH1D* h_tl_preFit  = new TH1D("h_tl_preFit","",45,145,205);
   TH1D* h_th_preFit  = new TH1D("h_th_preFit","",45,145,205);
   
   //Example of plots for saving the final output values, i.e. the post fit chi^2 and the MTop
   TH1D* h_tl_posFit  = new TH1D("h_tl_posFit","",45,145,205);
   TH1D* h_th_posFit  = new TH1D("h_th_posFit","",45,145,205);
   //...
   TH1D* h_mt   = new TH1D("h_mt","",45,145,205);
   TH1D* h_chi2 = new TH1D("h_chi2","",150,0,1500);

   
   //Actual loop over the entries
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      //   if( jentry > 1000) break;//Uncomment for debugging a few events at the time
      
      //Retrieval of b-jets, jets and leptons from the ttree
      vector<TLorentzVector> bjet;
      vector<TLorentzVector> ljet;
      vector<int> jettruth_decay;
      
      for (int i=0;i<njet;i++){
	        TLorentzVector tmp;

	        if (jetlabel[i]==5){//This is the b-tagging label 	    
	          tmp.SetPtEtaPhiM(jetpt[i],jeteta[i],jetphi[i],5);
	          bjet.push_back(tmp);
	          jettruth_decay.push_back(jettruth[i]);//is ==1 if this is a "truth" leptonic decay
	        } else {//not tagged jets (NB: always second 2 jets in the ntuples)
	          tmp.SetPtEtaPhiM(jetpt[i],jeteta[i],jetphi[i],0);
	          ljet.push_back(tmp);
	        }
      }

      ljet0=ljet[0];
      ljet1=ljet[1];
      mu.SetPtEtaPhiM(mupt,mueta,muphi,0.1);
      nu.SetPtEtaPhiM(metpt,0,metphi,0);

      //Now test the two possible assignation of the b-jet and decide for the best one
      TLorentzVector toplAssign0 = bjet[0]+mu+nu;
      TLorentzVector tophAssign0 = bjet[1]+ljet0+ljet1;
      TLorentzVector toplAssign1 = bjet[1]+mu+nu;
      TLorentzVector tophAssign1 = bjet[0]+ljet0+ljet1;
      //Examples of M from had decay with b assigned in the first case or second case
      h_tophAssign0->Fill(tophAssign0.M());	    
      h_tophAssign1->Fill(tophAssign1.M());

      //(tophAssign0.M()-tophAssign1.M()): Diff. between well assign. and misassign. b to had. decay
      //(toplAssign0.M()-toplAssign1.M()): Diff. between well assign. and misassing. b to lep. decay	  
      //NB: an initial test was using only the lepton assignment (a bit waker)
      double hadLep_assign_diff = (tophAssign0.M()-tophAssign1.M())+ (toplAssign0.M()-toplAssign1.M());
      int bjetl_assign_i, bjeth_assign_i;
      h_diff->Fill( hadLep_assign_diff );
      if(  hadLep_assign_diff > 0 ){//use type 0 assignation
	      bjeth=bjet[1];
	      bjetl=bjet[0];
	      bjetl_assign_i=0;
	      bjeth_assign_i=1;
	    } else {//use type 1 assignation (but I know from truth that this is wrong)
	      bjeth=bjet[0];
	      bjetl=bjet[1];
	      bjetl_assign_i=1;
	      bjeth_assign_i=0;
	      //NB: because of the simplified procedure to build the ntuple, this answer is "uncorrect", i.e.
	      // the "assignation 0" is what was done in "truth". However it is instructive to understand "assignation"
	      // on the basis of h_diff because the few erroneous cases are due to tails in the b-jet resolution.
	      // Uncomment the following "continue" if you want to test only cases where all M_top are correctly
	      // assigned as you may see small shifts in the final mTop value of the kinematic fit
	      //  continue;
      }
      
      
      ////////////////////////////////////////////////////////
      ///// NOW STEP 2 WITH THE DEFINITION OF MINIMIZED FCN
      ////////////////////////////////////////////////////////

      // MINUIT PARAMETER DEFINITION
      // A different number of parameters is possible
      // 6 parameters (4 jet pt, 1 pZnu, 1 mTop)
      const int nPar = 6;
      //Initialization of Miniuit fitter
      TMinuit fitter(nPar);
      // Which function to minimise (FCN for historical reasons)
      fitter.SetFCN(chi2_f);
      // Printout verbosity (-1 ==> silence; 0 ==> standard)
      fitter.SetPrintLevel(-1);
      // Shift from the minimum for estimating the errors on the parameters
      fitter.SetErrorDef(1.);


      // actual obecjts which are the parameters to the fit (NB init pZ =0)
      double par[nPar] = { bjetl.Pt(),bjeth.Pt(),ljet0.Pt(),ljet1.Pt(), 0, 172};
      
      // incremental fixed steps for the minimisation
      double step[nPar]={0.01,0.01,0.01,0.01,0.01,0.01};
      // 0 stands for no bounds on the minima
      double min[nPar] ={0.00,0.00,0.00,0.00,0.00,0};
      // 0 stands for no bounds on the maxima
      double max[nPar] ={0.00,0.00,0.00,0.00,0.00,0};
      // strings for the 6 parameters
      string parName[nPar]={"jbl","jbh","jl1","jl2","pnz","mtop"};
      
      for (int i=0; i<nPar; i++){
	      fitter.DefineParameter(i, parName[i].c_str(), par[i],step[i],min[i],max[i]); 
      }

      //Run MIGRAD and minimize the function and retrieve the parameters
      fitter.Migrad();

      // retrieve the output of Minuit parameters
      double chi2_val;
      double pbl, epbl, pbh, epbh, pl0, epl0, pl1, epl1, xt, ext;
      double pz, epz;      

      //extract post fit parameters
      //pT of the b-jet assigned to leptonic top decay
      fitter.GetParameter(0,pbl,epbl);
      fitter.GetParameter(0,pbh,epbh);
      fitter.GetParameter(0,pl0,epl0);
      fitter.GetParameter(0,pl1,epl1);
      
      //fitted top mass
      fitter.GetParameter(5,xt,ext);

      //Evaluate the chi2 at the minimum
      int status = fitter.Eval(nPar, nullptr, chi2_val, par, 3);
      if(verb) cout<<"chi2_val = "<<chi2_val<<endl;
      h_chi2->Fill(chi2_val);

      //Initial (pre-fit) 4-vector for hadronic and leptonic top-quarks
      TLorentzVector mtl_preFit(bjetl + mu + nu); //here nu is the pre-fit MET, no pZ!
      TLorentzVector mth_preFit(bjeth + ljet0 + ljet1);//hadronic decay
      h_th_preFit->Fill(mth_preFit.M());
      h_tl_preFit->Fill(mtl_preFit.M());

      //Example of  the b-jet pT resolution based on the truth
      h_blPtReso->Fill( (bjetl.Pt() - jettruthpt[bjetl_assign_i])/jettruthpt[bjetl_assign_i] );
      h_bhPtReso->Fill( (bjeth.Pt() - jettruthpt[bjeth_assign_i])/jettruthpt[bjetl_assign_i] );


      //Example of the b-jet pT resolution after the fit
      TLorentzVector  posFit_bjetl(bjetl);
      posFit_bjetl.SetPtEtaPhiM(pbl, bjetl.Eta(), bjetl.Phi(), 5);
      h_blPtPosFitReso->Fill( (posFit_bjetl.Pt() - jettruthpt[bjetl_assign_i])/jettruthpt[bjetl_assign_i] );
      TLorentzVector  posFit_bjeth(bjeth);
      posFit_bjeth.SetPtEtaPhiM(pbh, bjeth.Eta(), bjeth.Phi(), 5);
      h_bhPtPosFitReso->Fill( (posFit_bjeth.Pt() - jettruthpt[bjeth_assign_i])/jettruthpt[bjeth_assign_i] );

      TLorentzVector posFit_ljet0(ljet0);
      posFit_bjetl.SetPtEtaPhiM(pl0, ljet0.Eta(), ljet0.Phi(), 0);
      TLorentzVector posFit_ljet1(ljet1);
      posFit_bjetl.SetPtEtaPhiM(pl1, ljet1.Eta(), ljet1.Phi(), 0);

      TLorentzVector  posFit_nu(nu);//NB THIS IS STILL THE PREFIT ONE!!

      //Example of the the reconstructed m-top after the fit
      TLorentzVector mtl_posFit(posFit_bjetl + mu + posFit_nu);
      h_tl_posFit->Fill(mtl_posFit.M());
      TLorentzVector mth_posFit(posFit_bjeth + posFit_ljet0 + posFit_ljet1);
      h_th_posFit->Fill(mth_posFit.M());

      h_mt->Fill(xt);
   }

   TCanvas *c_diff = new TCanvas;
   h_diff->Draw();
   
   TCanvas *c_reso = new TCanvas("c_reso","c_reso",900, 600);
   c_reso->Divide(2,1);
   c_reso->cd(1);
   h_blPtReso->Draw();
   h_blPtPosFitReso->SetLineColor(1);
   h_blPtPosFitReso->Draw("SAME");
   
   c_reso->cd(2);
   h_bhPtReso->Draw();
   h_bhPtPosFitReso->SetLineColor(1);
   h_bhPtPosFitReso->Draw("SAME");

   //c_reso->cd(5);   
   //h_blPtPosFitReso->Draw();
   //...
   //...

   TCanvas *topM_c = new TCanvas("topM_c","topM_c",600, 600);
   h_mt->Draw();

   if(outFileStr != ""){
       TFile fOut(outFileStr.c_str(), "recreate");
       h_mt->Write();
       h_chi2->Write();
       fOut.Close();
   }

}