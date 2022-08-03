#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <list>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <typeinfo>
#include <TGraph.h>
#include <TCanvas.h>
#include <TNtuple.h> 
#include <TH1.h>
#include <TAttMarker.h>

void plotsaver(){
    // no muon dxy track dxy, not sure if trigger matched delta R delta Z no JPSI VERTEXPROB no jspi Lxy no cosalpha no jpsi pt, not sure if trigger matched
    TFile *f = TFile::Open("/eos/user/j/jodedra/Run3/NoTrigMatch/MergedFourStreamsNotrigMatchDoubleEle.root");
    TTree *events = (TTree*)f->Get("nano/tree");
    TCanvas *c1 = new TCanvas("c1"); 
    TH1D *myh = new TH1D("myh"," B_mass plot 2022 data HLTPhysics",40,4.5,6.0);
    myh->SetOption("E");
    myh->SetXTitle("Mass ee+k (GeV)");

    /*
    TCut hlt = "HLT_DoubleMu4_JpsiTrk_Displaced>0.99";

    TCut ptmu2 = "SkimBToKMuMu_fit_l2_pt>4";
    TCut ptmu1 = "SkimBToKMuMu_fit_l2_pt>4";
    TCut kpt = "SkimBToKMuMu_fit_k_pt>1.2";

    TCut dRmumu = "SkimBToKMuMu_l1l2_dr<1";
    TCut dzmumu = "SkimBToKMuMu_l1l2_dz<1";
    TCut mumumass = "SkimBToKMuMu_mll_fullfit<3.3&&SkimBToKMuMu_mll_fullfit>2.9";

    
    TCut ptbcand = "SkimBToKMuMu_fit_pt>6";
    //TCut dzlk = "SkimBToKMuMu_lk_dz<1";


    TCut displacementisg = "SkimBToKMuMu_l_xy_sig>3";
    TCut cos2d = "SkimBToKMuMu_fit_cos2D>0.90";
    TCut probbvertex = "SkimBToKMuMu_svprob>0.1";
    
   
    TCut final =hlt&&ptmu1&&ptmu2&&kpt&&dRmumu&&dzmumu&&mumumass&&ptbcand&&probbvertex&&cos2d&&displacementisg ;*/
    TCut jpsimassl = "JpsiKE_Jpsi_mass>2.9";
    TCut jpsimassh = "JpsiKE_Jpsi_mass<3.3";
    TCut jpsipdig = "TMath::Abs(JpsiKE_pi_pdg)==211";
    TCut vertexprob = "JpsiKE_B_vprob>0.01";
    TCut fls3d ="JpsiKE_B_fls3d>3.";
    TCut lower = "JpsiKE_B_mass>4.5";
    TCut higher ="JpsiKE_B_mass<6.0";
    TCut final = lower&&higher&&jpsimassl&&jpsimassh&&jpsipdig&&vertexprob&&fls3d;
    events->Draw("JpsiKE_B_mass>>myh",final);
    myh->SetTitle("JpsiKE_Jpsi_mass>2.9 , JpsiKE_Jpsi_mass<3.3 , TMath::Abs(JpsiKE_pi_pdg)==211 , JpsiKE_B_vprob>0.01 , JpsiKE_B_fls3d>3. , JpsiKE_B_mass>4.5 , JpsiKE_B_mass<6.0, INT LUMI average = 0.674 FB-1 ");
    c1->SaveAs("BMASS5GB.png");  
    c1->SaveAs("BMASS5GB.pdf");   
 
}
