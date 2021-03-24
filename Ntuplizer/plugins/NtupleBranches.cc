#include "../interface/NtupleBranches.h"

//===================================================================================================================        
NtupleBranches::NtupleBranches( std::map< std::string, bool >& runFlags, TTree* tree )
  : tree_( tree )
{
  branch( runFlags );
}

//===================================================================================================================
NtupleBranches::~NtupleBranches( void )
{
}

//===================================================================================================================      
void NtupleBranches::branch( std::map< std::string, bool >& runFlags ){


  if ( runFlags["runOnMC"] ){
    if ( runFlags["doGenParticles"] ){
      /** genParticles */
      tree_->Branch( "genParticle_N"	     , &genParticle_N	       );
      tree_->Branch( "genParticle_pt"	     , &genParticle_pt	       ); 
      //      tree_->Branch( "genParticle_px"	     , &genParticle_px	       ); 
      //      tree_->Branch( "genParticle_py"	     , &genParticle_py	       ); 
      //      tree_->Branch( "genParticle_pz"	     , &genParticle_pz	       ); 
      //      tree_->Branch( "genParticle_e" 	     , &genParticle_e	       ); 
      tree_->Branch( "genParticle_eta"	     , &genParticle_eta        ); 
      tree_->Branch( "genParticle_phi"	     , &genParticle_phi        ); 
      tree_->Branch( "genParticle_mass"	     , &genParticle_mass       ); 
      tree_->Branch( "genParticle_pdgId"     , &genParticle_pdgId      );
      tree_->Branch( "genParticle_status"    , &genParticle_status     );
      tree_->Branch( "genParticle_isPrompt"  , &genParticle_isPrompt   );
      tree_->Branch( "genParticle_isDirectPromptTauDecayProduct"  , &genParticle_isDirectPromptTauDecayProduct);
      tree_->Branch( "genParticle_isDirectHardProcessTauDecayProductFinalState"  , &genParticle_isDirectHardProcessTauDecayProductFinalState);
      tree_->Branch( "genParticle_fromHardProcessFinalState"  , &genParticle_fromHardProcessFinalState   );
      tree_->Branch( "genParticle_mother"    , &genParticle_mother     );
      tree_->Branch( "genParticle_mother_pt" , &genParticle_mother_pt  );
      tree_->Branch( "genParticle_nMoth"     , &genParticle_nMoth      );
      tree_->Branch( "genParticle_nDau"	     , &genParticle_nDau       ); 
      tree_->Branch( "genParticle_dau"	     , &genParticle_dau        );

      //      tree_->Branch( "genParticle_pmother"	     , &genParticle_pmother        );
      tree_->Branch( "genParticle_pdgs"	     , &genParticle_pdgs        );
      tree_->Branch( "genParticle_layers"	     , &genParticle_layers        );
      tree_->Branch( "genParticle_ppt"	     , &genParticle_ppt        );
      tree_->Branch( "genParticle_peta"	     , &genParticle_peta        );
      tree_->Branch( "genParticle_pphi"	     , &genParticle_pphi        );
      // ////Adding weight based on B decay chain for B generic bkg sample
      // tree_->Branch( "genWeightBkgB"             , &genWeightBkgB       );
	
    } //doGenParticles
      
    if ( runFlags["doGenEvent"] ){
      /** generator info */
      tree_->Branch( "lheV_pt"	             , &lheV_pt                ); 
      tree_->Branch( "lheHT"	             , &lheHT                  ); 
      tree_->Branch( "lheNj"	             , &lheNj                  );
      tree_->Branch( "lheNb"	             , &lheNb                  );
      tree_->Branch( "lheNl"	             , &lheNl                  );
      tree_->Branch( "lheV_mass"           , &lheV_mass              ); 
      tree_->Branch( "genWeight"	         , &genWeight              );
      tree_->Branch( "genFacWeightUp"	     , &genFacWeightUp         );
      tree_->Branch( "genFacWeightDown"	   , &genFacWeightDown       );
      tree_->Branch( "genRenWeightUp"	     , &genRenWeightUp         );
      tree_->Branch( "genRenWeightDown"	   , &genRenWeightDown       );
      tree_->Branch( "genFacRenWeightUp"	 , &genFacRenWeightUp      );
      tree_->Branch( "genFacRenWeightDown" , &genFacRenWeightDown    );
      tree_->Branch( "qScale"	             , &qScale                 );
      tree_->Branch( "PDF_rms"	           , &PDF_rms                );
      tree_->Branch( "PDF_x"	             , &PDF_x                  );
      tree_->Branch( "PDF_xPDF"	           , &PDF_xPDF               );
      tree_->Branch( "PDF_id"	             , &PDF_id                 );

    } //doGenEvent
  } //runOnMC
  
  
  if (runFlags["doTriggerDecisions"]) {
    /** HLT trigger decisions */
    tree_->Branch("HLT_isFired", &HLT_isFired );
  }

  
  if (runFlags["doTriggerObjects"]) {
    /** HLT trigger objects */
    tree_->Branch("triggerObject_pt"		, &triggerObject_pt		);
    tree_->Branch("triggerObject_eta"		, &triggerObject_eta		);
    tree_->Branch("triggerObject_phi"		, &triggerObject_phi	        );
    tree_->Branch("triggerObject_mass"		, &triggerObject_mass	        );
    tree_->Branch("triggerObject_lastname"	, &triggerObject_lastname	);
    tree_->Branch("triggerObject_filterLabels"	, &triggerObject_filterLabels	);
    tree_->Branch("triggerObject_firedTrigger"	, &triggerObject_firedTrigger	);
    tree_->Branch("triggerObject_filterIDs"	, &triggerObject_filterIDs	);

  } //doTriggerObjects
  
  if (runFlags["doHltFilters"]) {
    /** HLT filter decisions */
    tree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
    tree_->Branch("passFilter_HBHELoose"            ,&passFilter_HBHELoose_	          ,"passFilter_HBHELoose_/O");
    tree_->Branch("passFilter_HBHETight"            ,&passFilter_HBHETight_	          ,"passFilter_HBHETight_/O");
    tree_->Branch("passFilter_HBHEIso"              ,&passFilter_HBHEIso_	            ,"passFilter_HBHEIso_/O");
    tree_->Branch("passFilter_CSCHalo"              ,&passFilter_CSCHalo_             ,"passFilter_CSCHalo_/O");
    tree_->Branch("passFilter_CSCTightHalo2015"     ,&passFilter_CSCTightHalo2015_    ,"passFilter_CSCTightHalo2015_/O");
    tree_->Branch("passFilter_HCALlaser"            ,&passFilter_HCALlaser_           ,"passFilter_HCALlaser_/O");
    tree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
    tree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
    tree_->Branch("passFilter_TrkFailure"           ,&passFilter_TrkFailure_          ,"passFilter_TrkFailure_/O");
    tree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
    tree_->Branch("passFilter_ECALlaser"            ,&passFilter_ECALlaser_           ,"passFilter_ECALlaser_/O");
    tree_->Branch("passFilter_TrkPOG"               ,&passFilter_TrkPOG_              ,"passFilter_TrkPOG_/O");
    tree_->Branch("passFilter_TrkPOG_manystrip"     ,&passFilter_TrkPOG_manystrip_    ,"passFilter_TrkPOG_manystrip_/O");
    tree_->Branch("passFilter_TrkPOG_toomanystrip"  ,&passFilter_TrkPOG_toomanystrip_ ,"passFilter_TrkPOG_toomanystrip_/O");
    tree_->Branch("passFilter_TrkPOG_logError"      ,&passFilter_TrkPOG_logError_     ,"passFilter_TrkPOG_logError_/O");
    tree_->Branch("passFilter_METFilters"           ,&passFilter_METFilters_          ,"passFilter_METFilters_/O");
    
    //NEW FOR ICHEP
    tree_->Branch("passFilter_CSCTightHaloTrkMuUnvetoFilter", &passFilter_CSCTightHaloTrkMuUnvetoFilter_   ,"passFilter_CSCTightHaloTrkMuUnvetoFilter_/O");
    tree_->Branch("passFilter_globalTightHalo2016"          , &passFilter_globalTightHalo2016_             ,"passFilter_globalTightHalo2016_/O");
    tree_->Branch("passFilter_globalSuperTightHalo2016"          , &passFilter_globalSuperTightHalo2016_             ,"passFilter_globalSuperTightHalo2016_/O");
    tree_->Branch("passFilter_HcalStripHalo"                , &passFilter_HcalStripHalo_                   ,"passFilter_HcalStripHalo_/O");
    tree_->Branch("passFilter_chargedHadronTrackResolution" , &passFilter_chargedHadronTrackResolution_    ,"passFilter_chargedHadronTrackResolution_/O");
    tree_->Branch("passFilter_muonBadTrack"                 , &passFilter_muonBadTrack_                    ,"passFilter_muonBadTrack_/O");
    tree_->Branch("flag_badMuons"                 , &flag_badMuons_                    ,"flag_badMuons_/O");
    tree_->Branch("flag_duplicateMuons"                 , &flag_duplicateMuons_                    ,"flag_duplicateMuons_/O");
    tree_->Branch("flag_nobadMuons"                 , &flag_nobadMuons_                    ,"flag_nobadMuons_/O");
    tree_->Branch("passFilter_ecalBadCalib_"    ,&passFilter_ecalBadCalib_, "passFilter_ecalBadCalib_/O");
  } //do HltFilters

  if (runFlags["doMissingEt"]) {
    /** MET */
    tree_->Branch( "rho", &rho );
    tree_->Branch("METraw_et"		        , &METraw_et	     );
    tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 
    tree_->Branch("METraw_sumEt"		, &METraw_sumEt	     );   
    tree_->Branch("MET_corrPx"		        , &MET_corrPx	     ); 
    tree_->Branch("MET_corrPy"		        , &MET_corrPy	     );   
    tree_->Branch("MET_et"	                , &MET_et  	     ); 
    tree_->Branch("MET_phi"	                , &MET_phi           );
    tree_->Branch("MET_puppi_et"	        , &MET_puppi_et      ); 
    tree_->Branch("MET_puppi_phi"               , &MET_puppi_phi     );
    tree_->Branch("MET_sumEt"	                , &MET_sumEt 	     ); 
    tree_->Branch("MET_JetEnUp"	                , &MET_JetEnUp 	     ); 
    tree_->Branch("MET_JetEnDown"	                , &MET_JetEnDown 	     ); 
    tree_->Branch("MET_JetResUp"	                , &MET_JetResUp 	     ); 
    tree_->Branch("MET_JetResDown"	                , &MET_JetResDown 	     ); 
    tree_->Branch("MET_UnclusteredEnUp"	                , &MET_UnclusteredEnUp 	     ); 
    tree_->Branch("MET_UnclusteredEnDown"	                , &MET_UnclusteredEnDown 	     ); 
    
  } //doMissingEt

  if ( runFlags["doMVAMET"] ){
    /** MET SVift*/
    tree_->Branch("MET_Nmva"	                , &MET_Nmva 	     ); 
    tree_->Branch("MET_mva_et"	                , &MET_mva_et        ); 
    tree_->Branch("MET_mva_phi"                 , &MET_mva_phi       );
    tree_->Branch( "MET_mva_cov00"                                        , &MET_mva_cov00 );
    tree_->Branch( "MET_mva_cov10"                                        , &MET_mva_cov10 );
    tree_->Branch( "MET_mva_cov11"                                        , &MET_mva_cov11 );
    tree_->Branch( "MET_mva_recoil_pt"                                        , &MET_mva_recoil_pt );
    tree_->Branch( "MET_mva_recoil_eta"                                        , &MET_mva_recoil_eta );
    tree_->Branch( "MET_mva_recoil_phi"                                        , &MET_mva_recoil_phi );
    tree_->Branch( "MET_mva_recoil_pdgId"                                        , &MET_mva_recoil_pdgId );

  }

  
  /*------------- ------EVENT infos-----------------------------*/
  tree_->Branch("EVENT_event"	 , &EVENT_event     );
  tree_->Branch("EVENT_run"	 , &EVENT_run	    );
  tree_->Branch("EVENT_lumiBlock", &EVENT_lumiBlock );
  
  if (runFlags["runOnMC"]) {
    if (runFlags["doPileUp"]) {
      /*--------------------------PU infos--------------------------*/
      tree_->Branch("nPuVtxTrue", &nPuVtxTrue );	
      tree_->Branch("nPuVtx"    , &nPuVtx     );
      tree_->Branch("bX"	, &bX	      );
    } //doPileUp
  } //runOnMC
  
  if (runFlags["doVertices"]) {  
    /*--------------------------PV infos--------------------------*/
    tree_->Branch("PV_N"     , &PV_N      );
    tree_->Branch("PV_filter", &PV_filter );
    tree_->Branch("PV_chi2"  , &PV_chi2   );
    tree_->Branch("PV_ndof"  , &PV_ndof   );
    tree_->Branch("PV_rho"   , &PV_rho    );
    tree_->Branch("PV_z"     , &PV_z      );
    tree_->Branch("BeamSpot_x0" , &BeamSpot_x0     );
    tree_->Branch("BeamSpot_y0" , &BeamSpot_y0     );
    tree_->Branch("BeamSpot_z0" , &BeamSpot_z0     );


  }


  if (runFlags["doJpsiMu"]){
    tree_->Branch("IsJpsiMu" , &IsJpsiMu  );
    tree_->Branch("HLT_BPH_isFired", &HLT_BPH_isFired );

    tree_->Branch("JpsiMu_nCandidates", &JpsiMu_nCandidates );

    tree_->Branch("JpsiMu_mu1_pt", &JpsiMu_mu1_pt );
    tree_->Branch("JpsiMu_mu1_eta", &JpsiMu_mu1_eta );
    tree_->Branch("JpsiMu_mu1_phi", &JpsiMu_mu1_phi );
    tree_->Branch("JpsiMu_mu1_mass", &JpsiMu_mu1_mass );
    tree_->Branch("JpsiMu_mu1_unfit_pt", &JpsiMu_mu1_unfit_pt );
    tree_->Branch("JpsiMu_mu1_unfit_eta", &JpsiMu_mu1_unfit_eta );
    tree_->Branch("JpsiMu_mu1_unfit_phi", &JpsiMu_mu1_unfit_phi );
    tree_->Branch("JpsiMu_mu1_unfit_mass", &JpsiMu_mu1_unfit_mass );
    tree_->Branch("JpsiMu_mu1_q", &JpsiMu_mu1_q );
    tree_->Branch("JpsiMu_mu1_isLoose"  , &JpsiMu_mu1_isLoose   );
    tree_->Branch("JpsiMu_mu1_isTight"  , &JpsiMu_mu1_isTight   );
    tree_->Branch("JpsiMu_mu1_isPF"     , &JpsiMu_mu1_isPF      );
    tree_->Branch("JpsiMu_mu1_isGlobal" , &JpsiMu_mu1_isGlobal  );
    tree_->Branch("JpsiMu_mu1_isTracker", &JpsiMu_mu1_isTracker );
    tree_->Branch("JpsiMu_mu1_isSoft"   , &JpsiMu_mu1_isSoft    );
    tree_->Branch("JpsiMu_mu1_vx"   , &JpsiMu_mu1_vx    );
    tree_->Branch("JpsiMu_mu1_vy"   , &JpsiMu_mu1_vy    );
    tree_->Branch("JpsiMu_mu1_vz"   , &JpsiMu_mu1_vz    );
    tree_->Branch("JpsiMu_mu1_iso"   , &JpsiMu_mu1_iso    );
    tree_->Branch("JpsiMu_mu1_dbiso"   , &JpsiMu_mu1_dbiso    );

    tree_->Branch("JpsiMu_mu2_pt", &JpsiMu_mu2_pt );
    tree_->Branch("JpsiMu_mu2_eta", &JpsiMu_mu2_eta );
    tree_->Branch("JpsiMu_mu2_phi", &JpsiMu_mu2_phi );
    tree_->Branch("JpsiMu_mu2_mass", &JpsiMu_mu2_mass );
    tree_->Branch("JpsiMu_mu2_unfit_pt", &JpsiMu_mu2_unfit_pt );
    tree_->Branch("JpsiMu_mu2_unfit_eta", &JpsiMu_mu2_unfit_eta );
    tree_->Branch("JpsiMu_mu2_unfit_phi", &JpsiMu_mu2_unfit_phi );
    tree_->Branch("JpsiMu_mu2_unfit_mass", &JpsiMu_mu2_unfit_mass );
    tree_->Branch("JpsiMu_mu2_q", &JpsiMu_mu2_q );
    tree_->Branch("JpsiMu_mu2_isLoose"  , &JpsiMu_mu2_isLoose   );
    tree_->Branch("JpsiMu_mu2_isTight"  , &JpsiMu_mu2_isTight   );
    tree_->Branch("JpsiMu_mu2_isPF"     , &JpsiMu_mu2_isPF      );
    tree_->Branch("JpsiMu_mu2_isGlobal" , &JpsiMu_mu2_isGlobal  );
    tree_->Branch("JpsiMu_mu2_isTracker", &JpsiMu_mu2_isTracker );
    tree_->Branch("JpsiMu_mu2_isSoft"   , &JpsiMu_mu2_isSoft    );
    tree_->Branch("JpsiMu_mu2_vx"   , &JpsiMu_mu2_vx    );
    tree_->Branch("JpsiMu_mu2_vy"   , &JpsiMu_mu2_vy    );
    tree_->Branch("JpsiMu_mu2_vz"   , &JpsiMu_mu2_vz    );
    tree_->Branch("JpsiMu_mu2_iso"   , &JpsiMu_mu2_iso    );
    tree_->Branch("JpsiMu_mu2_dbiso"   , &JpsiMu_mu2_dbiso    );

    tree_->Branch("JpsiMu_mu3_pt", &JpsiMu_mu3_pt );
    tree_->Branch("JpsiMu_mu3_eta", &JpsiMu_mu3_eta );
    tree_->Branch("JpsiMu_mu3_phi", &JpsiMu_mu3_phi );
    tree_->Branch("JpsiMu_mu3_mass", &JpsiMu_mu3_mass );
    tree_->Branch("JpsiMu_mu3_unfit_pt", &JpsiMu_mu3_unfit_pt );
    tree_->Branch("JpsiMu_mu3_unfit_eta", &JpsiMu_mu3_unfit_eta );
    tree_->Branch("JpsiMu_mu3_unfit_phi", &JpsiMu_mu3_unfit_phi );
    tree_->Branch("JpsiMu_mu3_unfit_mass", &JpsiMu_mu3_unfit_mass );
    tree_->Branch("JpsiMu_mu3_doca2mu1", &JpsiMu_mu3_doca2mu1 );
    tree_->Branch("JpsiMu_mu3_doca2mu2", &JpsiMu_mu3_doca2mu2 );
    tree_->Branch("JpsiMu_mu3_q", &JpsiMu_mu3_q );
    tree_->Branch("JpsiMu_mu3_isLoose"  , &JpsiMu_mu3_isLoose   );
    tree_->Branch("JpsiMu_mu3_isTight"  , &JpsiMu_mu3_isTight   );
    tree_->Branch("JpsiMu_mu3_isPF"     , &JpsiMu_mu3_isPF      );
    tree_->Branch("JpsiMu_mu3_isGlobal" , &JpsiMu_mu3_isGlobal  );
    tree_->Branch("JpsiMu_mu3_isTracker", &JpsiMu_mu3_isTracker );
    tree_->Branch("JpsiMu_mu3_isSoft"   , &JpsiMu_mu3_isSoft    );
    tree_->Branch("JpsiMu_mu3_vx"   , &JpsiMu_mu3_vx    );
    tree_->Branch("JpsiMu_mu3_vy"   , &JpsiMu_mu3_vy    );
    tree_->Branch("JpsiMu_mu3_vz"   , &JpsiMu_mu3_vz    );
    tree_->Branch("JpsiMu_mu3_iso"   , &JpsiMu_mu3_iso    );
    tree_->Branch("JpsiMu_mu3_dbiso"   , &JpsiMu_mu3_dbiso    );

    tree_->Branch("JpsiMu_PV_vx", &JpsiMu_PV_vx );
    tree_->Branch("JpsiMu_PV_vy", &JpsiMu_PV_vy );
    tree_->Branch("JpsiMu_PV_vz", &JpsiMu_PV_vz );

    tree_->Branch("JpsiMu_bbPV_vx", &JpsiMu_bbPV_vx );
    tree_->Branch("JpsiMu_bbPV_vy", &JpsiMu_bbPV_vy );
    tree_->Branch("JpsiMu_bbPV_vz", &JpsiMu_bbPV_vz );

    //    tree_->Branch("JpsiMu_bbPV_refit_vx", &JpsiMu_bbPV_vx );
    //    tree_->Branch("JpsiMu_bbPV_refit_vy", &JpsiMu_bbPV_vy );
    //    tree_->Branch("JpsiMu_bbPV_refit_vz", &JpsiMu_bbPV_vz );


    tree_->Branch("JpsiMu_Jpsi_pt", &JpsiMu_Jpsi_pt );
    tree_->Branch("JpsiMu_Jpsi_eta", &JpsiMu_Jpsi_eta );
    tree_->Branch("JpsiMu_Jpsi_phi", &JpsiMu_Jpsi_phi );
    tree_->Branch("JpsiMu_Jpsi_mass", &JpsiMu_Jpsi_mass );
    tree_->Branch("JpsiMu_Jpsi_vprob", &JpsiMu_Jpsi_vprob );
    tree_->Branch("JpsiMu_Jpsi_lip", &JpsiMu_Jpsi_lip);
    tree_->Branch("JpsiMu_Jpsi_lips", &JpsiMu_Jpsi_lips);
    tree_->Branch("JpsiMu_Jpsi_pvip", &JpsiMu_Jpsi_pvip);
    tree_->Branch("JpsiMu_Jpsi_pvips", &JpsiMu_Jpsi_pvips);
    tree_->Branch("JpsiMu_Jpsi_fl3d", &JpsiMu_Jpsi_fl3d);
    tree_->Branch("JpsiMu_Jpsi_fls3d", &JpsiMu_Jpsi_fls3d);
    tree_->Branch("JpsiMu_Jpsi_alpha", &JpsiMu_Jpsi_alpha);
    tree_->Branch("JpsiMu_Jpsi_maxdoca", &JpsiMu_Jpsi_maxdoca);
    tree_->Branch("JpsiMu_Jpsi_mindoca", &JpsiMu_Jpsi_mindoca);
    tree_->Branch("JpsiMu_Jpsi_vx", &JpsiMu_Jpsi_vx );
    tree_->Branch("JpsiMu_Jpsi_vy", &JpsiMu_Jpsi_vy );
    tree_->Branch("JpsiMu_Jpsi_vz", &JpsiMu_Jpsi_vz );
    tree_->Branch("JpsiMu_Jpsi_unfit_pt", &JpsiMu_Jpsi_unfit_pt );
    tree_->Branch("JpsiMu_Jpsi_unfit_mass", &JpsiMu_Jpsi_unfit_mass );
    tree_->Branch("JpsiMu_Jpsi_unfit_vprob", &JpsiMu_Jpsi_unfit_vprob );
    tree_->Branch("JpsiMu_Jpsi_unfit_vx", &JpsiMu_Jpsi_unfit_vx );
    tree_->Branch("JpsiMu_Jpsi_unfit_vy", &JpsiMu_Jpsi_unfit_vy );
    tree_->Branch("JpsiMu_Jpsi_unfit_vz", &JpsiMu_Jpsi_unfit_vz );


    tree_->Branch("JpsiMu_B_pt", &JpsiMu_B_pt );
    tree_->Branch("JpsiMu_B_eta", &JpsiMu_B_eta );
    tree_->Branch("JpsiMu_B_phi", &JpsiMu_B_phi );
    tree_->Branch("JpsiMu_B_mass", &JpsiMu_B_mass );
    tree_->Branch("JpsiMu_B_mcorr", &JpsiMu_B_mcorr );
    tree_->Branch("JpsiMu_B_vprob", &JpsiMu_B_vprob );
    tree_->Branch("JpsiMu_B_lip", &JpsiMu_B_lip);
    tree_->Branch("JpsiMu_B_lips", &JpsiMu_B_lips);
    tree_->Branch("JpsiMu_B_pvip", &JpsiMu_B_pvip);
    tree_->Branch("JpsiMu_B_pvips", &JpsiMu_B_pvips);
    tree_->Branch("JpsiMu_B_fl3d", &JpsiMu_B_fl3d);
    tree_->Branch("JpsiMu_B_fls3d", &JpsiMu_B_fls3d);
    tree_->Branch("JpsiMu_B_alpha", &JpsiMu_B_alpha);
    tree_->Branch("JpsiMu_B_maxdoca", &JpsiMu_B_maxdoca);
    tree_->Branch("JpsiMu_B_mindoca", &JpsiMu_B_mindoca);
    tree_->Branch("JpsiMu_B_vx", &JpsiMu_B_vx );
    tree_->Branch("JpsiMu_B_vy", &JpsiMu_B_vy );
    tree_->Branch("JpsiMu_B_vz", &JpsiMu_B_vz );
    tree_->Branch("JpsiMu_B_iso", &JpsiMu_B_iso);
    tree_->Branch("JpsiMu_B_iso_ntracks", &JpsiMu_B_iso_ntracks );
    tree_->Branch("JpsiMu_B_iso_mindoca", &JpsiMu_B_iso_mindoca );
    tree_->Branch("JpsiMu_B_unfit_pt", &JpsiMu_B_unfit_pt );
    tree_->Branch("JpsiMu_B_unfit_mass", &JpsiMu_B_unfit_mass );
    tree_->Branch("JpsiMu_B_unfit_vprob", &JpsiMu_B_unfit_vprob );
    tree_->Branch("JpsiMu_B_unfit_vx", &JpsiMu_B_unfit_vx );
    tree_->Branch("JpsiMu_B_unfit_vy", &JpsiMu_B_unfit_vy );
    tree_->Branch("JpsiMu_B_unfit_vz", &JpsiMu_B_unfit_vz );
    tree_->Branch("JpsiMu_B_q2", &JpsiMu_B_q2 );
    tree_->Branch("JpsiMu_B_mm2", &JpsiMu_B_mm2 );
    tree_->Branch("JpsiMu_B_ptmiss", &JpsiMu_B_ptmiss );
    tree_->Branch("JpsiMu_B_Es", &JpsiMu_B_Es );
    tree_->Branch("JpsiMu_B_ptback", &JpsiMu_B_ptback );

  if (runFlags["runOnMC"] ){
    tree_->Branch("JpsiMu_genPV_vx", &JpsiMu_genPV_vx );
    tree_->Branch("JpsiMu_genPV_vy", &JpsiMu_genPV_vy );
    tree_->Branch("JpsiMu_genPV_vz", &JpsiMu_genPV_vz );

    tree_->Branch("JpsiMu_ngenmuons", &JpsiMu_ngenmuons);
    tree_->Branch("JpsiMu_isgenmatched", &JpsiMu_isgenmatched);
    tree_->Branch("JpsiMu_mu3_isgenmatched", &JpsiMu_mu3_isgenmatched);
    tree_->Branch("JpsiMu_q2_gen", &JpsiMu_q2_gen);
    tree_->Branch("JpsiMu_B_pt_gen", &JpsiMu_B_pt_gen);
    tree_->Branch("JpsiMu_B_eta_gen", &JpsiMu_B_eta_gen);
    tree_->Branch("JpsiMu_B_phi_gen", &JpsiMu_B_phi_gen);
    tree_->Branch("JpsiMu_B_mass_gen", &JpsiMu_B_mass_gen);

    tree_->Branch("JpsiMu_hammer_ff", &JpsiMu_hammer_ff);
    tree_->Branch("JpsiMu_hammer_ebe", &JpsiMu_hammer_ebe);
    tree_->Branch("JpsiMu_hammer_ebe_toy", &JpsiMu_hammer_ebe_toy);
    //    tree_->Branch("JpsiMu_hammer_ebe_up", &JpsiMu_hammer_ebe_up);
    //    tree_->Branch("JpsiMu_hammer_ebe_down", &JpsiMu_hammer_ebe_down);
    //    tree_->Branch("JpsiMu_hammer_ebe_rate_up", &JpsiMu_hammer_ebe_rate_up);
    //    tree_->Branch("JpsiMu_hammer_ebe_rate_down", &JpsiMu_hammer_ebe_rate_down);
///    tree_->Branch("JpsiMu_hammer_ebe_a0_up", &JpsiMu_hammer_ebe_a0_up);
///    tree_->Branch("JpsiMu_hammer_ebe_a0_down", &JpsiMu_hammer_ebe_a0_down);
///    tree_->Branch("JpsiMu_hammer_ebe_a1_up", &JpsiMu_hammer_ebe_a1_up);
///    tree_->Branch("JpsiMu_hammer_ebe_a1_down", &JpsiMu_hammer_ebe_a1_down);
///    tree_->Branch("JpsiMu_hammer_ebe_a2_up", &JpsiMu_hammer_ebe_a2_up);
///    tree_->Branch("JpsiMu_hammer_ebe_a2_down", &JpsiMu_hammer_ebe_a2_down);
///
///    tree_->Branch("JpsiMu_hammer_ebe_b0_up", &JpsiMu_hammer_ebe_b0_up);
///    tree_->Branch("JpsiMu_hammer_ebe_b0_down", &JpsiMu_hammer_ebe_b0_down);
///    tree_->Branch("JpsiMu_hammer_ebe_b1_up", &JpsiMu_hammer_ebe_b1_up);
///    tree_->Branch("JpsiMu_hammer_ebe_b1_down", &JpsiMu_hammer_ebe_b1_down);
///    tree_->Branch("JpsiMu_hammer_ebe_b2_up", &JpsiMu_hammer_ebe_b2_up);
///    tree_->Branch("JpsiMu_hammer_ebe_b2_down", &JpsiMu_hammer_ebe_b2_down);
///
///    tree_->Branch("JpsiMu_hammer_ebe_c1_up", &JpsiMu_hammer_ebe_c1_up);
///    tree_->Branch("JpsiMu_hammer_ebe_c1_down", &JpsiMu_hammer_ebe_c1_down);
///    tree_->Branch("JpsiMu_hammer_ebe_c2_up", &JpsiMu_hammer_ebe_c2_up);
///    tree_->Branch("JpsiMu_hammer_ebe_c2_down", &JpsiMu_hammer_ebe_c2_down);
///
///    tree_->Branch("JpsiMu_hammer_ebe_d0_up", &JpsiMu_hammer_ebe_d0_up);
///    tree_->Branch("JpsiMu_hammer_ebe_d0_down", &JpsiMu_hammer_ebe_d0_down);
///    tree_->Branch("JpsiMu_hammer_ebe_d1_up", &JpsiMu_hammer_ebe_d1_up);
///    tree_->Branch("JpsiMu_hammer_ebe_d1_down", &JpsiMu_hammer_ebe_d1_down);
///    tree_->Branch("JpsiMu_hammer_ebe_d2_up", &JpsiMu_hammer_ebe_d2_up);
///    tree_->Branch("JpsiMu_hammer_ebe_d2_down", &JpsiMu_hammer_ebe_d2_down);
  }


  }



  if (runFlags["doJpsiTau"]){
    tree_->Branch("HLT_BPH_isFired", &HLT_BPH_isFired );

    tree_->Branch("JpsiTau_nCandidates", &JpsiTau_nCandidates );


    tree_->Branch("JpsiTau_mu1_pt", &JpsiTau_mu1_pt );
    tree_->Branch("JpsiTau_mu1_eta", &JpsiTau_mu1_eta );
    tree_->Branch("JpsiTau_mu1_phi", &JpsiTau_mu1_phi );
    tree_->Branch("JpsiTau_mu1_mass", &JpsiTau_mu1_mass );
    tree_->Branch("JpsiTau_mu1_q", &JpsiTau_mu1_q );
    tree_->Branch("JpsiTau_mu1_isLoose"  , &JpsiTau_mu1_isLoose   );
    tree_->Branch("JpsiTau_mu1_isTight"  , &JpsiTau_mu1_isTight   );
    tree_->Branch("JpsiTau_mu1_isPF"     , &JpsiTau_mu1_isPF      );
    tree_->Branch("JpsiTau_mu1_isGlobal" , &JpsiTau_mu1_isGlobal  );
    tree_->Branch("JpsiTau_mu1_isTracker", &JpsiTau_mu1_isTracker );
    tree_->Branch("JpsiTau_mu1_isSoft"   , &JpsiTau_mu1_isSoft    );
    tree_->Branch("JpsiTau_mu1_vx"   , &JpsiTau_mu1_vx    );
    tree_->Branch("JpsiTau_mu1_vy"   , &JpsiTau_mu1_vy    );
    tree_->Branch("JpsiTau_mu1_vz"   , &JpsiTau_mu1_vz    );
    //    tree_->Branch("JpsiTau_mu1_iso"   , &JpsiTau_mu1_iso    );
    tree_->Branch("JpsiTau_mu1_dbiso"   , &JpsiTau_mu1_dbiso    );

    tree_->Branch("JpsiTau_mu2_pt", &JpsiTau_mu2_pt );
    tree_->Branch("JpsiTau_mu2_eta", &JpsiTau_mu2_eta );
    tree_->Branch("JpsiTau_mu2_phi", &JpsiTau_mu2_phi );
    tree_->Branch("JpsiTau_mu2_mass", &JpsiTau_mu2_mass );
    tree_->Branch("JpsiTau_mu2_q", &JpsiTau_mu2_q );
    tree_->Branch("JpsiTau_mu2_isLoose"  , &JpsiTau_mu2_isLoose   );
    tree_->Branch("JpsiTau_mu2_isTight"  , &JpsiTau_mu2_isTight   );
    tree_->Branch("JpsiTau_mu2_isPF"     , &JpsiTau_mu2_isPF      );
    tree_->Branch("JpsiTau_mu2_isGlobal" , &JpsiTau_mu2_isGlobal  );
    tree_->Branch("JpsiTau_mu2_isTracker", &JpsiTau_mu2_isTracker );
    tree_->Branch("JpsiTau_mu2_isSoft"   , &JpsiTau_mu2_isSoft    );
    tree_->Branch("JpsiTau_mu2_vx"   , &JpsiTau_mu2_vx    );
    tree_->Branch("JpsiTau_mu2_vy"   , &JpsiTau_mu2_vy    );
    tree_->Branch("JpsiTau_mu2_vz"   , &JpsiTau_mu2_vz    );
    //    tree_->Branch("JpsiTau_mu2_iso"   , &JpsiTau_mu2_iso    );
    tree_->Branch("JpsiTau_mu2_dbiso"   , &JpsiTau_mu2_dbiso    );

    //    tree_->Branch("JpsiTau_tau_fullfit_pt", &JpsiTau_tau_fullfit_pt );
    //    tree_->Branch("JpsiTau_tau_fullfit_eta", &JpsiTau_tau_fullfit_eta );
    //    tree_->Branch("JpsiTau_tau_fullfit_phi", &JpsiTau_tau_fullfit_phi );
    //    tree_->Branch("JpsiTau_tau_fullfit_mass", &JpsiTau_tau_fullfit_mass );
    tree_->Branch("JpsiTau_tau_pt", &JpsiTau_tau_pt );
    tree_->Branch("JpsiTau_tau_eta", &JpsiTau_tau_eta );
    tree_->Branch("JpsiTau_tau_phi", &JpsiTau_tau_phi );
    tree_->Branch("JpsiTau_tau_mass", &JpsiTau_tau_mass );
    tree_->Branch("JpsiTau_tau_q", &JpsiTau_tau_q );
    tree_->Branch("JpsiTau_tau_vx"   , &JpsiTau_tau_vx    );
    tree_->Branch("JpsiTau_tau_vy"   , &JpsiTau_tau_vy    );
    tree_->Branch("JpsiTau_tau_vz"   , &JpsiTau_tau_vz    );

    tree_->Branch("JpsiTau_tau_max_dr_3prong", &JpsiTau_tau_max_dr_3prong);
    tree_->Branch("JpsiTau_tau_lip", &JpsiTau_tau_lip);
    tree_->Branch("JpsiTau_tau_lips", &JpsiTau_tau_lips);
    tree_->Branch("JpsiTau_tau_pvip", &JpsiTau_tau_pvip);
    tree_->Branch("JpsiTau_tau_pvips", &JpsiTau_tau_pvips);
    tree_->Branch("JpsiTau_tau_fl3d", &JpsiTau_tau_fl3d);
    tree_->Branch("JpsiTau_tau_fls3d", &JpsiTau_tau_fls3d);
    tree_->Branch("JpsiTau_tau_alpha", &JpsiTau_tau_alpha);
    tree_->Branch("JpsiTau_tau_vprob", &JpsiTau_tau_vprob);

    tree_->Branch("JpsiTau_tau_fl3d_wjpsi", &JpsiTau_tau_fl3d_wjpsi);
    tree_->Branch("JpsiTau_tau_fls3d_wjpsi", &JpsiTau_tau_fls3d_wjpsi);

//    tree_->Branch("JpsiTau_tau_dr1", &JpsiTau_tau_dr1);
//    tree_->Branch("JpsiTau_tau_dr2", &JpsiTau_tau_dr2);
//    tree_->Branch("JpsiTau_tau_dr3", &JpsiTau_tau_dr3);
//    tree_->Branch("JpsiTau_tau_ptres1", &JpsiTau_tau_ptres1);
//    tree_->Branch("JpsiTau_tau_ptres2", &JpsiTau_tau_ptres2);
//    tree_->Branch("JpsiTau_tau_ptres3", &JpsiTau_tau_ptres3);
//    tree_->Branch("JpsiTau_tau_matched_ppdgId", &JpsiTau_tau_matched_ppdgId);
//    tree_->Branch("JpsiTau_tau_matched_gentaupt", &JpsiTau_tau_matched_gentaupt);
    //    tree_->Branch("JpsiTau_tau_gentaupt", &JpsiTau_tau_gentaupt);
    tree_->Branch("JpsiTau_tau_sumofdnn", &JpsiTau_tau_sumofdnn);
    tree_->Branch("JpsiTau_tau_sumofdnn_1prong", &JpsiTau_tau_sumofdnn_1prong);
    tree_->Branch("JpsiTau_tau_sumofdnn_otherB", &JpsiTau_tau_sumofdnn_otherB);
    tree_->Branch("JpsiTau_tau_sumofdnn_pu", &JpsiTau_tau_sumofdnn_pu);
    //    tree_->Branch("JpsiTau_tau_sumofdnn_old", &JpsiTau_tau_sumofdnn_old);
    //    tree_->Branch("JpsiTau_tau_sumofdnn_others", &JpsiTau_tau_sumofdnn_others);
    //    tree_->Branch("JpsiTau_tau_pi1_dnn", &JpsiTau_tau_pi1_dnn );
    //    tree_->Branch("JpsiTau_tau_pi2_dnn", &JpsiTau_tau_pi2_dnn );
    //    tree_->Branch("JpsiTau_tau_pi3_dnn", &JpsiTau_tau_pi3_dnn );

    //    tree_->Branch("JpsiTau_tau_pi1_doca", &JpsiTau_tau_pi1_doca );
    //    tree_->Branch("JpsiTau_tau_pi2_doca", &JpsiTau_tau_pi2_doca );
    //    tree_->Branch("JpsiTau_tau_pi3_doca", &JpsiTau_tau_pi3_doca );

    //    tree_->Branch("JpsiTau_tau_pi1_pv", &JpsiTau_tau_pi1_pv );
    //    tree_->Branch("JpsiTau_tau_pi2_pv", &JpsiTau_tau_pi2_pv );
    //    tree_->Branch("JpsiTau_tau_pi3_pv", &JpsiTau_tau_pi3_pv );


    tree_->Branch("JpsiTau_tau_rhomass1", &JpsiTau_tau_rhomass1 );
    tree_->Branch("JpsiTau_tau_rhomass2", &JpsiTau_tau_rhomass2 );

    tree_->Branch("JpsiTau_tau_pi1_pt", &JpsiTau_tau_pi1_pt );
    tree_->Branch("JpsiTau_tau_pi1_eta", &JpsiTau_tau_pi1_eta );
    tree_->Branch("JpsiTau_tau_pi1_phi", &JpsiTau_tau_pi1_phi );
    tree_->Branch("JpsiTau_tau_pi1_mass", &JpsiTau_tau_pi1_mass );
    tree_->Branch("JpsiTau_tau_pi1_q", &JpsiTau_tau_pi1_q );

    tree_->Branch("JpsiTau_tau_pi1_doca3d", &JpsiTau_tau_pi1_doca3d);
    tree_->Branch("JpsiTau_tau_pi1_doca3de", &JpsiTau_tau_pi1_doca3de);
    tree_->Branch("JpsiTau_tau_pi1_doca2d", &JpsiTau_tau_pi1_doca2d);
    tree_->Branch("JpsiTau_tau_pi1_doca2de", &JpsiTau_tau_pi1_doca2de);
    //    tree_->Branch("JpsiTau_tau_pi1_doca1d", &JpsiTau_tau_pi1_doca1d);
    //    tree_->Branch("JpsiTau_tau_pi1_doca1de", &JpsiTau_tau_pi1_doca1de);
    //    tree_->Branch("JpsiTau_tau_pi1_isRight", &JpsiTau_tau_pi1_isRight);
    tree_->Branch("JpsiTau_tau_pi1_dz", &JpsiTau_tau_pi1_dz);
    tree_->Branch("JpsiTau_tau_pi1_near_dz", &JpsiTau_tau_pi1_near_dz);
    tree_->Branch("JpsiTau_tau_pi1_isAssociate", &JpsiTau_tau_pi1_isAssociate);
    tree_->Branch("JpsiTau_tau_pi1_pvAssociationQuality", &JpsiTau_tau_pi1_pvAssociationQuality);
    tree_->Branch("JpsiTau_tau_pi1_isBdecay", &JpsiTau_tau_pi1_isBdecay);
    tree_->Branch("JpsiTau_tau_pi1_isBdecaypdg", &JpsiTau_tau_pi1_isBdecaypdg);
    tree_->Branch("JpsiTau_tau_pi1_isBdecayppdg", &JpsiTau_tau_pi1_isBdecayppdg);
    tree_->Branch("JpsiTau_tau_pi1_isSignal", &JpsiTau_tau_pi1_isSignal);
    tree_->Branch("JpsiTau_tau_pi1_nprong", &JpsiTau_tau_pi1_nprong);
    tree_->Branch("JpsiTau_tau_pi1_nprong_pi0", &JpsiTau_tau_pi1_nprong_pi0);
    tree_->Branch("JpsiTau_tau_pi1_dnn", &JpsiTau_tau_pi1_dnn);
    tree_->Branch("JpsiTau_tau_pi1_dnn_1prong", &JpsiTau_tau_pi1_dnn_1prong);
    tree_->Branch("JpsiTau_tau_pi1_dnn_otherB", &JpsiTau_tau_pi1_dnn_otherB);
    tree_->Branch("JpsiTau_tau_pi1_dnn_pu", &JpsiTau_tau_pi1_dnn_pu);
    //    tree_->Branch("JpsiTau_tau_pi1_dnn_old", &JpsiTau_tau_pi1_dnn_old);
    tree_->Branch("JpsiTau_tau_pi1_trigMatch", &JpsiTau_tau_pi1_trigMatch);
    tree_->Branch("JpsiTau_tau_pi1_trigMatch_dr", &JpsiTau_tau_pi1_trigMatch_dr);


    tree_->Branch("JpsiTau_tau_pi2_pt", &JpsiTau_tau_pi2_pt );
    tree_->Branch("JpsiTau_tau_pi2_eta", &JpsiTau_tau_pi2_eta );
    tree_->Branch("JpsiTau_tau_pi2_phi", &JpsiTau_tau_pi2_phi );
    tree_->Branch("JpsiTau_tau_pi2_mass", &JpsiTau_tau_pi2_mass );
    tree_->Branch("JpsiTau_tau_pi2_q", &JpsiTau_tau_pi2_q );

    tree_->Branch("JpsiTau_tau_pi2_doca3d", &JpsiTau_tau_pi2_doca3d);
    tree_->Branch("JpsiTau_tau_pi2_doca3de", &JpsiTau_tau_pi2_doca3de);
    tree_->Branch("JpsiTau_tau_pi2_doca2d", &JpsiTau_tau_pi2_doca2d);
    tree_->Branch("JpsiTau_tau_pi2_doca2de", &JpsiTau_tau_pi2_doca2de);
    //    tree_->Branch("JpsiTau_tau_pi2_doca1d", &JpsiTau_tau_pi2_doca1d);
    //    tree_->Branch("JpsiTau_tau_pi2_doca1de", &JpsiTau_tau_pi2_doca1de);
    //    tree_->Branch("JpsiTau_tau_pi2_isRight", &JpsiTau_tau_pi2_isRight);
    tree_->Branch("JpsiTau_tau_pi2_dz", &JpsiTau_tau_pi2_dz);
    tree_->Branch("JpsiTau_tau_pi2_near_dz", &JpsiTau_tau_pi2_near_dz);
    tree_->Branch("JpsiTau_tau_pi2_isAssociate", &JpsiTau_tau_pi2_isAssociate);
    tree_->Branch("JpsiTau_tau_pi2_pvAssociationQuality", &JpsiTau_tau_pi2_pvAssociationQuality);
    tree_->Branch("JpsiTau_tau_pi2_isBdecay", &JpsiTau_tau_pi2_isBdecay);
    tree_->Branch("JpsiTau_tau_pi2_isBdecaypdg", &JpsiTau_tau_pi2_isBdecaypdg);
    tree_->Branch("JpsiTau_tau_pi2_isBdecayppdg", &JpsiTau_tau_pi2_isBdecayppdg);
    tree_->Branch("JpsiTau_tau_pi2_isSignal", &JpsiTau_tau_pi2_isSignal);
    tree_->Branch("JpsiTau_tau_pi2_nprong", &JpsiTau_tau_pi2_nprong);
    tree_->Branch("JpsiTau_tau_pi2_nprong_pi0", &JpsiTau_tau_pi2_nprong_pi0);
    tree_->Branch("JpsiTau_tau_pi2_dnn", &JpsiTau_tau_pi2_dnn);
    tree_->Branch("JpsiTau_tau_pi2_dnn_1prong", &JpsiTau_tau_pi2_dnn_1prong);
    tree_->Branch("JpsiTau_tau_pi2_dnn_otherB", &JpsiTau_tau_pi2_dnn_otherB);
    tree_->Branch("JpsiTau_tau_pi2_dnn_pu", &JpsiTau_tau_pi2_dnn_pu);
    //    tree_->Branch("JpsiTau_tau_pi2_dnn_old", &JpsiTau_tau_pi2_dnn_old);
    tree_->Branch("JpsiTau_tau_pi2_trigMatch", &JpsiTau_tau_pi2_trigMatch);
    tree_->Branch("JpsiTau_tau_pi2_trigMatch_dr", &JpsiTau_tau_pi2_trigMatch_dr);


    tree_->Branch("JpsiTau_tau_pi3_pt", &JpsiTau_tau_pi3_pt );
    tree_->Branch("JpsiTau_tau_pi3_eta", &JpsiTau_tau_pi3_eta );
    tree_->Branch("JpsiTau_tau_pi3_phi", &JpsiTau_tau_pi3_phi );
    tree_->Branch("JpsiTau_tau_pi3_mass", &JpsiTau_tau_pi3_mass );
    tree_->Branch("JpsiTau_tau_pi3_q", &JpsiTau_tau_pi3_q );

    tree_->Branch("JpsiTau_tau_pi3_doca3d", &JpsiTau_tau_pi3_doca3d);
    tree_->Branch("JpsiTau_tau_pi3_doca3de", &JpsiTau_tau_pi3_doca3de);
    tree_->Branch("JpsiTau_tau_pi3_doca2d", &JpsiTau_tau_pi3_doca2d);
    tree_->Branch("JpsiTau_tau_pi3_doca2de", &JpsiTau_tau_pi3_doca2de);
    //    tree_->Branch("JpsiTau_tau_pi3_doca1d", &JpsiTau_tau_pi3_doca1d);
    //    tree_->Branch("JpsiTau_tau_pi3_doca1de", &JpsiTau_tau_pi3_doca1de);
    //    tree_->Branch("JpsiTau_tau_pi3_isRight", &JpsiTau_tau_pi3_isRight);
    tree_->Branch("JpsiTau_tau_pi3_dz", &JpsiTau_tau_pi3_dz);
    tree_->Branch("JpsiTau_tau_pi3_near_dz", &JpsiTau_tau_pi3_near_dz);
    tree_->Branch("JpsiTau_tau_pi3_isAssociate", &JpsiTau_tau_pi3_isAssociate);
    tree_->Branch("JpsiTau_tau_pi3_pvAssociationQuality", &JpsiTau_tau_pi3_pvAssociationQuality);
    tree_->Branch("JpsiTau_tau_pi3_isBdecay", &JpsiTau_tau_pi3_isBdecay);
    tree_->Branch("JpsiTau_tau_pi3_isBdecaypdg", &JpsiTau_tau_pi3_isBdecaypdg);
    tree_->Branch("JpsiTau_tau_pi3_isBdecayppdg", &JpsiTau_tau_pi3_isBdecayppdg);
    tree_->Branch("JpsiTau_tau_pi3_isSignal", &JpsiTau_tau_pi3_isSignal);
    tree_->Branch("JpsiTau_tau_pi3_nprong", &JpsiTau_tau_pi3_nprong);
    tree_->Branch("JpsiTau_tau_pi3_nprong_pi0", &JpsiTau_tau_pi3_nprong_pi0);

    tree_->Branch("JpsiTau_tau_pi3_dnn", &JpsiTau_tau_pi3_dnn);
    tree_->Branch("JpsiTau_tau_pi3_dnn_1prong", &JpsiTau_tau_pi3_dnn_1prong);
    tree_->Branch("JpsiTau_tau_pi3_dnn_otherB", &JpsiTau_tau_pi3_dnn_otherB);
    tree_->Branch("JpsiTau_tau_pi3_dnn_pu", &JpsiTau_tau_pi3_dnn_pu);
    //    tree_->Branch("JpsiTau_tau_pi3_dnn_old", &JpsiTau_tau_pi3_dnn_old);
    tree_->Branch("JpsiTau_tau_pi3_trigMatch", &JpsiTau_tau_pi3_trigMatch);
    tree_->Branch("JpsiTau_tau_pi3_trigMatch_dr", &JpsiTau_tau_pi3_trigMatch_dr);


    tree_->Branch("JpsiTau_tau_delta_chi2", &JpsiTau_tau_delta_chi2);
    tree_->Branch("JpsiTau_tau_delta_n_ch", &JpsiTau_tau_delta_n_ch);
    tree_->Branch("JpsiTau_tau_delta_n_mu", &JpsiTau_tau_delta_n_mu);
    tree_->Branch("JpsiTau_tau_vweight", &JpsiTau_tau_vweight);
    tree_->Branch("JpsiTau_tau_refit_vx", &JpsiTau_tau_refit_vx);
    tree_->Branch("JpsiTau_tau_refit_vy", &JpsiTau_tau_refit_vy);
    tree_->Branch("JpsiTau_tau_refit_vz", &JpsiTau_tau_refit_vz);
    tree_->Branch("JpsiTau_tau_refit_chi2", &JpsiTau_tau_refit_chi2);
    tree_->Branch("JpsiTau_tau_refit_ndof", &JpsiTau_tau_refit_ndof);
    tree_->Branch("JpsiTau_tau_refit_rho", &JpsiTau_tau_refit_rho);

    tree_->Branch("JpsiTau_tau_iso", &JpsiTau_tau_iso);
    tree_->Branch("JpsiTau_tau_iso_ntracks", &JpsiTau_tau_iso_ntracks );
    tree_->Branch("JpsiTau_tau_iso_mindoca", &JpsiTau_tau_iso_mindoca );

    //  std::vector<float> JpsiTau_tau_delta_chi2;
//  std::vector<float> JpsiTau_tau_delta_n;
//  std::vector<float> JpsiTau_tau_vweight;
//  std::vector<float> JpsiTau_tau_refit_vx;
//  std::vector<float> JpsiTau_tau_refit_vy;
//  std::vector<float> JpsiTau_tau_refit_vz;


    tree_->Branch("JpsiTau_ptbal", &JpsiTau_ptbal );
    tree_->Branch("JpsiTau_jpsi_tau_alpha", &JpsiTau_jpsi_tau_alpha );


    tree_->Branch("JpsiTau_PV_vx", &JpsiTau_PV_vx );
    tree_->Branch("JpsiTau_PV_vy", &JpsiTau_PV_vy );
    tree_->Branch("JpsiTau_PV_vz", &JpsiTau_PV_vz );

    tree_->Branch("JpsiTau_bbPV_vx", &JpsiTau_bbPV_vx );
    tree_->Branch("JpsiTau_bbPV_vy", &JpsiTau_bbPV_vy );
    tree_->Branch("JpsiTau_bbPV_vz", &JpsiTau_bbPV_vz );
    tree_->Branch("JpsiTau_bbPV_chi2", &JpsiTau_bbPV_chi2 );
    tree_->Branch("JpsiTau_bbPV_ndof", &JpsiTau_bbPV_ndof );
    tree_->Branch("JpsiTau_bbPV_rho", &JpsiTau_bbPV_rho );

    //    tree_->Branch("JpsiTau_bbPV_refit_vx", &JpsiTau_bbPV_vx );
    //    tree_->Branch("JpsiTau_bbPV_refit_vy", &JpsiTau_bbPV_vy );
    //    tree_->Branch("JpsiTau_bbPV_refit_vz", &JpsiTau_bbPV_vz );



    tree_->Branch("JpsiTau_Jpsi_pt", &JpsiTau_Jpsi_pt );
    tree_->Branch("JpsiTau_Jpsi_eta", &JpsiTau_Jpsi_eta );
    tree_->Branch("JpsiTau_Jpsi_phi", &JpsiTau_Jpsi_phi );
    tree_->Branch("JpsiTau_Jpsi_mass", &JpsiTau_Jpsi_mass );
    tree_->Branch("JpsiTau_Jpsi_vprob", &JpsiTau_Jpsi_vprob );
    tree_->Branch("JpsiTau_Jpsi_lip", &JpsiTau_Jpsi_lip);
    tree_->Branch("JpsiTau_Jpsi_lips", &JpsiTau_Jpsi_lips);
    tree_->Branch("JpsiTau_Jpsi_pvip", &JpsiTau_Jpsi_pvip);
    tree_->Branch("JpsiTau_Jpsi_pvips", &JpsiTau_Jpsi_pvips);
    tree_->Branch("JpsiTau_Jpsi_fl3d", &JpsiTau_Jpsi_fl3d);
    tree_->Branch("JpsiTau_Jpsi_fls3d", &JpsiTau_Jpsi_fls3d);
    tree_->Branch("JpsiTau_Jpsi_alpha", &JpsiTau_Jpsi_alpha);
    tree_->Branch("JpsiTau_Jpsi_maxdoca", &JpsiTau_Jpsi_maxdoca);
    tree_->Branch("JpsiTau_Jpsi_mindoca", &JpsiTau_Jpsi_mindoca);
    tree_->Branch("JpsiTau_Jpsi_vx", &JpsiTau_Jpsi_vx );
    tree_->Branch("JpsiTau_Jpsi_vy", &JpsiTau_Jpsi_vy );
    tree_->Branch("JpsiTau_Jpsi_vz", &JpsiTau_Jpsi_vz );
//    tree_->Branch("JpsiTau_Jpsi_unfit_pt", &JpsiTau_Jpsi_unfit_pt );
//    tree_->Branch("JpsiTau_Jpsi_unfit_mass", &JpsiTau_Jpsi_unfit_mass );
//    tree_->Branch("JpsiTau_Jpsi_unfit_vprob", &JpsiTau_Jpsi_unfit_vprob );
//    tree_->Branch("JpsiTau_Jpsi_unfit_vx", &JpsiTau_Jpsi_unfit_vx );
//    tree_->Branch("JpsiTau_Jpsi_unfit_vy", &JpsiTau_Jpsi_unfit_vy );
//    tree_->Branch("JpsiTau_Jpsi_unfit_vz", &JpsiTau_Jpsi_unfit_vz );


    tree_->Branch("JpsiTau_B_pt", &JpsiTau_B_pt );
    tree_->Branch("JpsiTau_B_eta", &JpsiTau_B_eta );
    tree_->Branch("JpsiTau_B_phi", &JpsiTau_B_phi );
    tree_->Branch("JpsiTau_B_mass", &JpsiTau_B_mass );
    tree_->Branch("JpsiTau_B_mcorr", &JpsiTau_B_mcorr );
    tree_->Branch("JpsiTau_B_vprob", &JpsiTau_B_vprob );
    tree_->Branch("JpsiTau_B_lip", &JpsiTau_B_lip);
    tree_->Branch("JpsiTau_B_lips", &JpsiTau_B_lips);
    tree_->Branch("JpsiTau_B_pvip", &JpsiTau_B_pvip);
    tree_->Branch("JpsiTau_B_pvips", &JpsiTau_B_pvips);
    tree_->Branch("JpsiTau_B_fl3d", &JpsiTau_B_fl3d);
    tree_->Branch("JpsiTau_B_fls3d", &JpsiTau_B_fls3d);
    tree_->Branch("JpsiTau_B_alpha", &JpsiTau_B_alpha);
    tree_->Branch("JpsiTau_B_maxdoca", &JpsiTau_B_maxdoca);
    tree_->Branch("JpsiTau_B_mindoca", &JpsiTau_B_mindoca);
    tree_->Branch("JpsiTau_B_vx", &JpsiTau_B_vx );
    tree_->Branch("JpsiTau_B_vy", &JpsiTau_B_vy );
    tree_->Branch("JpsiTau_B_vz", &JpsiTau_B_vz );

    //    tree_->Branch("JpsiTau_B_iso_nocut", &JpsiTau_B_iso_nocut);
    //    tree_->Branch("JpsiTau_B_iso_ntracks_nocut", &JpsiTau_B_iso_ntracks_nocut );
    //    tree_->Branch("JpsiTau_B_iso_mindoca_nocut", &JpsiTau_B_iso_mindoca_nocut );

    tree_->Branch("JpsiTau_B_q2", &JpsiTau_B_q2 );
    tree_->Branch("JpsiTau_B_mm2", &JpsiTau_B_mm2 );
    tree_->Branch("JpsiTau_B_ptmiss", &JpsiTau_B_ptmiss );
    tree_->Branch("JpsiTau_B_Es", &JpsiTau_B_Es );
    tree_->Branch("JpsiTau_B_ptback", &JpsiTau_B_ptback );

    tree_->Branch("JpsiTau_B_pt_simple", &JpsiTau_B_pt_simple );
    tree_->Branch("JpsiTau_B_eta_simple", &JpsiTau_B_eta_simple );
    tree_->Branch("JpsiTau_B_phi_simple", &JpsiTau_B_phi_simple );
    tree_->Branch("JpsiTau_B_mass_simple", &JpsiTau_B_mass_simple );
    tree_->Branch("JpsiTau_B_q2_simple", &JpsiTau_B_q2_simple );
    tree_->Branch("JpsiTau_B_mm2_simple", &JpsiTau_B_mm2_simple );
    tree_->Branch("JpsiTau_B_ptmiss_simple", &JpsiTau_B_ptmiss_simple );
    tree_->Branch("JpsiTau_B_Es_simple", &JpsiTau_B_Es_simple );
    tree_->Branch("JpsiTau_B_ptback_simple", &JpsiTau_B_ptback_simple );


  if (runFlags["runOnMC"] ){
      ////Adding weight based on B decay chain for B generic bkg sample                                                                                                                                                                            
      tree_->Branch( "genWeightBkgB"             , &genWeightBkgB       );
    tree_->Branch("JpsiTau_genPV_vx", &JpsiTau_genPV_vx );
    tree_->Branch("JpsiTau_genPV_vy", &JpsiTau_genPV_vy );
    tree_->Branch("JpsiTau_genPV_vz", &JpsiTau_genPV_vz );
    tree_->Branch("JpsiTau_genSV_vx", &JpsiTau_genSV_vx );
    tree_->Branch("JpsiTau_genSV_vy", &JpsiTau_genSV_vy );
    tree_->Branch("JpsiTau_genSV_vz", &JpsiTau_genSV_vz );
    tree_->Branch("JpsiTau_ngenmuons", &JpsiTau_ngenmuons);
    tree_->Branch("JpsiTau_isgenmatched", &JpsiTau_isgenmatched);
    //    tree_->Branch("JpsiTau_isgen3", &JpsiTau_isgen3);
    //    tree_->Branch("JpsiTau_isgen3matched", &JpsiTau_isgen3matched);
    //    tree_->Branch("JpsiTau_ngentau3", &JpsiTau_ngentau3); 
    //    tree_->Branch("JpsiTau_ngentau", &JpsiTau_ngentau);
    //    tree_->Branch("JpsiTau_gentaupt", &JpsiTau_gentaupt);
    //    tree_->Branch("JpsiTau_gentaueta", &JpsiTau_gentaueta);
    //    tree_->Branch("JpsiTau_gentauphi", &JpsiTau_gentauphi);
    //    tree_->Branch("JpsiTau_gentaumass", &JpsiTau_gentaumass);
    //    tree_->Branch("JpsiTau_gentaupt_bd", &JpsiTau_gentaupt_bd);
    //    tree_->Branch("JpsiTau_gentaueta_bd", &JpsiTau_gentaueta_bd);
    //    tree_->Branch("JpsiTau_gentauphi_bd", &JpsiTau_gentauphi_bd);
    //    tree_->Branch("JpsiTau_gentaumass_bd", &JpsiTau_gentaumass_bd);
    //    tree_->Branch("JpsiTau_gentaudm", &JpsiTau_gentaudm);
    tree_->Branch("JpsiTau_q2_gen", &JpsiTau_q2_gen);
    tree_->Branch("JpsiTau_B_pt_gen", &JpsiTau_B_pt_gen);
    tree_->Branch("JpsiTau_B_eta_gen", &JpsiTau_B_eta_gen);
    tree_->Branch("JpsiTau_B_phi_gen", &JpsiTau_B_phi_gen);
    tree_->Branch("JpsiTau_B_mass_gen", &JpsiTau_B_mass_gen);

    tree_->Branch("JpsiTau_hammer_ff", &JpsiTau_hammer_ff);
    tree_->Branch("JpsiTau_hammer_ebe", &JpsiTau_hammer_ebe);
    tree_->Branch("JpsiTau_hammer_ebe_toy", &JpsiTau_hammer_ebe_toy);
    //    tree_->Branch("JpsiTau_hammer_ebe_up", &JpsiTau_hammer_ebe_up);
    //    tree_->Branch("JpsiTau_hammer_ebe_down", &JpsiTau_hammer_ebe_down);
    //    tree_->Branch("JpsiTau_hammer_ebe_rate_up", &JpsiTau_hammer_ebe_rate_up);
    //    tree_->Branch("JpsiTau_hammer_ebe_rate_down", &JpsiTau_hammer_ebe_rate_down);
//    tree_->Branch("JpsiTau_hammer_ebe_a0_up", &JpsiTau_hammer_ebe_a0_up);
//    tree_->Branch("JpsiTau_hammer_ebe_a0_down", &JpsiTau_hammer_ebe_a0_down);
//    tree_->Branch("JpsiTau_hammer_ebe_a1_up", &JpsiTau_hammer_ebe_a1_up);
//    tree_->Branch("JpsiTau_hammer_ebe_a1_down", &JpsiTau_hammer_ebe_a1_down);
//    tree_->Branch("JpsiTau_hammer_ebe_a2_up", &JpsiTau_hammer_ebe_a2_up);
//    tree_->Branch("JpsiTau_hammer_ebe_a2_down", &JpsiTau_hammer_ebe_a2_down);
//
//    tree_->Branch("JpsiTau_hammer_ebe_b0_up", &JpsiTau_hammer_ebe_b0_up);
//    tree_->Branch("JpsiTau_hammer_ebe_b0_down", &JpsiTau_hammer_ebe_b0_down);
//    tree_->Branch("JpsiTau_hammer_ebe_b1_up", &JpsiTau_hammer_ebe_b1_up);
//    tree_->Branch("JpsiTau_hammer_ebe_b1_down", &JpsiTau_hammer_ebe_b1_down);
//    tree_->Branch("JpsiTau_hammer_ebe_b2_up", &JpsiTau_hammer_ebe_b2_up);
//    tree_->Branch("JpsiTau_hammer_ebe_b2_down", &JpsiTau_hammer_ebe_b2_down);
//
//    tree_->Branch("JpsiTau_hammer_ebe_c1_up", &JpsiTau_hammer_ebe_c1_up);
//    tree_->Branch("JpsiTau_hammer_ebe_c1_down", &JpsiTau_hammer_ebe_c1_down);
//    tree_->Branch("JpsiTau_hammer_ebe_c2_up", &JpsiTau_hammer_ebe_c2_up);
//    tree_->Branch("JpsiTau_hammer_ebe_c2_down", &JpsiTau_hammer_ebe_c2_down);
//
//    tree_->Branch("JpsiTau_hammer_ebe_d0_up", &JpsiTau_hammer_ebe_d0_up);
//    tree_->Branch("JpsiTau_hammer_ebe_d0_down", &JpsiTau_hammer_ebe_d0_down);
//    tree_->Branch("JpsiTau_hammer_ebe_d1_up", &JpsiTau_hammer_ebe_d1_up);
//    tree_->Branch("JpsiTau_hammer_ebe_d1_down", &JpsiTau_hammer_ebe_d1_down);
//    tree_->Branch("JpsiTau_hammer_ebe_d2_up", &JpsiTau_hammer_ebe_d2_up);
//    tree_->Branch("JpsiTau_hammer_ebe_d2_down", &JpsiTau_hammer_ebe_d2_down);

  }

  tree_->Branch("JpsiTau_nch", &JpsiTau_nch);
  //  tree_->Branch("JpsiTau_nch_after_dnn", &JpsiTau_nch_after_dnn);
  tree_->Branch("JpsiTau_nch_before", &JpsiTau_nch_before);
  //  tree_->Branch("JpsiTau_nch_pvipcut", &JpsiTau_nch_pvipcut);
  //  tree_->Branch("JpsiTau_nch_pvipcut_dz0", &JpsiTau_nch_pvipcut_dz0);
  //  tree_->Branch("JpsiTau_nch_pvipcut_dz2p5", &JpsiTau_nch_pvipcut_dz2p5);

  //  tree_->Branch("JpsiTau_nch_qr", &JpsiTau_nch_qr);

  tree_->Branch("JpsiTau_gen_pion_pt", &JpsiTau_gen_pion_pt);
  tree_->Branch("JpsiTau_gen_pion_eta", &JpsiTau_gen_pion_eta);
  tree_->Branch("JpsiTau_gen_pion_phi", &JpsiTau_gen_pion_phi);
  tree_->Branch("JpsiTau_gen_pion_matched", &JpsiTau_gen_pion_matched);

  tree_->Branch("JpsiTau_gen_tau_pt", &JpsiTau_gen_tau_pt);
  tree_->Branch("JpsiTau_gen_tau_eta", &JpsiTau_gen_tau_eta);
  tree_->Branch("JpsiTau_gen_tau_phi", &JpsiTau_gen_tau_phi);
  tree_->Branch("JpsiTau_gen_tau_nprong", &JpsiTau_gen_tau_nprong);
  tree_->Branch("JpsiTau_gen_tau_nmatched", &JpsiTau_gen_tau_nmatched);


  //  tree_->Branch("JpsiTau_st_doca3d_max", &JpsiTau_st_doca3d_max);
  //  tree_->Branch("JpsiTau_st_doca2d_max", &JpsiTau_st_doca2d_max);
  //  tree_->Branch("JpsiTau_st_doca1d_max", &JpsiTau_st_doca1d_max);
  //  tree_->Branch("JpsiTau_st_doca3ds_max", &JpsiTau_st_doca3ds_max);
  //  tree_->Branch("JpsiTau_st_doca2ds_max", &JpsiTau_st_doca2ds_max);
  //  tree_->Branch("JpsiTau_st_doca1ds_max", &JpsiTau_st_doca1ds_max);
  //  tree_->Branch("JpsiTau_st_doca3de_max", &JpsiTau_st_doca3de_max);
  //  tree_->Branch("JpsiTau_st_doca2de_max", &JpsiTau_st_doca2de_max);
  //  tree_->Branch("JpsiTau_st_doca1de_max", &JpsiTau_st_doca1de_max);

  //  tree_->Branch("JpsiTau_st_dz_max", &JpsiTau_st_dz_max);
  tree_->Branch("JpsiTau_st_nch", &JpsiTau_st_nch);
  tree_->Branch("JpsiTau_st_nch_matched", &JpsiTau_st_nch_matched);
  tree_->Branch("JpsiTau_st_n_charged_pions", &JpsiTau_st_n_charged_pions);
  tree_->Branch("JpsiTau_st_n_neutral_pions", &JpsiTau_st_n_neutral_pions);
  tree_->Branch("JpsiTau_st_n_mu_decay", &JpsiTau_st_n_mu_decay);
  tree_->Branch("JpsiTau_st_n_e_decay", &JpsiTau_st_n_e_decay);
  tree_->Branch("JpsiTau_st_n_occurance", &JpsiTau_st_n_occurance);
  tree_->Branch("JpsiTau_st_decayid", &JpsiTau_st_decayid);
  tree_->Branch("JpsiTau_st_gentau_pt", &JpsiTau_st_gentau_pt);
  tree_->Branch("JpsiTau_st_gentau_eta", &JpsiTau_st_gentau_eta);
  tree_->Branch("JpsiTau_st_gentau_phi", &JpsiTau_st_gentau_phi);

  tree_->Branch("JpsiTau_st_genjpsi_pt", &JpsiTau_st_genjpsi_pt);
  tree_->Branch("JpsiTau_st_genjpsi_eta", &JpsiTau_st_genjpsi_eta);
  tree_->Branch("JpsiTau_st_genjpsi_phi", &JpsiTau_st_genjpsi_phi);


  //  tree_->Branch("JpsiTau_st_isRight", &JpsiTau_st_isRight);
  tree_->Branch("JpsiTau_st_isBdecay", &JpsiTau_st_isBdecay);
  tree_->Branch("JpsiTau_st_isBdecaypdg", &JpsiTau_st_isBdecaypdg);
  tree_->Branch("JpsiTau_st_isBdecayppdg", &JpsiTau_st_isBdecayppdg);
  tree_->Branch("JpsiTau_st_isSignal", &JpsiTau_st_isSignal);
  tree_->Branch("JpsiTau_st_nprong", &JpsiTau_st_nprong);
  tree_->Branch("JpsiTau_st_nprong_pi0", &JpsiTau_st_nprong_pi0);

  tree_->Branch("JpsiTau_st_idx", &JpsiTau_st_idx);

  tree_->Branch("JpsiTau_st_doca3d", &JpsiTau_st_doca3d);
  tree_->Branch("JpsiTau_st_doca2d", &JpsiTau_st_doca2d);
  //  tree_->Branch("JpsiTau_st_doca1d", &JpsiTau_st_doca1d);
  tree_->Branch("JpsiTau_st_doca3ds", &JpsiTau_st_doca3ds);
  tree_->Branch("JpsiTau_st_doca2ds", &JpsiTau_st_doca2ds);
  //  tree_->Branch("JpsiTau_st_doca1ds", &JpsiTau_st_doca1ds);
  tree_->Branch("JpsiTau_st_doca3de", &JpsiTau_st_doca3de);
  tree_->Branch("JpsiTau_st_doca2de", &JpsiTau_st_doca2de);
  //  tree_->Branch("JpsiTau_st_doca1de", &JpsiTau_st_doca1de);

  tree_->Branch("JpsiTau_st_dz", &JpsiTau_st_dz);
  tree_->Branch("JpsiTau_st_isAssociate", &JpsiTau_st_isAssociate);
  tree_->Branch("JpsiTau_st_near_dz", &JpsiTau_st_near_dz);
  tree_->Branch("JpsiTau_st_dr_jpsi", &JpsiTau_st_dr_jpsi);
  tree_->Branch("JpsiTau_st_trigMatch", &JpsiTau_st_trigMatch);
  tree_->Branch("JpsiTau_st_trigMatch_dr", &JpsiTau_st_trigMatch_dr);
  tree_->Branch("JpsiTau_st_pvAssociationQuality", &JpsiTau_st_pvAssociationQuality);
  tree_->Branch("JpsiTau_st_pt", &JpsiTau_st_pt);
  tree_->Branch("JpsiTau_st_eta", &JpsiTau_st_eta);
  tree_->Branch("JpsiTau_st_phi", &JpsiTau_st_phi);
  tree_->Branch("JpsiTau_st_charge", &JpsiTau_st_charge);
  tree_->Branch("JpsiTau_st_mass", &JpsiTau_st_mass);
  tree_->Branch("JpsiTau_st_dnn", &JpsiTau_st_dnn);
  tree_->Branch("JpsiTau_st_dnn_1prong", &JpsiTau_st_dnn_1prong);
  tree_->Branch("JpsiTau_st_dnn_otherB", &JpsiTau_st_dnn_otherB);
  tree_->Branch("JpsiTau_st_dnn_pu", &JpsiTau_st_dnn_pu);
  //  tree_->Branch("JpsiTau_st_dnn_old", &JpsiTau_st_dnn_old);
  tree_->Branch("JpsiTau_st_matchidx", &JpsiTau_st_matchidx);

  //  tree_->Branch("JpsiTau_perEVT_old", &JpsiTau_perEVT_old);
  tree_->Branch("JpsiTau_perEVT_otherB", &JpsiTau_perEVT_otherB);
  tree_->Branch("JpsiTau_perEVT_sig", &JpsiTau_perEVT_sig);
  tree_->Branch("JpsiTau_perEVT_leptonic", &JpsiTau_perEVT_leptonic);
  tree_->Branch("JpsiTau_perEVT_1prong", &JpsiTau_perEVT_1prong);

  //  tree_->Branch("JpsiTau_st_isOtherPV", &JpsiTau_st_isOtherPV);



//    tree_->Branch("JpsiTau_ed_pfeta", &JpsiTau_ed_pfeta);
//    tree_->Branch("JpsiTau_ed_pfphi", &JpsiTau_ed_pfphi);
//    tree_->Branch("JpsiTau_ed_isRight", &JpsiTau_ed_isRight);
//    tree_->Branch("JpsiTau_ed_id", &JpsiTau_ed_id);
//    tree_->Branch("JpsiTau_ed_pfdnn", &JpsiTau_ed_pfdnn);
//    tree_->Branch("JpsiTau_ed_genpt", &JpsiTau_ed_genpt);


  }



//  if (runFlags["doBsTauTau"]){
//    tree_->Branch("IsBsTauTau", &IsBsTauTau );
//
//    tree_->Branch("BsTauTau_nCandidates", &BsTauTau_nCandidates );
//
//    tree_->Branch("BsTauTau_mu1_pt", &BsTauTau_mu1_pt );
//    tree_->Branch("BsTauTau_mu1_eta", &BsTauTau_mu1_eta );
//    tree_->Branch("BsTauTau_mu1_phi", &BsTauTau_mu1_phi );
//    tree_->Branch("BsTauTau_mu1_mass", &BsTauTau_mu1_mass );
//    tree_->Branch("BsTauTau_mu1_q", &BsTauTau_mu1_q );
//    tree_->Branch("BsTauTau_mu1_isLoose"  , &BsTauTau_mu1_isLoose   );
//    tree_->Branch("BsTauTau_mu1_isTight"  , &BsTauTau_mu1_isTight   );
//    tree_->Branch("BsTauTau_mu1_isPF"     , &BsTauTau_mu1_isPF      );
//    tree_->Branch("BsTauTau_mu1_isGlobal" , &BsTauTau_mu1_isGlobal  );
//    tree_->Branch("BsTauTau_mu1_isTracker", &BsTauTau_mu1_isTracker );
//    tree_->Branch("BsTauTau_mu1_isSoft"   , &BsTauTau_mu1_isSoft    );
//    tree_->Branch("BsTauTau_mu1_vx"   , &BsTauTau_mu1_vx    );
//    tree_->Branch("BsTauTau_mu1_vy"   , &BsTauTau_mu1_vy    );
//    tree_->Branch("BsTauTau_mu1_vz"   , &BsTauTau_mu1_vz    );
//    tree_->Branch("BsTauTau_mu1_iso"   , &BsTauTau_mu1_iso    );
//    tree_->Branch("BsTauTau_mu1_dbiso"   , &BsTauTau_mu1_dbiso    );
//
//
//    tree_->Branch("BsTauTau_tau_pt", &BsTauTau_tau_pt );
//    tree_->Branch("BsTauTau_tau_eta", &BsTauTau_tau_eta );
//    tree_->Branch("BsTauTau_tau_phi", &BsTauTau_tau_phi );
//    tree_->Branch("BsTauTau_tau_mass", &BsTauTau_tau_mass );
//    tree_->Branch("BsTauTau_tau_rhomass1", &BsTauTau_tau_rhomass1 );
//    tree_->Branch("BsTauTau_tau_rhomass2", &BsTauTau_tau_rhomass2 );
//    tree_->Branch("BsTauTau_tau_q", &BsTauTau_tau_q );
//    tree_->Branch("BsTauTau_tau_vx"   , &BsTauTau_tau_vx    );
//    tree_->Branch("BsTauTau_tau_vy"   , &BsTauTau_tau_vy    );
//    tree_->Branch("BsTauTau_tau_vz"   , &BsTauTau_tau_vz    );
//
//    tree_->Branch("BsTauTau_tau_max_dr_3prong", &BsTauTau_tau_max_dr_3prong);
//    tree_->Branch("BsTauTau_tau_lip", &BsTauTau_tau_lip);
//    tree_->Branch("BsTauTau_tau_lips", &BsTauTau_tau_lips);
//    tree_->Branch("BsTauTau_tau_pvip", &BsTauTau_tau_pvip);
//    tree_->Branch("BsTauTau_tau_pvips", &BsTauTau_tau_pvips);
//    tree_->Branch("BsTauTau_tau_fl3d", &BsTauTau_tau_fl3d);
//    tree_->Branch("BsTauTau_tau_fls3d", &BsTauTau_tau_fls3d);
//    tree_->Branch("BsTauTau_tau_alpha", &BsTauTau_tau_alpha);
//    tree_->Branch("BsTauTau_tau_vprob", &BsTauTau_tau_vprob);
//    tree_->Branch("BsTauTau_tau_isRight", &BsTauTau_tau_isRight);
//    tree_->Branch("BsTauTau_tau_isRight1", &BsTauTau_tau_isRight1);
//    tree_->Branch("BsTauTau_tau_isRight2", &BsTauTau_tau_isRight2);
//    tree_->Branch("BsTauTau_tau_isRight3", &BsTauTau_tau_isRight3);
////    tree_->Branch("BsTauTau_tau_dr1", &BsTauTau_tau_dr1);
////    tree_->Branch("BsTauTau_tau_dr2", &BsTauTau_tau_dr2);
////    tree_->Branch("BsTauTau_tau_dr3", &BsTauTau_tau_dr3);
////    tree_->Branch("BsTauTau_tau_ptres1", &BsTauTau_tau_ptres1);
////    tree_->Branch("BsTauTau_tau_ptres2", &BsTauTau_tau_ptres2);
////    tree_->Branch("BsTauTau_tau_ptres3", &BsTauTau_tau_ptres3);
//    tree_->Branch("BsTauTau_tau_matched_ppdgId", &BsTauTau_tau_matched_ppdgId);
//    tree_->Branch("BsTauTau_tau_matched_gentaupt", &BsTauTau_tau_matched_gentaupt);
//    //    tree_->Branch("BsTauTau_tau_gentaupt", &BsTauTau_tau_gentaupt);
//    tree_->Branch("BsTauTau_tau_sumofdnn", &BsTauTau_tau_sumofdnn);
//    tree_->Branch("BsTauTau_tau_pfidx1", &BsTauTau_tau_pfidx1);
//    tree_->Branch("BsTauTau_tau_pfidx2", &BsTauTau_tau_pfidx2);
//    tree_->Branch("BsTauTau_tau_pfidx3", &BsTauTau_tau_pfidx3);
//    tree_->Branch("BsTauTau_tau_pi1_dnn", &BsTauTau_tau_pi1_dnn );
//    tree_->Branch("BsTauTau_tau_pi2_dnn", &BsTauTau_tau_pi2_dnn );
//    tree_->Branch("BsTauTau_tau_pi3_dnn", &BsTauTau_tau_pi3_dnn );
//
//    tree_->Branch("BsTauTau_tau_pi1_pt", &BsTauTau_tau_pi1_pt );
//    tree_->Branch("BsTauTau_tau_pi1_eta", &BsTauTau_tau_pi1_eta );
//    tree_->Branch("BsTauTau_tau_pi1_phi", &BsTauTau_tau_pi1_phi );
//    tree_->Branch("BsTauTau_tau_pi1_mass", &BsTauTau_tau_pi1_mass );
//
//    tree_->Branch("BsTauTau_tau_pi2_pt", &BsTauTau_tau_pi2_pt );
//    tree_->Branch("BsTauTau_tau_pi2_eta", &BsTauTau_tau_pi2_eta );
//    tree_->Branch("BsTauTau_tau_pi2_phi", &BsTauTau_tau_pi2_phi );
//    tree_->Branch("BsTauTau_tau_pi2_mass", &BsTauTau_tau_pi2_mass );
//
//    tree_->Branch("BsTauTau_tau_pi3_pt", &BsTauTau_tau_pi3_pt );
//    tree_->Branch("BsTauTau_tau_pi3_eta", &BsTauTau_tau_pi3_eta );
//    tree_->Branch("BsTauTau_tau_pi3_phi", &BsTauTau_tau_pi3_phi );
//    tree_->Branch("BsTauTau_tau_pi3_mass", &BsTauTau_tau_pi3_mass );
//
//
//    tree_->Branch("BsTauTau_PV_vx", &BsTauTau_PV_vx );
//    tree_->Branch("BsTauTau_PV_vy", &BsTauTau_PV_vy );
//    tree_->Branch("BsTauTau_PV_vz", &BsTauTau_PV_vz );
//
//    tree_->Branch("BsTauTau_bbPV_vx", &BsTauTau_bbPV_vx );
//    tree_->Branch("BsTauTau_bbPV_vy", &BsTauTau_bbPV_vy );
//    tree_->Branch("BsTauTau_bbPV_vz", &BsTauTau_bbPV_vz );
//
//    tree_->Branch("BsTauTau_bbPV_refit_vx", &BsTauTau_bbPV_vx );
//    tree_->Branch("BsTauTau_bbPV_refit_vy", &BsTauTau_bbPV_vy );
//    tree_->Branch("BsTauTau_bbPV_refit_vz", &BsTauTau_bbPV_vz );
//
//    tree_->Branch("BsTauTau_genPV_vx", &BsTauTau_genPV_vx );
//    tree_->Branch("BsTauTau_genPV_vy", &BsTauTau_genPV_vy );
//    tree_->Branch("BsTauTau_genPV_vz", &BsTauTau_genPV_vz );
//
//    tree_->Branch("BsTauTau_B_pt", &BsTauTau_B_pt );
//    tree_->Branch("BsTauTau_B_eta", &BsTauTau_B_eta );
//    tree_->Branch("BsTauTau_B_phi", &BsTauTau_B_phi );
//    tree_->Branch("BsTauTau_B_mass", &BsTauTau_B_mass );
//    tree_->Branch("BsTauTau_B_vprob", &BsTauTau_B_vprob );
//    tree_->Branch("BsTauTau_B_lip", &BsTauTau_B_lip);
//    tree_->Branch("BsTauTau_B_lips", &BsTauTau_B_lips);
//    tree_->Branch("BsTauTau_B_pvip", &BsTauTau_B_pvip);
//    tree_->Branch("BsTauTau_B_pvips", &BsTauTau_B_pvips);
//    tree_->Branch("BsTauTau_B_fl3d", &BsTauTau_B_fl3d);
//    tree_->Branch("BsTauTau_B_fls3d", &BsTauTau_B_fls3d);
//    tree_->Branch("BsTauTau_B_alpha", &BsTauTau_B_alpha);
//    tree_->Branch("BsTauTau_B_maxdoca", &BsTauTau_B_maxdoca);
//    tree_->Branch("BsTauTau_B_mindoca", &BsTauTau_B_mindoca);
//    tree_->Branch("BsTauTau_B_vx", &BsTauTau_B_vx );
//    tree_->Branch("BsTauTau_B_vy", &BsTauTau_B_vy );
//    tree_->Branch("BsTauTau_B_vz", &BsTauTau_B_vz );
//    tree_->Branch("BsTauTau_B_iso", &BsTauTau_B_iso);
//    tree_->Branch("BsTauTau_B_iso_ntracks", &BsTauTau_B_iso_ntracks );
//    tree_->Branch("BsTauTau_B_iso_mindoca", &BsTauTau_B_iso_mindoca );
//
//    tree_->Branch("BsTauTau_ngenmuons", &BsTauTau_ngenmuons);
//    tree_->Branch("BsTauTau_isgen3", &BsTauTau_isgen3);
//    tree_->Branch("BsTauTau_isgen3matched", &BsTauTau_isgen3matched);
//    tree_->Branch("BsTauTau_nch", &BsTauTau_nch);
//    tree_->Branch("BsTauTau_nch_after_dnn", &BsTauTau_nch_after_dnn);
//    tree_->Branch("BsTauTau_nch_before_dnn", &BsTauTau_nch_before_dnn);
//    tree_->Branch("BsTauTau_nch_qr", &BsTauTau_nch_qr);
//    tree_->Branch("BsTauTau_ngentau3", &BsTauTau_ngentau3); 
//    tree_->Branch("BsTauTau_ngentau", &BsTauTau_ngentau);
//    tree_->Branch("BsTauTau_gentaupt", &BsTauTau_gentaupt);
//    tree_->Branch("BsTauTau_gentaudm", &BsTauTau_gentaudm);
//
//  }
//
//
//  if (runFlags["doBsTauTauFH"]){
//    tree_->Branch("IsBsTauTauFH", &IsBsTauTauFH );
//
//
//    tree_->Branch("BsTauTauFH_nCandidates", &BsTauTauFH_nCandidates );
//    tree_->Branch("BsTauTauFH_ntaus", &BsTauTauFH_ntaus );
//
//    tree_->Branch("BsTauTauFH_mu1_pt", &BsTauTauFH_mu1_pt );
//    tree_->Branch("BsTauTauFH_mu1_eta", &BsTauTauFH_mu1_eta );
//    tree_->Branch("BsTauTauFH_mu1_phi", &BsTauTauFH_mu1_phi );
//    tree_->Branch("BsTauTauFH_mu1_mass", &BsTauTauFH_mu1_mass );
//    tree_->Branch("BsTauTauFH_mu1_q", &BsTauTauFH_mu1_q );
//    tree_->Branch("BsTauTauFH_mu1_isLoose"  , &BsTauTauFH_mu1_isLoose   );
//    tree_->Branch("BsTauTauFH_mu1_isTight"  , &BsTauTauFH_mu1_isTight   );
//    tree_->Branch("BsTauTauFH_mu1_isPF"     , &BsTauTauFH_mu1_isPF      );
//    tree_->Branch("BsTauTauFH_mu1_isGlobal" , &BsTauTauFH_mu1_isGlobal  );
//    tree_->Branch("BsTauTauFH_mu1_isTracker", &BsTauTauFH_mu1_isTracker );
//    tree_->Branch("BsTauTauFH_mu1_isSoft"   , &BsTauTauFH_mu1_isSoft    );
//    tree_->Branch("BsTauTauFH_mu1_vx"   , &BsTauTauFH_mu1_vx    );
//    tree_->Branch("BsTauTauFH_mu1_vy"   , &BsTauTauFH_mu1_vy    );
//    tree_->Branch("BsTauTauFH_mu1_vz"   , &BsTauTauFH_mu1_vz    );
//    tree_->Branch("BsTauTauFH_mu1_iso"   , &BsTauTauFH_mu1_iso    );
//    tree_->Branch("BsTauTauFH_mu1_dbiso"   , &BsTauTauFH_mu1_dbiso    );
//
//
//    tree_->Branch("BsTauTauFH_tau1_pt", &BsTauTauFH_tau1_pt );
//    tree_->Branch("BsTauTauFH_tau1_eta", &BsTauTauFH_tau1_eta );
//    tree_->Branch("BsTauTauFH_tau1_phi", &BsTauTauFH_tau1_phi );
//    tree_->Branch("BsTauTauFH_tau1_mass", &BsTauTauFH_tau1_mass );
//    tree_->Branch("BsTauTauFH_tau1_rhomass1", &BsTauTauFH_tau1_rhomass1 );
//    tree_->Branch("BsTauTauFH_tau1_rhomass2", &BsTauTauFH_tau1_rhomass2 );
//    tree_->Branch("BsTauTauFH_tau1_q", &BsTauTauFH_tau1_q );
//    tree_->Branch("BsTauTauFH_tau1_vx"   , &BsTauTauFH_tau1_vx    );
//    tree_->Branch("BsTauTauFH_tau1_vy"   , &BsTauTauFH_tau1_vy    );
//    tree_->Branch("BsTauTauFH_tau1_vz"   , &BsTauTauFH_tau1_vz    );
//
//    tree_->Branch("BsTauTauFH_tau1_max_dr_3prong", &BsTauTauFH_tau1_max_dr_3prong);
//    tree_->Branch("BsTauTauFH_tau1_lip", &BsTauTauFH_tau1_lip);
//    tree_->Branch("BsTauTauFH_tau1_lips", &BsTauTauFH_tau1_lips);
//    tree_->Branch("BsTauTauFH_tau1_pvip", &BsTauTauFH_tau1_pvip);
//    tree_->Branch("BsTauTauFH_tau1_pvips", &BsTauTauFH_tau1_pvips);
//    tree_->Branch("BsTauTauFH_tau1_fl3d", &BsTauTauFH_tau1_fl3d);
//    tree_->Branch("BsTauTauFH_tau1_fls3d", &BsTauTauFH_tau1_fls3d);
//    tree_->Branch("BsTauTauFH_tau1_alpha", &BsTauTauFH_tau1_alpha);
//    tree_->Branch("BsTauTauFH_tau1_vprob", &BsTauTauFH_tau1_vprob);
//    tree_->Branch("BsTauTauFH_tau1_isRight", &BsTauTauFH_tau1_isRight);
//    tree_->Branch("BsTauTauFH_tau1_matched_ppdgId", &BsTauTauFH_tau1_matched_ppdgId);
//    tree_->Branch("BsTauTauFH_tau1_matched_gentaupt", &BsTauTauFH_tau1_matched_gentaupt);
//    //    tree_->Branch("BsTauTauFH_tau1_gentaupt", &BsTauTauFH_tau1_gentaupt);
//    tree_->Branch("BsTauTauFH_tau1_sumofdnn", &BsTauTauFH_tau1_sumofdnn);
//    tree_->Branch("BsTauTauFH_tau1_pfidx1", &BsTauTauFH_tau1_pfidx1);
//    tree_->Branch("BsTauTauFH_tau1_pfidx2", &BsTauTauFH_tau1_pfidx2);
//    tree_->Branch("BsTauTauFH_tau1_pfidx3", &BsTauTauFH_tau1_pfidx3);
//    tree_->Branch("BsTauTauFH_tau1_pi1_dnn", &BsTauTauFH_tau1_pi1_dnn );
//    tree_->Branch("BsTauTauFH_tau1_pi2_dnn", &BsTauTauFH_tau1_pi2_dnn );
//    tree_->Branch("BsTauTauFH_tau1_pi3_dnn", &BsTauTauFH_tau1_pi3_dnn );
//
//    tree_->Branch("BsTauTauFH_tau1_pi1_pt", &BsTauTauFH_tau1_pi1_pt );
//    tree_->Branch("BsTauTauFH_tau1_pi1_eta", &BsTauTauFH_tau1_pi1_eta );
//    tree_->Branch("BsTauTauFH_tau1_pi1_phi", &BsTauTauFH_tau1_pi1_phi );
//    tree_->Branch("BsTauTauFH_tau1_pi1_mass", &BsTauTauFH_tau1_pi1_mass );
//
//    tree_->Branch("BsTauTauFH_tau1_pi2_pt", &BsTauTauFH_tau1_pi2_pt );
//    tree_->Branch("BsTauTauFH_tau1_pi2_eta", &BsTauTauFH_tau1_pi2_eta );
//    tree_->Branch("BsTauTauFH_tau1_pi2_phi", &BsTauTauFH_tau1_pi2_phi );
//    tree_->Branch("BsTauTauFH_tau1_pi2_mass", &BsTauTauFH_tau1_pi2_mass );
//
//    tree_->Branch("BsTauTauFH_tau1_pi3_pt", &BsTauTauFH_tau1_pi3_pt );
//    tree_->Branch("BsTauTauFH_tau1_pi3_eta", &BsTauTauFH_tau1_pi3_eta );
//    tree_->Branch("BsTauTauFH_tau1_pi3_phi", &BsTauTauFH_tau1_pi3_phi );
//    tree_->Branch("BsTauTauFH_tau1_pi3_mass", &BsTauTauFH_tau1_pi3_mass );
//
//    
//    tree_->Branch("BsTauTauFH_tau2_pt", &BsTauTauFH_tau2_pt );
//    tree_->Branch("BsTauTauFH_tau2_eta", &BsTauTauFH_tau2_eta );
//    tree_->Branch("BsTauTauFH_tau2_phi", &BsTauTauFH_tau2_phi );
//    tree_->Branch("BsTauTauFH_tau2_mass", &BsTauTauFH_tau2_mass );
//    tree_->Branch("BsTauTauFH_tau2_rhomass1", &BsTauTauFH_tau2_rhomass1 );
//    tree_->Branch("BsTauTauFH_tau2_rhomass2", &BsTauTauFH_tau2_rhomass2 );
//    tree_->Branch("BsTauTauFH_tau2_q", &BsTauTauFH_tau2_q );
//    tree_->Branch("BsTauTauFH_tau2_vx"   , &BsTauTauFH_tau2_vx    );
//    tree_->Branch("BsTauTauFH_tau2_vy"   , &BsTauTauFH_tau2_vy    );
//    tree_->Branch("BsTauTauFH_tau2_vz"   , &BsTauTauFH_tau2_vz    );
//
//    tree_->Branch("BsTauTauFH_tau2_max_dr_3prong", &BsTauTauFH_tau2_max_dr_3prong);
//    tree_->Branch("BsTauTauFH_tau2_lip", &BsTauTauFH_tau2_lip);
//    tree_->Branch("BsTauTauFH_tau2_lips", &BsTauTauFH_tau2_lips);
//    tree_->Branch("BsTauTauFH_tau2_pvip", &BsTauTauFH_tau2_pvip);
//    tree_->Branch("BsTauTauFH_tau2_pvips", &BsTauTauFH_tau2_pvips);
//    tree_->Branch("BsTauTauFH_tau2_fl3d", &BsTauTauFH_tau2_fl3d);
//    tree_->Branch("BsTauTauFH_tau2_fls3d", &BsTauTauFH_tau2_fls3d);
//    tree_->Branch("BsTauTauFH_tau2_alpha", &BsTauTauFH_tau2_alpha);
//    tree_->Branch("BsTauTauFH_tau2_vprob", &BsTauTauFH_tau2_vprob);
//    tree_->Branch("BsTauTauFH_tau2_isRight", &BsTauTauFH_tau2_isRight);
//    tree_->Branch("BsTauTauFH_tau2_matched_ppdgId", &BsTauTauFH_tau2_matched_ppdgId);
//    tree_->Branch("BsTauTauFH_tau2_matched_gentaupt", &BsTauTauFH_tau2_matched_gentaupt);
//    //    tree_->Branch("BsTauTauFH_tau2_gentaupt", &BsTauTauFH_tau2_gentaupt);
//    tree_->Branch("BsTauTauFH_tau2_sumofdnn", &BsTauTauFH_tau2_sumofdnn);
//    tree_->Branch("BsTauTauFH_tau2_pfidx1", &BsTauTauFH_tau2_pfidx1);
//    tree_->Branch("BsTauTauFH_tau2_pfidx2", &BsTauTauFH_tau2_pfidx2);
//    tree_->Branch("BsTauTauFH_tau2_pfidx3", &BsTauTauFH_tau2_pfidx3);
//    tree_->Branch("BsTauTauFH_tau2_pi1_dnn", &BsTauTauFH_tau2_pi1_dnn );
//    tree_->Branch("BsTauTauFH_tau2_pi2_dnn", &BsTauTauFH_tau2_pi2_dnn );
//    tree_->Branch("BsTauTauFH_tau2_pi3_dnn", &BsTauTauFH_tau2_pi3_dnn );
//
//    tree_->Branch("BsTauTauFH_tau2_pi1_pt", &BsTauTauFH_tau2_pi1_pt );
//    tree_->Branch("BsTauTauFH_tau2_pi1_eta", &BsTauTauFH_tau2_pi1_eta );
//    tree_->Branch("BsTauTauFH_tau2_pi1_phi", &BsTauTauFH_tau2_pi1_phi );
//    tree_->Branch("BsTauTauFH_tau2_pi1_mass", &BsTauTauFH_tau2_pi1_mass );
//
//    tree_->Branch("BsTauTauFH_tau2_pi2_pt", &BsTauTauFH_tau2_pi2_pt );
//    tree_->Branch("BsTauTauFH_tau2_pi2_eta", &BsTauTauFH_tau2_pi2_eta );
//    tree_->Branch("BsTauTauFH_tau2_pi2_phi", &BsTauTauFH_tau2_pi2_phi );
//    tree_->Branch("BsTauTauFH_tau2_pi2_mass", &BsTauTauFH_tau2_pi2_mass );
//
//    tree_->Branch("BsTauTauFH_tau2_pi3_pt", &BsTauTauFH_tau2_pi3_pt );
//    tree_->Branch("BsTauTauFH_tau2_pi3_eta", &BsTauTauFH_tau2_pi3_eta );
//    tree_->Branch("BsTauTauFH_tau2_pi3_phi", &BsTauTauFH_tau2_pi3_phi );
//    tree_->Branch("BsTauTauFH_tau2_pi3_mass", &BsTauTauFH_tau2_pi3_mass );
//
//
//    tree_->Branch("BsTauTauFH_PV_vx", &BsTauTauFH_PV_vx );
//    tree_->Branch("BsTauTauFH_PV_vy", &BsTauTauFH_PV_vy );
//    tree_->Branch("BsTauTauFH_PV_vz", &BsTauTauFH_PV_vz );
//
//    tree_->Branch("BsTauTauFH_bbPV_vx", &BsTauTauFH_bbPV_vx );
//    tree_->Branch("BsTauTauFH_bbPV_vy", &BsTauTauFH_bbPV_vy );
//    tree_->Branch("BsTauTauFH_bbPV_vz", &BsTauTauFH_bbPV_vz );
//
//    tree_->Branch("BsTauTauFH_bbPV_refit_vx", &BsTauTauFH_bbPV_vx );
//    tree_->Branch("BsTauTauFH_bbPV_refit_vy", &BsTauTauFH_bbPV_vy );
//    tree_->Branch("BsTauTauFH_bbPV_refit_vz", &BsTauTauFH_bbPV_vz );
//
//    tree_->Branch("BsTauTauFH_genPV_vx", &BsTauTauFH_genPV_vx );
//    tree_->Branch("BsTauTauFH_genPV_vy", &BsTauTauFH_genPV_vy );
//    tree_->Branch("BsTauTauFH_genPV_vz", &BsTauTauFH_genPV_vz );
//
//    tree_->Branch("BsTauTauFH_B_pt", &BsTauTauFH_B_pt );
//    tree_->Branch("BsTauTauFH_B_eta", &BsTauTauFH_B_eta );
//    tree_->Branch("BsTauTauFH_B_phi", &BsTauTauFH_B_phi );
//    tree_->Branch("BsTauTauFH_B_mass", &BsTauTauFH_B_mass );
//    tree_->Branch("BsTauTauFH_B_vprob", &BsTauTauFH_B_vprob );
//    tree_->Branch("BsTauTauFH_B_lip", &BsTauTauFH_B_lip);
//    tree_->Branch("BsTauTauFH_B_lips", &BsTauTauFH_B_lips);
//    tree_->Branch("BsTauTauFH_B_pvip", &BsTauTauFH_B_pvip);
//    tree_->Branch("BsTauTauFH_B_pvips", &BsTauTauFH_B_pvips);
//    tree_->Branch("BsTauTauFH_B_fl3d", &BsTauTauFH_B_fl3d);
//    tree_->Branch("BsTauTauFH_B_fls3d", &BsTauTauFH_B_fls3d);
//    tree_->Branch("BsTauTauFH_B_alpha", &BsTauTauFH_B_alpha);
//    tree_->Branch("BsTauTauFH_B_maxdoca", &BsTauTauFH_B_maxdoca);
//    tree_->Branch("BsTauTauFH_B_mindoca", &BsTauTauFH_B_mindoca);
//    tree_->Branch("BsTauTauFH_B_vx", &BsTauTauFH_B_vx );
//    tree_->Branch("BsTauTauFH_B_vy", &BsTauTauFH_B_vy );
//    tree_->Branch("BsTauTauFH_B_vz", &BsTauTauFH_B_vz );
//    tree_->Branch("BsTauTauFH_B_iso", &BsTauTauFH_B_iso);
//    tree_->Branch("BsTauTauFH_B_iso_ntracks", &BsTauTauFH_B_iso_ntracks );
//    tree_->Branch("BsTauTauFH_B_iso_mindoca", &BsTauTauFH_B_iso_mindoca );
//
//    tree_->Branch("BsTauTauFH_ngenmuons", &BsTauTauFH_ngenmuons);
//    tree_->Branch("BsTauTauFH_isgen3", &BsTauTauFH_isgen3);
//    tree_->Branch("BsTauTauFH_isgen3matched", &BsTauTauFH_isgen3matched);
//    tree_->Branch("BsTauTauFH_nch", &BsTauTauFH_nch);
//    tree_->Branch("BsTauTauFH_nch_after_dnn", &BsTauTauFH_nch_after_dnn);
//    tree_->Branch("BsTauTauFH_nch_before_dnn", &BsTauTauFH_nch_before_dnn);
//    tree_->Branch("BsTauTauFH_nch_qr", &BsTauTauFH_nch_qr);
//    tree_->Branch("BsTauTauFH_ngentau3", &BsTauTauFH_ngentau3); 
//    tree_->Branch("BsTauTauFH_ngentau", &BsTauTauFH_ngentau);
//    tree_->Branch("BsTauTauFH_gentaupt", &BsTauTauFH_gentaupt);
//    tree_->Branch("BsTauTauFH_gentaudm", &BsTauTauFH_gentaudm);
//
//  }
//
//  if (runFlags["doBsTauTauFH_mr"]){
//    tree_->Branch("BsTauTauFH_mr_tau_pi1_pt", &BsTauTauFH_mr_tau_pi1_pt );
//    tree_->Branch("BsTauTauFH_mr_tau_pi1_eta", &BsTauTauFH_mr_tau_pi1_eta );
//    tree_->Branch("BsTauTauFH_mr_tau_pi1_phi", &BsTauTauFH_mr_tau_pi1_phi );
//    tree_->Branch("BsTauTauFH_mr_tau_pi1_mass", &BsTauTauFH_mr_tau_pi1_mass );
//    tree_->Branch("BsTauTauFH_mr_tau_pi2_pt", &BsTauTauFH_mr_tau_pi2_pt );
//    tree_->Branch("BsTauTauFH_mr_tau_pi2_eta", &BsTauTauFH_mr_tau_pi2_eta );
//    tree_->Branch("BsTauTauFH_mr_tau_pi2_phi", &BsTauTauFH_mr_tau_pi2_phi );
//    tree_->Branch("BsTauTauFH_mr_tau_pi2_mass", &BsTauTauFH_mr_tau_pi2_mass );
//    tree_->Branch("BsTauTauFH_mr_tau_pi3_pt", &BsTauTauFH_mr_tau_pi3_pt );
//    tree_->Branch("BsTauTauFH_mr_tau_pi3_eta", &BsTauTauFH_mr_tau_pi3_eta );
//    tree_->Branch("BsTauTauFH_mr_tau_pi3_phi", &BsTauTauFH_mr_tau_pi3_phi );
//    tree_->Branch("BsTauTauFH_mr_tau_pi3_mass", &BsTauTauFH_mr_tau_pi3_mass );
//
//    tree_->Branch("BsTauTauFH_mr_tau_genpt", &BsTauTauFH_mr_tau_genpt );
//    tree_->Branch("BsTauTauFH_mr_tau_geneta", &BsTauTauFH_mr_tau_geneta );
//    tree_->Branch("BsTauTauFH_mr_tau_genphi", &BsTauTauFH_mr_tau_genphi );
//    tree_->Branch("BsTauTauFH_mr_tau_genmass", &BsTauTauFH_mr_tau_genmass );
//
//    tree_->Branch("BsTauTauFH_mr_tau_genpt_bd", &BsTauTauFH_mr_tau_genpt_bd );
//    tree_->Branch("BsTauTauFH_mr_tau_geneta_bd", &BsTauTauFH_mr_tau_geneta_bd );
//    tree_->Branch("BsTauTauFH_mr_tau_genphi_bd", &BsTauTauFH_mr_tau_genphi_bd );
//    tree_->Branch("BsTauTauFH_mr_tau_genmass_bd", &BsTauTauFH_mr_tau_genmass_bd );
//    
//  }
//
//  if (runFlags["doBsDstarTauNu"]){
//
//    tree_->Branch("IsBsDstarTauNu", &IsBsDstarTauNu );
//
//    tree_->Branch("BsDstarTauNu_nCandidates", &BsDstarTauNu_nCandidates );
//
//    tree_->Branch("BsDstarTauNu_mu1_pt", &BsDstarTauNu_mu1_pt );
//    tree_->Branch("BsDstarTauNu_mu1_eta", &BsDstarTauNu_mu1_eta );
//    tree_->Branch("BsDstarTauNu_mu1_phi", &BsDstarTauNu_mu1_phi );
//    tree_->Branch("BsDstarTauNu_mu1_mass", &BsDstarTauNu_mu1_mass );
//    tree_->Branch("BsDstarTauNu_mu1_q", &BsDstarTauNu_mu1_q );
//    tree_->Branch("BsDstarTauNu_mu1_isLoose"  , &BsDstarTauNu_mu1_isLoose   );
//    tree_->Branch("BsDstarTauNu_mu1_isTight"  , &BsDstarTauNu_mu1_isTight   );
//    tree_->Branch("BsDstarTauNu_mu1_isPF"     , &BsDstarTauNu_mu1_isPF      );
//    tree_->Branch("BsDstarTauNu_mu1_isGlobal" , &BsDstarTauNu_mu1_isGlobal  );
//    tree_->Branch("BsDstarTauNu_mu1_isTracker", &BsDstarTauNu_mu1_isTracker );
//    tree_->Branch("BsDstarTauNu_mu1_isSoft"   , &BsDstarTauNu_mu1_isSoft    );
//    tree_->Branch("BsDstarTauNu_mu1_vx"   , &BsDstarTauNu_mu1_vx    );
//    tree_->Branch("BsDstarTauNu_mu1_vy"   , &BsDstarTauNu_mu1_vy    );
//    tree_->Branch("BsDstarTauNu_mu1_vz"   , &BsDstarTauNu_mu1_vz    );
//    tree_->Branch("BsDstarTauNu_mu1_iso"   , &BsDstarTauNu_mu1_iso    );
//    tree_->Branch("BsDstarTauNu_mu1_dbiso"   , &BsDstarTauNu_mu1_dbiso    );
//
//
//    tree_->Branch("BsDstarTauNu_tau_fullfit_pt", &BsDstarTauNu_tau_fullfit_pt );
//    tree_->Branch("BsDstarTauNu_tau_fullfit_eta", &BsDstarTauNu_tau_fullfit_eta );
//    tree_->Branch("BsDstarTauNu_tau_fullfit_phi", &BsDstarTauNu_tau_fullfit_phi );
//    tree_->Branch("BsDstarTauNu_tau_fullfit_mass", &BsDstarTauNu_tau_fullfit_mass );
//    tree_->Branch("BsDstarTauNu_tau_pt", &BsDstarTauNu_tau_pt );
//    tree_->Branch("BsDstarTauNu_tau_eta", &BsDstarTauNu_tau_eta );
//    tree_->Branch("BsDstarTauNu_tau_phi", &BsDstarTauNu_tau_phi );
//    tree_->Branch("BsDstarTauNu_tau_mass", &BsDstarTauNu_tau_mass );
//    tree_->Branch("BsDstarTauNu_tau_rhomass1", &BsDstarTauNu_tau_rhomass1 );
//    tree_->Branch("BsDstarTauNu_tau_rhomass2", &BsDstarTauNu_tau_rhomass2 );
//    tree_->Branch("BsDstarTauNu_tau_q", &BsDstarTauNu_tau_q );
//    tree_->Branch("BsDstarTauNu_tau_vx"   , &BsDstarTauNu_tau_vx    );
//    tree_->Branch("BsDstarTauNu_tau_vy"   , &BsDstarTauNu_tau_vy    );
//    tree_->Branch("BsDstarTauNu_tau_vz"   , &BsDstarTauNu_tau_vz    );
//
//    tree_->Branch("BsDstarTauNu_tau_max_dr_3prong", &BsDstarTauNu_tau_max_dr_3prong);
//    tree_->Branch("BsDstarTauNu_tau_lip", &BsDstarTauNu_tau_lip);
//    tree_->Branch("BsDstarTauNu_tau_lips", &BsDstarTauNu_tau_lips);
//    tree_->Branch("BsDstarTauNu_tau_pvip", &BsDstarTauNu_tau_pvip);
//    tree_->Branch("BsDstarTauNu_tau_pvips", &BsDstarTauNu_tau_pvips);
//    tree_->Branch("BsDstarTauNu_tau_fl3d", &BsDstarTauNu_tau_fl3d);
//    tree_->Branch("BsDstarTauNu_tau_fls3d", &BsDstarTauNu_tau_fls3d);
//    tree_->Branch("BsDstarTauNu_tau_alpha", &BsDstarTauNu_tau_alpha);
//    tree_->Branch("BsDstarTauNu_tau_vprob", &BsDstarTauNu_tau_vprob);
//    tree_->Branch("BsDstarTauNu_tau_isRight", &BsDstarTauNu_tau_isRight);
//    tree_->Branch("BsDstarTauNu_tau_matched_ppdgId", &BsDstarTauNu_tau_matched_ppdgId);
//    tree_->Branch("BsDstarTauNu_tau_matched_gentaupt", &BsDstarTauNu_tau_matched_gentaupt);
//    //    tree_->Branch("BsDstarTauNu_tau_gentaupt", &BsDstarTauNu_tau_gentaupt);
//    tree_->Branch("BsDstarTauNu_tau_sumofdnn", &BsDstarTauNu_tau_sumofdnn);
//    tree_->Branch("BsDstarTauNu_tau_pfidx1", &BsDstarTauNu_tau_pfidx1);
//    tree_->Branch("BsDstarTauNu_tau_pfidx2", &BsDstarTauNu_tau_pfidx2);
//    tree_->Branch("BsDstarTauNu_tau_pfidx3", &BsDstarTauNu_tau_pfidx3);
//
//    tree_->Branch("BsDstarTauNu_tau_pi1_pt", &BsDstarTauNu_tau_pi1_pt );
//    tree_->Branch("BsDstarTauNu_tau_pi1_eta", &BsDstarTauNu_tau_pi1_eta );
//    tree_->Branch("BsDstarTauNu_tau_pi1_phi", &BsDstarTauNu_tau_pi1_phi );
//    tree_->Branch("BsDstarTauNu_tau_pi1_mass", &BsDstarTauNu_tau_pi1_mass );
//
//    tree_->Branch("BsDstarTauNu_tau_pi2_pt", &BsDstarTauNu_tau_pi2_pt );
//    tree_->Branch("BsDstarTauNu_tau_pi2_eta", &BsDstarTauNu_tau_pi2_eta );
//    tree_->Branch("BsDstarTauNu_tau_pi2_phi", &BsDstarTauNu_tau_pi2_phi );
//    tree_->Branch("BsDstarTauNu_tau_pi2_mass", &BsDstarTauNu_tau_pi2_mass );
//
//    tree_->Branch("BsDstarTauNu_tau_pi3_pt", &BsDstarTauNu_tau_pi3_pt );
//    tree_->Branch("BsDstarTauNu_tau_pi3_eta", &BsDstarTauNu_tau_pi3_eta );
//    tree_->Branch("BsDstarTauNu_tau_pi3_phi", &BsDstarTauNu_tau_pi3_phi );
//    tree_->Branch("BsDstarTauNu_tau_pi3_mass", &BsDstarTauNu_tau_pi3_mass );
//
//
//    tree_->Branch("BsDstarTauNu_PV_vx", &BsDstarTauNu_PV_vx );
//    tree_->Branch("BsDstarTauNu_PV_vy", &BsDstarTauNu_PV_vy );
//    tree_->Branch("BsDstarTauNu_PV_vz", &BsDstarTauNu_PV_vz );
//
//    tree_->Branch("BsDstarTauNu_bbPV_vx", &BsDstarTauNu_bbPV_vx );
//    tree_->Branch("BsDstarTauNu_bbPV_vy", &BsDstarTauNu_bbPV_vy );
//    tree_->Branch("BsDstarTauNu_bbPV_vz", &BsDstarTauNu_bbPV_vz );
//
//    tree_->Branch("BsDstarTauNu_bbPV_refit_vx", &BsDstarTauNu_bbPV_vx );
//    tree_->Branch("BsDstarTauNu_bbPV_refit_vy", &BsDstarTauNu_bbPV_vy );
//    tree_->Branch("BsDstarTauNu_bbPV_refit_vz", &BsDstarTauNu_bbPV_vz );
//
//    tree_->Branch("BsDstarTauNu_genPV_vx", &BsDstarTauNu_genPV_vx );
//    tree_->Branch("BsDstarTauNu_genPV_vy", &BsDstarTauNu_genPV_vy );
//    tree_->Branch("BsDstarTauNu_genPV_vz", &BsDstarTauNu_genPV_vz );
//
//    tree_->Branch("BsDstarTauNu_Ds_pt", &BsDstarTauNu_Ds_pt );
//    tree_->Branch("BsDstarTauNu_Ds_eta", &BsDstarTauNu_Ds_eta );
//    tree_->Branch("BsDstarTauNu_Ds_phi", &BsDstarTauNu_Ds_phi );
//    tree_->Branch("BsDstarTauNu_Ds_mass", &BsDstarTauNu_Ds_mass );
//    tree_->Branch("BsDstarTauNu_Ds_vprob", &BsDstarTauNu_Ds_vprob );
//    tree_->Branch("BsDstarTauNu_Ds_lip", &BsDstarTauNu_Ds_lip);
//    tree_->Branch("BsDstarTauNu_Ds_lips", &BsDstarTauNu_Ds_lips);
//    tree_->Branch("BsDstarTauNu_Ds_pvip", &BsDstarTauNu_Ds_pvip);
//    tree_->Branch("BsDstarTauNu_Ds_pvips", &BsDstarTauNu_Ds_pvips);
//    tree_->Branch("BsDstarTauNu_Ds_fl3d", &BsDstarTauNu_Ds_fl3d);
//    tree_->Branch("BsDstarTauNu_Ds_fls3d", &BsDstarTauNu_Ds_fls3d);
//    tree_->Branch("BsDstarTauNu_Ds_alpha", &BsDstarTauNu_Ds_alpha);
//    tree_->Branch("BsDstarTauNu_Ds_vx", &BsDstarTauNu_Ds_vx );
//    tree_->Branch("BsDstarTauNu_Ds_vy", &BsDstarTauNu_Ds_vy );
//    tree_->Branch("BsDstarTauNu_Ds_vz", &BsDstarTauNu_Ds_vz );
//    tree_->Branch("BsDstarTauNu_Ds_unfit_pt", &BsDstarTauNu_Ds_unfit_pt );
//    tree_->Branch("BsDstarTauNu_Ds_unfit_mass", &BsDstarTauNu_Ds_unfit_mass );
//    tree_->Branch("BsDstarTauNu_Ds_ptfrac", &BsDstarTauNu_Ds_ptfrac );
//
//    tree_->Branch("BsDstarTauNu_D0_pt", &BsDstarTauNu_D0_pt );
//    tree_->Branch("BsDstarTauNu_D0_eta", &BsDstarTauNu_D0_eta );
//    tree_->Branch("BsDstarTauNu_D0_phi", &BsDstarTauNu_D0_phi );
//    tree_->Branch("BsDstarTauNu_D0_mass", &BsDstarTauNu_D0_mass );
//    tree_->Branch("BsDstarTauNu_D0_vprob", &BsDstarTauNu_D0_vprob );
//    tree_->Branch("BsDstarTauNu_D0_lip", &BsDstarTauNu_D0_lip);
//    tree_->Branch("BsDstarTauNu_D0_lips", &BsDstarTauNu_D0_lips);
//    tree_->Branch("BsDstarTauNu_D0_pvip", &BsDstarTauNu_D0_pvip);
//    tree_->Branch("BsDstarTauNu_D0_pvips", &BsDstarTauNu_D0_pvips);
//    tree_->Branch("BsDstarTauNu_D0_fl3d", &BsDstarTauNu_D0_fl3d);
//    tree_->Branch("BsDstarTauNu_D0_fls3d", &BsDstarTauNu_D0_fls3d);
//    tree_->Branch("BsDstarTauNu_D0_alpha", &BsDstarTauNu_D0_alpha);
//    tree_->Branch("BsDstarTauNu_D0_vx", &BsDstarTauNu_D0_vx );
//    tree_->Branch("BsDstarTauNu_D0_vy", &BsDstarTauNu_D0_vy );
//    tree_->Branch("BsDstarTauNu_D0_vz", &BsDstarTauNu_D0_vz );
//    tree_->Branch("BsDstarTauNu_D0_unfit_pt", &BsDstarTauNu_D0_unfit_pt );
//    tree_->Branch("BsDstarTauNu_D0_unfit_mass", &BsDstarTauNu_D0_unfit_mass );
//    tree_->Branch("BsDstarTauNu_D0_ptfrac", &BsDstarTauNu_D0_ptfrac );
//
//    tree_->Branch("BsDstarTauNu_k_charge", &BsDstarTauNu_k_charge );
//    tree_->Branch("BsDstarTauNu_pi_charge", &BsDstarTauNu_pi_charge );
//    tree_->Branch("BsDstarTauNu_spi_charge", &BsDstarTauNu_spi_charge );
//
//    tree_->Branch("BsDstarTauNu_B_pt", &BsDstarTauNu_B_pt );
//    tree_->Branch("BsDstarTauNu_B_eta", &BsDstarTauNu_B_eta );
//    tree_->Branch("BsDstarTauNu_B_phi", &BsDstarTauNu_B_phi );
//    tree_->Branch("BsDstarTauNu_B_mass", &BsDstarTauNu_B_mass );
//    tree_->Branch("BsDstarTauNu_B_vprob", &BsDstarTauNu_B_vprob );
//    tree_->Branch("BsDstarTauNu_B_lip", &BsDstarTauNu_B_lip);
//    tree_->Branch("BsDstarTauNu_B_lips", &BsDstarTauNu_B_lips);
//    tree_->Branch("BsDstarTauNu_B_pvip", &BsDstarTauNu_B_pvip);
//    tree_->Branch("BsDstarTauNu_B_pvips", &BsDstarTauNu_B_pvips);
//    tree_->Branch("BsDstarTauNu_B_fl3d", &BsDstarTauNu_B_fl3d);
//    tree_->Branch("BsDstarTauNu_B_fls3d", &BsDstarTauNu_B_fls3d);
//    tree_->Branch("BsDstarTauNu_B_alpha", &BsDstarTauNu_B_alpha);
//    tree_->Branch("BsDstarTauNu_B_maxdoca", &BsDstarTauNu_B_maxdoca);
//    tree_->Branch("BsDstarTauNu_B_mindoca", &BsDstarTauNu_B_mindoca);
//    tree_->Branch("BsDstarTauNu_B_vx", &BsDstarTauNu_B_vx );
//    tree_->Branch("BsDstarTauNu_B_vy", &BsDstarTauNu_B_vy );
//    tree_->Branch("BsDstarTauNu_B_vz", &BsDstarTauNu_B_vz );
//    tree_->Branch("BsDstarTauNu_B_iso", &BsDstarTauNu_B_iso);
//    tree_->Branch("BsDstarTauNu_B_iso_ntracks", &BsDstarTauNu_B_iso_ntracks );
//    tree_->Branch("BsDstarTauNu_B_iso_mindoca", &BsDstarTauNu_B_iso_mindoca );
//    tree_->Branch("BsDstarTauNu_B_mm2", &BsDstarTauNu_B_mm2 );
//    tree_->Branch("BsDstarTauNu_B_q2", &BsDstarTauNu_B_q2 );
//    tree_->Branch("BsDstarTauNu_B_Es", &BsDstarTauNu_B_Es );
//    tree_->Branch("BsDstarTauNu_B_ptback", &BsDstarTauNu_B_ptback );
//
//    tree_->Branch("BsDstarTauNu_ngenmuons", &BsDstarTauNu_ngenmuons);
//    tree_->Branch("BsDstarTauNu_isgen3", &BsDstarTauNu_isgen3);
//    tree_->Branch("BsDstarTauNu_isgen3matched", &BsDstarTauNu_isgen3matched);
//    tree_->Branch("BsDstarTauNu_nch", &BsDstarTauNu_nch);
//    tree_->Branch("BsDstarTauNu_nch_after_dnn", &BsDstarTauNu_nch_after_dnn);
//    tree_->Branch("BsDstarTauNu_nch_before_dnn", &BsDstarTauNu_nch_before_dnn);
//    tree_->Branch("BsDstarTauNu_nch_qr", &BsDstarTauNu_nch_qr);
//    tree_->Branch("BsDstarTauNu_ngentau3", &BsDstarTauNu_ngentau3); 
//    tree_->Branch("BsDstarTauNu_ngentau", &BsDstarTauNu_ngentau);
//    tree_->Branch("BsDstarTauNu_gentaupt", &BsDstarTauNu_gentaupt);
//    tree_->Branch("BsDstarTauNu_gentaudm", &BsDstarTauNu_gentaudm);
//    
//  }

}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  /** genParticle */
  genParticle_N = 0;
  genParticle_pt.clear();
  //  genParticle_px.clear();
  //  genParticle_py.clear();
  //  genParticle_pz.clear();
  //  genParticle_e.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_isPrompt.clear();
  genParticle_isDirectPromptTauDecayProduct.clear();
  genParticle_fromHardProcessFinalState.clear();
  genParticle_isDirectHardProcessTauDecayProductFinalState.clear();
  genParticle_status.clear();
  genParticle_mother.clear();
  genParticle_mother_pt.clear();
  genParticle_nMoth.clear();
  genParticle_nDau.clear();
  genParticle_dau.clear();
  genWeightBkgB = 0;

  //  genParticle_pmother.clear();
  genParticle_pdgs.clear();
  genParticle_layers.clear();
  genParticle_ppt.clear();
  genParticle_peta.clear();
  genParticle_pphi.clear();

  /** generator info */
  genWeight   = 0;
  qScale      = 0;
  genFacWeightUp       = 0;
  genFacWeightDown     = 0;
  genRenWeightUp       = 0;
  genRenWeightDown     = 0;
  genFacRenWeightUp    = 0;
  genFacRenWeightDown  = 0;
  PDF_rms = 0;
  PDF_id.clear();  
  PDF_x.clear();	
  PDF_xPDF.clear();
  lheV_pt = 0;
  lheHT = 0;
  lheNj = 0;
  lheNb = 0;
  lheV_mass = 0;
    
 

  /** HLT trigger decisions */
  HLT_isFired.clear();
  HLT_BPH_isFired.clear();

  /** HLT trigger objects */
  triggerObject_pt.clear();
  triggerObject_eta.clear();
  triggerObject_phi.clear();
  triggerObject_mass.clear();
  triggerObject_lastname.clear();
  triggerObject_filterIDs.clear();
  triggerObject_filterLabels.clear();
  triggerObject_firedTrigger.clear();

  /** HLT filter decisions */
  passFilter_HBHE_                  = false;
  passFilter_HBHELoose_             = false;
  passFilter_HBHETight_             = false;
  passFilter_HBHEIso_               = false;
  passFilter_CSCHalo_               = false;
  passFilter_CSCTightHalo2015_      = false;
  passFilter_HCALlaser_             = false;
  passFilter_ECALDeadCell_          = false;
  passFilter_GoodVtx_               = false;
  passFilter_TrkFailure_            = false;
  passFilter_EEBadSc_               = false;
  passFilter_ECALlaser_             = false;
  passFilter_TrkPOG_                = false;
  passFilter_TrkPOG_manystrip_      = false;
  passFilter_TrkPOG_toomanystrip_   = false;
  passFilter_TrkPOG_logError_       = false;
  passFilter_METFilters_            = false;
  //NEW FOR ICHEP
  passFilter_CSCTightHaloTrkMuUnvetoFilter_   = false;
  passFilter_globalTightHalo2016_             = false;
  passFilter_globalSuperTightHalo2016_             = false;
  passFilter_HcalStripHalo_                   = false;
  passFilter_chargedHadronTrackResolution_    = false;
  passFilter_muonBadTrack_                    = false;
  flag_badMuons_                    = false;
  flag_duplicateMuons_              = false;
  flag_nobadMuons_                  = false;

  /** energy density */
  rho = 0;
  
 

  /** MET */
  METraw_et.clear();	 
  METraw_phi.clear();
  METraw_sumEt.clear();
  MET_corrPx.clear();
  MET_corrPy.clear();
  MET_et.clear();
  MET_phi.clear();
  MET_puppi_et.clear();
  MET_puppi_phi.clear();

  MET_sumEt.clear();
  MET_T1Uncertainty.clear();
  
  MET_JetEnUp.clear();
  MET_JetEnDown.clear();
  MET_JetResUp.clear();
  MET_JetResDown.clear();
  MET_UnclusteredEnUp.clear();
  MET_UnclusteredEnDown.clear();

  /** MET SVift*/
  MET_significance.clear();
  MET_cov00.clear();
  MET_cov10.clear();
  MET_cov11.clear();
  MET_mva_et.clear();
  MET_mva_phi.clear();
  MET_mva_cov00.clear();
  MET_mva_cov10.clear();
  MET_mva_cov11.clear();
  MET_mva_recoil_pt.clear();
  MET_mva_recoil_eta.clear();
  MET_mva_recoil_phi.clear();
  MET_mva_recoil_pdgId.clear();
  MET_Nmva.clear();

  /*------------------------EVENT infos-------------------------*/    
  EVENT_event = 0;
  EVENT_run = 0;
  EVENT_lumiBlock = 0;

  /*--------------------------PV infos--------------------------*/
  PV_N = 0;
  PV_filter = true;
  PV_chi2.clear();
  PV_ndof.clear();
  PV_rho.clear();
  PV_z.clear();
  BeamSpot_x0.clear();
  BeamSpot_y0.clear();
  BeamSpot_z0.clear();
  /*--------------------------PU infos--------------------------*/  			       
  nPuVtxTrue.clear();
  nPuVtx.clear();
  bX.clear();

  /*-------------------------JPSI infos--------------------------*/ 
  IsJpsiMu.clear();
  //  IsBsTauTau.clear();
  //  IsBsTauTauFH.clear();
  //  IsBsDstarTauNu.clear();

  JpsiMu_nCandidates.clear();
  JpsiMu_mu1_pt.clear();
  JpsiMu_mu1_eta.clear();
  JpsiMu_mu1_phi.clear();
  JpsiMu_mu1_mass.clear();
  JpsiMu_mu1_unfit_pt.clear();
  JpsiMu_mu1_unfit_eta.clear();
  JpsiMu_mu1_unfit_phi.clear();
  JpsiMu_mu1_unfit_mass.clear();
  JpsiMu_mu1_q.clear();
  JpsiMu_mu1_isLoose.clear();
  JpsiMu_mu1_isTight.clear();
  JpsiMu_mu1_isPF.clear();
  JpsiMu_mu1_isGlobal.clear();
  JpsiMu_mu1_isTracker.clear();
  JpsiMu_mu1_isSoft.clear();
  JpsiMu_mu1_vx.clear();
  JpsiMu_mu1_vy.clear();
  JpsiMu_mu1_vz.clear();
  JpsiMu_mu1_iso.clear();
  JpsiMu_mu1_dbiso.clear();

  JpsiMu_mu2_pt.clear();
  JpsiMu_mu2_eta.clear();
  JpsiMu_mu2_phi.clear();
  JpsiMu_mu2_mass.clear();
  JpsiMu_mu2_unfit_pt.clear();
  JpsiMu_mu2_unfit_eta.clear();
  JpsiMu_mu2_unfit_phi.clear();
  JpsiMu_mu2_unfit_mass.clear();
  JpsiMu_mu2_q.clear();
  JpsiMu_mu2_isLoose.clear();
  JpsiMu_mu2_isTight.clear();
  JpsiMu_mu2_isPF.clear();
  JpsiMu_mu2_isGlobal.clear();
  JpsiMu_mu2_isTracker.clear();
  JpsiMu_mu2_isSoft.clear();
  JpsiMu_mu2_vx.clear();
  JpsiMu_mu2_vy.clear();
  JpsiMu_mu2_vz.clear();
  JpsiMu_mu2_iso.clear();
  JpsiMu_mu2_dbiso.clear();

  JpsiMu_mu3_pt.clear();
  JpsiMu_mu3_eta.clear();
  JpsiMu_mu3_phi.clear();
  JpsiMu_mu3_mass.clear();
  JpsiMu_mu3_unfit_pt.clear();
  JpsiMu_mu3_unfit_eta.clear();
  JpsiMu_mu3_unfit_phi.clear();
  JpsiMu_mu3_unfit_mass.clear();
  JpsiMu_mu3_doca2mu1.clear();
  JpsiMu_mu3_doca2mu2.clear();
  JpsiMu_mu3_q.clear();
  JpsiMu_mu3_isLoose.clear();
  JpsiMu_mu3_isTight.clear();
  JpsiMu_mu3_isPF.clear();
  JpsiMu_mu3_isGlobal.clear();
  JpsiMu_mu3_isTracker.clear();
  JpsiMu_mu3_isSoft.clear();
  JpsiMu_mu3_vx.clear();
  JpsiMu_mu3_vy.clear();
  JpsiMu_mu3_vz.clear();
  JpsiMu_mu3_iso.clear();
  JpsiMu_mu3_dbiso.clear();

  JpsiMu_Jpsi_pt.clear();
  JpsiMu_Jpsi_eta.clear();
  JpsiMu_Jpsi_phi.clear();
  JpsiMu_Jpsi_mass.clear();
  JpsiMu_Jpsi_vprob.clear();
  JpsiMu_Jpsi_lip.clear();
  JpsiMu_Jpsi_lips.clear();
  JpsiMu_Jpsi_pvip.clear();
  JpsiMu_Jpsi_pvips.clear();
  JpsiMu_Jpsi_fl3d.clear();
  JpsiMu_Jpsi_fls3d.clear();
  JpsiMu_Jpsi_alpha.clear();
  JpsiMu_Jpsi_maxdoca.clear();
  JpsiMu_Jpsi_mindoca.clear();
  JpsiMu_Jpsi_vx.clear();
  JpsiMu_Jpsi_vy.clear();
  JpsiMu_Jpsi_vz.clear();
  JpsiMu_Jpsi_unfit_pt.clear();
  JpsiMu_Jpsi_unfit_mass.clear();
  JpsiMu_Jpsi_unfit_vprob.clear();
  JpsiMu_Jpsi_unfit_vx.clear();
  JpsiMu_Jpsi_unfit_vy.clear();
  JpsiMu_Jpsi_unfit_vz.clear();

  JpsiMu_B_pt.clear();
  JpsiMu_B_eta.clear();
  JpsiMu_B_phi.clear();
  JpsiMu_B_mass.clear();
  JpsiMu_B_mcorr.clear();
  JpsiMu_B_vprob.clear();
  JpsiMu_B_lip.clear();
  JpsiMu_B_lips.clear();
  JpsiMu_B_pvip.clear();
  JpsiMu_B_pvips.clear();
  JpsiMu_B_fl3d.clear();
  JpsiMu_B_fls3d.clear();
  JpsiMu_B_alpha.clear();
  JpsiMu_B_maxdoca.clear();
  JpsiMu_B_mindoca.clear();
  JpsiMu_B_vx.clear();
  JpsiMu_B_vy.clear();
  JpsiMu_B_vz.clear();
  JpsiMu_B_iso.clear();
  JpsiMu_B_iso_ntracks.clear();
  JpsiMu_B_iso_mindoca.clear();
  JpsiMu_B_unfit_pt.clear();
  JpsiMu_B_unfit_mass.clear();
  JpsiMu_B_unfit_vprob.clear();
  JpsiMu_B_unfit_vx.clear();
  JpsiMu_B_unfit_vy.clear();
  JpsiMu_B_unfit_vz.clear();
  JpsiMu_B_q2.clear();
  JpsiMu_B_mm2.clear();
  JpsiMu_B_ptmiss.clear();
  JpsiMu_B_Es.clear();
  JpsiMu_B_ptback.clear();

  JpsiMu_PV_vx.clear();
  JpsiMu_PV_vy.clear();
  JpsiMu_PV_vz.clear();

  JpsiMu_bbPV_vx.clear();
  JpsiMu_bbPV_vy.clear();
  JpsiMu_bbPV_vz.clear();

  //  JpsiMu_bbPV_refit_vx.clear();
  //  JpsiMu_bbPV_refit_vy.clear();
  //  JpsiMu_bbPV_refit_vz.clear();

  JpsiMu_genPV_vx.clear();
  JpsiMu_genPV_vy.clear();
  JpsiMu_genPV_vz.clear();

  JpsiMu_ngenmuons.clear();
  JpsiMu_isgenmatched.clear();
  JpsiMu_mu3_isgenmatched.clear();
  JpsiMu_q2_gen.clear();
  JpsiMu_B_pt_gen.clear();
  JpsiMu_B_eta_gen.clear();
  JpsiMu_B_phi_gen.clear();
  JpsiMu_B_mass_gen.clear();

  JpsiMu_hammer_ff.clear();
  JpsiMu_hammer_ebe.clear();
  JpsiMu_hammer_ebe_toy.clear();
  //  JpsiMu_hammer_ebe_up.clear();
  //  JpsiMu_hammer_ebe_down.clear();
  //  JpsiMu_hammer_ebe_rate_up.clear();
//  JpsiMu_hammer_ebe_rate_down.clear();
//  JpsiMu_hammer_ebe_a0_up.clear();
//  JpsiMu_hammer_ebe_a0_down.clear();
//  JpsiMu_hammer_ebe_a1_up.clear();
//  JpsiMu_hammer_ebe_a1_down.clear();
//  JpsiMu_hammer_ebe_a2_up.clear();
//  JpsiMu_hammer_ebe_a2_down.clear();
//  
//  JpsiMu_hammer_ebe_b0_up.clear();
//  JpsiMu_hammer_ebe_b0_down.clear();
//  JpsiMu_hammer_ebe_b1_up.clear();
//  JpsiMu_hammer_ebe_b1_down.clear();
//  JpsiMu_hammer_ebe_b2_up.clear();
//  JpsiMu_hammer_ebe_b2_down.clear();
//  
//  JpsiMu_hammer_ebe_c1_up.clear();
//  JpsiMu_hammer_ebe_c1_down.clear();
//  JpsiMu_hammer_ebe_c2_up.clear();
//  JpsiMu_hammer_ebe_c2_down.clear();
//  
//  JpsiMu_hammer_ebe_d0_up.clear();
//  JpsiMu_hammer_ebe_d0_down.clear();
//  JpsiMu_hammer_ebe_d1_up.clear();
//  JpsiMu_hammer_ebe_d1_down.clear();
//  JpsiMu_hammer_ebe_d2_up.clear();
//  JpsiMu_hammer_ebe_d2_down.clear();




  JpsiTau_nCandidates = -99;

  JpsiTau_mu1_pt = -99;
  JpsiTau_mu1_eta = -99;
  JpsiTau_mu1_phi = -99;
  JpsiTau_mu1_mass = -99;
  JpsiTau_mu1_q = -99;
  JpsiTau_mu1_isLoose = -99;
  JpsiTau_mu1_isTight = -99;
  JpsiTau_mu1_isPF = -99;
  JpsiTau_mu1_isGlobal = -99;
  JpsiTau_mu1_isTracker = -99;
  JpsiTau_mu1_isSoft = -99;
  JpsiTau_mu1_vx = -99;
  JpsiTau_mu1_vy = -99;
  JpsiTau_mu1_vz = -99;
  //  JpsiTau_mu1_iso = -99;
  JpsiTau_mu1_dbiso = -99;

  JpsiTau_mu2_pt = -99;
  JpsiTau_mu2_eta = -99;
  JpsiTau_mu2_phi = -99;
  JpsiTau_mu2_mass = -99;
  JpsiTau_mu2_q = -99;
  JpsiTau_mu2_isLoose = -99;
  JpsiTau_mu2_isTight = -99;
  JpsiTau_mu2_isPF = -99;
  JpsiTau_mu2_isGlobal = -99;
  JpsiTau_mu2_isTracker = -99;
  JpsiTau_mu2_isSoft = -99;
  JpsiTau_mu2_vx = -99;
  JpsiTau_mu2_vy = -99;
  JpsiTau_mu2_vz = -99;
  //  JpsiTau_mu2_iso = -99;
  JpsiTau_mu2_dbiso = -99;

  JpsiTau_PV_vx = -99;
  JpsiTau_PV_vy = -99;
  JpsiTau_PV_vz = -99;

  JpsiTau_bbPV_vx = -99;
  JpsiTau_bbPV_vy = -99;
  JpsiTau_bbPV_vz = -99;
  JpsiTau_bbPV_chi2 = -99;
  JpsiTau_bbPV_ndof = -99;
  JpsiTau_bbPV_rho = -99;

  JpsiTau_Jpsi_pt = -99;
  JpsiTau_Jpsi_eta = -99;
  JpsiTau_Jpsi_phi = -99;
  JpsiTau_Jpsi_mass = -99;
  JpsiTau_Jpsi_vprob = -99;
  JpsiTau_Jpsi_lip = -99;
  JpsiTau_Jpsi_lips = -99;
  JpsiTau_Jpsi_pvip = -99;
  JpsiTau_Jpsi_pvips = -99;
  JpsiTau_Jpsi_fl3d = -99;
  JpsiTau_Jpsi_fls3d = -99;
  JpsiTau_Jpsi_alpha = -99;
  JpsiTau_Jpsi_maxdoca = -99;
  JpsiTau_Jpsi_mindoca = -99;
  JpsiTau_Jpsi_vx = -99;
  JpsiTau_Jpsi_vy = -99;
  JpsiTau_Jpsi_vz = -99;
//  JpsiTau_Jpsi_unfit_pt = -99;
//  JpsiTau_Jpsi_unfit_mass = -99;
//  JpsiTau_Jpsi_unfit_vprob = -99;
//  JpsiTau_Jpsi_unfit_vx = -99;
//  JpsiTau_Jpsi_unfit_vy = -99;
//  JpsiTau_Jpsi_unfit_vz = -99;



  //  JpsiTau_tau_fullfit_pt.clear();
  //  JpsiTau_tau_fullfit_eta.clear();
  //  JpsiTau_tau_fullfit_phi.clear();
  //  JpsiTau_tau_fullfit_mass.clear();
  JpsiTau_tau_pt.clear();
  JpsiTau_tau_eta.clear();
  JpsiTau_tau_phi.clear();
  JpsiTau_tau_mass.clear();
  JpsiTau_tau_q.clear();
  JpsiTau_tau_vx.clear();
  JpsiTau_tau_vy.clear();
  JpsiTau_tau_vz.clear();


  JpsiTau_tau_max_dr_3prong.clear();
  JpsiTau_tau_lip.clear();
  JpsiTau_tau_lips.clear();
  JpsiTau_tau_pvip.clear();
  JpsiTau_tau_pvips.clear();
  JpsiTau_tau_fl3d.clear();
  JpsiTau_tau_fls3d.clear();
  JpsiTau_tau_alpha.clear();
  JpsiTau_tau_vprob.clear();

  JpsiTau_tau_fl3d_wjpsi.clear();
  JpsiTau_tau_fls3d_wjpsi.clear();


  //  JpsiTau_tau_matched_isBdecay.clear();
  //  JpsiTau_tau_matched_gentaupt.clear();
  JpsiTau_tau_sumofdnn.clear();
  JpsiTau_tau_sumofdnn_1prong.clear();
  JpsiTau_tau_sumofdnn_otherB.clear();
  JpsiTau_tau_sumofdnn_pu.clear();
  //  JpsiTau_tau_sumofdnn_old.clear();
  //  JpsiTau_tau_sumofdnn_others.clear();
  //  JpsiTau_tau_pi1_dnn.clear();
  //  JpsiTau_tau_pi2_dnn.clear();
  //  JpsiTau_tau_pi3_dnn.clear();

  JpsiTau_tau_rhomass1.clear();
  JpsiTau_tau_rhomass2.clear();

  JpsiTau_tau_pi1_pt.clear();
  JpsiTau_tau_pi1_eta.clear();
  JpsiTau_tau_pi1_phi.clear();
  JpsiTau_tau_pi1_mass.clear();
  JpsiTau_tau_pi1_q.clear();  
  JpsiTau_tau_pi1_doca3d.clear();
  JpsiTau_tau_pi1_doca3de.clear();
  JpsiTau_tau_pi1_doca2d.clear();
  JpsiTau_tau_pi1_doca2de.clear();
  //  JpsiTau_tau_pi1_isRight.clear();
  JpsiTau_tau_pi1_dz.clear();
  JpsiTau_tau_pi1_near_dz.clear();
  JpsiTau_tau_pi1_isAssociate.clear();
  JpsiTau_tau_pi1_pvAssociationQuality.clear();
  JpsiTau_tau_pi1_isBdecay.clear();
  JpsiTau_tau_pi1_isBdecaypdg.clear();
  JpsiTau_tau_pi1_isBdecayppdg.clear();
  JpsiTau_tau_pi1_isSignal.clear();
  JpsiTau_tau_pi1_nprong.clear();
  JpsiTau_tau_pi1_nprong_pi0.clear();
  JpsiTau_tau_pi1_dnn.clear();
  JpsiTau_tau_pi1_dnn_1prong.clear();
  JpsiTau_tau_pi1_dnn_otherB.clear();
  JpsiTau_tau_pi1_dnn_pu.clear();
  //  JpsiTau_tau_pi1_dnn_old.clear();
  JpsiTau_tau_pi1_trigMatch.clear();
  JpsiTau_tau_pi1_trigMatch_dr.clear();




  JpsiTau_tau_pi2_pt.clear();
  JpsiTau_tau_pi2_eta.clear();
  JpsiTau_tau_pi2_phi.clear();
  JpsiTau_tau_pi2_mass.clear();
  JpsiTau_tau_pi2_q.clear();
  JpsiTau_tau_pi2_doca3d.clear();
  JpsiTau_tau_pi2_doca3de.clear();
  JpsiTau_tau_pi2_doca2d.clear();
  JpsiTau_tau_pi2_doca2de.clear();
  //  JpsiTau_tau_pi2_isRight.clear();
  JpsiTau_tau_pi2_dz.clear();
  JpsiTau_tau_pi2_near_dz.clear();
  JpsiTau_tau_pi2_isAssociate.clear();
  JpsiTau_tau_pi2_pvAssociationQuality.clear();
  JpsiTau_tau_pi2_isBdecay.clear();
  JpsiTau_tau_pi2_isBdecaypdg.clear();
  JpsiTau_tau_pi2_isBdecayppdg.clear();
  JpsiTau_tau_pi2_isSignal.clear();
  JpsiTau_tau_pi2_nprong.clear();
  JpsiTau_tau_pi2_nprong_pi0.clear();
  JpsiTau_tau_pi2_dnn.clear();
  JpsiTau_tau_pi2_dnn_1prong.clear();
  JpsiTau_tau_pi2_dnn_otherB.clear();
  JpsiTau_tau_pi2_dnn_pu.clear();
  //  JpsiTau_tau_pi2_dnn_old.clear();
  JpsiTau_tau_pi2_trigMatch.clear();
  JpsiTau_tau_pi2_trigMatch_dr.clear();


  JpsiTau_tau_pi3_pt.clear();
  JpsiTau_tau_pi3_eta.clear();
  JpsiTau_tau_pi3_phi.clear();
  JpsiTau_tau_pi3_mass.clear();
  JpsiTau_tau_pi3_q.clear();
  JpsiTau_tau_pi3_doca3d.clear();
  JpsiTau_tau_pi3_doca3de.clear();
  JpsiTau_tau_pi3_doca2d.clear();
  JpsiTau_tau_pi3_doca2de.clear();
  //  JpsiTau_tau_pi3_isRight.clear();
  JpsiTau_tau_pi3_dz.clear();
  JpsiTau_tau_pi3_near_dz.clear();
  JpsiTau_tau_pi3_isAssociate.clear();
  JpsiTau_tau_pi3_pvAssociationQuality.clear();
  JpsiTau_tau_pi3_isBdecay.clear();
  JpsiTau_tau_pi3_isBdecaypdg.clear();
  JpsiTau_tau_pi3_isBdecayppdg.clear();
  JpsiTau_tau_pi3_isSignal.clear();
  JpsiTau_tau_pi3_nprong.clear();
  JpsiTau_tau_pi3_nprong_pi0.clear();
  JpsiTau_tau_pi3_dnn.clear();
  JpsiTau_tau_pi3_dnn_1prong.clear();
  JpsiTau_tau_pi3_dnn_otherB.clear();
  JpsiTau_tau_pi3_dnn_pu.clear();
  //  JpsiTau_tau_pi3_dnn_old.clear();
  JpsiTau_tau_pi3_trigMatch.clear();
  JpsiTau_tau_pi3_trigMatch_dr.clear();


  JpsiTau_tau_delta_chi2.clear();
  JpsiTau_tau_delta_n_ch.clear();
  JpsiTau_tau_delta_n_mu.clear();
  JpsiTau_tau_vweight.clear();
  JpsiTau_tau_refit_vx.clear();
  JpsiTau_tau_refit_vy.clear();
  JpsiTau_tau_refit_vz.clear();
  JpsiTau_tau_refit_chi2.clear();
  JpsiTau_tau_refit_ndof.clear();
  JpsiTau_tau_refit_rho.clear();

  JpsiTau_tau_iso.clear();
  JpsiTau_tau_iso_ntracks.clear();
  JpsiTau_tau_iso_mindoca.clear();


  JpsiTau_ptbal.clear();
  JpsiTau_jpsi_tau_alpha.clear();


  JpsiTau_B_pt.clear();
  JpsiTau_B_eta.clear();
  JpsiTau_B_phi.clear();
  JpsiTau_B_mass.clear();
  JpsiTau_B_mcorr.clear();
  JpsiTau_B_vprob.clear();
  JpsiTau_B_lip.clear();
  JpsiTau_B_lips.clear();
  JpsiTau_B_pvip.clear();
  JpsiTau_B_pvips.clear();
  JpsiTau_B_fl3d.clear();
  JpsiTau_B_fls3d.clear();
  JpsiTau_B_alpha.clear();
  JpsiTau_B_maxdoca.clear();
  JpsiTau_B_mindoca.clear();
  JpsiTau_B_vx.clear();
  JpsiTau_B_vy.clear();
  JpsiTau_B_vz.clear();

  //  JpsiTau_B_iso_nocut.clear();
  //  JpsiTau_B_iso_ntracks_nocut.clear();
  //  JpsiTau_B_iso_mindoca_nocut.clear();

  JpsiTau_B_q2.clear();
  JpsiTau_B_mm2.clear();
  JpsiTau_B_ptmiss.clear();
  JpsiTau_B_Es.clear();
  JpsiTau_B_ptback.clear();

  JpsiTau_B_pt_simple.clear();
  JpsiTau_B_eta_simple.clear();
  JpsiTau_B_phi_simple.clear();
  JpsiTau_B_mass_simple.clear();
  JpsiTau_B_q2_simple.clear();
  JpsiTau_B_mm2_simple.clear();
  JpsiTau_B_ptmiss_simple.clear();
  JpsiTau_B_Es_simple.clear();
  JpsiTau_B_ptback_simple.clear();




  //  JpsiTau_bbPV_refit_vx.clear();
  //  JpsiTau_bbPV_refit_vy.clear();
  //  JpsiTau_bbPV_refit_vz.clear();

  JpsiTau_genPV_vx = -99;
  JpsiTau_genPV_vy = -99;
  JpsiTau_genPV_vz = -99;
  JpsiTau_genSV_vx = -99;
  JpsiTau_genSV_vy = -99;
  JpsiTau_genSV_vz = -99;

  JpsiTau_ngenmuons = -99;
  JpsiTau_isgenmatched = -99;

//  JpsiTau_isgen3 = -99;
//  JpsiTau_isgen3matched = -99;

  JpsiTau_nch = -99;
  //  JpsiTau_nch_after_dnn = -99;
  JpsiTau_nch_before = -99;
  //  JpsiTau_nch_qr = -99;

//  JpsiTau_ngentau3 = -99;
//  JpsiTau_ngentau = -99;
//  JpsiTau_gentaupt = -99;
//  JpsiTau_gentaueta = -99;
//  JpsiTau_gentauphi = -99;
//  JpsiTau_gentaumass = -99;
//  JpsiTau_gentaupt_bd = -99;
//  JpsiTau_gentaueta_bd = -99;
//  JpsiTau_gentauphi_bd = -99;
//  JpsiTau_gentaumass_bd = -99;
//  JpsiTau_gentaudm = -99;
  JpsiTau_q2_gen = -99;
  JpsiTau_B_pt_gen = -99;
  JpsiTau_B_eta_gen = -99;
  JpsiTau_B_phi_gen = -99;
  JpsiTau_B_mass_gen = -99;

  JpsiTau_st_idx.clear();

  JpsiTau_st_doca3d.clear();
  JpsiTau_st_doca2d.clear();
  JpsiTau_st_doca3ds.clear();
  JpsiTau_st_doca2ds.clear();
  JpsiTau_st_doca3de.clear();
  JpsiTau_st_doca2de.clear();
  //  JpsiTau_st_isRight.clear();
  JpsiTau_st_isBdecay.clear();
  JpsiTau_st_isBdecaypdg.clear();
  JpsiTau_st_isBdecayppdg.clear();
  JpsiTau_st_isSignal.clear();
  JpsiTau_st_nprong.clear();
  JpsiTau_st_nprong_pi0.clear();

  JpsiTau_st_dz.clear();
  JpsiTau_st_isAssociate.clear();
  JpsiTau_st_near_dz.clear();
  JpsiTau_st_dr_jpsi.clear();
  JpsiTau_st_trigMatch.clear();
  JpsiTau_st_trigMatch_dr.clear();
  JpsiTau_st_pvAssociationQuality.clear();
  JpsiTau_st_pt.clear();
  JpsiTau_st_eta.clear();
  JpsiTau_st_phi.clear();
  JpsiTau_st_charge.clear();
  JpsiTau_st_mass.clear();
  JpsiTau_st_dnn.clear();
  JpsiTau_st_dnn_1prong.clear();
  JpsiTau_st_dnn_otherB.clear();
  JpsiTau_st_dnn_pu.clear();
  //  JpsiTau_st_dnn_old.clear();
  JpsiTau_st_matchidx.clear();

  //  JpsiTau_st_doca3d_max.clear();
  //  JpsiTau_st_doca2d_max.clear();
  //  JpsiTau_st_doca3ds_max.clear();
  //  JpsiTau_st_doca2ds_max.clear();
  //  JpsiTau_st_doca3de_max.clear();
  //  JpsiTau_st_doca2de_max.clear();
  //  JpsiTau_st_dz_max.clear();

  JpsiTau_gen_pion_pt.clear();
  JpsiTau_gen_pion_eta.clear();
  JpsiTau_gen_pion_phi.clear();
  JpsiTau_gen_pion_matched.clear();
  
  JpsiTau_gen_tau_pt.clear();
  JpsiTau_gen_tau_eta.clear();
  JpsiTau_gen_tau_phi.clear();
  JpsiTau_gen_tau_nprong.clear();
  JpsiTau_gen_tau_nmatched.clear();



  JpsiTau_st_nch = 0;
  JpsiTau_st_nch_matched = 0;
  //  JpsiTau_st_decayid = -1;
  //  JpsiTau_st_npi0 = -1;

  JpsiTau_st_n_charged_pions = -99;
  JpsiTau_st_n_neutral_pions = -99;
  JpsiTau_st_n_mu_decay = -99;
  JpsiTau_st_n_e_decay = -99;
  JpsiTau_st_n_occurance = -99;
  JpsiTau_st_decayid = -99;
  JpsiTau_st_gentau_pt = -99;
  JpsiTau_st_gentau_eta = -99;
  JpsiTau_st_gentau_phi = -99;
  JpsiTau_st_genjpsi_pt = -99;
  JpsiTau_st_genjpsi_eta = -99;
  JpsiTau_st_genjpsi_phi = -99;

  //  JpsiTau_perEVT_old = -99;

  JpsiTau_perEVT_otherB = -99;
  JpsiTau_perEVT_sig = -99;
  JpsiTau_perEVT_leptonic = -99;
  JpsiTau_perEVT_1prong = -99;


//  JpsiTau_ed_pfeta.clear();
//  JpsiTau_ed_pfphi.clear();
//  JpsiTau_ed_isRight.clear();
//  JpsiTau_ed_id.clear();
//  JpsiTau_ed_pfdnn.clear();
//  JpsiTau_ed_genpt.clear();

  JpsiTau_hammer_ff.clear();
  JpsiTau_hammer_ebe.clear();
  JpsiTau_hammer_ebe_toy.clear();
  //  JpsiTau_hammer_ebe_up.clear();
  //  JpsiTau_hammer_ebe_down.clear();
  //  JpsiTau_hammer_ebe_rate_up.clear();
  //  JpsiTau_hammer_ebe_rate_down.clear();
//  JpsiTau_hammer_ebe_a0_up.clear();
//  JpsiTau_hammer_ebe_a0_down.clear();
//  JpsiTau_hammer_ebe_a1_up.clear();
//  JpsiTau_hammer_ebe_a1_down.clear();
//  JpsiTau_hammer_ebe_a2_up.clear();
//  JpsiTau_hammer_ebe_a2_down.clear();
//  
//  JpsiTau_hammer_ebe_b0_up.clear();
//  JpsiTau_hammer_ebe_b0_down.clear();
//  JpsiTau_hammer_ebe_b1_up.clear();
//  JpsiTau_hammer_ebe_b1_down.clear();
//  JpsiTau_hammer_ebe_b2_up.clear();
//  JpsiTau_hammer_ebe_b2_down.clear();
//  
//  JpsiTau_hammer_ebe_c1_up.clear();
//  JpsiTau_hammer_ebe_c1_down.clear();
//  JpsiTau_hammer_ebe_c2_up.clear();
//  JpsiTau_hammer_ebe_c2_down.clear();
//  
//  JpsiTau_hammer_ebe_d0_up.clear();
//  JpsiTau_hammer_ebe_d0_down.clear();
//  JpsiTau_hammer_ebe_d1_up.clear();
//  JpsiTau_hammer_ebe_d1_down.clear();
//  JpsiTau_hammer_ebe_d2_up.clear();
//  JpsiTau_hammer_ebe_d2_down.clear();


  /////////////////


//  BsTauTau_nCandidates.clear();
//
//  BsTauTau_mu1_pt.clear();
//  BsTauTau_mu1_eta.clear();
//  BsTauTau_mu1_phi.clear();
//  BsTauTau_mu1_mass.clear();
//  BsTauTau_mu1_q.clear();
//  BsTauTau_mu1_isLoose.clear();
//  BsTauTau_mu1_isTight.clear();
//  BsTauTau_mu1_isPF.clear();
//  BsTauTau_mu1_isGlobal.clear();
//  BsTauTau_mu1_isTracker.clear();
//  BsTauTau_mu1_isSoft.clear();
//  BsTauTau_mu1_vx.clear();
//  BsTauTau_mu1_vy.clear();
//  BsTauTau_mu1_vz.clear();
//  BsTauTau_mu1_iso.clear();
//  BsTauTau_mu1_dbiso.clear();
//
//  BsTauTau_tau_pt.clear();
//  BsTauTau_tau_eta.clear();
//  BsTauTau_tau_phi.clear();
//  BsTauTau_tau_mass.clear();
//  BsTauTau_tau_rhomass1.clear();
//  BsTauTau_tau_rhomass2.clear();
//  BsTauTau_tau_q.clear();
//  BsTauTau_tau_vx.clear();
//  BsTauTau_tau_vy.clear();
//  BsTauTau_tau_vz.clear();
//
//
//  BsTauTau_tau_max_dr_3prong.clear();
//  BsTauTau_tau_lip.clear();
//  BsTauTau_tau_lips.clear();
//  BsTauTau_tau_pvip.clear();
//  BsTauTau_tau_pvips.clear();
//  BsTauTau_tau_fl3d.clear();
//  BsTauTau_tau_fls3d.clear();
//  BsTauTau_tau_alpha.clear();
//  BsTauTau_tau_vprob.clear();
//  BsTauTau_tau_isRight.clear();
//  BsTauTau_tau_isRight1.clear();
//  BsTauTau_tau_isRight2.clear();
//  BsTauTau_tau_isRight3.clear();
////  BsTauTau_tau_dr1.clear();
////  BsTauTau_tau_dr2.clear();
////  BsTauTau_tau_dr3.clear();
////  BsTauTau_tau_ptres1.clear();
////  BsTauTau_tau_ptres2.clear();
////  BsTauTau_tau_ptres3.clear();
//  BsTauTau_tau_matched_ppdgId.clear();
//  BsTauTau_tau_matched_gentaupt.clear();
//  BsTauTau_tau_sumofdnn.clear();
//  BsTauTau_tau_pfidx1.clear();
//  BsTauTau_tau_pfidx2.clear();
//  BsTauTau_tau_pfidx3.clear();
//  BsTauTau_tau_pi1_dnn.clear();
//  BsTauTau_tau_pi2_dnn.clear();
//  BsTauTau_tau_pi3_dnn.clear();
//
//  BsTauTau_tau_pi1_pt.clear();
//  BsTauTau_tau_pi1_eta.clear();
//  BsTauTau_tau_pi1_phi.clear();
//  BsTauTau_tau_pi1_mass.clear();
//  BsTauTau_tau_pi2_pt.clear();
//  BsTauTau_tau_pi2_eta.clear();
//  BsTauTau_tau_pi2_phi.clear();
//  BsTauTau_tau_pi2_mass.clear();
//  BsTauTau_tau_pi3_pt.clear();
//  BsTauTau_tau_pi3_eta.clear();
//  BsTauTau_tau_pi3_phi.clear();
//  BsTauTau_tau_pi3_mass.clear();
//
//
//  BsTauTau_B_pt.clear();
//  BsTauTau_B_eta.clear();
//  BsTauTau_B_phi.clear();
//  BsTauTau_B_mass.clear();
//  BsTauTau_B_vprob.clear();
//  BsTauTau_B_lip.clear();
//  BsTauTau_B_lips.clear();
//  BsTauTau_B_pvip.clear();
//  BsTauTau_B_pvips.clear();
//  BsTauTau_B_fl3d.clear();
//  BsTauTau_B_fls3d.clear();
//  BsTauTau_B_alpha.clear();
//  BsTauTau_B_maxdoca.clear();
//  BsTauTau_B_mindoca.clear();
//  BsTauTau_B_vx.clear();
//  BsTauTau_B_vy.clear();
//  BsTauTau_B_vz.clear();
//  BsTauTau_B_iso.clear();
//  BsTauTau_B_iso_ntracks.clear();
//  BsTauTau_B_iso_mindoca.clear();
//
//  BsTauTau_PV_vx.clear();
//  BsTauTau_PV_vy.clear();
//  BsTauTau_PV_vz.clear();
//
//  BsTauTau_bbPV_vx.clear();
//  BsTauTau_bbPV_vy.clear();
//  BsTauTau_bbPV_vz.clear();
//
//  BsTauTau_bbPV_refit_vx.clear();
//  BsTauTau_bbPV_refit_vy.clear();
//  BsTauTau_bbPV_refit_vz.clear();
//
//  BsTauTau_genPV_vx.clear();
//  BsTauTau_genPV_vy.clear();
//  BsTauTau_genPV_vz.clear();
//
//  BsTauTau_ngenmuons.clear();
//
//  BsTauTau_isgen3.clear();
//  BsTauTau_isgen3matched.clear();
//  BsTauTau_nch.clear();
//  BsTauTau_nch_after_dnn.clear();
//  BsTauTau_nch_before_dnn.clear();
//  BsTauTau_nch_qr.clear();
//  BsTauTau_ngentau3.clear();
//  BsTauTau_ngentau.clear();
//  BsTauTau_gentaupt.clear();
//  BsTauTau_gentaudm.clear();
//
//
//
//  ///////////////////////////////
//
//
//
//  BsTauTauFH_nCandidates.clear();
//  BsTauTauFH_ntaus.clear();
//
//  BsTauTauFH_mu1_pt.clear();
//  BsTauTauFH_mu1_eta.clear();
//  BsTauTauFH_mu1_phi.clear();
//  BsTauTauFH_mu1_mass.clear();
//  BsTauTauFH_mu1_q.clear();
//  BsTauTauFH_mu1_isLoose.clear();
//  BsTauTauFH_mu1_isTight.clear();
//  BsTauTauFH_mu1_isPF.clear();
//  BsTauTauFH_mu1_isGlobal.clear();
//  BsTauTauFH_mu1_isTracker.clear();
//  BsTauTauFH_mu1_isSoft.clear();
//  BsTauTauFH_mu1_vx.clear();
//  BsTauTauFH_mu1_vy.clear();
//  BsTauTauFH_mu1_vz.clear();
//  BsTauTauFH_mu1_iso.clear();
//  BsTauTauFH_mu1_dbiso.clear();
//
//
//  BsTauTauFH_tau1_pt.clear();
//  BsTauTauFH_tau1_eta.clear();
//  BsTauTauFH_tau1_phi.clear();
//  BsTauTauFH_tau1_mass.clear();
//  BsTauTauFH_tau1_rhomass1.clear();
//  BsTauTauFH_tau1_rhomass2.clear();
//  BsTauTauFH_tau1_q.clear();
//  BsTauTauFH_tau1_vx.clear();
//  BsTauTauFH_tau1_vy.clear();
//  BsTauTauFH_tau1_vz.clear();
//
//
//  BsTauTauFH_tau1_max_dr_3prong.clear();
//  BsTauTauFH_tau1_lip.clear();
//  BsTauTauFH_tau1_lips.clear();
//  BsTauTauFH_tau1_pvip.clear();
//  BsTauTauFH_tau1_pvips.clear();
//  BsTauTauFH_tau1_fl3d.clear();
//  BsTauTauFH_tau1_fls3d.clear();
//  BsTauTauFH_tau1_alpha.clear();
//  BsTauTauFH_tau1_vprob.clear();
//  BsTauTauFH_tau1_isRight.clear();
//  BsTauTauFH_tau1_matched_ppdgId.clear();
//  BsTauTauFH_tau1_matched_gentaupt.clear();
//  BsTauTauFH_tau1_sumofdnn.clear();
//  BsTauTauFH_tau1_pfidx1.clear();
//  BsTauTauFH_tau1_pfidx2.clear();
//  BsTauTauFH_tau1_pfidx3.clear();
//  BsTauTauFH_tau1_pi1_dnn.clear();
//  BsTauTauFH_tau1_pi2_dnn.clear();
//  BsTauTauFH_tau1_pi3_dnn.clear();
//
//  BsTauTauFH_tau1_pi1_pt.clear();
//  BsTauTauFH_tau1_pi1_eta.clear();
//  BsTauTauFH_tau1_pi1_phi.clear();
//  BsTauTauFH_tau1_pi1_mass.clear();
//
//  BsTauTauFH_tau1_pi2_pt.clear();
//  BsTauTauFH_tau1_pi2_eta.clear();
//  BsTauTauFH_tau1_pi2_phi.clear();
//  BsTauTauFH_tau1_pi2_mass.clear();
//
//  BsTauTauFH_tau1_pi3_pt.clear();
//  BsTauTauFH_tau1_pi3_eta.clear();
//  BsTauTauFH_tau1_pi3_phi.clear();
//  BsTauTauFH_tau1_pi3_mass.clear();
//
//
//  BsTauTauFH_tau2_pt.clear();
//  BsTauTauFH_tau2_eta.clear();
//  BsTauTauFH_tau2_phi.clear();
//  BsTauTauFH_tau2_mass.clear();
//  BsTauTauFH_tau2_rhomass1.clear();
//  BsTauTauFH_tau2_rhomass2.clear();
//  BsTauTauFH_tau2_q.clear();
//  BsTauTauFH_tau2_vx.clear();
//  BsTauTauFH_tau2_vy.clear();
//  BsTauTauFH_tau2_vz.clear();
//
//
//  BsTauTauFH_tau2_max_dr_3prong.clear();
//  BsTauTauFH_tau2_lip.clear();
//  BsTauTauFH_tau2_lips.clear();
//  BsTauTauFH_tau2_pvip.clear();
//  BsTauTauFH_tau2_pvips.clear();
//  BsTauTauFH_tau2_fl3d.clear();
//  BsTauTauFH_tau2_fls3d.clear();
//  BsTauTauFH_tau2_alpha.clear();
//  BsTauTauFH_tau2_vprob.clear();
//  BsTauTauFH_tau2_isRight.clear();
//  BsTauTauFH_tau2_matched_ppdgId.clear();
//  BsTauTauFH_tau2_matched_gentaupt.clear();
//  BsTauTauFH_tau2_sumofdnn.clear();
//  BsTauTauFH_tau2_pfidx1.clear();
//  BsTauTauFH_tau2_pfidx2.clear();
//  BsTauTauFH_tau2_pfidx3.clear();
//  BsTauTauFH_tau2_pi1_dnn.clear();
//  BsTauTauFH_tau2_pi2_dnn.clear();
//  BsTauTauFH_tau2_pi3_dnn.clear();
//
//  BsTauTauFH_tau2_pi1_pt.clear();
//  BsTauTauFH_tau2_pi1_eta.clear();
//  BsTauTauFH_tau2_pi1_phi.clear();
//  BsTauTauFH_tau2_pi1_mass.clear();
//
//  BsTauTauFH_tau2_pi2_pt.clear();
//  BsTauTauFH_tau2_pi2_eta.clear();
//  BsTauTauFH_tau2_pi2_phi.clear();
//  BsTauTauFH_tau2_pi2_mass.clear();
//
//  BsTauTauFH_tau2_pi3_pt.clear();
//  BsTauTauFH_tau2_pi3_eta.clear();
//  BsTauTauFH_tau2_pi3_phi.clear();
//  BsTauTauFH_tau2_pi3_mass.clear();
//
//  BsTauTauFH_B_pt.clear();
//  BsTauTauFH_B_eta.clear();
//  BsTauTauFH_B_phi.clear();
//  BsTauTauFH_B_mass.clear();
//  BsTauTauFH_B_vprob.clear();
//  BsTauTauFH_B_lip.clear();
//  BsTauTauFH_B_lips.clear();
//  BsTauTauFH_B_pvip.clear();
//  BsTauTauFH_B_pvips.clear();
//  BsTauTauFH_B_fl3d.clear();
//  BsTauTauFH_B_fls3d.clear();
//  BsTauTauFH_B_alpha.clear();
//  BsTauTauFH_B_maxdoca.clear();
//  BsTauTauFH_B_mindoca.clear();
//  BsTauTauFH_B_vx.clear();
//  BsTauTauFH_B_vy.clear();
//  BsTauTauFH_B_vz.clear();
//  BsTauTauFH_B_iso.clear();
//  BsTauTauFH_B_iso_ntracks.clear();
//  BsTauTauFH_B_iso_mindoca.clear();
//
//  BsTauTauFH_PV_vx.clear();
//  BsTauTauFH_PV_vy.clear();
//  BsTauTauFH_PV_vz.clear();
//
//  BsTauTauFH_bbPV_vx.clear();
//  BsTauTauFH_bbPV_vy.clear();
//  BsTauTauFH_bbPV_vz.clear();
//
//  BsTauTauFH_bbPV_refit_vx.clear();
//  BsTauTauFH_bbPV_refit_vy.clear();
//  BsTauTauFH_bbPV_refit_vz.clear();
//
//  BsTauTauFH_genPV_vx.clear();
//  BsTauTauFH_genPV_vy.clear();
//  BsTauTauFH_genPV_vz.clear();
//
//  BsTauTauFH_ngenmuons.clear();
//
//  BsTauTauFH_isgen3.clear();
//  BsTauTauFH_isgen3matched.clear();
//  BsTauTauFH_nch.clear();
//  BsTauTauFH_nch_after_dnn.clear();
//  BsTauTauFH_nch_before_dnn.clear();
//  BsTauTauFH_nch_qr.clear();
//  BsTauTauFH_ngentau3.clear();
//  BsTauTauFH_ngentau.clear();
//  BsTauTauFH_gentaupt.clear();
//  BsTauTauFH_gentaudm.clear();
//
//
//  //////////////////////////////
//
//
//  BsTauTauFH_mr_tau_pi1_pt.clear();
//  BsTauTauFH_mr_tau_pi1_eta.clear();
//  BsTauTauFH_mr_tau_pi1_phi.clear();
//  BsTauTauFH_mr_tau_pi1_mass.clear();
//  BsTauTauFH_mr_tau_pi2_pt.clear();
//  BsTauTauFH_mr_tau_pi2_eta.clear();
//  BsTauTauFH_mr_tau_pi2_phi.clear();
//  BsTauTauFH_mr_tau_pi2_mass.clear();
//  BsTauTauFH_mr_tau_pi3_pt.clear();
//  BsTauTauFH_mr_tau_pi3_eta.clear();
//  BsTauTauFH_mr_tau_pi3_phi.clear();
//  BsTauTauFH_mr_tau_pi3_mass.clear();
//  
//  BsTauTauFH_mr_tau_genpt.clear();
//  BsTauTauFH_mr_tau_geneta.clear();
//  BsTauTauFH_mr_tau_genphi.clear();
//  BsTauTauFH_mr_tau_genmass.clear();
//
//  BsTauTauFH_mr_tau_genpt_bd.clear();
//  BsTauTauFH_mr_tau_geneta_bd.clear();
//  BsTauTauFH_mr_tau_genphi_bd.clear();
//  BsTauTauFH_mr_tau_genmass_bd.clear();
//  
//
//
//  BsDstarTauNu_nCandidates.clear();
//
//
//  BsDstarTauNu_mu1_pt.clear();
//  BsDstarTauNu_mu1_eta.clear();
//  BsDstarTauNu_mu1_phi.clear();
//  BsDstarTauNu_mu1_mass.clear();
//  BsDstarTauNu_mu1_q.clear();
//  BsDstarTauNu_mu1_isLoose.clear();
//  BsDstarTauNu_mu1_isTight.clear();
//  BsDstarTauNu_mu1_isPF.clear();
//  BsDstarTauNu_mu1_isGlobal.clear();
//  BsDstarTauNu_mu1_isTracker.clear();
//  BsDstarTauNu_mu1_isSoft.clear();
//  BsDstarTauNu_mu1_vx.clear();
//  BsDstarTauNu_mu1_vy.clear();
//  BsDstarTauNu_mu1_vz.clear();
//  BsDstarTauNu_mu1_iso.clear();
//  BsDstarTauNu_mu1_dbiso.clear();
//
//  BsDstarTauNu_tau_fullfit_pt.clear();
//  BsDstarTauNu_tau_fullfit_eta.clear();
//  BsDstarTauNu_tau_fullfit_phi.clear();
//  BsDstarTauNu_tau_fullfit_mass.clear();
//  BsDstarTauNu_tau_pt.clear();
//  BsDstarTauNu_tau_eta.clear();
//  BsDstarTauNu_tau_phi.clear();
//  BsDstarTauNu_tau_mass.clear();
//  BsDstarTauNu_tau_rhomass1.clear();
//  BsDstarTauNu_tau_rhomass2.clear();
//  BsDstarTauNu_tau_q.clear();
//  BsDstarTauNu_tau_vx.clear();
//  BsDstarTauNu_tau_vy.clear();
//  BsDstarTauNu_tau_vz.clear();
//
//
//  BsDstarTauNu_tau_max_dr_3prong.clear();
//  BsDstarTauNu_tau_lip.clear();
//  BsDstarTauNu_tau_lips.clear();
//  BsDstarTauNu_tau_pvip.clear();
//  BsDstarTauNu_tau_pvips.clear();
//  BsDstarTauNu_tau_fl3d.clear();
//  BsDstarTauNu_tau_fls3d.clear();
//  BsDstarTauNu_tau_alpha.clear();
//  BsDstarTauNu_tau_vprob.clear();
//  BsDstarTauNu_tau_isRight.clear();
//  BsDstarTauNu_tau_matched_ppdgId.clear();
//  BsDstarTauNu_tau_matched_gentaupt.clear();
//  BsDstarTauNu_tau_sumofdnn.clear();
//  BsDstarTauNu_tau_pfidx1.clear();
//  BsDstarTauNu_tau_pfidx2.clear();
//  BsDstarTauNu_tau_pfidx3.clear();
//
//  BsDstarTauNu_tau_pi1_pt.clear();
//  BsDstarTauNu_tau_pi1_eta.clear();
//  BsDstarTauNu_tau_pi1_phi.clear();
//  BsDstarTauNu_tau_pi1_mass.clear();
//  BsDstarTauNu_tau_pi2_pt.clear();
//  BsDstarTauNu_tau_pi2_eta.clear();
//  BsDstarTauNu_tau_pi2_phi.clear();
//  BsDstarTauNu_tau_pi2_mass.clear();
//  BsDstarTauNu_tau_pi3_pt.clear();
//  BsDstarTauNu_tau_pi3_eta.clear();
//  BsDstarTauNu_tau_pi3_phi.clear();
//  BsDstarTauNu_tau_pi3_mass.clear();
//
//  BsDstarTauNu_B_pt.clear();
//  BsDstarTauNu_B_eta.clear();
//  BsDstarTauNu_B_phi.clear();
//  BsDstarTauNu_B_mass.clear();
//  BsDstarTauNu_B_vprob.clear();
//  BsDstarTauNu_B_lip.clear();
//  BsDstarTauNu_B_lips.clear();
//  BsDstarTauNu_B_pvip.clear();
//  BsDstarTauNu_B_pvips.clear();
//  BsDstarTauNu_B_fl3d.clear();
//  BsDstarTauNu_B_fls3d.clear();
//  BsDstarTauNu_B_alpha.clear();
//  BsDstarTauNu_B_maxdoca.clear();
//  BsDstarTauNu_B_mindoca.clear();
//  BsDstarTauNu_B_vx.clear();
//  BsDstarTauNu_B_vy.clear();
//  BsDstarTauNu_B_vz.clear();
//  BsDstarTauNu_B_iso.clear();
//  BsDstarTauNu_B_iso_ntracks.clear();
//  BsDstarTauNu_B_iso_mindoca.clear();
//  BsDstarTauNu_B_mm2.clear();
//  BsDstarTauNu_B_q2.clear();
//  BsDstarTauNu_B_Es.clear();
//  BsDstarTauNu_B_ptback.clear();
//
//
//  BsDstarTauNu_PV_vx.clear();
//  BsDstarTauNu_PV_vy.clear();
//  BsDstarTauNu_PV_vz.clear();
//
//  BsDstarTauNu_bbPV_vx.clear();
//  BsDstarTauNu_bbPV_vy.clear();
//  BsDstarTauNu_bbPV_vz.clear();
//
//  BsDstarTauNu_bbPV_refit_vx.clear();
//  BsDstarTauNu_bbPV_refit_vy.clear();
//  BsDstarTauNu_bbPV_refit_vz.clear();
//
//  BsDstarTauNu_genPV_vx.clear();
//  BsDstarTauNu_genPV_vy.clear();
//  BsDstarTauNu_genPV_vz.clear();
//
//  BsDstarTauNu_Ds_pt.clear();
//  BsDstarTauNu_Ds_eta.clear();
//  BsDstarTauNu_Ds_phi.clear();
//  BsDstarTauNu_Ds_mass.clear();
//  BsDstarTauNu_Ds_vprob.clear();
//  BsDstarTauNu_Ds_lip.clear();
//  BsDstarTauNu_Ds_lips.clear();
//  BsDstarTauNu_Ds_pvip.clear();
//  BsDstarTauNu_Ds_pvips.clear();
//  BsDstarTauNu_Ds_fl3d.clear();
//  BsDstarTauNu_Ds_fls3d.clear();
//  BsDstarTauNu_Ds_alpha.clear();
//  BsDstarTauNu_Ds_vx.clear();
//  BsDstarTauNu_Ds_vy.clear();
//  BsDstarTauNu_Ds_vz.clear();
//  BsDstarTauNu_Ds_unfit_pt.clear();
//  BsDstarTauNu_Ds_unfit_mass.clear();
//  BsDstarTauNu_Ds_ptfrac.clear();
//
//  BsDstarTauNu_k_charge.clear();
//  BsDstarTauNu_pi_charge.clear();
//  BsDstarTauNu_spi_charge.clear();
//
//  BsDstarTauNu_D0_pt.clear();
//  BsDstarTauNu_D0_eta.clear();
//  BsDstarTauNu_D0_phi.clear();
//  BsDstarTauNu_D0_mass.clear();
//  BsDstarTauNu_D0_vprob.clear();
//  BsDstarTauNu_D0_lip.clear();
//  BsDstarTauNu_D0_lips.clear();
//  BsDstarTauNu_D0_pvip.clear();
//  BsDstarTauNu_D0_pvips.clear();
//  BsDstarTauNu_D0_fl3d.clear();
//  BsDstarTauNu_D0_fls3d.clear();
//  BsDstarTauNu_D0_alpha.clear();
//  BsDstarTauNu_D0_vx.clear();
//  BsDstarTauNu_D0_vy.clear();
//  BsDstarTauNu_D0_vz.clear();
//  BsDstarTauNu_D0_unfit_pt.clear();
//  BsDstarTauNu_D0_unfit_mass.clear();
//  BsDstarTauNu_D0_ptfrac.clear();
//
//  BsDstarTauNu_ngenmuons.clear();
//
//  BsDstarTauNu_isgen3.clear();
//  BsDstarTauNu_isgen3matched.clear();
//  BsDstarTauNu_nch.clear();
//  BsDstarTauNu_nch_after_dnn.clear();
//  BsDstarTauNu_nch_before_dnn.clear();
//  BsDstarTauNu_nch_qr.clear();
//  BsDstarTauNu_ngentau3.clear();
//  BsDstarTauNu_ngentau.clear();
//  BsDstarTauNu_gentaupt.clear();
//  BsDstarTauNu_gentaudm.clear();


  ////////////////////////


 
} 


void NtupleBranches::LabelHistograms( std::map< std::string, bool >& runFlags ){
  std::vector bins_string = {"Precut", "Trigger","muons", "J/#psi", "J/#psi fit","Tau  presence"};
  for(size_t i=0; i< bins_string.size(); i++){
    cutflow_perevt->GetXaxis()->SetBinLabel(i+1, bins_string[i]);
  }


  if ( runFlags["useHammer"] ){
    std::vector bins_hammer = {"den", "num", "num_a0_up", "num_a0_down", "num_a1_up", "num_a1_down","num_a2_up", "num_a2_down", "num_b0_up", "num_b0_down", "num_b1_up", "num_b1_down","num_b2_up", "num_b2_down", "num_c1_up", "num_c1_down", "num_c2_up", "num_c2_down", "num_d0_up", "num_d0_down", "num_d1_up", "num_d1_down","num_d2_up", "num_d2_down"};
    
    for(size_t i=0; i< bins_hammer.size(); i++){
      hammer_width->GetXaxis()->SetBinLabel(i+1, bins_hammer[i]);
    }
  }


  if ( runFlags["doGenHist"] ){
    std::vector bins_string = { "#mu", "#pi^{0}", "#pi^{#pm}","#rho^{0}","#rho^{#pm}","#eta","#eta^{`}","#omega","#phi","K^{0}","K^{#pm}","K^{*0}","K^{*#pm}","D^{#pm}","D^{0}","#eta_{c}","#eta_{b}","#Upsilon"};
    for(size_t i=0; i< bins_string.size(); i++){
      genParticle_Bdau_X_id->GetXaxis()->SetBinLabel(i+1, bins_string[i]);
    }       
  }
}
