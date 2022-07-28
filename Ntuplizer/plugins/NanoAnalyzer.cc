// **********************************************************************//
// Prototype implementation of NanoAnalyzer                              //
// for the creation of Run 1 nanoAOD-like ntuple from AOD or RECO input  //
// and corresp. Run 2 reference/validation ntuples from AOD or miniAOD   //
// **********************************************************************//

// ****************************************************
// for implementation history, see testnanoreadme.txt *
// ****************************************************
// Direct contributors:  A. Anuar, A. Bermudez, A. Geiser (coordinator), 
// N.Z. Jomhari, S. Wunsch, Q. Wang, H. Yang, Y. Yang, 2018-2021. 

// ******************************************
// automatically set appropriate CMSSW flag *
// ******************************************
// recognize automatically (no user action needed): 
// CMSSW is tied to particular compiler versions
// 42X taken from 4_2_8, 53X from 5_3_32, 7XX from 7_6_4
#define GCC_VERSION ( 10000 * __GNUC__ + 100 * __GNUC_MINOR__ + __GNUC_PATCHLEVEL__ )
#if GCC_VERSION < 40305
// Early Run 1 legacy CMSSW42X (e.g. 2010, 4_2_8)
#define CMSSW42X
#elif GCC_VERSION < 40703
// Main application is Run 1 legacy CMSSW53X (e.g. 2011/2, 5_3_32)
#define CMSSW53X
#elif GCC_VERSION > 40902
// instead in case CMSSW version is for Run 2, CMSSW7 or higher
// (e.g. 2015-18), for validation purposes
#define CMSSW7plus
#if GCC_VERSION < 50000
// for CMSSW 7_6_X  (GCC_VERSION might need to be changed/sharpened)
// (not clear whether this logic will still work for CMSSW 8 and 9)
#define CMSSW7XX
#endif
#endif
// for ultra-legacy AOD (10_6_4 or higher)
#if GCC_VERSION > 70400
// this one comes on top of the previous
#define CMSSW106plus
#endif
#if GCC_VERSION > 80300
// for Run 3 MC studies
// this one comes on top of the previous
// GCC_VERSION preliminary, might need to be changed/sharpened
#define CMSSW11plus
#endif
#if GCC_VERSION > 90299
// for 2021 pilot data 
// this one comes on top of the previous
#define CMSSW12plus
#endif

// ********************************************************************
// the following are flags which can be (de)activated manually
// defaults are set such that for normal use no manual action is needed
// ********************************************************************

// activate this to achieve maximum compatibility with official Run 2 nanoAOD
// (default is off -> better performance, e.g. use of beam spot constraint)
// (mainly for Run 2 validation - not recommended for Run 1)
// only relevant for Run 2 AOD, turn on for better consistency with miniAOD 
//                        (slightly worse performance)
// nanoext flag can also be steered/overruled via configuration file
//#define Compatibility

// activation flag to check and store compatibility with Golden JSON 
// default is on for 2010 data (not MC, automatic), off otherwise
#ifdef CMSSW42X
#define JSONcheck
#endif

#ifdef CMSSW42X
// activate this to check 2010 Golden JSON setting, i.e. abort upon nonJSON event
// activate this only when Golden JSON is activated in configuration
// (protection against accidental use of wrong or no JSON in configuration,
//  for check and validation, not strictly needed)
#define JSONcheckabort
#endif

// activate this only for data sets for which "plus" part of trigger treatment 
// is already implemented; checks and aborts in case of inconsistency;
// should be activated by default if trigger is implemented for dataset
// (protection against inconsistencies in NanoTrigger implementation)
//#define trigcheckabort

#ifdef CMSSW7plus
// activate this when you read from miniAOD for validation (Run 2/3 only!)
//#define miniAOD
#endif


// turn this on to deactivate code related to jet corrections
//  (e.g. in case of problems with the configuration, and only if you do 
//   *not* use jets in your analysis)
//  default should be flag off
#define noJetCor 

// turn this on to activate code related to charm final states
#define charm

#ifdef charm
// turn this on to activate D meson cuts optimized for large rapidities
// default (on) is for new looser cuts
#define beauty
#endif

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>

#include "../interface/helperfunc.h"
#include "../interface/MyStruct.h"

#ifdef CMSSW42X
#include <boost/unordered_map.hpp>
using boost::unordered_map;
#else
#include <unordered_map>
using std::unordered_map;
#endif

//*****************************
// general user include files *
//*****************************
#include "FWCore/Framework/interface/Frameworkfwd.h"
#ifndef CMSSW12plus
// Run 1 and 2
#include "FWCore/Framework/interface/EDAnalyzer.h"
#else
// Run 3
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#endif
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for fillDescriptions
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ------ EXTRA HEADER FILES--------------------//
#include "FWCore/Framework/interface/EventSetup.h"
#ifndef CMSSW12plus
#include "FWCore/Framework/interface/ESHandle.h"
#endif
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"
// the following does not seem to be needed (included through dEdx), but also not to hurt
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
// to get release version
#include "FWCore/Version/interface/GetReleaseVersion.h"
// header for conversion tools
#ifndef CMSSW11plus
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#else
#include "CommonTools/Egamma/interface/ConversionTools.h"
#endif
// effective area for rho
//https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/EffectiveAreas.h
//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

// for math and file handling
#include "TMath.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

// for ntuple output
#include "TLorentzVector.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#ifdef charm
// for control histograms
#include "TH1.h"
#include "TH2.h"
#endif

//**************************
// for trigger information *
//**************************
#include "FWCore/Common/interface/TriggerNames.h"
// not needed?
//#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/TriggerResults.h"
// for automatic trigger recognition
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
// for trigger objects
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <cassert> 
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#ifndef CMSSW11plus
// no longer exists in CMSSW_11_2_X, not needed??
#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"
#endif

//***************************
// for tracking information *
//***************************
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
// for dEdx
#include "DataFormats/TrackReco/interface/DeDxData.h" 	

//#ifdef miniAOD
// tracks from PATCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// there are three track collections in miniAOD (in addition to the muon 
// and electron collections), embedded into particleflow objects:
// packedPFCandidates allows to rebuild tracks (for pt>0.5 GeV) using
//                    pseudotrack() (object) or besttrack() (pointer)
// PackedPFCandidatesDiscarded presumably contains only discarded duplicate 
//                    muon candidates -> do not use
// lostTracks (high purity only) contains some of the non-vertex tracks 
//#endif

/// for track parametrization 
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"

//*************************
// for vertex information *
//*************************
// reconstructed primary is typically within 0.02 cm of true primary in z
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>

// for vertices refit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// for beamspot information (beam pipe radius: 5.8 -> 4.3 cm, 
//                           beam spot at x~0.2, y~0.4, z~0.3 cm,
//                           width x~0.002?, y~0.002?, z~5.5 cm,
// beam-spot constrained vertices: x~0.001 , y~0.001 , z~0.004 cm) 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//****************************
// for error matrix handling *
//****************************
#include "TMatrixT.h"
#include "Math/SMatrix.h"
#include "Math/StaticCheck.h"
//#include <MatrixRepresentationsStatic.h>
#include "Math/Expression.h"

// set namespaces
   using namespace edm;
   using namespace reco;
   using namespace std;

//***********************
// for muon information *
//***********************


#include "DataFormats/PatCandidates/interface/Muon.h"

//***************************
// for electron information *
//***************************


#include "DataFormats/PatCandidates/interface/Electron.h"



// for photon information
#include "DataFormats/PatCandidates/interface/Photon.h"

//**********************
// for MET information *
#include "DataFormats/PatCandidates/interface/MET.h"

//**********************
// for jet information *
//**********************
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


//*******************************
// for gen particle information *
//*******************************
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//********************************
// for particle flow information *
//********************************
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//***************************
// class member declaration *
//***************************

//***********************************
// main analyzer class (EDAnalyzer) *
//***********************************

#ifndef CMSSW12plus
// Run 1 and 2
class NanoAnalyzer : public edm::EDAnalyzer
#else
// Run 3
class NanoAnalyzer : public edm::one::EDAnalyzer<>
#endif
{
public:
  explicit NanoAnalyzer(const edm::ParameterSet&);
  ~NanoAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

  // this is the place to define global variables and parameters 

        // declare global trigger variables
  //#include "NanoTrigger.h"

private:

  virtual void beginJob(const edm::ParameterSet& iConfig);
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  void createBranch();

  void reset();
  
  // as it says on the tin
  // branch title creation 
  //  void eventDoc();
  //  void createTitle(const std::string &name, const std::string &title); 

  // HLT config for reading the table and its associated process name
  //HLTConfigProvider hlt_cfg;   // superseded above


  EDGetTokenT<pat::ElectronCollection>    electronToken_;
  EDGetTokenT<reco::VertexCollection> verticeToken_;
  EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken_;
  EDGetTokenT<edm::TriggerResults> triggerToken_;  
  EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobjectToken_;

  ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken;
  ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken;


  Handle<pat::ElectronCollection>      		       electrons_;
  Handle< reco::VertexCollection >                        vertices_;
  Handle< std::vector<pat::PackedCandidate> >             packedpfcandidates_   ;
  Handle< edm::TriggerResults> 			       HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;

  //  const MagneticField                 *fMagneticField;

/////////////////////////////////////////////////////////////////////////
////////////////////////// declare tree, file, //////////////////////////
/////////////////////////////////////////////////////////////////////////
  
  edm::Service<TFileService> fs;
  TTree* tree_;

  //  TFile *file;
  //  TTree *tree_;

  /// original nanoAOD ///
  UInt_t run;
  ULong64_t event;
  UInt_t luminosityBlock;


  float                JpsiKE_e1_pt      ;
  float                JpsiKE_e1_eta     ;
  float                JpsiKE_e1_phi     ;
  float                JpsiKE_e1_mass     ;
  int                  JpsiKE_e1_q      ;   
  float                JpsiKE_e1_vx       ;
  float                JpsiKE_e1_vy       ;
  float                JpsiKE_e1_vz       ;

  float                JpsiKE_e2_pt      ;
  float                JpsiKE_e2_eta     ;
  float                JpsiKE_e2_phi     ;
  float                JpsiKE_e2_mass     ;
  int                  JpsiKE_e2_q   ;   
  float                JpsiKE_e2_vx       ;
  float                JpsiKE_e2_vy       ;
  float                JpsiKE_e2_vz       ;

  float                JpsiKE_PV_vx       ;
  float                JpsiKE_PV_vy       ;
  float                JpsiKE_PV_vz       ;
  
  float                JpsiKE_bbPV_vx       ;
  float                JpsiKE_bbPV_vy       ;
  float                JpsiKE_bbPV_vz       ;
  float                JpsiKE_bbPV_chi2       ;
  float                JpsiKE_bbPV_ndof       ;
  float                JpsiKE_bbPV_rho       ;

  float                JpsiKE_Jpsi_pt      ;
  float                JpsiKE_Jpsi_eta     ;
  float                JpsiKE_Jpsi_phi     ;
  float                JpsiKE_Jpsi_mass       ;
  float                JpsiKE_Jpsi_vprob    ;
  float                JpsiKE_Jpsi_lip;
  float                JpsiKE_Jpsi_lips;
  float                JpsiKE_Jpsi_pvip;
  float                JpsiKE_Jpsi_pvips;
  float                JpsiKE_Jpsi_fl3d;
  float                JpsiKE_Jpsi_fls3d;
  float                JpsiKE_Jpsi_alpha;
  float                JpsiKE_Jpsi_maxdoca;
  float                JpsiKE_Jpsi_mindoca;
  float                JpsiKE_Jpsi_vx      ;
  float                JpsiKE_Jpsi_vy      ;
  float                JpsiKE_Jpsi_vz      ;

  std::vector<float>                JpsiKE_B_pt      ;
  std::vector<float>                JpsiKE_B_eta     ;
  std::vector<float>                JpsiKE_B_phi     ;
  std::vector<float>                JpsiKE_B_mass    ;
  std::vector<float>                JpsiKE_B_mass_nofit    ;
  std::vector<float>                JpsiKE_B_mcorr    ;
  std::vector<float>                JpsiKE_B_vprob ;
  std::vector<float>                JpsiKE_B_lip;
  std::vector<float>                JpsiKE_B_lips;
  std::vector<float>                JpsiKE_B_pvip;
  std::vector<float>                JpsiKE_B_pvips;
  std::vector<float>                JpsiKE_B_fl3d;
  std::vector<float>                JpsiKE_B_fls3d;
  std::vector<float>                JpsiKE_B_alpha;
  std::vector<float>                JpsiKE_B_maxdoca;
  std::vector<float>                JpsiKE_B_mindoca;
  std::vector<float>                JpsiKE_B_vx      ;
  std::vector<float>                JpsiKE_B_vy      ;
  std::vector<float>                JpsiKE_B_vz      ;


  std::vector<float>       JpsiKE_pi_pt;
  std::vector<float>       JpsiKE_pi_eta;
  std::vector<float>       JpsiKE_pi_phi;
  std::vector<float>       JpsiKE_pi_mass;
  std::vector<int>       JpsiKE_pi_q;
  std::vector<int>       JpsiKE_pi_pdg;

  std::vector<float> JpsiKE_pi_doca3d;
  std::vector<float> JpsiKE_pi_doca3de;
  std::vector<float> JpsiKE_pi_doca2d;
  std::vector<float> JpsiKE_pi_doca2de;
  std::vector<float> JpsiKE_pi_dz;
  std::vector<float> JpsiKE_pi_near_dz;
  std::vector<bool> JpsiKE_pi_isAssociate;
  std::vector<int> JpsiKE_pi_pvAssociationQuality;

  std::vector<float>                JpsiKE_B_q2;
  std::vector<float>                JpsiKE_B_mm2;
  std::vector<float>                JpsiKE_B_ptmiss;
  std::vector<float>                JpsiKE_B_Es;
  std::vector<float>                JpsiKE_B_ptback;

  
  int JpsiKE_nch;
  int JpsiKE_nch_filter;


  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  TH1F * hist; 

}; // end of class member



//////////////////////////////////////////////////////////////////////////////
//                        set analysis loop parameters                      //
//////////////////////////////////////////////////////////////////////////////

NanoAnalyzer::NanoAnalyzer(const edm::ParameterSet& iConfig)
{

  std::cout << "DC: NanoAnalyzer" << std::endl;

  electronToken_           = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));

  packedpfcandidatesToken_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates")); 


  verticeToken_            = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  triggerToken_	      	   = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_	   = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));


  ttkToken = esConsumes(edm::ESInputTag{"","TransientTrackBuilder"});
  fieldToken = esConsumes(edm::ESInputTag{"",""});

  //  aux = new helperfunc();

  //  file = new TFile("test.root", "recreate"); // check
  //  tree_ = new TTree("Events", "Events");
  //  tree_->SetImplicitMT(false);
  //tree_->Branch("run", &run, "run/i");

  std::cout << "test" << std::endl;
  //hist = new TH1F("cutflow", "cutflow", 10,0,10);
  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);

  std::cout << "test2" << std::endl;

  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();

  std::cout << "test3" << std::endl;

} // end of constructor

NanoAnalyzer::~NanoAnalyzer()
{
  std::cout << "Finish NanoAnalyzer" << std::endl;

  //  file->cd();
  //  tree_->Write();
  //  hist->Write();

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

///////////////////////////////////////////////////////////////////////////////
////////////// main analysis loop: method called for each event ///////////////
///////////////////////////////////////////////////////////////////////////////

void
NanoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  reset();

  iEvent.getByToken(triggerToken_, HLTtriggers_);

  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();

  
  hist->Fill(1);

  bool isTriggered = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string finalTriggerName="";
  //  std::string finalTriggerFilterObjName="";

  //  std::cout << "------------------------------" << std::endl;

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
    
    //    std::cout << "trig_name:" <<trigNames.triggerName(i) << std::endl;

    if(trigNames.triggerName(i).find("eta1p22_mMax6")!= std::string::npos){
      //      nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	finalTriggerName=trigNames.triggerName(i);  
	
	std::cout << "!!!!!!! This is fired: " << finalTriggerName << std::endl;
	//	finalTriggerFilterObjName="hltJpsiTkVertexFilter";
      }
    }
  }


  if(!isTriggered) return; 
  //  nBranches_->cutflow->Fill(1);

  hist->Fill(2);


  iEvent.getByToken(electronToken_	, electrons_    );
  iEvent.getByToken(triggerobjectToken_ , triggerObjects);
  
  std::vector<pat::Electron> electroncollection;

  electroncollection.clear();


  for(size_t ielectron = 0; ielectron < electrons_->size(); ++ ielectron){

    const pat::Electron & electron = (*electrons_)[ielectron];

    if(electron.pt() < 2.5) continue;
    if(TMath::Abs(electron.eta()) > 2.4) continue;
    if (!electron.passConversionVeto()) continue;
    const reco::GsfTrackRef gsfTrk = electron.gsfTrack();

    if(!gsfTrk.isNonnull()) continue;


    /// Trigger matching 


////    bool trigMatch = false;
////
////    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
////    
////      obj.unpackPathNames(trigNames);
////      obj.unpackFilterLabels(iEvent, *HLTtriggers_);
////
////      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
////
////      bool isPathExist = false;
////
////      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
////	if(pathNamesAll[h]==finalTriggerName) isPathExist = true;
////      }
////
////      //      std::cout << "isPathExist!" << std::endl;
////      if(!isPathExist) continue;
////      //      std::cout << "isPathExist passed!" << std::endl;
////
/////////      bool isFilterExist = false;
/////////    
/////////      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
/////////	
/////////	//	if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
/////////	  //	if(obj.filterLabels()[hh].find("hltDisplacedmumuFilterDoubleMu4Jpsi") != std::string::npos){
/////////	  isFilterExist = true;
/////////	}
/////////      }
/////////      
/////////      if(!isFilterExist) continue;
/////////      std::cout << "isFireExist passed!" << std::endl;
////      
////      //      if(TMath::Abs(obj.pdgId()) != 13) continue;
////      
////      Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
////					electron.eta(), electron.phi());
////      
////      if(trigger_dR < 0.015 &&
////	 obj.pt()/electron.pt() > 0.85 &&
////	 obj.pt()/electron.pt() < 1.15
////	 ){
////	trigMatch = true;
////	//	std::cout << "Muon" << imuon << " matches to the trigger object dR = " << trigger_dR << ", obj pT = " << obj.pt() << ", muon pT = " << muon.pt() << " " << obj.pdgId() << std::endl;
////      }
////    }
////
////    if(!trigMatch) continue;





    electroncollection.push_back(electron);
  }



  if(electroncollection.size() < 2) return;
  hist->Fill(3);


  const TransientTrackBuilder* builder = &iSetup.getData(ttkToken);

//  edm::ESHandle<MagneticField> fieldHandle;
//  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
//  fMagneticField = fieldHandle.product();


  const MagneticField *fMagneticField = &iSetup.getData(fieldToken);
  //  fMagneticField = fieldHandle.product();



  float jpsi_max_pt = -1;
  unsigned int mcidx_e1 = -1;
  unsigned int mcidx_e2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int ie = 0; ie < (int)electroncollection.size(); ie++){
    for(int je = ie+1; je < (int)electroncollection.size(); je++){

      const pat::Electron e1 = electroncollection[ie];
      const pat::Electron e2 = electroncollection[je];

      TLorentzVector tlv_e1;
      TLorentzVector tlv_e2;

      //      std::cout << e1.mass() << " " << e2.mass() << std::endl;
      tlv_e1.SetPtEtaPhiM(e1.pt(), e1.eta(), e1.phi(), aux.mass_electron);
      tlv_e2.SetPtEtaPhiM(e2.pt(), e2.eta(), e2.phi(), aux.mass_electron);

      //      std::cout << "e1 mass = " << e1.mass() << std::endl;

      TLorentzVector tlv_jpsi = (tlv_e1 + tlv_e2);

      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      if(e1.charge() + e2.charge() !=0) continue;
      if(jpsi_mass < 2.95) continue; // a little bit broad winder to take into account FSR ...
      if(jpsi_mass > 3.25) continue;

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_e1 = ie;
	mcidx_e2 = je;
	jpsi_tlv_highest = tlv_jpsi;
      }
    }
  }

  if(jpsi_max_pt == -1) return;
  std::cout << "J/psi found!" << std::endl;
  hist->Fill(4);

  const reco::GsfTrackRef track1_electron = electroncollection[mcidx_e1].gsfTrack();
  const reco::GsfTrackRef track2_electron = electroncollection[mcidx_e2].gsfTrack();
  reco::TransientTrack tt1_electron = (*builder).build(track1_electron);
  reco::TransientTrack tt2_electron = (*builder).build(track2_electron);


  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> electronParticles;

  electronParticles.push_back(pFactory.particle(tt1_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));
  electronParticles.push_back(pFactory.particle(tt2_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));

  
  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;

  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(electronParticles, -1, -1);


  if(!jpsifit_flag) return;


  hist->Fill(5);

  std::vector< RefCountedKinematicParticle > jpsi_children = jpTree->finalStateParticles();

  //  math::PtEtaPhiMLorentzVector e1_fit = aux.daughter_p4(jpsi_children, 0);
  //  math::PtEtaPhiMLorentzVector e2_fit = aux.daughter_p4(jpsi_children, 1);

  std::vector<pat::Electron> electroncollection_selected;
  electroncollection_selected.push_back(electroncollection[mcidx_e1]);
  electroncollection_selected.push_back(electroncollection[mcidx_e2]);

  iEvent.getByToken(verticeToken_   , vertices_     );

  // edm::ESHandle<MagneticField> fieldHandle;
  //  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  //  fMagneticField = fieldHandle.product();


  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TransverseImpactPointExtrapolator extrapolatort(fMagneticField);


  float max_criteria = 999;
  reco::Vertex closestVertex; 
  int counter = 0;

  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){

    particle_cand cand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, *vtx);
 
    if(TMath::Abs(cand.lip) < max_criteria){
      max_criteria = TMath::Abs(cand.lip);
      closestVertex = *vtx;
    }

    counter += 1;
  }


  iEvent.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  Int_t nfilter = 0; 


  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];

    
    if(pf.pt() < 0.5) continue;
    if(!pf.hasTrackDetails()) continue;

  
    // use the PF candidates that come from closestVertex    
    //  float precut_dz = pf.vz() - closestVertex.position().z();
    //  if(TMath::Abs(precut_dz) > c_dz) continue;
    
    Bool_t hpflag = pf.trackHighPurity();
    if(!hpflag) continue;

    if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
    if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
    if(pf.pseudoTrack().normalizedChi2() > 100) continue;

    //    if(TMath::Abs(pf.pdgId())!=211) continue;
    if(TMath::Abs(pf.eta()) > 2.5) continue; 
    


    float precut_dz = pf.vz() - closestVertex.position().z();
    if(TMath::Abs(precut_dz) > 0.12) continue;



    reco::TransientTrack  track = (*builder).build(pf.pseudoTrack());
  

    TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(track.impactPointState(), jpsi_vertex->position());


    TrajectoryStateOnSurface _tsost_pf = extrapolatort.extrapolate(track.impactPointState(), jpsi_vertex->position());
    

    std::pair<bool,Measurement1D> _cur3DIP_pf = aux.signedImpactParameter3D(_tsos_pf, jpsi_vertex, closestVertex);

    std::pair<bool,Measurement1D> _cur3DIPt_pf = aux.signedTransverseImpactParameter(_tsost_pf, jpsi_vertex, closestVertex);
    
    float doca3d = _cur3DIP_pf.second.value();
    float doca3de = _cur3DIP_pf.second.error();
    //    float doca3ds = _cur3DIP_pf.second.significance();

    float doca2d = _cur3DIPt_pf.second.value();
    float doca2de = _cur3DIPt_pf.second.error();
    //    float doca2ds = _cur3DIPt_pf.second.significance();


    float near_dz = -999;
    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
      
      if(pf.vertexRef()->z()==vtx->position().z()){
	near_dz = closestVertex.position().z() - vtx->position().z();
      }
    }


    Bool_t isAssociate = (bool)(pf.vertexRef()->z()==closestVertex.position().z());


    float dr_jpsi = reco::deltaR(pf.eta(),
				 pf.phi(),
				 jpsi_part->currentState().globalMomentum().eta(),
				 jpsi_part->currentState().globalMomentum().phi());


    //    if(doca3d < -0.04 || doca3d > 0.06) continue;
    if(dr_jpsi > 1.) continue;



    //    pat::PackedCandidate pf = pfcands[iii].pfcand;
    //    reco::TransientTrack track = pfcands[iii].track;
    //    attribute attr = pfcands[iii].pfaux;

    std::vector<RefCountedKinematicParticle> allParticles;
      
    allParticles.push_back(pFactory.particle(track, aux.kaon_mass, chi, ndf, aux.kaon_sigma));
    allParticles.push_back(pFactory.particle(tt1_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));
    allParticles.push_back(pFactory.particle(tt2_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));


    RefCountedKinematicParticle b_part;
    RefCountedKinematicVertex b_vertex;
    RefCountedKinematicTree bTree;
    Bool_t bfit_flag;
    std::tie(bfit_flag, b_part, b_vertex, bTree) = aux.KinematicFit(allParticles, -1, -1);
    if(!bfit_flag){
      std::cout <<"The Fit fails ..." << std::endl;
      continue;
    }


    //    if(TMath::Prob(b_vertex->chiSquared(), b_vertex->degreesOfFreedom()) <= 0.05) continue;
        


    particle_cand Bcand; 
    Bcand = aux.calculateIPvariables(extrapolator, b_part, b_vertex, closestVertex);


    TLorentzVector tlv_k;
    tlv_k.SetPtEtaPhiM(pf.pt(), pf.eta(), pf.phi(), aux.mass_kaon);

    


    nfilter+=1;

    JpsiKE_B_pt.push_back(b_part->currentState().globalMomentum().perp());
    JpsiKE_B_eta.push_back(b_part->currentState().globalMomentum().eta());
    JpsiKE_B_phi.push_back(b_part->currentState().globalMomentum().phi());
    JpsiKE_B_mass.push_back(b_part->currentState().mass());
    JpsiKE_B_mass_nofit.push_back((tlv_k + jpsi_tlv_highest).M());
    JpsiKE_B_vprob.push_back( TMath::Prob(b_vertex->chiSquared(), b_vertex->degreesOfFreedom()) ); //b_part->currentState().mass());
    JpsiKE_B_lip.push_back(Bcand.lip);
    JpsiKE_B_lips.push_back(Bcand.lips);
    JpsiKE_B_pvip.push_back(Bcand.pvip);
    JpsiKE_B_pvips.push_back(Bcand.pvips);
    JpsiKE_B_fls3d.push_back(Bcand.fls3d);
    JpsiKE_B_fl3d.push_back(Bcand.fl3d);
    JpsiKE_B_alpha.push_back(Bcand.alpha);    
    JpsiKE_B_vx.push_back(b_vertex->vertexState().position().x());
    JpsiKE_B_vy.push_back(b_vertex->vertexState().position().y());
    JpsiKE_B_vz.push_back(b_vertex->vertexState().position().z());

    JpsiKE_pi_pt.push_back(pf.pt());
    JpsiKE_pi_eta.push_back(pf.eta());
    JpsiKE_pi_phi.push_back(pf.phi());
    JpsiKE_pi_mass.push_back(pf.mass());
    JpsiKE_pi_q.push_back(pf.charge());
    JpsiKE_pi_doca3d.push_back(doca3d);
    JpsiKE_pi_doca3de.push_back(doca3de);
    JpsiKE_pi_doca2d.push_back(doca2d);
    JpsiKE_pi_doca2de.push_back(doca2de);
    JpsiKE_pi_dz.push_back(precut_dz);
    JpsiKE_pi_near_dz.push_back(near_dz);
    JpsiKE_pi_isAssociate.push_back(isAssociate);
    JpsiKE_pi_pvAssociationQuality.push_back(pf.pvAssociationQuality());
    JpsiKE_pi_pdg.push_back(pf.pdgId());



  }

  
  JpsiKE_e1_pt = electroncollection[mcidx_e1].pt();
  JpsiKE_e1_eta = electroncollection[mcidx_e1].eta();
  JpsiKE_e1_phi = electroncollection[mcidx_e1].phi();
  JpsiKE_e1_mass = electroncollection[mcidx_e1].mass();
  JpsiKE_e1_q = electroncollection[mcidx_e1].charge();
  JpsiKE_e1_vx = electroncollection[mcidx_e1].vx();
  JpsiKE_e1_vy = electroncollection[mcidx_e1].vy();
  JpsiKE_e1_vz = electroncollection[mcidx_e1].vz();
  
  JpsiKE_e2_pt = electroncollection[mcidx_e2].pt();
  JpsiKE_e2_eta = electroncollection[mcidx_e2].eta();
  JpsiKE_e2_phi = electroncollection[mcidx_e2].phi();
  JpsiKE_e2_mass = electroncollection[mcidx_e2].mass();
  JpsiKE_e2_q = electroncollection[mcidx_e2].charge();
  JpsiKE_e2_vx = electroncollection[mcidx_e2].vx();
  JpsiKE_e2_vy = electroncollection[mcidx_e2].vy();
  JpsiKE_e2_vz = electroncollection[mcidx_e2].vz();


  JpsiKE_PV_vx = vertices_->begin()->position().x();
  JpsiKE_PV_vy = vertices_->begin()->position().y();
  JpsiKE_PV_vz = vertices_->begin()->position().z();


  JpsiKE_bbPV_vx = closestVertex.position().x();
  JpsiKE_bbPV_vy = closestVertex.position().y();
  JpsiKE_bbPV_vz = closestVertex.position().z();
  JpsiKE_bbPV_chi2 = closestVertex.chi2();
  JpsiKE_bbPV_ndof = closestVertex.ndof();
  JpsiKE_bbPV_rho = closestVertex.position().Rho();


  JpsiKE_Jpsi_pt = jpsi_part->currentState().globalMomentum().perp();
  JpsiKE_Jpsi_eta = jpsi_part->currentState().globalMomentum().eta();
  JpsiKE_Jpsi_phi = jpsi_part->currentState().globalMomentum().phi();
  JpsiKE_Jpsi_mass = jpsi_part->currentState().mass();
  JpsiKE_Jpsi_vprob = TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom());


  particle_cand JPcand;
  JPcand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);


  JpsiKE_Jpsi_lip = JPcand.lip;
  JpsiKE_Jpsi_lips = JPcand.lips;
  JpsiKE_Jpsi_pvip = JPcand.pvip;
  JpsiKE_Jpsi_pvips = JPcand.pvips;
  JpsiKE_Jpsi_fl3d = JPcand.fl3d;
  JpsiKE_Jpsi_fls3d = JPcand.fls3d;
  JpsiKE_Jpsi_alpha = JPcand.alpha;
  JpsiKE_Jpsi_maxdoca = aux.getMaxDoca(electronParticles);
  JpsiKE_Jpsi_mindoca = aux.getMinDoca(electronParticles);
  JpsiKE_Jpsi_vx = jpsi_vertex->vertexState().position().x();
  JpsiKE_Jpsi_vy = jpsi_vertex->vertexState().position().y();
  JpsiKE_Jpsi_vz = jpsi_vertex->vertexState().position().z();  


  JpsiKE_nch = (int)packedpfcandidates_->size();
  JpsiKE_nch_filter = nfilter;
  tree_->Fill();

  return;

}//NanoAnalyzer::analyze ends



//**************************************************
//************* additional methods *****************
//**************************************************

// ------------ method called once each job just before starting event loop  ------------
void NanoAnalyzer::beginJob(const edm::ParameterSet& iConfig)
{


}

// ------------ method called once each job just before starting run  ------------
void NanoAnalyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp)
{
  std::cout << "beginRun" << std::endl;

  //If the hltConfig can be initialized, then the below is an example of how to extract the config information for the trigger from the so-called provenance.
  //The trigger configuration can change from run to run (during the run is the same), so it needs to be called here.
  //"init" return value indicates whether intitialisation has succeeded
  //"changed" parameter indicates whether the config has actually changed
  //cout << "hello beginRun end" << endl; 
}

void NanoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) 
{
  std::cout << "fillDescriptions!" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------  Qun below
void NanoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{

  std::cout << "endRun" << std::endl;

}

void 
NanoAnalyzer::endJob()
{
}

//include methods for special trigger variables
//#include "NanoTrigger.cc.forinclude"

//define this as a plug-in

// branch title creation
void NanoAnalyzer::createBranch() { 

  std::cout << "createBranch" << std::endl;

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");

  tree_->Branch("JpsiKE_e1_pt", &JpsiKE_e1_pt );
  tree_->Branch("JpsiKE_e1_eta", &JpsiKE_e1_eta );
  tree_->Branch("JpsiKE_e1_phi", &JpsiKE_e1_phi );
  tree_->Branch("JpsiKE_e1_mass", &JpsiKE_e1_mass );
  tree_->Branch("JpsiKE_e1_q", &JpsiKE_e1_q );
  tree_->Branch("JpsiKE_e1_vx"   , &JpsiKE_e1_vx    );
  tree_->Branch("JpsiKE_e1_vy"   , &JpsiKE_e1_vy    );
  tree_->Branch("JpsiKE_e1_vz"   , &JpsiKE_e1_vz    );

  tree_->Branch("JpsiKE_e2_pt", &JpsiKE_e2_pt );
  tree_->Branch("JpsiKE_e2_eta", &JpsiKE_e2_eta );
  tree_->Branch("JpsiKE_e2_phi", &JpsiKE_e2_phi );
  tree_->Branch("JpsiKE_e2_mass", &JpsiKE_e2_mass );
  tree_->Branch("JpsiKE_e2_q", &JpsiKE_e2_q );
  tree_->Branch("JpsiKE_e2_vx"   , &JpsiKE_e2_vx    );
  tree_->Branch("JpsiKE_e2_vy"   , &JpsiKE_e2_vy    );
  tree_->Branch("JpsiKE_e2_vz"   , &JpsiKE_e2_vz    );

  tree_->Branch("JpsiKE_PV_vx", &JpsiKE_PV_vx );
  tree_->Branch("JpsiKE_PV_vy", &JpsiKE_PV_vy );
  tree_->Branch("JpsiKE_PV_vz", &JpsiKE_PV_vz );

  tree_->Branch("JpsiKE_bbPV_vx", &JpsiKE_bbPV_vx );
  tree_->Branch("JpsiKE_bbPV_vy", &JpsiKE_bbPV_vy );
  tree_->Branch("JpsiKE_bbPV_vz", &JpsiKE_bbPV_vz );
  tree_->Branch("JpsiKE_bbPV_chi2", &JpsiKE_bbPV_chi2 );
  tree_->Branch("JpsiKE_bbPV_ndof", &JpsiKE_bbPV_ndof );
  tree_->Branch("JpsiKE_bbPV_rho", &JpsiKE_bbPV_rho );


  tree_->Branch("JpsiKE_Jpsi_pt", &JpsiKE_Jpsi_pt );
  tree_->Branch("JpsiKE_Jpsi_eta", &JpsiKE_Jpsi_eta );
  tree_->Branch("JpsiKE_Jpsi_phi", &JpsiKE_Jpsi_phi );
  tree_->Branch("JpsiKE_Jpsi_mass", &JpsiKE_Jpsi_mass );
  tree_->Branch("JpsiKE_Jpsi_vprob", &JpsiKE_Jpsi_vprob );
  tree_->Branch("JpsiKE_Jpsi_lip", &JpsiKE_Jpsi_lip);
  tree_->Branch("JpsiKE_Jpsi_lips", &JpsiKE_Jpsi_lips);
  tree_->Branch("JpsiKE_Jpsi_pvip", &JpsiKE_Jpsi_pvip);
  tree_->Branch("JpsiKE_Jpsi_pvips", &JpsiKE_Jpsi_pvips);
  tree_->Branch("JpsiKE_Jpsi_fl3d", &JpsiKE_Jpsi_fl3d);
  tree_->Branch("JpsiKE_Jpsi_fls3d", &JpsiKE_Jpsi_fls3d);
  tree_->Branch("JpsiKE_Jpsi_alpha", &JpsiKE_Jpsi_alpha);
  tree_->Branch("JpsiKE_Jpsi_maxdoca", &JpsiKE_Jpsi_maxdoca);
  tree_->Branch("JpsiKE_Jpsi_mindoca", &JpsiKE_Jpsi_mindoca);
  tree_->Branch("JpsiKE_Jpsi_vx", &JpsiKE_Jpsi_vx );
  tree_->Branch("JpsiKE_Jpsi_vy", &JpsiKE_Jpsi_vy );
  tree_->Branch("JpsiKE_Jpsi_vz", &JpsiKE_Jpsi_vz );


  tree_->Branch("JpsiKE_B_pt", &JpsiKE_B_pt );
  tree_->Branch("JpsiKE_B_eta", &JpsiKE_B_eta );
  tree_->Branch("JpsiKE_B_phi", &JpsiKE_B_phi );
  tree_->Branch("JpsiKE_B_mass", &JpsiKE_B_mass );
  tree_->Branch("JpsiKE_B_mass_nofit", &JpsiKE_B_mass_nofit );
  tree_->Branch("JpsiKE_B_mcorr", &JpsiKE_B_mcorr );
  tree_->Branch("JpsiKE_B_vprob", &JpsiKE_B_vprob );
  tree_->Branch("JpsiKE_B_lip", &JpsiKE_B_lip);
  tree_->Branch("JpsiKE_B_lips", &JpsiKE_B_lips);
  tree_->Branch("JpsiKE_B_pvip", &JpsiKE_B_pvip);
  tree_->Branch("JpsiKE_B_pvips", &JpsiKE_B_pvips);
  tree_->Branch("JpsiKE_B_fl3d", &JpsiKE_B_fl3d);
  tree_->Branch("JpsiKE_B_fls3d", &JpsiKE_B_fls3d);
  tree_->Branch("JpsiKE_B_alpha", &JpsiKE_B_alpha);
  tree_->Branch("JpsiKE_B_maxdoca", &JpsiKE_B_maxdoca);
  tree_->Branch("JpsiKE_B_mindoca", &JpsiKE_B_mindoca);
  tree_->Branch("JpsiKE_B_vx", &JpsiKE_B_vx );
  tree_->Branch("JpsiKE_B_vy", &JpsiKE_B_vy );
  tree_->Branch("JpsiKE_B_vz", &JpsiKE_B_vz );



  tree_->Branch("JpsiKE_pi_pt", &JpsiKE_pi_pt );
  tree_->Branch("JpsiKE_pi_eta", &JpsiKE_pi_eta );
  tree_->Branch("JpsiKE_pi_phi", &JpsiKE_pi_phi );
  tree_->Branch("JpsiKE_pi_mass", &JpsiKE_pi_mass );
  tree_->Branch("JpsiKE_pi_q", &JpsiKE_pi_q );

  tree_->Branch("JpsiKE_pi_doca3d", &JpsiKE_pi_doca3d);
  tree_->Branch("JpsiKE_pi_doca3de", &JpsiKE_pi_doca3de);
  tree_->Branch("JpsiKE_pi_doca2d", &JpsiKE_pi_doca2d);
  tree_->Branch("JpsiKE_pi_doca2de", &JpsiKE_pi_doca2de);
  tree_->Branch("JpsiKE_pi_dz", &JpsiKE_pi_dz);
  tree_->Branch("JpsiKE_pi_near_dz", &JpsiKE_pi_near_dz);
  tree_->Branch("JpsiKE_pi_isAssociate", &JpsiKE_pi_isAssociate);
  tree_->Branch("JpsiKE_pi_pvAssociationQuality", &JpsiKE_pi_pvAssociationQuality);

  tree_->Branch("JpsiKE_pi_pdg", &JpsiKE_pi_pdg );

  tree_->Branch("JpsiKE_B_q2", &JpsiKE_B_q2 );
  tree_->Branch("JpsiKE_B_mm2", &JpsiKE_B_mm2 );
  tree_->Branch("JpsiKE_B_ptmiss", &JpsiKE_B_ptmiss );
  tree_->Branch("JpsiKE_B_Es", &JpsiKE_B_Es );
  tree_->Branch("JpsiKE_B_ptback", &JpsiKE_B_ptback );
  tree_->Branch("JpsiKE_nch", &JpsiKE_nch );
  tree_->Branch("JpsiKE_nch_filter", &JpsiKE_nch_filter );



  std::cout << "createBranch finishes" << std::endl;







}

//void NanoAnalyzer::createTitle(const std::string &name, const std::string &title) { 
// TBranch *branch = tree_->GetBranch(name.c_str()); branch->SetTitle(title.c_str());
//} 


void NanoAnalyzer::reset(void){

  std::cout << "RESET" << std::endl;

  run = -1;
  event = -1;
  luminosityBlock = -1;

  JpsiKE_e1_pt = -99;
  JpsiKE_e1_eta = -99;
  JpsiKE_e1_phi = -99;
  JpsiKE_e1_mass = -99;
  JpsiKE_e1_q = -99;
  JpsiKE_e1_vx = -99;
  JpsiKE_e1_vy = -99;
  JpsiKE_e1_vz = -99;

  JpsiKE_e2_pt = -99;
  JpsiKE_e2_eta = -99;
  JpsiKE_e2_phi = -99;
  JpsiKE_e2_mass = -99;
  JpsiKE_e2_q = -99;
  JpsiKE_e2_vx = -99;
  JpsiKE_e2_vy = -99;
  JpsiKE_e2_vz = -99;

  JpsiKE_PV_vx = -99;
  JpsiKE_PV_vy = -99;
  JpsiKE_PV_vz = -99;

  JpsiKE_bbPV_vx = -99;
  JpsiKE_bbPV_vy = -99;
  JpsiKE_bbPV_vz = -99;
  JpsiKE_bbPV_chi2 = -99;
  JpsiKE_bbPV_ndof = -99;
  JpsiKE_bbPV_rho = -99;

  JpsiKE_Jpsi_pt = -99;
  JpsiKE_Jpsi_eta = -99;
  JpsiKE_Jpsi_phi = -99;
  JpsiKE_Jpsi_mass = -99;
  JpsiKE_Jpsi_vprob = -99;
  JpsiKE_Jpsi_lip = -99;
  JpsiKE_Jpsi_lips = -99;
  JpsiKE_Jpsi_pvip = -99;
  JpsiKE_Jpsi_pvips = -99;
  JpsiKE_Jpsi_fl3d = -99;
  JpsiKE_Jpsi_fls3d = -99;
  JpsiKE_Jpsi_alpha = -99;
  JpsiKE_Jpsi_maxdoca = -99;
  JpsiKE_Jpsi_mindoca = -99;
  JpsiKE_Jpsi_vx = -99;
  JpsiKE_Jpsi_vy = -99;
  JpsiKE_Jpsi_vz = -99;

  JpsiKE_B_pt.clear();
  JpsiKE_B_eta.clear();
  JpsiKE_B_phi.clear();
  JpsiKE_B_mass.clear();
  JpsiKE_B_mcorr.clear();
  JpsiKE_B_vprob.clear();
  JpsiKE_B_lip.clear();
  JpsiKE_B_lips.clear();
  JpsiKE_B_pvip.clear();
  JpsiKE_B_pvips.clear();
  JpsiKE_B_fl3d.clear();
  JpsiKE_B_fls3d.clear();
  JpsiKE_B_alpha.clear();
  JpsiKE_B_maxdoca.clear();
  JpsiKE_B_mindoca.clear();
  JpsiKE_B_vx.clear();
  JpsiKE_B_vy.clear();
  JpsiKE_B_vz.clear();


  JpsiKE_pi_pt.clear();
  JpsiKE_pi_eta.clear();
  JpsiKE_pi_phi.clear();
  JpsiKE_pi_mass.clear();
  JpsiKE_pi_q.clear();  
  JpsiKE_pi_doca3d.clear();
  JpsiKE_pi_doca3de.clear();
  JpsiKE_pi_doca2d.clear();
  JpsiKE_pi_doca2de.clear();
  JpsiKE_pi_dz.clear();
  JpsiKE_pi_near_dz.clear();
  JpsiKE_pi_isAssociate.clear();
  JpsiKE_pi_pvAssociationQuality.clear();
  JpsiKE_pi_pdg.clear();

  JpsiKE_B_q2.clear();
  JpsiKE_B_mm2.clear();
  JpsiKE_B_ptmiss.clear();
  JpsiKE_B_Es.clear();
  JpsiKE_B_ptback.clear();

  JpsiKE_nch = -99;
  JpsiKE_nch_filter = -99;
  std::cout << "RESET end" << std::endl;

}

DEFINE_FWK_MODULE(NanoAnalyzer);
