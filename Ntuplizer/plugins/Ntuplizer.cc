#include "../interface/Ntuplizer.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/METsNtuplizer.h"
#include "../interface/PileUpNtuplizer.h"
#include "../interface/GenEventNtuplizer.h"
#include "../interface/GenParticlesNtuplizer.h"
#include "../interface/VerticesNtuplizer.h"
#include "../interface/BsTauTauNtuplizer.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"

// #include "DataFormats/METReco/interface/PFMET.h"


///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
 

  beamToken_                  (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  vtxToken_             	    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_             	    (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  packedpfcandidatesToken_    (consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates"))), 
  svToken_                    (consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("SecondaryVertices"))), 
  puinfoToken_          	    (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  geneventToken_        	    (consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),     
  lheEventProductToken_       (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("externallheProducer"))),     
  genparticleToken_     	    (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
  gentauToken_     	          (consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("gentaus"))),
  
  muonToken_	      	        (consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  CaloTowerCollection_        (consumes<edm::SortedCollection<CaloTower>>(edm::InputTag("towerMaker"))),
  
  
  triggerToken_	      	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"))),
  triggerObjects_	      	    (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects")))
{


  nevents = 0;
  nevents_all = 0;

  /*=======================================================================================*/
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
 
  //std::map< std::string, bool > runFlags;
  runFlags["runOnMC"] = iConfig.getParameter<bool>("runOnMC");
//  runFlags["useHammer"] = iConfig.getParameter<bool>("useHammer");
  runFlags["doGenParticles"] = iConfig.getParameter<bool>("doGenParticles");
  runFlags["doGenEvent"] = iConfig.getParameter<bool>("doGenEvent");
  runFlags["doPileUp"] = iConfig.getParameter<bool>("doPileUp");
  runFlags["doVertices"] = iConfig.getParameter<bool>("doVertices");
  runFlags["doBsTauTau"] = iConfig.getParameter<bool>("doBsTauTau");
  runFlags["doGenHist"] = iConfig.getParameter<bool>("doGenHist");
  runFlags["isTruth"] = iConfig.getParameter<bool>("isTruth");
  runFlags["verbose"] = iConfig.getParameter<bool>("verbose");
  runValues["dzcut"] = iConfig.getParameter<double>("dzcut");
  runValues["fsigcut"] = iConfig.getParameter<double>("fsigcut");
  runValues["vprobcut"] = iConfig.getParameter<double>("vprobcut");
  runValues["tau_charge"] = iConfig.getParameter<unsigned int>("tau_charge");

  std::cout << "Ntuplizer: (dzcut, fsigcut, vprobcut, tau_charge) = " << runValues["dzcut"] << " " << runValues["fsigcut"] << " " << runValues["vprobcut"] << " " << runValues["tau_charge"] << std::endl;
  
//  std::cout << "useammer" << runFlags["useHammer"] << std::endl;

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  jecpath = "EXOVVNtuplizerRunII/Ntuplizer/data/" + jecpath;
  //  std::cout << "jecpath  "<< jecpath  <<std::endl;
 
  nBranches_ = new NtupleBranches( runFlags, tree );
  
  /*=======================================================================================*/
  /* Histogram buildinng, definition in NtupleBrances */
  /* Histogram for cutflow */


//  if(runFlags["useHammer"]){
//    nBranches_->hammer_width = fs->make<TH1F>("hammer_width", "Hammer width", 24, 0, 24);  
//  }

  /* Histogram for genParticles */ 
  if (runFlags["doGenHist"]){
  } 
  //  if (runFlags["doGenHist"]) {
  nBranches_-> LabelHistograms( runFlags );
      //  }

  /*=======================================================================================*/


 
  /*=======================================================================================*/  
  std::vector<edm::EDGetTokenT<reco::VertexCollection>> vtxTokens;
  vtxTokens.push_back( vtxToken_  );  


  /*=======================================================================================*/  

						      
  if (runFlags["doVertices"]) {
    nTuplizers_["vertices"] = new VerticesNtuplizer( vtxTokens   , 
						     beamToken_,
                                                     nBranches_  ,
						     runFlags    );
  }

  if (runFlags["doBsTauTau"]) {
    std::cout<<"\n\n --->GETTING INSIDE doBsTauTau<---\n\n"<<std::endl;
    nTuplizers_["BsTauTau"] = new BsTauTauNtuplizer( muonToken_   , 
						     CaloTowerCollection_ ,
						     vtxToken_   , 
						     packedpfcandidatesToken_,
						     triggerToken_,
						     genparticleToken_,
						     gentauToken_,
						     runFlags,
						     runValues,
						     runStrings,
						     nBranches_ );
  }


 
  /*=======================================================================================*/    
  if ( runFlags["runOnMC"] ){

     
    if (runFlags["doGenParticles"]) {
      std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> genpTokens;
      genpTokens.push_back( genparticleToken_ );

      nTuplizers_["genParticles"] = new GenParticlesNtuplizer( genpTokens, nBranches_, runFlags );
    }

    if (runFlags["doPileUp"]) {
      std::vector<edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > puTokens;
      puTokens.push_back( puinfoToken_ );
      nTuplizers_["PU"] = new PileUpNtuplizer( puTokens, nBranches_, runFlags );
    }

    if (runFlags["doGenEvent"]) {
      std::vector<edm::EDGetTokenT< GenEventInfoProduct > > geneTokens;
      geneTokens.push_back( geneventToken_ );
      std::vector<edm::EDGetTokenT<  LHEEventProduct > > lheTokens;
      lheTokens.push_back( lheEventProductToken_);
      nTuplizers_["genEvent"] = new GenEventNtuplizer( geneTokens, nBranches_ , lheTokens, runFlags);
    }
  }


  // if (runFlags["doGenHist"] || runFlags["doGenHist"]){
  //     LabelHistogram(nBranches_, runFlags );


  // }
}

///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::~Ntuplizer()
{
	  
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){
    //std::cout << "deconstructor: Branches: " << it->first << std::endl;
    delete it->second;
  }
   nTuplizers_.clear();
   
   delete nBranches_;
   
   std::cout << "Total number of filled events = " << nevents << "/" << nevents_all << ", eff = " << Float_t(nevents)/Float_t(nevents_all) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  nBranches_->reset();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if( vertices->empty() ) return; // skip the event if no PV found
  
  nBranches_->EVENT_event     = iEvent.id().event();
  nBranches_->EVENT_run       = iEvent.id().run();
  nBranches_->EVENT_lumiBlock = iEvent.id().luminosityBlock();  

  //  std::cout<<" ----------------- before the branches loop"<<std::endl; 

  bool isSave = true;
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){

    isSave = (it->second)->fillBranches( iEvent, iSetup );
    //    std::cout << "Fill Branchines: " << it->first << ", result = " << isSave << std::endl;
    if(!isSave) break;
  }

  //  std::cout << "isSave =  " << isSave << std::endl;
  if(isSave){
    //    std::cout << "-------------- save -----------" << std::endl;
    nBranches_->fillTree();
    nevents++;
  }
  
  nBranches_->reset();    
 
  nevents_all++;
}






///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginJob(){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endJob( ) {
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
