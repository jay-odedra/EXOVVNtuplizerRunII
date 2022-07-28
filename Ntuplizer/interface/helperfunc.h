#ifndef HELPERFUNC_H
#define HELPERFUNC_H

//#include "../interface/CandidateNtuplizer.h"
#include "../interface/MyStruct.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

/////#include "DataFormats/VertexReco/interface/VertexFwd.h"
/////#include "DataFormats/VertexReco/interface/Vertex.h"
///////#include "DataFormats/PatCandidates/interface/Jet.h"
/////#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
/////#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
/////#include "DataFormats/Candidate/interface/Candidate.h"
/////
/////
/////#include "JetMETCorrections/Modules/interface/JetResolution.h"
///////#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
/////
/////#include "DataFormats/BeamSpot/interface/BeamSpot.h"
/////
/////// for vertexing
/////#include "PhysicsTools/PatUtils/interface/TriggerHelperfunc.h"
/////#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
/////
///////#include "FWCore/Framework/interface/ESHandle.h"
/////#include "FWCore/MessageLogger/interface/MessageLogger.h"
/////#include "FWCore/Common/interface/TriggerNames.h"
/////
/////#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
/////#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"   
/////#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
/////#include "DataFormats/TrackReco/interface/DeDxData.h" 
/////#include "DataFormats/MuonReco/interface/MuonSelectors.h"
/////#include "DataFormats/MuonReco/interface/MuonQuality.h"
/////
/////#include "DataFormats/Math/interface/deltaR.h"
///////#include "DataFormats/Common/interface/TriggerResults.h"
///////#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
///////#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
/////
/////
/////
/////#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
/////#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
/////#include "DataFormats/L1Trigger/interface/BXVector.h"
/////#include "DataFormats/L1Trigger/interface/Jet.h"
/////#include "DataFormats/L1Trigger/interface/EGamma.h"
/////#include "DataFormats/L1Trigger/interface/Muon.h"
/////#include "DataFormats/L1Trigger/interface/Tau.h"
/////#include "DataFormats/L1Trigger/interface/L1Candidate.h"
/////#include "DataFormats/HLTReco/interface/TriggerEvent.h"
/////#include "DataFormats/HLTReco/interface/TriggerObject.h"
/////#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
/////#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
/////#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
/////#include "DataFormats/HLTReco/interface/EgammaObject.h"
/////#include "DataFormats/HLTReco/interface/EgammaObjectFwd.h"
/////
/////
/////
/////
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"             
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
/////
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
/////
/////
/////// kinematic fit package
/////#include "DataFormats/TrackReco/interface/Track.h"
/////#include "DataFormats/TrackReco/interface/TrackFwd.h"
/////#include "FWCore/Framework/interface/MakerMacros.h"
/////
//////// Run3 
/////#include "FWCore/Framework/interface/ConsumesCollector.h"
/////#include "FWCore/Framework/interface/one/EDAnalyzer.h"
///////#include "FWCore/Framework/interface/Event.h"
/////#include "FWCore/Framework/interface/ESConsumesCollector.h"
/////#include "FWCore/Framework/interface/Frameworkfwd.h"
/////
/////#include "FWCore/Framework/interface/Event.h"
/////#include "FWCore/Framework/interface/MakerMacros.h"
/////
/////#include "FWCore/ParameterSet/interface/ParameterSet.h"
/////
/////#include "FWCore/Framework/interface/EventSetup.h"
///////#include "FWCore/Framework/interface/ESHandle.h"
/////
/////
/////
/////
/////
/////
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
/////
/////
/////#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
/////
/////#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
/////
/////#include "FWCore/ParameterSet/interface/FileInPath.h"
/////#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
/////
/////
/////#include <tuple>
/////#include <sstream>
/////#include <iostream>
/////#include <cstdlib>
/////#include <fstream>
/////#include <string>
/////#include <algorithm> 
/////#include <stdio.h>
/////#include <stdlib.h>
/////#include <vector>
/////#include <map>
/////
///////#include "Hammer/Hammer.hh"
///////#include "Hammer/Process.hh"
///////#include "Hammer/Particle.hh"
/////
/////#include "FWCore/Framework/interface/global/EDAnalyzer.h"
/////#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
/////#include "DataFormats/Provenance/interface/ModuleDescription.h"
/////#include "DataFormats/TestObjects/interface/ToyProducts.h"
/////#include "FWCore/Framework/interface/Event.h"
/////#include "FWCore/Framework/interface/EventSetup.h"
/////#include "FWCore/Framework/interface/ESHandle.h"
/////#include "FWCore/Integration/interface/ESTestData.h"
/////#include "FWCore/Integration/interface/ESTestRecords.h"
/////#include "FWCore/MessageLogger/interface/MessageLogger.h"
/////#include "FWCore/Framework/interface/MakerMacros.h"
/////#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
/////#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
/////#include "FWCore/Utilities/interface/ESGetToken.h"



//#include "../interface/MyStruct.h"
//#include <iostream>
//#include <string.h>
#include <TVector3.h>
//#include <TRandom.h>

using namespace std;

// mass


class helperfunc{
  
 public:

  ParticleMass electron_mass = 0.000511;
  ParticleMass muon_mass = 0.1056583;
  ParticleMass jpsi_mass = 3.09687;
  ParticleMass pion_mass = 0.139571;
  ParticleMass kaon_mass = 0.493677;
  ParticleMass ds_mass = 2.01026;
  ParticleMass d0_mass = 1.86483;
  
  float electron_sigma = 0.0000001;
  float muon_sigma = 0.0000001;
  float jp_m_sigma = 0.00004;
  float pion_sigma = 0.000016;
  float kaon_sigma = 0.000016;
  float phi_sigma = 0.000016;
  float ds_sigma = 0.00005;
  float d0_sigma = 0.00005;
  
  float mass_kaon = 0.493677;
  float mass_pion = 0.139571;
  float mass_D0 = 1.86483;
  float mass_Dstar = 2.01026;
  float mass_B0 = 5.27963;
  float mass_Bc = 6.2749;
  float mass_electron = 0.000511;

  


  // return vertex x, y, and z
  // return tau decay mode
  // max. or min. of distance of closest approach
  float getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  float getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  
  // nagivator of the kinematic fit class
  //  void printout(const RefCountedKinematicVertex& myVertex);
  //  void printout(const RefCountedKinematicParticle& myParticle);
  //  void printout(const RefCountedKinematicTree& myTree);
  

  // calculate IP related variables
  particle_cand calculateIPvariables(AnalyticalImpactPointExtrapolator extrapolator,
				     RefCountedKinematicParticle particle,
				     RefCountedKinematicVertex vertex,
				     reco::Vertex wrtVertex);

  std::pair<float, float> calculateIPvariables(RefCountedKinematicVertex tauVertex,
					       RefCountedKinematicVertex jpsiVertex,
					       reco::Vertex refitVertex);

  // absolute impact parameter
  std::pair<bool, Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
							 RefCountedKinematicVertex vertex,
							 VertexDistance& distanceComputer);

  std::pair<bool, Measurement1D> absoluteImpactParameter3D(const TrajectoryStateOnSurface& tsos,
							   RefCountedKinematicVertex vertex);


  std::pair<bool, Measurement1D> absoluteTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
								   RefCountedKinematicVertex vertex);
  

  std::pair<bool, Measurement1D> signedTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
								 RefCountedKinematicVertex vertex,
								 reco::Vertex wrtVertex);

  std::pair<bool, Measurement1D> signedImpactParameter3D(const TrajectoryStateOnSurface& tsos,
							 RefCountedKinematicVertex vertex,
							 reco::Vertex wrtVertex);


  math::PtEtaPhiMLorentzVector daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i);
  
  // Vertex probability calculater
  std::tuple<float, TransientVertex> vertexProb( const std::vector<reco::TransientTrack>& tracks);


  std::tuple<Bool_t, RefCountedKinematicParticle, RefCountedKinematicVertex, RefCountedKinematicTree> KinematicFit(std::vector<RefCountedKinematicParticle> particles, float constrain_mass, float constrain_error);

  //  bool basicPFcut(pat::PackedCandidate pf);

};





#endif
