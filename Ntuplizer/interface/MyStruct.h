#ifndef MyStruct_H
#define MyStruct_H

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

// kinematic fit package
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"



#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" 
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


struct particle_cand 
{
  
  // Longitudinal impact parameter and its significance
  Float_t lip;
  Float_t lips;
  
  // Impact parameter for the PV and its significance
  Float_t pvip;
  Float_t pvips;
  
  // Flight length and its significance
  Float_t fl3d;
  Float_t fls3d;
  
  // opening angle
  Float_t alpha;
  
};

struct attribute
{

  Float_t doca3d;
  Float_t doca3de;
  Float_t doca3ds;

  Float_t doca2d;
  Float_t doca2de;
  Float_t doca2ds;

//  Float_t doca1d;
//  Float_t doca1de;
//  Float_t doca1ds;

//  Bool_t isRight;

  Float_t dz;
  
  Bool_t isAssociate;
  Int_t pvAssociationQuality;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Int_t charge;
  Float_t mass;
  
  Bool_t isBdecay;
  Int_t isBdecaypdg;
  Int_t isBdecayppdg;
  Bool_t isSignal;
  Int_t nprong;
  Int_t nprong_pi0;

  Float_t near_dz; 
  Bool_t trigMatch;
  Float_t trigMatch_dr;
  Float_t dr_jpsi;
};






#endif
