#include "../interface/helperfunc.h"






float helperfunc::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){

    double maxDoca = -1.0;

    TwoTrackMinimumDistance md;
    std::vector<RefCountedKinematicParticle>::iterator in_it, out_it;

    for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
        for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
            md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
            if (md.distance() > maxDoca)
                maxDoca = md.distance();
        }
    }

    return maxDoca;
}



float helperfunc::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {

    double minDoca = 99999.9;

    TwoTrackMinimumDistance md;
    unsigned j,k,n;

    n = kinParticles.size();
    for (j = 0; j < n; j++) {
        for (k = j+1; k < n; k++) {
            md.calculate(kinParticles[j]->currentState().freeTrajectoryState(),kinParticles[k]->currentState().freeTrajectoryState());
            if (md.distance() < minDoca)
                minDoca = md.distance();
        }
    }

    return minDoca;
}




particle_cand helperfunc::calculateIPvariables(
					   AnalyticalImpactPointExtrapolator extrapolator,
					   RefCountedKinematicParticle particle,
					   RefCountedKinematicVertex vertex,
					   reco::Vertex wrtVertex
					   ){
  
    TrajectoryStateOnSurface tsos = extrapolator.extrapolate(particle->currentState().freeTrajectoryState(),
                                                             RecoVertex::convertPos(wrtVertex.position()));


    VertexDistance3D a3d;  

    std::pair<bool,Measurement1D> currentIp = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), wrtVertex);
    std::pair<bool,Measurement1D> cur3DIP = IPTools::absoluteImpactParameter(tsos, wrtVertex, a3d);
  
    // flight length
    float fl3d = a3d.distance(wrtVertex, vertex->vertexState()).value();
    float fl3de = a3d.distance(wrtVertex, vertex->vertexState()).error();
    float fls3d = -1;

    if(fl3de!=0) fls3d = fl3d/fl3de;

    // longitudinal impact parameters
    float lip = currentIp.second.value();
    float lipe = currentIp.second.error();
    float lips = -1;
  
    if(lipe!=0) lips = lip/lipe;

    // impact parameter to the PV
    float pvip = cur3DIP.second.value();
    float pvipe = cur3DIP.second.error();
    float pvips = -1;
  
    if(pvipe!=0) pvips = pvip/pvipe;

    // opening angle
    TVector3 plab = TVector3(particle->currentState().globalMomentum().x(),
                             particle->currentState().globalMomentum().y(),
                             particle->currentState().globalMomentum().z());

    const TVector3 tv3diff = TVector3(vertex->vertexState().position().x() - wrtVertex.position().x(),
                                      vertex->vertexState().position().y() - wrtVertex.position().y(),
                                      vertex->vertexState().position().z() - wrtVertex.position().z()
                                      );

    float alpha = -1;

    if(plab.Mag() != 0. && tv3diff.Mag()!=0){
        alpha = plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag());
    }

    particle_cand cand = {
        lip,
        lips,
        pvip, 
        pvips,
        fl3d,
        fls3d,
        alpha
    };


    return cand;
}




std::pair<float, float> helperfunc::calculateIPvariables(
						     RefCountedKinematicVertex tauVertex,
						     RefCountedKinematicVertex jpsiVertex,
						     reco::Vertex refitVertex
						     ){
  
  //  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(particle->currentState().freeTrajectoryState(),
  //							   RecoVertex::convertPos(jpsiVertex->vertexState().position()));
  
  VertexDistance3D a3d;   
  // flight length
  float fl3d = a3d.distance(jpsiVertex->vertexState(), tauVertex->vertexState()).value();
  float fl3de = a3d.distance(jpsiVertex->vertexState(), tauVertex->vertexState()).error();
  float fls3d = -1;
  
  if(fl3de!=0) fls3d = fl3d/fl3de;
  
  GlobalVector IPVec(tauVertex->vertexState().position().x() - jpsiVertex->vertexState().position().x(), 
		     tauVertex->vertexState().position().y() - jpsiVertex->vertexState().position().y(), 
		     tauVertex->vertexState().position().z() - jpsiVertex->vertexState().position().z());

  GlobalVector direction(jpsiVertex->vertexState().position().x() - refitVertex.position().x(), 
			 jpsiVertex->vertexState().position().y() - refitVertex.position().y(), 
			 jpsiVertex->vertexState().position().z() - refitVertex.position().z());

  double prod = IPVec.dot(direction);
  double sign = (prod >= 0) ? 1. : -1.;

  return pair<float, float>(sign*fl3d, sign*fls3d);
}


std::pair<bool, Measurement1D> helperfunc::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
							       RefCountedKinematicVertex vertex,
							       VertexDistance& distanceComputer){
  if (!tsos.isValid()) {
    return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
  }
  GlobalPoint refPoint = tsos.globalPosition();
  GlobalError refPointErr = tsos.cartesianError().position();
  GlobalPoint vertexPosition = vertex->vertexState().position();
  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
  return std::pair<bool, Measurement1D>(
					true,
					distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
}


std::pair<bool, Measurement1D> helperfunc::absoluteImpactParameter3D(const TrajectoryStateOnSurface& tsos,
								 RefCountedKinematicVertex vertex){

  VertexDistance3D dist;
  
  return absoluteImpactParameter(tsos, vertex, dist);
}


std::pair<bool, Measurement1D> helperfunc::absoluteTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
									 RefCountedKinematicVertex vertex){

  VertexDistanceXY dist;
  
  return absoluteImpactParameter(tsos, vertex, dist);
}





std::pair<bool, Measurement1D> helperfunc::signedTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
								       RefCountedKinematicVertex vertex,
								       reco::Vertex wrtVertex){
//  if (!tsos.isValid()) {
//    return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
//    }
//  GlobalPoint refPoint = tsos.globalPosition();
//  GlobalError refPointErr = tsos.cartesianError().position();
//  GlobalPoint vertexPosition = vertex->vertexState().position();
//  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
//  return std::pair<bool, Measurement1D>(
//					true,
//					distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
  
  VertexDistanceXY dist;
  
  std::pair<bool,Measurement1D> result = absoluteImpactParameter(tsos, vertex, dist);
  if (!result.first)
    return result;

  //Compute Sign
  GlobalPoint impactPoint = tsos.globalPosition();
  GlobalVector IPVec(impactPoint.x() - vertex->vertexState().position().x(), impactPoint.y() - vertex->vertexState().position().y(), 0.);
  GlobalVector direction(vertex->vertexState().position().x() - wrtVertex.position().x(), 
			 vertex->vertexState().position().y() - wrtVertex.position().y(), 0);

  double prod = IPVec.dot(direction);
  double sign = (prod >= 0) ? 1. : -1.;
  
  //Apply sign to the result
  return pair<bool, Measurement1D>(result.first, Measurement1D(sign * result.second.value(), result.second.error()));
  
}


std::pair<bool, Measurement1D> helperfunc::signedImpactParameter3D(const TrajectoryStateOnSurface& tsos,
							       RefCountedKinematicVertex vertex,
							       reco::Vertex wrtVertex){
//  if (!tsos.isValid()) {
//    return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
//    }
//  GlobalPoint refPoint = tsos.globalPosition();
//  GlobalError refPointErr = tsos.cartesianError().position();
//  GlobalPoint vertexPosition = vertex->vertexState().position();
//  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
//  return std::pair<bool, Measurement1D>(
//					true,
//					distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
  
  VertexDistance3D dist;
  
  std::pair<bool,Measurement1D> result = absoluteImpactParameter(tsos, vertex, dist);
  if (!result.first)
    return result;
  
  //Compute Sign
  GlobalPoint impactPoint = tsos.globalPosition();
  GlobalVector IPVec(impactPoint.x() - vertex->vertexState().position().x(), 
		     impactPoint.y() - vertex->vertexState().position().y(),  
		     impactPoint.z() - vertex->vertexState().position().z());

  GlobalVector direction(vertex->vertexState().position().x() - wrtVertex.position().x(), 
			 vertex->vertexState().position().y() - wrtVertex.position().y(), 
			 vertex->vertexState().position().z() - wrtVertex.position().z());

  double prod = IPVec.dot(direction);
  double sign = (prod >= 0) ? 1. : -1.;
  
  //Apply sign to the result
  return pair<bool, Measurement1D>(result.first, Measurement1D(sign * result.second.value(), result.second.error()));
  
}









math::PtEtaPhiMLorentzVector helperfunc::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
  const auto& state = fitted_children.at(i)->currentState();

  return math::PtEtaPhiMLorentzVector(
				      state.globalMomentum().perp(), 
				      state.globalMomentum().eta() ,
				      state.globalMomentum().phi() ,
				      state.mass()
				      );
}


std::tuple<float, TransientVertex> helperfunc::vertexProb( const std::vector<reco::TransientTrack>& tracks){

    float vprob = -1;
  
    KalmanVertexFitter kalman_fitter;
    TransientVertex vertex;

    try{
        vertex = kalman_fitter.vertex(tracks);
    }catch(std::exception e){
      std::cout << "No vertex found ... return" << std::endl;
      return std::forward_as_tuple(-9, vertex);
    }

    if(vertex.isValid()){

        vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());

        //    vx = vertex.position().x();
        //    vy = vertex.position().y();
        //    vz = vertex.position().z();
    
        return std::forward_as_tuple(vprob, vertex);

    }else{

        return std::forward_as_tuple(-9, vertex);

    }
}




std::tuple<Bool_t, RefCountedKinematicParticle, RefCountedKinematicVertex, RefCountedKinematicTree> helperfunc::KinematicFit(std::vector<RefCountedKinematicParticle> particles, float constrain_mass, float constrain_error){
  
  //creating the vertex fitter
  KinematicParticleVertexFitter kpvFitter;
   
  //reconstructing a J/Psi decay
  RefCountedKinematicTree tree = kpvFitter.fit(particles);
  RefCountedKinematicParticle part; // = tree->currentParticle();
  RefCountedKinematicVertex vertex; // = tree->currentDecayVertex();

  if(!tree->isEmpty() && tree->isValid() && tree->isConsistent()){

    //creating the particle fitter
    KinematicParticleFitter csFitter;
    
    // creating the constraint

    if(constrain_mass!=-1){
      //      std::cout << "Constrained fit with mass = " << constrain_mass << " error = " <<  constrain_error << std::endl;
      KinematicConstraint* constraint = new MassKinematicConstraint(constrain_mass, constrain_error);
      //the constrained fit
      tree = csFitter.fit(constraint, tree);
    } //else{
      //      std::cout << "No mass constrained fit" << std::endl;
    //    }


    //getting the J/Psi KinematicParticle
    //    std::cout <<"check" <<  tree->isEmpty() << std::endl;
    tree->movePointerToTheTop();
    part = tree->currentParticle();

    if(part->currentState().isValid()){
    
      vertex = tree->currentDecayVertex();

      if(vertex->vertexIsValid()){
      
	if(TMath::Prob(vertex->chiSquared(), vertex->degreesOfFreedom()) > 0){

	  return std::forward_as_tuple(true, part, vertex, tree);

	}
      }
    }
  }

  
  return std::forward_as_tuple(false, part, vertex, tree);

}

