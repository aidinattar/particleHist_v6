#include "ProperTime.h"
#include "../AnalysisFramework/Event.h"
#include "../AnalysisUtilities/Utilities.h"
#include "../AnalysisUtilities/Constants.h"
#include "ParticleReco.h"

#include <iostream>
#include <math.h>

using namespace std;



ProperTime::ProperTime() {
}


ProperTime::~ProperTime() {
}


// recompute tag informations for new event
void ProperTime::update( const Event& ev ) {
  Constants cost;
  
  // compute distance from the origin
  double d = sqrt( pow( ev.getX(), 2 ) +
                   pow( ev.getY(), 2 ) +
                   pow( ev.getZ(), 2 ) );

  // get invariant mass and energy
  static ParticleReco* pr = ParticleReco::instance();
  double m = pr->invariantMass();
  double e = pr->  totalEnergy();

  // get momenta
  double p = sqrt( e * e - m * m );

  // compute decay time
  time = d * m / ( p * cost.lightVelocity );
  
  return;

}


double ProperTime::decayTime(){
  // check for new event and return result
  check();
  return time;

}