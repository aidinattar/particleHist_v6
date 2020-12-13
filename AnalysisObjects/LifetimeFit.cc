// Created by Aidin Attar. 2019

// LifetimeFit.cc

#include "LifetimeFit.h"
#include "../AnalysisFramework/Event.h"
#include "ParticleReco.h"
#include "ProperTime.h"

#include <cmath>
#include <vector>

#include "../AnalysisUtilities/QuadraticFitter.h"


LifetimeFit::LifetimeFit( float  mmin, float  mmax,
                          double tmin, double tmax,
                          double smin, double smax,
                          double step ):
  // initializations
  massMin( mmin ), massMax( mmax ),
  timeMin( tmin ), timeMax( tmax ),
  scanMin( smin ), scanMax( smax ),
  scanstp( step ),
  timemea(    0 ), timeerr(    0 ) {
}


LifetimeFit::~LifetimeFit() {
}


// add data from a new event
bool LifetimeFit::add( const Event& ev ) {
  // get invariant mass
  static ParticleReco* pr = ParticleReco::instance();
  double             mass = pr->invariantMass();

  // get invariant massnEvt;
  static ProperTime* pt = ProperTime::instance();
  double           time = pt->decayTime();

  // check for mass being in range
  if ( mass < massMin || mass > massMax ) 
    return false;

  // update decay times
  decay_times.push_back( time );

  return true;

}


void LifetimeFit::compute() {
  QuadraticFitter* qf = new QuadraticFitter;
  double t_s;
  int i;
  int n = decay_times.size();
  double L_t;
  double t_i;

  for ( t_s = scanMin; t_s <= scanMax; t_s += scanstp ) {
    L_t = 0;
    
    for ( i = 0; i < n; ++i ) {
      t_i  = decay_times[i];
      L_t += t_i / t_s                  + log( t_s )             +
             log( exp( -timeMin / t_s ) - exp( -timeMax / t_s ) );
    }
    qf->add( t_s, L_t );
  }

//  double a = qf->a();
  double b = qf->b();
  double c = qf->c();

  timemea = -b /     ( 2 * c );
  timeerr =  1 / sqrt( 2 * c );

  return;
}


// return number of selected events
int LifetimeFit::nEvents() {
  return decay_times.size();
}


// return mean of decay times
double LifetimeFit::lifeTime() {
  return timemea;
}


// return error of decay time mean
double LifetimeFit::lifeTimeError() {
  return timeerr;
}