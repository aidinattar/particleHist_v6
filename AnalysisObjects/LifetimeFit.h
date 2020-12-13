// Created by Aidin Attar. 2019

// LifetimeFit.h

#ifndef LifetimeFit_h
#define LifetimeFit_h

#include <vector>

class Event;

class LifetimeFit {

 public:

  LifetimeFit( float  mmin, float  mmax, // mass range
               double tmin, double tmax, // time range
               double smin, double smax, // scan range
               double step);             // scan  step
  ~LifetimeFit();

  bool add( const Event& ev ); 
  void compute(); // compute mean and rms

  int    nEvents(); // return number of accepted events

  double lifeTime(); // return lifetime
  double lifeTimeError(); // return lifetime error

 private:

  double massMin; // min mass
  double massMax; // max mass

  double timeMin; // min time
  double timeMax; // max time

  double scanMin; // min scan
  double scanMax; // max scan
  double scanstp; // scan step

  double timemea;
  double timeerr;

  std::vector <double> decay_times; // all decay times
  

};

#endif