#include "ParticleLifetime.h"

#include "../AnalysisFramework/Event.h"
#include "../AnalysisObjects/LifetimeFit.h"
#include "../AnalysisFramework/AnalysisInfo.h"
#include "../AnalysisFramework/AnalysisFactory.h"
#include "../AnalysisObjects/ProperTime.h"

#include "TH1F.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>


// arbitrary bin number
int nBin = 128;

// concrete factory to create an ParticleLifetime analyzer
class ParticleLifetimeFactory: public AnalysisFactory::AbsFactory {
 public:
  // assign "time" as name for this analyzer and factory
  ParticleLifetimeFactory(): AnalysisFactory::AbsFactory( "time" ) {}
  // create an ParticleLifetime when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new ParticleLifetime( info );
  }
};
// create a global ParticleLifetimeFactory, so that it is created and registered 
// before main execution start:
// when the AnalysisFactory::create function is run,
// an ParticleLifetimeFactory will be available with name "time".
static ParticleLifetimeFactory er;

ParticleLifetime::ParticleLifetime( const AnalysisInfo* info ):
 AnalysisSteering( info ) {
}


ParticleLifetime::~ParticleLifetime() {
}


// function to be called at execution start
void ParticleLifetime::beginJob() {

  // create mass distributions for different particles
  pList.reserve( 10 );
  std::ifstream file( "particleFitters" );
  std::string name;
  float   min,  max;
  double tmin, tmax;
  double smin, smax,
               step;
  while ( file >> name 
               >> min  >>  max 
               >> tmin >> tmax
               >> smin >> smax 
               >> step         )
    pCreate( name,  min,  max,
                   tmin, tmax,
                   smin, smax,
                         step );

  return;

}


// function to be called at execution end
void ParticleLifetime::endJob() {
  // save current working area
  TDirectory* currentDir = gDirectory;
  // open histogram file
  TFile* file = new TFile( aInfo->value( "time" ).c_str(), "RECREATE" );

  // loop over elements
  int n = pList.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    LifetimeFit* pMean = pList[i]->pMean;
    TH1F*        hMean = pList[i]->hMean;
    // compute results
    pMean->compute();

    double mean = pMean->lifeTime     ();
    double  err = pMean->lifeTimeError();

    // print dacay time
    std::cout << pMean->nEvents() << ' '
              << mean             << ' '
              << err              << std::endl;

    // save distribution
    hMean->Write();
  }

  // close file
  file->Close();
  delete file;
  // restore working area
  currentDir->cd();

  return;

}


// function to be called for each event
void ParticleLifetime::update( const Event& ev ) {
  // loop over mass distributions and pass event to each of them
  unsigned int n = pList.size();
  unsigned int i;

  // get invariant mass
  static ProperTime*   pt = ProperTime::instance();
  double             time =        pt->decayTime();

  for ( i = 0; i < n; ++i ){
    // get mass informations
    LifetimeFit* pMean = pList[i]->pMean;
    TH1F*        hMean = pList[i]->hMean;

    if( pMean->add( ev ) == true ){
      
      // set center values in the graph
      // by using SetBinContent, bin count starts from 1
      hMean->Fill( time );
    }
  }
  return;
}


void ParticleLifetime::pCreate( const std::string&        name,
                                float      min, float      max,
                                double timeMin, double timeMax,
                                double scanMin, double scanMax,
                                double scanstp ){

  // create mass distribution for events with mass in given range

  // create name for TH1F object
  const char* hName = ( "time" + name ).c_str();

  // create TH1F and statistic objects and store their pointers
  Particle* pm = new Particle;
  pm-> name = name;
  pm->pMean = new LifetimeFit(     min,     max,
                               timeMin, timeMax,
                               scanMin, scanMax,
                               scanstp          );
  pm->hMean = new TH1F( hName, hName, nBin, timeMin, timeMax );
  pList.push_back( pm );

  return;

}