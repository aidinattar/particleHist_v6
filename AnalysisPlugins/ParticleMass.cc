#include "ParticleMass.h"

#include "../AnalysisFramework/Event.h"
#include "../AnalysisObjects/MassMean.h"
#include "../AnalysisFramework/AnalysisInfo.h"
#include "../AnalysisFramework/AnalysisFactory.h"
#include "../AnalysisObjects/ParticleReco.h"

#include "TH1F.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

// arbitrary bin number
int nBinD = 128;

// concrete factory to create an ParticleMass analyzer
class ParticleMassFactory: public AnalysisFactory::AbsFactory {
 public:
  // assign "plot" as name for this analyzer and factory
  ParticleMassFactory(): AnalysisFactory::AbsFactory( "plot" ) {}
  // create an ParticleMass when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new ParticleMass( info );
  }
};
// create a global ParticleMassFactory, so that it is created and registered 
// before main execution start:
// when the AnalysisFactory::create function is run,
// an ParticleMassFactory will be available with name "plot".
static ParticleMassFactory er;

ParticleMass::ParticleMass( const AnalysisInfo* info ):
 AnalysisSteering( info ) {
}


ParticleMass::~ParticleMass() {
}


// function to be called at execution start
void ParticleMass::beginJob() {

  // create mass distributions for different particles
  pList.reserve( 10 );

  std::ifstream file( aInfo->value( "ranges" ).c_str() );
  std::string name;
  float min;
  float max;
  while ( file >> name >> min >> max )
    pCreate( name, min, max );

  return;

}


// function to be called at execution end
void ParticleMass::endJob() {
  // save current working area
  TDirectory* currentDir = gDirectory;
  // open histogram file
  TFile* file = new TFile( aInfo->value( "plot" ).c_str(), "RECREATE" );

  // loop over elements
  int n = pList.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    MassMean* pMean = pList[i]->pMean;
    TH1F*     hMean = pList[i]->hMean;
    // compute results
    pMean->compute();
    // get mean and rms mass
    double mean = pMean->mMeanM();
    double rms  = pMean->mRMS  ();

    // print results
    std::cout << pMean->nEvents() << ' ' 
              << mean             << ' ' 
              << rms              << std::endl;
    
    // print number of events
//    std::cout << pMean->nEvents() << std::endl;
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
void ParticleMass::update( const Event& ev ) {
  // loop over mass distributions and pass event to each of them
  unsigned int n = pList.size();
  unsigned int i;

  // get invariant mass
  static ParticleReco* pr = ParticleReco::instance();
  double Mass = pr->invariantMass();

  for ( i = 0; i < n; ++i ){
    // get mass informations
    MassMean* pMean = pList[i]->pMean;
    TH1F*     hMean = pList[i]->hMean;

    if( pMean->add( ev ) == true )
      // set center values in the graph
      // by using SetBinContent, bin count starts from 1
      hMean->Fill( Mass );
  }
  return;
}


void ParticleMass::pCreate( const std::string& name, float min, float max ) {

  // create mass distribution for events with mass in given range

  // create name for TH1F object
  const char* hName = ( "mass" + name ).c_str();

  // create TH1F and statistic objects and store their pointers
  Particle* pm = new Particle;
  pm-> name = name;
  pm->pMean = new MassMean( min, max );
  pm->hMean = new TH1F( hName, hName, nBinD, min, max );
  pList.push_back( pm );

  return;

}