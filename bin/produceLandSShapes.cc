/**                                                                                                                                                                              
  \class    produceLandSShapes produceLandSShapes.cc "UserCode/LIP/TopTaus/bin/produceLandSShapes.cc"                                                                     
  \brief    executable for producing the shapes for limit estimation with landS, following the instructions at http://cmsdoc.cern.ch/~anikiten/cms-higgs/LandsInstructions1.txt
  
  \author   Pietro Vischia

  \version  $Id: produceLandSShapes.cc,v 1.5 2013/04/19 14:46:45 vischia Exp $                                                                                                       
*/

#include "LIP/TopTaus/interface/LandSShapesProducer.hh"

// System includes
#include <string>

// ROOT includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TSystem.h"
/*

@run me like this

tauDileptonPDFBuilderFitter mutau HIGGS2BKG
tauDileptonPDFBuilderFitter mutau HIGGS2BKG

tauDileptonPDFBuilderFitter mutau HIGGS2BKG
tauDileptonPDFBuilderFitter mutau HIGGS3BKG

*/

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 4 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py --produceOnly [true|false]" << std::endl;
    return 0;
  }
  
  string sProduceOnly("");
  bool produceOnly(false);
  // Get input arguments
  for(int i=2;i<argc;i++){
    string arg(argv[i]);
    //    if(arg.find("--help")        !=string::npos) { printHelp(); return -1;}   
    if(arg.find("--produceOnly") !=string::npos) { sProduceOnly = argv[i+1];}
    //check arguments // FIXME: implement --blah bih
  }
  
  if(sProduceOnly == "true")
    produceOnly=true;
  else if(sProduceOnly == "false")
    produceOnly=false;
  else{
    cout << "Error. ProduceOnly value not defined. Defaulting to true" << endl;
    produceOnly=true;
  }

  string parSet(argv[1]);

  LandSShapesProducer* myProducer = new LandSShapesProducer(parSet, produceOnly);
  myProducer->Produce();
  cout << "Shapes producer reached its natural end" << endl;

  return 0;
  
}
