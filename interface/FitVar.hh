#ifndef _FitVar_hh
#define _FitVar_hh
/**                                                                                                                                                                              
																						 \class    likelihoodFitter::tauDileptonPDFBuilderFitter tauDileptonPDFBuilderFitter.cc "UserCode/LIP/TopTaus/tauDileptonPDFBuilderFitter.cc"                                                                     
  \brief    executable for performing multivariable likelihood fit in order to improve estimation of N_fakes
  
  \author   Pietro Vischia

  \version  $Id: TauDileptonPDFBuilderFitter.hh,v 1.6 2012/09/18 17:12:33 vischia Exp $                                                                                                       
*/

//#if !defined(__CINT__) || defined(__MAKECINT__)
//
//#ifndef __CINT__
////#include "RooGlobalFunc.h"
//#endif

// System headers
#include <string>


// ROOT headers
#include "TString.h"

//#endif

class FitVar{
  
public:
  FitVar(std::string, double, double, double, double, double, Int_t, Int_t);
  std::string getVarName();
  double getMin();
  double getMax();
  int getBins();
  double getHmin();
  double getHmax();
  Int_t getUnbinned();
  Int_t getSmoothOrder();
private:
  std::string varName_;
  double min_;
  double max_;
  double bins_;
  double hmin_;
  double hmax_;
  Int_t unbinned_;
  Int_t smoothOrder_;
};

#endif //_FitVar_hh

