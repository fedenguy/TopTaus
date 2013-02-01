#ifndef _AnalysisMonitoring_hh
#define _AnalysisMonitoring_hh
/**                                                                                                                                                                              
																			 \class    AnalysisMonitoring AnalysisMonitoring.cc "UserCode/LIP/TopTaus/src/AnalysisMonitoring.cc"                                                                     
																			 \brief    Class for monitoring the different steps of the cutflow and storing monitors into the root output file
  
  \author   Pietro Vischia, pietro.vischia@gmail.com

  \version  $Id: AnalysisMonitoring.hh,v 0.1 2013/02/01 17:53:22 vischia Exp $                                                                                                       
*/

// System headers
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string>
#include <iostream>
#include <iomanip>
#endif

// ROOT headers


#include "LIP/TopTaus/interface/MiniEvent.hh"

#include "LIP/TopTaus/interface/SampleProcessor.hh"
#include "LIP/TopTaus/interface/TauDileptonTableBuilder.hh"
#include "LIP/TopTaus/interface/HistogramBuilder.hh"

// #include <TopTauDileptonUtils.hh>



// Forward declarations

using namespace std;


class AnalysisMonitoring :  public TauDileptonTableBuilder , public HistogramBuilder, public SampleProcessor {
  
public:
  
  AnalysisMonitoring(double tauPtCut) :  TauDileptonTableBuilder(), HistogramBuilder(), SampleProcessor(tauPtCut) { outFile_=0; }
  
  void buildStructures();
  void buildTrees(TFile*);
  
  virtual ~AnalysisMonitoring(){  infoFile_.close(); }

  // FIXME: check whether that is true, and in this case add systematics on tau id  
  ////////////////////////////////////////////////////////////////
  // WARNING : Important note :
  // This script does not take into account systematics on tau id  
  ////////////////////////////////////////////////////////////////
  TGraphAsymmErrors*  scriptForPlotWithSystematics( TString, int, string);

  //TauDileptonHistogramBuilder * tauDileptonHisto_;
  
  map<TString,TTree*> treeMap_;
  
  map<TString,TDirectory*> dirsMaps_;                           // map with dirs to store
  map< TString, SelectionMonitor* >  tauDileptonStepsMons_;     // object monitors at different steps of tau dilepton selection 
  map< TString, SelectionMonitor* >  tauDileptonYieldsMons_;    // yield monitors at different steps of tau dilepton selection 
  map< TString, SelectionMonitor* >  tauDileptonMCYieldsMons_;  // yield monitors at different steps of tau dilepton selection 
  map< TString, SelectionMonitor* >  tauDileptonBkgMons_;       // yield monitors at different steps of tau dilepton selection
  map< TString, SelectionMonitor* >  tauDileptonDebugMons_;     // yield monitors at different steps of tau dilepton selection 

};

#endif
