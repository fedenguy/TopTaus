/**   
      \class	    physicsAnalysis physicsAnalysis.cc "UserCode/LIP/TopTaus/bin/physicsAnalysis.cc"                                                                     
      \brief    executable for performing the cutflow analysis and control plots for the various channels and data samples. Must expand to calling other classes
      
      \author   Pietro Vischia
      
      \version  $Id: physicsAnalysis.cc,v 1.22 2013/04/19 15:13:00 vischia Exp $                                                                                                       
*/

#include "LIP/TopTaus/interface/CutflowAnalyzer.hh"
#include "LIP/TopTaus/interface/HistogramPlotter.hh"

// System includes	
#include <string>

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// ROOT includes
#include "TSystem.h"



/*
  
@run me like this

physicsAnalysis test/physicsAnalysisParSets_cfg.py sample --> analyze sample

physicsAnalysis test/physicsAnalysisParSets_cfg.py doPlots --> draw plots

physicsAnalysis test/physicsAnalysisParSets_cfg.py doTables --> draw tables

*/

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //check arguments
  if ( argc < 3 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py action" << std::endl;
    return 0;
  }
  
  
  //  // Get input arguments
  //  string sDoDatacards("");
  //  bool doDatacards(false);
  
  //  for(int i=3;i<argc;i++){
  //    string arg(argv[i]);
  //    //    if(arg.find("--help")        !=string::npos) { printHelp(); return -1;}                                                                                                                                  
  //    if(arg.find("--doDatacards") !=string::npos) { sDoDatacards = argv[i+1];}
  //    //check arguments // FIXME: implement --blah bih                             
  //  }
  //  if(sDoDatacards == "true")
  //    doDatacards=true;
  //  else if(sDoDatacards == "false")
  //    doDatacards=false;
  //  else{
  //    cout << "Error. DoDatacards value not defined. Defaulting to false (tables mode)" << endl;
  //    doDatacards=false;
  //  }
  
  
  
  string parSet(argv[1]);
  string runOn(argv[2]);
  
  const edm::ParameterSet &pSet = edm::readPSetsFrom(parSet)->getParameter<edm::ParameterSet>("PhysicsAnalysisParSet");
  
  double tauPtCut = pSet.getParameter<double>("tauPtCut");
  bool noUncertainties = pSet.getParameter<bool>("noUncertainties");
  bool doWPlusJetsAnalysis = pSet.getParameter<bool>("doWPlusJetsAnalysis");
  TString inputArea = TString(pSet.getParameter<string>("inputArea"));
  TString outputArea = TString(pSet.getParameter<string>("outputArea"));
  TString puFileName = TString(pSet.getParameter<string>("puFileName"));
  TString runRange   = TString(pSet.getParameter<string>("runRange"));
  std::vector<double> brHtaunu = pSet.getParameter<std::vector<double> >("brHtaunu");
  std::vector<double> brHtb    = pSet.getParameter<std::vector<double> >("brHtb");
  
  CutflowAnalyzer* analyzer = new CutflowAnalyzer( tauPtCut, noUncertainties, doWPlusJetsAnalysis, inputArea, outputArea, puFileName, runRange, brHtaunu, brHtb /*parSet*/ );
  
  std::cout << "Analyzer has been set with a cut on tau pt of " << tauPtCut << " GeV/c " << std::endl;
  
  
  if(runOn ==      "data_muonA")          analyzer->process_data_RunA()         ;
  else if(runOn == "data_muonAr")         analyzer->process_data_RunAr()        ;
  else if(runOn == "data_muonB_1")        analyzer->process_data_RunB(1)         ;
  else if(runOn == "data_muonB_2")        analyzer->process_data_RunB(2)         ;
  else if(runOn == "data_muonB_3")        analyzer->process_data_RunB(3)         ;
  else if(runOn == "data_muonB_4")        analyzer->process_data_RunB(4)         ;
  else if(runOn == "data_muonB_5")        analyzer->process_data_RunB(5)         ;
  else if(runOn == "data_muonB_6")        analyzer->process_data_RunB(6)         ;
  else if(runOn == "data_muonB_7")        analyzer->process_data_RunB(7)         ;
  else if(runOn == "data_muonB_8")        analyzer->process_data_RunB(8)         ;
  else if(runOn == "data_muonB_9")        analyzer->process_data_RunB(9)         ;
  else if(runOn == "data_muonB_10")       analyzer->process_data_RunB(10)         ;
  else if(runOn == "data_muonB_11")       analyzer->process_data_RunB(11)         ;
  else if(runOn == "data_muonB_12")       analyzer->process_data_RunB(12)         ;
  else if(runOn == "data_muonB_13")       analyzer->process_data_RunB(13)         ;
  else if(runOn == "data_muonC1")         analyzer->process_data_RunC1()        ;
  else if(runOn == "data_muonC2_1")       analyzer->process_data_RunC2(1)        ;
  else if(runOn == "data_muonC2_2")       analyzer->process_data_RunC2(2)        ;
  else if(runOn == "data_muonC2_3")       analyzer->process_data_RunC2(3)        ;
  else if(runOn == "data_muonC2_4")       analyzer->process_data_RunC2(4)        ;
  else if(runOn == "data_muonC2_5")       analyzer->process_data_RunC2(5)        ;
  else if(runOn == "data_muonC2_6")       analyzer->process_data_RunC2(6)        ;
  else if(runOn == "data_muonC2_7")       analyzer->process_data_RunC2(7)        ;
  else if(runOn == "data_muonC2_8")       analyzer->process_data_RunC2(8)        ;
  else if(runOn == "data_muonC2_9")       analyzer->process_data_RunC2(9)        ;
  else if(runOn == "data_muonC2_10")      analyzer->process_data_RunC2(10)        ;
  else if(runOn == "data_muonC2_11")      analyzer->process_data_RunC2(11)        ;
  else if(runOn == "data_muonC2_12")      analyzer->process_data_RunC2(12)        ;
  else if(runOn == "data_muonC2_13")      analyzer->process_data_RunC2(13)        ;
  else if(runOn == "data_muonD_1")        analyzer->process_data_RunD(1)         ;
  else if(runOn == "data_muonD_2")        analyzer->process_data_RunD(2)         ;
  else if(runOn == "data_muonD_3")        analyzer->process_data_RunD(3)         ;
  else if(runOn == "data_muonD_4")        analyzer->process_data_RunD(4)         ;
  else if(runOn == "data_muonD_5")        analyzer->process_data_RunD(5)         ;
  else if(runOn == "data_muonD_6")        analyzer->process_data_RunD(6)         ;
  else if(runOn == "data_muonD_7")        analyzer->process_data_RunD(7)         ;
  else if(runOn == "data_muonD_8")        analyzer->process_data_RunD(8)         ;
  else if(runOn == "data_muonD_9")        analyzer->process_data_RunD(9)         ;
  else if(runOn == "data_muonD_10")        analyzer->process_data_RunD(10)         ;
  else if(runOn == "data_muonD_11")        analyzer->process_data_RunD(11)         ;
  else if(runOn == "data_muonD_12")        analyzer->process_data_RunD(12)         ;
  else if(runOn == "data_muonD_13")        analyzer->process_data_RunD(13)         ;
  else if(runOn == "dibosons")            analyzer->process_dibosons()           ;
  else if(runOn == "hh_higgs_bychannel")  analyzer->process_hh_higgs_bychannel() ;
  else if(runOn == "hh_higgs")            analyzer->process_hh_higgs()           ;
  else if(runOn == "qcd")                 analyzer->process_qcd()                ;
  else if(runOn == "singletop")           analyzer->process_singletop()          ;
  //  else if(runOn == "trigger_mht")         analyzer->process_trigger_mht()        ;
  //  else if(runOn == "trigger_mhtjets")     analyzer->process_trigger_mhtjets()    ;
  else if(runOn == "ttbar_mutau")         analyzer->process_ttbar_mutau()    ;
  else if(runOn == "ttbar_mumu")          analyzer->process_ttbar_mumu()    ;
  else if(runOn == "ttbar_emu")           analyzer->process_ttbar_emu()    ;
  else if(runOn == "ttbar_ddbkg")         analyzer->process_ttbar_ddbkg()    ;
  else if(runOn == "ttbar_mcbkg")         analyzer->process_ttbar_mcbkg()    ;
  else if(runOn == "ttbar_unc")           analyzer->process_ttbar_unc()          ;
  else if(runOn == "ttbar")               analyzer->process_ttbar()              ;
  else if(runOn == "wh_higgs_bychannel")  analyzer->process_wh_higgs_bychannel() ;
  else if(runOn == "wh_higgs")            analyzer->process_wh_higgs()           ;
  else if(runOn == "tbh_higgs_bychannel") analyzer->process_tbh_higgs_bychannel() ;
  else if(runOn == "tbh_higgs")           analyzer->process_tbh_higgs()           ;
  else if(runOn == "htb_higgs_bychannel") analyzer->process_htb_higgs_bychannel() ;
  else if(runOn == "htb_higgs")           analyzer->process_htb_higgs()           ;
  else if(runOn == "wjets")               analyzer->process_wjets()              ;
  else if(runOn == "zjets_from50")        analyzer->process_zjets_from50()       ;
  else if(runOn == "zjets_10to50")        analyzer->process_zjets_10to50()       ;
  // FIXME: manage with --blah for having parameters like "withHiggs" and so on
  else if(runOn == "doTables"){
    cout << "Doing tables" << endl;
    bool onlyhiggs(true), sm(false), doNotPrintAllErrors(false), printAllErrors(true), includeSoverB(true), doNotincludeSoverB(false), withShapes(true), withStatShapes(false);
    int detailed(2), notDetailed(1);      
    analyzer->mcTable(notDetailed, includeSoverB, printAllErrors, sm, "PFlow", "yields-mc-", false, false, false); 
    cout << "Done SM table" << endl;
    analyzer->mcTable(notDetailed, includeSoverB, printAllErrors, onlyhiggs, "PFlow", "yields-mc-", false, false, false); 
    cout << "Done TBH table" << endl;
    analyzer->summaryTable( notDetailed, true, false, false, false, withShapes, withStatShapes);
    cout << "Done summary table" << endl;
  }
  else if(runOn == "doPlots"){
//    PlotStyle sty();
//    sty.setTDRStyle();
      
    TString samples("data/plotter/samples.xml");
    TString outFolder("plots/"); // move to input line
    TString cmd = "mkdir -p "+outFolder+"; cp data/plotter/index.html+"+outFolder+";";
    gSystem->Exec(cmd);
    TString limits("data/plotter/limits.xml");
    TString leptons("data/plotter/leptons.xml");
    TString met("data/plotter/met.xml");
    TString mt("data/plotter/mt.xml");
    TString wjets("data/plotter/wjets.xml");
    TString wtmass("data/plotter/wtmass.xml");
    TString jets("data/plotter/jets.xml");
    TString yields("data/plotter/yields.xml");
    TString taus("data/plotter/taus.xml");
    TString vertex("data/plotter/vertex.xml");
    TString afterR("data/plotter/afterR.xml");
    TString test("data/plotter/test.xml");
    TString debug("data/plotter/debug.xml");

    HistogramPlotter a; // Move to input line or cfg file the choice of what to plot
    a.parse(samples,vertex,outFolder);
    a.parse(samples,met,outFolder);      
    //    a.parse(samples,met,outFolder);      
    a.parse(samples,leptons,outFolder);  
    //a.parse(samples,limits,outFolder);      
    a.parse(samples,mt,outFolder);      
    a.parse(samples,yields,outFolder);
    a.parse(samples,jets,outFolder); 
    //a.parse(samples,debug,outFolder);  
    //a.parse(samples,afterR,outFolder);
    //
    //  //a.parse(samples,wjets,outFolder);
    //  //a.parse(samples,wtmass,outFolder);
    //  //a.parse(higgssamples,taus,outFolder);
    //  //a.parse(samples,taus,outFolder);
    //  //a.parse(samples,test,outFolder);
    
  }
  else
    cout << "Sample does not exist" << endl;
  


  
  cout << "physicsAnalysis reached its natural end" << endl;
  
  return 0;
  
}
//  // creates the physics analyzer
//  PhysicsAnalyzer analyzer(20);
//
//  
//  //analyzer.process_ttbar();  
//  //analyzer.process_data();
//  //analyzer.process_ttbar_bychannel();
//
//  // Table processing starts here
//  // include all errors -> (xx,xx,false,false)
//  bool onlyhiggs(true), sm(false), doNotPrintAllErrors(false), printAllErrors(true), includeSoverB(true), doNotincludeSoverB(false); 
//  int detailed(2), notDetailed(1);      
//analyzer.mcTable(notDetailed, includeSoverB, printAllErrors, sm, "PFlow", "yields-mc-", false, false, false); 
//analyzer.mcTable(notDetailed, includeSoverB, printAllErrors, onlyhiggs, "PFlow", "yields-mc-", false, false, false); 
////analyzer.mcTable(notDetailed, includeSoverB, doNotPrintAllErrors, onlyhiggs, "PFlow", "yields-mc-", false, false, false); 
//  //DEBUG (include all systematics )
//  //analyzer.mcTable(notDetailed, includeSoverB, doNotPrintAllErrors, sm, "PFlow", "yields-mc-", false, false, false);  
//  
//  //
// analyzer.summaryTable( notDetailed, false, false, false, false);
// analyzer.summaryTable( notDetailed, true, false, false, false);
//
// analyzer.summaryTable( notDetailed, false, true, false, false);
// analyzer.summaryTable( notDetailed, false, false, true, false);
// analyzer.summaryTable( notDetailed, false, true, true, false);
// analyzer.summaryTable( notDetailed, false, true, true, true);
//
// analyzer.summaryTable( notDetailed, true, true, false, false);
// analyzer.summaryTable( notDetailed, true, false, true, false);
// analyzer.summaryTable( notDetailed, true, true, true, false);
// analyzer.summaryTable( notDetailed, true, true, true, true);
//
//// analyzer.summaryTable( notDetailed, true, false, false, false);
//
//*/




