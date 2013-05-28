#include "LIP/TopTaus/interface/LandSShapesProducer.hh"


// System includes
#include <sstream>
#include <algorithm>

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT includes
#include <THStack.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TText.h>
#include <TPaveText.h>

using namespace std;
using namespace RooFit;

LandSShapesProducer::LandSShapesProducer(string parSet):
  parSet_(parSet)
{
  
  using namespace std; 
  using namespace RooFit ;
  
  Init();

  SetOptions();

}


void LandSShapesProducer::Init(){

  displayMin_ = 0.001;
  displayMax_ = 1.001;

//  displayMin_ = 0.;
//  displayMax_ = 200.;
  
  // Clear vectors
  nVars_ = 0;
  //  uncSources_.clear();
  systComponents_.clear();
  systFancyComponents_.clear();
  nSamples_ = 0;
  nMassPoints_ = 0;
  currentMassPoint_ = 999;
  fitVars_.clear();
  vars_.clear();
  mins_.clear(); 
  maxs_.clear();
  bins_.clear();
  hmin_.clear();
  hmax_.clear();
  unbinned_.clear(); 
  smoothOrder_.clear();
  
  minitree_.clear();
  baseDir_.clear();

  signalShapesToCompare_.clear();
  signalShapesToCompareHH_.clear();
  signalShapesToCompareWH_.clear();
  
  // Init RooFit variables
  // unused here. It's from TauDileptonPDFBuilderFitter


  // Reset identifier string
  identifier_ = "";

  cout << "Getting Parameter sets" << endl;

  // Get ParameterSet from cfg file
  const edm::ParameterSet &mFitPars = edm::readPSetsFrom(parSet_)->getParameter<edm::ParameterSet>("LandSShapesProducerParSet");

  outFolder_        = mFitPars.getParameter<std::string>("outFolder");
  outputFileName_  = mFitPars.getParameter<std::string>("outputFileName");
  massPointName_  = mFitPars.getParameter<vector<std::string> >("massPointName");
  nMassPoints_ = massPointName_.size();

  vector<string>baseDirTemp = mFitPars.getParameter<vector<std::string> >("baseDir");
  for(size_t f=0; f<baseDirTemp.size(); f++)
    baseDir_.push_back( baseDirTemp[f] );

 
  vector<string>inputFileNameTemp = mFitPars.getParameter<vector<std::string> >("inputFileName");
  for(size_t f=0; f<inputFileNameTemp.size(); f++)
    inputFileName_.push_back( inputFileNameTemp[f] );
  
  sampleName_      = mFitPars.getParameter<vector<std::string> >("sampleName");
  
  fancySampleName_      = mFitPars.getParameter<vector<std::string> >("fancySampleName");
  
  sampleColour_     = mFitPars.getParameter<vector<Int_t> >("sampleColour");
  sampleFillStyle_ = mFitPars.getParameter<vector<Int_t> >("sampleFillStyle");
  
  produceOnly_ = mFitPars.getParameter<bool>("produceOnly");

  useOS_            = mFitPars.getParameter<bool>("useOS");  
  standaloneTTbar_  = mFitPars.getParameter<bool>("standaloneTTbar");  
  catNBtags_        = mFitPars.getParameter<int>("NBtags");
  
  doMultiDimensionalShapes_       = mFitPars.getParameter<bool>("doMultiDimensionalShapes"); 

  unsplitUncertainties_       = mFitPars.getParameter<bool>("unsplitUncertainties"); // The finnish guys for some reason stack all the variations simultaneously.

  vector<string>minitreeTemp = mFitPars.getParameter<vector<std::string> >("minitree");
  for(size_t f=0; f<minitreeTemp.size(); f++)
    minitree_.push_back( minitreeTemp[f] );
  
  isFromData_       = mFitPars.getParameter<vector<Int_t> >("isFromData"); // Protect data and DD-bkg from syst
  isDDbkg_       = mFitPars.getParameter<vector<Int_t> >("isDDbkg"); // DD
  isSignal_       = mFitPars.getParameter<vector<Int_t> >("isSignal"); // DD


  
  vector<string>systComponentsTemp = mFitPars.getParameter<vector<std::string> >("systComponents");
  for(size_t f=0; f<systComponentsTemp.size(); f++)
    systComponents_.push_back( systComponentsTemp[f] );
  
  vector<string>systFancyComponentsTemp = mFitPars.getParameter<vector<std::string> >("systFancyComponents");
  for(size_t f=0; f<systFancyComponentsTemp.size(); f++)
    systFancyComponents_.push_back( systFancyComponentsTemp[f] );
  
  
  osCutEff_    = mFitPars.getParameter<double>("osCutEff");
  
  vars_        = mFitPars.getParameter<vector<string> >("vars");
  mins_        = mFitPars.getParameter<vector<double> >("mins");
  maxs_        = mFitPars.getParameter<vector<double> >("maxs");
  bins_        = mFitPars.getParameter<vector<int> >("bins");
  hmin_        = mFitPars.getParameter<vector<double> >("hmin");
  hmax_        = mFitPars.getParameter<vector<double> >("hmax");
  unbinned_    = mFitPars.getParameter<vector<Int_t> >("unbinned");
  smoothOrder_ = mFitPars.getParameter<vector<Int_t> >("smoothOrder");
  
  cout << "Got all parameter sets" << endl;
  
  // Open files and get trees
  // ddBkg is the only to be taken from data driven estimation (tree)
  nSamples_ = inputFileName_.size();
  nSysts_     = systFancyComponents_.size();
  
  cout << "Files opened" << endl;
  // Set variables
  //
  nVars_ = vars_.size();
  
  
  for(size_t i=0; i<nVars_; i++)
    fitVars_.push_back( new FitVar(vars_[i], mins_[i], maxs_[i], bins_[i], hmin_[i], hmax_[i], unbinned_[i], smoothOrder_[i]));
  
  
  // Set canvas
  SetTDRStyle();
  canvas_ = new TCanvas("canvas","My plots ",0,0,1000,1000);
  canvas_->cd();
  // set logy
  
  resultsFile_.open ((outFolder_+resultsFileName_).c_str());
  //  Uncomment following line in order to redirect stdout to file
  //  streamToFile_ = std::cout.rdbuf(resultsFile_.rdbuf());
  
  cout << "Init process complete" << endl;
}

void LandSShapesProducer::SetOptions(){
  //  myStyle_->SetOptStat(0);
}

void LandSShapesProducer::InitMassPoint(size_t s){

  cout << "Initializing mass point " << endl;

  th2_         .clear();
  th2StatUp_   .clear();
  th2StatDown_ .clear();
  th2Syst_     .clear();

  unrolled_         .clear();
  unrolledStatUp_   .clear();
  unrolledStatDown_ .clear();
  unrolledSyst_     .clear();
  

  currentMassPoint_ = s;
  inputFile_.clear();
  inputTree_.clear();
  
  systTree_.clear();

  masterHist_.clear();
  masterHistNames_.clear();
  masterShapes_.clear();

  perMassPointSignalShapesToCompareHH_.clear();
  perMassPointSignalShapesToCompareWH_.clear();
  perMassPointSignalShapesToCompare_  .clear();

  // Get files
  for(size_t f=0; f<nSamples_; f++){
    TString myTempFileName = inputFileName_[f];
    //if( f==1 || f == 2) myTempFileName = inputFileName_[f] + massPointName_[s] + TString(".root"); // Signals name manipulation
    cout << "Opening file " << baseDir_[f]   + myTempFileName << " opened." << endl;
    inputFile_.push_back( TFile::Open(baseDir_[f]   + myTempFileName  ) );
    cout << "File " << baseDir_[f]   + myTempFileName << " opened." << endl;
  }
  cout << "Got files" << endl;

  // Get base trees
  for(size_t f=0; f<nSamples_; f++)
    inputTree_.push_back( (TTree*) inputFile_[f]     ->Get(minitree_[f]) );
  
  cout << "Got base trees " << endl;
  // Get syst trees
  for(size_t a=0; a<nSysts_; a++){
    cout << "Getting tree " << systComponents_[a] << endl;
    vector<TTree*> temp;
    temp.clear();
    for(size_t f=0; f<nSamples_; f++){
      if(isFromData_[f]) temp.push_back( (TTree*) inputFile_[f] ->Get(minitree_[f]) ); // dummy fill with base minitree. Will not be used, just needed for vectors size
      else temp.push_back( (TTree*) inputFile_[f]     ->Get(minitree_[f] + systComponents_[a]) );
    }
    systTree_.push_back(temp);
  }
  
  cout << "Mass points set" << endl;
}
 
void LandSShapesProducer::InitPerVariableAmbient(size_t i){
  
  // Base trees
  for(size_t f=0; f<nSamples_; f++) inputTree_[f] ->ResetBranchAddresses()  ;

  // Syst trees
  for(size_t a=0; a<nSysts_; a++){
    for(size_t f=0; f<nSamples_; f++){
      cout << "variable " << i << ", syst " << a << ", sample " << f << ", entries " << systTree_[a][f]->GetEntries() << endl;
      systTree_[a][f]->ResetBranchAddresses();
    }
  }
  
  histo_.clear();

  myDSName_.clear();
  myDS_.clear();

  for(size_t a=0; a<nSysts_; a++){
    mySystDSName_.clear();
  }


  myDS_.clear();

  mySystDS_    .clear();

  hist_.clear();
  histStatUp_.clear();
  histStatDown_.clear();
  systHist_.clear();
  
  histo_.clear();
  systHisto_.clear();

  leg_            = 0;
  

  cout << "------------------------------------------------------" << endl;
  cout << "----------------- index " << i << "..................." << endl;
  cout << "Variable: " << fitVars_[i]->getVarName().c_str() << endl;
  cout << "Minimum: " << fitVars_[i]->getMin() << endl;
  cout << "Maximum: " << fitVars_[i]->getMax() << endl;
  cout << "Bins: " << fitVars_[i]->getBins() << endl;
  cout << "SmoothOrder: " << fitVars_[i]->getSmoothOrder() << endl;
  cout << "------------------------------------------------------" << endl;

  // Binned fit variable
  myvar_           = new RooRealVar(fitVars_[i]->getVarName().c_str(), fitVars_[i]->getVarName().c_str(), fitVars_[i]->getMin(), fitVars_[i]->getMax());  myvar_->setBins(fitVars_[i]->getBins()); 
  myvar_weights_   = new RooRealVar("weight","weight",0,1000);
  isOSvar_         = new RooRealVar("is_os","is_os",0,2);

  // In multidimensional case, add the next variable too
  if(doMultiDimensionalShapes_){
    size_t nextVar;
    if(i <nVars_-1)
      nextVar = i+1;
    else if(i== nVars_-1)
      nextVar = i-1;
    myvarMultiD_ = new RooRealVar(fitVars_[nextVar]->getVarName().c_str(), fitVars_[nextVar]->getVarName().c_str(), fitVars_[nextVar]->getMin(), fitVars_[nextVar]->getMax());
    myvarMultiD_->setBins(fitVars_[nextVar]->getBins());
  }


  //Define data sets /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(size_t f=0; f<nSamples_; f++){
    //fn    myDSName_.push_back( fitVars_[i]->getVarName().c_str() + sampleName_[f] + massPointName_[currentMassPoint_] );
    myDSName_.push_back( fitVars_[i]->getVarName().c_str() + sampleName_[f]);
  }

  // Define syst data sets
  for(size_t a=0; a<nSysts_; a++){
    vector<string> temp;
    temp.clear();
    for(size_t f=0; f<nSamples_; f++){
      temp.push_back( myDSName_[f] + string(systComponents_[a]) );
     }
     mySystDSName_.push_back(temp);
  }
  
  sumWeights_ = 0;
  cout << "per-variable ambient initialized" << endl;
}

void LandSShapesProducer::BuildDatasets(size_t i){

  
  // Temp variables for setting branch addresses
  double myVarAllocator, myVarWeightAllocator;
  double isOSsig, btagmult;


  // In multidimensional case, add the next variable too
  if(doMultiDimensionalShapes_){
    size_t nextVar;
    cout << "i is " << i << ", and nVars_-1 is " << nVars_-1 << endl;
    if(i <nVars_-1){
      cout << "I am inside" << endl;
      nextVar = i+1;
      //      else if(i== nVars_-1)
      //	nextVar = i-1;
      
      
      double myNextVarAllocator;
      
      for(size_t f=0; f<nSamples_; f++){
	
	th2_.push_back(new TH2D(massPointName_[currentMassPoint_] + sampleName_[f] + fitVars_[i]->getVarName().c_str() + TString("th2") , fitVars_[i]->getVarName().c_str() + TString("th2"), fitVars_[i]->getBins(), fitVars_[i]->getMin(), fitVars_[i]->getMax(), fitVars_[nextVar]->getBins(), fitVars_[nextVar]->getMin(), fitVars_[nextVar]->getMax())    );
	th2StatUp_.push_back(new TH2D(massPointName_[currentMassPoint_] + sampleName_[f] + fitVars_[i]->getVarName().c_str() + TString("th2StatUp") , fitVars_[i]->getVarName().c_str() + TString("th2StatUp"), fitVars_[i]->getBins(), fitVars_[i]->getMin(), fitVars_[i]->getMax(), fitVars_[nextVar]->getBins(), fitVars_[nextVar]->getMin(), fitVars_[nextVar]->getMax())    );
	th2StatDown_.push_back(new TH2D(massPointName_[currentMassPoint_] + sampleName_[f] + fitVars_[i]->getVarName().c_str() + TString("th2StatDown") , fitVars_[i]->getVarName().c_str() + TString("th2StatDown"), fitVars_[i]->getBins(), fitVars_[i]->getMin(), fitVars_[i]->getMax(), fitVars_[nextVar]->getBins(), fitVars_[nextVar]->getMin(), fitVars_[nextVar]->getMax())    );
	th2_[f]->Sumw2();
	th2StatUp_[f]->Sumw2();
	th2StatDown_[f]->Sumw2();
	
	myDS_.push_back( new RooDataSet(myDSName_[f].c_str(),myDSName_[f].c_str(), RooArgSet(*myvar_,*myvarMultiD_,*myvar_weights_), "weight") ); // This constructor does not accept the cut parameter
	// Get events
	inputTree_[f]->SetBranchAddress(fitVars_[i]->getVarName().c_str(), &myVarAllocator);
	inputTree_[f]->SetBranchAddress(fitVars_[nextVar]->getVarName().c_str(), &myNextVarAllocator);
	inputTree_[f]->SetBranchAddress("weight", &myVarWeightAllocator);
	inputTree_[f]->SetBranchAddress("is_os", &isOSsig);
	inputTree_[f]->SetBranchAddress("btagmultiplicity_j", &btagmult);
	//    cout << "getIsoS      ";
	for(unsigned int ev=0; ev<inputTree_[f]->GetEntries(); ev++){
	  inputTree_[f]->GetEntry(ev);
	  if( myVarAllocator < fitVars_[i]->getMin() /* || myVarAllocator > fitVars_[i]->getMax() */
	      ) continue;
	  if (useOS_ && !isDDbkg_[f] ) {
	    if(isOSsig < 0.5 ) continue;
	  }
	  if(catNBtags_ == 1 && btagmult != 1) continue;
	  if(catNBtags_ > 1 && btagmult < 2) continue;

	  if(nextVar== i-1 ) continue;
	  myvar_->setVal(myVarAllocator);
	  myvarMultiD_->setVal(myNextVarAllocator);
	  sumWeights_ += myVarWeightAllocator;
	  if(isDDbkg_[f]){
	    myvar_weights_->setVal(osCutEff_*myVarWeightAllocator);
	    myDS_[f]->add(RooArgSet(*myvar_,*myvar_weights_),osCutEff_*myVarWeightAllocator);
	  }
	  else{
	    myvar_weights_->setVal(myVarWeightAllocator);
	    myDS_[f]->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
	  }
	  
	  th2_[f]->Fill(myVarAllocator, myNextVarAllocator, myVarWeightAllocator);
	  //	  cout << "Filled th2 with " << myVarAllocator << ", " << myNextVarAllocator << ", " << myVarWeightAllocator << endl;
	}
      }
      
      
      cout << "Got base datasets" << endl;
      
      ///// Syst get
      for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
	
	vector<TH2D*> tempTh2;
	tempTh2.clear();
	
	vector<RooDataSet*> temp;
	temp.clear();
	for(size_t f=0; f<nSamples_; f++){
	  tempTh2.push_back(new TH2D(massPointName_[currentMassPoint_] + sampleName_[f] + fitVars_[i]->getVarName().c_str() + systComponents_[a] + TString("th2") , fitVars_[i]->getVarName().c_str() + systComponents_[a] + TString("th2"), fitVars_[i]->getBins(), fitVars_[i]->getMin(), fitVars_[i]->getMax(), fitVars_[nextVar]->getBins(), fitVars_[nextVar]->getMin(), fitVars_[nextVar]->getMax())    );
	  tempTh2[f]->Sumw2();
	  temp.push_back( new RooDataSet(mySystDSName_[a][f].c_str(),mySystDSName_[a][f].c_str(), RooArgSet(*myvar_,*myvarMultiD_,*myvar_weights_), "weight") ); // This constructor does not accept the cut parameter
	}
	mySystDS_.push_back(temp);
	th2Syst_.push_back(tempTh2);
	for(size_t f=0; f<nSamples_; f++){
	  // Get MCBkg events
	  systTree_[a][f]->SetBranchAddress(fitVars_[i]->getVarName().c_str(), &myVarAllocator);
	  systTree_[a][f]->SetBranchAddress(fitVars_[nextVar]->getVarName().c_str(), &myNextVarAllocator);
	  systTree_[a][f]->SetBranchAddress("weight", &myVarWeightAllocator);
	  systTree_[a][f]->SetBranchAddress("is_os", &isOSsig);
	  systTree_[a][f]->SetBranchAddress("btagmultiplicity_j", &btagmult);
	  //    cout << "getIsoS      ";
	  for(unsigned int ev=0; ev<systTree_[a][f]->GetEntries(); ev++){
	    systTree_[a][f]->GetEntry(ev);
	    if( myVarAllocator < fitVars_[i]->getMin() /* || myVarAllocator > fitVars_[i]->getMax() */
		) continue;
	    if (useOS_ && !isDDbkg_[f]) {
	      if(isOSsig < 0.5 ) continue;
	    }
	    if(catNBtags_ == 1 && btagmult != 1) continue;
	    if(catNBtags_ > 1 && btagmult < 2) continue;	    

	    if(nextVar== i-1 ) continue;
	    myvar_->setVal(myVarAllocator);
	    myvarMultiD_->setVal(myNextVarAllocator);
	    sumWeights_ += myVarWeightAllocator;
	    if(isDDbkg_[f]){
	      myvar_weights_->setVal(osCutEff_*myVarWeightAllocator);
	      mySystDS_[a][f]->add(RooArgSet(*myvar_,*myvar_weights_),osCutEff_*myVarWeightAllocator);
	    }
	    else{
	      myvar_weights_->setVal(myVarWeightAllocator);
	      mySystDS_[a][f]->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
	    }
	    th2Syst_[a][f]->Fill(myVarAllocator, myNextVarAllocator, myVarWeightAllocator);
	    //	    cout << "Filled th2 with " << myVarAllocator << ", " << myNextVarAllocator << ", " << myVarWeightAllocator << endl;
	  }
	}
	
      } // End loop on syst components
      ///// End syst get
    }
    cout << "Done multiD stuff" << endl;
  }
  else {
  
    for(size_t f=0; f<nSamples_; f++){
      myDS_.push_back( new RooDataSet(myDSName_[f].c_str(),myDSName_[f].c_str(), RooArgSet(*myvar_,*myvar_weights_), "weight") ); // This constructor does not accept the cut parameter
      // Get events
      inputTree_[f]->SetBranchAddress(fitVars_[i]->getVarName().c_str(), &myVarAllocator);
      inputTree_[f]->SetBranchAddress("weight", &myVarWeightAllocator);
      inputTree_[f]->SetBranchAddress("is_os", &isOSsig);
      inputTree_[f]->SetBranchAddress("btagmultiplicity_j", &btagmult);
      //    cout << "getIsoS      ";
      for(unsigned int ev=0; ev<inputTree_[f]->GetEntries(); ev++){
	inputTree_[f]->GetEntry(ev);
	if( myVarAllocator < fitVars_[i]->getMin() /* || myVarAllocator > fitVars_[i]->getMax() */
	    ) continue;
	if (useOS_ && !isDDbkg_[f]) {
	  if(isOSsig < 0.5) continue;
	}
	if(catNBtags_ == 1 && btagmult != 1) continue;
	if(catNBtags_ > 1 && btagmult < 2) continue;
	
	myvar_->setVal(myVarAllocator);
	sumWeights_ += myVarWeightAllocator;
	if(isDDbkg_[f]){
	  myvar_weights_->setVal(osCutEff_*myVarWeightAllocator);
	  myDS_[f]->add(RooArgSet(*myvar_,*myvar_weights_),osCutEff_*myVarWeightAllocator);
	}
	else{
	  myvar_weights_->setVal(myVarWeightAllocator);
	  myDS_[f]->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
	}
      }
    }

    cout << "Got base datasets" << endl;
    
    ///// Syst get
    for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
      vector<RooDataSet*> temp;
      temp.clear();
      for(size_t f=0; f<nSamples_; f++){
	temp.push_back( new RooDataSet(mySystDSName_[a][f].c_str(),mySystDSName_[a][f].c_str(), RooArgSet(*myvar_,*myvar_weights_), "weight") ); // This constructor does not accept the cut parameter
      }
      mySystDS_.push_back(temp);
      for(size_t f=0; f<nSamples_; f++){
	// Get MCBkg events
	systTree_[a][f]->SetBranchAddress(fitVars_[i]->getVarName().c_str(), &myVarAllocator);
	systTree_[a][f]->SetBranchAddress("weight", &myVarWeightAllocator);
	systTree_[a][f]->SetBranchAddress("is_os", &isOSsig);
	systTree_[a][f]->SetBranchAddress("btagmultiplicity_j", &btagmult);
	//    cout << "getIsoS      ";
	for(unsigned int ev=0; ev<systTree_[a][f]->GetEntries(); ev++){
	  systTree_[a][f]->GetEntry(ev);
	  if( myVarAllocator < fitVars_[i]->getMin() /* || myVarAllocator > fitVars_[i]->getMax() */
	      ) continue;
	  if (useOS_ && !isDDbkg_[f]) {
	    if(isOSsig < 0.5) continue;
	  }
	  if(catNBtags_ == 1 && btagmult != 1) continue;
	  if(catNBtags_ > 1 && btagmult < 2) continue;	  

	  myvar_->setVal(myVarAllocator);
	  sumWeights_ += myVarWeightAllocator;
	  if(isDDbkg_[f]){
	    myvar_weights_->setVal(osCutEff_*myVarWeightAllocator);
	    mySystDS_[a][f]->add(RooArgSet(*myvar_,*myvar_weights_),osCutEff_*myVarWeightAllocator);
	  }
	  else{
	    myvar_weights_->setVal(myVarWeightAllocator);
	    mySystDS_[a][f]->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
	  }
	}
      }
      
    } // End loop on syst components
    ///// End syst get
    

  cout << "Got syst datasets" << endl;
  
  
  cout << "All datasets built" << endl;
  
  // // Legacy: seems not to give correct inputs
  //  string myOsCut = "is_os>0.5";
  //  unrMyDDBkgDS_      = new RooDataSet(myDDBkgDSName_.c_str(), myDDBkgDSName_.c_str(),  ddBkgTree_,  RooArgSet(*myvar_,*myvar_weights_,*isOSvar_),myOsCut.c_str(),"weight" );
  //  myDDBkgDS_ = (RooDataSet*) unrMyDDBkgDS_->reduce(RooArgSet(*myvar_,*myvar_weights_));
  //  if(standaloneTTbar_){
  //    unrMyTTBARMCBkgDS_ = new RooDataSet(myTTBARMCBkgDSName_.c_str(), myTTBARMCBkgDSName_.c_str(),  ttbarmcBkgTree_,  RooArgSet(*myvar_,*myvar_weights_,*isOSvar_),myOsCut.c_str(),"weight" );
  //    myTTBARMCBkgDS_ = (RooDataSet*) unrMyTTBARMCBkgDS_->reduce(RooArgSet(*myvar_,*myvar_weights_));
  //  }
  //  unrMyMCBkgDS_      = new RooDataSet(myMCBkgDSName_.c_str(), myMCBkgDSName_.c_str(),  mcBkgTree_,  RooArgSet(*myvar_,*myvar_weights_,*isOSvar_),myOsCut.c_str(),"weight" );
  //  myMCBkgDS_ = (RooDataSet*) unrMyMCBkgDS_->reduce(RooArgSet(*myvar_,*myvar_weights_));
  //  unrMyDataDS_       = new RooDataSet(myDataDSName_.c_str(),  myDataDSName_.c_str(),   dataTree_,   RooArgSet(*myvar_,*isOSvar_), myOsCut.c_str() );
  //  myDataDS_ = (RooDataSet*) unrMyDataDS_->reduce(RooArgSet(*myvar_,*myvar_weights_));
  //  
  
  // Build binned clones
  for(size_t f=0; f<nSamples_; f++) histo_.push_back( myDS_[f]->binnedClone());

  cout << "Binned clones built" << endl;

  for(size_t f=0; f<nSamples_; f++)
    cout << myDSName_[f] << " unbinned entries: " << myDS_[f]->numEntries() << ". weighted entries: " << histo_[f]->sum(kFALSE) << endl;
  
  cout << "BuildDatasets successful" << endl;
  


  /// Build binned clones - syst case
  for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
    vector<RooDataHist*> temp;
    temp.clear();
    for(size_t f=0; f<nSamples_; f++){
      temp.push_back( mySystDS_[a][f]->binnedClone() );
    }
    systHisto_.push_back(temp);
    cout << "Syst binned clones built" << endl;
    
    for(size_t f=0; f<nSamples_; f++)
      cout << mySystDSName_[a][f] << " unbinned entries: " << mySystDS_[a][f]->numEntries() << ". weighted entries: " << systHisto_[a][f]->sum(kFALSE) << endl;
    
  } // End syst loop
    /// End build binned clones - syst case
  }  
  cout << "BuildDatasets successful" << endl;
}

//void LandSShapesProducer::BuildPDFs(size_t i){
//  
//  switch(fitVars_[i]->getUnbinned()){
//  case 1 : // Unbinned
//    if(includeSignal_) u_signalModel_ = new RooKeysPdf(signalModelName_.c_str(), signalModelName_.c_str(), *myvar_, *mySignalDS_, RooKeysPdf::MirrorBoth,2);
//    u_ddbkgModel_  = new RooKeysPdf(ddbkgModelName_.c_str(),  ddbkgModelName_.c_str(),  *myvar_, *myDDBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    if(standaloneTTbar_) u_ttbarmcbkgModel_  = new RooKeysPdf(ttbarmcbkgModelName_.c_str(),  ttbarmcbkgModelName_.c_str(),  *myvar_, *myTTBARMCBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    u_mcbkgModel_  = new RooKeysPdf(mcbkgModelName_.c_str(),  mcbkgModelName_.c_str(),  *myvar_, *myMCBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    break;
//  case 0:  // Binned (w/ or w/out smoothing)
//    if(includeSignal_) b_signalModel_ = new RooHistPdf(signalModelName_.c_str(), signalModelName_.c_str(), *myvar_, *signalHisto_, fitVars_[i]->getSmoothOrder());
//    b_ddbkgModel_  = new RooHistPdf(ddbkgModelName_.c_str(),  ddbkgModelName_.c_str(), *myvar_, *ddbkgHisto_ , fitVars_[i]->getSmoothOrder());
//    if(standaloneTTbar_) b_ttbarmcbkgModel_ = new RooHistPdf(ttbarmcbkgModelName_.c_str(),  ttbarmcbkgModelName_.c_str(),  *myvar_, *ttbarmcbkgHisto_  , fitVars_[i]->getSmoothOrder()); 
//    b_mcbkgModel_  = new RooHistPdf(mcbkgModelName_.c_str(),  mcbkgModelName_.c_str(),  *myvar_, *mcbkgHisto_  , fitVars_[i]->getSmoothOrder()); 
//    break;
//  default : // Dummy - should never arrive here
//    cout<<"Neither binned not unbinned. Check your options motherfucker."<<endl;
//  }
//}


void LandSShapesProducer::DrawTemplates(size_t i){
  
  // Open file for storing shapes
  string outputFileName = outputFileName_ + string("_") + massPointName_[currentMassPoint_] + string("_") + fitVars_[i]->getVarName();
  TFile* outputFile;
  if(produceOnly_ && !doMultiDimensionalShapes_)
    outputFile = new TFile((outFolder_+outputFileName+string(".root")).c_str(), "RECREATE");
  
  // Reset canvas
  canvas_->cd();
  canvas_->Clear();
  


  //  /////////////////////////////////////////////////////////////////////
  
  // mc bkg histogram ////////////////////////////////////////////////
  for(size_t f=0; f<nSamples_; f++){
    hist_.push_back( histo_[f]->createHistogram(fitVars_[i]->getVarName().c_str(),fitVars_[i]->getBins() ) );
    hist_[f]->Sumw2();
    if(isSignal_[f]){
      hist_[f]->SetOption("0000"); // What did that do, again? LoL
      hist_[f]->SetLineColor(sampleColour_[f]);
    }
    else
      hist_[f]->SetFillColor(sampleColour_[f]);
    hist_[f]->SetLineWidth(3);
    hist_[f]->SetFillStyle(sampleFillStyle_[f]); //1001 for background and data, 0 for signal // 3017);
    
    
    cout << "sono nel punto giusto" << endl;
    cout << "sample " << sampleName_[f] << ", integral " << hist_[f]->Integral() << endl;
    if(isDDbkg_[f]){
      hist_[f]->Scale(osCutEff_*2933.21/hist_[f]->Integral());
      //hist_[f]->Scale(1.0458/hist_[f]->Integral());
      ddbkgHistUp_ =   (TH1*) hist_[f]->Clone(hist_[f]->GetName() + TString("Up") );
      ddbkgHistDown_ = (TH1*) hist_[f]->Clone(hist_[f]->GetName() + TString("Down") );
      
      ddbkgHistUp_->SetFillColor(sampleColour_[f]);
      ddbkgHistUp_->SetLineWidth(3);
      ddbkgHistUp_->SetFillStyle(sampleFillStyle_[f]);//3017);
      
      ddbkgHistDown_->SetFillColor(sampleColour_[f]);
      ddbkgHistDown_->SetLineWidth(3);
      ddbkgHistDown_->SetFillStyle(sampleFillStyle_[f]);//3017);
      
      ddbkgHistUp_->Scale(osCutEff_*1.088*2933.21/ddbkgHistUp_->Integral());
      ddbkgHistDown_->Scale(osCutEff_*0.912*2933.21/ddbkgHistDown_->Integral());
      cout << "sample " << sampleName_[f] << ", integral " << hist_[f]->Integral() << " after rescaling " << endl;
    }
  }
  //    ///////////////////////////////////////////////////////////////////
  
  
  // hist_.push_back( (TH1*)hist_[nSamples_-1]->Clone( sampleName_[nSamples_].c_str())   ); // nSamples-1 Must always be di-bosons
  // hist_.push_back( (TH1*)hist_[nSamples_-1]->Clone( sampleName_[nSamples_+1].c_str()) ); // nSamples-1 Must always be di-bosons
  // hist_[nSamples_]  ->Sumw2();
  // hist_[nSamples_+1]->Sumw2();

  // hist_[nSamples_]->Scale(50.5/hist_[nSamples_]->Integral());      // Normalize to sample
  // hist_[nSamples_+1]->Scale(0.4/hist_[nSamples_+1]->Integral()); // Normalize to sample

  // cout << "sample " << sampleName_[nSamples_] << ", integral " << hist_[nSamples_]->Integral() << endl;
  // cout << "sample " << sampleName_[nSamples_+1] << ", integral " << hist_[nSamples_+1]->Integral() << endl;

  ddbkgHistUp_ =   (TH1*) hist_[1]->Clone(hist_[1]->GetName() + TString("Up") );
  ddbkgHistDown_ = (TH1*) hist_[1]->Clone(hist_[1]->GetName() + TString("Down") );

  if ( nSamples_ > 9) { //fn To be fixed seriously...

    for(size_t f=0; f<nSamples_+2; f++){
      if(f>nSamples_-1){ // Repeat colours for newly cloned histos
	if(isSignal_[f]){
	  hist_[f]->SetOption("0000"); // What did that do, again? LoL
	  hist_[f]->SetLineColor(sampleColour_[f]);
	}
	else
	  hist_[f]->SetFillColor(sampleColour_[f]);
	hist_[f]->SetLineWidth(3);
	hist_[f]->SetFillStyle(sampleFillStyle_[f]); //1001 for background and data, 0 for signal // 3017);
      }
      
      //    hist_[f]->Sumw2();
      histStatUp_  .push_back( (TH1*)hist_[f]->Clone(hist_[f]->GetName() + TString("_StatUp")));
      histStatDown_.push_back( (TH1*)hist_[f]->Clone(hist_[f]->GetName() + TString("_StatDown")));
      histStatUp_[f]->Sumw2();
      histStatDown_[f]->Sumw2();
    }
    
    
    /// Produce plots - syst case
    for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
      
      // mc bkg histogram ////////////////////////////////////////////////
      vector<TH1*> temp;
      temp.clear();
      for(size_t f=0; f<nSamples_; f++){
	temp.push_back( systHisto_[a][f]->createHistogram(fitVars_[i]->getVarName().c_str(),fitVars_[i]->getBins() ) );
      }
      systHist_.push_back(temp);
      
      for(size_t f=0; f<nSamples_; f++){
	systHist_[a][f]->Sumw2();
	systHist_[a][f]->SetFillColor(sampleColour_[f]);
	systHist_[a][f]->SetFillStyle(sampleFillStyle_[f]);
	if(isDDbkg_[f]){    
	  if(a == 0 || a == 2 || a == 4) //fn I fill DataDriven syst shapes, but I don't use them...
	    systHist_[a][f]->Scale((222+11.4/5)/systHist_[a][f]->Integral());
	  else
	    systHist_[a][f]->Scale((222-11.4/5)/systHist_[a][f]->Integral());
	}
      }
    //    ///////////////////////////////////////////////////////////////////
    
      systHist_[a].push_back( (TH1*)systHist_[a][nSamples_-1]->Clone(  sampleName_[nSamples_].c_str()) );
      systHist_[a].push_back( (TH1*)systHist_[a][nSamples_-1]->Clone(  sampleName_[nSamples_+1].c_str()) );
      systHist_[a][nSamples_]  ->Sumw2();
      systHist_[a][nSamples_+1]->Sumw2();
      
      systHist_[a][nSamples_]  ->Scale( (50.5/systHist_[a][nSamples_]->Integral()  )*( systHist_[a][nSamples_-1]->Integral() / hist_[nSamples_-1]->Integral()  )    );
      systHist_[a][nSamples_+1]->Scale( (0.4/systHist_[a][nSamples_+1]->Integral() )*( systHist_[a][nSamples_-1]->Integral() / hist_[nSamples_-1]->Integral()  )  );
      
      for(size_t f=0; f<nSamples_+2; f++){
	if(f>nSamples_-1){ // Repeat colours for newly cloned histos
	  systHist_[a][f]->Sumw2();
	  if(isSignal_[f]){
	    systHist_[a][f]->SetOption("0000"); // What did that do, again? LoL
	    systHist_[a][f]->SetLineColor(sampleColour_[f]);
	  }
	  else
	    systHist_[a][f]->SetFillColor(sampleColour_[f]);
	  systHist_[a][f]->SetLineWidth(3);
	  systHist_[a][f]->SetFillStyle(sampleFillStyle_[f]); //1001 for background and data, 0 for signal // 3017);
	}
	cout << "sample " << sampleName_[f] << "" << systComponents_[a] << ", integral " << systHist_[a][f]->Integral() << endl;
      }
    } // End loop on systs
    // End produce plots - syst case
   
  // Stat filling and histogram putting nonzeroes instead of zeroes (RooFit mojo)
  // Perhaps do a single loop with systs.
  
    for(size_t f=0; f<nSamples_+2; f++){
      for(int ibin=1; ibin<hist_[f]->GetNbinsX(); ibin++){ // <= is for overflow
	
	hist_[nSamples_]->Scale(50.5/hist_[nSamples_]->Integral());      // Normalize to sample
	hist_[nSamples_+1]->Scale(0.4/hist_[nSamples_+1]->Integral()); // Normalize to sample
	
	//      if(f == 1 | f == 2)
	if(ibin==1 && hist_[f]->GetBinContent(ibin) != 0){
	  hist_[f]->SetBinError(ibin,  hist_[f]->GetBinContent(ibin) * hist_[f]->GetBinError(ibin+1)/hist_[f]->GetBinContent(ibin+1)  );
	  cout << "bin " << ibin << ", bin error: " << hist_[f]->GetBinError(ibin) << endl;
	  //	if(f<nSamples_){
	  //	  cout << "sample " << hist_[f]->GetName();
	  //	  cout << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	  //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral()     ) );
	  //	  cout << "done sample " << hist_[f]->GetName() << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	  //	} else{
	  //	  cout << "sample " << hist_[f]->GetName();
	  //	  cout << ", bincontent scaled " <<    myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral() << endl;
	  //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral()     ) );
	  //	}
	}
	
	if(hist_[f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	  hist_[f]->SetBinContent(ibin, 1e-6);
	  cout << hist_[f]->GetName() << " has 0 bin with content " << hist_[f]->GetBinContent(ibin) << " and error " << hist_[f]->GetBinError(ibin) << endl;
	  //	hist_[f]->SetBinError(ibin,0);
	}
	// Fill stat uncertainty histograms
	histStatUp_[f]  ->SetBinContent(ibin, std::max( hist_[f]->GetBinContent(ibin) + hist_[f]->GetBinError(ibin) , 1e-6 ) );
	histStatDown_[f]->SetBinContent(ibin, std::max( hist_[f]->GetBinContent(ibin) - hist_[f]->GetBinError(ibin) , 1e-6 ) );
	cout << "==================================" << hist_[f]->GetName() << ": statUp " << hist_[f]->GetBinContent(ibin) + hist_[f]->GetBinError(ibin) << ", statDown " << hist_[f]->GetBinContent(ibin) - hist_[f]->GetBinError(ibin) << endl;
      }
    }

  //  hist_[nSamples_]->Scale(50.5/hist_[nSamples_]->Integral());      // Normalize to sample
  //  hist_[nSamples_+1]->Scale(0.4/hist_[nSamples_+1]->Integral()); // Normalize to sample
  // Rescale stat plots - they must have the same integral as the base ones
    for(size_t f=0; f<nSamples_+2; f++){
      histStatUp_[f]->Scale(hist_[f]->Integral()/histStatUp_[f]->Integral());
      histStatDown_[f]->Scale(hist_[f]->Integral()/histStatDown_[f]->Integral());
    }

  }
  else { //CLOSED by the "End of the finnish way" line... 

    for(size_t f=0; f<nSamples_; f++){
      hist_[f]->SetFillColor(sampleColour_[f]);
      hist_[f]->SetLineWidth(3);
      hist_[f]->SetFillStyle(sampleFillStyle_[f]); //1001 for background and data, 0 for signal // 3017);
      
      //    hist_[f]->Sumw2();
      histStatUp_  .push_back( (TH1*)hist_[f]->Clone(hist_[f]->GetName() + TString("_StatUp")));
      histStatDown_.push_back( (TH1*)hist_[f]->Clone(hist_[f]->GetName() + TString("_StatDown")));
      histStatUp_[f]->Sumw2();
      histStatDown_[f]->Sumw2();
    }

    /// Produce plots - syst case
    for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
      
      // mc bkg histogram ////////////////////////////////////////////////
      vector<TH1*> temp;
      temp.clear();
      for(size_t f=0; f<nSamples_; f++){
	temp.push_back( systHisto_[a][f]->createHistogram(fitVars_[i]->getVarName().c_str(),fitVars_[i]->getBins() ) );
      }
      systHist_.push_back(temp);
      
      for(size_t f=0; f<nSamples_; f++){
	systHist_[a][f]->Sumw2();
	systHist_[a][f]->SetFillColor(sampleColour_[f]);
	systHist_[a][f]->SetFillStyle(sampleFillStyle_[f]);
	systHist_[a][f]->SetLineWidth(3);

	cout << "sample " << sampleName_[f] << "" << systComponents_[a] << ", integral " << systHist_[a][f]->Integral() << endl; 
      }
    } // End loop on systs
    // End produce plots - syst case


    // Stat filling and histogram putting nonzeroes instead of zeroes (RooFit mojo)
    // Perhaps do a single loop with systs.
     for(size_t f=0; f<nSamples_; f++){
      for(int ibin=1; ibin<hist_[f]->GetNbinsX(); ibin++){ // <= is for overflow	
	
	//      if(f == 1 | f == 2)
	if(ibin==1 && hist_[f]->GetBinContent(ibin) != 0){
	  hist_[f]->SetBinError(ibin,  hist_[f]->GetBinContent(ibin) * hist_[f]->GetBinError(ibin+1)/hist_[f]->GetBinContent(ibin+1)  );
	  cout << "bin " << ibin << ", bin error: " << hist_[f]->GetBinError(ibin) << endl;
	}
	
	if(hist_[f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	  hist_[f]->SetBinContent(ibin, 1e-6);
	  cout << hist_[f]->GetName() << " has 0 bin with content " << hist_[f]->GetBinContent(ibin) << " and error " << hist_[f]->GetBinError(ibin) << endl;
	  //	hist_[f]->SetBinError(ibin,0);
	}
	// Fill stat uncertainty histograms
	histStatUp_[f]  ->SetBinContent(ibin, std::max( hist_[f]->GetBinContent(ibin) + hist_[f]->GetBinError(ibin) , 1e-6 ) );
	histStatDown_[f]->SetBinContent(ibin, std::max( hist_[f]->GetBinContent(ibin) - hist_[f]->GetBinError(ibin) , 1e-6 ) );
	cout << "==================================" << hist_[f]->GetName() << ": statUp " << hist_[f]->GetBinContent(ibin) + hist_[f]->GetBinError(ibin) << ", statDown " << hist_[f]->GetBinContent(ibin) - hist_[f]->GetBinError(ibin) << endl;
      }
     }

  //  hist_[nSamples_]->Scale(50.5/hist_[nSamples_]->Integral());      // Normalize to sample
  //  hist_[nSamples_+1]->Scale(0.4/hist_[nSamples_+1]->Integral()); // Normalize to sample
  // Rescale stat plots - they must have the same integral as the base ones
//fn     for(size_t f=0; f<nSamples_; f++){
//fn       histStatUp_[f]->Scale(hist_[f]->Integral()/histStatUp_[f]->Integral());
//fn       histStatDown_[f]->Scale(hist_[f]->Integral()/histStatDown_[f]->Integral());
//fn     }
    
    if(!unsplitUncertainties_){ // the correct thing to do
      for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
	for(size_t f=0; f<nSamples_; f++){
	  for(int ibin=1; ibin<systHist_[a][f]->GetNbinsX(); ibin++){ // <= is for overflow
	    
	    if(ibin==1 && systHist_[a][f]->GetBinContent(ibin) != 0){
	      systHist_[a][f]->SetBinError(ibin,  systHist_[a][f]->GetBinContent(ibin) * systHist_[a][f]->GetBinError(ibin+1)/systHist_[a][f]->GetBinContent(ibin+1)  );
	      //	if(f<nSamples_){
	      //	  cout << "sample " << hist_[f]->GetName();
	      //	  cout << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	      //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral()     ) );
	      //	  cout << "done sample " << hist_[f]->GetName() << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	      //	} else{
	      //	  cout << "sample " << hist_[f]->GetName();
	      //	  cout << ", bincontent scaled " <<    myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral() << endl;
	      //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral()     ) );
	      //	}
	    }
	    
	    if(systHist_[a][f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	      systHist_[a][f]->SetBinContent(ibin, 1e-6);
	      cout << systHist_[a][f]->GetName() << " has 0 bin error " << systHist_[a][f]->GetBinError(ibin) << endl;
	      //	hist_[f]->SetBinError(ibin,0);
	    }
	  }
	}
      } // end loop on syst components
    } // End of the correct way
    else{ // The finnish way: JES, JER, MET are stacked one over the other like if they were correlated
      
      // Stack them
      for(size_t f=0; f<nSamples_; f++){
	systHist_[0][f]->Add(systHist_[2][f],1);
	systHist_[0][f]->Add(systHist_[4][f],1);
	systHist_[0][f]->Add(hist_[f],-2);
	
	systHist_[1][f]->Add(systHist_[3][f],1);
	systHist_[1][f]->Add(systHist_[5][f],1);
	systHist_[1][f]->Add(hist_[f],-2);
	
	systHist_[0][f]->Scale(hist_[f]->Integral()/systHist_[0][f]->Integral());
	systHist_[1][f]->Scale(hist_[f]->Integral()/systHist_[1][f]->Integral());
	
      }
      
      
      for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
	for(size_t f=0; f<nSamples_; f++){
	  for(int ibin=1; ibin<systHist_[a][f]->GetNbinsX(); ibin++){ // <= is for overflow
	    
	    if(ibin==1 && systHist_[a][f]->GetBinContent(ibin) != 0){
	      systHist_[a][f]->SetBinError(ibin,  systHist_[a][f]->GetBinContent(ibin) * systHist_[a][f]->GetBinError(ibin+1)/systHist_[a][f]->GetBinContent(ibin+1)  );
	      //	if(f<nSamples_){
	      //	  cout << "sample " << hist_[f]->GetName();
	      //	  cout << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	      //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral()     ) );
	      //	  cout << "done sample " << hist_[f]->GetName() << ", bincontent scaled " <<    myDS_[f]->numEntries()*hist_[f]->GetBinContent(ibin)/hist_[f]->Integral() << endl;
	      //	} else{
	      //	  cout << "sample " << hist_[f]->GetName();
	      //	  cout << ", bincontent scaled " <<    myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral() << endl;
	      //	  hist_[f]->SetBinError(ibin, sqrt( myDS_[nSamples_-1]->numEntries()*hist_[nSamples_-1]->GetBinContent(ibin)/hist_[nSamples_-1]->Integral()     ) );
	      //	}
	    }
	    
	    
	    
	    
	    //	  if(systHist_[a][f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	    //	    systHist_[a][f]->SetBinContent(ibin, 1e-6);
	    //	    cout << systHist_[a][f]->GetName() << " has 0 bin error " << systHist_[a][f]->GetBinError(ibin) << endl;
	    //	    //	hist_[f]->SetBinError(ibin,0);
	    //	  }
	  }
	}
      } // end loop on syst components
      
      for(size_t f=0; f<nSamples_; f++){
	systHist_[0][f]->Add(systHist_[2][f], 1);
	systHist_[0][f]->Add(systHist_[4][f], 1);
	systHist_[0][f]->Add(hist_[f], -2);
	
	systHist_[1][f]->Add(systHist_[3][f], 1);
	systHist_[1][f]->Add(systHist_[5][f], 1);
	systHist_[1][f]->Add(hist_[f], -2);
	
	for(int ibin=1; ibin<systHist_[0][f]->GetNbinsX(); ibin++){ // <= is for overflow
	  if(systHist_[0][f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	    systHist_[0][f]->SetBinContent(ibin, 1e-6);
	  }
	}
	
	for(int ibin=1; ibin<systHist_[1][f]->GetNbinsX(); ibin++){ // <= is for overflow
	  if(systHist_[1][f]->GetBinContent(ibin) == 0){ // Protection from RooStat/Fit breaking in combine
	    systHist_[1][f]->SetBinContent(ibin, 1e-6);
	  }
	}
	
      }    
      
      
    } // End of the finnish way 
    

    // Set names
    if(!doMultiDimensionalShapes_){
      
      for(size_t f=0; f<nSamples_; f++){
	hist_[f]->SetName(sampleName_[f].c_str());
	histStatUp_[f]->SetName(sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatUp")); // Double name because of datacard syntax (indipendent lines for each stat)
	histStatDown_[f]->SetName(sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatDown")); // Double name because of datacard syntax (indipendent lines for each stat)
      
      //    hist_[f]->Sumw2();
      //    histStatUp_[f]->Sumw2();
      //    histStatDown_[f]->Sumw2();
    }
      
      cout << "nSysts_ =" << nSysts_ << endl; 
      // Syst case
      if(!unsplitUncertainties_){ // Correct way
	
	//fn	ddbkgHistUp_ =   (TH1*) histStatUp_[2]->Clone();
	//fn	ddbkgHistDown_ =   (TH1*) histStatDown_[2]->Clone();
	//fn	ddbkgHistUp_->SetName(sampleName_[2].c_str()+systFancyComponents_[10]);
	//fn	ddbkgHistDown_->SetName(sampleName_[2].c_str()+systFancyComponents_[11]);
	
	for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
	  for(size_t f=0; f<nSamples_; f++){
	    systHist_[a][f]->SetName(sampleName_[f].c_str()+systFancyComponents_[a]);
	  }
	} // End syst loop
      } else{ // Finnish way (all stacked)
	for(size_t f=0; f<nSamples_; f++){
	  systHist_[0][f]->SetName(sampleName_[f].c_str()+TString("_7Up"));
	  systHist_[1][f]->SetName(sampleName_[f].c_str()+TString("_7Down"));
	}
      }
      // End syst case
    } // End if do multidimensional shapes
 
  }
   
  if(produceOnly_ && !doMultiDimensionalShapes_)
    outputFile->Write();
  
  if(doMultiDimensionalShapes_){
    cout << "Now filling multiDim list of histos" << endl;
    
    vector<TH1*> tempMasterHist;
    vector<TString> tempMasterHistNames;
    tempMasterHist.clear();
    tempMasterHistNames.clear();
    // Store histograms
    for(size_t f=0; f<nSamples_+2; f++){
      tempMasterHist.push_back(  (TH1*) hist_[f]        ->Clone( TString( hist_[f]        ->GetName() + fitVars_[i]->getVarName())  ));
      tempMasterHist.push_back(  (TH1*) histStatUp_[f]  ->Clone( TString( histStatUp_[f]  ->GetName() + fitVars_[i]->getVarName())  ));
      tempMasterHist.push_back(  (TH1*) histStatDown_[f]->Clone( TString( histStatDown_[f]->GetName() + fitVars_[i]->getVarName())  ));
      
      tempMasterHistNames.push_back( sampleName_[f].c_str() 							     );
      tempMasterHistNames.push_back( sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatUp")   );
      tempMasterHistNames.push_back( sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatDown") );
    }
    
    // Syst case
    for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
      for(size_t f=0; f<nSamples_+2; f++){
	tempMasterHist.push_back( (TH1*) systHist_[a][f]->Clone( TString( systHist_[a][f]->GetName() + fitVars_[i]->getVarName() ) ));
	tempMasterHistNames.push_back( sampleName_[f].c_str()+systFancyComponents_[a] );
      }
    }
    // End syst loop
    cout << "Stored temp histograms" << endl;
    masterHist_.push_back( tempMasterHist );
    masterHistNames_.push_back( tempMasterHistNames );
    cout << "Stored master histograms" << endl;
  }
  // ---




if(!produceOnly_){ // For now raw indexes. This will improve.
    // 0: data
    // 1: WH
    // 2: HH
    // 3: DD
    // 4: other MCs
    TGraphAsymmErrors ddbkgBands;
    getErrorBands(*(hist_[3]), *ddbkgHistUp_, *ddbkgHistDown_, ddbkgBands);

    ddbkgBands.SetFillColor(1);
    ddbkgBands.SetFillStyle(3004);
    ddbkgBands.GetYaxis()->SetTitle("a.u.");
    ddbkgBands.GetYaxis()->SetTitleOffset(0.85);
    ddbkgBands.GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    ddbkgBands.GetXaxis()->SetTitleOffset(0.85);
    ddbkgBands.GetXaxis()->SetTitleSize(5);
    ddbkgBands.GetXaxis()->SetRange(displayMin_,displayMax_);    
    ddbkgBands.GetXaxis()->SetRangeUser(displayMin_,displayMax_);    


    TLegend* leg2_ = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg2_->SetTextFont(62);
    leg2_->SetBorderSize(0);
    leg2_->SetLineColor(1);
    leg2_->SetLineStyle(1);
    leg2_->SetLineWidth(1);
    leg2_->SetFillColor(0);
    leg2_->SetFillStyle(0);
    leg2_->AddEntry(hist_[1],"Base","l");
    leg2_->AddEntry(systHist_[0][1],"JESup","l");
    leg2_->AddEntry(systHist_[1][1],"JESdown","l");
    leg2_->AddEntry(&ddbkgBands, "stat+method","f");

    leg_ = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg_->SetTextFont(62);
    leg_->SetBorderSize(0);
    leg_->SetLineColor(1);
    leg_->SetLineStyle(1);
    leg_->SetLineWidth(1);
    leg_->SetFillColor(0);
    leg_->SetFillStyle(0);
    leg_->AddEntry(hist_[1],"Base","l");
    leg_->AddEntry(systHist_[0][1],"JESup","l");
    leg_->AddEntry(systHist_[1][1],"JESdown","l");

    hist_[1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[0][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[0][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[1][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[1][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    
    systHist_[0][1]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[1]->SetLineColor(1);
    systHist_[0][1]->SetLineColor(2);
    systHist_[1][1]->SetLineColor(3);
    systHist_[0][1]->Draw("hist");
    systHist_[1][1]->Draw("histsame");
    hist_[1]->Draw("histsame");
    leg_->Draw();
    canvas_->cd(); 
    // Order chosen to have good Y axis boundaries

    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsWH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsWH.png")).c_str());
    canvas_->Clear();
    hist_[1]->SetLineColor(sampleColour_[1]); // Reset line color from black to final one

    hist_[2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[0][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[0][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[1][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[1][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

    systHist_[0][2]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[2]->SetLineColor(1);
    systHist_[0][2]->SetLineColor(2);
    systHist_[1][2]->SetLineColor(3);
    systHist_[0][2]->Draw("hist");
    systHist_[1][2]->Draw("histsame");
    hist_[2]->Draw("histsame");
    leg_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsHH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsHH.png")).c_str());
    canvas_->Clear();
    hist_[2]->SetLineColor(sampleColour_[2]); // Reset line color from black to final one
    
    systHist_[0][3]->SetFillStyle(0);
    systHist_[1][3]->SetFillStyle(0);
    hist_[3]->SetFillStyle(0);

    hist_[3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[0][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[0][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[1][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[1][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

    systHist_[0][3]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[3]->SetMarkerColor(1);
    systHist_[0][3]->SetLineColor(2);
    systHist_[1][3]->SetLineColor(3);
    systHist_[0][3]->Draw("hist");
    systHist_[1][3]->Draw("histsame");
    hist_[3]->Draw("same");
    ddbkgBands.Draw("2same");
    leg2_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsDD.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariationsDD.png")).c_str());
    canvas_->Clear();
    hist_[3]->SetLineColor(sampleColour_[3]); // Reset line color from black to final one

    systHist_[0][3]->SetFillStyle(1001);
    systHist_[1][3]->SetFillStyle(1001);
    hist_[3]->SetFillStyle(1001);
    
    for(size_t f=4; f<nSamples_+2; f++){
      
      systHist_[0][f]->SetFillStyle(0);
      systHist_[1][f]->SetFillStyle(0);
      hist_[f]->SetFillStyle(0);
      
      hist_[f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      hist_[f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[0][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[0][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[1][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[1][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
 
      systHist_[0][f]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
      hist_[f]->SetLineColor(1);
      systHist_[0][f]->SetLineColor(2);
      systHist_[1][f]->SetLineColor(3);
      systHist_[0][f]->Draw("hist");
      systHist_[1][f]->Draw("histsame");
      hist_[f]->Draw("histsame");
      leg_->Draw();
      canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariations")+sampleName_[f].c_str()+string(".pdf")).c_str());
      canvas_->SaveAs((outFolder_+outputFileName+string("_jes_sigVariations")+sampleName_[f].c_str()+string(".png")).c_str());
      canvas_->Clear();
      
      hist_[f]->SetLineColor(sampleColour_[f]);
      systHist_[0][f]->SetFillStyle(1001);
      systHist_[1][f]->SetFillStyle(1001);
      hist_[f]->SetFillStyle(1001);

    }

    
    ////////////////////////////
    leg2_  = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg2_->SetTextFont(62);
    leg2_->SetBorderSize(0);
    leg2_->SetLineColor(1);
    leg2_->SetLineStyle(1);
    leg2_->SetLineWidth(1);
    leg2_->SetFillColor(0);
    leg2_->SetFillStyle(0);
    leg2_->AddEntry(hist_[1],"Base","l");
    leg2_->AddEntry(systHist_[0][1],"METup","l");
    leg2_->AddEntry(systHist_[1][1],"METdown","l");
    leg2_->AddEntry(&ddbkgBands, "stat+method","f"); // Dummy style show


    leg_ = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg_->SetTextFont(62);
    leg_->SetBorderSize(0);
    leg_->SetLineColor(1);
    leg_->SetLineStyle(1);
    leg_->SetLineWidth(1);
    leg_->SetFillColor(0);
    leg_->SetFillStyle(0);
    leg_->AddEntry(hist_[1],"Base","l");
    leg_->AddEntry(systHist_[2][1],"METup","l");
    leg_->AddEntry(systHist_[3][1],"METdown","l");

    hist_[1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[2][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[2][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[3][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[3][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

//    for(int mm=0; mm<signalHistWH_->GetNbinsX(); mm++)
//      cout << "bin: " << mm << ", diff up: " << signalHistWH_->GetBinContent(mm) - signalSystHistWH_[2]->GetBinContent(mm) << ", diff down: " << signalHistWH_->GetBinContent(mm) - signalSystHistWH_[3]->GetBinContent(mm) << endl;
//
    systHist_[2][1]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[1]->SetLineColor(1);
    systHist_[2][1]->SetLineColor(2);
    systHist_[3][1]->SetLineColor(3);
    systHist_[2][1]->Draw("hist");
    systHist_[3][1]->Draw("histsame");
    hist_[1]->Draw("histsame");
    leg_->Draw();
    canvas_->cd(); 
    // Order chosen to have good Y axis boundaries

    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsWH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsWH.png")).c_str());
    canvas_->Clear();
    hist_[1]->SetLineColor(sampleColour_[1]);

    hist_[2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[2][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[2][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[3][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[3][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

    systHist_[2][2]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[2]->SetLineColor(1);
    systHist_[2][2]->SetLineColor(2);
    systHist_[3][2]->SetLineColor(3);
    systHist_[2][2]->Draw("hist");
    systHist_[3][2]->Draw("histsame");
    hist_[2]->Draw("histsame");
    leg_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsHH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsHH.png")).c_str());
    canvas_->Clear();
    hist_[2]->SetLineColor(sampleColour_[2]);
    
    systHist_[2][3]->SetFillStyle(0);
    systHist_[3][3]->SetFillStyle(0);
    hist_[3]->SetFillStyle(0);

    hist_[3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[2][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[2][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[3][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[3][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

    systHist_[2][3]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[3]->SetMarkerColor(1);
    systHist_[2][3]->SetLineColor(2);
    systHist_[3][3]->SetLineColor(3);
    systHist_[2][3]->Draw("hist");
    systHist_[3][3]->Draw("histsame");
    hist_[3]->Draw("same");
    ddbkgBands.Draw("2same");
    leg2_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsDD.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariationsDD.png")).c_str());
    canvas_->Clear();
    hist_[3]->SetLineColor(sampleColour_[3]);
    
    systHist_[2][3]->SetFillStyle(1001);
    systHist_[3][3]->SetFillStyle(1001);
    hist_[3]->SetFillStyle(1001);


    for(size_t f=4; f<nSamples_+2; f++){
      
      systHist_[2][f]->SetFillStyle(0);
      systHist_[3][f]->SetFillStyle(0);
      hist_[f]->SetFillStyle(0);

      hist_[f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      hist_[f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[2][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[2][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[3][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[3][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
            
      systHist_[2][f]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
      hist_[f]->SetLineColor(1);
      systHist_[2][f]->SetLineColor(2);
      systHist_[3][f]->SetLineColor(3);
      systHist_[2][f]->Draw("hist");
      systHist_[3][f]->Draw("histsame");
      hist_[f]->Draw("histsame");
      leg_->Draw();
      canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariations")+sampleName_[f].c_str()+string(".pdf")).c_str());
      canvas_->SaveAs((outFolder_+outputFileName+string("_met_sigVariations")+sampleName_[f].c_str()+string(".png")).c_str());
      canvas_->Clear();
      
      hist_[f]->SetLineColor(sampleColour_[f]);
      systHist_[2][f]->SetFillStyle(1001);
      systHist_[3][f]->SetFillStyle(1001);
      hist_[f]->SetFillStyle(1001);

    }

    /////////////////////////

    leg2_  = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg2_->SetTextFont(62);
    leg2_->SetBorderSize(0);
    leg2_->SetLineColor(1);
    leg2_->SetLineStyle(1);
    leg2_->SetLineWidth(1);
    leg2_->SetFillColor(0);
    leg2_->SetFillStyle(0);
    leg2_->AddEntry(hist_[1],"Base","l");
    leg2_->AddEntry(systHist_[0][1],"JERup","l");
    leg2_->AddEntry(systHist_[1][1],"JERdown","l");
    leg2_->AddEntry(&ddbkgBands, "stat+method","f");

    leg_ = new TLegend(0.23,0.65,0.62,0.80,NULL,"brNDC");
    leg_->SetTextFont(62);
    leg_->SetBorderSize(0);
    leg_->SetLineColor(1);
    leg_->SetLineStyle(1);
    leg_->SetLineWidth(1);
    leg_->SetFillColor(0);
    leg_->SetFillStyle(0);
    leg_->AddEntry(hist_[1],"Base","l");
    leg_->AddEntry(systHist_[4][1],"JERup","l");
    leg_->AddEntry(systHist_[5][1],"JERdown","l");

    hist_[1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[4][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[4][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[5][1]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[5][1]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    
    systHist_[4][1]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[1]->SetLineColor(1);
    systHist_[4][1]->SetLineColor(2);
    systHist_[5][1]->SetLineColor(3);
    systHist_[4][1]->Draw("hist");
    systHist_[5][1]->Draw("histsame");
    hist_[1]->Draw("histsame");
    leg_->Draw();
    canvas_->cd(); 
    // Order chosen to have good Y axis boundaries

    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsWH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsWH.png")).c_str());
    canvas_->Clear();
    hist_[1]->SetLineColor(sampleColour_[1]);
    
    hist_[2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[4][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[4][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[5][2]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[5][2]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    

    systHist_[4][2]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[2]->SetLineColor(1);
    systHist_[4][2]->SetLineColor(2);
    systHist_[5][2]->SetLineColor(3);
    systHist_[4][2]->Draw("hist");
    systHist_[5][2]->Draw("histsame");
    hist_[2]->Draw("histsame");
    leg_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsHH.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsHH.png")).c_str());
    canvas_->Clear();
    hist_[2]->SetLineColor(sampleColour_[2]);
    
    systHist_[4][3]->SetFillStyle(0);
    systHist_[5][3]->SetFillStyle(0);
    hist_[3]->SetFillStyle(0);

    hist_[3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[4][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[4][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    systHist_[5][3]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    systHist_[5][3]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    
    systHist_[4][3]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
    hist_[3]->SetMarkerColor(1);
    systHist_[4][3]->SetLineColor(2);
    systHist_[5][3]->SetLineColor(3);
    systHist_[4][3]->Draw("hist");
    systHist_[5][3]->Draw("histsame");
    hist_[3]->Draw("same");
    ddbkgBands.Draw("2same");
    leg2_->Draw();
    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsDD.pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariationsDD.png")).c_str());
    canvas_->Clear();
    hist_[3]->SetLineColor(sampleColour_[3]);

    
    systHist_[4][3]->SetFillStyle(1001);
    systHist_[5][3]->SetFillStyle(1001);
    hist_[3]->SetFillStyle(1001);


    for(size_t f=4; f<nSamples_+2; f++){

      systHist_[4][f]->SetFillStyle(0);
      systHist_[5][f]->SetFillStyle(0);
      hist_[f]->SetFillStyle(0);

      hist_[f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      hist_[f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[4][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[4][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      systHist_[5][f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      systHist_[5][f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      
      systHist_[4][f]->GetXaxis()->SetTitle("R = p_{T}^{lead.track}/E^{#tau}");
      hist_[f]->SetLineColor(1);
      systHist_[4][f]->SetLineColor(2);
      systHist_[5][f]->SetLineColor(3);
      systHist_[4][f]->Draw("hist");
      systHist_[5][f]->Draw("histsame");
      hist_[f]->Draw("histsame");
      leg_->Draw();
      canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariations")+sampleName_[f].c_str()+string(".pdf")).c_str());
      canvas_->SaveAs((outFolder_+outputFileName+string("_jer_sigVariations")+sampleName_[f].c_str()+string(".png")).c_str());
      canvas_->Clear();
      
      hist_[f]->SetLineColor(sampleColour_[f]);
      systHist_[4][f]->SetFillStyle(1001);
      systHist_[5][f]->SetFillStyle(1001);
      hist_[f]->SetFillStyle(1001);
      
    }
    
    //// END SYST

    //  leg_ = new TLegend(0.3,0.635,0.63,0.93,NULL,"brNDC");
    leg_ = new TLegend(0.23,0.535,0.62,0.93,NULL,"brNDC");
    //  leg_ = new TLegend(0.7147651,0.6346154,0.9446309,0.9353147,NULL,"brNDC");
    //leg_ = new TLegend(0.75,0.6346154,1.,1.,NULL,"brNDC");
    
    leg_->SetTextFont(62);
    leg_->SetBorderSize(0);
    leg_->SetLineColor(1);
    leg_->SetLineStyle(1);
    leg_->SetLineWidth(1);
    leg_->SetFillColor(0);
    leg_->SetFillStyle(0);
    leg_->AddEntry(hist_[0],(fancySampleName_[0]).c_str(),"lep");
    leg_->AddEntry(hist_[1],fancySampleName_[1].c_str(),"l");
    //  leg_->AddEntry(signalHistHH_,signalSampleNameHH_.c_str(),"f");


    // Because colour printer is a nobrainer
    leg_->AddEntry(hist_[3],fancySampleName_[3].c_str(),"f");
    leg_->AddEntry(hist_[8],fancySampleName_[8].c_str(),"f");
    leg_->AddEntry(hist_[7],fancySampleName_[7].c_str(),"f");
    leg_->AddEntry(hist_[6],fancySampleName_[6].c_str(),"f");
    leg_->AddEntry(hist_[5],fancySampleName_[5].c_str(),"f");
    leg_->AddEntry(hist_[4],fancySampleName_[4].c_str(),"f");

//    leg_->AddEntry(hist_[3],fancySampleName_[3].c_str(),"f");
//    for(size_t f=4; f<nSamples_+2; f++) leg_->AddEntry(hist_[f],fancySampleName_[f].c_str(),"f");
    // for fig7 // for(size_t f=0; f<nMcSamples_; f++) leg_->AddEntry(mcbkgHist_[f],mcBkgFancySampleName_[f].c_str(),"f");
    canvas_->cd(); 
    // Order chosen to have good Y axis boundaries
   
    perMassPointSignalShapesToCompareHH_.push_back((TH1*)hist_[2]->Clone(hist_[2]->GetName() +TString("comparison") + massPointName_[currentMassPoint_].c_str()) );
    perMassPointSignalShapesToCompareWH_.push_back((TH1*)hist_[1]->Clone(hist_[1]->GetName() +TString("comparison") + massPointName_[currentMassPoint_].c_str()) );
    
    TH1* higgsH_ = 0;
    double cHiggsBR_ = 0.05; // Perhaps move to cfg file.
    double fhh(cHiggsBR_*cHiggsBR_) , fhw( 2*(1-cHiggsBR_)*cHiggsBR_), ftt(1-fhh-fhw);
    //  signalHistWH_->Scale(fhw/signalHistWH_->Integral());
    //  signalHistHH_->Scale(fhh/signalHistHH_->Integral());
    
    //  signalHistWH_->Add(signalHistHH_, fhh);
    higgsH_ = hist_[1];
    higgsH_->Scale(fhw/higgsH_->Integral());
    higgsH_->Add(hist_[2],fhh);
    
    perMassPointSignalShapesToCompare_.push_back((TH1*)higgsH_->Clone(higgsH_->GetName() + TString("totComparison") + massPointName_[currentMassPoint_].c_str() ) );    
    // Stacked drawing -------------------------------------------
    
    
    
    hist_[0]->GetYaxis()->SetTitle("a.u.");
    hist_[0]->GetYaxis()->SetTitleOffset(0.85);
    hist_[0]->GetXaxis()->SetTitle("p_{T}^{lead.track}/E^{#tau}");
    hist_[0]->GetXaxis()->SetTitleOffset(0.85);
    
    cout << " DATA HIST BIN 0: " << hist_[0]->GetBinContent(0) << " +/- " << hist_[0]->GetBinError(0) << endl;
    cout << " DATA HIST BIN 1: " << hist_[0]->GetBinContent(1) << " +/- " << hist_[0]->GetBinError(1) << endl;
    
    hist_[1]->GetYaxis()->SetTitle("a.u.");
    hist_[1]->GetYaxis()->SetTitleOffset(0.85);
    hist_[1]->GetXaxis()->SetTitle("p_{T}^{lead.track}/E^{#tau}");
    hist_[1]->GetXaxis()->SetTitleOffset(0.85);
    
    

    THStack hs("hs","stacked");
    
    hist_[2]->DrawNormalized("histsame");
    for(size_t f=3; f<nSamples_+2; f++){
      // for fig7 // for(size_t f=0; f<nMcSamples_; f++){
      hist_[f]->GetYaxis()->SetTitle("a.u.");
      hist_[f]->GetYaxis()->SetTitleOffset(0.85);
      hist_[f]->GetXaxis()->SetTitle("p_{T}^{lead.track}/E^{#tau}");
      hist_[f]->GetXaxis()->SetTitleOffset(0.85);
      hist_[f]->GetXaxis()->SetTitleSize(5);
      hist_[f]->GetXaxis()->SetRange(displayMin_,displayMax_);    
      hist_[f]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
      if(isDDbkg_[f] == 0) hs.Add(hist_[f],"hist");
    }
    
    for(size_t f=3; f<nSamples_+2; f++){
      if(isDDbkg_[f]) hs.Add(hist_[f],"hist");
    }
    
    //  dataHist_->SetMarkerStyle(1);
    //  dataHist_->SetMarkerSize(0.8);
    
 
    cout << "dd integral: " << hist_[3]->Integral() << endl;
    normalize(hs, 1.);
    hs.SetMaximum(0.4);
    hs.Draw("hist");
    hs.GetXaxis()->SetRange(displayMin_,displayMax_);    
    hs.GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    hs.GetYaxis()->SetTitle("a.u.");
    //  hs.GetYaxis()->SetTitleOffset(1.5);
    hs.GetXaxis()->SetTitle("p^{lead.track}/E^{#tau}");
    //  hs.GetXaxis()->SetTitleOffset(1.5);
    
    hist_[0]->Scale(1./hist_[0]->Integral());
    hist_[0]->GetXaxis()->SetRange(displayMin_,displayMax_);    
    hist_[0]->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
  
    hist_[0]->Draw("same");
    higgsH_->Scale(1./higgsH_->Integral());    /// ??? was signalHistWH_->Integral()); instead of higgsH->Integral());
    higgsH_->GetXaxis()->SetRange(displayMin_,displayMax_);    
    higgsH_->GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    higgsH_->Draw("histsame");
    
    TGraphErrors myBkgError;
    getErrorBands(hs, myBkgError);
    myBkgError.SetFillColor(1);
    myBkgError.SetFillStyle(3004);
    myBkgError.GetYaxis()->SetTitle("a.u.");
    myBkgError.GetYaxis()->SetTitleOffset(0.85);
    myBkgError.GetXaxis()->SetTitle("p_{T}^{lead.track}/E^{#tau}");
    myBkgError.GetXaxis()->SetTitleOffset(0.85);
    myBkgError.GetXaxis()->SetTitleSize(5);
    myBkgError.GetXaxis()->SetRange(displayMin_,displayMax_);    
    myBkgError.GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    myBkgError.Draw("2");
    
    TGraphErrors mySignalError;
    getErrorBands(*higgsH_, mySignalError);
    mySignalError.SetName("blahSignalError");
    mySignalError.SetFillColor(616);
    mySignalError.SetFillStyle(3005);
    mySignalError.GetYaxis()->SetTitle("a.u.");
    mySignalError.GetYaxis()->SetTitleOffset(0.85);
    mySignalError.GetXaxis()->SetTitle("p_{T}^{lead.track}/E^{#tau}");
    mySignalError.GetXaxis()->SetTitleOffset(0.85);
    mySignalError.GetXaxis()->SetTitleSize(5);
    mySignalError.GetXaxis()->SetRange(displayMin_,displayMax_);    
    mySignalError.GetXaxis()->SetRangeUser(displayMin_,displayMax_);    
    mySignalError.Draw("2");
    
    leg_->AddEntry(myBkgError.GetName(),"bkg total unc.","f");
    leg_->AddEntry(mySignalError.GetName(),"signal total unc.","f");
    
    TPaveText *pt1 = new TPaveText(0.17,0.45,0.65,0.5, "brNDC");
    pt1->SetBorderSize(1);
    pt1->SetFillColor(19);
    pt1->SetFillStyle(0);
    pt1->SetLineColor(0);
    pt1->SetTextFont(132);
    pt1->SetTextSize(0.055);
    //  TText *text = pt1->AddText("#splitline{m_{H^{#pm}} = 120 GeV/c^{2},}{BR(t #rightarrow H^{+}b) = 0.05}");
    TText *text = pt1->AddText("B(t #rightarrow H^{+}b) = 0.05");
    text->SetTextAlign(11);
    pt1->Draw("same");
    
    TPaveText *pt = new TPaveText(0.15,0.93,0.9,1.0, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillColor(19);
    pt->SetFillStyle(0);
    pt->SetLineColor(0);
    pt->SetTextFont(132);
    pt->SetTextSize(0.045);
    //TText *text = pt->AddText("#sqrt{s} = 7 TeV,  2.1 fb^{-1}  CMS Preliminary");
    //TText *text = pt->AddText("#sqrt{s} = 7 TeV,  2.1 fb^{-1} CMS ");
    //TText *text = pt->AddText("#sqrt{s} = 7 TeV,  2.2 fb^{-1} CMS Preliminary");
    //TText *text = pt->AddText("#sqrt{s} = 7 TeV,  1.9 fb^{-1}  CMS Preliminary");
    //   TText *text = pt->AddText("#sqrt{s} = 7 TeV,  1.9 fb^{-1}  CMS ");
    //   TText *text = pt->AddText("#sqrt{s} = 7 TeV,  4.0 fb^{-1} CMS Preliminary");
    //   TText *text = pt->AddText("#sqrt{s} = 7 TeV,  4.7 fb^{-1} CMS Preliminary");
    TText *textPrel = pt->AddText("#sqrt{s} = 7 TeV,  4.9 fb^{-1} CMS Preliminary");
    textPrel->SetTextAlign(11);
    pt->Draw("same");
    
    
    // End stacked drawing ---------------------------------------------
    
    // Normal drawing ----------------------------------------
    //  signalHistWH_->Draw("hist");
    //  //  signalHistHH_->DrawNormalized("histsame");
    //  
    //  ddbkgHist_->Draw("histsame");
    //  for(size_t f=0; f<nMcSamples_; f++)
    //    mcbkgHist_[f]->Draw("histsame");
    //  dataHist_->Draw("histsame");    
    // End normal drawing -----------------------------------
    
    leg_->Draw();
    
    canvas_->SaveAs((outFolder_+outputFileName+string(".pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string(".png")).c_str());
    canvas_->cd();
    canvas_->Clear();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    
    
  } // End if !produceOnly_
  

  
  if(produceOnly_ && !doMultiDimensionalShapes_)
     outputFile->Close();
  
  cout << "File closed" << endl;
  
}

void LandSShapesProducer::UnrollMultiDimensionalShape(){
  
  // Is 2D. Must translate into genericD

  masterShapes_.clear(); // Just to be sure
  
  cout << "Unrolling shapes" << endl;

  // Open file for storing shapes
  string outputFileName = outputFileName_ + string("_") + massPointName_[currentMassPoint_] + string("_") /* + string("multiDimShape_")*/ + fitVars_[0]->getVarName();// + string("_") + fitVars_[1]->getVarName();
  TFile* outputFile;
  if(produceOnly_ && doMultiDimensionalShapes_)
    outputFile = new TFile((outFolder_+outputFileName+string(".root")).c_str(), "RECREATE");
  cout << "File opened" << endl;


  
//  unrolled_         .clear();
//  unrolledStatUp_   .clear();
//  unrolledStatDown_ .clear();
//  unrolledSyst_     .clear();

  th2_.push_back( (TH2D*)th2_[nSamples_-1]->Clone( sampleName_[nSamples_].c_str()   + TString(th2_[nSamples_-1]->GetName() ))   ); // nSamples-1 Must always be di-bosons
  th2_.push_back( (TH2D*)th2_[nSamples_-1]->Clone( sampleName_[nSamples_+1].c_str() + TString(th2_[nSamples_-1]->GetName() )) ); // nSamples-1 Must always be di-bosons
  th2_[nSamples_]  ->Sumw2();
  th2_[nSamples_+1]->Sumw2();
  
  // Create unrolled stuff
  for(size_t f=0; f<nSamples_+2; f++){
    int nBins = th2_[f]->GetXaxis()->GetNbins() * th2_[f]->GetYaxis()->GetNbins();
    unrolled_.push_back( new TH1D(sampleName_[f].c_str(), sampleName_[f].c_str(),  nBins, 0 , nBins   ));
    unrolledStatUp_.push_back( new TH1D(sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatUp"),sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatUp"), nBins, 0 , nBins ) );
    unrolledStatDown_.push_back( new TH1D(sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatDown"), sampleName_[f].c_str()+TString("_")+sampleName_[f].c_str()+TString("_StatDown"),  nBins, 0, nBins  ) );

    
    unrolled_[f]->Sumw2();
    unrolledStatUp_[f]->Sumw2();
    unrolledStatDown_[f]->Sumw2();
  }
  
  // Syst case
  for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
    vector<TH1D*> unrolledSystTemp;
    unrolledSystTemp.clear();
    for(size_t f=0; f<nSamples_+2; f++){
      int nBins = th2_[f]->GetXaxis()->GetNbins() * th2_[f]->GetYaxis()->GetNbins();
      unrolledSystTemp.push_back( new TH1D(sampleName_[f].c_str()+systFancyComponents_[a], sampleName_[f].c_str()+systFancyComponents_[a], nBins, 0,  nBins  ));
      unrolledSystTemp[f]->Sumw2();

    }
    unrolledSyst_.push_back(unrolledSystTemp);
  } // End syst loop
    // End syst case
  
  // End create unrolled stuff
  
  // Fill statUp/Down

  // Stat filling and histogram putting nonzeroes instead of zeroes (RooFit mojo)
  // Perhaps do a single loop with systs.

 
  
  for(size_t f=0; f<nSamples_+2; f++){
    int ibin = 1;
    for(int xbin=1; xbin<th2_[f]->GetXaxis()->GetNbins(); xbin++){ // <= is for overflow
      for(int ybin=1; ybin<th2_[f]->GetYaxis()->GetNbins(); ybin++){ // <= is for overflow
	//	int ibin = th2_[f]->GetBin(xbin,ybin);
	Double_t val =th2_[f]->GetBinContent(xbin,ybin);
	unrolled_[f]->SetBinContent(ibin, val);
	unrolled_[f]->SetBinError(ibin, th2_[f]->GetBinError(xbin,ybin));
	cout << "histogram " << unrolled_[f]->GetName() << ": th2 has bin " << ibin << " - ["<<xbin<<","<<ybin<<"] th2:" <<  val << ", th1: " << unrolled_[f]->GetBinContent(ibin) << endl;
	if(unrolled_[f]->GetBinContent(ibin) == 0)
	  unrolled_[f]->SetBinContent(ibin, 1e-6); // Protection from RooStat/Fit breaking in combine
	unrolledStatUp_[f]  ->SetBinContent(ibin, std::max( th2_[f]->GetBinContent(xbin,ybin) + th2_[f]->GetBinError(xbin,ybin), 1e-6 ) );
	unrolledStatDown_[f]->SetBinContent(ibin, std::max( th2_[f]->GetBinContent(xbin,ybin) - th2_[f]->GetBinError(xbin,ybin), 1e-6 ) );

	ibin++;
      }
    }
  }

    
  unrolled_[nSamples_]->Scale(50.5/unrolled_[nSamples_]->Integral());      // Normalize to sample
  unrolled_[nSamples_+1]->Scale(0.4/unrolled_[nSamples_+1]->Integral()); // Normalize to sample

  
  // Rescale stat plots - they must have the same integral as the base ones
  for(size_t f=0; f<nSamples_+2; f++){
    if(isDDbkg_[f])
      unrolled_[f]->Scale(222./unrolled_[f]->Integral());
    
    unrolledStatUp_[f]  ->Scale(unrolled_[f]->Integral()/unrolledStatUp_[f]->Integral());
    unrolledStatDown_[f]->Scale(unrolled_[f]->Integral()/unrolledStatDown_[f]->Integral());
  }
  
  
  for(size_t a=0; a<nSysts_; a++){ // Loop on syst components
    th2Syst_[a].push_back( (TH2D*)th2Syst_[a][nSamples_-1]->Clone( sampleName_[nSamples_].c_str()  + systComponents_[a])   ); // nSamples-1 Must always be di-bosons
    th2Syst_[a].push_back( (TH2D*)th2Syst_[a][nSamples_-1]->Clone( sampleName_[nSamples_+1].c_str()+ systComponents_[a]) ); // nSamples-1 Must always be di-bosons
    th2Syst_[a][nSamples_]  ->Sumw2();
    th2Syst_[a][nSamples_+1]->Sumw2();

    for(size_t f=0; f<nSamples_+2; f++){
      int ibin=1;
      for(int xbin=1; xbin<th2_[f]->GetXaxis()->GetNbins(); xbin++){ // <= is for overflow
	for(int ybin=1; ybin<th2_[f]->GetYaxis()->GetNbins(); ybin++){ // <= is for overflow
	  unrolledSyst_[a][f]->SetBinContent(ibin, th2Syst_[a][f]->GetBinContent(xbin,ybin));
	  if(unrolledSyst_[a][f]->GetBinContent(ibin) == 0) // Protection from RooStat/Fit breaking in combine
	    unrolledSyst_[a][f]->SetBinContent(ibin, 1e-6);
	  ibin++;
	}
      }
    }
    unrolledSyst_[a][nSamples_]->Scale(50.5/unrolledSyst_[a][nSamples_]->Integral());      // Normalize to sample
    unrolledSyst_[a][nSamples_+1]->Scale(0.4/unrolledSyst_[a][nSamples_+1]->Integral()); // Normalize to sample

  }

  //////  // First cycle, to create histos
  //////  for(size_t iHist=0; iHist<masterHist_[0].size(); iHist++){
  //////    masterShapes_.push_back( (TH1*) masterHist_[0][iHist]->Clone() );
  //////    double nBins = masterHist_[0][iHist]->GetNbinsX() + masterHist_[1][iHist]->GetNbinsX();
  //////    cout << masterHist_[0][iHist]->GetName() << " has " << masterHist_[0][iHist]->GetNbinsX() << " bins, " << masterHist_[1][iHist]->GetName() << " has " << masterHist_[1][iHist]->GetNbinsX() << " bins. BinStep is " << fitVars_[0]->getBinStep() << endl;
  //////    masterShapes_[iHist]->SetBins(nBins, fitVars_[0]->getMin() , nBins * fitVars_[0]->getBinStep()); // 10 is personalized reset of binstep to 1 from 0.1
  //////    masterShapes_[iHist]->SetName( masterHistNames_[0][iHist] ); // the others are the same
  //////  }
  //////  
  //////  cout << "Master shapes created" << endl;
  //////  // Fill the bins (it is necessary to fill with the cloned one too, because SetBins destroys contents)
  //////  for(size_t iHist=0; iHist<masterHist_[0].size(); iHist++){
  //////    for (int ibin=0; ibin <= masterHist_[0][iHist]->GetNbinsX(); ibin++){
  //////      masterShapes_[iHist]->SetBinContent(ibin,       masterHist_[0][iHist]->GetBinContent(ibin) );
  //////      masterShapes_[iHist]->SetBinError(    ibin,     masterHist_[0][iHist]->GetBinError(ibin) );
  //////    }
  //////    for (int ibin=0; ibin <= masterHist_[1][iHist]->GetNbinsX(); ibin++){
  //////      masterShapes_[iHist]->SetBinContent(ibin + masterHist_[0][iHist]->GetNbinsX() , masterHist_[1][iHist]->GetBinContent(ibin) );
  //////      masterShapes_[iHist]->SetBinError(  ibin + masterHist_[0][iHist]->GetNbinsX() , masterHist_[1][iHist]->GetBinError(ibin) );
  //////    }
  //////
  //////    masterShapes_[iHist]->Scale(masterHist_[0][iHist]->Integral() / masterShapes_[iHist]->Integral() );
//////  }
//////  cout << "Master shapes bins are ok" << endl;
  
 //  masterHist_.clear(); // In order not to have that into the rootfile // See that





  if(produceOnly_ && doMultiDimensionalShapes_)
    outputFile->Write();
  
  outputFile->Close();
}

 void LandSShapesProducer::StorePerMassSignalShapes(){
  for(size_t i=0; i<perMassPointSignalShapesToCompare_.size(); i++){
    perMassPointSignalShapesToCompare_[i]  ->SetLineColor(currentMassPoint_);
    perMassPointSignalShapesToCompareHH_[i]->SetLineColor(currentMassPoint_);
    perMassPointSignalShapesToCompareWH_[i]->SetLineColor(currentMassPoint_);
  }
  signalShapesToCompare_  .push_back(perMassPointSignalShapesToCompare_  );
  signalShapesToCompareHH_.push_back(perMassPointSignalShapesToCompareHH_);
  signalShapesToCompareWH_.push_back(perMassPointSignalShapesToCompareWH_);

  cout << "Stored per mass signal shapes" << endl;
}

void LandSShapesProducer::DrawSignalShapesComparison(){
  cout << "Drawing signal shapes comparisons" << endl;
  canvas_->cd();
  
  leg_ = new TLegend(0.23,0.535,0.62,0.93,NULL,"brNDC");
  //  leg_ = new TLegend(0.7147651,0.6346154,0.9446309,0.9353147,NULL,"brNDC");
  //leg_ = new TLegend(0.75,0.6346154,1.,1.,NULL,"brNDC");
  leg_->SetTextFont(62);
  leg_->SetBorderSize(0);
  leg_->SetLineColor(1);
  leg_->SetLineStyle(1);
  leg_->SetLineWidth(1);
  leg_->SetFillColor(0);
  leg_->SetFillStyle(0);
  
  for(size_t s=0; s<signalShapesToCompare_.size(); s++){
    //    for(size_t i=0; i<signalShapesToCompare_[s].size(); i++)
    leg_->AddEntry(signalShapesToCompare_[s][0],(massPointName_[s]).c_str(),"l");
  }  
  
  // Draw them
  for(size_t i=0; i<nVars_; i++){
    canvas_->cd();
    for(size_t s=0; s<signalShapesToCompare_.size(); s++){     
      cout << "DEBUG: signalShapesToCompare: " << signalShapesToCompare_[s].size() << endl;
      signalShapesToCompare_[s][i]->SetLineColor(s);
      cout << "DEBUG: line color has been set" << endl;
      if(s=0) signalShapesToCompare_[s][i]->Draw("hist");
      else signalShapesToCompare_[s][i]->Draw("histsame");
      cout << "DEBUG: histograms have been drawn" << endl;
    }
    leg_->Draw();
    string outputFileName = outputFileName_ + string("_") + fitVars_[i]->getVarName();
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalShapesComparison")+string(".pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalShapesComparison")+string(".png")).c_str());
    canvas_->Clear();
    cout << "DEBUG: canvas saved and cleared" << endl;
  }

  for(size_t i=0; i<nVars_; i++){
    canvas_->cd();
    for(size_t s=0; s<signalShapesToCompareHH_.size(); s++){     
      signalShapesToCompareHH_[s][i]->SetLineColor(s);
      if(s=0) signalShapesToCompareHH_[s][i]->Draw("hist");
      else signalShapesToCompareHH_[s][i]->Draw("histsame");
    }
    leg_->Draw();
    string outputFileName = outputFileName_ + string("_") + fitVars_[i]->getVarName();
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalHHShapesComparison")+string(".pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalHHShapesComparison")+string(".png")).c_str());
    canvas_->Clear();
  }

  for(size_t i=0; i<nVars_; i++){
    canvas_->cd();
    for(size_t s=0; s<signalShapesToCompareWH_.size(); s++){     
      signalShapesToCompareWH_[s][i]->SetLineColor(s);
      if(s=0) signalShapesToCompareWH_[s][i]->Draw("hist");
      else signalShapesToCompareWH_[s][i]->Draw("histsame");
    }
    leg_->Draw();
    string outputFileName = outputFileName_ + string("_") + fitVars_[i]->getVarName();
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalWHShapesComparison")+string(".pdf")).c_str());
    canvas_->SaveAs((outFolder_+outputFileName+string("_signalWHShapesComparison")+string(".png")).c_str());
    canvas_->Clear();
  }

    
}

void LandSShapesProducer::Produce(){
  
  //cout << "INIT: signal tree entries:"<<signalTree_->GetEntries()<<endl;

  //fn for(size_t s=0; s< nMassPoints_; s++){
 
  size_t s = 1;
  InitMassPoint(s);
    
  cout << "=============================================================================" << endl;
  cout << "============================= MASS POINT: " << s << " ===============================" << endl;
  cout << "=============================================================================" << endl;
    for(size_t i = 0; i< nVars_; i++){



      InitPerVariableAmbient(i);

      cout << "============================= VARIABLE: " << i << ": init ambient completed ===============================" << endl;            
      BuildDatasets(i);
      cout << "============================= VARIABLE: " << i << ": datasets built ===============================" << endl;      
      //    BuildPDFs(i);
      
      if(!doMultiDimensionalShapes_)  DrawTemplates(i);


      //    BuildConstrainedModels(i);
      //    
      //    DoPerVariableFit(i);
      //    
      //    DrawPerVariableFit(i);
      //    
      //    DoPerVariableLikelihoodFit(i);
      
    }

   //fn  cout << "============================= MASS POINT: " << s << " ended processing of variables ===============================" << endl;
  //fn   if(doMultiDimensionalShapes_)
  //fn     UnrollMultiDimensionalShape();
  //fn   cout << "============================= MASS POINT: " << s << " unrolled stuff ===============================" << endl;
  //fn   //    StorePerMassSignalShapes();
  //fn }


  //  DoCombinedLikelihoodFit();
  //  DrawSignalShapesComparison();
  //}
  //  cout.rdbuf(old); // restore   
}
