#include "LIP/TopTaus/interface/LandSShapesProducer.hh"

// System includes
#include <sstream>

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT includes


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

  // Clear vectors
  nVars_ = 0;
  nMcSamples_ = 0;
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
  
  // Init RooFit variables
  sigVar_             = 0;              
  sigMeanVar_         = 0;              
  sigSigmaVar_        = 0;              
  ttbarmcbkgVar_      = 0;
  ttbarmcbkgMeanVar_  = 0;
  ttbarmcbkgSigmaVar_ = 0;
  mcbkgVar_           = 0;               
  mcbkgMeanVar_       = 0;         
  mcbkgSigmaVar_      = 0;         
  ddbkgVar_           = 0;              
  ddbkgMeanVar_       = 0;         
  ddbkgSigmaVar_      = 0;        

  identifier_ = "";

  cout << "Getting Parameter sets" << endl;

  // Get ParameterSet from cfg file
  const edm::ParameterSet &mFitPars = edm::readPSetsFrom(parSet_)->getParameter<edm::ParameterSet>("LandSShapesProducerParSet");

  outFolder_        = mFitPars.getParameter<std::string>("outFolder");
  outputFileName_  = mFitPars.getParameter<std::string>("outputFileName");
  outputFileNameSuffix_  = mFitPars.getParameter<vector<std::string> >("outputFileNameSuffix");
  nMassPoints_ = outputFileNameSuffix_.size();
  baseMCDir_        = mFitPars.getParameter<std::string>("baseMCDir");
  baseDataDir_      = mFitPars.getParameter<std::string>("baseDataDir");
  
  signalFileNameWH_   = mFitPars.getParameter<vector<std::string> >("signalFileNameWH");
  signalFileNameHH_   = mFitPars.getParameter<vector<std::string> >("signalFileNameHH");
  ddBkgFileName_      = mFitPars.getParameter<std::string>("ddBkgFileName");
  vector<string>mcBkgFileNameTemp = mFitPars.getParameter<vector<std::string> >("mcBkgFileName");
  for(size_t f=0; f<mcBkgFileNameTemp.size(); f++)
    mcBkgFileName_.push_back( mcBkgFileNameTemp[f] );
  dataFileName_       = mFitPars.getParameter<std::string>("dataFileName");

  signalSampleNameWH_   = mFitPars.getParameter<std::string>("signalSampleNameWH");
  signalSampleNameHH_   = mFitPars.getParameter<std::string>("signalSampleNameHH");
  ddBkgSampleName_      = mFitPars.getParameter<std::string>("ddBkgSampleName");
  mcBkgSampleName_      = mFitPars.getParameter<vector<std::string> >("mcBkgSampleName");
//  vector<string>mcBkgSampleNameTemp = mFitPars.getParameter<vector<std::string> >("mcBkgSampleName");
//  for(size_t f=0; f<mcBkgSampleNameTemp.size(); f++)
//    mcBkgSampleName_.push_back( mcBkgSampleNameTemp[f] );
  dataSampleName_       = mFitPars.getParameter<std::string>("dataSampleName");

  
  minitreeSelected_   = mFitPars.getParameter<std::string>("minitreeSelected");
  minitreeDataDriven_ = mFitPars.getParameter<std::string>("minitreeDataDriven");
  
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
  nMcSamples_ = mcBkgFileName_.size();
  
  
  cout << "Files opened" << endl;
  // Set variables
  //
  nVars_ = vars_.size();
  
  
  for(size_t i=0; i<nVars_; i++)
    {
      FitVar myVar(vars_[i], mins_[i], maxs_[i], bins_[i], hmin_[i], hmax_[i], unbinned_[i], smoothOrder_[i]);
      fitVars_.push_back(myVar);
    }

  // Set canvas
  canvas_ = new TCanvas("canvas","My plots ",0,0,1000,500);
  canvas_->cd();
  
  resultsFile_.open ((outFolder_+resultsFileName_).c_str());
  
  //  Uncomment following line in order to redirect stdout to file
  //  streamToFile_ = std::cout.rdbuf(resultsFile_.rdbuf());
  
  cout << "Init process complete" << endl;
}

void LandSShapesProducer::SetOptions(){
  myStyle_->SetOptStat(0);
}

void LandSShapesProducer::InitMassPoint(size_t s){

  mcBkgFile_.clear();

  currentMassPoint_ = s;

  signalFileWH_   = TFile::Open(baseMCDir_   + signalFileNameWH_[currentMassPoint_]  ); signalTreeWH_   = (TTree*) signalFileWH_  ->Get(minitreeSelected_);
  signalFileHH_   = TFile::Open(baseMCDir_   + signalFileNameHH_[currentMassPoint_]  ); signalTreeHH_   = (TTree*) signalFileHH_  ->Get(minitreeSelected_);

  ddBkgFile_      = TFile::Open(baseDataDir_   + ddBkgFileName_     ); ddBkgTree_      = (TTree*) ddBkgFile_     ->Get(minitreeDataDriven_);
  for(size_t f=0; f<nMcSamples_; f++){ mcBkgFile_.push_back( TFile::Open(baseMCDir_   + mcBkgFileName_[f]     ) ); mcBkgTree_.push_back( (TTree*) mcBkgFile_[f]     ->Get(minitreeSelected_) );}
  dataFile_       = TFile::Open(baseDataDir_ + dataFileName_      ); dataTree_       = (TTree*) dataFile_      ->Get(minitreeSelected_);
  

}

void LandSShapesProducer::InitPerVariableAmbient(size_t i){
  
  myMCBkgDSName_.clear();
  
  // Binned fit variable
  myvar_           = new RooRealVar(fitVars_[i].getVarName().c_str(), fitVars_[i].getVarName().c_str(), fitVars_[i].getMin(), fitVars_[i].getMax());  myvar_->setBins(fitVars_[i].getBins()); 
  myvar_weights_   = new RooRealVar("weight","weight",0,1000);
  isOSvar_         = new RooRealVar("is_os","is_os",0,2);
  
  
  //Define data sets /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  mySignalDSNameWH_= fitVars_[i].getVarName().c_str() + signalSampleNameWH_;
  mySignalDSNameHH_= fitVars_[i].getVarName().c_str() + signalSampleNameHH_;
  myDDBkgDSName_ = fitVars_[i].getVarName().c_str() + ddBkgSampleName_;
  for(size_t f=0; f<nMcSamples_; f++){
    myMCBkgDSName_.push_back( fitVars_[i].getVarName().c_str() + mcBkgSampleName_[f] );
  }
  myDataDSName_   = fitVars_[i].getVarName().c_str() + dataSampleName_;
  
  mySignalDSWH_        = 0; //= new RooDataSet(mySignalDSName.c_str(),mySignalDSName.c_str(), signalTree_, RooArgSet(*myvar,*myvar_weights),0,"weight" );
  mySignalDSHH_        = 0; //= new RooDataSet(mySignalDSName.c_str(),mySignalDSName.c_str(), signalTree_, RooArgSet(*myvar,*myvar_weights),0,"weight" );
  myDDBkgDS_      = 0; 
  myMCBkgDS_.clear();
  myDataDS_       = 0; 
  
  signalHistWH_     = 0;
  signalHistHH_     = 0;
  ddbkgHist_      = 0;
  mcbkgHist_.clear();

  signalHistoWH_     = 0;
  signalHistoHH_     = 0;
  ddbkgHisto_      = 0;
  mcbkgHisto_.clear();

  leg_            = 0;
  
  
  sumWeights_ = 0;
}

void LandSShapesProducer::BuildDatasets(size_t i){
  
  // Temp variables for setting branch addresses
  double myVarAllocator, myVarWeightAllocator;
  double isOSsig;
  
  mySignalDSWH_         = new RooDataSet(mySignalDSNameWH_.c_str(),mySignalDSNameWH_.c_str(),              RooArgSet(*myvar_,*myvar_weights_), "weight"); // This constructor does not accept the cut parameter  
  // Get WH events
  signalTreeWH_->SetBranchAddress(fitVars_[i].getVarName().c_str(), &myVarAllocator);
  signalTreeWH_->SetBranchAddress("weight", &myVarWeightAllocator);
  signalTreeWH_->SetBranchAddress("is_os", &isOSsig);
  for(unsigned int ev=0; ev<signalTreeWH_->GetEntries(); ev++){
    signalTreeWH_->GetEntry(ev);
    if(isOSsig<0.5) continue;
    myvar_->setVal(myVarAllocator);
    sumWeights_ += myVarWeightAllocator;
    myvar_weights_->setVal(myVarWeightAllocator);
    mySignalDSWH_->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
  }
  mySignalDSHH_         = new RooDataSet(mySignalDSNameHH_.c_str(),mySignalDSNameHH_.c_str(),              RooArgSet(*myvar_,*myvar_weights_), "weight"); // This constructor does not accept the cut parameter
  // Get HH events
  signalTreeHH_->SetBranchAddress(fitVars_[i].getVarName().c_str(), &myVarAllocator);
  signalTreeHH_->SetBranchAddress("weight", &myVarWeightAllocator);
  signalTreeHH_->SetBranchAddress("is_os", &isOSsig);
  //    cout << "getIsoS      ";
  for(unsigned int ev=0; ev<signalTreeHH_->GetEntries(); ev++){
    signalTreeHH_->GetEntry(ev);
    if(isOSsig < 0.5) continue;
    myvar_->setVal(myVarAllocator);
    sumWeights_ += myVarWeightAllocator;
    myvar_weights_->setVal(myVarWeightAllocator);
    mySignalDSHH_->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
  }
  
  myDDBkgDS_         = new RooDataSet(myDDBkgDSName_.c_str(),myDDBkgDSName_.c_str(),              RooArgSet(*myvar_,*myvar_weights_), "weight"); // This constructor does not accept the cut parameter
  // Get DD events
  ddBkgTree_->SetBranchAddress(fitVars_[i].getVarName().c_str(), &myVarAllocator);
  ddBkgTree_->SetBranchAddress("weight", &myVarWeightAllocator);
  ddBkgTree_->SetBranchAddress("is_os", &isOSsig);
  //    cout << "getIsoS      ";
  for(unsigned int ev=0; ev<ddBkgTree_->GetEntries(); ev++){
    ddBkgTree_->GetEntry(ev);
    if(isOSsig < 0.5) continue;
    myvar_->setVal(myVarAllocator);
    sumWeights_ += myVarWeightAllocator;
    myvar_weights_->setVal(myVarWeightAllocator);
    myDDBkgDS_->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
  }

  for(size_t f=0; f<nMcSamples_; f++){
    myMCBkgDS_.push_back( new RooDataSet(myMCBkgDSName_[f].c_str(),myMCBkgDSName_[f].c_str(),              RooArgSet(*myvar_,*myvar_weights_), "weight") ); // This constructor does not accept the cut parameter
    // Get MCBkg events
    mcBkgTree_[f]->SetBranchAddress(fitVars_[i].getVarName().c_str(), &myVarAllocator);
    mcBkgTree_[f]->SetBranchAddress("weight", &myVarWeightAllocator);
    mcBkgTree_[f]->SetBranchAddress("is_os", &isOSsig);
    //    cout << "getIsoS      ";
    for(unsigned int ev=0; ev<mcBkgTree_[f]->GetEntries(); ev++){
      mcBkgTree_[f]->GetEntry(ev);
      if(isOSsig < 0.5) continue;
      myvar_->setVal(myVarAllocator);
      sumWeights_ += myVarWeightAllocator;
      myvar_weights_->setVal(myVarWeightAllocator);
      myMCBkgDS_[f]->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
    }
  }
  
  myDataDS_         = new RooDataSet(myDataDSName_.c_str(),myDataDSName_.c_str(),              RooArgSet(*myvar_,*myvar_weights_), "weight"); // This constructor does not accept the cut parameter
  // Get Data events
  dataTree_->SetBranchAddress(fitVars_[i].getVarName().c_str(), &myVarAllocator);
  dataTree_->SetBranchAddress("weight", &myVarWeightAllocator);
  dataTree_->SetBranchAddress("is_os", &isOSsig);
  //    cout << "getIsoS      ";
  for(unsigned int ev=0; ev<dataTree_->GetEntries(); ev++){
    dataTree_->GetEntry(ev);
    if(isOSsig < 0.5) continue;
    myvar_->setVal(myVarAllocator);
    sumWeights_ += myVarWeightAllocator;
    myvar_weights_->setVal(myVarWeightAllocator);
    myDataDS_->add(RooArgSet(*myvar_,*myvar_weights_),myVarWeightAllocator);
  }
  
  cout << "Datasets built" << endl;
  
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
  signalHistoWH_ = mySignalDSWH_->binnedClone();
  signalHistoHH_ = mySignalDSHH_->binnedClone();
  ddbkgHisto_    = myDDBkgDS_->binnedClone();
  for(size_t f=0; f<nMcSamples_; f++) mcbkgHisto_.push_back( myMCBkgDS_[f]->binnedClone());
  dataHisto_  = myDataDS_ ->binnedClone();

  cout << "Binned clones built" << endl;

  cout << mySignalDSNameWH_ << " unbinned entries: " << mySignalDSWH_->numEntries() << ". weighted entries: " << signalHistoWH_->sum(kFALSE) << endl;
  cout << mySignalDSNameHH_ << " unbinned entries: " << mySignalDSHH_->numEntries() << ". weighted entries: " << signalHistoHH_->sum(kFALSE) << endl;
  cout << myDDBkgDSName_ << " unbinned entries: " << myDDBkgDS_->numEntries() << ". weighted entries: " << ddbkgHisto_->sum(kFALSE) << endl;
  for(size_t f=0; f<nMcSamples_; f++)  cout << myMCBkgDSName_[f] << " unbinned entries: " << myMCBkgDS_[f]->numEntries() << ". weighted entries: " << mcbkgHisto_[f]->sum(kFALSE) << endl;
  cout << myDataDSName_ << " unbinned entries: " << myDataDS_->numEntries() << ". weighted entries: " << dataHisto_->sum(kFALSE) << endl;
  
  cout << "BuildDatasets successful" << endl;
}

//void LandSShapesProducer::BuildPDFs(size_t i){
//  
//  switch(fitVars_[i].getUnbinned()){
//  case 1 : // Unbinned
//    if(includeSignal_) u_signalModel_ = new RooKeysPdf(signalModelName_.c_str(), signalModelName_.c_str(), *myvar_, *mySignalDS_, RooKeysPdf::MirrorBoth,2);
//    u_ddbkgModel_  = new RooKeysPdf(ddbkgModelName_.c_str(),  ddbkgModelName_.c_str(),  *myvar_, *myDDBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    if(standaloneTTbar_) u_ttbarmcbkgModel_  = new RooKeysPdf(ttbarmcbkgModelName_.c_str(),  ttbarmcbkgModelName_.c_str(),  *myvar_, *myTTBARMCBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    u_mcbkgModel_  = new RooKeysPdf(mcbkgModelName_.c_str(),  mcbkgModelName_.c_str(),  *myvar_, *myMCBkgDS_ , RooKeysPdf::MirrorBoth,2);
//    break;
//  case 0:  // Binned (w/ or w/out smoothing)
//    if(includeSignal_) b_signalModel_ = new RooHistPdf(signalModelName_.c_str(), signalModelName_.c_str(), *myvar_, *signalHisto_, fitVars_[i].getSmoothOrder());
//    b_ddbkgModel_  = new RooHistPdf(ddbkgModelName_.c_str(),  ddbkgModelName_.c_str(), *myvar_, *ddbkgHisto_ , fitVars_[i].getSmoothOrder());
//    if(standaloneTTbar_) b_ttbarmcbkgModel_ = new RooHistPdf(ttbarmcbkgModelName_.c_str(),  ttbarmcbkgModelName_.c_str(),  *myvar_, *ttbarmcbkgHisto_  , fitVars_[i].getSmoothOrder()); 
//    b_mcbkgModel_  = new RooHistPdf(mcbkgModelName_.c_str(),  mcbkgModelName_.c_str(),  *myvar_, *mcbkgHisto_  , fitVars_[i].getSmoothOrder()); 
//    break;
//  default : // Dummy - should never arrive here
//    cout<<"Neither binned not unbinned. Check your options motherfucker."<<endl;
//  }
//}

void LandSShapesProducer::DrawTemplates(size_t i){
  
  string outputFileName = outputFileName_ + string("_") + outputFileNameSuffix_[currentMassPoint_] + fitVars_[i].getVarName() + string(".root");
  TFile* outputFile = new TFile((outFolder_+outputFileName).c_str(), "RECREATE");
  
  canvas_->cd();
  canvas_->Clear();
  
  //signal histogram ////////////////////////////////////////////////////////////
  signalHistWH_ = signalHistoWH_->createHistogram(fitVars_[i].getVarName().c_str(),fitVars_[i].getBins() );
  signalHistHH_ = signalHistoHH_->createHistogram(fitVars_[i].getVarName().c_str(),fitVars_[i].getBins() );

  //  signalHist_->SetOption("0000");
  //  signalHist_->SetLineWidth(3);
  //  signalHist_->SetTitle("");
  //  signalHist_->GetYaxis()->SetTitle("a.u.");
  //  signalHist_->GetYaxis()->SetTitleOffset(1.5);
  //  signalHist_->SetLineColor(kGreen);
  //  signalHist_->SetFillColor(kGreen);
  //  ///////////////////////////////////////////////////////////////////////////////
  
  
  // dd bkg histogram /////////////////////////////////////////////////
  ddbkgHist_ = ddbkgHisto_->createHistogram(fitVars_[i].getVarName().c_str(),fitVars_[i].getBins() );
//  ddbkgHist_->SetLineColor(kRed);
//  ddbkgHist_->SetFillColor(kRed);
//  ddbkgHist_->SetLineWidth(3);
//  ddbkgHist_->SetFillStyle(3017);
//  /////////////////////////////////////////////////////////////////////
  
  
  // mc bkg histogram ////////////////////////////////////////////////
  for(size_t f=0; f<nMcSamples_; f++)
    mcbkgHist_.push_back( mcbkgHisto_[f]->createHistogram(fitVars_[i].getVarName().c_str(),fitVars_[i].getBins() ) );
  //    mcbkgHist_->SetLineColor(kBlack);
  //    mcbkgHist_->SetFillColor(kBlack);
  //    mcbkgHist_->SetLineWidth(3);
  //    mcbkgHist_->SetFillStyle(3017);
  //    ///////////////////////////////////////////////////////////////////
    
  mcbkgHist_.push_back( (TH1*)mcbkgHist_[0]->Clone(  mcBkgSampleName_[nMcSamples_].c_str()) );
  mcbkgHist_.push_back( (TH1*)mcbkgHist_[0]->Clone(  mcBkgSampleName_[nMcSamples_+1].c_str()) );
  mcbkgHist_[nMcSamples_]->Reset();
  mcbkgHist_[nMcSamples_+1]->Reset();





  dataHist_ = dataHisto_->createHistogram(fitVars_[i].getVarName().c_str(),fitVars_[i].getBins() );
  
  leg_ = new TLegend(0.3,0.665,0.85,0.86,NULL,"brNDC");
  leg_->SetTextFont(62);
  leg_->SetBorderSize(0);
  leg_->SetLineColor(1);
  leg_->SetLineStyle(1);
  leg_->SetLineWidth(1);
  leg_->SetFillColor(0);
  leg_->SetFillStyle(1001);
  leg_->AddEntry(signalHistWH_,signalSampleNameWH_.c_str(),"f");
  leg_->AddEntry(signalHistHH_,signalSampleNameHH_.c_str(),"f");
  leg_->AddEntry(ddbkgHist_,ddBkgSampleName_.c_str(),"f");
  for(size_t f=0; f<nMcSamples_; f++) leg_->AddEntry(mcbkgHist_[f],mcBkgSampleName_[f].c_str(),"f");
  leg_->AddEntry(dataHist_,dataSampleName_.c_str(),"f");
  canvas_->cd(); 
  // Order chosen to have good Y axis boundaries

  signalHistWH_->DrawNormalized("hist");
  signalHistHH_->DrawNormalized("histsame");
  ddbkgHist_->DrawNormalized("histsame");
  for(size_t f=0; f<nMcSamples_; f++)
    mcbkgHist_[f]->DrawNormalized("histsame");
  dataHist_->DrawNormalized("histsame");    
  leg_->Draw();
  canvas_->SaveAs((outFolder_+identifier_+string("_shapes")+string(".pdf")).c_str());
  canvas_->SaveAs((outFolder_+identifier_+string("_shapes")+string(".png")).c_str());
  canvas_->cd();
  canvas_->Clear();
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  signalHistWH_->SetName(signalSampleNameWH_.c_str());
  signalHistHH_->SetName(signalSampleNameHH_.c_str());

  ddbkgHist_->SetName(ddBkgSampleName_.c_str());
  for(size_t f=0; f<nMcSamples_; f++)
    mcbkgHist_[f]->SetName(mcBkgSampleName_[f].c_str());
  mcbkgHist_[nMcSamples_]->SetName(mcBkgSampleName_[nMcSamples_].c_str());
  mcbkgHist_[nMcSamples_+1]->SetName(mcBkgSampleName_[nMcSamples_+1].c_str());
  
  dataHist_->SetName(dataSampleName_.c_str());
  
  outputFile->Write();
  outputFile->Close();
  
}


void LandSShapesProducer::Produce(){
  
  //cout << "INIT: signal tree entries:"<<signalTree_->GetEntries()<<endl;

  for(size_t s=0; s< nMassPoints_; s++){
    
    InitMassPoint(s);
    
    for(size_t i = 0; i< nVars_; i++){
      
      InitPerVariableAmbient(i);
      
      BuildDatasets(i);
      
      //    BuildPDFs(i);
      
      DrawTemplates(i);
      
      //    BuildConstrainedModels(i);
      //    
      //    DoPerVariableFit(i);
      //    
      //    DrawPerVariableFit(i);
      //    
      //    DoPerVariableLikelihoodFit(i);
      
    }
  }
  //  DoCombinedLikelihoodFit();
  
  //}
  //  cout.rdbuf(old); // restore   
}
