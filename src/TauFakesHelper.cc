
/**   
      \class	    TauFakesHelper TauFakesHelper.cc "UserCode/LIP/TopTaus/src/TauFakesHelper.cc"                                                                     
      \brief    Class for managing the operations needed in order to estimate tauFakes :)
      
      \author   Pietro Vischia
      
      \version  $Id: TauFakesHelper.cc,v 1.5 2013/04/23 10:51:16 vischia Exp $                                                                                                       
*/




/*

   FRDataQCDTree q(20);
   q.processEvents();

*/
 
#include "LIP/TopTaus/interface/TauFakesHelper.hh"


// TMVA
//#include "PhysicsTools/MVAComputer/interface/MVAModuleHelper.h"
//#include "RecoTauTag/TauTagTools/interface/TauMVADBConfiguration.h"
//
//#include "PhysicsTools/MVAComputer/interface/MVAComputerRecord.h"
//#include "CondFormats/PhysicsToolsObjects/interface/MVAComputer.h"
#include "PhysicsTools/MVAComputer/interface/Variable.h"
#include "PhysicsTools/MVAComputer/interface/MVAComputer.h"
#include "PhysicsTools/MVATrainer/interface/MVATrainer.h"
#include "PhysicsTools/MVAComputer/interface/TreeReader.h"
//#include "PhysicsTools/MVATrainer/interface/MVATrainerLooper.h"
//#include "PhysicsTools/MVATrainer/interface/MVATrainerContainerSave.h"

//#include "FWCore/Framework/interface/IOVSyncValue.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"


// ROOT includes
#include <TLegend.h>
#include <TLegendEntry.h>

using namespace event; // FIXME: define inside the code to avoid ambiguities
using namespace PhysicsTools;

TauFakesHelper::TauFakesHelper(double tauPtCut,
			       TString inputArea,
			       TString trainingTreesArea,
			       TString configArea,
			       TString trainingOutputArea,
			       TString outputArea,
			       TString puFileName,
			       TString ntuplesArea):
  ObjectSelector(tauPtCut),
  tauPtCut_            (tauPtCut	    ),   
  inputArea_		(inputArea	    ),   
  trainingTreesArea_	(trainingTreesArea  ),
  configArea_		(configArea	    ),   
  trainingOutputArea_	(trainingOutputArea ),
  outputArea_		(outputArea	    ),   
  puFileName_         	(puFileName         ),
  ntuplesArea_          (ntuplesArea        ),
  PlotStyle()
{

  setTDRStyle();
  
}



void TauFakesHelper::SetParameters(){
  
  switch(qualifier_){
  case WMUDATA:
    infix_ = "WMuData";
    break;
  case WMUMC:
    infix_ = "WMuMC";
    break;
  case QCDDATA:
    infix_ = "QCDData";
    break;
  case QCDMC:
    infix_ = "QCDMC";
    break;
  default: 
    cout << "No valid qualifier specified."<< endl;
    exit(-1);
    break;
  }


  // FRWJetsTree
  if(qualifier_ == WMUDATA || qualifier_ == WMUMC){
    
    DRMIN_JET_E_ = 0.5;
    DRMIN_JET_M_ = 0.5;
    DRMIN_T_E_   = 0.5;
    DRMIN_T_M_   = 0.5;
    DRMIN_T_J_   = 0.3;
    
    MET_CUT_   = 50.0;
    
    std::vector<float> MCPUDist, DataPUDist;
    float dataDist_2011A[50] = {0.00592951, 0.0255428, 0.0591468, 0.097016, 0.126287, 0.138848, 0.134117, 0.11691, 0.0937398, 0.0700927, 0.0493627, 0.0329741, 0.0209976, 0.0127917, 0.00747402, 0.00419649, 0.00226774, 0.00118102, 0.000593481, 0.000288109, 0.000135272, 6.14976e-05, 2.71017e-05, 1.15906e-05, 4.81577e-06};

    TFile* dataWeightsFile = new TFile(puFileName_);
    if(!dataWeightsFile){ cout<<endl<<"File : "<< puFileName_ << " not found .. "<<endl; exit(0);}
    
    TH1D * hist = (TH1D *) dataWeightsFile->Get("pileup");
    
    for(int i=0; i<50;i++){ dataDist_2011A[i]=hist->GetBinContent(i+1); DataPUDist.push_back(dataDist_2011A[i]);}    
    
    TFile* mcwjetsFile = new TFile(inputArea_+TString("/WJetsToLNuForFakes.root"));
    TH1 * mcPileupHisto = (TH1*) mcwjetsFile->Get("myEventSelector/Selection/pileup");
    int histoEntries(51);
    if( mcPileupHisto ){  for( int s=1; s < mcPileupHisto->GetNbinsX() && s<histoEntries; s++ ){ MCPUDist.push_back( mcPileupHisto->GetBinContent(s)); } }
    LumiWeights_ = reweight::LumiReWeighting(MCPUDist, DataPUDist);
  } // if is WMu
  else
    {
      
      DRMIN_JET_E_ = 0.3;
      DRMIN_JET_M_ = 0.3;
      DRMIN_T_E_   = 0.3;
      DRMIN_T_M_   = 0.3;
      DRMIN_T_J_   = 0.3;
      MIN_JET_PT_CUT_ = 20.0;
      MIN_TAU_PT_CUT_ = 20.0;
      
      TAUPRONGS_ = 0;    
      
      JES_=0.0;
    }
}



void TauFakesHelper::ProcessEvents(unsigned int qualifier){
  qualifier_ = qualifier;
  SetParameters();

  switch(qualifier_){
  case WMUDATA:
    ComputeFakeRate("PFlow", true, true);
    ComputeFakeRate("PFlow", false, true);
    break;
  case WMUMC:
    ComputeFakeRate("PFlow", true, false);
    ComputeFakeRate("PFlow", false, false);
    break;
  case QCDDATA:
    ComputeFakeRate("PFlow", true);
    ComputeFakeRate("PFlow", false);
    break;
  case QCDMC:
    ComputeFakeRate("PFlow", true);
    ComputeFakeRate("PFlow", false);
    break;
  default: 
    cout << "No valid qualifier specified."<< endl;
    exit(-1);
    break;
  }


}

//  switch(qualifier_){
//  case WMUDATA:
//    break;
//  case WMUMC:
//    break;
//  case QCDDATA:
//    break;
//  case QCDMC:
//    break;
//  default: 
//    cout << "No valid qualifier specified."<< endl;
//    exit(-1);
//    break;
//  }

void TauFakesHelper::ComputeFakeRate(TString myKey, bool passing, bool isAntiBTag){
  TString url("");
  bool isDATA(true);
  TString histoFiles_("");
  if(passing) histoFiles_ = trainingTreesArea_+"/TauFR"+infix_+"Passing";
  else        histoFiles_ = trainingTreesArea_+"/TauFR"+infix_+"Failing";
  
  switch(qualifier_){
  case WMUDATA:
    url = inputArea_+"/MuonForFakes.root";
    isDATA = true;
    break;
  case WMUMC:
    url = inputArea_+"/WJetsToLNuForFakes.root";
    isDATA = false;
    break;
  case QCDDATA:
    url = inputArea_+TString("JetsForFakes.root");
    isDATA = true;
    break;
  case QCDMC:
    isDATA = false;
    break;
  default: 
    cout << "No valid qualifier specified."<< endl;
    exit(-1);
    break;
  }
  
  // Output txt file
  ofstream txtFile;
  txtFile.open(histoFiles_+".txt");
  
  

  // Analysis
  //objects to use
  unsigned int jetAlgo(event::AK5),metAlgo(event::CALOAK5), tauType(CALOTAU), leptonType(event::STDLEPTON);
  if(myKey=="TC") { jetAlgo=event::AK5JPT, metAlgo=event::TC, tauType = CALOTAU;}
  else if(myKey=="PF") { jetAlgo=event::AK5PF, metAlgo=event::PF, tauType = PFTAU;}
  else if(myKey.Contains("TaNC")) { jetAlgo=event::AK5PF, metAlgo=event::PF, tauType = PFTAU;}
  else if(myKey.Contains("HPS")) { jetAlgo=event::AK5PF, metAlgo=event::PF, tauType = HPSTAU;}
  //else if(myKey=="PFlow") { jetAlgo=event::AK5PFLOW, metAlgo=event::PFLOWMET, tauType = PFLOWTAU; leptonType = event::PFLOWLEPTON;}
  else if(myKey.Contains("PFlow")) { jetAlgo=event::AK5PFLOW, metAlgo=event::PFLOWMET, tauType = PFLOWTAU; leptonType = event::PFLOWLEPTON;}
  
  cout << "Events will be taken from "  << url << endl
       << "and analyzed with jet-algo="<< jetAlgo << " met-algo=" << metAlgo << endl;
  

  switch(qualifier_){
  case QCDMC:
    {
      //create the file ////////////////////////////////
      TString Filename_ = histoFiles_ +"_"+myKey+".root";
      TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
      //////////////////////////////////////////////////
      
      TH1D* wplusjets_pt          = new TH1D("wplusjets_pt"         ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_unknown  = new TH1D("wplusjets_pt_unknown" ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_quark    = new TH1D("wplusjets_pt_quark"   ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_bquark   = new TH1D("wplusjets_pt_bquark"  ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_gluon    = new TH1D("wplusjets_pt_gluon"   ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_allquark = new TH1D("wplusjets_pt_allquark",TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      
      TTree *tree = new TTree("tree", "tree");
      tree->Branch("JetPt", &JetPt, "JetPt/D");
      tree->Branch("AbsJetEta", &AbsJetEta, "AbsJetEta/D");
      tree->Branch("JetWidth", &JetWidth, "JetWidth/D");
      tree->Branch("__WEIGHT__", &__WEIGHT__, "__WEIGHT__/D");
      tree->Branch("__TARGET__", &__TARGET__, "__TARGET__/I");
      
      ///////////////////////
      
      // qcd multijet
      TString urlpath = inputArea_;
      //TString urlpath = "/lustre/data3/cmslocal/vischia/tau_dilepton/fakerate2011/inputSamples/";
      //double CrossBin[6] = {36750000000., 815900000., 53120000., 6359000., 784300., 115100.};  //cross section
      //      double CrossBin[5] = {815900000., 53120000., 6359000., 784300., 115100.};  //cross section // OLD
      double CrossBin[5] = {988287420., 66285328., 8148778., 1033680., 156293.3};  //cross section
      //Fall11
      ///  double CrossBin[5] = {816000000., 53100000., 6360000., 784000., 115000.};
      //for low stat files
      // find number of run events in file 
      double EventsBin[5] = {8213600., 1500000., 1500000., 1500000., 1500000.};
      TString DataSets[5] = {"qcd_Pt-15.root", "qcd_Pt-30.root", "qcd_Pt-50.root", "qcd_Pt-80.root", "qcd_Pt-120.root"};
      /*
      //for all stat files
      double EventsBin[5] = {8213600., 6529320., 4301392., 6407738., 6090400.};
      TString DataSets[5] = {"qcd_Pt-15.root", "qcd_Pt-30-B.root", "qcd_Pt-50-B.root", "qcd_Pt-80-B.root", "qcd_Pt-120-B.root"};
      */
      int nFiles = 5;
      //  double Lumi = 0.007653; //in pb-1 (for HLT_Jet30_v*)
      double Lumi = 0.00485234; //in pb-1 (effective lumi for HLT_Jet30_v*)
      //double ProcEvents[5] = {-1, 1000000., 1000000., 1000000., 1000000.}; //events to process, scale accordingly
      //double ProcEvents[5] = {-1, 1000000, 500000, 500000, 100000};
      double ProcEvents[5] = {-1, -1, -1, -1, -1};
      
      evReader = new event::Reader();
      for(int ifile = 1; ifile <= nFiles; ifile++){
	TString urlname = DataSets[ifile-1];
	TString url = urlpath+urlname;
	
	//open file with mini events
	TFile *f = TFile::Open(url);
	if(f==0) return ;
	if(f->IsZombie()) { f->Close(); return; }  
	TString tag("data");
	int nEntries = evReader->AssignMiniEventTreeFrom(f,tag,"myEventSelector"); 
	cout << nEntries << " events are available" << endl;
	if( nEntries == 0) { f->Close(); return; }
	
	double weight = Lumi*CrossBin[ifile-1]/EventsBin[ifile-1];
	//double weight = CrossBin[ifile-1]/double(nInitEvent);
	//weight = 1.;
	if(ProcEvents[ifile-1] != -1 && nEntries > ProcEvents[ifile-1]){
	  weight = (weight*nEntries/ProcEvents[ifile-1]);
	}
	
	int nTriggEvent = 0;
	int nProcEvents = 0;
	int nToProcess = nEntries;
	if(nToProcess > 10000000) nToProcess = 10000000;
	for(int i=0; i<nToProcess; ++i)
	  //    for(int i=0; i<nEntries; ++i)
	  {
	    nProcEvents++;
	    if(ProcEvents[ifile-1] != -1 && nProcEvents > ProcEvents[ifile-1])break;
	    //
	    // GET OBJECTS FOR NEW EVENT (nb. if kin results are available include them)
	    //
	    event::MiniEvent_t *ev = evReader->GetNewMiniEvent(i,tag);      
	    if(ev==0) continue;
	    
	    //trigger
	    TVectorD *trig = (TVectorD *)ev->triggerColl->At(0);
	    bool hasJettrig = ((*trig)[2]>0);
	    //if(hasJettrig)cout<<"hasJettrig"<<hasJettrig<<endl;
	    if(!hasJettrig) continue;
	    nTriggEvent++;  
	    
	    //get MC weight for pileup
	    double evtWeight = weight;
	    
	    TVectorD *classif = (TVectorD *)ev->eventClassif->At(0);
	    if(classif!=0){
	      if(classif->GetNoElements() >= 6){
		evtWeight *= (double)(*classif)[5];
	      }
	    }
	    
	    //Collect Physics Objects
	    std::vector<PhysicsObject> vertices                 = evReader->GetVertexCollection(ev);
	    std::vector<PhysicsObject> jets_without_arbitration = evReader->GetPhysicsObjectsFrom(ev,event::JET, jetAlgo  ); 
	    std::vector<PhysicsObject> jets = evReader->GetJetsWithArbitration( jets_without_arbitration, event::Reader::VTXCONSTRAINED );
	    std::vector<PhysicsObject> muons = evReader->GetPhysicsObjectsFrom(ev,event::MUON);
	    std::vector<PhysicsObject> electrons = evReader->GetPhysicsObjectsFrom(ev,event::ELECTRON);
	    std::vector<PhysicsObject> tausColl = evReader->GetPhysicsObjectsFrom(ev,event::TAU);
	    std::vector<PhysicsObject> taus; 
	    for(size_t iorigtau=0; iorigtau<tausColl.size(); iorigtau++){
	      if(tausColl[iorigtau][17] == tauType){
		taus.push_back(tausColl[iorigtau]);
	      }
	    }
	    
	    if(vertices.size()==0){ cout<<endl<<" Vertex was zero ???????"<<endl; continue; }
	    PhysicsObject & primaryVertex = vertices[0];
	    
	    JetCorrectionUncertainty * junc(0);
	    vector<double> jerFactors; jerFactors.clear();
	    vector<PhysicsObject> newJets; newJets.clear();
	    for(unsigned int i=0;i<jets.size();++i){ 
	      double corr_jer(1.);
	      /*if(!isData_)*/ newJets.push_back( smearedJet(jets[i], jets[i][34], 0/* 0=genpt, 1=random */, 0 /* 0=base, 1=jerup, 2=jerdown*/, corr_jer) );
	      jerFactors.push_back(1.);
	      //    std::cout << " jerfac: " << corr_jer << std::endl;
	    } 
	    /*if(!isData_)*/ jets=newJets;

	    // preselect objects /////////////////////////////////
	    DisableLtkCutOnJets(); Pt_Jet(MIN_JET_PT_CUT_); // select hard jets
	    //pt_Tau(MIN_TAU_PT_CUT_);  //select tau pT
	    vector<int> e_init, m_init, j_init, t_init;
	    PreSelectMuons( evReader,&m_init,muons, primaryVertex ); 
	    PreSelectElectrons(  evReader,&e_init,electrons, primaryVertex );
	    PreSelectJets(false, jerFactors, JES_, junc, jetAlgo,&j_init,jets);
	    PreSelectTaus( &t_init,taus, TAUPRONGS_, myKey, primaryVertex);
	    //////////////////////////////////////////////////////
	    
	    // only accept jets if dr > drmin in respect to electrons and muons /////////////////////
	    vector<int> emptyColl, j_toRemove; 
	    vector<int> j_afterLeptonRemoval;
	    ProcessCleaning(&j_init, &j_toRemove, &e_init, &emptyColl, jets, electrons, DRMIN_JET_E_ );
	    ProcessCleaning(&j_init, &j_toRemove, &m_init, &emptyColl, jets, muons,     DRMIN_JET_M_ );
	    ApplyCleaning(&j_init, &j_toRemove, &j_afterLeptonRemoval);
	    /////////////////////////////////////////////////////////////////////////////////////
	    
	    ///only accept taus if dr > drmin in respect to electrons and muons /////////////////////
	    vector<int> t_toRemove;
	    vector<int> t_final; 
	    ProcessCleaning(&t_init, &t_toRemove, &e_init, &emptyColl, taus, electrons, DRMIN_T_E_ );
	    ProcessCleaning(&t_init, &t_toRemove, &m_init, &emptyColl, taus, muons,     DRMIN_T_M_ );
	    ApplyCleaning(&t_init, &t_toRemove, &t_final);
	    /////////////////////////////////////////////////////////////////////////////////////////
	    
	    //remove jet matching to the trigger////////////////////////////////////////////////////// 
	    vector<int> j_final;
	    int nTrigJet = 0;
	    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	      int ind =  j_afterLeptonRemoval[ijet];
	      if(jets[ind][4] > 30)nTrigJet++;
	    }
	    //if(nTrigJet)cout<<"nTrig Jet "<<nTrigJet<<endl;
	    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	      int ind =  j_afterLeptonRemoval[ijet];
	      if(nTrigJet == 1 && jets[ind][4] > 30) continue; //not to use trigger matched Jet
	      j_final.push_back(ind);
	    }
	    //////////////////////////////////////////////////////////////////////////////////////////
	    
	    
	    //Fill Histograms
	    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
	      int ind_jet = j_final[ijet];
	      
	      

	      

	      if(isAntiBTag && (jets[ind_jet][33] > 0.679)) continue; //anti-btag CSVM
	      
	      double Rjet = ((jets[ind_jet][29]+jets[ind_jet][30]) > 0) ? sqrt(jets[ind_jet][29]+jets[ind_jet][30]) : 0.;
	      
	      int pgid(commondefinitions::PGID_UNKNOWN); if(!isDATA) pgid = TMath::Abs(jets[ind_jet][commondefinitions::JETPGID_OLD]); 	      
	      if( ! isDATA  ){       
		wplusjets_pt->Fill(jets[ind_jet].Pt(),evtWeight);
		// if we have at least one btag add also soft jet contribution to the collection of fakable taus   /////////////////////////////////////////////////////////
		
		if     ( pgid==commondefinitions::PGID_UNKNOWN )                                                               { wplusjets_pt_unknown ->Fill(jets[ind_jet].Pt(),evtWeight);}
		else if( pgid == commondefinitions::PGID_D || pgid == commondefinitions::PGID_U || pgid == commondefinitions::PGID_S || pgid == commondefinitions::PGID_C )             { wplusjets_pt_quark   ->Fill(jets[ind_jet].Pt(),evtWeight);}
		else if( pgid == commondefinitions::PGID_B     )                                                               { wplusjets_pt_bquark  ->Fill(jets[ind_jet].Pt(),evtWeight);}
		else if( pgid == commondefinitions::PGID_G     )                                                               { wplusjets_pt_gluon   ->Fill(jets[ind_jet].Pt(),evtWeight);}
		if( pgid == commondefinitions::PGID_D || pgid == commondefinitions::PGID_U || pgid == commondefinitions::PGID_S || pgid == commondefinitions::PGID_C || pgid == commondefinitions::PGID_B ){ wplusjets_pt_allquark->Fill(jets[ind_jet].Pt(),evtWeight);}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	      }
	      
	      

	      //fill for selected taus
	      double mindr = 0.5; int jtau = -1;
	      for(size_t itau = 0; itau < t_final.size(); itau++){
		int ind_tau = t_final[itau];
		double idr=taus[ind_tau].DeltaR(jets[ind_jet]);
		if(idr < mindr){
		  mindr = idr;
		  jtau = ind_tau;
		}
	      }
	      if(passing && mindr <  DRMIN_T_J_ ){
		JetPt = jets[ind_jet].Pt();
		AbsJetEta = fabs(jets[ind_jet].Eta());
		JetWidth = Rjet;
		__WEIGHT__ = evtWeight;
		__TARGET__ = 1;
		
		txtFile<< JetPt << " " << AbsJetEta << " " << JetWidth << " " << __WEIGHT__ << endl;
		tree->Fill();
	      }
	      else if(!passing && mindr >  DRMIN_T_J_ ){
		JetPt = jets[ind_jet].Pt();
		AbsJetEta = fabs(jets[ind_jet].Eta());
		JetWidth = Rjet;
		__WEIGHT__ = evtWeight;
		__TARGET__ = 0;
		txtFile<< JetPt << " " << AbsJetEta << " " << JetWidth << " " << __WEIGHT__ << endl;
		tree->Fill();
	      }
	    }
	    
	  } //end of each file
	cout<<" no. of triggered Events "<<nTriggEvent<<endl;
	txtFile.close();
	f->Close();
      } //end of all files
      outFile_->Write();
      outFile_->Close();
      break;
    }
  default:
    {
      //open file with mini events
      evReader = new event::Reader();
      TFile *f = TFile::Open(url);
      if(f==0) return ;
      if(f->IsZombie()) { f->Close(); return; }  
      TString tag("data");
      int nEntries = evReader->AssignMiniEventTreeFrom(f,tag,"myEventSelector"); 
      cout << nEntries << " events are available" << endl;
      if( nEntries == 0) { f->Close(); return; }
      
      //create the file ////////////////////////////////
      TString Filename_ = histoFiles_ +"_"+myKey+".root";
      TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
      //////////////////////////////////////////////////
      TH1D* wplusjets_pt          = new TH1D("wplusjets_pt"         ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_unknown  = new TH1D("wplusjets_pt_unknown" ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_quark    = new TH1D("wplusjets_pt_quark"   ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_bquark   = new TH1D("wplusjets_pt_bquark"  ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_gluon    = new TH1D("wplusjets_pt_gluon"   ,TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      TH1D* wplusjets_pt_allquark = new TH1D("wplusjets_pt_allquark",TString("j;jet p_{T} [GeV/c]; Entries"),25,0,250); 
      
      TTree *tree = new TTree("tree", "tree");
      tree->Branch("JetPt", &JetPt, "JetPt/D");
      tree->Branch("AbsJetEta", &AbsJetEta, "AbsJetEta/D");
      tree->Branch("JetWidth", &JetWidth, "JetWidth/D");
      tree->Branch("__WEIGHT__", &__WEIGHT__, "__WEIGHT__/D");
      tree->Branch("__TARGET__", &__TARGET__, "__TARGET__/I");
      
      ///////////////////////
      
      int nTriggEvent = 0, nSelEvents = 0;
      int nToProcess = nEntries;
      if(nToProcess > 2000000) nToProcess = 2000000; // FIXME: set this limit in another way
      for(int i=0; i<nToProcess; ++i)
	{
	  //
	  // GET OBJECTS FOR NEW EVENT (nb. if kin results are available include them)
	  //
	  if( ! (i% 10000 ) ){cout<<" -> Processing ... : "<<(i+1)<<endl;}
	  
	  event::MiniEvent_t *ev = evReader->GetNewMiniEvent(i,tag);      
	  if(ev==0) continue;
	  
	  //trigger
	  TVectorD *trig = (TVectorD *)ev->triggerColl->At(0);
	  bool hasMutrig = ((*trig)[0]>0); //bool hasEGtrig = ((*trig)[1]>0);
	  bool hasJettrig = ((*trig)[2]>0);
	  //if(!hasMutrig && !hasEGtrig) continue;
	  if(qualifier_ != QCDDATA && !hasMutrig) continue;
	  if(qualifier_ == QCDDATA && !hasJettrig) continue;
	  nTriggEvent++;
	  //////////////////////////////////////////////////////////
	  
	  double evtWeight = 1.0;
	  if(!isDATA){
	    TVectorD *classif = (TVectorD *)ev->eventClassif->At(0);
	    if(classif!=0){
	      /*
		if(classif->GetNoElements() >= 6){
		evtWeight *= (double)(*classif)[5];
		}
	      */
	      int npv = (*classif)[2]; //intime pu vertices // FIXME: this will be trueNVX in new nutples
	      float correctPUWeight_ = LumiWeights_.ITweight(npv);
	      evtWeight *= correctPUWeight_;
	    }
	  }
	  
	  //Collect Physics Objects
	  std::vector<PhysicsObject> vertices                 = evReader->GetVertexCollection(ev);
	  std::vector<PhysicsObject> jets_without_arbitration = evReader->GetPhysicsObjectsFrom(ev,event::JET, jetAlgo  ); 
	  std::vector<PhysicsObject> jets = evReader->GetJetsWithArbitration( jets_without_arbitration, event::Reader::VTXCONSTRAINED );
	  std::vector<PhysicsObject> muons = evReader->GetPhysicsObjectsFrom(ev,event::MUON, leptonType);
	  std::vector<PhysicsObject> electrons = evReader->GetPhysicsObjectsFrom(ev,event::ELECTRON);
	  std::vector<PhysicsObject> tausColl = evReader->GetPhysicsObjectsFrom(ev,event::TAU);
	  std::vector<PhysicsObject> taus; 
	  //      std::cout << "DEBUG: muons " << muons.size() << std::endl;
	  for(size_t iorigtau=0; iorigtau<tausColl.size(); iorigtau++){
	    if(tausColl[iorigtau][17] == tauType){
	      taus.push_back(tausColl[iorigtau]);
	    }
	  }
	  
	  if(vertices.size()==0){ cout<<endl<<" Vertex was zero ???????"<<endl; continue; }
	  PhysicsObject & primaryVertex = vertices[0];
	  
	  PhysicsObjectCollection mets   = evReader->GetPhysicsObjectsFrom(ev,event::MET,metAlgo);
	  if(qualifier_ != QCDDATA && mets.size()==0) { cout << "No met available"<<endl; continue;}
	  PhysicsObject met = mets[0];
	  ////////////////////////////////////////////////////////////////////////////////
	  
	  JetCorrectionUncertainty * junc(0);
	  vector<double> jerFactors; jerFactors.clear();
	  vector<PhysicsObject> newJets; newJets.clear();
	  for(unsigned int i=0;i<jets.size();++i){ 
	    double corr_jer(1.);
	    if(qualifier_==WMUMC)
	      newJets.push_back( smearedJet(jets[i], jets[i][34], 0/* 0=genpt, 1=random */, 0 /* 0=base, 1=jerup, 2=jerdown*/, corr_jer) );
	    jerFactors.push_back(1.);
	    //    std::cout << " jerfac: " << corr_jer << std::endl;
	  } 
	  if(qualifier_==WMUMC) jets=newJets;
	  // preselect objects /////////////////////////////////
	  DisableLtkCutOnJets(); Pt_Jet(MIN_JET_PT_CUT_); // select hard jets
	  //pt_Tau(MIN_TAU_PT_CUT_);  //select tau pT
	  vector<int> e_init, m_init, j_init, t_init;
	  PreSelectMuons( evReader,&m_init,muons, primaryVertex ); 
	  PreSelectElectrons(  evReader,&e_init,electrons, primaryVertex );
	  PreSelectJets(isDATA, jerFactors, JES_, junc, jetAlgo,&j_init,jets);
	  PreSelectTaus( &t_init,taus, TAUPRONGS_, myKey, primaryVertex);
	  //////////////////////////////////////////////////////
	  
	  if(qualifier_ != QCDDATA){
	    //Apply Lepton selection//////////////////////////////
	    //int nLepton = e_init.size() + m_init.size();
	    int nLepton = m_init.size();  //only muon +jet events
	    //      std::cout << "nLepton: " << nLepton << std::endl;
	    if(nLepton != 1)continue;
	    //std::cout << "DEBUG: 1lep" << std::endl;
	    if(e_init.size() > 0)continue;
	    //std::cout << "DEBUG: electron veto" << std::endl;
	    if( LooseMuonVeto( m_init[0], muons ))continue;
	    //std::cout << "DEBUG: loose mu veto" << std::endl;
	    if( LooseElectronVeto(evReader,-1,electrons))continue;
	    //std::cout << "DEBUG: loose e veto" << std::endl;
	    //Apply MET cut for W+jet selection
	    //if(met.Pt() < MET_CUT_) continue;
	    //Reconstruct W transverse mass
	    double leptonPt(0), deltaPhi(0);
	    double metPt = met.Pt();
	    if( m_init.size() == 1 ){ int m_i = m_init[0]; leptonPt = TMath::Abs(muons[m_i].Pt());     deltaPhi = muons[m_i].DeltaPhi( met );     }
	    //if( e_init.size() == 1 ){ int e_i = e_init[0]; leptonPt = TMath::Abs(electrons[e_i].Pt()); deltaPhi = electrons[e_i].DeltaPhi( met);  }
	    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
	    // Disabled MT cut
	    //if(mt < MET_CUT_) continue;
	    nSelEvents++;
	    ///////////////////////////////////////////////////////
	  }
	  
	  
	  // only accept jets if dr > drmin in respect to electrons and muons /////////////////////
	  vector<int> emptyColl, j_toRemove; 
	  vector<int> j_afterLeptonRemoval;
	  ProcessCleaning(&j_init, &j_toRemove, &e_init, &emptyColl, jets, electrons, DRMIN_JET_E_ );
	  ProcessCleaning(&j_init, &j_toRemove, &m_init, &emptyColl, jets, muons,     DRMIN_JET_M_ );
	  ApplyCleaning(&j_init, &j_toRemove, &j_afterLeptonRemoval);
	  /////////////////////////////////////////////////////////////////////////////////////////
	  
	  // only accept taus if dr > drmin in respect to electrons and muons /////////////////////
	  vector<int> t_toRemove;
	  vector<int> t_final; 
	  ProcessCleaning(&t_init, &t_toRemove, &e_init, &emptyColl, taus, electrons, DRMIN_T_E_ );
	  ProcessCleaning(&t_init, &t_toRemove, &m_init, &emptyColl, taus, muons,     DRMIN_T_M_ );
	  ApplyCleaning(&t_init, &t_toRemove, &t_final);
	  /////////////////////////////////////////////////////////////////////////////////////////
	  
	  vector<int> j_final;
	  if(qualifier_ == QCDDATA){
	    //remove jet matching to the trigger////////////////////////////////////////////////////// 
	    int nTrigJet = 0;
	    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	      int ind =  j_afterLeptonRemoval[ijet];
	      if(jets[ind][4] > 30)nTrigJet++;
	    }
	    //if(nTrigJet)cout<<"nTrig Jet "<<nTrigJet<<endl;
	    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	      int ind =  j_afterLeptonRemoval[ijet];
	      if(nTrigJet == 1 && jets[ind][4] > 30) continue; //not to use trigger matched Jet
	      j_final.push_back(ind);
	    }
	    //////////////////////////////////////////////////////////////////////////////////////////
	  }
	  else
	    j_final = j_afterLeptonRemoval;
	  
	  
	  //Fill Histograms
	  for(size_t ijet = 0; ijet < j_final.size(); ijet++){
	    int ind_jet = j_final[ijet];
	    
	    if(isAntiBTag && (jets[ind_jet][33] > 0.679)) continue; //anti-btag CSVM
	    
	    double jetPt = jets[ind_jet].Pt();
	    if(isDATA)jetPt = GetJetResidualPt(jets[ind_jet]);
	    double Rjet = ((jets[ind_jet][29]+jets[ind_jet][30]) > 0) ? sqrt(jets[ind_jet][29]+jets[ind_jet][30]) : 0.;
	    
	    
	    int pgid(commondefinitions::PGID_UNKNOWN); if(!isDATA) pgid = TMath::Abs(jets[ind_jet][commondefinitions::JETPGID_OLD]); 	      
	    if( ! isDATA  ){       
	      wplusjets_pt->Fill(jets[ind_jet].Pt(),evtWeight);
	      // if we have at least one btag add also soft jet contribution to the collection of fakable taus   /////////////////////////////////////////////////////////
	      if     ( pgid==commondefinitions::PGID_UNKNOWN )                                                               { wplusjets_pt_unknown ->Fill(jets[ind_jet].Pt(),evtWeight);}
	      else if( pgid == commondefinitions::PGID_D || pgid == commondefinitions::PGID_U || pgid == commondefinitions::PGID_S || pgid == commondefinitions::PGID_C )             { wplusjets_pt_quark   ->Fill(jets[ind_jet].Pt(),evtWeight);}
	      else if( pgid == commondefinitions::PGID_B     )                                                               { wplusjets_pt_bquark  ->Fill(jets[ind_jet].Pt(),evtWeight);}
	      else if( pgid == commondefinitions::PGID_G     )                                                               { wplusjets_pt_gluon   ->Fill(jets[ind_jet].Pt(),evtWeight);}
	      if( pgid == commondefinitions::PGID_D || pgid == commondefinitions::PGID_U || pgid == commondefinitions::PGID_S || pgid == commondefinitions::PGID_C || pgid == commondefinitions::PGID_B ){ wplusjets_pt_allquark->Fill(jets[ind_jet].Pt(),evtWeight);}
	      
	      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	      }
	      
	    
	    
	    
	    //fill for selected taus
	    double mindr = 0.5; int jtau = -1;
	    for(size_t itau = 0; itau < t_final.size(); itau++){
	      int ind_tau = t_final[itau];
	      double idr=taus[ind_tau].DeltaR(jets[ind_jet]);
	      if(idr < mindr){
		mindr = idr;
		jtau = ind_tau;
	      }
	    }
	    if(passing && mindr < DRMIN_T_J_){
	      JetPt = jetPt;
	      AbsJetEta = fabs(jets[ind_jet].Eta());
	      JetWidth = Rjet;
	      __WEIGHT__ = evtWeight;
	      __TARGET__ = 1;
	      
	      txtFile<< JetPt << " " << AbsJetEta << " " << JetWidth << " " << __WEIGHT__ << endl;
	      tree->Fill();
	    }
	    else if(!passing && mindr > DRMIN_T_J_){
	      JetPt = jetPt;
	      AbsJetEta = fabs(jets[ind_jet].Eta());
	      JetWidth = Rjet;
	      __WEIGHT__ = evtWeight;
	      __TARGET__ = 0;
	      
	      txtFile<< JetPt << " " << AbsJetEta << " " << JetWidth << " " << __WEIGHT__ << endl;
	      tree->Fill();
	    }
	    
	  }
	  
	}
      cout<<" no. of triggered Events "<<nTriggEvent<<endl;
      if(qualifier_ != QCDDATA) cout<<" no. of selected Events "<<nSelEvents<<endl;
      
      txtFile.close();
      f->Close();
      outFile_->Write();
      outFile_->Close();
      break;
    }
  }
}

void TauFakesHelper::Trainer(unsigned int qualifier)
{
  qualifier_ = qualifier;
  SetParameters();
  
  
  // Obtain signal and background training trees;
  
  fPassing_ = new TFile(trainingTreesArea_+TString("/TauFR")+infix_+TString("Passing_PFlow.root")); // FIXME: move myKey to global variable
  if(fPassing_==0) return ;
  if(fPassing_->IsZombie()) { fPassing_->Close(); return; }
  tPassing_ = (TTree *) fPassing_->Get("tree");
  
  fFailing_ = new TFile(trainingTreesArea_+TString("/TauFR")+infix_+TString("Failing_PFlow.root")); // FIXME: move myKey to global variable
  if(fFailing_==0) return ;
  if(fFailing_->IsZombie()) { fFailing_->Close(); return; }
  tFailing_ = (TTree *) fFailing_->Get("tree");
  
  std::cout << "Training with " << tPassing_->GetEntries()
	    << " signal events." <<  std::endl;
  std::cout << "Training with " << tFailing_->GetEntries()
	    << " background events." << std::endl;
  
  ///  // Note: one tree argument -> tree has to contain a branch __TARGET__
  ///  //       two tree arguments -> signal and background tree
  ///  
  ///  PhysicsTools::TreeTrainer trainer(tPassing_); // FIXME: set monitoring
  ///  trainer.addTree(tFailing_);
  ///  
  ///  Calibration::MVAComputer *calib = trainer.train(configArea_+TString("fakerateMVADef")+infix_+TString(".xml"));
  ///  
  ///  MVAComputer::writeCalibration(outputArea_+TString("TrainedWMu")+infix_+TString("PFlow_pt20.mva"), calib);
  ///  
  ///  cout << "Trained.mva written." << endl;
  ///  
  ///  //  delete calib;
  
  MVATrainer trainer(string( configArea_+TString("fakerateMVADef")+infix_+TString(".xml") ) );
  trainer.setMonitoring(true);// ROOT file with histograms // FIXME: investigate why it seems not to be working
  
  // obtain signal and background training trees;
  TreeReader sigReader(tPassing_);
  TreeReader bkgReader(tFailing_);
  
  // __TARGET__ branch is already into trees
  //  bool sigTarget = true;
  //  sigReader.addSingle("__TARGET__", &sigTarget);
  //
  //  bool bkgTarget = false;
  //  bkgReader.addSingle("__TARGET__", &bkgTarget);
  
  
  // looping over dataset until trainer is satisfied
  // Note: In case of ROOT tree, TreeTrainer can do it for you
  //       (see testFWLiteTreeTrainer.C)
  
  for(;;) {
    Calibration::MVAComputer *calib = 
      trainer.getTrainCalibration();
    
    if (!calib)
      break;
    
    MVAComputer computer(calib, true);
    
    sigReader.loop(&computer);
    bkgReader.loop(&computer);
    
    trainer.doneTraining(calib);
  }
  
  // Get calibration after training
  Calibration::MVAComputer *calib = trainer.getCalibration();
  
  // Save calibration to .mva file
  MVAComputer::writeCalibration(trainingOutputArea_+TString("Trained_")+infix_+TString("_PFlow.mva"), calib);


  std::cout << trainingOutputArea_+TString("Trained_")+infix_+TString("_PFlow.mva") << " written" << std::endl;

  delete calib;
  
}

void TauFakesHelper::Monitor(unsigned int qualifier)
{
  qualifier_ = qualifier;
  SetParameters();
 
  TString calibrationFile = trainingOutputArea_+"/Trained_"+infix_+"_PFlow.mva";

  // FIXME: tdrstyle
  // SetTDRStyle();

  TString inputFileName = trainingOutputArea_+TString("/train_")+infix_+TString("/train_FakeRateMethod_output.root");
  TFile* file = TFile::Open(inputFileName);
  TString testTreeName = "TestTree";
  TTree* testTree = (TTree*)file->Get(testTreeName);
  Float_t JetPt, JetWidth, AbsJetEta, weight, KNN;
  Int_t classID; 

  testTree->SetBranchAddress("JetPt", &JetPt);
  testTree->SetBranchAddress("JetWidth", &JetWidth);
  testTree->SetBranchAddress("AbsJetEta", &AbsJetEta);
  testTree->SetBranchAddress("weight", &weight);
  testTree->SetBranchAddress("KNN", &KNN);
  testTree->SetBranchAddress("classID", &classID);

  Int_t nentries = (Int_t)testTree->GetEntries();
  
  double ptBin[16] = {0, 10, 20, 22.5, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 250};
  TH1D *ptJetAll = new TH1D("ptJetAll", " ", 15, ptBin);
  TH1D *ptJetPassMVA = new TH1D("ptJetPassMVA", " ", 15, ptBin);
  TH1D *ptJetPass = new TH1D("ptJetPass", " ", 15, ptBin);
  TH1D *ptJetFR = new TH1D("ptJetFR", " ", 15, ptBin);
  TH1D *ptJetFRMVA = new TH1D("ptJetFRMVA", " ", 15, ptBin);
  TH1D *absEtaJetAll = new TH1D("absEtaJetAll", " ", 6, 0, 3.0);
  TH1D *absEtaJetPassMVA = new TH1D("absEtaJetPassMVA", " ", 6, 0, 3.0);
  TH1D *absEtaJetPass = new TH1D("absEtaJetPass", " ", 6, 0, 3.0);
  TH1D *absEtaJetFR = new TH1D("absEtaJetFR", " ", 6, 0, 3.0);
  TH1D *absEtaJetFRMVA = new TH1D("absEtaJetFRMVA", " ", 6, 0, 3.0);
  TH1D *rJetAll = new TH1D("rJetAll", " ", 12, 0, 0.48);
  TH1D *rJetPassMVA = new TH1D("rJetPassMVA", " ", 12, 0, 0.48);
  TH1D *rJetPass = new TH1D("rJetPass", " ", 12, 0, 0.48);
  TH1D *rJetFR = new TH1D("rJetFR", " ", 12, 0, 0.48);
  TH1D *rJetFRMVA = new TH1D("rJetFRMVA", " ", 12, 0, 0.48);

  for(Int_t ev = 0; ev < nentries; ev++){
    testTree->GetEntry(ev);
    double frWeight = 0;
    if(classID == 0)frWeight = 1.0;
    ptJetAll->Fill(JetPt, weight);
    ptJetPassMVA->Fill(JetPt, KNN*weight);
    ptJetPass->Fill(JetPt, frWeight*weight);
    absEtaJetAll->Fill(AbsJetEta, weight);
    absEtaJetPassMVA->Fill(AbsJetEta, KNN*weight);
    absEtaJetPass->Fill(AbsJetEta, frWeight*weight);
    rJetAll->Fill(JetWidth, weight);
    rJetPassMVA->Fill(JetWidth, KNN*weight);
    rJetPass->Fill(JetWidth, frWeight*weight);
  }
  AddFakeRate(ptJetPass, ptJetAll, ptJetFR);
  AddFakeRate(absEtaJetPass, absEtaJetAll, absEtaJetFR);
  AddFakeRate(rJetPass, rJetAll, rJetFR); 
  AddFakeRate(ptJetPassMVA, ptJetAll, ptJetFRMVA);
  AddFakeRate(absEtaJetPassMVA, absEtaJetAll, absEtaJetFRMVA);
  AddFakeRate(rJetPassMVA, rJetAll, rJetFRMVA);

  //Draw histograms
  GetReweightedPlots(ptJetPass, ptJetPassMVA, TString("ptDistribution_") + infix_,           "Jet p_{T} [GeV]");
  GetReweightedPlots(absEtaJetPass, absEtaJetPassMVA, TString("etaDistribution_") + infix_,  "Jet |#eta|");
  GetReweightedPlots(rJetPass, rJetPassMVA, TString("rDistribution_") + infix_,               "Jet Width");
  
  GetPlotsFRPt(ptJetFRMVA, ptJetFR, TString("pt_fr_diff_mc_") + infix_ );

  GetPlotsFRAbsEta(absEtaJetFRMVA, absEtaJetFR, TString("eta_fr_diff_mc_") + infix_ );

  GetPlotsFRR(rJetFRMVA, rJetFR, TString("r_fr_diff_mc_") + infix_ );
  
  
}

void TauFakesHelper::AddFakeRate(TH1* hn, TH1* hd, TH1* hf){
  TH1D * n; TH1 *d;
  n = (TH1D*)hn->Clone();
  d = (TH1D*)hd->Clone();

  int nBinsX1 = n->GetNbinsX(); int nBinsX2 = d->GetNbinsX();

  if( nBinsX1 != nBinsX2 ){
    //cout<<endl<<endl<<"H1 N bins is = "<<nBinsX1<<" but H2 N bin is = "<<nBinsX2<<" ... quit! "<<endl;
    return;
  }

  for(int i=1; i<= nBinsX1;i++){
    double numContent = n->GetBinContent(i); 
    double denContent = d->GetBinContent(i);
    double err = GetErrorFraction(numContent,denContent);
    //cout<<endl<<" i "<<i<<" num = "<<numContent<<" den = "<<denContent<<endl;
    
    if(denContent){ hf->SetBinContent(i, numContent/denContent); hf->SetBinError(i,err); } 
    
  }

  delete n; delete d;

}


double TauFakesHelper::GetErrorFraction( double& a,double& b){
  double ret(0);
  if (b){
    double temp = fabs(a)/(b*b)+ fabs( (a*a)/(b*b*b)) ;
    ret = sqrt(temp);
  }
  return ret;
}

void TauFakesHelper::GetReweightedPlots(TH1D* trueDist, TH1D* mvaDist, TString name, TString xaxis_title)
{
  TCanvas* c1 = new TCanvas();
  c1->SetLogy();
  c1->cd();
  
  mvaDist->SetMarkerStyle(20);
  mvaDist->SetMarkerSize(1.0);
  mvaDist->SetMarkerColor(2);
  mvaDist->SetLineColor(2);
  trueDist->SetMarkerStyle(21);
  trueDist->SetMarkerSize(1.0);
  trueDist->SetMarkerColor(4);
  trueDist->SetLineColor(4);
  
  //  mvaDist->GetYaxis()->SetRangeUser(ymin, ymax);
  //  mvaDist->GetXaxis()->SetRangeUser(xmin, xmax);
  mvaDist->GetXaxis()->SetTitle(xaxis_title);
  mvaDist->GetYaxis()->SetTitleOffset(1.4);
  mvaDist->GetYaxis()->SetTitle("Events/bin");

  setXerrorBars();
  mvaDist->Draw();
  trueDist->Draw("same");
  
  TLegend *pl = new TLegend(0.6,0.80,0.8,0.9);
  pl->SetTextSize(0.04);
  pl->SetFillColor(0);
  pl->SetTextFont(132);
  pl->SetBorderSize(0);
  TLegendEntry *ple;
  ple = pl->AddEntry(mvaDist, "kNN output", "lp");
  ple = pl->AddEntry(trueDist, "True", "lp");
  ple->SetMarkerSize(0.1);
  pl->Draw();
  
  c1->Update();
  //  gPad->SetTopMargin(0);
  //  gPad->SetBottomMargin(0.3);
  //  gPad->SetGridy();
  //  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
  
  //  h1->GetYaxis()->SetTitleSize(0.1);
  //  h1->GetXaxis()->SetTitleSize(0.1);
  //  h1->GetXaxis()->SetLabelSize(0.1);
  //  h1->GetYaxis()->SetLabelSize(0.06);
  
  c1->Update();
  
  TString filename = outputArea_+"/plotsMVA/"+name;
  c1->SaveAs(filename+".png");
  c1->SaveAs(filename+".pdf");
  c1->SaveAs(filename+".root");
  
  c1->Close();
}

void TauFakesHelper::GetPlotsFRPt(TH1D* data, TH1D* mc, TString name, const char* xaxis_title, float ymin, float ymax, float xmin, float xmax)
{
  TCanvas* c1 = new TCanvas();
  //c1->SetLogy();
  //Pad positioning stuff
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.3,0.3,1.0};
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
  gPad->SetLogy();

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.0);
  data->SetMarkerColor(2);
  data->SetLineColor(2);
  mc->SetMarkerStyle(21);
  mc->SetMarkerSize(1.0);
  mc->SetMarkerColor(4);
  mc->SetLineColor(4);

  data->GetYaxis()->SetRangeUser(ymin, ymax);
  data->GetXaxis()->SetRangeUser(xmin, xmax);
  data->GetXaxis()->SetTitle(xaxis_title);
  data->GetYaxis()->SetTitleOffset(1.4);
  data->GetYaxis()->SetTitle("Fake Rate");

  data->Draw();
  mc->Draw("same");
  
  TLegend *pl = new TLegend(0.6,0.80,0.8,0.9);
  pl->SetTextSize(0.04);
  pl->SetFillColor(0);
  pl->SetTextFont(132);
  pl->SetBorderSize(0);
  TLegendEntry *ple;
  ple = pl->AddEntry(data, "kNN output", "lp");
  ple = pl->AddEntry(mc, "True", "lp");
  ple->SetMarkerSize(0.1);
  pl->Draw();
  
  c1->Update();
  c1->cd(2);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridy();
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);

  TH1D *h1 = new TH1D(*data);
  h1->Add(mc, -1.0);
  h1->Divide(mc);
  h1->SetMarkerStyle(24);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h1->GetYaxis()->SetRangeUser(-1, 1);
  h1->GetXaxis()->SetRangeUser(xmin, xmax);
  h1->GetXaxis()->SetTitle(xaxis_title);
  h1->GetYaxis()->SetTitleOffset(0.5);
  h1->GetYaxis()->SetTitle("Ratio");
  h1->GetYaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetLabelSize(0.1);
  h1->GetYaxis()->SetLabelSize(0.06);
  
  h1->Draw();
  c1->Update();
  
  TString filename = outputArea_+"/plotsMVA/"+name;
  c1->SaveAs(filename+".png");
  c1->SaveAs(filename+".pdf");
  c1->SaveAs(filename+".root");

  c1->Close();
}

void TauFakesHelper::GetPlotsFRAbsEta(TH1D* data, TH1D* mc, TString name, const char* xaxis_title, float ymin, float ymax, float xmin, float xmax)
{
  TCanvas* c1 = new TCanvas();
  //c1->SetLogy();
  //Pad positioning stuff
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.3,0.3,1.0};
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
  gPad->SetLogy();

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.0);
  data->SetMarkerColor(2);
  data->SetLineColor(2);
  mc->SetMarkerStyle(21);
  mc->SetMarkerSize(1.0);
  mc->SetMarkerColor(4);
  mc->SetLineColor(4);

  data->GetYaxis()->SetRangeUser(ymin, ymax);
  data->GetXaxis()->SetRangeUser(xmin, xmax);
  data->GetXaxis()->SetTitle(xaxis_title);
  data->GetYaxis()->SetTitleOffset(1.4);
  data->GetYaxis()->SetTitle("Fake Rate");

  data->Draw();
  mc->Draw("same");
  
  TLegend *pl = new TLegend(0.6,0.80,0.8,0.9);
  pl->SetTextSize(0.04);
  pl->SetFillColor(0);
  pl->SetTextFont(132);
  pl->SetBorderSize(0);
  TLegendEntry *ple;
  ple = pl->AddEntry(data, "kNN output", "lp");
  ple = pl->AddEntry(mc, "True", "lp");
  ple->SetMarkerSize(0.1);
  pl->Draw();
  
  c1->Update();
  c1->cd(2);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridy();
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);

  TH1D *h1 = new TH1D(*data);
  h1->Add(mc, -1.0);
  h1->Divide(mc);
  h1->SetMarkerStyle(24);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h1->GetYaxis()->SetRangeUser(-1, 1);
  h1->GetXaxis()->SetRangeUser(xmin, xmax);
  h1->GetXaxis()->SetTitle(xaxis_title);
  h1->GetYaxis()->SetTitleOffset(0.5);
  h1->GetYaxis()->SetTitle("Ratio");
  h1->GetYaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetLabelSize(0.1);
  h1->GetYaxis()->SetLabelSize(0.06);
  
  h1->Draw();
  c1->Update();
  
  TString filename = outputArea_+"/plotsMVA/"+name;
  c1->SaveAs(filename+".png");
  c1->SaveAs(filename+".pdf");
  c1->SaveAs(filename+".root");

  c1->Close();
}


void TauFakesHelper::GetPlotsFRR(TH1D* data, TH1D* mc, TString name,  const char* xaxis_title, float ymin, float ymax, float xmin, float xmax)
{
  TCanvas* c1 = new TCanvas();
  //c1->SetLogy();
  //Pad positioning stuff
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.3,0.3,1.0};
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
  gPad->SetLogy();

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.0);
  data->SetMarkerColor(2);
  data->SetLineColor(2);
  mc->SetMarkerStyle(21);
  mc->SetMarkerSize(1.0);
  mc->SetMarkerColor(4);
  mc->SetLineColor(4);

  data->GetYaxis()->SetRangeUser(ymin, ymax);
  data->GetXaxis()->SetRangeUser(xmin, xmax);
  data->GetXaxis()->SetTitle(xaxis_title);
  data->GetYaxis()->SetTitleOffset(1.4);
  data->GetYaxis()->SetTitle("Fake Rate");

  data->Draw();
  mc->Draw("same");
  
  TLegend *pl = new TLegend(0.6,0.80,0.8,0.9);
  pl->SetTextSize(0.04);
  pl->SetFillColor(0);
  pl->SetTextFont(132);
  pl->SetBorderSize(0);
  TLegendEntry *ple;
  ple = pl->AddEntry(data, "kNN output", "lp");
  ple = pl->AddEntry(mc, "True", "lp");
  ple->SetMarkerSize(0.1);
  pl->Draw();
  
  c1->Update();
  c1->cd(2);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridy();
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);

  TH1D *h1 = new TH1D(*data);
  h1->Add(mc, -1.0);
  h1->Divide(mc);
  h1->SetMarkerStyle(24);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h1->GetYaxis()->SetRangeUser(-1, 1);
  h1->GetXaxis()->SetRangeUser(xmin, xmax);
  h1->GetXaxis()->SetTitle(xaxis_title);
  h1->GetYaxis()->SetTitleOffset(0.5);
  h1->GetYaxis()->SetTitle("Ratio");
  h1->GetYaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetTitleSize(0.1);
  h1->GetXaxis()->SetLabelSize(0.1);
  h1->GetYaxis()->SetLabelSize(0.06);
  
  h1->Draw();
  c1->Update();

  TString filename = outputArea_+"/plotsMVA/"+name;
  c1->SaveAs(filename+".png");
  c1->SaveAs(filename+".pdf");
  c1->SaveAs(filename+".root");
  
  c1->Close();
}




void TauFakesHelper::PrepareFiles(unsigned int qualifier){
  qualifier_ = qualifier;
  SetParameters();

  using namespace PhysicsTools;

  //  TFile *f = TFile::Open("/lustre/lip.pt/data/cmslocal/anayak/TMVA_CMSSW_4_2_4/FakeRateFiles/TauFRMCFailing_pt20_30To170_PFlow.root");
  //  TFile *f = TFile::Open("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_8_patch4b/mc/processed-muon-Summer11-oldTau/out-data.root");

  ///  TFile *f = TFile::Open("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_8_patch4b/data/processed-muon-RunA11-oldTau/out-data.root");
  ///  TFile *f = TFile::Open("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_8_patch4b/mc/processed-muon-Summer11-oldTau/out-wjets.root");

  //------------

//  TString inputFilePath = "./QCDDiJet/train/DiJetData_Pt20/";
//  MVAComputer mva("tempTraiedFiles/TrainedDiJetDataPFlow_pt20_new.mva");  
//  TFile *f = TFile::Open(inputFilePath.Append("train_FakeRateMethod_output.root"));
//  TFile *outFile_ = TFile::Open( "histosDiJetQCDJetData.root", "RECREATE" );

///////  TString inputFilePath = "./QCDDiJet/train/DiJetMC_Pt20/";
///////  MVAComputer mva("tempTraiedFiles/TrainedDiJetMCPFlow_pt20_new.mva");
///////  TFile *f = TFile::Open(inputFilePath.Append("train_FakeRateMethod_output.root"));
///////  TFile *outFile_ = TFile::Open( "histosDiJetQCDJetMC.root", "RECREATE" );


  TString inputFileName = trainingOutputArea_ + TString("/train_") + infix_ + TString("/train_FakeRateMethod_output.root");
  std::cout << "Opening input filename " << inputFileName << std::endl;
  TFile* f = TFile::Open(inputFileName);
  std::cout << "Writing to filename " << inputFileName << std::endl;
  TFile* outFile_ = TFile::Open(outputArea_+TString("histos")+infix_+TString(".root"), "RECREATE");
  
//  TString inputFilePath = "./WMuJet/train/WMuData_pt20/";
//  //  MVAComputer mva("/lustre/data3/cmslocal/vischia/MyTMVA_temporaneo/tempTraiedFiles/TrainedWMuDataPFlow_pt20_againstInvalid2.mva");
//  MVAComputer mva("./TrainedWMuDataPFlow_pt20.mva");
//  TFile *f = TFile::Open(inputFilePath.Append("train_FakeRateMethod_output.root"));
//  TFile *outFile_ = TFile::Open( "histosDiJetWMuData.root", "RECREATE" );

//  TString inputFilePath = "./WMuJet/train/WMuMC_pt20/";
//  //  MVAComputer mva("tempTraiedFiles/TrainedWMuMCPFlow_pt20.mva");
//  //  MVAComputer mva("/lustre/data3/cmslocal/vischia/MyTMVA_temporaneo/tempTraiedFiles/TrainedWMuMCPFlow_pt20_againstInvalid.mva");
//  MVAComputer mva("./TrainedWMuMCPFlow_pt20.mva");
//  TFile *f = TFile::Open(inputFilePath.Append("train_FakeRateMethod_output.root"));
//  TFile *outFile_ = TFile::Open( "histosDiJetWMuMC.root", "RECREATE" );



  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  ///  TTree *tree = (TTree *) f->Get("m_tau_DataDriven");
  ///  TTree *tree = (TTree *) f->Get("tree");
  TTree *tree = (TTree *) f->Get("TestTree");
  Float_t JetPt, JetWidth, AbsJetEta, weight, KNN;
  Int_t classID;

  Float_t RcT, LepPt, MetPt, LepMetDPhi, JetMulti, JetBMulti, TauMetDPhi;
 
  tree->SetBranchAddress("JetPt", &JetPt);
  tree->SetBranchAddress("JetWidth", &JetWidth);
  tree->SetBranchAddress("AbsJetEta", &AbsJetEta);
  //  tree->SetBranchAddress("__WEIGHT__", &Weight);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("KNN", &KNN);

  
  /* TODO: rc_t as spectator variable
      tree->SetBranchAddress("pt_t", &JetPt);
      tree->SetBranchAddress("radius_t", &JetWidth);
      tree->SetBranchAddress("eta_t", &AbsJetEta);
      tree->SetBranchAddress("weight", &Weight);
      tree->SetBranchAddress("is_os", &IsOS);
      tree->SetBranchAddress("rc_t", &RcT);
      
      
      //    tree->SetBranchAddress(nJetPt;
      //    tree->SetBranchAddress("pt_t");               "pt_t", &nJetPt);    
      tree->SetBranchAddress("pt_l",               &LepPt);
      tree->SetBranchAddress("pt_met",             &MetPt);
      tree->SetBranchAddress("Dphi_lepton_met",    &LepMetDPhi);
      tree->SetBranchAddress("multiplicity_j",     &JetMulti);
      tree->SetBranchAddress("btagmultiplicity_j", &JetBMulti);
      tree->SetBranchAddress("Dphi_tau_met",       &TauMetDPhi);
      
  */



///  tree->SetBranchAddress("pt_t", &JetPt);
///  tree->SetBranchAddress("radius_t", &JetWidth);
///  tree->SetBranchAddress("eta_t", &AbsJetEta);
///  tree->SetBranchAddress("weight", &Weight);

  Int_t nentries = (Int_t)tree->GetEntries();
  cout<<"entries "<<nentries<<endl;
  //TreeReader reader(tree);
  

  double ptBin[16] = {0, 10, 20, 22.5, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 250};
  //20, 0, 100.

  TH1D *ptJetAll = new TH1D("ptJetAll", " ", 15, ptBin);
  TH1D *ptJetPass = new TH1D("ptJetPass", " ", 15, ptBin);
  TH1D *ptJetFR = new TH1D("ptJetFR", " ", 15, ptBin);
  TH1D *absEtaJetAll = new TH1D("absEtaJetAll", " ", 6, 0, 3.0);
  TH1D *absEtaJetPass = new TH1D("absEtaJetPass", " ", 6, 0, 3.0);
  TH1D *absEtaJetFR = new TH1D("absEtaJetFR", " ", 6, 0, 3.0);
  TH1D *rJetAll = new TH1D("rJetAll", " ", 12, 0, 0.48);
  TH1D *rJetPass = new TH1D("rJetPass", " ", 12, 0, 0.48);
  TH1D *rJetFR = new TH1D("rJetFR", " ", 12, 0, 0.48);

  double sumEvent = 0;
  //int nev = (nentries < 2000000) ? nentries : 2000000;
  int nev = nentries;
  cout<<"events to process "<<nev<<endl;

  for(Int_t ev = 0; ev < nev; ev++){
    tree->GetEntry(ev);
    ///    Variable::ValueList vars;
    ///    vars.add("JetPt", JetPt);
    ///    vars.add("AbsJetEta", AbsJetEta);
    ///    vars.add("JetWidth", JetWidth);
    ///    
    ///    double value = mva.eval(vars);
    
    //    sumEvent += value;
    sumEvent +=KNN;
    //cout << value << endl;
    ptJetAll->Fill(JetPt, weight);
    ptJetPass->Fill(JetPt, KNN*weight);
    absEtaJetAll->Fill(AbsJetEta, weight);
    absEtaJetPass->Fill(AbsJetEta, KNN*weight);
    rJetAll->Fill(JetWidth, weight);
    rJetPass->Fill(JetWidth, KNN*weight);
  }

  cout<<"events pass "<<sumEvent<<endl;

  AddFakeRate(ptJetPass, ptJetAll, ptJetFR);
  AddFakeRate(absEtaJetPass, absEtaJetAll, absEtaJetFR);
  AddFakeRate(rJetPass, rJetAll, rJetFR);
  
  
  f->Close();
  outFile_->Write();
  outFile_->Close();

  // note that there is also a TTree interface which can read entries
  // directly from a ROOT tree. This interface is found in the
  // class "TreeReader" in this package.
  //
  // Note that the TreeReader can also be used for ROOT-like
  // interfacing of the MVAComputer::eval method, like the
  // TTree::Fill() method. This might come in handy and might be
  // simpler to use than the ValueList interface.
  //
  // See "MVATrainer/test/testFWLiteEvaluation.C" for an example.
  
  
  //Variable::ValueList vars;
  //vars.add("JetPt", 50.0);
  //vars.add("AbsJetEta", 2.0);
  //vars.add("JetWidth", 0.2);
  //
  //cout << mva.eval(vars) << endl;
  std::cout << "Files prepared" << std::endl;
}



void TauFakesHelper::ComputeTauFake(string type, vector<double>& finalValues, double TauPt, int option)
{
  string kNNfileData_, kNNfileMC_; 
  string frFileNameData_, frFileNameMC_;
  string frHistoName_ = "ptJetFR";
  
  if(type == "DiJet"){
    kNNfileData_ = trainingOutputArea_+"/Trained_QCDData_PFlow.mva";
    //kNNfileMC_   = trainingOutputArea_+"/Trained_QCDMC_PFlow.mva";
    kNNfileMC_   = trainingOutputArea_+"/Trained_QCDData_PFlow.mva";
    frFileNameData_ = outputArea_+"/histosQCDData.root";
    frFileNameMC_   = outputArea_+"/histosQCDMC.root";
  }
  if(type == "WMu"){
    kNNfileData_ = trainingOutputArea_+"/Trained_WMuData_PFlow.mva";
    kNNfileMC_   = trainingOutputArea_+"/Trained_WMuMC_PFlow.mva";
    //kNNfileMC_   = trainingOutputArea_+"/Trained_WMuData_PFlow.mva";
    frFileNameData_ = outputArea_+"/histosWMuData.root";
    frFileNameMC_   = outputArea_+"/histosWMuMC.root";
  }


  string DataMVA_(kNNfileData_);
  string mcMVA_(kNNfileMC_);
  
  
  TFile *fr_data_file = new TFile(frFileNameData_.c_str());
  TH1F* fr_data = (TH1F*)fr_data_file->Get(frHistoName_.c_str())->Clone("fr_data");
  
  TFile *fr_mc_file = new TFile(frFileNameMC_.c_str());
  TH1F* fr_mc = (TH1F*)fr_mc_file->Get(frHistoName_.c_str())->Clone("fr_mc");

  
  cout<<DataMVA_<<" "<<mcMVA_<<endl;
  
  vector<TString> listOfurls_; listOfurls_.clear();
  TString dataFolder, mcFolder;
  //mu-tau
  //dataFolder  = "/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/data/processed-v2/v2-data-20GeV-noweights/";
  //mcFolder    = "/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/processed-v2/v2-mc-20GeV-noweights/";
  
  cout << "Ntuples area: " << ntuplesArea_ << endl;

  dataFolder = ntuplesArea_+"/nomt-2012-V1-data-MU-20GeV/";
  mcFolder   = ntuplesArea_+"/nomt-2012-V1-mc-MU-20GeV/";
  
  listOfurls_.push_back(dataFolder+TString("out-data.root"));
  listOfurls_.push_back(mcFolder+TString("out-ttbar.root"));
   
  listOfurls_.push_back(mcFolder+TString("out-wjets.root"));  
  listOfurls_.push_back(mcFolder+TString("out-zjets.root"));
  listOfurls_.push_back(mcFolder+TString("out-singletop.root"));
  listOfurls_.push_back(mcFolder+TString("out-ww.root"));
  listOfurls_.push_back(mcFolder+TString("out-wz.root"));
  listOfurls_.push_back(mcFolder+TString("out-zz.root"));
  listOfurls_.push_back(mcFolder+TString("out-qcd.root"));

///  listOfurls_.push_back(dataFolder+TString("out-data.root"));
///  //  listOfurls_.push_back(mcFolder+TString("out-ttbar.root"));
///  listOfurls_.push_back(mcFolder+TString("out-ttbar_mutau.root"));
///
///  //  listOfurls_.push_back(mcFolder+TString("out-ddbkg.root"));
///  //  listOfurls_.push_back(mcFolder+TString("out-mcbkg.root"));
///
///  listOfurls_.push_back(mcFolder+TString("out-wjets.root"));  
///  listOfurls_.push_back(mcFolder+TString("out-zjets.root"));
///  listOfurls_.push_back(mcFolder+TString("out-singletop.root"));
///  listOfurls_.push_back(mcFolder+TString("out-ww.root"));
///  listOfurls_.push_back(mcFolder+TString("out-wz.root"));
///
///  //listOfurls_.push_back(mcFolder+TString("out-zz.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-qcd_30to80_BCtoE.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-qcd_80to170_BCtoE.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-qcd_EM_Pt30to80.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-qcd_EM_Pt80to170.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-qcdmu10_20toinf.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-PhotonJets_30to50.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-PhotonJets_50to80.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-PhotonJets_80to120.root"));
///  //listOfurls_.push_back(mcFolder+TString("out-PhotonJets_120to170.root"));
  
  /*
  vector<double>listOfWeights_; listOfWeights_.clear();
  listOfWeights_.push_back(1); //data
  listOfWeights_.push_back(0.02769897);  //ttbar
  listOfWeights_.push_back(0.40701022);  //WJets
  listOfWeights_.push_back(0.25768951);  //zjets
  listOfWeights_.push_back(0.00057319);  //singleTop_s
  listOfWeights_.push_back(0.00823326);  //singleTop_t
  listOfWeights_.push_back(0.00412381);  //singleTop_tW
  listOfWeights_.push_back(0.00402021);  //WW
  listOfWeights_.push_back(0.00164376);  //WZ
  listOfWeights_.push_back(0.00053283);  //ZZ
  listOfWeights_.push_back(13.05346354); //qcd_BCToE_30_80
  listOfWeights_.push_back(1.70806357);  //qcd_BCToE_80_170
  listOfWeights_.push_back(184.03691716); //qcd_EM_30_80
  listOfWeights_.push_back(3.67397446);  //qcd_EM_80_170
  listOfWeights_.push_back(8.41145800);  //qcdMu
  listOfWeights_.push_back(3.09952236);  //qcd_Ph_30_50
  listOfWeights_.push_back(0.50589264);  //qcd_Ph_50_80
  listOfWeights_.push_back(0.08123749);  //qcd_Ph_80_120
  listOfWeights_.push_back(0.01612019);  //qcd_Ph_120_170
  
  //mc scaling for lepton trigger
  double muScale(0.97), eleScale(0.99);
  */

//  TString jetObjectNoBTag_("PFlow/Bkg/btag0/lep_tau/PFlow lep_tau bkg 3dp clean");
  TString jetObjectNoBTag_("m_tau_DataDriven");
//  //TString jetObjectNoBTagExtra_("PFlow/Bkg/lep_tau/PFlow lep_tau bkg 3dp extra");
//  TString jetObjectBTag_("PFlow/Bkg/btag1/lep_tau/PFlow lep_tau2 bkg 3dp clean");
//  //TString jetObjectBTagExtra_("PFlow/Bkg/btag1/lep_tau/PFlow lep_tau2 bkg 3dp extra");
//
//  TString jetObjectNoBTagEle_("PFlow/Bkg/btag0/e_tau/PFlow e_tau bkg 3dp clean");
//  TString jetObjectNoBTagMu_("PFlow/Bkg/btag0/m_tau/PFlow m_tau bkg 3dp clean");
  TString jetObjectNoBTagMu_("m_tau_DataDriven");
//  TString jetObjectNoBTagExtraEle_("PFlow/Bkg/btag0/e_tau/PFlow e_tau bkg 3dp extra");
//  TString jetObjectNoBTagExtraMu_("PFlow/Bkg/btag0/m_tau/PFlow m_tau bkg 3dp extra");
  TString jetObjectNoBTagExtraMu_("m_tau_DataDriven");
//  TString jetObjectBTagEle_("PFlow/Bkg/btag1/e_tau/PFlow e_tau2 bkg 3dp clean");
//  TString jetObjectBTagMu_("PFlow/Bkg/btag1/m_tau/PFlow m_tau2 bkg 3dp clean");
//  TString jetObjectBTagExtraEle_("PFlow/Bkg/btag1/e_tau/PFlow e_tau2 bkg 3dp extra");
//  TString jetObjectBTagExtraMu_("PFlow/Bkg/btag1/m_tau/PFlow m_tau2 bkg 3dp extra");
  
  //Get Estimate from Data
  MVAComputer mva1(DataMVA_.c_str( ));
  ///  double fake_data_btag = TauFakeEstimate(listOfurls_[0], jetObjectBTag_, &mva1);
  double fake_data = 0, fake_data_err = 0;
  TauFakeEstimate(listOfurls_[0], jetObjectNoBTag_, &mva1, fr_data, &fake_data, &fake_data_err);
  cout<<"estimation from data : no-btag "<<fake_data<<" +/- "<<fake_data_err<<endl;///" 1 b-tag "<<fake_data_btag<<endl;
  
  //Get Estimation from MC
  MVAComputer mva2(mcMVA_.c_str());
  double fake_mc = 0, fake_mc_err = 0, fake_mc_btag = 0;
  for(size_t i = 1; i < listOfurls_.size(); i++){
    std::cout << "Processing file " << listOfurls_[i] << std::endl;
    //electron
    ///    double fake_btag = TauFakeEstimate(listOfurls_[i], jetObjectBTagEle_, &mva2);
    ///    fake_mc_btag += fake_btag; //(fake_btag*listOfWeights_[i]*eleScale);
    //muon
    ///    double fake_btag = TauFakeEstimate(listOfurls_[i], jetObjectBTagMu_, &mva2);
    ///    fake_mc_btag += fake_btag; //(fake_btag*listOfWeights_[i]*muScale);
    
    //electron
    ///    double fake = TauFakeEstimate(listOfurls_[i], jetObjectNoBTagEle_, &mva2);
    ///    fake_mc += fake; //(fake*listOfWeights_[i]*eleScale);
    //muon
    double fake = 0, fake_err = 0;
    TauFakeEstimate(listOfurls_[i], jetObjectNoBTagMu_, &mva2, fr_mc, &fake, &fake_err);
    fake_mc += fake; //(fake*listOfWeights_[i]*muScale);
    fake_mc_err += pow(fake_err,2);
  }
  
  cout<<"estimation from MC : no b-tag "<<fake_mc<<" +/- "<<sqrt(fake_mc_err)<<"    1 b-tag "<<fake_mc_btag<<endl;
  
  //Get Residual Estimation from MC
  double fake_mc_res = 0, fake_mc_res_err = 0, fake_mc_res_btag = 0;
  
  listOfurls_[1] = mcFolder+TString("out-ttbar_mcbkg.root");
  listOfurls_.push_back( mcFolder+TString("out-ttbar_mutau.root"));
//  listOfurls_.clear();
//  listOfurls_.push_back( mcFolder+TString("out-ttbar_mcbkg.root") );
  
//  for(size_t i = 0; i < listOfurls_.size(); i++){
  for(size_t i = 1; i < listOfurls_.size(); i++){
    //electron
    ///    double fake_btag = TauFakeEstimate(listOfurls_[i], jetObjectBTagExtraEle_, &mva2);
    ///    fake_mc_res_btag += fake_btag; //(fake_btag*listOfWeights_[i]*eleScale);
    //muon
    ///   fake_btag = TauFakeEstimate(listOfurls_[i], jetObjectBTagExtraMu_, &mva2);
    ///   fake_mc_res_btag += fake_btag; //(fake_btag*listOfWeights_[i]*muScale);
    
    //electron
    /// double fake = TauFakeEstimate(listOfurls_[i], jetObjectNoBTagExtraEle_, &mva2);
    ///    fake_mc_res += fake; //(fake*listOfWeights_[i]*eleScale);
    //muon
    if(i!=2){
      double fake = 0, fake_err = 0;
      TauFakeEstimate(listOfurls_[i], jetObjectNoBTagExtraMu_, &mva2, fr_mc, &fake, &fake_err);
          fake_mc_res += fake; //(fake*listOfWeights_[i]*muScale);
	  fake_mc_res_err += pow(fake_err,2);
    }
  }
  
  cout<<"Residual estimation from MC : no b-tag "<<fake_mc_res<<" +/- "<<sqrt(fake_mc_res_err) <<"    1 b-tag "<<fake_mc_res_btag<<endl;

  finalValues.push_back(fake_data);
  finalValues.push_back(fake_data_err);
  finalValues.push_back(fake_mc);
  finalValues.push_back(sqrt(fake_mc_err));
  finalValues.push_back(fake_mc_res);
  finalValues.push_back(sqrt(fake_mc_res_err));
   
}

void TauFakesHelper::TauFakeEstimate(TString filename_, TString jetObject_, PhysicsTools::MVAComputer *mva, TH1F* frHisto, double* sumFakes, double* sumFakesError)
{
  using namespace PhysicsTools;
  
  TFile *f = TFile::Open(filename_);
  if(f==0) return;
  if(f->IsZombie()) { f->Close(); return ; }


///  TPolyMarker3D *data_ = (TPolyMarker3D*) f->Get(jetObject_);
///  TString weightObject_ = jetObject_+" weights";
///  TPolyMarker3D *weight_ = (TPolyMarker3D*) f->Get(weightObject_)->Clone("weight_");
///
///  int nev = data_->Size();
///  //cout<<"events to process "<<nev<<endl;
///  if(nev <= 0) return 0;
///
///  double JetPt, JetWidth, JetEta;
  double wtMCScale, wtPileup, wtLepEff;

// -----
  TTree *tree = (TTree *) f->Get(jetObject_);
  double JetPt, JetWidth, AbsJetEta, Weight;
  //  tree->SetBranchAddress("JetPt", &JetPt);
  //  tree->SetBranchAddress("JetWidth", &JetWidth);
  //  tree->SetBranchAddress("AbsJetEta", &AbsJetEta);
  //  tree->SetBranchAddress("__WEIGHT__", &Weight);
  tree->SetBranchAddress("pt_t", &JetPt);
  tree->SetBranchAddress("radius_t", &JetWidth);
  tree->SetBranchAddress("eta_t", &AbsJetEta);
  tree->SetBranchAddress("weight", &Weight);
  Int_t nentries = (Int_t)tree->GetEntries();
  cout<<"entries "<<nentries<<endl;
  //int nev = (nentries < 2000000) ? nentries : 2000000;

  int nev = nentries;
  // -----

  int index = 0;
  double sumEvent = 0;
  double sumEvent_err = 0;
  double sumUnw = 0;
  double checkPt=0, checkEta=0, checkWidth=0;
  for(Int_t ev = 0; ev < nev; ev++){
///    data_->GetPoint(ev, JetPt, JetEta, JetWidth);
///    weight_->GetPoint(ev, wtMCScale, wtPileup, wtLepEff);
    tree->GetEntry(ev);
    double wtTotal = Weight; ///wtMCScale*wtPileup*wtLepEff;
    //cout<<"pt "<<JetPt<<" eta "<<JetEta<<" r "<<JetWidth<<endl;
    //cout<<wtMCScale<<" "<<wtPileup<<" "<<wtLepEff<<" "<<wtTotal<<endl;
    if(JetPt==checkPt || AbsJetEta == checkEta || JetWidth == checkWidth)cout<<"repeated"<<endl;

    Variable::ValueList vars;
    vars.add("JetPt", JetPt);
    ///    vars.add("AbsJetEta", fabs(JetEta));
    vars.add("AbsJetEta", AbsJetEta);
    vars.add("JetWidth", JetWidth);
    
    double value = mva->eval(vars);

    Int_t ibin = frHisto->FindBin(JetPt);
    if(JetPt > 200) continue; // FIXME: test this

    sumEvent += value*wtTotal;
    sumUnw += wtTotal;
    double Err = wtTotal * value*(      frHisto->GetBinError(ibin) / frHisto->GetBinContent(ibin)     );
    sumEvent_err += pow(Err,2);
    //cout << value << endl;
    index++;
  }
  
  double finErr = sqrt(sumEvent_err);
  (*sumFakes) = sumEvent;
  (*sumFakesError) = finErr;
  cout<<"events pass "<<sumEvent<<", squaredError " << sumEvent_err << " unweighted "<<sumUnw<<endl;
  
  f->Close();

  return;// sumEvent;
}


void TauFakesHelper::ProduceDataDrivenDistributions(){
  
  // Example of Root macro to copy a subset of a Tree to a new Tree
  // Only selected entries are copied to the new Tree.
  // The input file has been generated by the program in $ROOTSYS/test/Event
  // with   Event 1000 1 99 1
  //Author: Rene Brun
  //   
  //   gSystem->Load("$ROOTSYS/test/libEvent");
  //
  //   //Get old file, old tree and set top branch address
  //   TFile *oldfile = new TFile("$ROOTSYS/test/Event.root");
  //   TTree *oldtree = (TTree*)oldfile->Get("T");
  //   Long64_t nentries = oldtree->GetEntries();
  //   Event *event   = 0;
  //   oldtree->SetBranchAddress("event",&event);
  //
  //   //Create a new file + a clone of old tree in new file
  //   TFile *newfile = new TFile("small.root","recreate");
  //   TTree *newtree = oldtree->CloneTree(0);

  using namespace PhysicsTools;
  
  string WMuMVA(trainingOutputArea_+"/Trained_QCDData_PFlow.mva");
  string QCDMVA(trainingOutputArea_+"/Trained_WMuData_PFlow.mva");
  
  
  cout << "Rescaling data using:\n \t\t" << WMuMVA<<"\n \t\t "<<QCDMVA<<endl;
  
  vector<TString> listOfurls_; listOfurls_.clear();
  vector<TString> rescaledListOfurls_; rescaledListOfurls_.clear();
  
  listOfurls_.push_back(ntuplesArea_+"/nomt-2012-V1-data-MU-20GeV/"+TString("out-data.root"));  
  
  rescaledListOfurls_.push_back(ntuplesArea_+"/nomt-2012-V1-data-MU-20GeV/"+TString("out-data_rescaled.root"));  
  
  //  TString jetObject_("m_tau_DataDriven");
  std::vector<TString>jetObject_;
  jetObject_.push_back("m_tau_Selected");
  jetObject_.push_back("m_tau_DataDriven");

  //  // Syst // FIXME: look for additional ones (tauenergy?)
  //  jetObject_.push_back("m_tau_Selected_plus"     );
  //  jetObject_.push_back("m_tau_Selected_minus"    );
  //  jetObject_.push_back("m_tau_Selected_uncplus"  );
  //  jetObject_.push_back("m_tau_Selected_uncminus" );
  //  jetObject_.push_back("m_tau_Selected_jerplus"  );
  //  jetObject_.push_back("m_tau_Selected_jerminus" );
  
  
  //Get Estimate from WMu
  MVAComputer mvaWMu(WMuMVA.c_str());
  //Get Estimation from QCD
  MVAComputer mvaQCD(QCDMVA.c_str());
  
  for(size_t i = 0; i < listOfurls_.size(); i++){
    
    // get the two, scale each
    
    TFile *f = TFile::Open(listOfurls_[i]);
    if(f==0) return;
    if(f->IsZombie()) { f->Close(); return; }
    
    TFile *newfile = new TFile(rescaledListOfurls_[i],"recreate");
    
    double wtMCScale, wtPileup, wtLepEff;
    
    // -----
    for(size_t tr=0; tr<jetObject_.size(); tr++){
      
      
      TTree *tree = (TTree *) f->Get(jetObject_[tr]);
      //    TFile *newfile = new TFile(rescaledListOfurls_,"recreate");
      cout << "New file opening...";

      cout << "succeeded." << endl;
      TTree *newtree = tree->CloneTree(0);
      
      double JetPt, JetWidth, AbsJetEta, Weight, IsOS, RcT;
      
      double LepPt, MetPt, LepMetDPhi, JetMulti, JetBMulti, TauMetDPhi;
      
      //  tree->SetBranchAddress("JetPt", &JetPt);
      //  tree->SetBranchAddress("JetWidth", &JetWidth);
      //  tree->SetBranchAddress("AbsJetEta", &AbsJetEta);
      //  tree->SetBranchAddress("__WEIGHT__", &Weight);
      tree->SetBranchAddress("pt_t", &JetPt);
      tree->SetBranchAddress("radius_t", &JetWidth);
      tree->SetBranchAddress("eta_t", &AbsJetEta);
      tree->SetBranchAddress("weight", &Weight);
      tree->SetBranchAddress("is_os", &IsOS);
      tree->SetBranchAddress("rc_t", &RcT);
      
      
      //    tree->SetBranchAddress(nJetPt;
      //    tree->SetBranchAddress("pt_t");               "pt_t", &nJetPt);    
      tree->SetBranchAddress("pt_l",               &LepPt);
      tree->SetBranchAddress("pt_met",             &MetPt);
      tree->SetBranchAddress("Dphi_lepton_met",    &LepMetDPhi);
      tree->SetBranchAddress("multiplicity_j",     &JetMulti);
      tree->SetBranchAddress("btagmultiplicity_j", &JetBMulti);
      tree->SetBranchAddress("Dphi_tau_met",       &TauMetDPhi);
      
      
      Int_t nentries = (Int_t)tree->GetEntries();
      cout<<"entries "<<nentries<<endl;
      int nev = nentries; // Always rescale everything, otherwise yields get fucked upB
      // -----
      
      int index = 0;
      double sumEvent = 0;
      double sumUnw = 0;
      for(Int_t ev = 0; ev < nev; ev++){
	if (ev%1000==0) cout<<listOfurls_[i]<<": entry "<<ev<<"/"<<nev<<endl;
	
	///    data_->GetPoint(ev, JetPt, JetEta, JetWidth);
	///    weight_->GetPoint(ev, wtMCScale, wtPileup, wtLepEff);
	tree->GetEntry(ev);
	double wtTotal = Weight; ///wtMCScale*wtPileup*wtLepEff;      
	//cout<<"pt "<<JetPt<<" eta "<<JetEta<<" r "<<JetWidth<<endl;
	//cout<<wtMCScale<<" "<<wtPileup<<" "<<wtLepEff<<" "<<wtTotal<<endl;
	Variable::ValueList vars;
	vars.add("JetPt", JetPt);
	///    vars.add("AbsJetEta", fabs(JetEta));
	vars.add("AbsJetEta", AbsJetEta);
	vars.add("JetWidth", JetWidth);
	
	double wmuValue =mvaWMu.eval(vars);
	double qcdValue =mvaQCD.eval(vars);
	
	double value = (wmuValue + qcdValue )/2.;
	
	//JetPt       *=value;
	//AbsJetEta   *=value;
	//JetWidth    *=value;
	//LepPt       *=value;
	//MetPt       *=value;
	//LepMetDPhi  *=value;
	//JetMulti    *=value;
	//JetBMulti   *=value;
	//TauMetDPhi  *=value;       
	
	Weight *=value;
	
	sumEvent += Weight;
	sumUnw += wtTotal;
	//      cout << "WMu: " << wmuValue << ", QCD: " << qcdValue << ", TotaL: " << value << ", wtbase: " << wtTotal << ", wtFinal: " << value*wtTotal << endl;
	index++;
	
	//if(value != 0) newtree->Fill();
	newtree->Fill();
	//event->Clear();
	
      }
      
    
      newtree->Print();
      newtree->AutoSave();
      
      
     

      tree->ResetBranchAddresses();
      //    return sumEvent;
      cout<< "SumEvent: " << sumEvent <<", sumUnw: " << sumUnw << endl;
    }    
     f->Close();
  }

}




