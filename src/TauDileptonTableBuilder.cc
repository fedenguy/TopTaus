#include "LIP/TopTaus/interface/TauDileptonTableBuilder.hh"

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "LIP/TopTaus/interface/CommonDefinitions.hh"
// System headers
#include <cmath>
#include <map>
#include <fstream>

#endif

// ROOT headers
#include <TLegendEntry.h>
#include <TFile.h>
#include <TH2D.h>

using namespace commondefinitions;  

// FIXME: clean code
namespace tableutils{

  enum dataPrecision { D1=1, D2=2, D3=3 };
  
  
  // systset1 -> if enabled include only JES + MET +JER ,  systset2 -> if enabled include only btag unc ,  systset3 -> if enabled only trigger
  // detailed -> select between 1 or 2 floating points
  void TauDileptonTableBuilder::mcTable( int detailed, bool includeSoverB , bool printAllErrors, bool higgs, TString key, TString name, bool systset1, bool systset2, bool systset3){ 


    //    if(!doTheDatacards) // For datacards we want 234., for tables we want 1.1 // This does not need conditional, because of the default values
      for(size_t i=0; i<brHtaunu_.size(); ++i)
	brHtaunu_[i] = 5./234.;
    
    /*
      cout<<endl;
      cout<<"RATIO on data   : "<<306./455<<" +- "<<getErrorFraction(306.,432.);
      cout<<"RATIO on MC     : "<<279./432.<<" +- "<<getErrorFraction(279.,432.,17.,22);
      cout<<"RATIO on w+jets : "<<207.48/330.77<<" +- "<<getErrorFraction(207.48,330.77,16.77,21.18);
    */
    
    
    
    
    bool triggerunc(true);  //enable/disable syst on trigger
    
    int STARTINGPOINT = STARTING_AT_LJETS_;         // starts selection after requiring lepton+3jets
    //int STARTINGPOINT = STARTING_AT_LJETSPLUSMET_;    // starts selection after requiring lepton+3jets+MHT
    
    // fw and fhw used in H+ analysis to get out the hh and hw contribution for a given BR
    // hh contribution : x2 , hw contribution :2.(1-x)x , whh and fhw can be set as funcion of x (ex 5%)
    // 5% contribution //
    //  double fhh(0.05*0.05) , fhw( 2*(1-0.05)*0.05) ;
    double fhh(1),fhw(1); // this is for the datacards (so that they are normalized to 234 (xsecttbar)
    // for ftbh see below (search for "vector<double> ftbh")
    //    double ftbh(1); // This is for the datacards, so that they are normalized to 234


     
    bool XSECMEASUREMENT( true );  if( XSECMEASUREMENT ){ TTBAR_CS_ERR = XSEC_EXP_ERR_; } else { SIGNAL_EFF_ERR = 0; }
    
    // select wich runing mode we are in ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int TAUSTEP,BTAG1STEP; int COL1,COL2,COL3,COL4,COL5,COL6,COL7,COL8,COL9;
    
    TAUSTEP = TAU_STEP2; BTAG1STEP = BTAG1_STEP2;
    // nomt //  if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;  COL4=BTAG1_STEP2; COL5=TAU_STEP2; COL6=OS_STEP2;    COL7=BTAG2_STEP2; COL8=R_STEP1; }
    // nomt //  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;   COL3=BTAG1_STEP2 ;  COL4=TAU_STEP2;   COL5=OS_STEP2;  COL6=BTAG2_STEP2; COL7=R_STEP2;                   }
    
    if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;  COL4=MT_STEP2;      COL5=BTAG1_STEP2; COL6=TAU_STEP2; COL7=OS_STEP2;    COL8=BTAG2_STEP2; COL9=R_STEP1; }
    else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;  COL3=MT_STEP2;   COL4=BTAG1_STEP2 ;  COL5=TAU_STEP2;   COL6=OS_STEP2;  COL7=BTAG2_STEP2; COL8=R_STEP2;                   }
    
    
    // table options ////////////////////////////
    bool incDocument(true);
    bool includeSM,hh_contrib,hw_contrib,tbh_contrib(true);
  ////////////////////////////////////////////
    
    
    if(higgs){ 
      hh_contrib = false; // now doing heavy
      hw_contrib = false; // now doing heavy 
      tbh_contrib = true;
      includeSM = false; }
    else     { hh_contrib = false; hw_contrib = false; tbh_contrib = false;
      includeSM = true; }
    
    
    vector<TString> hwEffHitoNames;
    hwEffHitoNames.push_back("WH80");  hwEffHitoNames.push_back("WH100"); hwEffHitoNames.push_back("WH120"); hwEffHitoNames.push_back("WH140"); 
    hwEffHitoNames.push_back("WH150"); hwEffHitoNames.push_back("WH155"); hwEffHitoNames.push_back("WH160"); 
    
    vector<TString> hhEffHitoNames;
    hhEffHitoNames.push_back("HH80");  hhEffHitoNames.push_back("HH100"); hhEffHitoNames.push_back("HH120"); hhEffHitoNames.push_back("HH140"); 
    hhEffHitoNames.push_back("HH150"); hhEffHitoNames.push_back("HH155");
    hhEffHitoNames.push_back("HH160"); 

    vector<TString> tbhEffHitoNames;
    hwEffHitoNames.push_back("TBH180"); 
    hwEffHitoNames.push_back("TBH190"); 
    hwEffHitoNames.push_back("TBH200"); 
    hwEffHitoNames.push_back("TBH220"); 
    hwEffHitoNames.push_back("TBH250");
    hwEffHitoNames.push_back("TBH300");


    vector<TString> hhFiles; vector<TString> hwFiles;
    vector<TString> tbhFiles; 
    if(hh_contrib){
      hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m80.root"));   hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m100.root")); hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m120.root"));   
      hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m140.root"));  hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m150.root")); hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m155.root"));   
      hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m160.root"));
  }
  
  if(hw_contrib){
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m80.root"));  hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m100.root")); hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m120.root")); 
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m140.root")); hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m150.root")); hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m155.root"));   
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m160.root"));
  }

  if(tbh_contrib){
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m180.root")); 
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m190.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m200.root")); 
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m220.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m250.root")); 
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m300.root"));   
  }


  TString pFile(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-mc.root"));       TFile * processedFile     = TFile::Open(pFile); 
  TString dFile(outputArea_+TString("nomt-2012-V1-data-MU-20GeV/out-data.root"));   TFile * processedDataFile = TFile::Open(dFile); 

  vector<TFile *> processedHWFiles; vector<TFile *> processedHHFiles;
  vector<TFile*> processedTBHFiles;
  for(size_t i = 0; i<hhFiles.size(); i++ ){  processedHHFiles.push_back( TFile::Open( hhFiles[i] ) ); }
  for(size_t i = 0; i<hwFiles.size(); i++ ){  processedHWFiles.push_back( TFile::Open( hwFiles[i] ) ); }
  for(size_t i = 0; i<tbhFiles.size(); i++ ){  processedTBHFiles.push_back( TFile::Open( tbhFiles[i] ) ); }

  // Names for data histograms 
  TString opString,leptauString,etauString,mtauString, leptauStringMC;

  leptauStringMC=TString("1");   leptauString=TString("1"); etauString=TString("2"); mtauString=TString("3"); 
  TString hNameData;    hNameData   = key + TString("/Yields")+ TString("/lep_tau/cutflow_yields_") + leptauString;
  TString hNameData_e;  hNameData_e = key + TString("/Yields")+ TString("/e_tau/cutflow_yields_")   + etauString;
  TString hNameData_m;  hNameData_m = key + TString("/Yields")+ TString("/m_tau/cutflow_yields_")   + mtauString ;

  //Names for total MC yields histograms (2D histograms)
  TString hName;             hName          = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_")          + leptauStringMC; 
  TString hName_plus;        hName_plus     = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_plus_")     + leptauStringMC; 
  TString hName_minus;       hName_minus    = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_minus_")    + leptauStringMC; 
  TString hName_uncplus;     hName_uncplus  = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_uncplus_")  + leptauStringMC;
  TString hName_uncminus;    hName_uncminus = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_uncminus_") + leptauStringMC;
  TString hName_jerplus;     hName_jerplus  = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_jerplus_")  + leptauStringMC; 
  TString hName_jerminus;    hName_jerminus = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_jerminus_") + leptauStringMC;
  TString hName_bminus;      hName_bminus   = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_bminus_")   + leptauStringMC; 
  TString hName_bplus;       hName_bplus    = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_bplus_")    + leptauStringMC; 
  TString hName_trigger;     hName_trigger  = key+TString("/YieldsMC")+TString("/cutflow_mc_triggErr_")        + leptauStringMC; 

  //Names for total MC electron/muon yields histograms
  TString hName_higgs_e_tau;           hName_higgs_e_tau          = key+TString("/Yields") + TString("/e_tau/cutflow_yields_")          + etauString;
  TString hName_higgs_e_tau_plus;      hName_higgs_e_tau_plus     = key+TString("/Yields") + TString("/e_tau/cutflow_yields_plus_")     + etauString;
  TString hName_higgs_e_tau_minus;     hName_higgs_e_tau_minus    = key+TString("/Yields") + TString("/e_tau/cutflow_yields_minus_")    + etauString;
  TString hName_higgs_e_tau_uncplus;   hName_higgs_e_tau_uncplus  = key+TString("/Yields") + TString("/e_tau/cutflow_yields_uncplus_")  + etauString;
  TString hName_higgs_e_tau_uncminus;  hName_higgs_e_tau_uncminus = key+TString("/Yields") + TString("/e_tau/cutflow_yields_uncminus_") + etauString;
  TString hName_higgs_e_tau_jerplus;   hName_higgs_e_tau_jerplus  = key+TString("/Yields") + TString("/e_tau/cutflow_yields_jerplus_")  + etauString; 
  TString hName_higgs_e_tau_jerminus;  hName_higgs_e_tau_jerminus = key+TString("/Yields") + TString("/e_tau/cutflow_yields_jerminus_") + etauString;
  TString hName_higgs_e_tau_bplus;     hName_higgs_e_tau_bplus    = key+TString("/Yields") + TString("/e_tau/cutflow_yields_bplus_")    + etauString; 
  TString hName_higgs_e_tau_bminus;    hName_higgs_e_tau_bminus   = key+TString("/Yields") + TString("/e_tau/cutflow_yields_bminus_")   + etauString;
  TString hName_higgs_e_tau_trigger;   hName_higgs_e_tau_trigger  = key+TString("/Yields") + TString("/e_tau/cutflow_triggErr_")        + etauString;
  TString hName_higgs_m_tau;           hName_higgs_m_tau          = key+TString("/Yields") + TString("/m_tau/cutflow_yields_")          + mtauString;
  TString hName_higgs_m_tau_plus;      hName_higgs_m_tau_plus     = key+TString("/Yields") + TString("/m_tau/cutflow_yields_plus_")     + mtauString;
  TString hName_higgs_m_tau_minus;     hName_higgs_m_tau_minus    = key+TString("/Yields") + TString("/m_tau/cutflow_yields_minus_")    + mtauString;     
  TString hName_higgs_m_tau_uncplus;   hName_higgs_m_tau_uncplus  = key+TString("/Yields") + TString("/m_tau/cutflow_yields_uncplus_")  + mtauString;  
  TString hName_higgs_m_tau_uncminus;  hName_higgs_m_tau_uncminus = key+TString("/Yields") + TString("/m_tau/cutflow_yields_uncminus_") + mtauString; 
  TString hName_higgs_m_tau_jerplus;   hName_higgs_m_tau_jerplus  = key+TString("/Yields") + TString("/m_tau/cutflow_yields_jerplus_")  + mtauString;  
  TString hName_higgs_m_tau_jerminus;  hName_higgs_m_tau_jerminus = key+TString("/Yields") + TString("/m_tau/cutflow_yields_jerminus_") + mtauString;
  TString hName_higgs_m_tau_bplus;     hName_higgs_m_tau_bplus    = key+TString("/Yields") + TString("/m_tau/cutflow_yields_bplus_")    + mtauString;  
  TString hName_higgs_m_tau_bminus;    hName_higgs_m_tau_bminus   = key+TString("/Yields") + TString("/m_tau/cutflow_yields_bminus_")   + mtauString;
  TString hName_higgs_m_tau_trigger;   hName_higgs_m_tau_trigger  = key+TString("/Yields") + TString("/m_tau/cutflow_triggErr_")        + mtauString;

  // MC split channel plots
  TH2D * h(0); TH2D * h_plus(0); TH2D * h_minus(0); 
  TH2D * h_uncplus(0); TH2D * h_uncminus(0); TH2D * h_jerplus(0); TH2D * h_jerminus(0); TH2D * h_bplus(0); TH2D * h_bminus(0);TH2D * h_trigger(0);

  // data plots
  TH1D * hData(0) ; TH1D * hData_e(0); TH1D * hData_m(0);
  // MC e/m channels
  TH1D * h_e(0); TH1D * h_m(0);


  vector< TH1D * > h_hh_e, h_hh_e_plus, h_hh_e_minus, h_hh_e_uncplus, h_hh_e_uncminus, h_hh_e_jerplus, h_hh_e_jerminus, h_hh_e_bplus, h_hh_e_bminus, h_hh_e_trigger;
  vector< TH1D * > h_hh_m, h_hh_m_plus, h_hh_m_minus, h_hh_m_uncplus, h_hh_m_uncminus, h_hh_m_jerplus, h_hh_m_jerminus, h_hh_m_bplus, h_hh_m_bminus, h_hh_m_trigger;
    
  vector< TH1D * > h_hw_e, h_hw_e_plus, h_hw_e_minus, h_hw_e_uncplus, h_hw_e_uncminus, h_hw_e_jerplus, h_hw_e_jerminus, h_hw_e_bplus, h_hw_e_bminus, h_hw_e_trigger; 
  vector< TH1D * > h_hw_m, h_hw_m_plus, h_hw_m_minus, h_hw_m_uncplus, h_hw_m_uncminus, h_hw_m_jerplus, h_hw_m_jerminus, h_hw_m_bplus, h_hw_m_bminus, h_hw_m_trigger;

  vector< TH1D * > h_tbh_e,
    h_tbh_e_plus,
    h_tbh_e_minus,
    h_tbh_e_uncplus,
    h_tbh_e_uncminus,
    h_tbh_e_jerplus,
    h_tbh_e_jerminus,
    h_tbh_e_bplus,
    h_tbh_e_bminus,
    h_tbh_e_trigger; 
  vector< TH1D * > h_tbh_m, 
    h_tbh_m_plus, 
    h_tbh_m_minus,
    h_tbh_m_uncplus,
    h_tbh_m_uncminus,
    h_tbh_m_jerplus,
    h_tbh_m_jerminus, 
    h_tbh_m_bplus, 
    h_tbh_m_bminus, 
    h_tbh_m_trigger;

  if(processedFile){ 
    h          = (TH2D*) processedFile->Get(hName);  
    h_plus     = (TH2D*) processedFile->Get(hName_plus); 
    h_minus    = (TH2D*) processedFile->Get(hName_minus); 
    h_uncplus  = (TH2D*) processedFile->Get(hName_uncplus);
    h_uncminus = (TH2D*) processedFile->Get(hName_uncminus);
    h_jerplus  = (TH2D*) processedFile->Get(hName_jerplus);
    h_jerminus = (TH2D*) processedFile->Get(hName_jerminus);
    h_bplus    = (TH2D*) processedFile->Get(hName_bplus);
    h_bminus   = (TH2D*) processedFile->Get(hName_bminus);
    if (triggerunc )h_trigger  = (TH2D*) processedFile->Get(hName_trigger); 
    h_e        = (TH1D*) processedFile->Get(hName_higgs_e_tau); // use the same strings for higgs
    h_m        = (TH1D*) processedFile->Get(hName_higgs_m_tau); // use the same strings for higgs
  } 

  if(processedDataFile){ hData =(TH1D *) processedDataFile->Get(hNameData); hData_e=(TH1D *) processedDataFile->Get(hNameData_e); hData_m=(TH1D *) processedDataFile->Get(hNameData_m);}
 
  for(uint i = 0; i< hhFiles.size(); i++ ){ 

    if(processedHHFiles[i]){ 
      h_hh_e.push_back(         (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau)         ); h_hh_m.push_back(         (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau)         );
      h_hh_e_plus.push_back(    (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_plus)    ); h_hh_m_plus.push_back(    (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_plus)    );  
      h_hh_e_minus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_minus)   ); h_hh_m_minus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_minus)   );  
      h_hh_e_uncplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_uncplus) ); h_hh_m_uncplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_uncplus) );   
      h_hh_e_uncminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_uncminus)); h_hh_m_uncminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_uncminus));
      h_hh_e_jerplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_jerplus) ); h_hh_m_jerplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_jerplus) );   
      h_hh_e_jerminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_jerminus)); h_hh_m_jerminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_jerminus));  
      h_hh_e_bplus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_bplus) );   h_hh_m_bplus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_bplus)   );   
      h_hh_e_bminus.push_back(  (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_bminus));   h_hh_m_bminus.push_back(  (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_bminus)  );  
      if (triggerunc ){ 
        h_hh_e_trigger.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_trigger));  h_hh_m_trigger.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_trigger) );  
      }else{ h_hh_e_trigger.push_back(0);  h_hh_m_trigger.push_back(0);  }
    }
    else { 
      h_hh_e.push_back(0);  h_hh_e_plus.push_back(0); h_hh_e_minus.push_back(0); h_hh_e_trigger.push_back(0); 
      h_hh_e_uncplus.push_back(0); h_hh_e_uncminus.push_back(0); h_hh_e_jerplus.push_back(0); h_hh_e_jerminus.push_back(0); h_hh_e_bplus.push_back(0);   h_hh_e_bminus.push_back(0);

      h_hh_m.push_back(0);  h_hh_m_plus.push_back(0); h_hh_m_minus.push_back(0); h_hh_m_trigger.push_back(0);
      h_hh_m_uncplus.push_back(0); h_hh_m_uncminus.push_back(0); h_hh_m_jerplus.push_back(0); h_hh_m_jerminus.push_back(0); h_hh_m_bplus.push_back(0); h_hh_m_bminus.push_back(0); 

    }
  }
  for(uint i = 0; i< hwFiles.size(); i++ ){ 


    if(processedHWFiles[i]){ 
      h_hw_e.push_back(         (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau)         ); h_hw_m.push_back(         (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau)         );  
      h_hw_e_plus.push_back(    (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_plus)    ); h_hw_m_plus.push_back(    (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_plus)    ); 
      h_hw_e_minus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_minus)   ); h_hw_m_minus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_minus)   ); 
      h_hw_e_uncplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_uncplus) ); h_hw_m_uncplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_uncplus) ); 
      h_hw_e_uncminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_uncminus)); h_hw_m_uncminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_uncminus)); 
      h_hw_e_jerplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_jerplus) ); h_hw_m_jerplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_jerplus) ); 
      h_hw_e_jerminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_jerminus)); h_hw_m_jerminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_jerminus)); 
      h_hw_e_bplus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_bplus) );   h_hw_m_bplus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_bplus)   ); 
      h_hw_e_bminus.push_back(  (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_bminus));   h_hw_m_bminus.push_back(  (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_bminus)  ); 
      if(triggerunc ){
        h_hw_e_trigger.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_trigger));  h_hw_m_trigger.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_trigger));
      }else{ h_hw_e_trigger.push_back(0);  h_hw_m_trigger.push_back(0);}
    }
    else { 
      h_hw_e.push_back(0); h_hw_e_plus.push_back(0); h_hw_e_minus.push_back(0); h_hw_e_trigger.push_back(0);
      h_hw_e_uncplus.push_back(0); h_hw_e_uncminus.push_back(0); h_hw_e_jerplus.push_back(0); h_hw_e_jerminus.push_back(0); h_hw_e_bplus.push_back(0); h_hw_e_bminus.push_back(0);

      h_hw_m.push_back(0); h_hw_m_plus.push_back(0); h_hw_m_minus.push_back(0); h_hw_m_trigger.push_back(0);
      h_hw_m_uncplus.push_back(0); h_hw_m_uncminus.push_back(0); h_hw_m_jerplus.push_back(0); h_hw_m_jerminus.push_back(0); h_hw_m_bplus.push_back(0); h_hw_m_bminus.push_back(0);

    }
  }


  for(uint i = 0; i< tbhFiles.size(); i++ ){ 


    if(processedTBHFiles[i]){ 
      h_tbh_e.push_back(         (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau)         ); h_tbh_m.push_back(         (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau)         );  
      h_tbh_e_plus.push_back(    (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_plus)    ); h_tbh_m_plus.push_back(    (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_plus)    ); 
      h_tbh_e_minus.push_back(   (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_minus)   ); h_tbh_m_minus.push_back(   (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_minus)   ); 
      h_tbh_e_uncplus.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_uncplus) ); h_tbh_m_uncplus.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_uncplus) ); 
      h_tbh_e_uncminus.push_back((TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_uncminus)); h_tbh_m_uncminus.push_back((TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_uncminus)); 
      h_tbh_e_jerplus.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_jerplus) ); h_tbh_m_jerplus.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_jerplus) ); 
      h_tbh_e_jerminus.push_back((TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_jerminus)); h_tbh_m_jerminus.push_back((TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_jerminus)); 
      h_tbh_e_bplus.push_back(   (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_bplus) );   h_tbh_m_bplus.push_back(   (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_bplus)   ); 
      h_tbh_e_bminus.push_back(  (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_bminus));   h_tbh_m_bminus.push_back(  (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_bminus)  ); 
      if(triggerunc ){
        h_tbh_e_trigger.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_e_tau_trigger));  h_tbh_m_trigger.push_back( (TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m_tau_trigger));
      }else{ h_tbh_e_trigger.push_back(0);  h_tbh_m_trigger.push_back(0);}
    }
    else { 
      h_tbh_e.push_back(0); h_tbh_e_plus.push_back(0); h_tbh_e_minus.push_back(0); h_tbh_e_trigger.push_back(0);
      h_tbh_e_uncplus.push_back(0); h_tbh_e_uncminus.push_back(0); h_tbh_e_jerplus.push_back(0); h_tbh_e_jerminus.push_back(0); h_tbh_e_bplus.push_back(0); h_tbh_e_bminus.push_back(0);

      h_tbh_m.push_back(0); h_tbh_m_plus.push_back(0); h_tbh_m_minus.push_back(0); h_tbh_m_trigger.push_back(0);
      h_tbh_m_uncplus.push_back(0); h_tbh_m_uncminus.push_back(0); h_tbh_m_jerplus.push_back(0); h_tbh_m_jerminus.push_back(0); h_tbh_m_bplus.push_back(0); h_tbh_m_bminus.push_back(0);

    }
  }


  double dataTotal[14], dataTotal_e[14], dataTotal_m[14], syst_error_plus[14], syst_error_minus[14];
   
  // copy all data to matrix ( x = 14 steps of selection ,y =16 = 13 ev types + 1 total + 2 for S/B and S/sqrt B )
  double d[14][16], e[14][16], d_plus[14][16], d_minus[14][16], 
         d_uncplus[14][16],    d_uncminus[14][16], 
         d_jerplus[14][16],    d_jerminus[14][16], 
         d_bplus[14][16],      d_bminus[14][16], 
         d_trigger[14][16],
         systErrPlus[14][16],  systErrMinus[14][16];
         
  // split by e channel and mu channel /////////////////////////////////////////////////////////////////////////
  double d_e[14], d_m[14], d_e_systErrMinus[14], d_e_systErrPlus[14], d_m_systErrPlus[14], d_m_systErrMinus[14];
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////



  // fill HH samples and WH matrix with data and errors //////////////////////////////////////////////////////////////////////////////////////////////////////
  map< int, vector<double> > hh_e, hh_m, hw_e, hw_m, hh_e_err, hh_m_err, hw_e_err, hw_m_err;
  map< int, vector<double> > hh_e_syst_plus, hh_e_syst_minus, hh_m_syst_plus, hh_m_syst_minus, hw_e_syst_plus, hw_e_syst_minus, hw_m_syst_plus, hw_m_syst_minus;


  for(uint samp=0; samp<h_hh_e.size(); samp++){
    
    vector<double> data_hh_e, data_hh_m, data_hh_e_err, data_hh_m_err, data_syst_e_plus, data_syst_e_minus, data_syst_m_plus, data_syst_m_minus; 
    for(int x=0; x<14;x++){
      
      // get data and stats errors ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double data_hh_e_plus,    data_hh_e_minus,    data_hh_m_plus,    data_hh_m_minus;
      double data_hh_e_uncplus, data_hh_e_uncminus, data_hh_m_uncplus, data_hh_m_uncminus;
      double data_hh_e_jerplus, data_hh_e_jerminus, data_hh_m_jerplus, data_hh_m_jerminus;
      double data_hh_e_bplus,   data_hh_e_bminus,   data_hh_m_bplus,   data_hh_m_bminus;
      double data_hh_e_trigger, data_hh_m_trigger;
      
      if( h_hh_e[samp]          !=0 ){ data_hh_e.push_back(     h_hh_e[samp]->GetBinContent(x+1) );     } else { data_hh_e.push_back( 0 );     }
      if( h_hh_e[samp]          !=0 ){ data_hh_e_err.push_back( h_hh_e[samp]->GetBinError(x+1)   );     } else { data_hh_e_err.push_back( 0 ); }
      if( h_hh_m[samp]          !=0 ){ data_hh_m.push_back(     h_hh_m[samp]->GetBinContent(x+1) );     } else { data_hh_m.push_back( 0 );     }
      if( h_hh_m[samp]          !=0 ){ data_hh_m_err.push_back( h_hh_m[samp]->GetBinError(x+1)   );     } else { data_hh_m_err.push_back( 0 ); }
      if( h_hh_e_plus[samp]     !=0 ){ data_hh_e_plus     = h_hh_e_plus[samp]->GetBinContent(x+1);      } else { data_hh_e_plus = 0 ;          }
      if( h_hh_e_minus[samp]    !=0 ){ data_hh_e_minus    = h_hh_e_minus[samp]->GetBinContent(x+1) ;    } else { data_hh_e_minus= 0 ;          }  
      if( h_hh_m_plus[samp]     !=0 ){ data_hh_m_plus     = h_hh_m_plus[samp]->GetBinContent(x+1);      } else { data_hh_m_plus = 0 ;          }
      if( h_hh_m_minus[samp]    !=0 ){ data_hh_m_minus    = h_hh_m_minus[samp]->GetBinContent(x+1);     } else { data_hh_m_minus= 0 ;          }
      if( h_hh_e_uncplus[samp]  !=0 ){ data_hh_e_uncplus  = h_hh_e_uncplus[samp]->GetBinContent(x+1);   } else { data_hh_e_uncplus = 0 ;       }
      if( h_hh_e_uncminus[samp] !=0 ){ data_hh_e_uncminus = h_hh_e_uncminus[samp]->GetBinContent(x+1) ; } else { data_hh_e_uncminus= 0 ;       } 
      if( h_hh_m_uncplus[samp]  !=0 ){ data_hh_m_uncplus  = h_hh_m_uncplus[samp]->GetBinContent(x+1);   } else { data_hh_m_uncplus = 0 ;       }
      if( h_hh_m_uncminus[samp] !=0 ){ data_hh_m_uncminus = h_hh_m_uncminus[samp]->GetBinContent(x+1);  } else { data_hh_m_uncminus= 0 ;       }
      if( h_hh_e_jerplus[samp]  !=0 ){ data_hh_e_jerplus  = h_hh_e_jerplus[samp]->GetBinContent(x+1);   } else { data_hh_e_jerplus = 0 ;       }
      if( h_hh_e_jerminus[samp] !=0 ){ data_hh_e_jerminus = h_hh_e_jerminus[samp]->GetBinContent(x+1) ; } else { data_hh_e_jerminus= 0 ;       } 
      if( h_hh_m_jerplus[samp]  !=0 ){ data_hh_m_jerplus  = h_hh_m_jerplus[samp]->GetBinContent(x+1);   } else { data_hh_m_jerplus = 0 ;       }
      if( h_hh_m_jerminus[samp] !=0 ){ data_hh_m_jerminus = h_hh_m_jerminus[samp]->GetBinContent(x+1);  } else { data_hh_m_jerminus= 0 ;       }
      if( h_hh_e_bplus[samp]    !=0 ){ data_hh_e_bplus    = h_hh_e_bplus[samp]->GetBinContent(x+1);     } else { data_hh_e_bplus = 0 ;         }
      if( h_hh_e_bminus[samp]   !=0 ){ data_hh_e_bminus   = h_hh_e_bminus[samp]->GetBinContent(x+1);    } else { data_hh_e_bminus= 0 ;         }  
      if( h_hh_m_bplus[samp]    !=0 ){ data_hh_m_bplus    = h_hh_m_bplus[samp]->GetBinContent(x+1);     } else { data_hh_m_bplus = 0 ;         }
      if( h_hh_m_bminus[samp]   !=0 ){ data_hh_m_bminus   = h_hh_m_bminus[samp]->GetBinContent(x+1);    } else { data_hh_m_bminus= 0 ;         }
      if(triggerunc ){
        if( h_hh_e_trigger[samp]  !=0 ){ data_hh_e_trigger  = h_hh_e_trigger[samp]->GetBinContent(x+1);   } else { data_hh_e_trigger = 0;      }
      if( h_hh_m_trigger[samp]  !=0 ){ data_hh_m_trigger  = h_hh_m_trigger[samp]->GetBinContent(x+1);   } else { data_hh_m_trigger = 0;      }
    }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      double syst_plus_e     = fabs( data_hh_e[x]-data_hh_e_plus     );  double syst_plus_m     = fabs( data_hh_m[x]-data_hh_m_plus );
      double syst_minus_e    = fabs( data_hh_e[x]-data_hh_e_minus    );  double syst_minus_m    = fabs( data_hh_m[x]-data_hh_m_minus );
      double syst_uncplus_e  = fabs( data_hh_e[x]-data_hh_e_uncplus  );  double syst_uncplus_m  = fabs( data_hh_m[x]-data_hh_m_uncplus );
      double syst_uncminus_e = fabs( data_hh_e[x]-data_hh_e_uncminus );  double syst_uncminus_m = fabs( data_hh_m[x]-data_hh_m_uncminus );
      double syst_jerplus_e  = fabs( data_hh_e[x]-data_hh_e_jerplus  );  double syst_jerplus_m  = fabs( data_hh_m[x]-data_hh_m_jerplus );
      double syst_jerminus_e = fabs( data_hh_e[x]-data_hh_e_jerminus );  double syst_jerminus_m = fabs( data_hh_m[x]-data_hh_m_jerminus );
      double syst_bplus_e    = fabs( data_hh_e[x]-data_hh_e_bplus  );    double syst_bplus_m    = fabs( data_hh_m[x]-data_hh_m_bplus );
      double syst_bminus_e   = fabs( data_hh_e[x]-data_hh_e_bminus );    double syst_bminus_m   = fabs( data_hh_m[x]-data_hh_m_bminus );

      double syst_lum_e      = LUM_ERR*data_hh_e[x];                     double syst_lum_m      = LUM_ERR*data_hh_m[x];
      double syst_eff_e      = LEP_EFF_ERR*data_hh_e[x];                 double syst_eff_m      = LEP_EFF_ERR*data_hh_m[x];
      double syst_ttbar_cs_e = TTBAR_CS_ERR*data_hh_e[x];                double syst_ttbar_cs_m = TTBAR_CS_ERR*data_hh_m[x];  
      double syst_trigger_e  = data_hh_e_trigger;                        double syst_trigger_m  = data_hh_m_trigger;
      double syst_tauid_e(0);  double syst_tauid_m(0);

      if( x >= TAUSTEP){ syst_tauid_e = TAU_ID_ERR*data_hh_e[x]; syst_tauid_m = TAU_ID_ERR*data_hh_m[x]; } 
      
      

      if(systset1){ 
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset2){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;            
      }
      if(systset3){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
      }


      double temp_e = pow(syst_tauid_e,2) +  pow(syst_eff_e,2) + pow(syst_ttbar_cs_e,2) +pow(syst_lum_e,2);
      double temp_m = pow(syst_tauid_m,2) +  pow(syst_eff_m,2) + pow(syst_ttbar_cs_m,2) +pow(syst_lum_m,2);

      if(triggerunc){ temp_e += pow(syst_trigger_e,2 ); temp_m += pow(syst_trigger_m,2 ); }

      data_syst_e_plus.push_back(  sqrt( pow(syst_plus_e,2)   + pow(syst_uncplus_e,2)  + pow(syst_jerplus_e,2)  + pow(syst_bplus_e,2)  + temp_e )); 
      data_syst_e_minus.push_back( sqrt( pow(syst_minus_e,2)  + pow(syst_uncminus_e,2) + pow(syst_jerminus_e,2) + pow(syst_bminus_e,2) + temp_e )); 
      data_syst_m_plus.push_back(  sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + temp_m )); 
      data_syst_m_minus.push_back( sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + temp_m ));

      // DEBUG
      //if(samp==1 && x == OS_STEP2 ){
      //  cout<<"\n value is "<<data_hh_m[x];
      //  cout<<"\n syst plus  = "<<(sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + temp_m ))<<endl;
      //  cout<<"\n syst minus = "<<(sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + temp_m ))<<endl;
      //  cout<<"\n detailed minus : "
      //     <<"\n syst_minus_m    : "<<syst_minus_m<<" minus :"<<data_hh_m_minus
      //     <<"\n syst_uncminus_m : "<<syst_uncminus_m<<" minus :"<<data_hh_m_uncminus
      //     <<"\n syst_jerminus_m : "<<syst_jerminus_m
      //     <<"\n syst_bminus_m   : "<<syst_bminus_m
      //     <<"\n temp m          : "<<temp_m<<endl;
      //  cout<<"\n detailed plus  : "
      //     <<"\n syst_plus_m     : "<<syst_plus_m<<" plus : "<<data_hh_m_plus
      //     <<"\n syst_uncplus_m  : "<<syst_uncplus_m<<" plus :"<<data_hh_m_uncplus
      //     <<"\n syst_jerplus_m  : "<<syst_jerplus_m
      //     <<"\n syst_bplus_m    : "<<syst_bplus_m
      //     <<"\n temp m          : "<<temp_m<<endl;
      //}
     
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    }
    hh_e[samp] = data_hh_e; hh_m[samp] = data_hh_m; hh_e_err[samp] = data_hh_e_err; hh_m_err[samp] = data_hh_m_err;
    hh_e_syst_plus[samp] = data_syst_e_plus; hh_e_syst_minus[samp] = data_syst_e_minus; 
    hh_m_syst_plus[samp] = data_syst_m_plus; hh_m_syst_minus[samp] = data_syst_m_minus; 

  }
  for(uint samp=0; samp<h_hw_e.size(); samp++){


    vector<double> data_hw_e, data_hw_m, data_hw_e_err, data_hw_m_err,data_syst_e_plus, data_syst_e_minus, data_syst_m_plus, data_syst_m_minus;
    for(int x=0; x<14;x++){

      // get data and stats erros ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double data_hw_e_plus,    data_hw_e_minus,    data_hw_m_plus,    data_hw_m_minus;
      double data_hw_e_uncplus, data_hw_e_uncminus, data_hw_m_uncplus, data_hw_m_uncminus;
      double data_hw_e_jerplus, data_hw_e_jerminus, data_hw_m_jerplus, data_hw_m_jerminus;
      double data_hw_e_bplus,   data_hw_e_bminus,   data_hw_m_bplus,   data_hw_m_bminus;
      double data_hw_e_trigger, data_hw_m_trigger;

      if( h_hw_e[samp]          !=0 ){ data_hw_e.push_back(     h_hw_e[samp]->GetBinContent(x+1) );   } else { data_hw_e.push_back( 0 );    }
      if( h_hw_e[samp]          !=0 ){ data_hw_e_err.push_back( h_hw_e[samp]->GetBinError(x+1)   );   } else { data_hw_e_err.push_back( 0 );}
      if( h_hw_m[samp]          !=0 ){ data_hw_m.push_back(     h_hw_m[samp]->GetBinContent(x+1) );   } else { data_hw_m.push_back( 0 );    }
      if( h_hw_m[samp]          !=0 ){ data_hw_m_err.push_back( h_hw_m[samp]->GetBinError(x+1)   );   } else { data_hw_m_err.push_back( 0 );}
      if( h_hw_e_plus[samp]     !=0 ){ data_hw_e_plus  = h_hw_e_plus[samp]->GetBinContent(x+1);       } else { data_hw_e_plus = 0 ;         }
      if( h_hw_e_minus[samp]    !=0 ){ data_hw_e_minus = h_hw_e_minus[samp]->GetBinContent(x+1);      } else { data_hw_e_minus= 0 ;         }
      if( h_hw_m_plus[samp]     !=0 ){ data_hw_m_plus  = h_hw_m_plus[samp]->GetBinContent(x+1);       } else { data_hw_m_plus = 0 ;         }
      if( h_hw_m_minus[samp]    !=0 ){ data_hw_m_minus = h_hw_m_minus[samp]->GetBinContent(x+1);      } else { data_hw_m_minus= 0 ;         }
      if( h_hw_e_uncplus[samp]  !=0 ){ data_hw_e_uncplus  = h_hw_e_uncplus[samp]->GetBinContent(x+1); } else { data_hw_e_uncplus = 0 ;      }
      if( h_hw_e_uncminus[samp] !=0 ){ data_hw_e_uncminus = h_hw_e_uncminus[samp]->GetBinContent(x+1);} else { data_hw_e_uncminus= 0 ;      } 
      if( h_hw_m_uncplus[samp]  !=0 ){ data_hw_m_uncplus  = h_hw_m_uncplus[samp]->GetBinContent(x+1); } else { data_hw_m_uncplus = 0 ;      }
      if( h_hw_m_uncminus[samp] !=0 ){ data_hw_m_uncminus = h_hw_m_uncminus[samp]->GetBinContent(x+1);} else { data_hw_m_uncminus= 0 ;      }
      if( h_hw_e_jerplus[samp]  !=0 ){ data_hw_e_jerplus  = h_hw_e_jerplus[samp]->GetBinContent(x+1); } else { data_hw_e_jerplus = 0 ;      }
      if( h_hw_e_jerminus[samp] !=0 ){ data_hw_e_jerminus = h_hw_e_jerminus[samp]->GetBinContent(x+1);} else { data_hw_e_jerminus= 0 ;      } 
      if( h_hw_m_jerplus[samp]  !=0 ){ data_hw_m_jerplus  = h_hw_m_jerplus[samp]->GetBinContent(x+1); } else { data_hw_m_jerplus = 0 ;      }
      if( h_hw_m_jerminus[samp] !=0 ){ data_hw_m_jerminus = h_hw_m_jerminus[samp]->GetBinContent(x+1);} else { data_hw_m_jerminus= 0 ;      }
      if( h_hw_e_bplus[samp]    !=0 ){ data_hw_e_bplus    = h_hw_e_bplus[samp]->GetBinContent(x+1);   } else { data_hw_e_bplus = 0 ;        }
      if( h_hw_e_bminus[samp]   !=0 ){ data_hw_e_bminus   = h_hw_e_bminus[samp]->GetBinContent(x+1);  } else { data_hw_e_bminus= 0 ;        } 
      if( h_hw_m_bplus[samp]    !=0 ){ data_hw_m_bplus    = h_hw_m_bplus[samp]->GetBinContent(x+1);   } else { data_hw_m_bplus = 0 ;        }
      if( h_hw_m_bminus[samp]   !=0 ){ data_hw_m_bminus   = h_hw_m_bminus[samp]->GetBinContent(x+1);  } else { data_hw_m_bminus= 0 ;        }
      if(triggerunc ){
        if( h_hw_e_trigger[samp]  !=0 ){ data_hw_e_trigger  = h_hw_e_trigger[samp]->GetBinContent(x+1); } else { data_hw_e_trigger = 0;       }
        if( h_hw_m_trigger[samp]  !=0 ){ data_hw_m_trigger  = h_hw_m_trigger[samp]->GetBinContent(x+1); } else { data_hw_m_trigger = 0;       }    
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      double syst_plus_e  = fabs( data_hw_e[x]-data_hw_e_plus );       double syst_plus_m  = fabs( data_hw_m[x]-data_hw_m_plus );
      double syst_minus_e = fabs( data_hw_e[x]-data_hw_e_minus );      double syst_minus_m = fabs( data_hw_m[x]-data_hw_m_minus );

      double syst_uncplus_e  = fabs( data_hw_e[x]-data_hw_e_uncplus );  double syst_uncplus_m  = fabs( data_hw_m[x]-data_hw_m_uncplus );
      double syst_uncminus_e = fabs( data_hw_e[x]-data_hw_e_uncminus ); double syst_uncminus_m = fabs( data_hw_m[x]-data_hw_m_uncminus );

      double syst_jerplus_e  = fabs( data_hw_e[x]-data_hw_e_jerplus );  double syst_jerplus_m  = fabs( data_hw_m[x]-data_hw_m_jerplus );
      double syst_jerminus_e = fabs( data_hw_e[x]-data_hw_e_jerminus ); double syst_jerminus_m = fabs( data_hw_m[x]-data_hw_m_jerminus );

      double syst_bplus_e  = fabs( data_hw_e[x]-data_hw_e_bplus );      double syst_bplus_m  = fabs( data_hw_m[x]-data_hw_m_bplus );
      double syst_bminus_e = fabs( data_hw_e[x]-data_hw_e_bminus );     double syst_bminus_m = fabs( data_hw_m[x]-data_hw_m_bminus );

      double syst_lum_e      = LUM_ERR*data_hw_e[x];                    double syst_lum_m      = LUM_ERR*data_hw_m[x];
      double syst_eff_e      = LEP_EFF_ERR*data_hw_e[x];                double syst_eff_m      = LEP_EFF_ERR*data_hw_m[x];
      double syst_ttbar_cs_e = TTBAR_CS_ERR*data_hw_e[x];               double syst_ttbar_cs_m = TTBAR_CS_ERR*data_hw_m[x];  
      double syst_trigger_e  = data_hw_e_trigger;                       double syst_trigger_m  = data_hw_m_trigger;
      double syst_tauid_e(0);                                           double syst_tauid_m(0);

      if( x >= TAUSTEP){ syst_tauid_e = TAU_ID_ERR*data_hw_e[x]; syst_tauid_m = TAU_ID_ERR*data_hw_m[x]; }  

      if(systset1){ 
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset2){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset3){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
      }

      double temp_e = pow(syst_tauid_e,2) +  pow(syst_eff_e,2) + pow(syst_ttbar_cs_e,2) +pow(syst_lum_e,2);
      double temp_m = pow(syst_tauid_m,2) +  pow(syst_eff_m,2) + pow(syst_ttbar_cs_m,2) +pow(syst_lum_m,2);

      if(triggerunc){ temp_e += pow(syst_trigger_e,2 ); temp_m += pow(syst_trigger_m,2 ); }


      data_syst_e_plus.push_back(  sqrt( pow(syst_plus_e,2)   + pow(syst_uncplus_e,2)  + pow(syst_jerplus_e,2)  + pow(syst_bplus_e,2)  + temp_e )); 
      data_syst_e_minus.push_back( sqrt( pow(syst_minus_e,2)  + pow(syst_uncminus_e,2) + pow(syst_jerminus_e,2) + pow(syst_bminus_e,2) + temp_e )); 
      data_syst_m_plus.push_back(  sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + temp_m )); 
      data_syst_m_minus.push_back( sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + temp_m ));
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    }


    hw_e[samp] = data_hw_e; hw_m[samp] = data_hw_m; hw_e_err[samp] = data_hw_e_err; hw_m_err[samp] = data_hw_m_err;
    hw_e_syst_plus[samp] = data_syst_e_plus; hw_e_syst_minus[samp] = data_syst_e_minus; 
    hw_m_syst_plus[samp] = data_syst_m_plus; hw_m_syst_minus[samp] = data_syst_m_minus; 


  }

  // TBH
  map< int, vector<double> > tbh_e, tbh_m, tbh_e_err, tbh_m_err;
  map< int, vector<double> > tbh_e_syst_plus, tbh_e_syst_minus, tbh_m_syst_plus, tbh_m_syst_minus;


  for(uint samp=0; samp<h_tbh_e.size(); samp++){
    
    vector<double> data_tbh_e, data_tbh_m, data_tbh_e_err, data_tbh_m_err, data_syst_e_plus, data_syst_e_minus, data_syst_m_plus, data_syst_m_minus; 
    for(int x=0; x<14;x++){
      
      // get data and stats errors ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double data_tbh_e_plus,    data_tbh_e_minus,    data_tbh_m_plus,    data_tbh_m_minus;
      double data_tbh_e_uncplus, data_tbh_e_uncminus, data_tbh_m_uncplus, data_tbh_m_uncminus;
      double data_tbh_e_jerplus, data_tbh_e_jerminus, data_tbh_m_jerplus, data_tbh_m_jerminus;
      double data_tbh_e_bplus,   data_tbh_e_bminus,   data_tbh_m_bplus,   data_tbh_m_bminus;
      double data_tbh_e_trigger, data_tbh_m_trigger;
      
      if( h_tbh_e[samp]          !=0 ){ data_tbh_e.push_back(     h_tbh_e[samp]->GetBinContent(x+1) );     } else { data_tbh_e.push_back( 0 );     }
      if( h_tbh_e[samp]          !=0 ){ data_tbh_e_err.push_back( h_tbh_e[samp]->GetBinError(x+1)   );     } else { data_tbh_e_err.push_back( 0 ); }
      if( h_tbh_m[samp]          !=0 ){ data_tbh_m.push_back(     h_tbh_m[samp]->GetBinContent(x+1) );     } else { data_tbh_m.push_back( 0 );     }
      if( h_tbh_m[samp]          !=0 ){ data_tbh_m_err.push_back( h_tbh_m[samp]->GetBinError(x+1)   );     } else { data_tbh_m_err.push_back( 0 ); }
      if( h_tbh_e_plus[samp]     !=0 ){ data_tbh_e_plus     = h_tbh_e_plus[samp]->GetBinContent(x+1);      } else { data_tbh_e_plus = 0 ;          }
      if( h_tbh_e_minus[samp]    !=0 ){ data_tbh_e_minus    = h_tbh_e_minus[samp]->GetBinContent(x+1) ;    } else { data_tbh_e_minus= 0 ;          }  
      if( h_tbh_m_plus[samp]     !=0 ){ data_tbh_m_plus     = h_tbh_m_plus[samp]->GetBinContent(x+1);      } else { data_tbh_m_plus = 0 ;          }
      if( h_tbh_m_minus[samp]    !=0 ){ data_tbh_m_minus    = h_tbh_m_minus[samp]->GetBinContent(x+1);     } else { data_tbh_m_minus= 0 ;          }
      if( h_tbh_e_uncplus[samp]  !=0 ){ data_tbh_e_uncplus  = h_tbh_e_uncplus[samp]->GetBinContent(x+1);   } else { data_tbh_e_uncplus = 0 ;       }
      if( h_tbh_e_uncminus[samp] !=0 ){ data_tbh_e_uncminus = h_tbh_e_uncminus[samp]->GetBinContent(x+1) ; } else { data_tbh_e_uncminus= 0 ;       } 
      if( h_tbh_m_uncplus[samp]  !=0 ){ data_tbh_m_uncplus  = h_tbh_m_uncplus[samp]->GetBinContent(x+1);   } else { data_tbh_m_uncplus = 0 ;       }
      if( h_tbh_m_uncminus[samp] !=0 ){ data_tbh_m_uncminus = h_tbh_m_uncminus[samp]->GetBinContent(x+1);  } else { data_tbh_m_uncminus= 0 ;       }
      if( h_tbh_e_jerplus[samp]  !=0 ){ data_tbh_e_jerplus  = h_tbh_e_jerplus[samp]->GetBinContent(x+1);   } else { data_tbh_e_jerplus = 0 ;       }
      if( h_tbh_e_jerminus[samp] !=0 ){ data_tbh_e_jerminus = h_tbh_e_jerminus[samp]->GetBinContent(x+1) ; } else { data_tbh_e_jerminus= 0 ;       } 
      if( h_tbh_m_jerplus[samp]  !=0 ){ data_tbh_m_jerplus  = h_tbh_m_jerplus[samp]->GetBinContent(x+1);   } else { data_tbh_m_jerplus = 0 ;       }
      if( h_tbh_m_jerminus[samp] !=0 ){ data_tbh_m_jerminus = h_tbh_m_jerminus[samp]->GetBinContent(x+1);  } else { data_tbh_m_jerminus= 0 ;       }
      if( h_tbh_e_bplus[samp]    !=0 ){ data_tbh_e_bplus    = h_tbh_e_bplus[samp]->GetBinContent(x+1);     } else { data_tbh_e_bplus = 0 ;         }
      if( h_tbh_e_bminus[samp]   !=0 ){ data_tbh_e_bminus   = h_tbh_e_bminus[samp]->GetBinContent(x+1);    } else { data_tbh_e_bminus= 0 ;         }  
      if( h_tbh_m_bplus[samp]    !=0 ){ data_tbh_m_bplus    = h_tbh_m_bplus[samp]->GetBinContent(x+1);     } else { data_tbh_m_bplus = 0 ;         }
      if( h_tbh_m_bminus[samp]   !=0 ){ data_tbh_m_bminus   = h_tbh_m_bminus[samp]->GetBinContent(x+1);    } else { data_tbh_m_bminus= 0 ;         }
      if(triggerunc ){
        if( h_tbh_e_trigger[samp]  !=0 ){ data_tbh_e_trigger  = h_tbh_e_trigger[samp]->GetBinContent(x+1);   } else { data_tbh_e_trigger = 0;      }
      if( h_tbh_m_trigger[samp]  !=0 ){ data_tbh_m_trigger  = h_tbh_m_trigger[samp]->GetBinContent(x+1);   } else { data_tbh_m_trigger = 0;      }
    }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      double syst_plus_e     = fabs( data_tbh_e[x]-data_tbh_e_plus     );  double syst_plus_m     = fabs( data_tbh_m[x]-data_tbh_m_plus );
      double syst_minus_e    = fabs( data_tbh_e[x]-data_tbh_e_minus    );  double syst_minus_m    = fabs( data_tbh_m[x]-data_tbh_m_minus );
      double syst_uncplus_e  = fabs( data_tbh_e[x]-data_tbh_e_uncplus  );  double syst_uncplus_m  = fabs( data_tbh_m[x]-data_tbh_m_uncplus );
      double syst_uncminus_e = fabs( data_tbh_e[x]-data_tbh_e_uncminus );  double syst_uncminus_m = fabs( data_tbh_m[x]-data_tbh_m_uncminus );
      double syst_jerplus_e  = fabs( data_tbh_e[x]-data_tbh_e_jerplus  );  double syst_jerplus_m  = fabs( data_tbh_m[x]-data_tbh_m_jerplus );
      double syst_jerminus_e = fabs( data_tbh_e[x]-data_tbh_e_jerminus );  double syst_jerminus_m = fabs( data_tbh_m[x]-data_tbh_m_jerminus );
      double syst_bplus_e    = fabs( data_tbh_e[x]-data_tbh_e_bplus  );    double syst_bplus_m    = fabs( data_tbh_m[x]-data_tbh_m_bplus );
      double syst_bminus_e   = fabs( data_tbh_e[x]-data_tbh_e_bminus );    double syst_bminus_m   = fabs( data_tbh_m[x]-data_tbh_m_bminus );

      double syst_lum_e      = LUM_ERR*data_tbh_e[x];                     double syst_lum_m      = LUM_ERR*data_tbh_m[x];
      double syst_eff_e      = LEP_EFF_ERR*data_tbh_e[x];                 double syst_eff_m      = LEP_EFF_ERR*data_tbh_m[x];
      double syst_ttbar_cs_e = TTBAR_CS_ERR*data_tbh_e[x];                double syst_ttbar_cs_m = TTBAR_CS_ERR*data_tbh_m[x];  
      double syst_trigger_e  = data_tbh_e_trigger;                        double syst_trigger_m  = data_tbh_m_trigger;
      double syst_tauid_e(0);  double syst_tauid_m(0);

      if( x >= TAUSTEP){ syst_tauid_e = TAU_ID_ERR*data_tbh_e[x]; syst_tauid_m = TAU_ID_ERR*data_tbh_m[x]; } 

      cout << "DEBUG TBH ME " << "-----------------------------------------------------------------------__" << endl;
      
      cout << " syst_plus_e     : " << syst_plus_e     << ", syst_plus_m     : " << syst_plus_m     << endl;    
      cout << " syst_minus_e    : " << syst_minus_e    << ", syst_minus_m    : " << syst_minus_m    << endl;
      cout << " syst_uncplus_e  : " << syst_uncplus_e  << ", syst_uncplus_m  : " << syst_uncplus_m  << endl;
      cout << " syst_uncminus_e : " << syst_uncminus_e << ", syst_uncminus_m : " << syst_uncminus_m << endl;
      cout << " syst_jerplus_e  : " << syst_jerplus_e  << ", syst_jerplus_m  : " << syst_jerplus_m  << endl;
      cout << " syst_jerminus_e : " << syst_jerminus_e << ", syst_jerminus_m : " << syst_jerminus_m << endl;
      cout << " syst_bplus_e    : " << syst_bplus_e    << ", syst_bplus_m    : " << syst_bplus_m    << endl;
      cout << " syst_bminus_e   : " << syst_bminus_e   << ", syst_bminus_m   : " << syst_bminus_m   << endl;
      cout << " syst_lum_e      : " << syst_lum_e      << ", syst_lum_m      : " << syst_lum_m      << endl;
      cout << " syst_eff_e      : " << syst_eff_e      << ", syst_eff_m      : " << syst_eff_m      << endl;
      cout << " syst_ttbar_cs_e : " << syst_ttbar_cs_e << ", syst_ttbar_cs_m : " << syst_ttbar_cs_m << endl;  
      cout << " syst_trigger_e  : " << syst_trigger_e  << ", syst_trigger_m  : " << syst_trigger_m  << endl;
      cout << " syst_tauid_e    : " << syst_tauid_e    << ", syst_tauid_m    : " << syst_tauid_m    << endl;

      cout << " END DEBUG " << endl;
      
      

      if(systset1){ 
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset2){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;            
      }
      if(systset3){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
      }


      double temp_e = pow(syst_tauid_e,2) +  pow(syst_eff_e,2) + pow(syst_ttbar_cs_e,2) +pow(syst_lum_e,2);
      double temp_m = pow(syst_tauid_m,2) +  pow(syst_eff_m,2) + pow(syst_ttbar_cs_m,2) +pow(syst_lum_m,2);

      if(triggerunc){ temp_e += pow(syst_trigger_e,2 ); temp_m += pow(syst_trigger_m,2 ); }

      data_syst_e_plus.push_back(  sqrt( pow(syst_plus_e,2)   + pow(syst_uncplus_e,2)  + pow(syst_jerplus_e,2)  + pow(syst_bplus_e,2)  + temp_e )); 
      data_syst_e_minus.push_back( sqrt( pow(syst_minus_e,2)  + pow(syst_uncminus_e,2) + pow(syst_jerminus_e,2) + pow(syst_bminus_e,2) + temp_e )); 
      data_syst_m_plus.push_back(  sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + temp_m )); 
      data_syst_m_minus.push_back( sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + temp_m ));

      // DEBUG
      //if(samp==1 && x == OS_STEP2 ){
      //  cout<<"\n value is "<<data_tbh_m[x];
      //  cout<<"\n syst plus  = "<<(sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + temp_m ))<<endl;
      //  cout<<"\n syst minus = "<<(sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + temp_m ))<<endl;
      //  cout<<"\n detailed minus : "
      //     <<"\n syst_minus_m    : "<<syst_minus_m<<" minus :"<<data_tbh_m_minus
      //     <<"\n syst_uncminus_m : "<<syst_uncminus_m<<" minus :"<<data_tbh_m_uncminus
      //     <<"\n syst_jerminus_m : "<<syst_jerminus_m
      //     <<"\n syst_bminus_m   : "<<syst_bminus_m
      //     <<"\n temp m          : "<<temp_m<<endl;
      //  cout<<"\n detailed plus  : "
      //     <<"\n syst_plus_m     : "<<syst_plus_m<<" plus : "<<data_tbh_m_plus
      //     <<"\n syst_uncplus_m  : "<<syst_uncplus_m<<" plus :"<<data_tbh_m_uncplus
      //     <<"\n syst_jerplus_m  : "<<syst_jerplus_m
      //     <<"\n syst_bplus_m    : "<<syst_bplus_m
      //     <<"\n temp m          : "<<temp_m<<endl;
      //}
     
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    }
    tbh_e[samp] = data_tbh_e; tbh_m[samp] = data_tbh_m; tbh_e_err[samp] = data_tbh_e_err; tbh_m_err[samp] = data_tbh_m_err;
    tbh_e_syst_plus[samp] = data_syst_e_plus; tbh_e_syst_minus[samp] = data_syst_e_minus; 
    tbh_m_syst_plus[samp] = data_syst_m_plus; tbh_m_syst_minus[samp] = data_syst_m_minus; 

  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  

  // reset matrix /////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int x=0;x<14;x++){ 
    //hh_e[x]= 0; hh_m[x] = 0; hw_e[x]= 0; hw_m[x] = 0;
    dataTotal[x]=0; syst_error_plus[x]=0; syst_error_minus[x]=0; 
    dataTotal_e[x]=0; dataTotal_m[x]=0;

    d_e[x]=0; d_m[x]=0; d_e_systErrPlus[x]=0; d_e_systErrMinus[x]=0; d_m_systErrPlus[x]=0; d_m_systErrMinus[x]=0;

    for(int y=0;y<15;y++){ 
      d[x][y]=0; e[x][y]=0; d_trigger[x][y]=0; d_plus[x][y]=0; d_minus[x][y]=0;
      d_uncplus[x][y]=0; d_uncminus[x][y]=0; d_jerplus[x][y]=0; d_jerminus[x][y]=0;  d_bplus[x][y]=0; d_bminus[x][y]=0; 
      systErrPlus[x][y]=0; systErrMinus[x][y]=0; 
    } 
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //vector<double> fakes; vector<double> fakes2_stat; vector<double> fakes2_syst;
  //for(int fff=0;fff<100;fff++){ fakes.push_back(0); fakes2_stat.push_back(0); fakes2_syst.push_back(0);}


  // fill data and compute S/B and S/sqrt(B) ////////////////////////////////////////////////////////////////
  for(int x=0;x<14;x++){
    
    // fill data ////////////////////////////////////////
    if(hData)   dataTotal[x]   = hData->GetBinContent(x+1);
    if(hData_e) dataTotal_e[x] = hData_e->GetBinContent(x+1);
    if(hData_m) dataTotal_m[x] = hData_m->GetBinContent(x+1);
    if(h_e)     d_e[x]         = h_e->GetBinContent(x+1);
    if(h_m)     d_m[x]         = h_m->GetBinContent(x+1);
    //////////////////////////////////////////////////////

    // fill values with histogram data //////////////////////////////////////////////////////////////////////
    for(int y=0;y<13;y++){  
   
      // yields //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(h)          { d[x][y]          = ( h->GetBinContent(x+1,y+1)         ) ;  e[x][y] = ( h->GetBinError(x+1,y+1) );            } 
      if(h_plus)     { d_plus[x][y]     = ( h_plus->GetBinContent(x+1,y+1)    ) ;                                                    } 
      if(h_minus)    { d_minus[x][y]    = ( h_minus->GetBinContent(x+1,y+1)   ) ;                                                    }
      if(h_uncplus)  { d_uncplus[x][y]  = ( h_uncplus->GetBinContent(x+1,y+1) ) ;                                                    } 
      if(h_uncminus) { d_uncminus[x][y] = ( h_uncminus->GetBinContent(x+1,y+1)) ;                                                    }
      if(h_jerplus)  { d_jerplus[x][y]  = ( h_jerplus->GetBinContent(x+1,y+1) ) ;                                                    } 
      if(h_jerminus) { d_jerminus[x][y] = ( h_jerminus->GetBinContent(x+1,y+1)) ;                                                    }
      if(h_bplus)    { d_bplus[x][y]    = ( h_bplus->GetBinContent(x+1,y+1)   ) ;                                                    } 
      if(h_bminus)   { d_bminus[x][y]   = ( h_bminus->GetBinContent(x+1,y+1)  ) ;                                                    }
      if (triggerunc && h_trigger)  { d_trigger[x][y]  = ( h_trigger->GetBinContent(x+1,y+1) ) ;                                     }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      cout << "Name: " << hName << ", value: " << d[x][y] << ", stat error: " << e[x][y] <<endl;

      if( XSECMEASUREMENT && ( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH)  ){  
        double scalefactor( XSEC_EXP_ / 234. );  //note it should be 234
        d[x][y]          *= scalefactor ; 
        d_plus[x][y]     *= scalefactor ;
        d_minus[x][y]    *= scalefactor ;     
        d_uncplus[x][y]  *= scalefactor ; 
        d_uncminus[x][y] *= scalefactor ;
        d_jerplus[x][y]  *= scalefactor ;
        d_jerminus[x][y] *= scalefactor ;
        d_bplus[x][y]    *= scalefactor ;
        d_bminus[x][y]   *= scalefactor ;  
        d_trigger[x][y]  *= scalefactor ;    
      }

      //compute systematics ///////////////////////////////////////////////////////////////////////////////////////////////
      double syst_plus     = fabs(d[x][y]-d_plus[x][y]);    
      double syst_minus    = fabs(d[x][y]-d_minus[x][y]);   
      double syst_uncplus  = fabs(d[x][y]-d_uncplus[x][y]);    
      double syst_uncminus = fabs(d[x][y]-d_uncminus[x][y]);   
      double syst_jerplus  = fabs(d[x][y]-d_jerplus[x][y]);    
      double syst_jerminus = fabs(d[x][y]-d_jerminus[x][y]); 
      double syst_bplus    = fabs(d[x][y]-d_bplus[x][y]);    
      double syst_bminus   = fabs(d[x][y]-d_bminus[x][y]); 
      double syst_lum      = LUM_ERR*d[x][y];
      double syst_eff_lep  = LEP_EFF_ERR*d[x][y];
      double syst_trigger  = d_trigger[x][y];


      double sample_cs(0);

      if(      y==QCD_TDCH )                                                                                                               { sample_cs = QCD_CS_ERR ;      }
      else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
      else if( y==SINGLETOP_TDCH )                                                                                                         { sample_cs = SINGLETOP_CS_ERR; }              
      else if( y==WJETS_TDCH     )                                                                                                         { sample_cs = WJETS_CS_ERR;     }              
      else                                                                                                                                 { sample_cs = OTHER_CS_ERR;     }

      double syst_cs = sample_cs * d[x][y];
     
      double syst_tauid(0);                                             
 
      if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=LEPLEP_TDCH     ){ syst_tauid = TAU_ID_ERR*d[x][y];     }  
      else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==LEPLEP_TDCH )  ){ syst_tauid = TAU_MISSID_ERR*d[x][y]; }




      // signal theory unc.
      double syst_theoryeff(0); if( XSECMEASUREMENT && ( y==ETAU_TDCH || y==MUTAU_TDCH )) syst_theoryeff = SIGNAL_EFF_ERR*d[x][y];

      cout << "DEBUG SM ME for sample " << x << " -----------------------------------------------------------------------__" << endl;
      cout << "  syst_plus      : " <<   syst_plus     << endl;
      cout << "  syst_minus     : " <<   syst_minus    << endl;
      cout << "  syst_uncplus   : " <<   syst_uncplus  << endl;
      cout << "  syst_uncminus  : " <<   syst_uncminus << endl;
      cout << "  syst_jerplus   : " <<   syst_jerplus  << endl;
      cout << "  syst_jerminus  : " <<   syst_jerminus << endl;
      cout << "  syst_bplus     : " <<   syst_bplus    << endl;
      cout << "  syst_bminus    : " <<   syst_bminus   << endl;
      cout << "  syst_lum       : " <<   syst_lum      << endl;
      cout << "  syst_eff_lep   : " <<   syst_eff_lep  << endl;
      cout << "  syst_trigger   : " <<   syst_trigger  << endl;
      cout << "  syst_theoryeff : " <<   syst_theoryeff<< endl;
      cout << " END DEBUG " << endl;


      if(systset1){ 
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ;
        syst_bplus = 0;  syst_bminus = 0; syst_theoryeff=0;
        syst_trigger=0;
      }
      if(systset2){
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ; 
        syst_uncplus = 0; syst_uncminus= 0;  syst_jerplus = 0;  syst_jerminus = 0; syst_plus=0; syst_minus=0; syst_theoryeff=0;
        syst_trigger=0;
      }
      if(systset3){
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ;
        syst_bplus = 0;  syst_bminus = 0; syst_theoryeff=0;
        syst_uncplus = 0; syst_uncminus= 0;  syst_jerplus = 0;  syst_jerminus = 0; syst_plus=0; syst_minus=0; syst_theoryeff=0;
        
      }

     
      // GOOD ERRORS
      double temp;
      temp = pow(syst_lum,2)+ pow(syst_tauid,2) + pow(syst_cs,2) + pow(syst_eff_lep,2)+pow(syst_theoryeff,2);
      if(triggerunc){ temp += pow(syst_trigger,2); }
      systErrPlus[x][y]  = sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2)  + temp ); 
      systErrMinus[x][y] = sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp );
      
   
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

      //cout<<endl<<" x ( channel code is "<<x<<" , selection step is "<<y<<" )  value is : "<<(d[x][y])<<endl;

      //compute fakes...
      //if(y==ZJETS_TDCH||y==LEPJETS_TDCH|| y== QCD_TDCH|| y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH || y==WJETS_TDCH ){
      //   fakes[x] += d[x][y]; fakes2_stat[x] += pow(e[x][y],2); fakes2_syst[x] += temp;
      //}

    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  

    // etau and mutau SM signal 

    //cout<<"\n debug "<<x<<" , etau is : "<<(d[x][0])<<" mutau is : "<<(d[x][1])<<endl;

    double signal(d[x][0] + d[x][1]);
    double signal_plus(   d_plus[x][0]    + d_plus[x][1]   ), signal_minus(    d_minus[x][0]    + d_minus[x][1]); 
    double signal_uncplus(d_uncplus[x][0] + d_uncplus[x][1]), signal_uncminus( d_uncminus[x][0] + d_uncminus[x][1]);
    double signal_jerplus(d_jerplus[x][0] + d_jerplus[x][1]), signal_jerminus( d_jerminus[x][0] + d_jerminus[x][1]);
    double signal_bplus(  d_bplus[x][0]   + d_bplus[x][1]  ), signal_bminus(   d_bminus[x][0]   + d_bminus[x][1]);
    double signal_trigger = pow( d_trigger[x][0],2) + pow(d_trigger[x][1],2) ;
    double err_signal2 = (e[x][0])*(e[x][0]) + (e[x][1])*(e[x][1]) ;

    
    // standard model background 
    double bkg(0),bkg_plus(0), bkg_minus(0), bkg_uncplus(0), bkg_uncminus(0), bkg_jerplus(0), bkg_jerminus(0), bkg_bplus(0), bkg_bminus(0), err_bkg2(0), bkg_trigger(0); 
    for(int y2=2;y2<13;y2++){ 
      bkg += d[x][y2]; err_bkg2 += (e[x][y2])*(e[x][y2]);
      bkg_plus     += d_plus[x][y2];     bkg_minus   += d_minus[x][y2];  
      bkg_uncminus += d_uncminus[x][y2]; bkg_uncplus += d_uncplus[x][y2];       
      bkg_jerminus += d_jerminus[x][y2]; bkg_jerplus += d_jerplus[x][y2];
      bkg_bminus   += d_bminus[x][y2];   bkg_bplus   += d_bplus[x][y2];
      bkg_trigger  += pow(d_trigger[x][y2],2); 
    }
    

    // S/B and S/sqrtB
    double sob(0), err_sob(0), sosqrtb(0), err_sosqrtb(0);
    if(bkg){ 
      sob     = signal/bkg;       err_sob     = getErrorFraction(signal,bkg);         // WARNING this should include the errors
      sosqrtb = signal/sqrt(bkg); err_sosqrtb = getErrorFractionWithSqrt(signal,bkg); // WARNING this should include the errors
    } 
   

    //cout<<" step : "<<x<<" sob "<<sob<<" sqrt b "<<sosqrtb<<" signal "<<signal<<" bkg "<<bkg<<endl;

    d_plus[x][13]    = signal_plus  + bkg_plus;    d_uncplus[x][13]  = signal_uncplus  + bkg_uncplus;   d_jerplus[x][13]  = signal_jerplus  + bkg_jerplus; 
    d_bplus[x][13]   = signal_bplus + bkg_bplus;
    
    d_minus[x][13]   = signal_minus + bkg_minus;   d_uncminus[x][13] = signal_uncminus + bkg_uncminus;  d_jerminus[x][13] = signal_jerminus + bkg_jerminus;
    d_bminus[x][13]  = signal_bminus  + bkg_bminus;

    d_trigger[x][13] = sqrt( signal_trigger +  bkg_trigger );


    d[x][13]=signal+bkg; e[x][13] = sqrt( err_signal2 + err_bkg2 );
    // S/B
    d[x][14] = sob;      e[x][14] = err_sob;
    // S/Sqrt(B)
    d[x][15] = sosqrtb;  e[x][15] = err_sosqrtb;


    double temp_syst_plus      = fabs(d_plus[x][13]    - d[x][13]);   
    double temp_syst_minus     = fabs(d_minus[x][13]   - d[x][13]);
    double temp_syst_uncplus   = fabs(d_uncplus[x][13] - d[x][13]);   
    double temp_syst_uncminus  = fabs(d_uncminus[x][13]- d[x][13]);
    double temp_syst_jerplus   = fabs(d_jerplus[x][13] - d[x][13]);   
    double temp_syst_jerminus  = fabs(d_jerminus[x][13]- d[x][13]);
    double temp_syst_bplus     = fabs(d_bplus[x][13]   - d[x][13]);
    double temp_syst_bminus    = fabs(d_bminus[x][13]  - d[x][13]);

    double temp_syst_ttbar      = TTBAR_CS_ERR     * ( d[x][ETAU_TDCH]+d[x][MUTAU_TDCH]+d[x][LEPJETS_TDCH]+d[x][LEPLEP_TDCH]+d[x][TAUJETS_TDCH]+d[x][TAUTAU_TDCH]+d[x][ALLJETS_TDCH]);
    double temp_syst_wjets      = WJETS_CS_ERR     * ( d[x][WJETS_TDCH] );
    double temp_syst_singletop  = SINGLETOP_CS_ERR * ( d[x][SINGLETOP_TDCH]);
    double temp_syst_diboson    = OTHER_CS_ERR     * ( d[x][EWKDI_TDCH] );
    double temp_syst_zjets      = OTHER_CS_ERR     * ( d[x][ZJETS_TDCH] + d[x][ZTAUTAU_TDCH]);
    double temp_syst_qcd        = QCD_CS_ERR       * ( d[x][QCD_TDCH] );

    double temp_syst_tauid      = 0; if( x>=TAUSTEP ){ temp_syst_tauid     = TAU_ID_ERR     *(d[x][ETAU_TDCH]+d[x][MUTAU_TDCH]+d[x][TAUJETS_TDCH]+d[x][TAUTAU_TDCH]);}
    double temp_syst_misstauid  = 0; if( x>=TAUSTEP ){ temp_syst_misstauid = TAU_MISSID_ERR *(d[x][LEPJETS_TDCH]+d[x][LEPLEP_TDCH]+d[x][ALLJETS_TDCH]);              }

    double temp_syst_lum      = LUM_ERR*d[x][13];
    double temp_syst_eff_lep  = LEP_EFF_ERR*d[x][13];

    double syst_theoryeff = 0;
    double temp_syst_trigger  = d_trigger[x][13] ;

    if (XSECMEASUREMENT) syst_theoryeff = SIGNAL_EFF_ERR*signal;
     
    if(systset1){ 
      temp_syst_misstauid = 0;
      temp_syst_eff_lep   = 0;
      temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
      temp_syst_bplus =0; temp_syst_bminus = 0; 
      syst_theoryeff=0;
      temp_syst_trigger=0;
    }
    if(systset2){
      temp_syst_misstauid = 0;
      temp_syst_eff_lep   = 0;
      temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
      temp_syst_plus = 0; temp_syst_minus = 0;
      temp_syst_uncplus = 0; temp_syst_uncminus = 0;
      temp_syst_jerplus = 0; temp_syst_jerminus = 0;
      syst_theoryeff=0;
      temp_syst_trigger=0;
    }
   if(systset3){
    
     temp_syst_misstauid = 0;
     temp_syst_eff_lep   = 0;
     temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
     temp_syst_bplus =0; temp_syst_bminus = 0; 
     temp_syst_plus = 0; temp_syst_minus = 0;
     temp_syst_uncplus = 0; temp_syst_uncminus = 0;
     temp_syst_jerplus = 0; temp_syst_jerminus = 0;
     syst_theoryeff=0;
   }


   // Good errors
   double temp(0);
   temp = pow(temp_syst_lum,2) + pow(temp_syst_tauid,2)+pow(temp_syst_misstauid,2)+pow(temp_syst_eff_lep,2)+ 
   pow(temp_syst_wjets,2)+pow(temp_syst_qcd,2)+pow(temp_syst_ttbar,2)+pow(temp_syst_singletop,2)+pow(temp_syst_diboson,2)+pow(temp_syst_zjets,2) + pow(syst_theoryeff,2);

   //cout<<"\n debug : temp should be zero : "<<temp<<endl;

   if( triggerunc ){ temp += pow(temp_syst_trigger,2);  }

   syst_error_plus[x]  = sqrt( pow(temp_syst_plus,2)  + pow(temp_syst_uncplus,2) +  pow(temp_syst_jerplus,2)  + pow(temp_syst_bplus,2) + temp );
   syst_error_minus[x] = sqrt( pow(temp_syst_minus,2) + pow(temp_syst_uncminus,2) + pow(temp_syst_jerminus,2) + pow(temp_syst_bminus,2) + temp );

   cout<<"\n debug : trigger unc = "<<temp_syst_trigger<<" syst plus = "<<syst_error_plus[x]<<endl;
  

    //BULL SHIT //////////////////////////////////////////////////////////
    //temp = pow(temp_syst_lum,2) + pow(temp_syst_wjets,2)+pow(temp_syst_qcd,2)+pow(temp_syst_ttbar,2)+pow(temp_syst_singletop,2)+pow(temp_syst_diboson,2)+pow(temp_syst_zjets,2) + pow(syst_theoryeff,2);
    //syst_error_plus[x]  = sqrt( temp );
    //syst_error_minus[x] = sqrt( temp );
    //////////////////////////////////////////////////////////////////////

  }
  /////////////////////////////////////////////////////////////////////////////////////////////

  TString fName = name+key+TString(".tex");
  TString vName("-v2"), systName("");


  if(systset1)             systName=TString("-JESMETJER");
  if(systset2)             systName=TString("-BTAG");
  if(systset3)             systName=TString("-TRIGGER");
  if( (systset1 && systset2) || (systset1 && systset3) || (systset2 && systset3)  ) systName=TString("-ERROR"); // not valid condition
  
  
  
  if(!hh_contrib && !hw_contrib){ fName = name+key+TString("-SM")+systName+vName+TString(".tex");}
  if(hh_contrib && hw_contrib ) { fName = name+key+TString("-HIGGS")+systName+vName+TString(".tex");}
  if(tbh_contrib)               { fName = name+key+TString("-TBH")+systName+vName+TString(".tex");}

  FILE * f = fopen(fName.Data(),"w"); 

  if(incDocument)  prepareDocument(f);     
  if(incDocument)  fprintf(f,"\\begin{center} \n"); 
  if(printAllErrors) fprintf(f,"\\begin{footnotesize} \n");

 
  if(      STARTINGPOINT == STARTING_AT_LJETS_ ){        fprintf(f,"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c|} \n");  fprintf(f,"\\multicolumn{10}{c}{ } \\\\ \n"); }
  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ fprintf(f,"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| } \n");    fprintf(f,"\\multicolumn{9}{c}{ } \\\\ \n"); }
  
  fprintf(f,"\\hline \n");

 
  if( STARTINGPOINT == STARTING_AT_LJETS_ ){ fprintf(f,"\\multicolumn{1}{| c|}{ } & \\multicolumn{1}{ c|}{ channel } &  \\multicolumn{1}{ c|}{ $N_{lep}=1$ + $N_{jets}\\ge 3$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 40GeV$ }  & \\multicolumn{1}{ c|}{ $M_{T}\\ge 40 GeV/c^{2}$ } & \\multicolumn{1}{ c |}{ $\\ge$ 1btag } & \\multicolumn{1}{ c| }{ $N_{\\tau}=1$ } & \\multicolumn{1}{ c |}{ OS }  & \\multicolumn{1}{ c |}{ $\\ge$ 2btags } &  \\multicolumn{1}{ c|}{ $R \\ge 0.9$ } \\\\ \n"); }
  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ fprintf(f,"\\multicolumn{1}{| c|}{ } & \\multicolumn{1}{ c|}{ channel } & \\multicolumn{1}{ c|}{ $N_{lep}=1$ + $\\not\\!\\!E_T$+ $N_{jets}\\ge 3$ } & \\multicolumn{1}{ c|}{ $M_{T}\\ge 40 GeV/c^{2}$ } & \\multicolumn{1}{ c |}{ $\\ge$ 1btag } & \\multicolumn{1}{ c| }{ $N_{\\tau}=1$ } & \\multicolumn{1}{ c |}{ OS } & \\multicolumn{1}{ c |}{ $\\ge$ 2btags } &  \\multicolumn{1}{ c|}{ $R \\ge 0.9$ } \\\\ \n");}
  


  TString Tcolumn_mc[] = {
    "$e\\tau_h$","$\\mu\\tau_h$", 
    "lepton+jets", "dilepton", "tau+jets", "$\\tau\\tau$", "full had.", 
    "W+jets", "Z(ll)+jets", "$Z(\\tau\\tau)$+jets", "Single Top", "QCD", "WW,WZ,ZZ",
    "total ", "S/B", "S/$\\sqrt{B}$"
  };


  TString T_HW_mc[] = {
    "HW, $M_{H}=80 GeV/c^{2}$ &", "HW, $M_{H}=100 GeV/c^{2}$ &","HW, $M_{H}=120 GeV/c^{2}$ &", "HW, $M_{H}=140 GeV/c^{2}$ &",
    "HW, $M_{H}=150 GeV/c^{2}$ &",   // comment if needed
    "HW, $M_{H}=155 GeV/c^{2}$ & ",  // comment if needed
    "HW, $M_{H}=160 GeV/c^{2}$ &",
  };  

  TString T_HH_mc[]      = {
   "HH $M_{H}=80 GeV/c^{2}$ & ", "HH $M_{H}=100 GeV/c^{2}$ & ", "HH $M_{H}=120 GeV/c^{2}$ & ", "HH $M_{H}=140 GeV/c^{2}$ & ",
   "HH $M_{H}=150 GeV/c^{2}$ & ", // comment if needed 
   "HH $M_{H}=155 GeV/c^{2}$ & ", // comment if needed
   "HH $M_{H}=160 GeV/c^{2}$ & ",
  }; 

  TString T_TBH_mc[]      = {
   "TBH $M_{H}=180 GeV/c^{2}$ & ",
   "TBH $M_{H}=190 GeV/c^{2}$ & ",
   "TBH $M_{H}=200 GeV/c^{2}$ & ", 
   "TBH $M_{H}=220 GeV/c^{2}$ & ",
   "TBH $M_{H}=250 GeV/c^{2}$ & ", // comment if needed 
   "TBH $M_{H}=300 GeV/c^{2}$ & ", // comment if needed
  }; 

  fprintf(f,"\\hline \n");
  //fprintf(f,"\\hline \n");
 
  // size of extra samples
  int numbHWSamples = processedHWFiles.size();
  int numbHHSamples = processedHHFiles.size();
  int numbTBHSamples = processedTBHFiles.size();

  int extraLines_hw = 2*numbHWSamples;
  int extraLines_hh = 2*numbHHSamples;
  int extraLines_tbh = 2*numbTBHSamples;
  int extraLines    = extraLines_hw + extraLines_hh;
  if(tbh_contrib) extraLines = extraLines_tbh;
  int totalLines    = 17 + extraLines;


  for(int l=1;l<totalLines;l++){
    
    //cout<<endl<<" Debug => Line => "<<l<<endl;

    //hh samples
    if( l ==1 && numbHHSamples != 0 ){
      fprintf(f,"\n"); fprintf(f,"\\hline \n");

      for( int k=0; k<numbHHSamples ; k++){ 

        // Electron contribution ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        TString line(""), dtemp(""); line.Append(T_HH_mc[k]); line.Append(Tcolumn_mc[0]);  
        
        if( STARTINGPOINT == STARTING_AT_LJETS_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n"); 
          if(        detailed == D2 )  {dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");}
          else  if(  detailed == D3 )  {dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n");}
          line.Append(dtemp); const char * data = line.Data();
          fprintf(f,data,fhh*(hh_e[k])[COL2],fhh*(hh_e[k])[COL3],fhh*(hh_e[k])[COL4],fhh*(hh_e[k])[COL5],fhh*(hh_e[k])[COL6],fhh*(hh_e[k])[COL7],fhh*(hh_e[k])[COL8],fhh*(hh_e[k])[COL9]);
        }
        else if ( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n"); 
          if     ( detailed == D2 )  {dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n");  }
          else if( detailed == D3 )  {dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n");  }
          line.Append(dtemp); const char * data = line.Data();
          fprintf(f, data, fhh*(hh_e[k])[COL2] , fhh*(hh_e[k])[COL3] ,fhh*(hh_e[k])[COL4], fhh*(hh_e[k])[COL5], fhh*(hh_e[k])[COL6], fhh*(hh_e[k])[COL7], fhh*(hh_e[k])[COL8]);
        }

        if(printAllErrors){
                 
          TString lx(""); TString dx(""); const char * datax;
          
          if( STARTINGPOINT == STARTING_AT_LJETS_ ){
            if      ( detailed == D2 ){ dx= TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$\\\\ \n"); }
            else if ( detailed == D2 ){ dx= TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else        { dx= TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); }
           
            lx.Append(dx); datax= lx.Data();
          
            fprintf(f, datax, 
              fhh*hh_e_err[k][COL2], fhh*hh_e_syst_plus[k][COL2], fhh*hh_e_syst_minus[k][COL2],
              fhh*hh_e_err[k][COL3], fhh*hh_e_syst_plus[k][COL3], fhh*hh_e_syst_minus[k][COL3],    
              fhh*hh_e_err[k][COL4], fhh*hh_e_syst_plus[k][COL4], fhh*hh_e_syst_minus[k][COL4],
              fhh*hh_e_err[k][COL5], fhh*hh_e_syst_plus[k][COL5], fhh*hh_e_syst_minus[k][COL5],
              fhh*hh_e_err[k][COL6], fhh*hh_e_syst_plus[k][COL6], fhh*hh_e_syst_minus[k][COL6],
              fhh*hh_e_err[k][COL7], fhh*hh_e_syst_plus[k][COL7], fhh*hh_e_syst_minus[k][COL7],              
              fhh*hh_e_err[k][COL8], fhh*hh_e_syst_plus[k][COL8], fhh*hh_e_syst_minus[k][COL8],
              fhh*hh_e_err[k][COL9], fhh*hh_e_syst_plus[k][COL9], fhh*hh_e_syst_minus[k][COL9]); 
          } 
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
            if      ( detailed == D2 ){ dx = TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if ( detailed == D3 ){ dx = TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
           else          { dx = TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }

            lx.Append(dx); datax= lx.Data();
          
            fprintf(f, datax, 
              fhh*hh_e_err[k][COL2], fhh*hh_e_syst_plus[k][COL2], fhh*hh_e_syst_minus[k][COL2],
              fhh*hh_e_err[k][COL3], fhh*hh_e_syst_plus[k][COL3], fhh*hh_e_syst_minus[k][COL3],    
              fhh*hh_e_err[k][COL4], fhh*hh_e_syst_plus[k][COL4], fhh*hh_e_syst_minus[k][COL4],
              fhh*hh_e_err[k][COL5], fhh*hh_e_syst_plus[k][COL5], fhh*hh_e_syst_minus[k][COL5],
              fhh*hh_e_err[k][COL6], fhh*hh_e_syst_plus[k][COL6], fhh*hh_e_syst_minus[k][COL6],
              fhh*hh_e_err[k][COL7], fhh*hh_e_syst_plus[k][COL7], fhh*hh_e_syst_minus[k][COL7],
              fhh*hh_e_err[k][COL8], fhh*hh_e_syst_plus[k][COL8], fhh*hh_e_syst_minus[k][COL8]);
          } 
        }
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        

        // Muon contribution //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        TString line2(" & "); line2.Append(Tcolumn_mc[1]); line2.Append(dtemp); const char * data2 = line2.Data();
    
        if(      STARTINGPOINT == STARTING_AT_LJETS_        ){ fprintf(f, data2,  fhh*hh_m[k][COL2], fhh*hh_m[k][COL3], fhh*hh_m[k][COL4], fhh*hh_m[k][COL5], fhh*hh_m[k][COL6], fhh*hh_m[k][COL7],fhh*hh_m[k][COL8], fhh*hh_m[k][COL9] ); }
        else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ fprintf(f, data2,  fhh*hh_m[k][COL2], fhh*hh_m[k][COL3], fhh*hh_m[k][COL4], fhh*hh_m[k][COL5], fhh*hh_m[k][COL6], fhh*hh_m[k][COL7], fhh*hh_m[k][COL8]); }
        

        if(printAllErrors){

          TString lx(""); TString dx(""); const char * datax;

          if( STARTINGPOINT      == STARTING_AT_LJETS_        ){
            if     ( detailed == D2){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if( detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{         dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); }
        
            lx.Append(dx);  datax = lx.Data();
            fprintf(f, datax, 
              fhh*hh_m_err[k][COL2], fhh*hh_m_syst_plus[k][COL2], fhh*hh_m_syst_minus[k][COL2],
              fhh*hh_m_err[k][COL3], fhh*hh_m_syst_plus[k][COL3], fhh*hh_m_syst_minus[k][COL3],    
              fhh*hh_m_err[k][COL4], fhh*hh_m_syst_plus[k][COL4], fhh*hh_m_syst_minus[k][COL4], 
              fhh*hh_m_err[k][COL5], fhh*hh_m_syst_plus[k][COL5], fhh*hh_m_syst_minus[k][COL5],
              fhh*hh_m_err[k][COL6], fhh*hh_m_syst_plus[k][COL6], fhh*hh_m_syst_minus[k][COL6],
              fhh*hh_m_err[k][COL7], fhh*hh_m_syst_plus[k][COL7], fhh*hh_m_syst_minus[k][COL7],
              fhh*hh_m_err[k][COL8], fhh*hh_m_syst_plus[k][COL8], fhh*hh_m_syst_minus[k][COL8],
              fhh*hh_m_err[k][COL9], fhh*hh_m_syst_plus[k][COL9], fhh*hh_m_syst_minus[k][COL9]);
          }
          else if(STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){        
            if     (detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if(detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{           dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }
   
            lx.Append(dx);  datax = lx.Data();

            fprintf(f, datax, 
              fhh*hh_m_err[k][COL2], fhh*hh_m_syst_plus[k][COL2], fhh*hh_m_syst_minus[k][COL2],
              fhh*hh_m_err[k][COL3], fhh*hh_m_syst_plus[k][COL3], fhh*hh_m_syst_minus[k][COL3],    
              fhh*hh_m_err[k][COL4], fhh*hh_m_syst_plus[k][COL4], fhh*hh_m_syst_minus[k][COL4], 
              fhh*hh_m_err[k][COL5], fhh*hh_m_syst_plus[k][COL5], fhh*hh_m_syst_minus[k][COL5],
              fhh*hh_m_err[k][COL6], fhh*hh_m_syst_plus[k][COL6], fhh*hh_m_syst_minus[k][COL6], 
              fhh*hh_m_err[k][COL7], fhh*hh_m_syst_plus[k][COL7], fhh*hh_m_syst_minus[k][COL7],
              fhh*hh_m_err[k][COL8], fhh*hh_m_syst_plus[k][COL8], fhh*hh_m_syst_minus[k][COL8]);
           }

           //DEBUG 
	  //          if(k==1){ // OS step
	  cout<<"\n DEBUG HH : k="<<k<<" ->"<<(hh_m[k][COL5])<<" stat error after OS +-"<<(hh_m_err[k][COL5])<<" syst + "<<(hh_m_syst_plus[k][COL5])<<" syst - "<<(hh_m_syst_minus[k][COL5])<<endl;
	  //}
	  
         

         } 
         /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         fprintf(f,"\\hline \n\n"); 
      }
    }
           
    //hw samples
    if( l ==1 && numbHWSamples != 0 ){
      
      fprintf(f,"\\hline \n");
      for( int k=0; k<numbHWSamples ; k++){ 
        
        // Electron contribution ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        TString line(""); line.Append(T_HW_mc[k]); line.Append(Tcolumn_mc[0]);  
        TString dtemp;

        // old selection
        if( STARTINGPOINT == STARTING_AT_LJETS_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n"); 
          if     (detailed==D2 ){dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");  }
          else if(detailed==D3 ){dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n");  }
        }
        else if(  STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n");         
          if( detailed == D2 ){dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");         }
          if( detailed == D3 ){dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n");         }
        }

        line.Append(dtemp); const char * data = line.Data();

        if( STARTINGPOINT == STARTING_AT_LJETS_ ){
          fprintf(f, data, fhw*hw_e[k][COL2], fhw*hw_e[k][COL3], fhw*hw_e[k][COL4],fhw*hw_e[k][COL5], fhw*hw_e[k][COL6], fhw*hw_e[k][COL7], fhw*hw_e[k][COL8], fhw*hw_e[k][COL9]);
        }
        else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
          fprintf(f, data, fhw*hw_e[k][COL2], fhw*hw_e[k][COL3], fhw*hw_e[k][COL4],fhw*hw_e[k][COL5], fhw*hw_e[k][COL6], fhw*hw_e[k][COL7], fhw*hw_e[k][COL8]  );
        }


        if(printAllErrors){
          
          TString lx(""); TString dx(""); const char * datax;

          if( STARTINGPOINT == STARTING_AT_LJETS_ ){

            if     (detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if(detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  \\\\ \n"); }
            else{ dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$   & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$\\\\ \n"); }
            lx.Append(dx);  datax = lx.Data();
            fprintf(f, datax, 
              fhw*hw_e_err[k][COL2], fhw*hw_e_syst_plus[k][COL2], fhw*hw_e_syst_minus[k][COL2],
              fhw*hw_e_err[k][COL3], fhw*hw_e_syst_plus[k][COL3], fhw*hw_e_syst_minus[k][COL3],    
              fhw*hw_e_err[k][COL4], fhw*hw_e_syst_plus[k][COL4], fhw*hw_e_syst_minus[k][COL4],   
              fhw*hw_e_err[k][COL5], fhw*hw_e_syst_plus[k][COL5], fhw*hw_e_syst_minus[k][COL5],
              fhw*hw_e_err[k][COL6], fhw*hw_e_syst_plus[k][COL6], fhw*hw_e_syst_minus[k][COL6],
              fhw*hw_e_err[k][COL7], fhw*hw_e_syst_plus[k][COL7], fhw*hw_e_syst_minus[k][COL7],
              fhw*hw_e_err[k][COL8], fhw*hw_e_syst_plus[k][COL8], fhw*hw_e_syst_minus[k][COL8],
              fhw*hw_e_err[k][COL9], fhw*hw_e_syst_plus[k][COL9], fhw*hw_e_syst_minus[k][COL9]);

          }
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
    
            if     ( detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if( detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{ dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }

            lx.Append(dx);  datax = lx.Data();
            fprintf(f, datax, 
              fhw*hw_e_err[k][COL2], fhw*hw_e_syst_plus[k][COL2], fhw*hw_e_syst_minus[k][COL2],
              fhw*hw_e_err[k][COL3], fhw*hw_e_syst_plus[k][COL3], fhw*hw_e_syst_minus[k][COL3],    
              fhw*hw_e_err[k][COL4], fhw*hw_e_syst_plus[k][COL4], fhw*hw_e_syst_minus[k][COL4],   
              fhw*hw_e_err[k][COL5], fhw*hw_e_syst_plus[k][COL5], fhw*hw_e_syst_minus[k][COL5],
              fhw*hw_e_err[k][COL6], fhw*hw_e_syst_plus[k][COL6], fhw*hw_e_syst_minus[k][COL6],
              fhw*hw_e_err[k][COL7], fhw*hw_e_syst_plus[k][COL7], fhw*hw_e_syst_minus[k][COL7],
              fhw*hw_e_err[k][COL8], fhw*hw_e_syst_plus[k][COL8], fhw*hw_e_syst_minus[k][COL8]);           
         }
	  //	  cout<<"\n DEBUG : k="<<k<<" ->"<<(hw_m[k][COL5])<<" stat error after OS +-"<<(hw_m_err[k][COL5])<<" syst + "<<(hhwm_syst_plus[k][COL5])<<" syst - "<<(hw_m_syst_minus[k][COL5])<<endl;
       }

         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        // Muon contribution //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        TString line2(" & "); line2.Append(Tcolumn_mc[1]); line2.Append(dtemp); const char * data2 = line2.Data();
        
        if(      STARTINGPOINT == STARTING_AT_LJETS_ ){        fprintf(f, data2,  fhw*hw_m[k][COL2], fhw*hw_m[k][COL3], fhw*hw_m[k][COL4], fhw*hw_m[k][COL5], fhw*hw_m[k][COL6], fhw*hw_m[k][COL7],fhw*hw_m[k][COL8], fhw*hw_m[k][COL9]);}
        else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ fprintf(f, data2,  fhw*hw_m[k][COL2], fhw*hw_m[k][COL3], fhw*hw_m[k][COL4], fhw*hw_m[k][COL5], fhw*hw_m[k][COL6],fhw*hw_m[k][COL7], fhw*hw_m[k][COL8]); }

        if(printAllErrors){

          TString lx("");  TString dx(""); const char * datax ;
 
          if( STARTINGPOINT == STARTING_AT_LJETS_ ){
            if(detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$   \\\\ \n"); }
            else if(detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  \\\\ \n"); }  
            else{         dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); }
            
            lx.Append(dx);  datax = lx.Data();
            fprintf(f, datax, 
              fhw*hw_m_err[k][COL2], fhw*hw_m_syst_plus[k][COL2], fhw*hw_m_syst_minus[k][COL2],
              fhw*hw_m_err[k][COL3], fhw*hw_m_syst_plus[k][COL3], fhw*hw_m_syst_minus[k][COL3],    
              fhw*hw_m_err[k][COL4], fhw*hw_m_syst_plus[k][COL4], fhw*hw_m_syst_minus[k][COL4], 
              fhw*hw_m_err[k][COL5], fhw*hw_m_syst_plus[k][COL5], fhw*hw_m_syst_minus[k][COL5],
              fhw*hw_m_err[k][COL6], fhw*hw_m_syst_plus[k][COL6], fhw*hw_m_syst_minus[k][COL6],
              fhw*hw_m_err[k][COL7], fhw*hw_m_syst_plus[k][COL7], fhw*hw_m_syst_minus[k][COL7],
              fhw*hw_m_err[k][COL8], fhw*hw_m_syst_plus[k][COL8], fhw*hw_m_syst_minus[k][COL8],
              fhw*hw_m_err[k][COL9], fhw*hw_m_syst_plus[k][COL9], fhw*hw_m_syst_minus[k][COL9]);
          }
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
            if     ( detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if( detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{           dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }

            lx.Append(dx); datax = lx.Data();
            fprintf(f, datax, 
              fhw*hw_m_err[k][COL2], fhw*hw_m_syst_plus[k][COL2], fhw*hw_m_syst_minus[k][COL2],
              fhw*hw_m_err[k][COL3], fhw*hw_m_syst_plus[k][COL3], fhw*hw_m_syst_minus[k][COL3],    
              fhw*hw_m_err[k][COL4], fhw*hw_m_syst_plus[k][COL4], fhw*hw_m_syst_minus[k][COL4], 
              fhw*hw_m_err[k][COL5], fhw*hw_m_syst_plus[k][COL5], fhw*hw_m_syst_minus[k][COL5],
              fhw*hw_m_err[k][COL6], fhw*hw_m_syst_plus[k][COL6], fhw*hw_m_syst_minus[k][COL6],
              fhw*hw_m_err[k][COL7], fhw*hw_m_syst_plus[k][COL7], fhw*hw_m_syst_minus[k][COL7],
              fhw*hw_m_err[k][COL8], fhw*hw_m_syst_plus[k][COL8], fhw*hw_m_syst_minus[k][COL8]
            ); 
          }

	  cout<<"\n DEBUG HW : k="<<k<<" ->"<<(hw_m[k][COL5])<<" stat error after OS +-"<<(hw_m_err[k][COL5])<<" syst + "<<(hw_m_syst_plus[k][COL5])<<" syst - "<<(hw_m_syst_minus[k][COL5])<<endl;
        } 
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fprintf(f,"\\hline \n\n"); 
         
      }
   
    }
   
    //tbh samples
    if( l ==1 && numbTBHSamples != 0 ){
      fprintf(f,"\n"); fprintf(f,"\\hline \n");

      for( int k=0; k<numbTBHSamples ; k++){ 

        // Electron contribution ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        TString line(""), dtemp(""); line.Append(T_TBH_mc[k]); line.Append(Tcolumn_mc[0]);  
        
        if( STARTINGPOINT == STARTING_AT_LJETS_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n"); 
          if(        detailed == D2 )  {dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");}
          else  if(  detailed == D3 )  {dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n");}
          line.Append(dtemp); const char * data = line.Data();
          fprintf(f,data,brHtaunu_[k]*(tbh_e[k])[COL2],brHtaunu_[k]*(tbh_e[k])[COL3],brHtaunu_[k]*(tbh_e[k])[COL4],brHtaunu_[k]*(tbh_e[k])[COL5],brHtaunu_[k]*(tbh_e[k])[COL6],brHtaunu_[k]*(tbh_e[k])[COL7],brHtaunu_[k]*(tbh_e[k])[COL8],brHtaunu_[k]*(tbh_e[k])[COL9]);
        }
        else if ( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
          dtemp = TString(" & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n"); 
          if     ( detailed == D2 )  {dtemp = TString(" & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n");  }
          else if( detailed == D3 )  {dtemp = TString(" & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n");  }
          line.Append(dtemp); const char * data = line.Data();
          fprintf(f, data, brHtaunu_[k]*(tbh_e[k])[COL2] , brHtaunu_[k]*(tbh_e[k])[COL3] ,brHtaunu_[k]*(tbh_e[k])[COL4], brHtaunu_[k]*(tbh_e[k])[COL5], brHtaunu_[k]*(tbh_e[k])[COL6], brHtaunu_[k]*(tbh_e[k])[COL7], brHtaunu_[k]*(tbh_e[k])[COL8]);
        }

        if(printAllErrors){
                 
          TString lx(""); TString dx(""); const char * datax;
          
          if( STARTINGPOINT == STARTING_AT_LJETS_ ){
            if      ( detailed == D2 ){ dx= TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$\\\\ \n"); }
            else if ( detailed == D2 ){ dx= TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else        { dx= TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); }
           
            lx.Append(dx); datax= lx.Data();
          
            fprintf(f, datax, 
              brHtaunu_[k]*tbh_e_err[k][COL2], brHtaunu_[k]*tbh_e_syst_plus[k][COL2], brHtaunu_[k]*tbh_e_syst_minus[k][COL2],
              brHtaunu_[k]*tbh_e_err[k][COL3], brHtaunu_[k]*tbh_e_syst_plus[k][COL3], brHtaunu_[k]*tbh_e_syst_minus[k][COL3],    
              brHtaunu_[k]*tbh_e_err[k][COL4], brHtaunu_[k]*tbh_e_syst_plus[k][COL4], brHtaunu_[k]*tbh_e_syst_minus[k][COL4],
              brHtaunu_[k]*tbh_e_err[k][COL5], brHtaunu_[k]*tbh_e_syst_plus[k][COL5], brHtaunu_[k]*tbh_e_syst_minus[k][COL5],
              brHtaunu_[k]*tbh_e_err[k][COL6], brHtaunu_[k]*tbh_e_syst_plus[k][COL6], brHtaunu_[k]*tbh_e_syst_minus[k][COL6],
              brHtaunu_[k]*tbh_e_err[k][COL7], brHtaunu_[k]*tbh_e_syst_plus[k][COL7], brHtaunu_[k]*tbh_e_syst_minus[k][COL7],              
              brHtaunu_[k]*tbh_e_err[k][COL8], brHtaunu_[k]*tbh_e_syst_plus[k][COL8], brHtaunu_[k]*tbh_e_syst_minus[k][COL8],
              brHtaunu_[k]*tbh_e_err[k][COL9], brHtaunu_[k]*tbh_e_syst_plus[k][COL9], brHtaunu_[k]*tbh_e_syst_minus[k][COL9]); 
          } 
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
            if      ( detailed == D2 ){ dx = TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if ( detailed == D3 ){ dx = TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
           else          { dx = TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }

            lx.Append(dx); datax= lx.Data();
          
            fprintf(f, datax, 
              brHtaunu_[k]*tbh_e_err[k][COL2], brHtaunu_[k]*tbh_e_syst_plus[k][COL2], brHtaunu_[k]*tbh_e_syst_minus[k][COL2],
              brHtaunu_[k]*tbh_e_err[k][COL3], brHtaunu_[k]*tbh_e_syst_plus[k][COL3], brHtaunu_[k]*tbh_e_syst_minus[k][COL3],    
              brHtaunu_[k]*tbh_e_err[k][COL4], brHtaunu_[k]*tbh_e_syst_plus[k][COL4], brHtaunu_[k]*tbh_e_syst_minus[k][COL4],
              brHtaunu_[k]*tbh_e_err[k][COL5], brHtaunu_[k]*tbh_e_syst_plus[k][COL5], brHtaunu_[k]*tbh_e_syst_minus[k][COL5],
              brHtaunu_[k]*tbh_e_err[k][COL6], brHtaunu_[k]*tbh_e_syst_plus[k][COL6], brHtaunu_[k]*tbh_e_syst_minus[k][COL6],
              brHtaunu_[k]*tbh_e_err[k][COL7], brHtaunu_[k]*tbh_e_syst_plus[k][COL7], brHtaunu_[k]*tbh_e_syst_minus[k][COL7],
              brHtaunu_[k]*tbh_e_err[k][COL8], brHtaunu_[k]*tbh_e_syst_plus[k][COL8], brHtaunu_[k]*tbh_e_syst_minus[k][COL8]);
          } 
        }
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        

        // Muon contribution //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        TString line2(" & "); line2.Append(Tcolumn_mc[1]); line2.Append(dtemp); const char * data2 = line2.Data();
    
        if(      STARTINGPOINT == STARTING_AT_LJETS_        ){ fprintf(f, data2,  brHtaunu_[k]*tbh_m[k][COL2], brHtaunu_[k]*tbh_m[k][COL3], brHtaunu_[k]*tbh_m[k][COL4], brHtaunu_[k]*tbh_m[k][COL5], brHtaunu_[k]*tbh_m[k][COL6], brHtaunu_[k]*tbh_m[k][COL7],brHtaunu_[k]*tbh_m[k][COL8], brHtaunu_[k]*tbh_m[k][COL9] ); }
        else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ fprintf(f, data2,  brHtaunu_[k]*tbh_m[k][COL2], brHtaunu_[k]*tbh_m[k][COL3], brHtaunu_[k]*tbh_m[k][COL4], brHtaunu_[k]*tbh_m[k][COL5], brHtaunu_[k]*tbh_m[k][COL6], brHtaunu_[k]*tbh_m[k][COL7], brHtaunu_[k]*tbh_m[k][COL8]); }
        

        if(printAllErrors){

          TString lx(""); TString dx(""); const char * datax;

          if( STARTINGPOINT      == STARTING_AT_LJETS_        ){
            if     ( detailed == D2){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if( detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{         dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); }
        
            lx.Append(dx);  datax = lx.Data();
            fprintf(f, datax, 
              brHtaunu_[k]*tbh_m_err[k][COL2], brHtaunu_[k]*tbh_m_syst_plus[k][COL2], brHtaunu_[k]*tbh_m_syst_minus[k][COL2],
              brHtaunu_[k]*tbh_m_err[k][COL3], brHtaunu_[k]*tbh_m_syst_plus[k][COL3], brHtaunu_[k]*tbh_m_syst_minus[k][COL3],    
              brHtaunu_[k]*tbh_m_err[k][COL4], brHtaunu_[k]*tbh_m_syst_plus[k][COL4], brHtaunu_[k]*tbh_m_syst_minus[k][COL4], 
              brHtaunu_[k]*tbh_m_err[k][COL5], brHtaunu_[k]*tbh_m_syst_plus[k][COL5], brHtaunu_[k]*tbh_m_syst_minus[k][COL5],
              brHtaunu_[k]*tbh_m_err[k][COL6], brHtaunu_[k]*tbh_m_syst_plus[k][COL6], brHtaunu_[k]*tbh_m_syst_minus[k][COL6],
              brHtaunu_[k]*tbh_m_err[k][COL7], brHtaunu_[k]*tbh_m_syst_plus[k][COL7], brHtaunu_[k]*tbh_m_syst_minus[k][COL7],
              brHtaunu_[k]*tbh_m_err[k][COL8], brHtaunu_[k]*tbh_m_syst_plus[k][COL8], brHtaunu_[k]*tbh_m_syst_minus[k][COL8],
              brHtaunu_[k]*tbh_m_err[k][COL9], brHtaunu_[k]*tbh_m_syst_plus[k][COL9], brHtaunu_[k]*tbh_m_syst_minus[k][COL9]);
          }
          else if(STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){        
            if     (detailed == D2 ){ dx=TString(" & & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if(detailed == D3 ){ dx=TString(" & & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$  &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n"); }
            else{           dx=TString(" & & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n"); }
   
            lx.Append(dx);  datax = lx.Data();

            fprintf(f, datax, 
              brHtaunu_[k]*tbh_m_err[k][COL2], brHtaunu_[k]*tbh_m_syst_plus[k][COL2], brHtaunu_[k]*tbh_m_syst_minus[k][COL2],
              brHtaunu_[k]*tbh_m_err[k][COL3], brHtaunu_[k]*tbh_m_syst_plus[k][COL3], brHtaunu_[k]*tbh_m_syst_minus[k][COL3],    
              brHtaunu_[k]*tbh_m_err[k][COL4], brHtaunu_[k]*tbh_m_syst_plus[k][COL4], brHtaunu_[k]*tbh_m_syst_minus[k][COL4], 
              brHtaunu_[k]*tbh_m_err[k][COL5], brHtaunu_[k]*tbh_m_syst_plus[k][COL5], brHtaunu_[k]*tbh_m_syst_minus[k][COL5],
              brHtaunu_[k]*tbh_m_err[k][COL6], brHtaunu_[k]*tbh_m_syst_plus[k][COL6], brHtaunu_[k]*tbh_m_syst_minus[k][COL6], 
              brHtaunu_[k]*tbh_m_err[k][COL7], brHtaunu_[k]*tbh_m_syst_plus[k][COL7], brHtaunu_[k]*tbh_m_syst_minus[k][COL7],
              brHtaunu_[k]*tbh_m_err[k][COL8], brHtaunu_[k]*tbh_m_syst_plus[k][COL8], brHtaunu_[k]*tbh_m_syst_minus[k][COL8]);
           }

           //DEBUG 
	  //          if(k==1){ // OS step
	  cout<<"\n DEBUG TBH : k="<<k<<" ->"<<(tbh_m[k][COL5])<<" stat error after OS +-"<<(tbh_m_err[k][COL5])<<" syst + "<<(tbh_m_syst_plus[k][COL5])<<" syst - "<<(tbh_m_syst_minus[k][COL5])<<endl;
	  //}
	  
         

         } 
         /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         fprintf(f,"\\hline \n\n"); 
      }
    }
    // End tbh


    if(includeSM){
    // PRINT extra background sources /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( l ==1 ){
    
      if(hh_contrib || hw_contrib || tbh_contrib) fprintf(f,"\\hline \n\n"); 

      for(int k=0; k<13;k++){ // 13 sources of background

        //if(k== 3 || k ==2 || k==5 || k == 6 || k==7 || k==11 ) continue;

        TString d2;
        if     (k ==0 ) { d2=TString("tau dilepton from $t\\bar{t}$ & $e\\tau_h$ & ");}
        else if(k == 1) { d2=TString(" & $\\mu\\tau_h$ & ");}
        else if(k == 2) { fprintf(f,"\\hline ");  d2=TString(" & lepton+jets & ");}
        else if(k == 3) { d2=TString(" & dilepton & ");}
        else if(k == 4) { d2=TString(" other $t\\bar{t}$ & tau+jets & ");}
        else if(k == 5) { d2=TString(" & $\\tau\\tau$ & ");}
        else if(k == 6) { d2=TString(" & full had. & ");}
        else if(k == 7) { fprintf(f,"\\hline");  d2=TString(" & W+jets & ");}
        else if(k == 8) { d2=TString(" & Z(ll)+jets & ");}
        else if(k == 9) { d2=TString("  other bkg & Z($\\tau\\tau$)+jets & ");}
        else if(k == 10){ d2=TString(" & single top & ");}
        else if(k == 11){ d2=TString(" & QCD & "); }
        else if(k == 12){ d2=TString(" & WW,WZ,ZZ & ");}

        if      ( detailed == D1 ){ 
          if(      STARTINGPOINT == STARTING_AT_LJETS_        ){ d2 = d2 + TString(" %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n");}  
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ d2 = d2 + TString(" %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ \n"); }
        }else if( detailed == D2 ){ 
          if(      STARTINGPOINT == STARTING_AT_LJETS_        ){ d2 = d2 + TString(" %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");}  
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ d2 = d2 + TString(" %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f\\\\ \n"); }
        }
        else if( detailed == D3 ){ 
          if(      STARTINGPOINT == STARTING_AT_LJETS_        ){ d2 = d2 + TString(" %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n");}  
          else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ d2 = d2 + TString(" %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n"); }
        }

        TString line(""); line.Append(d2);  const char * data = line.Data();

        //if(k==1){ cout<<" debug... "<<" Mutau k=1 col2 is = "<<COL2<<" and value is "<<d[COL2][k]<<endl;   }       
         

        // old selection WHAT FOR ?
        //if(k==4){fprintf(f, data, fakes[COL2], fakes[COL3], fakes[COL4], fakes[COL5], fakes[COL6], fakes[COL7]);}
        //else{ fprintf(f, data, d[COL2][k], d[COL3][k], d[COL4][k],d[COL5][k], d[COL6][k], d[COL7][k]); }
        //
        if(      STARTINGPOINT == STARTING_AT_LJETS_         ){ fprintf(f, data, d[COL2][k], d[COL3][k], d[COL4][k],d[COL5][k], d[COL6][k], d[COL7][k], d[COL8][k], d[COL9][k]); }
        else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_  ){ fprintf(f, data, d[COL2][k], d[COL3][k], d[COL4][k],d[COL5][k], d[COL6][k], d[COL7][k], d[COL8][k]); }

        if(printAllErrors){
          TString lx(""); TString dx(""); const char * datax;
        

          //if(k!=4) dx =TString(" &  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); 
          //else {
          //  dx =TString(" &  & $\\pm$ %8.1f $\\pm$ %8.1f &  $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f  & $\\pm$ %8.1f $\\pm$ %8.1f }$  \\\\ \n"); 
          
          
          ////if(k!=4) dx =TString(" &  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$  \\\\ \n"); 
          ////else {
          ////  dx =TString(" &  & $\\pm$ %8.1f $\\pm$ %8.1f &  $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f & $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); 
          ////}
          

           ////if(k==4){  NOTA ISTO VEM DEPOIS ....
         
           //old selection
           // fprintf(f, datax, sqrt( fakes2_stat[COL2]) , sqrt(fakes2_syst[COL2]),
           //   sqrt( fakes2_stat[COL3]) , sqrt(fakes2_syst[COL3]),
           //   sqrt( fakes2_stat[COL4]) , sqrt(fakes2_syst[COL4]),
           //   sqrt( fakes2_stat[COL5]) , sqrt(fakes2_syst[COL5]),
           //   sqrt( fakes2_stat[COL6]) , sqrt(fakes2_syst[COL6]),
           //   sqrt( fakes2_stat[COL7]) , sqrt(fakes2_syst[COL7]));
           // BULL SHIT
           // fprintf(f, datax, sqrt( fakes2_stat[COL2]) , sqrt(fakes2_syst[COL2]),
           //   sqrt( fakes2_stat[COL3]) , sqrt(fakes2_syst[COL3]),
           //   sqrt( fakes2_stat[COL4]) , sqrt(fakes2_syst[COL4]),
           //   sqrt( fakes2_stat[COL5]) , sqrt(fakes2_syst[COL5]),
           //   sqrt( fakes2_stat[COL6]) , sqrt(fakes2_syst[COL6]));
          //}
          //else{
          // old selection
          //  fprintf(f, datax, e[COL2][k], systErrPlus[COL2][k], systErrMinus[COL2][k],
          //    e[COL3][k], systErrPlus[COL3][k], systErrMinus[COL3][k],   
          //    e[COL4][k], systErrPlus[COL4][k], systErrMinus[COL4][k], 
          //    e[COL5][k], systErrPlus[COL5][k], systErrMinus[COL5][k],
          //    e[COL6][k], systErrPlus[COL6][k], systErrMinus[COL6][k],
          //    e[COL7][k], systErrPlus[COL7][k], systErrMinus[COL7][k]);
           

          if(      STARTINGPOINT == STARTING_AT_LJETS_   ){
            if      ( detailed == D2 ){ dx =TString(" &  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ \\\\ \n"); }
            else if ( detailed == D1 ){         dx =TString(" &  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n");  }
            else if ( detailed == D3 ){         dx =TString(" &  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n");  }

            lx.Append(dx); datax = lx.Data();
            
            fprintf(f, datax, e[COL2][k], systErrPlus[COL2][k], systErrMinus[COL2][k],
                              e[COL3][k], systErrPlus[COL3][k], systErrMinus[COL3][k],   
                              e[COL4][k], systErrPlus[COL4][k], systErrMinus[COL4][k], 
                              e[COL5][k], systErrPlus[COL5][k], systErrMinus[COL5][k],
                              e[COL6][k], systErrPlus[COL6][k], systErrMinus[COL6][k],
                              e[COL7][k], systErrPlus[COL7][k], systErrMinus[COL7][k],                              
                              e[COL8][k], systErrPlus[COL8][k], systErrMinus[COL8][k],
                              e[COL9][k], systErrPlus[COL9][k], systErrMinus[COL9][k]);
          }
          
          if(      STARTINGPOINT == STARTING_AT_LJETSPLUSMET_   ){
            if     ( detailed == D2 ){ dx =TString(" &  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ &  $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$ & $\\pm$ %8.2f $^{+ %8.2f }_{- %8.2f }$  \\\\ \n"); }
            else if( detailed == D1 ){         dx =TString(" &  & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ &  $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ & $\\pm$ %8.1f $^{+ %8.1f }_{- %8.1f }$ \\\\ \n");  }
            else if( detailed == D3 ){         dx =TString(" &  & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ &  $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ & $\\pm$ %8.3f $^{+ %8.3f }_{- %8.3f }$ \\\\ \n");  }

            lx.Append(dx); datax = lx.Data();
            
            fprintf(f, datax, e[COL2][k], systErrPlus[COL2][k], systErrMinus[COL2][k],
                              e[COL3][k], systErrPlus[COL3][k], systErrMinus[COL3][k],   
                              e[COL4][k], systErrPlus[COL4][k], systErrMinus[COL4][k], 
                              e[COL5][k], systErrPlus[COL5][k], systErrMinus[COL5][k],
                              e[COL6][k], systErrPlus[COL6][k], systErrMinus[COL6][k],
                              e[COL7][k], systErrPlus[COL7][k], systErrMinus[COL7][k],
                              e[COL8][k], systErrPlus[COL8][k], systErrMinus[COL8][k]);}
            
          }
	cout<<"\n DEBUG SM: k="<<k<<" ->"<<(d[COL5][k])<<" stat error after OS +-"<<(e[k][COL5])<<" syst + "<<(systErrPlus[COL5][k])<<" syst - "<<(systErrMinus[COL5][k])<<endl;
      
        }         
      }
    //WARNING}
    // PRINT etau, mutau and total
    if( l ==1 ){
      fprintf(f,"\n \\hline");
      
      TString d2(""),d3(""),d4(""),d5(""),line_e(""),line_m(""), line(""),line4("");
      const char * data;  const char * data_e = line_e.Data(); const char * data_m = line_m.Data();

      if(      STARTINGPOINT == STARTING_AT_LJETS_   ){
        d2=TString("  & $e\\tau$  & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f  & %8.1f & %8.1f\\\\ \n");
        d3=TString("  & $\\mu\\tau$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n");
      }
      if(      STARTINGPOINT == STARTING_AT_LJETSPLUSMET_   ){
        d2=TString("  & $e\\tau$  & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ \n");
        d3=TString("  & $\\mu\\tau$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n");
      }
      line_e.Append(d2);     line_m.Append(d3);

      if( STARTINGPOINT == STARTING_AT_LJETS_ ){
        fprintf(f, data_e, d_e[COL2], d_e[COL3], d_e[COL4], d_e[COL5], d_e[COL6], d_e[COL7], d_e[COL8], d_e[COL9]);
        fprintf(f, data_m, d_m[COL2], d_m[COL3], d_m[COL4], d_m[COL5], d_m[COL6], d_m[COL7], d_m[COL8], d_m[COL9]);
        d4= TString(" MC(SM)  &  total & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f  & %8.0f & %8.0f  & %8.0f \\\\ \n");
        line4.Append(d4);  const char * data4 = line4.Data();
        fprintf(f, data4, d[COL2][13], d[COL3][13], d[COL4][13], d[COL5][13], d[COL6][13], d[COL7][13], d[COL8][13], d[COL9][13]);
        d5=TString(" & stat $\\pm$ syst & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ &  $\\pm$ %8.1f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ \\\\ \n"); 
        line.Append(d5);  data = line.Data();
        fprintf(f, data, e[COL2][13], syst_error_plus[COL2], syst_error_minus[COL2],
                         e[COL3][13], syst_error_plus[COL3], syst_error_minus[COL3],    
                         e[COL4][13], syst_error_plus[COL4], syst_error_minus[COL4], 
                         e[COL5][13], syst_error_plus[COL5], syst_error_minus[COL5],
                         e[COL6][13], syst_error_plus[COL6], syst_error_minus[COL6],
                         e[COL7][13], syst_error_plus[COL7], syst_error_minus[COL7],
                         e[COL8][13], syst_error_plus[COL8], syst_error_minus[COL8],
                         e[COL9][13], syst_error_plus[COL9], syst_error_minus[COL9]);
      } 


      if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
        fprintf(f, data_e, d_e[COL2], d_e[COL3], d_e[COL4], d_e[COL5], d_e[COL6], d_e[COL7], d_e[COL8]);
        fprintf(f, data_m, d_m[COL2], d_m[COL3], d_m[COL4], d_m[COL5], d_m[COL6], d_m[COL7], d_m[COL8]);
        d4= TString(" MC(SM)  &  total & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f\\\\ \n");
        line4.Append(d4);  const char * data4 = line4.Data();
        fprintf(f, data4, d[COL2][13], d[COL3][13], d[COL4][13], d[COL5][13], d[COL6][13], d[COL7][13], d[COL8][13]);
        d5=TString(" & stat $\\pm$ syst & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ &  $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ \\\\ \n"); 
        line.Append(d5);  data = line.Data();
        fprintf(f, data, e[COL2][13], syst_error_plus[COL2], syst_error_minus[COL2],
          e[COL3][13], syst_error_plus[COL3], syst_error_minus[COL3],    
          e[COL4][13], syst_error_plus[COL4], syst_error_minus[COL4], 
          e[COL5][13], syst_error_plus[COL5], syst_error_minus[COL5],
          e[COL6][13], syst_error_plus[COL6], syst_error_minus[COL6],
          e[COL7][13], syst_error_plus[COL7], syst_error_minus[COL7],
          e[COL8][13], syst_error_plus[COL8], syst_error_minus[COL8]
        );
      } 

    }
    if(l==1 && includeSoverB ){
      fprintf(f,"\n \\hline"); fprintf(f,"\n \\hline"); 
      TString d1(""),d2(""),line1(""),line2(""); const char * data1; const char * data2;
      
      if( STARTINGPOINT == STARTING_AT_LJETS_ ){
        d1=TString(" & S/B &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n");
        line1.Append(d1); data1 = line1.Data();
        fprintf(f, data1, d[COL2][14], d[COL3][14], d[COL4][14], d[COL5][14], d[COL6][14], d[COL7][14] , d[COL8][14], d[COL9][14]);

        d2=TString(" & S/$\\sqrt{B}$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");
        line2.Append(d2); data2 = line2.Data();
        fprintf(f, data2, d[COL2][15], d[COL3][15], d[COL4][15], d[COL5][15],  d[COL6][15],  d[COL7][15],  d[COL8][15],  d[COL9][15]);
      }
      if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
        d1=TString(" & S/B &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n");
        line1.Append(d1); data1 = line1.Data();
        fprintf(f, data1, d[COL2][14], d[COL3][14], d[COL4][14], d[COL5][14], d[COL6][14], d[COL7][14], d[COL8][14]);

        d2=TString(" & S/$\\sqrt{B}$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");
        line2.Append(d2); data2 = line2.Data();
        fprintf(f, data2, d[COL2][15], d[COL3][15], d[COL4][15], d[COL5][15],  d[COL6][15], d[COL7][15],  d[COL8][15]);
      }
    }  
  } 
  }//----->WARNING
  // SM backgrounds
  if(includeSM){
    if(!includeSoverB) fprintf(f,"\\hline \n"); 
    fprintf(f,"\\hline \n"); 
    TString lastLine_e(""),lastLine_m,lastLine_total(""); const char * dataBuffer_e; const char * dataBuffer_m; const char * dataBuffer_total;
   


    if( STARTINGPOINT == STARTING_AT_LJETS_ ){
      lastLine_e=TString("   & e$\\tau$ & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f \\\\ \n"); dataBuffer_e = lastLine_e.Data(); 
      fprintf(f,dataBuffer_e, dataTotal_e[COL2],dataTotal_e[COL3],dataTotal_e[COL4],dataTotal_e[COL5],dataTotal_e[COL6],dataTotal_m[COL7],dataTotal_e[COL8],dataTotal_m[COL9] );
      lastLine_m=TString(" data & $\\mu\\tau$ & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f \\\\ \n"); dataBuffer_m = lastLine_m.Data(); 
      fprintf(f,dataBuffer_m,  dataTotal_m[COL2],dataTotal_m[COL3],dataTotal_m[COL4],dataTotal_m[COL5],dataTotal_m[COL6],dataTotal_m[COL7],dataTotal_m[COL8],dataTotal_m[COL9]  );
      lastLine_total=TString(" & total & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f  & %8.0f & %8.0f  \\\\ \n"); dataBuffer_total = lastLine_total.Data(); 
      fprintf(f,dataBuffer_total, dataTotal[COL2],dataTotal[COL3],dataTotal[COL4],dataTotal[COL5],dataTotal[COL6],dataTotal[COL7],dataTotal[COL8],dataTotal[COL9] );
    }
    else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
      lastLine_e=TString("   & e$\\tau$ & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f  & %8.0f & %8.0f \\\\ \n"); dataBuffer_e = lastLine_e.Data(); 
      fprintf(f,dataBuffer_e, dataTotal_e[COL2],dataTotal_e[COL3],dataTotal_e[COL4],dataTotal_e[COL5],dataTotal_e[COL6],dataTotal_e[COL7],dataTotal_e[COL8]);
      lastLine_m=TString(" data & $\\mu\\tau$ & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f \\\\ \n"); dataBuffer_m = lastLine_m.Data(); 
      fprintf(f,dataBuffer_m,  dataTotal_m[COL2],dataTotal_m[COL3],dataTotal_m[COL4],dataTotal_m[COL5],dataTotal_m[COL6],dataTotal_m[COL7],dataTotal_m[COL8] );
      lastLine_total=TString(" & total & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f  & %8.0f & %8.0f  \\\\ \n"); dataBuffer_total = lastLine_total.Data(); 
      fprintf(f,dataBuffer_total, dataTotal[COL2],dataTotal[COL3],dataTotal[COL4],dataTotal[COL5],dataTotal[COL6] ,dataTotal[COL7],dataTotal[COL8]);
    }

    fprintf(f,"\\hline \n"); 

  }
  //fprintf(f,"\\end{tabular*}");
  fprintf(f,"\\end{tabular}");
  if(printAllErrors) fprintf(f,"\\end{footnotesize} \n");
  if(incDocument) fprintf(f,"\\end{center} \n");   //COMMENT
  if(incDocument) fprintf(f,"\\end{document} \n"); //COMMENT
  fclose (f);
  processedFile->Close();
  
    
}

















// Summary Table


















// systset1 -> if enabled include only JES + MET +JER
// systset2 -> if enabled include only btag unc
// systset3 -> if enabled include only trigger
  void TauDileptonTableBuilder::summaryTable( bool detailed, bool higgs, bool systset1, bool systset2, bool systset3, bool withShapes, bool withStatShapes){ 

    //    if(!doTheDatacards) // For datacards we want 234., for tables we want 1.1 // This does not need conditional, because of the default values
      for(size_t i=0; i<brHtaunu_.size(); ++i)
	brHtaunu_[i] = 5./234.;
      for(size_t i=0; i<brHtb_.size(); ++i)
	brHtb_[i] = 5./234.;
      
       
  //  double ttbar_init(20416081); //SIMPLE evt processed
  double ttbar_init(8228517); // evt processed from debug.txt
   
  
  // For datacards
  vector<vector<double> > tbh_datacards;
  vector<vector<double> > htb_datacards;
  vector<vector<double> > sm_datacards;
  vector<double>          data_datacards;
  vector<double>          taufakes_datacards;
  //MUON-TAU channel ( Bkg after OS :  89.3 +/- 9.1 +/- 10.9 )
  //double scale(0.09585155); double lum(2136);    double lumFactor(1);                                               // NO SCALE REWEIGHTING
  //double lum(2203.);   double lumFactor(1); double scalefactor(2203./2136.);  double scale(0.09585155*scalefactor); // WITH SCALE REWEIGHTING :: WARNING DATA SHOULD NOT BE SCALED!!!


                              
  // HIGGS ===>
  //double tau_fakes(89.3);     double tau_fakes_stat2(9.2*9.2); double tau_fakes_syst2(10.9*10.9);
  // XSEC  ===>
  //double tau_fakes(89.3);     double tau_fakes_stat2(9.1*9.1); double tau_fakes_syst2(10.9*10.9);   //no mt
  //double tau_fakes(88.5);     double tau_fakes_stat2(0); double tau_fakes_syst2(10.8*10.8+8.9*8.9);   //no mt, new single top xsec (all syst included)
  //double tau_fakes(88.5);     double tau_fakes_stat2(8.9*8.9); double tau_fakes_syst2(10.8*10.8);   //no mt, new single top xsec (all syst included)



  //ttbar_init = 930480; scale = 0.37877225;   //SCALE UP    VALUES
  //ttbar_init = 967053; scale = 0.36444745;   //SCALE DOWN  VALUES
  //ttbar_init = 1062788; scale = 0.33161835;  //PDF UP   
  //ttbar_init = 1065321; scale = 0.33082986;  //PDF DOWN
  //double tau_fakes(73.6);     double tau_fakes_stat2(8.1*8.1); double tau_fakes_syst2(8.8*8.8); //MT fakes


  //ELECTRON-TAU channel inclusive , Bkg after OS : 19.9 +/- 3.4 +/- 3.2
  //double scale(0.08494740); double lum(xxx);   double lumFactor(1);  

 
  //ELECTRON-TAU channel (electron+2jets+met trigger) ,  Bkg after OS : 54.8.3 +/- 6.5 +/- 8.2
  //double scale(0.08562379); double lum(1908.1); double lumFactor(1); //NO SCALE REWEIGHTING
  
  //double lum(1961.33); double lumFactor(1); double scalefactor(1961.33/1908.1); double scale(0.08562379*scalefactor); //WITH SCALE REWEIGHTING

  // FIXME: hardcoded
  //  double lum(4904); double lumFactor(1); double scalefactor(1.); double scale(0.02184544); //    double scale(0.03284470); // double scale(0.08562379); // NEW VALUES
  double lum(commondefinitions::LUM_); double lumFactor(1); double scalefactor(1.); double scale(0.51393073); //    double scale(0.03284470); // double scale(0.08562379); // NEW VALUES

  //  sqrt( (eA.B)^2 + (A.eB)^2 ), /B^2

  //  double tau_fakes(2939.9*0.699);     double tau_fakes_stat2( 0 ); double tau_fakes_syst2(261.6*261.6*0.699*0.699 + 0.02*0.02*2939.9*2939.9); // Preliminar from kNN

  //  $\pm$8.55099\% 

  double f_1tau(2762.19);
  double fe_1tau(263.891);

  double tau_fakes(f_1tau*0.699);     double tau_fakes_stat2( 0 ); double tau_fakes_syst2(fe_1tau*fe_1tau*0.699*0.699 + 0.02*0.02*f_1tau*f_1tau); // Preliminar from kNN

  taufakes_datacards.push_back(tau_fakes);
  taufakes_datacards.push_back(tau_fakes_stat2);
  taufakes_datacards.push_back(sqrt(tau_fakes_syst2));


  // tau fakes: otherttbarlepton+jets + otherbkgw+jets (+otherbckZlljets?)
  //scale = 0.0856237;
  //ttbar_init = 930480; scale = 0.33835569; //SCALE UP VALUES
  //ttbar_init = 967053; scale = 0.32555941; //SCALE DOWN  VALUES
  //ttbar_init = 1062788; scale = 0.29623330;  //SCALE PDF UP VALUES
  //ttbar_init = 1065321;  scale = 0.29552895 ;  //SCALE PDF SOWN VALUES


  //int STARTINGPOINT = STARTING_AT_LJETS_;
  int STARTINGPOINT = STARTING_AT_LJETSPLUSMET_;

  bool XSECMEASUREMENT(true);

  // WARNING : I have set the signal eff error to zero on not xsec measurement 
  // Note Lum error is set to zero this will be deriver afterwards
  if( XSECMEASUREMENT ){ TTBAR_CS_ERR = 0.; LUM_ERR=0; /*SIGNAL_EFF_ERR=0;*/ }   
  else                 { SIGNAL_EFF_ERR=0;                                   }

  TString version("-v2");


   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int STEP,TAUSTEP ,BTAG1STEP, BTAG2STEP; int COL1,COL2,COL3,COL4,COL5,COL6,COL7,COL8,COL9;
  
  TAUSTEP = TAU_STEP2; BTAG1STEP = BTAG1_STEP2; /*STEP=BTAG1_STEP2; */ STEP = OS_STEP2;
// nomt //  if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;    COL4=BTAG1_STEP2; COL5=TAU_STEP2; COL6=OS_STEP2;    COL7=BTAG2_STEP2; COL8=R_STEP1; }
// nomt //  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;   COL3=BTAG1_STEP2 ;  COL4=TAU_STEP2;   COL5=OS_STEP2;  COL6=BTAG2_STEP2; COL7=R_STEP2; }                 

  if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;   COL4=MT_STEP2;  COL5=BTAG1_STEP2; COL6=TAU_STEP2; COL7=OS_STEP2;    COL8=BTAG2_STEP2; COL9=R_STEP1; }
  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;  COL3=MT_STEP2; COL4=BTAG1_STEP2 ;  COL5=TAU_STEP2;   COL6=OS_STEP2;  COL7=BTAG2_STEP2; COL8=R_STEP2; }                 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // table options /////////////////////////////////////////////////////////////////
  bool incDocument(true);
  bool showAllHiggs(true), processHH(false),processHW(false),processTBH(true),processHTB(true);
  if(!higgs){ showAllHiggs=false; processTBH=false; processHTB=false;}
  bool addXsection; if(XSECMEASUREMENT){ addXsection=true; } else addXsection=false;
  //////////////////////////////////////////////////////////////////////////////////
  
  // yields for xsection measurement ////////////////////////////////////////////////////////////////////////////////
  double data_yields(0);   
  double acceptance(0.);      double acceptance_stat(0.);      double acceptance_syst(0.);                           
  double other_bkg(0);        double other_bkg_stat2(0);       double other_bkg_syst2(0);   
  double tau_dilepton(0);     double tau_dilepton_stat(0);     double tau_dilepton_syst(0);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   
  vector<TString> keys; TString fName; keys.push_back("PFlow");
  TString vName("-v2");
  TString systName("");
  
  if( systset1 )             systName=TString("-JESMETJER");
  if( systset2 )             systName=TString("-BTAG");
  if( systset1 && systset2 ) systName=TString("-ERROR");
  
  
  if    (!addXsection){ fName=TString("summary-tab-")+keys[0]+ systName + vName +TString(".tex");}
  else if(addXsection){ fName=TString("summary-tab-xsec-")+keys[0]+ systName + vName +TString(".tex");}
  
  TFile * processedMCFile = TFile::Open(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-mc.root")); TFile * processedDataFile = TFile::Open(outputArea_+TString("nomt-2012-V1-data-MU-20GeV/out-data.root")); 

  map<TString,TString> hNameData_m, hNameMC_m, 
    hNameMC_plus_m, hNameMC_minus_m, 
    hNameMC_uncplus_m, hNameMC_uncminus_m, 
    hNameMC_jerplus_m, hNameMC_jerminus_m,
    hNameMC_bplus_m, hNameMC_bminus_m,
    hNameMC_trigger_m;
  
  map<TString,TString> hName_higgs_m, hName_higgs_plus_m, hName_higgs_minus_m, 
    hName_higgs_uncplus_m, hName_higgs_uncminus_m, 
    hName_higgs_jerminus_m, hName_higgs_jerplus_m,
    hName_higgs_bminus_m, hName_higgs_bplus_m, hName_higgs_trigger_m;
  
  

  TString leptauString,etauString,mtauString, leptauStringMC;
  
  
  leptauStringMC=TString("1");   leptauString=TString("1"); etauString=TString("2"); mtauString=TString("3"); 
  
  for(uint i=0; i<keys.size(); i++){
    
    TString algo=keys[i];
    
    TString hNameData( algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_") + leptauString); hNameData_m[algo]            = hNameData;
    
    cout<<endl<<" DATA NAME : "<<hNameData<<endl;
    
    TString hNameMC(         algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_")         + leptauStringMC); hNameMC_m[algo]           = hNameMC;    
    TString hNameMC_plus(    algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_plus_")    + leptauStringMC); hNameMC_plus_m[algo]      = hNameMC_plus;
    TString hNameMC_minus(   algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_minus_")   + leptauStringMC); hNameMC_minus_m[algo]     = hNameMC_minus;
    TString hNameMC_uncplus( algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_uncplus_") + leptauStringMC); hNameMC_uncplus_m[algo]   = hNameMC_uncplus;
    TString hNameMC_uncminus(algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_uncminus_")+ leptauStringMC); hNameMC_uncminus_m[algo]  = hNameMC_uncminus;
    TString hNameMC_jerplus( algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_jerplus_") + leptauStringMC); hNameMC_jerplus_m[algo]   = hNameMC_jerplus;
    TString hNameMC_jerminus(algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_jerminus_")+ leptauStringMC); hNameMC_jerminus_m[algo]  = hNameMC_jerminus;
    TString hNameMC_bplus(   algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_bplus_")   + leptauStringMC); hNameMC_bplus_m[algo]     = hNameMC_bplus;
    TString hNameMC_bminus(  algo+TString("/YieldsMC")+TString("/cutflow_yields_mc_bminus_")  + leptauStringMC); hNameMC_bminus_m[algo]    = hNameMC_bminus;
    TString hNameMC_trigger( algo+TString("/YieldsMC")+TString("/cutflow_mc_triggErr_")       + leptauStringMC); hNameMC_trigger_m[algo]   = hNameMC_trigger;
     

    TString hName_higgs(          algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_")         +leptauString); hName_higgs_m[algo]          = hName_higgs;
    TString hName_higgs_plus(     algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_plus_")    +leptauString); hName_higgs_plus_m[algo]     = hName_higgs_plus;
    TString hName_higgs_minus(    algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_minus_")   +leptauString); hName_higgs_minus_m[algo]    = hName_higgs_minus;
    TString hName_higgs_uncplus(  algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_uncplus_") +leptauString); hName_higgs_uncplus_m[algo]  = hName_higgs_uncplus;
    TString hName_higgs_uncminus( algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_uncminus_")+leptauString); hName_higgs_uncminus_m[algo] = hName_higgs_uncminus;
    TString hName_higgs_jerplus(  algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_jerplus_") +leptauString); hName_higgs_jerplus_m[algo]  = hName_higgs_jerplus;
    TString hName_higgs_jerminus( algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_jerminus_")+leptauString); hName_higgs_jerminus_m[algo] = hName_higgs_jerminus;
    TString hName_higgs_bplus(    algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_bplus_")   +leptauString); hName_higgs_bplus_m[algo]    = hName_higgs_bplus;
    TString hName_higgs_bminus(   algo + TString("/Yields")+ TString("/lep_tau/cutflow_yields_bminus_")  +leptauString); hName_higgs_bminus_m[algo]   = hName_higgs_bminus;
    TString hName_higgs_trigger(  algo + TString("/Yields")+ TString("/lep_tau/cutflow_triggErr_")+leptauString); hName_higgs_trigger_m[algo]  = hName_higgs_trigger;
  
  }
  
  
  // map histograms for data
  map< TString, TH1D * > hData;
  // map histograms for MC 
  map< TString, TH2D * > hMC, hMC_plus, hMC_minus, hMC_uncplus, hMC_uncminus, hMC_jerplus, hMC_jerminus, hMC_bplus, hMC_bminus, hMC_trigger;      
  // map histograms for higgs
  map< TString, vector< TH1D * > > h_hh, h_hh_plus, h_hh_minus, h_hh_uncplus, h_hh_uncminus, h_hh_jerplus, h_hh_jerminus, h_hh_bplus, h_hh_bminus, h_hh_trigger;
  map< TString, vector< TH1D * > > h_hw, h_hw_plus, h_hw_minus, h_hw_uncplus, h_hw_uncminus, h_hw_jerplus, h_hw_jerminus, h_hw_bplus, h_hw_bminus, h_hw_trigger;
  map< TString, vector< TH1D * > > h_tbh, h_tbh_plus, h_tbh_minus, h_tbh_uncplus, h_tbh_uncminus, h_tbh_jerplus, h_tbh_jerminus, h_tbh_bplus, h_tbh_bminus, h_tbh_trigger;
  map< TString, vector< TH1D * > > h_htb, h_htb_plus, h_htb_minus, h_htb_uncplus, h_htb_uncminus, h_htb_jerplus, h_htb_jerminus, h_htb_bplus, h_htb_bminus, h_htb_trigger;
  

  // Get Data Histogram //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(processedDataFile){
    for(uint i=0; i<keys.size(); i++){ TString algo = keys[i];  hData[algo] = ( TH1D * ) processedDataFile->Get(hNameData_m[algo]);  (hData[algo])->Scale(scalefactor);}
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Get MC histograms ///////////////////////////////////////////////////////////
  if(processedMCFile){ 
    for(uint i=0; i<keys.size(); i++){
      TString algo = keys[i];
      hMC[algo]          = (TH2D*) processedMCFile->Get(hNameMC_m[algo]);                   hMC[algo]->Scale(scalefactor);
      hMC_plus[algo]     = (TH2D*) processedMCFile->Get(hNameMC_plus_m[algo]);              hMC_plus[algo]->Scale(scalefactor);
      hMC_minus[algo]    = (TH2D*) processedMCFile->Get(hNameMC_minus_m[algo]);             hMC_minus[algo]->Scale(scalefactor);
      hMC_uncplus[algo]  = (TH2D*) processedMCFile->Get(hNameMC_uncplus_m[algo]);           hMC_uncplus[algo]->Scale(scalefactor);
      hMC_uncminus[algo] = (TH2D*) processedMCFile->Get(hNameMC_uncminus_m[algo]);          hMC_uncminus[algo]->Scale(scalefactor);
      hMC_jerplus[algo]  = (TH2D*) processedMCFile->Get(hNameMC_jerplus_m[algo]);           hMC_jerplus[algo]->Scale(scalefactor);
      hMC_jerminus[algo] = (TH2D*) processedMCFile->Get(hNameMC_jerminus_m[algo]);          hMC_jerminus[algo]->Scale(scalefactor);
      hMC_bplus[algo]    = (TH2D*) processedMCFile->Get(hNameMC_bplus_m[algo]);             hMC_bplus[algo]->Scale(scalefactor);
      hMC_bminus[algo]   = (TH2D*) processedMCFile->Get(hNameMC_bminus_m[algo]);            hMC_bminus[algo]->Scale(scalefactor);
      hMC_trigger[algo]  = (TH2D*) processedMCFile->Get(hNameMC_trigger_m[algo]);           hMC_trigger[algo]->Scale(scalefactor);
    }
  } 
  ////////////////////////////////////////////////////////////////////////////////


  //get Higgs histograms////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<TString> hhFiles, hwFiles, tbhFiles, htbFiles;
  if( processHH ){
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m80.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m100.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m120.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m140.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m150.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m155.root"));
    hhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-hh-pythia-m160.root"));
  }
  if( processHW ){
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m80.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m100.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m120.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m140.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m150.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m155.root"));
    hwFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-wh-pythia-m160.root"));
  }
  if( processTBH ){
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m180.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m190.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m200.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m220.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m250.root"));
    tbhFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-tbh-pythia-m300.root"));
  }
  if( processHTB ){
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m180.root"));
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m190.root"));
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m200.root"));
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m220.root"));
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m250.root"));
    htbFiles.push_back(outputArea_+TString("nomt-2012-V1-mc-MU-20GeV/out-htb-pythia-m300.root"));
  }
  
 
  vector<TFile *> processedHWFiles, processedHHFiles, processedTBHFiles, processedHTBFiles;
  for(uint i = 0; i<hhFiles.size(); i++ ){  processedHHFiles.push_back( TFile::Open( hhFiles[i] ) ); }
  for(uint i = 0; i<hwFiles.size(); i++ ){  processedHWFiles.push_back( TFile::Open( hwFiles[i] ) ); }
  for(uint i = 0; i<tbhFiles.size(); i++ ){  processedTBHFiles.push_back( TFile::Open( tbhFiles[i] ) ); }
  for(uint i = 0; i<htbFiles.size(); i++ ){  processedHTBFiles.push_back( TFile::Open( htbFiles[i] ) ); }
  cout << "Open hh and wh and tbh and htb files if any" << endl;
  for(uint j=0; j<keys.size(); j++){
    vector<TH1D * > v_hh, v_hh_plus, v_hh_minus, v_hh_uncplus, v_hh_uncminus, v_hh_jerplus, v_hh_jerminus, v_hh_bplus, v_hh_bminus, v_hh_trigger;
    TString algo = keys[j];
    for(uint i = 0; i< hhFiles.size(); i++ ){ 
      if(processedHHFiles[i]){ 
        TH1D *h;     
	cout << "Now opening histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_m[algo]);          h->Scale(scalefactor); v_hh.push_back(h); 
	cout << "open first histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_plus_m[algo]);     h->Scale(scalefactor); v_hh_plus.push_back(h);
	cout << "open second histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_minus_m[algo]);    h->Scale(scalefactor); v_hh_minus.push_back(h);
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_uncplus_m[algo]);  h->Scale(scalefactor); v_hh_uncplus.push_back(h);
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_uncminus_m[algo]); h->Scale(scalefactor); v_hh_uncminus.push_back(h);
	cout << "open fifth histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_jerplus_m[algo]);  h->Scale(scalefactor); v_hh_jerplus.push_back(h);
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_jerminus_m[algo]); h->Scale(scalefactor); v_hh_jerminus.push_back(h);
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_bplus_m[algo]);    h->Scale(scalefactor); v_hh_bplus.push_back(h);
	cout << "open ninth histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_bminus_m[algo]);   h->Scale(scalefactor); v_hh_bminus.push_back(h);
	cout << "open tenth histo" << endl;
        h=(TH1D *)(processedHHFiles[i])->Get(hName_higgs_trigger_m[algo]);  h->Scale(scalefactor); v_hh_trigger.push_back(h);
	cout << "open last histo" << endl;
      }
      else { 
        v_hh.push_back(0);         v_hh_plus.push_back(0);       v_hh_minus.push_back(0); 
        v_hh_uncplus.push_back(0); v_hh_uncminus.push_back(0); 
        v_hh_jerplus.push_back(0); v_hh_jerminus.push_back(0);
        v_hh_bplus.push_back(0);   v_hh_bminus.push_back(0);
        v_hh_trigger.push_back(0);
      }
    }
    h_hh[algo] = v_hh;  h_hh_plus[algo] = v_hh_plus;  h_hh_minus[algo] = v_hh_minus;
    h_hh_uncplus[algo] = v_hh_uncplus;  h_hh_uncminus[algo] = v_hh_uncminus;
    h_hh_jerplus[algo] = v_hh_jerplus;  h_hh_jerminus[algo] = v_hh_jerminus;
    h_hh_bplus[algo]   = v_hh_bplus;    h_hh_bminus[algo] = v_hh_bminus;
    h_hh_trigger[algo] = v_hh_trigger;
  }

  cout << "Got hh histograms, if any" << endl;

  for(uint j=0; j<keys.size(); j++){
    vector<TH1D * > v_hw, v_hw_plus, v_hw_minus, v_hw_uncplus, v_hw_uncminus, v_hw_jerplus, v_hw_jerminus, v_hw_bplus, v_hw_bminus, v_hw_trigger;
    TString algo = keys[j];
    for(uint i = 0; i< hwFiles.size(); i++ ){ 
      if(processedHWFiles[i]){ 
        TH1D *h;     
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_m[algo]);          h->Scale(scalefactor); v_hw.push_back(h); 
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_plus_m[algo]);     h->Scale(scalefactor); v_hw_plus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_minus_m[algo]);    h->Scale(scalefactor); v_hw_minus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_uncplus_m[algo]);  h->Scale(scalefactor); v_hw_uncplus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_uncminus_m[algo]); h->Scale(scalefactor); v_hw_uncminus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_jerplus_m[algo]);  h->Scale(scalefactor); v_hw_jerplus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_jerminus_m[algo]); h->Scale(scalefactor); v_hw_jerminus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_bplus_m[algo]);    h->Scale(scalefactor); v_hw_bplus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_bminus_m[algo]);   h->Scale(scalefactor); v_hw_bminus.push_back(h);
        h=(TH1D *)(processedHWFiles[i])->Get(hName_higgs_trigger_m[algo]);  h->Scale(scalefactor); v_hw_trigger.push_back(h);

      }
      else { 
        v_hw.push_back(0);         v_hw_plus.push_back(0);     v_hw_minus.push_back(0); 
        v_hw_uncplus.push_back(0); v_hw_uncminus.push_back(0); 
        v_hw_jerplus.push_back(0); v_hw_jerminus.push_back(0); 
        v_hw_bplus.push_back(0);   v_hw_bminus.push_back(0);
        v_hw_trigger.push_back(0);  
      }
    }
    h_hw[algo] = v_hw;  h_hw_plus[algo] = v_hw_plus;  h_hw_minus[algo] = v_hw_minus;
    h_hw_uncplus[algo] = v_hw_uncplus;  h_hw_uncminus[algo] = v_hw_uncminus;
    h_hw_jerplus[algo] = v_hw_jerplus;  h_hw_jerminus[algo] = v_hw_jerminus;
    h_hw_bplus[algo]   = v_hw_bplus;    h_hw_bminus[algo] = v_hw_bminus;
    h_hw_trigger[algo] = v_hw_trigger;
  } 
  cout << "Got wh  histograms, if any" << endl;

  for(uint j=0; j<keys.size(); j++){
    vector<TH1D * > v_tbh, v_tbh_plus, v_tbh_minus, v_tbh_uncplus, v_tbh_uncminus, v_tbh_jerplus, v_tbh_jerminus, v_tbh_bplus, v_tbh_bminus, v_tbh_trigger;
    TString algo = keys[j];
    for(uint i = 0; i< tbhFiles.size(); i++ ){ 
      if(processedTBHFiles[i]){ 
        TH1D *h;     
	cout << "Now opening histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_m[algo]);          h->Scale(scalefactor); v_tbh.push_back(h); 
	cout << "open first histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_plus_m[algo]);     h->Scale(scalefactor); v_tbh_plus.push_back(h);
	cout << "open second histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_minus_m[algo]);    h->Scale(scalefactor); v_tbh_minus.push_back(h);
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_uncplus_m[algo]);  h->Scale(scalefactor); v_tbh_uncplus.push_back(h);
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_uncminus_m[algo]); h->Scale(scalefactor); v_tbh_uncminus.push_back(h);
	cout << "open fifth histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_jerplus_m[algo]);  h->Scale(scalefactor); v_tbh_jerplus.push_back(h);
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_jerminus_m[algo]); h->Scale(scalefactor); v_tbh_jerminus.push_back(h);
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_bplus_m[algo]);    h->Scale(scalefactor); v_tbh_bplus.push_back(h);
	cout << "open ninth histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_bminus_m[algo]);   h->Scale(scalefactor); v_tbh_bminus.push_back(h);
	cout << "open tenth histo" << endl;
        h=(TH1D *)(processedTBHFiles[i])->Get(hName_higgs_trigger_m[algo]);  h->Scale(scalefactor); v_tbh_trigger.push_back(h);
	cout << "open last histo" << endl;
      }
      else { 
        v_tbh.push_back(0);         v_tbh_plus.push_back(0);       v_tbh_minus.push_back(0); 
        v_tbh_uncplus.push_back(0); v_tbh_uncminus.push_back(0); 
        v_tbh_jerplus.push_back(0); v_tbh_jerminus.push_back(0);
        v_tbh_bplus.push_back(0);   v_tbh_bminus.push_back(0);
        v_tbh_trigger.push_back(0);
      }
    }
    h_tbh[algo] = v_tbh;  h_tbh_plus[algo] = v_tbh_plus;  h_tbh_minus[algo] = v_tbh_minus;
    h_tbh_uncplus[algo] = v_tbh_uncplus;  h_tbh_uncminus[algo] = v_tbh_uncminus;
    h_tbh_jerplus[algo] = v_tbh_jerplus;  h_tbh_jerminus[algo] = v_tbh_jerminus;
    h_tbh_bplus[algo]   = v_tbh_bplus;    h_tbh_bminus[algo] = v_tbh_bminus;
    h_tbh_trigger[algo] = v_tbh_trigger;
  }
  cout << "got TBH histograms, if any" << endl;

  for(uint j=0; j<keys.size(); j++){
    vector<TH1D * > v_htb, v_htb_plus, v_htb_minus, v_htb_uncplus, v_htb_uncminus, v_htb_jerplus, v_htb_jerminus, v_htb_bplus, v_htb_bminus, v_htb_trigger;
    TString algo = keys[j];
    for(uint i = 0; i< htbFiles.size(); i++ ){ 
      if(processedHTBFiles[i]){ 
        TH1D *h;     
	cout << "Now opening histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_m[algo]);          h->Scale(scalefactor); v_htb.push_back(h); 
	cout << "open first histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_plus_m[algo]);     h->Scale(scalefactor); v_htb_plus.push_back(h);
	cout << "open second histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_minus_m[algo]);    h->Scale(scalefactor); v_htb_minus.push_back(h);
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_uncplus_m[algo]);  h->Scale(scalefactor); v_htb_uncplus.push_back(h);
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_uncminus_m[algo]); h->Scale(scalefactor); v_htb_uncminus.push_back(h);
	cout << "open fifth histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_jerplus_m[algo]);  h->Scale(scalefactor); v_htb_jerplus.push_back(h);
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_jerminus_m[algo]); h->Scale(scalefactor); v_htb_jerminus.push_back(h);
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_bplus_m[algo]);    h->Scale(scalefactor); v_htb_bplus.push_back(h);
	cout << "open ninth histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_bminus_m[algo]);   h->Scale(scalefactor); v_htb_bminus.push_back(h);
	cout << "open tenth histo" << endl;
        h=(TH1D *)(processedHTBFiles[i])->Get(hName_higgs_trigger_m[algo]);  h->Scale(scalefactor); v_htb_trigger.push_back(h);
	cout << "open last histo" << endl;
      }
      else { 
        v_htb.push_back(0);         v_htb_plus.push_back(0);       v_htb_minus.push_back(0); 
        v_htb_uncplus.push_back(0); v_htb_uncminus.push_back(0); 
        v_htb_jerplus.push_back(0); v_htb_jerminus.push_back(0);
        v_htb_bplus.push_back(0);   v_htb_bminus.push_back(0);
        v_htb_trigger.push_back(0);
      }
    }
    h_htb[algo] = v_htb;  h_htb_plus[algo] = v_htb_plus;  h_htb_minus[algo] = v_htb_minus;
    h_htb_uncplus[algo] = v_htb_uncplus;  h_htb_uncminus[algo] = v_htb_uncminus;
    h_htb_jerplus[algo] = v_htb_jerplus;  h_htb_jerminus[algo] = v_htb_jerminus;
    h_htb_bplus[algo]   = v_htb_bplus;    h_htb_bminus[algo] = v_htb_bminus;
    h_htb_trigger[algo] = v_htb_trigger;
  }
  cout << "got HTB histograms, if any" << endl;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





  // HH / WH data/error holders
  map< pair< TString, int > , vector<double> > hh, hw, hh_err, hw_err;
  map< pair< TString, int > , vector<double> > hh_syst_plus, hh_syst_minus, hw_syst_plus, hw_syst_minus;

  map< pair< TString, int > , vector<double> > tbh, tbh_err;
  map< pair< TString, int > , vector<double> > tbh_syst_plus, tbh_syst_minus;

  map< pair< TString, int > , vector<double> > htb, htb_err;
  map< pair< TString, int > , vector<double> > htb_syst_plus, htb_syst_minus;

  // SM MC data/errors holders
  map< pair< TString, pair<int,int> > , double > dMC, dMC_err, dMC_plus, dMC_minus, dMC_uncplus, dMC_uncminus, dMC_trigger,
                                                 dMC_jerplus, dMC_jerminus, dMC_bplus, dMC_bminus,dMC_syst_plus, dMC_syst_minus;

  // DATA holder
  map< pair<TString,int> , double> dataTotal;
 

  //algorithm, sample de cada algo, vector de dados
  for(uint j=0; j<keys.size(); j++){

    TString algo = keys[j];

    // hh samples //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(uint samp=0; samp < h_hh[algo].size(); samp++){
      pair< TString, int> key(algo,samp);      
      vector<double> data_hh, data_hh_err, data_syst_plus, data_syst_minus; 
      for(int x=0; x<12;x++){
        // get data and stats errors //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double data_hh_plus, data_hh_minus;
        double data_hh_uncplus, data_hh_uncminus, data_hh_jerplus, data_hh_jerminus, data_hh_bplus, data_hh_bminus, data_hh_trigger;
        if( h_hh[algo][samp]          != 0 ){ data_hh.push_back(      h_hh[algo][samp]->GetBinContent(x+1) );  } else { data_hh.push_back( 0 );   }
        if( h_hh[algo][samp]          != 0 ){ data_hh_err.push_back(  h_hh[algo][samp]->GetBinError(x+1)   );  } else { data_hh_err.push_back(0); }
        if( h_hh_plus[algo][samp]     != 0 ){ data_hh_plus  =    h_hh_plus[algo][samp]->GetBinContent(x+1);    } else { data_hh_plus       = 0 ;  }
        if( h_hh_minus[algo][samp]    != 0 ){ data_hh_minus =    h_hh_minus[algo][samp]->GetBinContent(x+1);   } else { data_hh_minus      = 0 ;  }  
        if( h_hh_uncplus[algo][samp]  != 0 ){ data_hh_uncplus  = h_hh_uncplus[algo][samp]->GetBinContent(x+1); } else { data_hh_uncplus    = 0 ;  }
        if( h_hh_uncminus[algo][samp] != 0 ){ data_hh_uncminus = h_hh_uncminus[algo][samp]->GetBinContent(x+1);} else { data_hh_uncminus   = 0 ;  }  
        if( h_hh_jerplus[algo][samp]  != 0 ){ data_hh_jerplus  = h_hh_jerplus[algo][samp]->GetBinContent(x+1); } else { data_hh_jerplus    = 0 ;  }
        if( h_hh_jerminus[algo][samp] != 0 ){ data_hh_jerminus = h_hh_jerminus[algo][samp]->GetBinContent(x+1);} else { data_hh_jerminus   = 0 ;  }  
        if( h_hh_bplus[algo][samp]    != 0 ){ data_hh_bplus    = h_hh_bplus[algo][samp]->GetBinContent(x+1); }   else { data_hh_bplus      = 0 ;  }
        if( h_hh_bminus[algo][samp]   != 0 ){ data_hh_bminus   = h_hh_bminus[algo][samp]->GetBinContent(x+1);}   else { data_hh_bminus     = 0 ;  }
        if( h_hh_trigger[algo][samp]  != 0 ){ data_hh_trigger  = h_hh_trigger[algo][samp]->GetBinContent(x+1);}  else { data_hh_trigger    = 0 ;  }  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double syst_plus     = fabs( data_hh[x]-data_hh_plus );         double syst_minus    = fabs( data_hh[x]-data_hh_minus );   
        double syst_uncplus  = fabs( data_hh[x]-data_hh_uncplus );      double syst_uncminus = fabs( data_hh[x]-data_hh_uncminus );   
        double syst_jerplus  = fabs( data_hh[x]-data_hh_jerplus );      double syst_jerminus = fabs( data_hh[x]-data_hh_jerminus );
        double syst_bplus    = fabs( data_hh[x]-data_hh_bplus );        double syst_bminus   = fabs( data_hh[x]-data_hh_bminus );      

        double syst_lum      = LUM_ERR*data_hh[x];         
        double syst_eff_lep  = LEP_EFF_ERR*data_hh[x];             
        double syst_ttbar_cs = TTBAR_CS_ERR*data_hh[x];
        
        double syst_tauid(0);  if( x >= TAUSTEP) { syst_tauid = TAU_ID_ERR*data_hh[x];   }  
        double syst_trigger = data_hh_trigger;
        
        if(systset1){ syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_trigger=0; } 
        if(systset2){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_trigger=0;
        } 
        if(systset3){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; 
        } 
     
        double temp =  pow(syst_lum,2) + pow(syst_eff_lep,2) + pow(syst_ttbar_cs,2) + pow(syst_tauid,2) + pow(syst_trigger,2);

        data_syst_plus.push_back(  sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2)  + temp) );   
        data_syst_minus.push_back( sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp) ); 
  
        //cout<<endl<<" HH algo : "<<algo<<" data syst plus "<<(data_syst_plus[x])<<" data syst minus "<<(data_syst_minus[x])<<" temp "<<temp<<endl;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      }

      hh[key] = data_hh; hh_err[key] = data_hh_err;  hh_syst_minus[key] = data_syst_minus; hh_syst_plus[key] = data_syst_plus; 
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // hw samples ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(uint samp=0; samp < h_hw[algo].size(); samp++){
      pair< TString, int> key(algo,samp);      
      vector<double> data_hw, data_hw_err, data_syst_plus, data_syst_minus, data_hw_blpus; 
      for(int x=0; x<12;x++){
        // get data and stats errors ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double data_hw_plus, data_hw_minus; double data_hw_uncplus, data_hw_uncminus, data_hw_jerplus, data_hw_jerminus, data_hw_bplus, data_hw_bminus, data_hw_trigger; 
     
        if( h_hw[algo][samp]          !=0 ){ data_hw.push_back(       h_hw[algo][samp]->GetBinContent(x+1) );  } else { data_hw.push_back( 0 );     }
        if( h_hw[algo][samp]          !=0 ){ data_hw_err.push_back(   h_hw[algo][samp]->GetBinError(x+1)   );  } else { data_hw_err.push_back( 0 ); }
        if( h_hw_plus[algo][samp]     !=0 ){ data_hw_plus     = h_hw_plus[algo][samp]->GetBinContent(x+1);     } else { data_hw_plus       = 0 ;    }
        if( h_hw_minus[algo][samp]    !=0 ){ data_hw_minus    = h_hw_minus[algo][samp]->GetBinContent(x+1);    } else { data_hw_minus      = 0 ;    }  
        if( h_hw_uncplus[algo][samp]  !=0 ){ data_hw_uncplus  = h_hw_uncplus[algo][samp]->GetBinContent(x+1);  } else { data_hw_uncplus    = 0 ;    }
        if( h_hw_uncminus[algo][samp] !=0 ){ data_hw_uncminus = h_hw_uncminus[algo][samp]->GetBinContent(x+1); } else { data_hw_uncminus   = 0 ;    }
        if( h_hw_jerplus[algo][samp]  !=0 ){ data_hw_jerplus  = h_hw_jerplus[algo][samp]->GetBinContent(x+1);  } else { data_hw_jerplus    = 0 ;    }
        if( h_hw_jerminus[algo][samp] !=0 ){ data_hw_jerminus = h_hw_jerminus[algo][samp]->GetBinContent(x+1); } else { data_hw_jerminus   = 0 ;    }  
        if( h_hw_bplus[algo][samp]    !=0 ){ data_hw_bplus    = h_hw_bplus[algo][samp]->GetBinContent(x+1);    } else { data_hw_bplus      = 0 ;    }
        if( h_hw_bminus[algo][samp]   !=0 ){ data_hw_bminus   = h_hw_bminus[algo][samp]->GetBinContent(x+1);   } else { data_hw_bminus     = 0 ;    }  
        if( h_hw_trigger[algo][samp]  != 0){ data_hw_trigger  = h_hw_trigger[algo][samp]->GetBinContent(x+1);  } else { data_hw_trigger    = 0 ;    }  
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double syst_plus     = fabs( data_hw[x]-data_hw_plus );      double syst_minus    = fabs( data_hw[x]-data_hw_minus );   
        double syst_uncplus  = fabs( data_hw[x]-data_hw_uncplus );   double syst_uncminus = fabs( data_hw[x]-data_hw_uncminus );   
        double syst_jerplus  = fabs( data_hw[x]-data_hw_jerplus );   double syst_jerminus = fabs( data_hw[x]-data_hw_jerminus );   
        double syst_bplus    = fabs( data_hw[x]-data_hw_bplus );     double syst_bminus   = fabs( data_hw[x]-data_hw_bminus );   
        double syst_trigger  = data_hw_trigger;   
   
       //if(samp==0 ){
       //  cout<<endl<<" x "<<x<<" minus    = "<<data_hw_minus<<" plus = "<<data_hw_plus<<" uncplus  = "<<data_hw_uncplus<<" unc minus = "<<data_hw_uncminus
       //                      <<" jer plus = "<<data_hw_jerplus<<" jer minus ="<<data_hw_jerminus<<endl;
       //}
  

        double syst_lum      = LUM_ERR*data_hw[x];         
        double syst_eff_lep  = LEP_EFF_ERR*data_hw[x];             
        double syst_ttbar_cs = TTBAR_CS_ERR*data_hw[x];
        
        double syst_tauid(0);  if( x >= TAUSTEP   ){ syst_tauid = TAU_ID_ERR*data_hw[x];   }  
        
        if(systset1){ syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0;  syst_trigger=0;} 
        if(systset2){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_trigger=0;
        } 
        if(systset3){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; 
        }       
   
        double temp = pow(syst_lum,2) + pow(syst_eff_lep,2) + pow(syst_ttbar_cs,2) + pow(syst_tauid,2);

        data_syst_plus.push_back(  sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2) + pow(syst_bplus,2)   + temp));  
        data_syst_minus.push_back( sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp)); 
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

      }

      hw[key] = data_hw; hw_err[key] = data_hw_err;  hw_syst_plus[key]  = data_syst_plus;  hw_syst_minus[key] = data_syst_minus;  
 
 
    }

    // tbh samples //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(uint samp=0; samp < h_tbh[algo].size(); samp++){
      pair< TString, int> key(algo,samp);      
      vector<double> data_tbh, data_tbh_err, data_syst_plus, data_syst_minus; 
      for(int x=0; x<12;x++){
        // get data and stats errors //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double data_tbh_plus, data_tbh_minus;
        double data_tbh_uncplus, data_tbh_uncminus, data_tbh_jerplus, data_tbh_jerminus, data_tbh_bplus, data_tbh_bminus, data_tbh_trigger;
        if( h_tbh[algo][samp]          != 0 ){ data_tbh.push_back(      h_tbh[algo][samp]->GetBinContent(x+1) );  } else { data_tbh.push_back( 0 );   }
        if( h_tbh[algo][samp]          != 0 ){ data_tbh_err.push_back(  h_tbh[algo][samp]->GetBinError(x+1)   );  } else { data_tbh_err.push_back(0); }
        if( h_tbh_plus[algo][samp]     != 0 ){ data_tbh_plus  =    h_tbh_plus[algo][samp]->GetBinContent(x+1);    } else { data_tbh_plus       = 0 ;  }
        if( h_tbh_minus[algo][samp]    != 0 ){ data_tbh_minus =    h_tbh_minus[algo][samp]->GetBinContent(x+1);   } else { data_tbh_minus      = 0 ;  }  
        if( h_tbh_uncplus[algo][samp]  != 0 ){ data_tbh_uncplus  = h_tbh_uncplus[algo][samp]->GetBinContent(x+1); } else { data_tbh_uncplus    = 0 ;  }
        if( h_tbh_uncminus[algo][samp] != 0 ){ data_tbh_uncminus = h_tbh_uncminus[algo][samp]->GetBinContent(x+1);} else { data_tbh_uncminus   = 0 ;  }  
        if( h_tbh_jerplus[algo][samp]  != 0 ){ data_tbh_jerplus  = h_tbh_jerplus[algo][samp]->GetBinContent(x+1); } else { data_tbh_jerplus    = 0 ;  }
        if( h_tbh_jerminus[algo][samp] != 0 ){ data_tbh_jerminus = h_tbh_jerminus[algo][samp]->GetBinContent(x+1);} else { data_tbh_jerminus   = 0 ;  }  
        if( h_tbh_bplus[algo][samp]    != 0 ){ data_tbh_bplus    = h_tbh_bplus[algo][samp]->GetBinContent(x+1); }   else { data_tbh_bplus      = 0 ;  }
        if( h_tbh_bminus[algo][samp]   != 0 ){ data_tbh_bminus   = h_tbh_bminus[algo][samp]->GetBinContent(x+1);}   else { data_tbh_bminus     = 0 ;  }
        if( h_tbh_trigger[algo][samp]  != 0 ){ data_tbh_trigger  = h_tbh_trigger[algo][samp]->GetBinContent(x+1);}  else { data_tbh_trigger    = 0 ;  }  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double syst_plus     = fabs( data_tbh[x]-data_tbh_plus );         double syst_minus    = fabs( data_tbh[x]-data_tbh_minus );   
        double syst_uncplus  = fabs( data_tbh[x]-data_tbh_uncplus );      double syst_uncminus = fabs( data_tbh[x]-data_tbh_uncminus );   
        double syst_jerplus  = fabs( data_tbh[x]-data_tbh_jerplus );      double syst_jerminus = fabs( data_tbh[x]-data_tbh_jerminus );
        double syst_bplus    = fabs( data_tbh[x]-data_tbh_bplus );        double syst_bminus   = fabs( data_tbh[x]-data_tbh_bminus );      

        double syst_lum      = LUM_ERR*data_tbh[x];         
        double syst_eff_lep  = LEP_EFF_ERR*data_tbh[x];             
        double syst_ttbar_cs = TTBAR_CS_ERR*data_tbh[x];
        
        double syst_tauid(0);  if( x >= TAUSTEP) { syst_tauid = TAU_ID_ERR*data_tbh[x];   }  
        double syst_trigger = data_tbh_trigger;
        
        if(systset1){ syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_trigger=0; } 
        if(systset2){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_trigger=0;
        } 
        if(systset3){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; 
        } 
     
        double temp =  pow(syst_lum,2) + pow(syst_eff_lep,2) + pow(syst_ttbar_cs,2) + pow(syst_tauid,2) + pow(syst_trigger,2);

        data_syst_plus.push_back(  sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2)  + temp) );   
        data_syst_minus.push_back( sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp) ); 
  
        //cout<<endl<<" TBH algo : "<<algo<<" data syst plus "<<(data_syst_plus[x])<<" data syst minus "<<(data_syst_minus[x])<<" temp "<<temp<<endl;

	// WORKING HERE
	
	if(x==11){ // OS

	  cout << "=========================== DEBUG: samp: " << samp << ", BRHTAUNU[samp]: " << brHtaunu_[samp] << endl;
	  vector<double> temp_tbh_datacards;
	  temp_tbh_datacards.push_back( data_tbh[x] );
	  temp_tbh_datacards.push_back( data_tbh_err[x]  );
	  temp_tbh_datacards.push_back( syst_plus     );
	  temp_tbh_datacards.push_back( syst_minus    );
	  temp_tbh_datacards.push_back( syst_uncplus  );
	  temp_tbh_datacards.push_back( syst_uncminus );
	  temp_tbh_datacards.push_back( syst_jerplus  );
	  temp_tbh_datacards.push_back( syst_jerminus );
	  temp_tbh_datacards.push_back( syst_bplus    );
	  temp_tbh_datacards.push_back( syst_bminus   );
	  temp_tbh_datacards.push_back( syst_lum      );
	  temp_tbh_datacards.push_back( syst_eff_lep  );
	  //        temp_tbh_datacards.push_back( syst_ttbar_cs );
	  temp_tbh_datacards.push_back( syst_tauid    );
	  temp_tbh_datacards.push_back( syst_trigger  );

	  
	  tbh_datacards.push_back(temp_tbh_datacards);
	}
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      }
      

      tbh[key] = data_tbh; tbh_err[key] = data_tbh_err;  tbh_syst_minus[key] = data_syst_minus; tbh_syst_plus[key] = data_syst_plus; 
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // htb samples //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(uint samp=0; samp < h_htb[algo].size(); samp++){
      pair< TString, int> key(algo,samp);      
      vector<double> data_htb, data_htb_err, data_syst_plus, data_syst_minus; 
      for(int x=0; x<12;x++){
        // get data and stats errors //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double data_htb_plus, data_htb_minus;
        double data_htb_uncplus, data_htb_uncminus, data_htb_jerplus, data_htb_jerminus, data_htb_bplus, data_htb_bminus, data_htb_trigger;
        if( h_htb[algo][samp]          != 0 ){ data_htb.push_back(      h_htb[algo][samp]->GetBinContent(x+1) );  } else { data_htb.push_back( 0 );   }
        if( h_htb[algo][samp]          != 0 ){ data_htb_err.push_back(  h_htb[algo][samp]->GetBinError(x+1)   );  } else { data_htb_err.push_back(0); }
        if( h_htb_plus[algo][samp]     != 0 ){ data_htb_plus  =    h_htb_plus[algo][samp]->GetBinContent(x+1);    } else { data_htb_plus       = 0 ;  }
        if( h_htb_minus[algo][samp]    != 0 ){ data_htb_minus =    h_htb_minus[algo][samp]->GetBinContent(x+1);   } else { data_htb_minus      = 0 ;  }  
        if( h_htb_uncplus[algo][samp]  != 0 ){ data_htb_uncplus  = h_htb_uncplus[algo][samp]->GetBinContent(x+1); } else { data_htb_uncplus    = 0 ;  }
        if( h_htb_uncminus[algo][samp] != 0 ){ data_htb_uncminus = h_htb_uncminus[algo][samp]->GetBinContent(x+1);} else { data_htb_uncminus   = 0 ;  }  
        if( h_htb_jerplus[algo][samp]  != 0 ){ data_htb_jerplus  = h_htb_jerplus[algo][samp]->GetBinContent(x+1); } else { data_htb_jerplus    = 0 ;  }
        if( h_htb_jerminus[algo][samp] != 0 ){ data_htb_jerminus = h_htb_jerminus[algo][samp]->GetBinContent(x+1);} else { data_htb_jerminus   = 0 ;  }  
        if( h_htb_bplus[algo][samp]    != 0 ){ data_htb_bplus    = h_htb_bplus[algo][samp]->GetBinContent(x+1); }   else { data_htb_bplus      = 0 ;  }
        if( h_htb_bminus[algo][samp]   != 0 ){ data_htb_bminus   = h_htb_bminus[algo][samp]->GetBinContent(x+1);}   else { data_htb_bminus     = 0 ;  }
        if( h_htb_trigger[algo][samp]  != 0 ){ data_htb_trigger  = h_htb_trigger[algo][samp]->GetBinContent(x+1);}  else { data_htb_trigger    = 0 ;  }  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double syst_plus     = fabs( data_htb[x]-data_htb_plus );         double syst_minus    = fabs( data_htb[x]-data_htb_minus );   
        double syst_uncplus  = fabs( data_htb[x]-data_htb_uncplus );      double syst_uncminus = fabs( data_htb[x]-data_htb_uncminus );   
        double syst_jerplus  = fabs( data_htb[x]-data_htb_jerplus );      double syst_jerminus = fabs( data_htb[x]-data_htb_jerminus );
        double syst_bplus    = fabs( data_htb[x]-data_htb_bplus );        double syst_bminus   = fabs( data_htb[x]-data_htb_bminus );      

        double syst_lum      = LUM_ERR*data_htb[x];         
        double syst_eff_lep  = LEP_EFF_ERR*data_htb[x];             
        double syst_ttbar_cs = TTBAR_CS_ERR*data_htb[x];
        
        double syst_tauid(0);  if( x >= TAUSTEP) { syst_tauid = TAU_ID_ERR*data_htb[x];   }  
        double syst_trigger = data_htb_trigger;
        
        if(systset1){ syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_trigger=0; } 
        if(systset2){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_trigger=0;
        } 
        if(systset3){ 
          syst_lum=0; syst_tauid=0; syst_ttbar_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; 
        } 
     
        double temp =  pow(syst_lum,2) + pow(syst_eff_lep,2) + pow(syst_ttbar_cs,2) + pow(syst_tauid,2) + pow(syst_trigger,2);

        data_syst_plus.push_back(  sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2)  + temp) );   
        data_syst_minus.push_back( sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp) ); 
  
        //cout<<endl<<" HTB algo : "<<algo<<" data syst plus "<<(data_syst_plus[x])<<" data syst minus "<<(data_syst_minus[x])<<" temp "<<temp<<endl;

	// WORKING HERE
	
	if(x==11){ // OS

	  cout << "=========================== DEBUG: samp: " << samp << ", BRHTAUNU[samp]: " << brHtaunu_[samp] << endl;
	  vector<double> temp_htb_datacards;
	  temp_htb_datacards.push_back( data_htb[x] );
	  temp_htb_datacards.push_back( data_htb_err[x]  );
	  temp_htb_datacards.push_back( syst_plus     );
	  temp_htb_datacards.push_back( syst_minus    );
	  temp_htb_datacards.push_back( syst_uncplus  );
	  temp_htb_datacards.push_back( syst_uncminus );
	  temp_htb_datacards.push_back( syst_jerplus  );
	  temp_htb_datacards.push_back( syst_jerminus );
	  temp_htb_datacards.push_back( syst_bplus    );
	  temp_htb_datacards.push_back( syst_bminus   );
	  temp_htb_datacards.push_back( syst_lum      );
	  temp_htb_datacards.push_back( syst_eff_lep  );
	  //        temp_htb_datacards.push_back( syst_ttbar_cs );
	  temp_htb_datacards.push_back( syst_tauid    );
	  temp_htb_datacards.push_back( syst_trigger  );

	  
	  htb_datacards.push_back(temp_htb_datacards);
	}
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      }
      

      htb[key] = data_htb; htb_err[key] = data_htb_err;  htb_syst_minus[key] = data_syst_minus; htb_syst_plus[key] = data_syst_plus; 
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << "Processed hw, hh, tbh, htb histos, if any" << endl;



    // fill data points////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int x=0;x<12;x++){
    
      // fill data /////////////////////////////////////////////////////////////////////////////////
      pair<TString,int> dk(algo,x); if( hData[algo] ) dataTotal[dk]= hData[algo]->GetBinContent(x+1);
      if(x==11) // OS
	data_datacards.push_back( dataTotal[dk] );
      //////////////////////////////////////////////////////////////////////////////////////////////
      
      // fill values with histogram data /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int y=0;y<13;y++){  

        //map< pair< TString, pair<int,int>, double> 
        pair<int,int> keytemp(x,y);
        pair< TString, pair<int,int> > k(algo,keytemp);

        if(hMC[algo])          { dMC[k]          = (hMC[algo]->GetBinContent(x+1,y+1)); dMC_err[k] = (hMC[algo]->GetBinError(x+1,y+1));              } 
        if(hMC_plus[algo])     { dMC_plus[k]     = (hMC_plus[algo]->GetBinContent(x+1,y+1));                                                         } 
        if(hMC_minus[algo])    { dMC_minus[k]    = (hMC_minus[algo]->GetBinContent(x+1,y+1));                                                        }
        if(hMC_uncplus[algo])  { dMC_uncplus[k]  = (hMC_uncplus[algo]->GetBinContent(x+1,y+1));                                                      } 
        if(hMC_uncminus[algo]) { dMC_uncminus[k] = (hMC_uncminus[algo]->GetBinContent(x+1,y+1));                                                     }
        if(hMC_jerplus[algo])  { dMC_jerplus[k]  = (hMC_jerplus[algo]->GetBinContent(x+1,y+1));                                                      }   
        if(hMC_jerminus[algo]) { dMC_jerminus[k] = (hMC_jerminus[algo]->GetBinContent(x+1,y+1));                                                     }
        if(hMC_bplus[algo])    { dMC_bplus[k]    = (hMC_bplus[algo]->GetBinContent(x+1,y+1));                                                        } 
        if(hMC_bminus[algo])   { dMC_bminus[k]   = (hMC_bminus[algo]->GetBinContent(x+1,y+1));                                                       } 
        if(hMC_trigger[algo])  { dMC_trigger[k]  = (hMC_trigger[algo]->GetBinContent(x+1,y+1));                                                      } 

        if( XSECMEASUREMENT && y==LEPLEP_TDCH ){   
          // Note : we re-scale ttbar lepton+lepton ( we get it from MC)
          double scalefactor( XSEC_EXP_/234.); //note it should be 234
          dMC[k]          *= scalefactor;
          dMC_plus[k]     *= scalefactor; 
          dMC_minus[k]    *= scalefactor;
          dMC_uncplus[k]  *= scalefactor;
          dMC_uncminus[k] *= scalefactor;
          dMC_jerplus[k]  *= scalefactor;
          dMC_jerminus[k] *= scalefactor;
          dMC_bplus[k]    *= scalefactor;
          dMC_trigger[k]  *= scalefactor;
        }

        //compute systematics //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double syst_plus  = fabs(dMC[k]-dMC_plus[k]);    double syst_uncplus  = fabs(dMC[k]-dMC_uncplus[k]);    double syst_jerplus  = fabs(dMC[k]-dMC_jerplus[k]); 
        double syst_minus = fabs(dMC[k]-dMC_minus[k]);   double syst_uncminus = fabs(dMC[k]-dMC_uncminus[k]);   double syst_jerminus = fabs(dMC[k]-dMC_jerminus[k]);

        double syst_bplus    = fabs(dMC[k]-dMC_bplus[k]);  double syst_bminus = fabs(dMC[k]-dMC_bminus[k]); 
        double syst_lum      = LUM_ERR*dMC[k];
        double syst_eff_lep  = LEP_EFF_ERR*dMC[k];
        double syst_trigger  = dMC_trigger[k];
        double syst_tauid(0);
        double syst_cs(0), sample_cs(0);
  
        
        if(      y==QCD_TDCH )                                                                                                               { sample_cs = QCD_CS_ERR ;      }
        else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
        else if( y == SINGLETOP_TDCH )                                                                                                       { sample_cs = SINGLETOP_CS_ERR; }              
        else if( y == WJETS_TDCH     )                                                                                                       { sample_cs = WJETS_CS_ERR;     }              
        else                                                                                                                                 { sample_cs = OTHER_CS_ERR;     }
        

        if(XSECMEASUREMENT && y==LEPLEP_TDCH ){ sample_cs = XSEC_EXP_ERR_; }

        syst_cs = sample_cs * dMC[k];
 
        double syst_theory(0);
        if( XSECMEASUREMENT && (y==ETAU_TDCH || y==MUTAU_TDCH) ){ syst_theory=SIGNAL_EFF_ERR*dMC[k];}


        if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=LEPLEP_TDCH   ){ syst_tauid = TAU_ID_ERR*dMC[k];     }  
        else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==LEPLEP_TDCH )){ syst_tauid = TAU_MISSID_ERR*dMC[k]; }

	if(x==11){ //OS
	  vector<double> temp_sm_datacards;

	  temp_sm_datacards.push_back(dMC[k] );
	  temp_sm_datacards.push_back(dMC_err[k] );	  
	  temp_sm_datacards.push_back( syst_plus     );
	  temp_sm_datacards.push_back( syst_minus    );
	  temp_sm_datacards.push_back( syst_uncplus  );
	  temp_sm_datacards.push_back( syst_uncminus );
	  temp_sm_datacards.push_back( syst_jerplus  );
	  temp_sm_datacards.push_back( syst_jerminus );
	  temp_sm_datacards.push_back( syst_bplus    );
	  temp_sm_datacards.push_back( syst_bminus   );
	  temp_sm_datacards.push_back( syst_lum      );
	  temp_sm_datacards.push_back( syst_eff_lep  );
	  // syst_theory? syst_cs? sample_cs? // FIXME: fixme
	  temp_sm_datacards.push_back( syst_tauid    );
	  temp_sm_datacards.push_back( syst_trigger  );
	  sm_datacards.push_back(temp_sm_datacards);
	  cout<<endl<< "KNOW BETTER WHAT IS WHAT " << sm_datacards.size() << " " << temp_sm_datacards[0] << endl; 
	}
	
        if(systset1){ syst_lum=0; syst_tauid=0; syst_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_theory=0; syst_trigger=0; } 
        if(systset2){ syst_lum=0; syst_tauid=0; syst_cs=0; syst_eff_lep=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_theory=0; syst_trigger=0;} 
        if(systset3){ syst_lum=0; syst_tauid=0; syst_cs=0; syst_eff_lep=0; syst_bplus=0; syst_bminus=0; syst_plus=0; syst_uncplus=0; syst_jerplus=0; syst_minus=0; syst_uncminus=0; syst_jerminus=0; syst_theory=0; } 


        double temp = pow(syst_lum,2) + pow(syst_tauid,2) + pow(syst_cs,2) + pow(syst_eff_lep,2) + pow(syst_theory,2);

        dMC_syst_plus[k] = sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2) + temp ); 
        dMC_syst_minus[k]= sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + temp ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

      }

    }
  }

  // write document...
  FILE * f = fopen(fName.Data(),"w");

  if(incDocument) prepareDocument(f);     

  //note size of table is restriced to 2 or 3 elements
  //if(allSysErrors)fprintf(f,"\\begin{footnotesize} \n");
  if(incDocument)  fprintf(f,"\\begin{center} \n"); 
  fprintf(f,"\\begin{tabular}{ c c } \n");     fprintf(f,"\\multicolumn{2}{c}{ } \\\\ \n"); 
  fprintf(f,"\\hline \n");
  fprintf(f,"\\multicolumn{1}{ c }{ Source } & \\multicolumn{1}{ c }{ $N_{events}$ $\\pm$ stat $\\pm$ syst } \\\\ \n"); 
 

  TString  T_HW_mc[] = {
    "HW, $M_{H^{+}}=80GeV/c^{2}$ &",  "HW, $M_{H^{+}}=100GeV/c^{2}$ &", "HW, $M_{H^{+}}=120GeV/c^{2}$ &",  "HW, $M_{H^{+}}=140GeV/c^{2}$ &",
    "HW, $M_{H^{+}}=150GeV/c^{2}$ &", "HW, $M_{H^{+}}=155GeV/c^{2}$ &", "HW, $M_{H^{+}}=160GeV/c^{2}$ &"
  };

  TString  T_HH_mc[] = {
    "HH, $M_{H^{+}}=80GeV/c^{2}$ &",  "HH, $M_{H^{+}}=100GeV/c^{2}$ &", "HH, $M_{H^{+}}=120GeV/c^{2}$ &", "HH, $M_{H^{+}}=140GeV/c^{2}$ &",
    "HH, $M_{H^{+}}=150GeV/c^{2}$ &", "HH, $M_{H^{+}}=155GeV/c^{2}$ &", "HH, $M_{H^{+}}=160GeV/c^{2}$ &"
  }; 

  TString  T_TBH_mc[] = {
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=180GeV/c^{2}$ &",
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=190GeV/c^{2}$ &",
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=200GeV/c^{2}$ &",
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=220GeV/c^{2}$ &",
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=250GeV/c^{2}$ &",
    "H^{+}#rightarrow#tau_{#mathrm{h}}#nu, $M_{H^{+}}=300GeV/c^{2}$ &"
  };

  TString  T_HTB_mc[] = {
    "H^{+}#rightarrow tb, $M_{H^{+}}=180GeV/c^{2}$ &",
    "H^{+}#rightarrow tb, $M_{H^{+}}=190GeV/c^{2}$ &",
    "H^{+}#rightarrow tb, $M_{H^{+}}=200GeV/c^{2}$ &",
    "H^{+}#rightarrow tb, $M_{H^{+}}=220GeV/c^{2}$ &",
    "H^{+}#rightarrow tb, $M_{H^{+}}=250GeV/c^{2}$ &",
    "H^{+}#rightarrow tb, $M_{H^{+}}=300GeV/c^{2}$ &"
  };

  TString    other[] = {
    "$t\\bar{t}$ tau-dilepton &","$\\tau$ fakes &", "$t\\bar{t}\\to\\ell\\ell$ &", "$Z/\\gamma\\rightarrow ee,\\mu\\mu$ &", "$Z/\\gamma\\rightarrow \\tau\\tau$ &",
    "single top &","WW,WZ,ZZ &"
  };

  fprintf(f,"\\hline \n");
  fprintf(f,"\\hline \n");

  TString algo0(keys[0]); 


  // HH contribution
  for(uint l=0; l<hhFiles.size(); l++){
    if( !showAllHiggs && l!= 2 ){ continue; }  // only show 120 GeV/c2
    TString line(""); line.Append(T_HH_mc[l]);
    
    if(detailed){ TString d(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); line.Append(d); }
    else{         TString d(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); line.Append(d); }
    
    int massPoint(l);
    pair<TString, int> k0(algo0,massPoint); 
    const char * data = line.Data();
    double systerr0(hh_syst_plus[k0][STEP]); if( hh_syst_minus[k0][STEP] > hh_syst_plus[k0][STEP] ) systerr0 = hh_syst_minus[k0][STEP];
    fprintf(f, data, hh[k0][STEP], hh_err[k0][STEP], systerr0 );
    
  }

  // HW contribution
  for(uint l=0; l<hhFiles.size(); l++){
    if( !showAllHiggs && l!= 2  ){ continue;}
     
    TString line(""); line.Append(T_HW_mc[l]); 

    if(detailed){ TString d(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); line.Append(d);}
    else        { TString d(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); line.Append(d);}                                    
    
    int massPoint(l);
    pair<TString, int> k0(algo0,massPoint); 
    const char * data = line.Data();
    double systerr0(hw_syst_plus[k0][STEP]); if( hw_syst_minus[k0][STEP] > hw_syst_plus[k0][STEP] ) systerr0 = hw_syst_minus[k0][STEP];
   
    fprintf(f, data, hw[k0][STEP], hw_err[k0][STEP], systerr0 );
   
    
  }


  // TBH contribution
  for(uint l=0; l<tbhFiles.size(); l++){
    if( !showAllHiggs && l!= 2  ){ continue;}
    
    TString line(""); line.Append(T_TBH_mc[l]); 
    
    if(detailed){ TString d(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); line.Append(d);}
    else        { TString d(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); line.Append(d);}                                    
    
    int massPoint(l);
    pair<TString, int> k0(algo0,massPoint); 
    const char * data = line.Data();
    double systerr0(tbh_syst_plus[k0][STEP]); if( tbh_syst_minus[k0][STEP] > tbh_syst_plus[k0][STEP] ) systerr0 = tbh_syst_minus[k0][STEP];
   
    fprintf(f, data, brHtaunu_[l]*tbh[k0][STEP], brHtaunu_[l]*tbh_err[k0][STEP], brHtaunu_[l]*systerr0 );
   
    
  }

  // HTB contribution
  for(uint l=0; l<htbFiles.size(); l++){
    if( !showAllHiggs && l!= 2  ){ continue;}
    
    TString line(""); line.Append(T_HTB_mc[l]); 
    
    if(detailed){ TString d(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); line.Append(d);}
    else        { TString d(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); line.Append(d);}                                    
    
    int massPoint(l);
    pair<TString, int> k0(algo0,massPoint); 
    const char * data = line.Data();
    double systerr0(tbh_syst_plus[k0][STEP]); if( htb_syst_minus[k0][STEP] > htb_syst_plus[k0][STEP] ) systerr0 = htb_syst_minus[k0][STEP];
   
    fprintf(f, data, brHtb_[l]*htb[k0][STEP], brHtb_[l]*htb_err[k0][STEP], brHtb_[l]*systerr0 );
   
    
  }


  //tt tau-dilepton
  {    
    fprintf(f,"\\hline \n");        
    pair<int,int> ktemp0_e(STEP,0);  pair<TString ,pair<int,int> > e0(algo0,ktemp0_e); pair<int,int> ktemp0_mu(STEP,1); pair<TString ,pair<int,int> > mu0(algo0,ktemp0_mu);
    


    // first algorithm /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double signal_0= (dMC[e0]+dMC[mu0]); 
    double signal_plus_0=(dMC_plus[e0]+dMC_plus[mu0]);             double signal_minus_0= (dMC_minus[e0]+dMC_minus[mu0]);
    double signal_uncplus_0(dMC_uncplus[e0]+dMC_uncplus[mu0]);     double signal_uncminus_0(dMC_uncminus[e0]+dMC_uncminus[mu0]);
    double signal_jerplus_0(dMC_jerplus[e0]+dMC_jerplus[mu0]);     double signal_jerminus_0(dMC_jerminus[e0]+dMC_jerminus[mu0]);
    double signal_bplus_0(dMC_bplus[e0]+dMC_bplus[mu0]);           double signal_bminus_0(dMC_bminus[e0]+dMC_bminus[mu0]);
    double signal_trigger = sqrt( pow(dMC_trigger[e0],2)+pow(dMC_trigger[mu0],2) );
 
    double err_0          = sqrt( dMC_err[e0]*dMC_err[e0] + dMC_err[mu0]*dMC_err[mu0] );      
    double syst_jes_0     = fabs ( signal_0 - signal_plus_0);    if( fabs( signal_0 - signal_minus_0    ) > syst_jes_0 )  syst_jes_0  = fabs( signal_0 - signal_minus_0 );
    double syst_unc_0     = fabs ( signal_0 - signal_uncplus_0); if( fabs( signal_0 - signal_uncminus_0 ) > syst_unc_0 )  syst_unc_0  = fabs( signal_0 - signal_uncminus_0 );
    double syst_jer_0     = fabs ( signal_0 - signal_jerplus_0); if( fabs( signal_0 - signal_jerminus_0 ) > syst_jer_0 )  syst_jer_0  = fabs( signal_0 - signal_jerminus_0 );
    double syst_btag_0    = fabs ( signal_0 - signal_bplus_0);   if( fabs( signal_0 - signal_bminus_0   ) > syst_btag_0 ) syst_btag_0 = fabs( signal_0 - signal_bminus_0 );

    double syst_tauid_0   = TAU_ID_ERR*signal_0; 
    double syst_ttbar_0   = TTBAR_CS_ERR*signal_0;
    double syst_lum_0     = LUM_ERR*signal_0;
    double syst_eff_lep_0 = LEP_EFF_ERR*signal_0;

    double syst_theory_0  = SIGNAL_EFF_ERR*signal_0;

    if(systset1){ syst_lum_0=0; syst_tauid_0=0; syst_ttbar_0=0; syst_eff_lep_0=0; syst_btag_0=0;  syst_theory_0=0; signal_trigger=0;} 
    if(systset2){
      syst_lum_0=0; syst_tauid_0=0; syst_ttbar_0=0; syst_eff_lep_0=0; 
      syst_jes_0=0; syst_unc_0=0; syst_jer_0=0;
      syst_theory_0=0; signal_trigger=0;
    }
    if(systset3){
      syst_lum_0=0; syst_tauid_0=0; syst_ttbar_0=0; syst_eff_lep_0=0; syst_btag_0=0;
      syst_jes_0=0; syst_unc_0=0; syst_jer_0=0;
      syst_theory_0=0; 
    }


    double temp_0 = pow(syst_tauid_0,2) + pow(syst_ttbar_0,2) + pow(syst_lum_0,2) + pow(syst_eff_lep_0,2)+pow(syst_theory_0,2);

    double syst_0  = sqrt( pow(syst_jes_0,2)+ pow(syst_unc_0,2) + pow(syst_jer_0,2) + pow(syst_btag_0,2)+ temp_0 );
    
    //FORCED :: Need to recheck
    //syst_0=11.9;
    //////////////////////////

    
    cout<<endl<<"\n systematics on tau-dilepton : "
              <<"\n jes        : "<<syst_jes_0
              <<"\n unc        : "<<syst_unc_0
              <<"\n btag       : "<<syst_btag_0
              <<"\n tau id     : "<<syst_tauid_0
              <<"\n ttbar xsec : "<<syst_ttbar_0
              <<"\n lum        : "<<syst_lum_0
              <<"\n lepton eff : "<<syst_eff_lep_0
              <<"\n theory unc : "<<syst_theory_0<<endl;



    tau_dilepton = signal_0; tau_dilepton_stat = err_0;  tau_dilepton_syst = syst_0; // acceptance computation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    TString line(""); line.Append(other[0]); TString d;


    if(detailed) d = TString(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); 
    else         d = TString(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); 
    line.Append(d); 

    

    const char * datatemp = line.Data(); fprintf(f, datatemp, signal_0, err_0, syst_0);  
    
    
  }

  
  cout<<"\n tau dilepton -> "<<tau_dilepton<<" stat +/- "<<tau_dilepton_stat<<" syst +/- "<<tau_dilepton_syst<<endl;



  // Now we can compute the acceptance ///////////////////////////////////////////////////////////////

  // acceptance        = (tau_dilepton/scale)/ttbar_init;
  // acceptance_stat   = getErrorFraction((tau_dilepton/scale),ttbar_init);

  acceptance        = tau_dilepton/(ttbar_init*scale);
  acceptance_stat   = getErrorFraction( tau_dilepton ,ttbar_init*scale);

  double true_error = getErrorFraction( tau_dilepton, ttbar_init*scale, tau_dilepton_stat , sqrt(ttbar_init)*scale );

  acceptance_syst = tau_dilepton_syst/(ttbar_init*scale);


//  cout<<endl<<" Acceptance   : "<<acceptance<<" +- "<<acceptance_stat<<" +- "<<acceptance_syst<<endl;
  cout<<endl<<" Acceptance   : "<<acceptance<<" +- "<<true_error<<" +- "<<acceptance_syst<<endl;
  cout<<endl<<" tau_dilepton : "<<tau_dilepton<<endl;

  cout<<endl<<" Numerator for acceptance : "<<(tau_dilepton)<<endl;
  cout<<endl<<" Initial is               : "<<ttbar_init<<endl<<endl;
  cout<<endl<<" scale                    : "<<scale<<endl<<endl;

  // new
  acceptance_stat = true_error;

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  // tau fakes
  { 

    TString line(""); line.Append(other[1]); TString d("");

    if(detailed)             d = TString(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); 
    else                     d = TString(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); 
    line.Append(d); 
    const char * data = line.Data(); fprintf(f, data, tau_fakes, sqrt(tau_fakes_stat2), sqrt(tau_fakes_syst2) );        
  }



  //LEPLEP_TDCH ;ZJETS_TDCH; ZTAUTAU_TDCH; SINGLETOP_TDCH; EWKDI_TDCH
  { 
    vector<int> codes; codes.push_back(LEPLEP_TDCH); codes.push_back(ZJETS_TDCH); codes.push_back(ZTAUTAU_TDCH); codes.push_back(SINGLETOP_TDCH); codes.push_back(EWKDI_TDCH); 
    for(uint k=0;k<codes.size();k++){
      TString line(""); line.Append(other[2+k]); 
      TString d("");

      if(detailed) d = TString(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); 
      else d = TString(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); 
      
      line.Append(d);

      pair<int,int> ktemp( STEP, codes[k] ); 
      pair< TString , pair<int,int> > entry_0( algo0, ktemp ); 
             
      double syst_0( dMC_syst_plus[entry_0] ); if( dMC_syst_minus[entry_0] > syst_0 ){ syst_0 = dMC_syst_minus[entry_0]; }
      
      
      other_bkg     += dMC[entry_0]; other_bkg_stat2 += pow(dMC_err[entry_0],2); other_bkg_syst2 += pow(syst_0,2);
      const char * datatemp2 = line.Data(); 


      //if(codes[k]==LEPLEP_TDCH    ){ syst_0 = 1.8; }
      //if(codes[k]==ZTAUTAU_TDCH   ){ syst_0 = 2.0; }
      //if(codes[k]==ZJETS_TDCH     ){ syst_0 = 1.0; }
      //if(codes[k]==LEPLEP_TDCH    ){ syst_0 = 1.8; }
      //if(codes[k]==SINGLETOP_TDCH ){ syst_0 = 0.8; }
      //if(codes[k]==EWKDI_TDCH     ){ syst_0 = 0.2; }

      cout<<endl<< "KNOW WHAT IS WHAT " << k << "     " << line << endl;
      fprintf(f, datatemp2, dMC[entry_0], dMC_err[entry_0], syst_0);   
      
    }
  }


  // total expected
  { 
    fprintf(f,"\\hline \n");
    double total(0);
    pair<int,int> ktemp0_e(STEP,0);  pair<TString ,pair<int,int> > e0(algo0,ktemp0_e); pair<int,int> ktemp0_mu(STEP,1); pair<TString ,pair<int,int> > mu0(algo0,ktemp0_mu);
   
    // first algorithm //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    double total_0(0);           
    double total_err_0(0);
    double total_plus_0(0);
    double total_minus_0(0);
    double total_uncplus_0(0);
    double total_uncminus_0(0);
    double total_jerplus_0(0);
    double total_jerminus_0(0);
    double total_bplus_0(0);
    double total_bminus_0(0);
    double total_trigger_0(0);


    double total_ttbaronly_0(0);
    double total_wjetsonly_0(0);
    double total_singletoponly_0(0);
    double total_zjetsonly_0(0);
    double total_dibosononly_0(0);
    double total_qcdonly_0(0);
    double total_signal_0(0);
    

    //double total_signalEfftheory_0(0);

     
    for(int y=0;y<13;y++){

      if(y == WJETS_TDCH || y== QCD_TDCH || y ==  LEPJETS_TDCH ) { continue;} 

      pair<int,int> k(STEP,y); pair<TString ,pair<int,int> > r_0(algo0,k); 
      total_0          += dMC[r_0];                    
      total_err_0      += dMC_err[r_0]*dMC_err[r_0];   
      total_plus_0     += dMC_plus[r_0];       total_minus_0    += dMC_minus[r_0];             
      total_uncplus_0  += dMC_uncplus[r_0];    total_uncminus_0 += dMC_uncminus[r_0];           
      total_jerplus_0  += dMC_jerplus[r_0];    total_jerminus_0 += dMC_jerminus[r_0];   
      total_bplus_0    += dMC_bplus[r_0];      total_bminus_0   += dMC_bminus[r_0];   
      total_trigger_0  += dMC_trigger[r_0];

      if( ETAU_TDCH || y==MUTAU_TDCH) total_signal_0 += dMC[r_0];

      if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ 
                                                    total_ttbaronly_0      += dMC[r_0]; 
      }
      else if( y==WJETS_TDCH                      ){ total_wjetsonly_0     += dMC[r_0];}
      else if( y==QCD_TDCH                        ){ total_qcdonly_0       += dMC[r_0];}
      else if( y==ZJETS_TDCH || y == ZTAUTAU_TDCH ){ total_zjetsonly_0     += dMC[r_0];}
      else if( y==SINGLETOP_TDCH                  ){ total_singletoponly_0 += dMC[r_0];}
      else if( y==EWKDI_TDCH                      ){ total_dibosononly_0   += dMC[r_0];}
      
    }
     
    total_err_0   = sqrt(total_err_0); 

    double syst_0       = fabs( total_0-total_plus_0);     if(  fabs( total_0 - total_minus_0    ) > syst_0     ){ syst_0     = fabs( total_0 - total_minus_0 ); }  
    double syst_unc0    = fabs( total_0-total_uncplus_0);  if(  fabs( total_0 - total_uncminus_0 ) > syst_unc0  ){ syst_unc0  = fabs( total_0 - total_uncminus_0 ); }
    double syst_jer0    = fabs( total_0-total_jerplus_0);  if(  fabs( total_0 - total_jerminus_0 ) > syst_jer0  ){ syst_jer0  = fabs( total_0 - total_jerminus_0 ); }
    double syst_btag0   = fabs( total_0-total_bplus_0);    if(  fabs( total_0 - total_bminus_0   ) > syst_btag0 ){ syst_btag0 = fabs( total_0 - total_bminus_0 ); }
    double syst_trigger = total_trigger_0;
   
    double temp_0 = pow(TAU_ID_ERR*total_0,2) + pow(LUM_ERR*total_0,2)+ pow(LEP_EFF_ERR*total_0,2)  +
      pow(TTBAR_CS_ERR*total_ttbaronly_0,2) + pow(WJETS_CS_ERR*total_wjetsonly_0,2)  + pow(SINGLETOP_CS_ERR*total_singletoponly_0,2)+ pow(QCD_CS_ERR*total_qcdonly_0,2)  +
      pow(OTHER_CS_ERR*total_zjetsonly_0,2) + pow(OTHER_CS_ERR*total_dibosononly_0,2)+ pow(SIGNAL_EFF_ERR*total_signal_0,2);
                    
                    
   
    if(systset1){ temp_0=0; syst_btag0 =0; syst_trigger=0;} 
    if(systset2){ temp_0=0; syst_0 = 0; syst_unc0=0; syst_jer0=0; syst_trigger=0;}
    if(systset3){ temp_0=0; syst_btag0 =0; syst_0 = 0; syst_unc0=0; syst_jer0=0;}


    syst_0 = sqrt( pow(syst_0,2)+pow(syst_unc0,2)+pow(syst_jer0,2)+pow(syst_btag0,2) +pow(syst_trigger,2)+ temp_0); 


    //NEED TO RECHECK ============================================>
    // note systematics need to supress w+jets and qcde erro
    //syst_0 = 
    //syst_0 =14.6;
    ///////////////////////////////////////////////////////////////



    TString line(""); line.Append("Total expected & "); TString d("");
     
    if(detailed) d=TString(" %8.2f $\\pm$ %8.2f $\\pm$ %8.2f \\\\ \n"); 
    else d=TString(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n");
    line.Append(d);

    //FORCED
    //total_0=124.3;
        
    //

     const char * datatemp3 = line.Data();  fprintf(f, datatemp3, total_0 + tau_fakes, sqrt( pow(total_err_0,2)+ tau_fakes_stat2 ), sqrt( pow(syst_0,2) + tau_fakes_syst2 ));  
    //const char * datatemp3 = line.Data();  fprintf(f, datatemp3, 124.3, sqrt( pow(total_err_0,2)+ tau_fakes_stat2 ), 14.6);  

  }

  // total data
  { 
    fprintf(f,"\\hline \n");
    pair< TString, int > k_0( algo0, STEP );          
    pair< TString, int > k_0_BTAG( algo0, BTAG1STEP ); 
    double total_0      = dataTotal[k_0];              
    double total_0_BTAG = dataTotal[k_0_BTAG];          

    data_yields=total_0;

   // TString line(""); line.Append("Data (b-tagged) &"); TString d("");
   TString line(""); line.Append("Data &"); TString d("");
  
  
    
    //d=TString(" %8.0f (%8.0f)  \\\\ \n"); line.Append(d);
    //const char * data = line.Data();  fprintf(f, data, total_0, total_0_BTAG);  
    d=TString(" %8.0f  \\\\ \n"); line.Append(d);
    const char * datatemp4 = line.Data();  fprintf(f, datatemp4, total_0 );  
    

  }
  fprintf(f,"\\hline \n"); 


  if(addXsection){

    // acceptance line //////////////////////////////////////////////////////
    TString line(""); line.Append("Acceptance &"); TString d("");
    d=TString(" %.2e $\\pm$ %.1e $\\pm$ %.1e \\\\ \n"); line.Append(d);
    const char * data = line.Data();  fprintf(f, data, acceptance, acceptance_stat, acceptance_syst);  
    fprintf(f,"\\hline \n"); 
    ////////////////////////////////////////////////////////////////////////

    // xsect computation ///////////////////////////////////////////////////////////////////////
        
    
    //    data_yields=176;
    //    data_yields=369;
    //    data_yields=434;
    //    data_yields=625;
    double bkg       = tau_fakes       + other_bkg*lumFactor;
    double bkg_stat2 = tau_fakes_stat2 + other_bkg_stat2*pow(lumFactor,2);
    double bkg_syst2 = tau_fakes_syst2 + other_bkg_syst2*pow(lumFactor,2);
   
    double num = data_yields - bkg; 

    
    double num_stat = sqrt( bkg_stat2 + data_yields );
   


    double num_syst = sqrt(bkg_syst2);

    double den      = acceptance*lum*lumFactor;  
    double den_stat = acceptance_stat*lum*lumFactor;  
    double den_syst = acceptance_syst*lum*lumFactor;

    double xsec = num/den ; 
    double xsec_stat = getErrorFraction( num,den,num_stat,den_stat);
    double xsec_syst = getErrorFraction( num,den,num_syst,den_syst);

    cout<<endl<<" debug  XSECTION "<<endl;
    cout<<"---------------------------------------------------------------------------------------------"<<endl;
    cout<<" tau fakes   : "<<tau_fakes<<" stat2 : "<<tau_fakes_stat2<<" syst2 : "<<tau_fakes_syst2<<endl;
    cout<<" other bkg   : "<<other_bkg<<" stat2 : "<<other_bkg_stat2<<" syst2 : "<<other_bkg_syst2<<endl;
    cout<<"---------------------------------------------------------------------------------------------"<<endl;
    cout<<" total bkg   : "<<bkg<<" stat2 : "<<bkg_stat2<<" syst2 : "<<bkg_syst2<<endl;
    cout<<"---------------------------------------------------------------------------------------------"<<endl;
    cout<<" data yields : "<<data_yields<<endl;
    cout<<"---------------------------------------------------------------------------------------------"<<endl;
    cout<<" XSEC "<<endl;
    cout<<" num  "<<num<<" +- "<<num_stat<<" +- "<<num_syst<<endl;
    cout<<" den  "<<den<<" +- "<<den_stat<<" +- "<<den_syst<<endl;
    cout<<" acceptance  : "<<acceptance<<" +- "<<acceptance_stat<<" +- "<<acceptance_syst<<endl;
    cout<<" xsec        : "<<xsec<<" +- "<<xsec_stat<<" +- "<<xsec_syst<<endl;



    /////////////////////////////////////////////////////////////////////////////////////////////

    // xsection line /////////////////////////////////////////////////////////
    TString line2(""); line2.Append(" $\\sigma_{l\\tau}$ &"); TString d2("");
    d2=TString(" %8.1f $\\pm$ %8.1f $\\pm$ %8.1f \\\\ \n"); line2.Append(d2);
    const char * data2 = line2.Data();  fprintf(f, data2,xsec, xsec_stat,xsec_syst); 
    /////////////////////////////////////////////////////////////////////////// 

    fprintf(f,"\\hline \n"); 
  }


  fprintf(f,"\\end{tabular}");
  //f(allSysErrors) fprintf(f,"\\end{footnotesize} \n");
  if(incDocument) fprintf(f,"\\end{center} \n");   //COMMENT
  if(incDocument) fprintf(f,"\\end{document} \n"); //COMMENT
  fclose(f);

  if(higgs) doDatacards(data_datacards, tbh_datacards, htb_datacards, sm_datacards, taufakes_datacards, withShapes, withStatShapes, string("PFlow"));


  processedMCFile->Close();
  processedDataFile->Close();

 

  }
  

  void TauDileptonTableBuilder::doDatacards(vector<double> data, vector<vector<double> > tbh, vector<vector<double> > htb, vector<vector<double> > sm, vector<double> fakes, bool withShapes, bool withStatShapes, string tauType){
    
    
    
    
    /// light charged ///     double vHHm80, sErHHm80, jesPHHm80, jesNHHm80, buPHHm80, buNHHm80, trPHHm80, trNHHm80, puHHm80; //yield, stat error, jer uncert. 
    /// light charged ///     double vHHm100, sErHHm100, jesPHHm100, jesNHHm100, buPHHm100, buNHHm100, trPHHm100, trNHHm100, puHHm100;
    /// light charged ///     double vHHm120, sErHHm120, jesPHHm120, jesNHHm120, buPHHm120, buNHHm120, trPHHm120, trNHHm120, puHHm120;
    /// light charged ///     double vHHm140, sErHHm140, jesPHHm140, jesNHHm140, buPHHm140, buNHHm140, trPHHm140, trNHHm140, puHHm140;
    /// light charged ///     double vHHm150, sErHHm150, jesPHHm150, jesNHHm150, buPHHm150, buNHHm150, trPHHm150, trNHHm150, puHHm150;
    /// light charged ///     double vHHm155, sErHHm155, jesPHHm155, jesNHHm155, buPHHm155, buNHHm155, trPHHm155, trNHHm155, puHHm155;
    /// light charged ///     double vHHm160, sErHHm160, jesPHHm160, jesNHHm160, buPHHm160, buNHHm160, trPHHm160, trNHHm160, puHHm160;
    /// light charged ///     double vWHm80, sErWHm80, jesPWHm80, jesNWHm80, buPWHm80, buNWHm80, trPWHm80, trNWHm80, puWHm80; //yield, stat error, jer uncert. 
    /// light charged ///     double vWHm100, sErWHm100, jesPWHm100, jesNWHm100, buPWHm100, buNWHm100, trPWHm100, trNWHm100, puWHm100;
    /// light charged ///     double vWHm120, sErWHm120, jesPWHm120, jesNWHm120, buPWHm120, buNWHm120, trPWHm120, trNWHm120, puWHm120;
    /// light charged ///     double vWHm140, sErWHm140, jesPWHm140, jesNWHm140, buPWHm140, buNWHm140, trPWHm140, trNWHm140, puWHm140;
    /// light charged ///     double vWHm150, sErWHm150, jesPWHm150, jesNWHm150, buPWHm150, buNWHm150, trPWHm150, trNWHm150, puWHm150;
    /// light charged ///     double vWHm155, sErWHm155, jesPWHm155, jesNWHm155, buPWHm155, buNWHm155, trPWHm155, trNWHm155, puWHm155;
    /// light charged ///     double vWHm160, sErWHm160, jesPWHm160, jesNWHm160, buPWHm160, buNWHm160, trPWHm160, trNWHm160, puWHm160;
    

    // TBH
    //yield, stat error, jer uncert. 
    double vTBHm180(tbh[0][0]); double sErTBHm180(tbh[0][1]); double jesPTBHm180(sqrt(tbh[0][2]*tbh[0][2] + tbh[0][4]*tbh[0][4] + tbh[0][6]*tbh[0][6])); double jesNTBHm180(sqrt(tbh[0][3]*tbh[0][3] + tbh[0][5]*tbh[0][5] + tbh[0][7]*tbh[0][7])); double buPTBHm180(tbh[0][8]); double buNTBHm180(tbh[0][9]); double trPTBHm180(tbh[0][13]); double trNTBHm180(tbh[0][13]); double puTBHm180(0); 
    double vTBHm190(tbh[1][0]); double sErTBHm190(tbh[1][1]); double jesPTBHm190(sqrt(tbh[1][2]*tbh[1][2] + tbh[1][4]*tbh[1][4] + tbh[1][6]*tbh[1][6])); double jesNTBHm190(sqrt(tbh[1][3]*tbh[1][3] + tbh[1][5]*tbh[1][5] + tbh[1][7]*tbh[1][7])); double buPTBHm190(tbh[1][8]); double buNTBHm190(tbh[1][9]); double trPTBHm190(tbh[1][13]); double trNTBHm190(tbh[1][13]); double puTBHm190(0);
    double vTBHm200(tbh[2][0]); double sErTBHm200(tbh[2][1]); double jesPTBHm200(sqrt(tbh[2][2]*tbh[2][2] + tbh[2][4]*tbh[2][4] + tbh[2][6]*tbh[2][6])); double jesNTBHm200(sqrt(tbh[2][3]*tbh[2][3] + tbh[2][5]*tbh[2][5] + tbh[2][7]*tbh[2][7])); double buPTBHm200(tbh[2][8]); double buNTBHm200(tbh[2][9]); double trPTBHm200(tbh[2][13]); double trNTBHm200(tbh[2][13]); double puTBHm200(0);
    double vTBHm220(tbh[3][0]); double sErTBHm220(tbh[3][1]); double jesPTBHm220(sqrt(tbh[3][2]*tbh[3][2] + tbh[3][4]*tbh[3][4] + tbh[3][6]*tbh[3][6])); double jesNTBHm220(sqrt(tbh[3][3]*tbh[3][3] + tbh[3][5]*tbh[3][5] + tbh[3][7]*tbh[3][7])); double buPTBHm220(tbh[3][8]); double buNTBHm220(tbh[3][9]); double trPTBHm220(tbh[3][13]); double trNTBHm220(tbh[3][13]); double puTBHm220(0);
    double vTBHm250(tbh[4][0]); double sErTBHm250(tbh[4][1]); double jesPTBHm250(sqrt(tbh[4][2]*tbh[4][2] + tbh[4][4]*tbh[4][4] + tbh[4][6]*tbh[4][6])); double jesNTBHm250(sqrt(tbh[4][3]*tbh[4][3] + tbh[4][5]*tbh[4][5] + tbh[4][7]*tbh[4][7])); double buPTBHm250(tbh[4][8]); double buNTBHm250(tbh[4][9]); double trPTBHm250(tbh[4][13]); double trNTBHm250(tbh[4][13]); double puTBHm250(0);
    double vTBHm300(tbh[5][0]); double sErTBHm300(tbh[5][1]); double jesPTBHm300(sqrt(tbh[5][2]*tbh[5][2] + tbh[5][4]*tbh[5][4] + tbh[5][6]*tbh[5][6])); double jesNTBHm300(sqrt(tbh[5][3]*tbh[5][3] + tbh[5][5]*tbh[5][5] + tbh[5][7]*tbh[5][7])); double buPTBHm300(tbh[5][8]); double buNTBHm300(tbh[5][9]); double trPTBHm300(tbh[5][13]); double trNTBHm300(tbh[5][13]); double puTBHm300(0);


    // HTB
    //yield, stat error, jer uncert. 
    double vHTBm180(htb[0][0]); double sErHTBm180(htb[0][1]); double jesPHTBm180(sqrt(htb[0][2]*htb[0][2] + htb[0][4]*htb[0][4] + htb[0][6]*htb[0][6])); double jesNHTBm180(sqrt(htb[0][3]*htb[0][3] + htb[0][5]*htb[0][5] + htb[0][7]*htb[0][7])); double buPHTBm180(htb[0][8]); double buNHTBm180(htb[0][9]); double trPHTBm180(htb[0][13]); double trNHTBm180(htb[0][13]); double puHTBm180(0); 
    double vHTBm190(htb[1][0]); double sErHTBm190(htb[1][1]); double jesPHTBm190(sqrt(htb[1][2]*htb[1][2] + htb[1][4]*htb[1][4] + htb[1][6]*htb[1][6])); double jesNHTBm190(sqrt(htb[1][3]*htb[1][3] + htb[1][5]*htb[1][5] + htb[1][7]*htb[1][7])); double buPHTBm190(htb[1][8]); double buNHTBm190(htb[1][9]); double trPHTBm190(htb[1][13]); double trNHTBm190(htb[1][13]); double puHTBm190(0);
    double vHTBm200(htb[2][0]); double sErHTBm200(htb[2][1]); double jesPHTBm200(sqrt(htb[2][2]*htb[2][2] + htb[2][4]*htb[2][4] + htb[2][6]*htb[2][6])); double jesNHTBm200(sqrt(htb[2][3]*htb[2][3] + htb[2][5]*htb[2][5] + htb[2][7]*htb[2][7])); double buPHTBm200(htb[2][8]); double buNHTBm200(htb[2][9]); double trPHTBm200(htb[2][13]); double trNHTBm200(htb[2][13]); double puHTBm200(0);
    double vHTBm220(htb[3][0]); double sErHTBm220(htb[3][1]); double jesPHTBm220(sqrt(htb[3][2]*htb[3][2] + htb[3][4]*htb[3][4] + htb[3][6]*htb[3][6])); double jesNHTBm220(sqrt(htb[3][3]*htb[3][3] + htb[3][5]*htb[3][5] + htb[3][7]*htb[3][7])); double buPHTBm220(htb[3][8]); double buNHTBm220(htb[3][9]); double trPHTBm220(htb[3][13]); double trNHTBm220(htb[3][13]); double puHTBm220(0);
    double vHTBm250(htb[4][0]); double sErHTBm250(htb[4][1]); double jesPHTBm250(sqrt(htb[4][2]*htb[4][2] + htb[4][4]*htb[4][4] + htb[4][6]*htb[4][6])); double jesNHTBm250(sqrt(htb[4][3]*htb[4][3] + htb[4][5]*htb[4][5] + htb[4][7]*htb[4][7])); double buPHTBm250(htb[4][8]); double buNHTBm250(htb[4][9]); double trPHTBm250(htb[4][13]); double trNHTBm250(htb[4][13]); double puHTBm250(0);
    double vHTBm300(htb[5][0]); double sErHTBm300(htb[5][1]); double jesPHTBm300(sqrt(htb[5][2]*htb[5][2] + htb[5][4]*htb[5][4] + htb[5][6]*htb[5][6])); double jesNHTBm300(sqrt(htb[5][3]*htb[5][3] + htb[5][5]*htb[5][5] + htb[5][7]*htb[5][7])); double buPHTBm300(htb[5][8]); double buNHTBm300(htb[5][9]); double trPHTBm300(htb[5][13]); double trNHTBm300(htb[5][13]); double puHTBm300(0);


    
    //tau fakes, stat err, sys err
    double tauF    (fakes[0]);
    double sErtauF (fakes[1]);
    double sytauF  (fakes[2]); 

    //tt tau-dilepton, stat err, sys err    
    double tt     (sm[1][0]);							   
    double sErtt  (sm[1][1]); 							   
    double jesPtt (sqrt(sm[1][2]*sm[1][2] + sm[1][4]*sm[1][4] +  sm[1][6]*sm[1][6]));   
    double jesNtt (sqrt(sm[1][3]*sm[1][3] + sm[1][5]*sm[1][5] + sm[1][7]*sm[1][7]));   
    double buPtt  (sm[1][8]);							   
    double buNtt  (sm[1][9]);							   
    double trPtt  (sm[1][13]);							   
    double trNtt  (sm[1][13]);							   
    double putt   (0);                                                                 

    double ttll     (sm[3][0]);							   
    double sErttll  (sm[3][1]); 							
    double jesPttll (sqrt(sm[3][2]*sm[3][2] + sm[3][4]*sm[3][4] + sm[3][6]*sm[3][6]));   
    double jesNttll (sqrt(sm[3][3]*sm[3][3] + sm[3][5]*sm[3][5] + sm[3][7]*sm[3][7]));   
    double buPttll  (sm[3][8]);							   
    double buNttll  (sm[3][9]);							   
    double trPttll  (sm[3][13]);						   
    double trNttll  (sm[3][13]);						   
    double puttll   (0);                                                           

    double Zll     (sm[8][0]);							   
    double sErZll  (sm[8][1]); 							   
    double jesPZll (sqrt(sm[8][2]*sm[8][2] + sm[8][4]*sm[8][4] +  sm[8][6]*sm[8][6]))   ;
    double jesNZll (sqrt(sm[8][3]*sm[8][3] + sm[8][5]*sm[8][5] + sm[8][7]*sm[8][7]))   ;
    double buPZll  (sm[8][8]);							   
    double buNZll  (sm[8][9]);							   
    double trPZll  (sm[8][13]);							   
    double trNZll  (sm[8][13]);							   
    double puZll   (0);                                                            

    double Ztau     (sm[9][0]);							   
    double sErZtau  (sm[9][1]) 							   ;
    double jesPZtau (sqrt(sm[9][2]*sm[9][2] + sm[9][4]*sm[9][4] +  sm[9][6]*sm[9][6]))   ;
    double jesNZtau (sqrt(sm[9][3]*sm[9][3] + sm[9][5]*sm[9][5] + sm[9][7]*sm[9][7]))   ;
    double buPZtau  (sm[9][8])							   ;
    double buNZtau  (sm[9][9])							   ;
    double trPZtau  (sm[9][13])							   ;
    double trNZtau  (sm[9][13])							   ;
    double puZtau   (0)                                                                 ;

    double sTop     (sm[10][0])							   ;
    double sErsTop  (sm[10][1]) 							   ;
    double jesPsTop (sqrt(sm[10][2]*sm[10][2] + sm[10][4]*sm[10][4] +  sm[10][6]*sm[10][6]))   ;
    double jesNsTop (sqrt(sm[10][3]*sm[10][3] + sm[10][5]*sm[10][5] + sm[10][7]*sm[10][7]))   ;
    double buPsTop  (sm[10][8])							   ;
    double buNsTop  (sm[10][9])							   ;
    double trPsTop  (sm[10][13])							   ;
    double trNsTop  (sm[10][13])							   ;
    double pusTop   (0)                                                                 ;

    cout << "SIZE: " << sm.size() << endl;
    double VV     (sm[12][0])							   ;
    double sErVV  (sm[12][1]) 							   ;
    double jesPVV (sqrt(sm[12][2]*sm[12][2] + sm[12][4]*sm[12][4] +  sm[12][6]*sm[12][6]))   ;
    double jesNVV (sqrt(sm[12][3]*sm[12][3] + sm[12][5]*sm[12][5] + sm[12][7]*sm[12][7]))   ;
    double buPVV  (sm[12][8])							   ;
    double buNVV  (sm[12][9])							   ;
    double trPVV  (sm[12][13])							   ;
    double trNVV  (sm[12][13])							   ;
    double puVV   (0)                                                                 ;

    double SM, sErSM, sySM; //SM Total

    double ndata(data[0]);

  
    //make for all mass points
    int HMass[6] = {180,
		    190,
		    200,
		    220,
		    250,
		    300};
    for(int im = 0; im < 6; im++){
      double vTBH, sErTBH, jesPTBH, jesNTBH, buTBH, trTBH, puTBH;
      if(im == 0)     {vTBH = vTBHm180; sErTBH = sErTBHm180; jesPTBH = jesPTBHm180; jesNTBH = jesNTBHm180; buTBH = buPTBHm180; trTBH = trPTBHm180; puTBH = puTBHm180;}
      else if(im == 1){vTBH = vTBHm190; sErTBH = sErTBHm190; jesPTBH = jesPTBHm190; jesNTBH = jesNTBHm190; buTBH = buPTBHm190; trTBH = trPTBHm190; puTBH = puTBHm190;}
      else if(im == 2){vTBH = vTBHm200; sErTBH = sErTBHm200; jesPTBH = jesPTBHm200; jesNTBH = jesNTBHm200; buTBH = buPTBHm200; trTBH = trPTBHm200; puTBH = puTBHm200;}
      else if(im == 3){vTBH = vTBHm220; sErTBH = sErTBHm220; jesPTBH = jesPTBHm220; jesNTBH = jesNTBHm220; buTBH = buPTBHm220; trTBH = trPTBHm220; puTBH = puTBHm220;}
      else if(im == 4){vTBH = vTBHm250; sErTBH = sErTBHm250; jesPTBH = jesPTBHm250; jesNTBH = jesNTBHm250; buTBH = buPTBHm250; trTBH = trPTBHm250; puTBH = puTBHm250;}
      else if(im == 5){vTBH = vTBHm300; sErTBH = sErTBHm300; jesPTBH = jesPTBHm300; jesNTBH = jesNTBHm300; buTBH = buPTBHm300; trTBH = trPTBHm300; puTBH = puTBHm300;}
      
      double vHTB, sErHTB, jesPHTB, jesNHTB, buHTB, trHTB, puHTB;
      if(im == 0)     {vHTB = vHTBm180; sErHTB = sErHTBm180; jesPHTB = jesPHTBm180; jesNHTB = jesNHTBm180; buHTB = buPHTBm180; trHTB = trPHTBm180; puHTB = puHTBm180;}
      else if(im == 1){vHTB = vHTBm190; sErHTB = sErHTBm190; jesPHTB = jesPHTBm190; jesNHTB = jesNHTBm190; buHTB = buPHTBm190; trHTB = trPHTBm190; puHTB = puHTBm190;}
      else if(im == 2){vHTB = vHTBm200; sErHTB = sErHTBm200; jesPHTB = jesPHTBm200; jesNHTB = jesNHTBm200; buHTB = buPHTBm200; trHTB = trPHTBm200; puHTB = puHTBm200;}
      else if(im == 3){vHTB = vHTBm220; sErHTB = sErHTBm220; jesPHTB = jesPHTBm220; jesNHTB = jesNHTBm220; buHTB = buPHTBm220; trHTB = trPHTBm220; puHTB = puHTBm220;}
      else if(im == 4){vHTB = vHTBm250; sErHTB = sErHTBm250; jesPHTB = jesPHTBm250; jesNHTB = jesNHTBm250; buHTB = buPHTBm250; trHTB = trPHTBm250; puHTB = puHTBm250;}
      else if(im == 5){vHTB = vHTBm300; sErHTB = sErHTBm300; jesPHTB = jesPHTBm300; jesNHTB = jesNHTBm300; buHTB = buPHTBm300; trHTB = trPHTBm300; puHTB = puHTBm300;}
      

      // Full Charged Higgs
      if(im == 4){vTBH = vHTBm250+vTBHm250; sErTBH = sErHTBm250+sErTBHm250; jesPTBH = jesPHTBm250+jesPTBHm250; jesNTBH = jesNHTBm250+jesNTBHm250; buTBH = buPHTBm250+buPTBHm250; trTBH = trPHTBm250+trPTBHm250; puTBH = puHTBm250+puTBHm250;}


      int iw = 0;
      string filename("datacard_m");
      if(im == 0)filename = filename+"180_";
      else if(im == 1)filename = filename+"190_";
      else if(im == 2)filename = filename+"200_";
      else if(im == 3)filename = filename+"220_";
      else if(im == 4)filename = filename+"250_";
      else if(im == 5)filename = filename+"300_";
      
      //  int HMass[7] = {80, 100, 120, 140, 150, 155, 160};
      //  for(int im = 0; im < 7; im++){
      //    double vHH, sErHH, jesPHH, jesNHH, buHH, trHH, puHH;
      //    if(im == 0){vHH = vHHm80; sErHH = sErHHm80; jesPHH = jesPHHm80; jesNHH = jesNHHm80; buHH = buPHHm80; trHH = trPHHm80; puHH = puHHm80;}
      //    else if(im == 1){vHH = vHHm100; sErHH = sErHHm100; jesPHH = jesPHHm100; jesNHH = jesNHHm100; buHH = buPHHm100; trHH = trPHHm100; puHH = puHHm100;}
      //    else if(im == 2){vHH = vHHm120; sErHH = sErHHm120; jesPHH = jesPHHm120; jesNHH = jesNHHm120; buHH = buPHHm120; trHH = trPHHm120; puHH = puHHm120;}
      //    else if(im == 3){vHH = vHHm140; sErHH = sErHHm140; jesPHH = jesPHHm140; jesNHH = jesNHHm140; buHH = buPHHm140; trHH = trPHHm140; puHH = puHHm140;}
      //    else if(im == 4){vHH = vHHm150; sErHH = sErHHm150; jesPHH = jesPHHm150; jesNHH = jesNHHm150; buHH = buPHHm150; trHH = trPHHm150; puHH = puHHm150;}
      //    else if(im == 5){vHH = vHHm155; sErHH = sErHHm155; jesPHH = jesPHHm155; jesNHH = jesNHHm155; buHH = buPHHm155; trHH = trPHHm155; puHH = puHHm155;}
      //    else if(im == 6){vHH = vHHm160; sErHH = sErHHm160; jesPHH = jesPHHm160; jesNHH = jesNHHm160; buHH = buPHHm160; trHH = trPHHm160; puHH = puHHm160;}
      //
      //    double vWH, sErWH, jesPWH, jesNWH, buWH, trWH, puWH;
      //    if(im == 0){vWH = vWHm80; sErWH = sErWHm80; jesPWH = jesPWHm80; jesNWH = jesNWHm80; buWH = buPWHm80; trWH = trPWHm80; puWH = puWHm80;}
      //    else if(im == 1){vWH = vWHm100; sErWH = sErWHm100; jesPWH = jesPWHm100; jesNWH = jesNWHm100; buWH = buPWHm100; trWH = trPWHm100; puWH = puWHm100;}
      //    else if(im == 2){vWH = vWHm120; sErWH = sErWHm120; jesPWH = jesPWHm120; jesNWH = jesNWHm120; buWH = buPWHm120; trWH = trPWHm120; puWH = puWHm120;}
      //    else if(im == 3){vWH = vWHm140; sErWH = sErWHm140; jesPWH = jesPWHm140; jesNWH = jesNWHm140; buWH = buPWHm140; trWH = trPWHm140; puWH = puWHm140;}
      //    else if(im == 4){vWH = vWHm150; sErWH = sErWHm150; jesPWH = jesPWHm150; jesNWH = jesNWHm150; buWH = buPWHm150; trWH = trPWHm150; puWH = puWHm150;}
      //    else if(im == 5){vWH = vWHm155; sErWH = sErWHm155; jesPWH = jesPWHm155; jesNWH = jesNWHm155; buWH = buPWHm155; trWH = trPWHm155; puWH = puWHm155;}
      //    else if(im == 6){vWH = vWHm160; sErWH = sErWHm160; jesPWH = jesPWHm160; jesNWH = jesNWHm160; buWH = buPWHm160; trWH = trPWHm160; puWH = puWHm160;}
      //    
      //    int iw = 0;
      //    string filename("datacard_m");
      //    if(im == 0)filename = filename+"80_";
      //    else if(im == 1)filename = filename+"100_";
      //    else if(im == 2)filename = filename+"120_";
      //    else if(im == 3)filename = filename+"140_";
      //    else if(im == 4)filename = filename+"150_";
      //    else if(im == 5)filename = filename+"155_";
      //    else if(im == 6)filename = filename+"160_";
    
      
      double nobs = ndata;
      filename = "datacards/"+filename + tauType + ".txt";
      char * file = new char[filename.size() + 1];
      strcpy(file, filename.c_str());
      
      //cout<< fixed << showpoint <<setprecision(2);
      time_t secs=time(0);
      tm *t=localtime(&secs);
      
      ofstream outfile;
      outfile.open(file);
      outfile<< fixed << showpoint <<setprecision(3);
      outfile<<"Data: "<<t->tm_mday<<"/"<<t->tm_mon+1<<"/"<<t->tm_year+1900<<endl;
      //outfile<<"Description: H+, mu-tau, mass "<<HMass[im]<<" GeV, lumi=2140 pb-1"<<endl;
      outfile<<"Description: H+, e-tau, mass "<<HMass[im]<<" GeV, lumi=" << commondefinitions::LUM_<<" pb-1"<<endl;
      outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
      outfile<<"imax   1  number of channels"<<endl;
      outfile<<"jmax   *  number of backgrounds"<<endl;
      outfile<<"kmax   *  number of nuisance parameters"<<endl;
      if(withShapes){
	outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	outfile<<"shapes * * shapes_m"<<HMass[im]<<"_rc_t.root $PROCESS $PROCESS_$SYSTEMATIC"<<endl;
      }
      outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
      outfile<<"bin a"<<endl;
      outfile<<"observation    "<<nobs<<endl;
      outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
      
      
      outfile<<"bin              a          a          a           a           a           a           a         a         "<<endl;
      outfile<<"process         TBH    tt_ltau    tt_ll     tau_fake    Z_eemumu   Z_tautau   singleTop  di_boson"<<endl;
      outfile<<"process          0          1          2           3           4           5           6         7         "<<endl;
      outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
      outfile<<"rate      "<<vTBH<<setw(10)<<tt<<setw(10)<<ttll<<setw(10)<<tauF<<setw(10)<<Zll<<setw(10)<<Ztau<<setw(10)<<sTop<<setw(10)<<VV<<endl;//"       Projected event rates"<<endl; 
      outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
      outfile<<endl;
      outfile<<" tauMetTrg      lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tau+MET trg"<<endl;
      outfile<<" eJetMHTTrg     lnN"<<setw(7)<<1+trTBH/vTBH<<setw(10)<<1+trPtt/tt<<setw(10)<<1+trPttll/ttll<<setw(10)<<1.00<<setw(10)<<1+trPZll/Zll<<setw(10)<<1+trPZtau/Ztau<<setw(10)<<1+trPsTop/sTop<<setw(10)<<1+trPVV/VV<<endl;//"    e+jet+MHT trig"<<endl;
      outfile<<" tauId  lnN"<<setw(7)<<1.06<<setw(10)<<1.06<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.06<<setw(10)<<1.06<<setw(10)<<1.06<<endl;//"    tau ID"<<endl;
      outfile<<" jetTauMisId    lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.15<<setw(10)<<1.00<<setw(10)<<1.15<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    jet-> tau miss ID"<<endl;
      outfile<<" fakesSyst      lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sytauF/tauF<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tau-fakes syst"<<endl;
      outfile<<" fakesStat      lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErtauF/tauF<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tau-fakes stat"<<endl;
      outfile<<" jesJerMet      lnN"<<setw(7)<<1-jesNTBH/vTBH<<"/"<<1+jesPTBH/vTBH<<setw(7)<<1-jesNtt/tt<<"/"<<1+jesPtt/tt<<setw(7)<<1-jesNttll/ttll<<"/"<<1+jesPttll/ttll<<setw(7)<<1.00<<setw(10)<<1-jesNZll/Zll<<"/"<<1+jesPZll/Zll<<setw(7)<<1-jesNZtau/Ztau<<"/"<<1+jesPZtau/Ztau<<setw(7)<<1-jesNsTop/sTop<<"/"<<1+jesPsTop/sTop<<setw(7)<<1-jesNVV/VV<<"/"<<1+jesPVV/VV<<endl;//"    JES/JER/MET scale"<<endl;
      outfile<<" leptEff        lnN"<<setw(7)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.00<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<endl;//"    Lepton Efficiency"<<endl;
      //    outfile<<" leptonVeto     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    Lepton Veto"<<endl;
      outfile<<" btagging       lnN"<<setw(7)<<1+buTBH/vTBH<<setw(10)<<1+buPtt/tt<<setw(10)<<1+buPttll/ttll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+buPsTop/sTop<<setw(10)<<1.00<<endl;//"    b-tagging"<<endl;
      outfile<<"bmistagging     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+buPZll/Zll<<setw(10)<<1+buPZtau/Ztau<<setw(10)<<1.00<<setw(10)<<1+buPVV/VV<<endl;//"    b-mis-tagging"<<endl;
      //    outfile<<17<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    HH     MC stat"<<endl;
      //    outfile<<17<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    HW     MC stat"<<endl;
      outfile<<" tbhStatistics     lnN"<<setw(7)<<1+sErTBH/vTBH<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    TBH     MC stat"<<endl;
      // outfile<<19<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    embedding MC stat"<<endl;
      outfile<<" ttltauStatistics  lnN"<<setw(7)<<1.00<<setw(10)<<1+sErtt/tt<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tt->ltau MC stat "<<endl;
      outfile<<" ttllStatistics    lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErttll/ttll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tt->ll  MC stat"<<endl;
      //    outfile<<22<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    tt->jjtau MC stat"<<endl;
      outfile<<" zeemumuStatistics     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErZll/Zll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<endl;//"    Z->ee,mumu MC stat"<<endl;
      outfile<<" ztautauStatistics     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErZtau/Ztau<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    Z->tautau  MC stat"<<endl;
      outfile<<" wjetsStatistics       lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    W+jets MC stat"<<endl;
      outfile<<" singleTopStatistics   lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErsTop/sTop<<setw(10)<<1.00<<endl;//"    single top MC stat"<<endl;
      outfile<<" dibosonStatistics     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErVV/VV<<endl;//"    diboson MC stat"<<endl;
      outfile<<" ttbarCrossSection     lnN"<<setw(7)<<0.90<<"/"<<1.07<<setw(10)<<0.90<<"/"<<1.07<<setw(7)<<0.90<<"/"<<1.07<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    ttbar cross-section"<<endl;
      outfile<<" wjetsCrossSection     lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    W+jets cross section"<<endl;
      outfile<<" singletopCrossSection lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.08<<setw(10)<<1.00<<setw(10)<<endl;//"    singleTop cross-section"<<endl; 
      outfile<<" zllCrossSection       lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.04<<setw(10)<<1.04<<setw(10)<<1.00<<setw(10)<<1.00<<endl;//"    Z->ll cross section"<<endl;
      outfile<<" dibosonCrossSection   lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.04<<endl;//"    diboson cross-section"<<endl;
      outfile<<" lumiErr                lnN"<<setw(7)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.00<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<endl;//"    Luminosity Error"<<endl;
      outfile<<" pileupErr               lnN"<<setw(7)<<1+puTBH/vTBH<<setw(10)<<1+putt/tt<<setw(10)<<1+puttll/ttll<<setw(10)<<1.00<<setw(10)<<1+puZll/Zll<<setw(10)<<1+puZtau/Ztau<<setw(10)<<1+pusTop/sTop<<setw(10)<<1+puVV/VV<<endl;//"    pileup"<<endl;
      if(withShapes){
	outfile<<"jes           "<<"  shape        1        1        1           1           1           1           1           1            1     "<<endl;  //        JES_effect_on_shape                             
	outfile<<"met           "<<"  shape        1        1        1           1           1           1           1           1            1     "<<endl; //        MET_effect_on_shape                             
	outfile<<"jer           "<<"  shape        1        1        1           1           1           1           1           1            1     "<<endl;  //       JER_effect_on_shape                
	if(!withStatShapes) outfile<<"#";outfile<<"TBH_Stat      "<<"  shape        1        -        -           -           -           -           -           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"tt_ltau_Stat  "<<"  shape        -        -        1           -           -           -           -           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"tt_ll_Stat    "<<"  shape        -        -        -           1           -           -           -           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"tau_fake_Stat "<<"  shape        -        -        -           -           1           -           -           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"Z_eemumu_Stat "<<"  shape        -        -        -           -           -           1           -           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"Z_tautau_Stat "<<"  shape        -        -        -           -           -           -           1           -            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"singleTop_Stat"<<"  shape        -        -        -           -           -           -           -           1            -     "<<endl;
	if(!withStatShapes) outfile<<"#";outfile<<"di_boson_Stat "<<"  shape        -        -        -           -           -           -           -           -            1     "<<endl;
	//             
	//outfile<< HH3_Stat        shape  1         -        -           -           -           -           -           -            -
	//outfile<< WH3_Stat        shape  -         1        -           -           -           -           -           -            -
      }
      outfile.close();
      
      //}
    }
    return ;
    
  }
  
  /*

  // NOTE : THIS PEACE OF CODE CAN BE USED FOR


  vector<int> colors;        
  vector<int> markerStyles;   
  markerStyles.push_back(21); colors.push_back(1); //black  square
  markerStyles.push_back(21); colors.push_back(2); //red    square
  markerStyles.push_back(21); colors.push_back(4); //blue   square

  markerStyles.push_back(20); colors.push_back(1); //black  circle
  markerStyles.push_back(20); colors.push_back(2); //red    circle
  markerStyles.push_back(20); colors.push_back(4); //blue   circle

  markerStyles.push_back(22); colors.push_back(1); //black triangle



  //Muon case
  if(numbHWSamples || numbHHSamples){

    //get plot with ratios
    //TH1F * efficiencies_hh = new TH1F ("efficiencies_hh", "Cut Flow; Step;N",3,0,3);
    //efficiencies_hh->GetXaxis()->SetBinLabel(1, "btag/pre-selection");
    //efficiencies_hh->GetXaxis()->SetBinLabel(2, "1tau/btag");
    //efficiencies_hh->GetXaxis()->SetBinLabel(3, "OS/1tau");

    vector< TGraphAsymmErrors * > histograms_hw, histograms_hh;

    for( int k=0; k<numbHWSamples ; k++){ 

      char name[100]; sprintf(name,"hh_%d",k); 
      TString histoName(name);
      
      TGraphAsymmErrors * hist = new TGraphAsymmErrors(3);
      hist->SetName(histoName);
      hist->GetXaxis()->SetRangeUser(0,4);
      hist->GetYaxis()->SetRangeUser(0.6,1);

      //hist->GetXaxis()->SetBinLabel(1, "btag/pre-selection");
      //hist->GetXaxis()->SetBinLabel(2, "1tau/btag");
      //hist->GetXaxis()->SetBinLabel(3, "OS/1tau");
      hist->SetMarkerColor( colors[k] );
      hist->SetMarkerStyle( markerStyles[k] );
      hist->SetLineColor(colors[k]);
      hist->SetLineStyle(3);

      double r1 = hw_m[k][COL3]/hw_m[k][COL2];      double errY1 =    getErrorFraction( hw_m[k][COL3],hw_m[k][COL2],hw_m_err[k][COL3],hw_m_err[k][COL2]);
       // for tau id we apply a factor on eff of 6
      double r2 = 6*(hw_m[k][COL4]/hw_m[k][COL3]);  double errY2 = 6*(getErrorFraction( hw_m[k][COL4],hw_m[k][COL3],hw_m_err[k][COL4],hw_m_err[k][COL3]));
      double r3 = hw_m[k][COL5]/hw_m[k][COL4];      double errY3 =    getErrorFraction( hw_m[k][COL5],hw_m[k][COL4],hw_m_err[k][COL5],hw_m_err[k][COL4]);

      if( hw_m[k][COL2] ) hist->SetPoint(0,1,r1 ); hist->SetPointError(0,0,0,errY1,errY1);
      if( hw_m[k][COL3] ) hist->SetPoint(1,2,r2);  hist->SetPointError(1,0,0,errY2,errY2);
      if( hw_m[k][COL4] ) hist->SetPoint(2,3,r3 ); hist->SetPointError(2,0,0,errY3,errY3);

      histograms_hw.push_back(hist);

    }
    
    TCanvas * c_hw = new TCanvas("c_hw", "title",4,46,600,600); c_hw->cd();    
    TLegend * leg_hw = new TLegend(0.6895973,0.2045455,0.8708054,0.5017483,NULL,"brNDC");
    leg_hw->SetFillColor(0);
    leg_hw->SetBorderSize(0);
  
    for(int i=0;i<histograms_hw.size();i++){

      TGraphAsymmErrors * h = histograms_hw[i];
      h->SetTitle("");
      h->GetXaxis()->SetTitle("selection step i ");
      //h->GetYaxis()->SetTitle("step[i]/step[i-1]");
      TLegendEntry * entry = leg_hw->AddEntry(h->GetName(),hwEffHitoNames[i],"LP");  
      entry->SetLineColor(colors[i]);
      entry->SetLineStyle(3);
      entry->SetMarkerColor(colors[i]);
      entry->SetMarkerStyle(markerStyles[i]);

      if( i==0 ) h->Draw("ALP");
      else h->Draw("LP");
    }  
    leg_hw->Draw();
    c_hw->SaveAs("eff_hw.png");



   for( int k=0; k<numbHHSamples ; k++){ 

      char name[100]; sprintf(name,"hh_%d",k); 
      TString histoName(name);
      
      TGraphAsymmErrors * hist = new TGraphAsymmErrors(3);
      hist->SetName(histoName);
      hist->GetXaxis()->SetRangeUser(0,4);
      hist->GetYaxis()->SetRangeUser(0.6,1);

      //hist->GetXaxis()->SetBinLabel(1, "btag/pre-selection");
      //hist->GetXaxis()->SetBinLabel(2, "1tau/btag");
      //hist->GetXaxis()->SetBinLabel(3, "OS/1tau");
      hist->SetMarkerColor( colors[k] );
      hist->SetMarkerStyle( markerStyles[k] );
      hist->SetLineColor(colors[k]);
      hist->SetLineStyle(3);

      double r1 = hh_m[k][COL3]/hh_m[k][COL2];      double errY1 =    getErrorFraction( hh_m[k][COL3],hh_m[k][COL2],hh_m_err[k][COL3],hh_m_err[k][COL2]);
       // for tau id we apply a factor on eff of 6
      double r2 = 4*(hh_m[k][COL4]/hh_m[k][COL3]);  double errY2 = 4*(getErrorFraction( hh_m[k][COL4],hh_m[k][COL3],hh_m_err[k][COL4],hh_m_err[k][COL3]));
      double r3 = hh_m[k][COL5]/hh_m[k][COL4];      double errY3 =    getErrorFraction( hw_m[k][COL5],hh_m[k][COL4],hh_m_err[k][COL5],hh_m_err[k][COL4]);

      if( hh_m[k][COL2] ) hist->SetPoint(0,1,r1 ); hist->SetPointError(0,0,0,errY1,errY1);
      if( hh_m[k][COL3] ) hist->SetPoint(1,2,r2);  hist->SetPointError(1,0,0,errY2,errY2);
      if( hh_m[k][COL4] ) hist->SetPoint(2,3,r3 ); hist->SetPointError(2,0,0,errY3,errY3);

      histograms_hh.push_back(hist);

    }
    
    TCanvas * c_hh = new TCanvas("c_hh", "title",4,46,600,600); c_hh->cd();    
    TLegend * leg_hh = new TLegend(0.6895973,0.2045455,0.8708054,0.5017483,NULL,"brNDC");
    leg_hh->SetFillColor(0);
    leg_hh->SetBorderSize(0);
  
    for(int i=0;i<histograms_hh.size();i++){

      TGraphAsymmErrors * h = histograms_hh[i];
      h->SetTitle("");
      h->GetXaxis()->SetTitle("selection step i ");
      //h->GetYaxis()->SetTitle("step[i]/step[i-1]");
      TLegendEntry * entry = leg_hh->AddEntry(h->GetName(),hhEffHitoNames[i],"LP");  
      entry->SetLineColor(colors[i]);
      entry->SetLineStyle(3);
      entry->SetMarkerColor(colors[i]);
      entry->SetMarkerStyle(markerStyles[i]);

      if( i==0 ) h->Draw("ALP");
      else h->Draw("LP");
    }  
    leg_hh->Draw();
    c_hh->SaveAs("eff_hw.png");

  }

*/



// Results Muon

/*

root [36]   analyzer.summaryTable(notDetailed, sm, 2,false,false,false);

 DATA NAME : PFlow/Yields/btag1/lep_tau/cutflow_yields_4


 systematics on tau-dilepton : 
 jes        : 8.20766
 unc        : 2.47702
 btag       : 7.69411
 tau id     : 9.36181
 ttbar xsec : 0
 lum        : 0
 lepton eff : 3.1206
 theory unc : 6.24121

 tau dilepton -> 156.03 stat +/- 3.86726 syst +/- 16.4023

 Acceptance   : 0.000442714 +- 1.09753e-05 +- 4.65392e-05

 tau_dilepton : 156.03

 Numerator for acceptance : 156.03

 Initial is               : 3.67694e+06


 scale                    : 0.0958515


 debug  XSECTION 
---------------------------------------------------------------------------------------------
 tau fakes   : 88.5 stat2 : 0 syst2 : 195.85
 other bkg   : 53.3716 stat2 : 19.585 syst2 : 48.68
---------------------------------------------------------------------------------------------
 total bkg   : 141.872 stat2 : 19.585 syst2 : 244.53
---------------------------------------------------------------------------------------------
 data yields : 288
---------------------------------------------------------------------------------------------
 XSEC 
 num  146.128 +- 17.5381 +- 15.6375
 den  0.945637 +- 0.0234432 +- 0.0994078
 acceptance  : 0.000442714 +- 1.09753e-05 +- 4.65392e-05
 xsec        : 154.529 +- 18.9379 +- 23.1805
root [37] 


*/



/*

 systematics on tau-dilepton : 
 jes        : 4.15531
 unc        : 0.477044
 btag       : 4.80082
 tau id     : 5.73621
 ttbar xsec : 0
 lum        : 0
 lepton eff : 2.86811
 theory unc : 3.82414

 tau dilepton -> 95.6035 stat +/- 2.86111 syst +/- 9.81307

 Acceptance   : 0.000303664 +- 9.08907e-06 +- 3.11691e-05

 tau_dilepton : 95.6035

 Numerator for acceptance : 95.6035

 Initial is               : 3.67694e+06


 scale                    : 0.0856238


 debug  XSECTION 
---------------------------------------------------------------------------------------------
 tau fakes   : 54.3 stat2 : 0 syst2 : 106.57
 other bkg   : 37.8657 stat2 : 13.788 syst2 : 13.5051
---------------------------------------------------------------------------------------------
 total bkg   : 92.1657 stat2 : 13.788 syst2 : 120.075
---------------------------------------------------------------------------------------------
 data yields : 176
---------------------------------------------------------------------------------------------
 XSEC 
 num  83.8343 +- 13.7764 +- 10.9579
 den  0.579421 +- 0.0173429 +- 0.0594738
 acceptance  : 0.000303664 +- 9.08907e-06 +- 3.11691e-05
 xsec        : 144.686 +- 24.1672 +- 24.046

*/



/* NEW LUMINOSITY ***************************************************************************************************************************************************************

 debug  XSECTION XSEC
---------------------------------------------------------------------------------------------
 tau fakes   : 88.5 stat2 : 0 syst2 : 195.85
 other bkg   : 55.0457 stat2 : 20.8329 syst2 : 46.6174
---------------------------------------------------------------------------------------------
 total bkg   : 143.546 stat2 : 20.8329 syst2 : 242.467
---------------------------------------------------------------------------------------------
 data yields : 288
---------------------------------------------------------------------------------------------
 XSEC 
 num  144.454 +- 17.5736 +- 15.5714
 den  0.975299 +- 0.0241785 +- 0.100476
 acceptance  : 0.000442714 +- 1.09753e-05 +- 4.56086e-05
 xsec        : 148.113 +- 18.389 +- 22.0846


 debug  XSECTION 
---------------------------------------------------------------------------------------------
 tau fakes   : 54.3 stat2 : 0 syst2 : 106.57
 other bkg   : 38.9221 stat2 : 14.568 syst2 : 14.1407
---------------------------------------------------------------------------------------------
 total bkg   : 93.2221 stat2 : 14.568 syst2 : 120.711
---------------------------------------------------------------------------------------------
 data yields : 176
---------------------------------------------------------------------------------------------
 XSEC 
 num  82.7779 +- 13.8046 +- 10.9868
 den  0.595585 +- 0.0178267 +- 0.0599399
 acceptance  : 0.000303664 +- 9.08907e-06 +- 3.05608e-05
 xsec        : 138.986 +- 23.5486 +- 23.1506



*/
/*
    ofstream outfile;
    outfile.open(file);
    outfile<< fixed << showpoint <<setprecision(3);
    outfile<<"Data: "<<t->tm_mday<<"/"<<t->tm_mon+1<<"/"<<t->tm_year+1900<<endl;
    //outfile<<"Description: H+, mu-tau, mass "<<HMass[im]<<" GeV, lumi=2140 pb-1"<<endl;
    outfile<<"Description: H+, e-tau, mass "<<HMass[im]<<" GeV, lumi=1900 pb-1"<<endl;
    outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    outfile<<"imax   1  number of channels"<<endl;
    outfile<<"jmax   *  number of backgrounds"<<endl;
    outfile<<"kmax  34  number of nuisance parameters"<<endl;
    outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    outfile<<"Observation    "<<nobs<<endl;
    outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    outfile<<"bin       1         1          1          1           1           1           1           1         1         "<<endl;
    outfile<<"process   HH       WH    tt->ltau    tt~->ll     tau-fake    Z->ee,mumu   Z->tautau   singleTop  di-boson"<<endl;
    outfile<<"process   -1        0          1          2           3           4           5           6         7         "<<endl;
    outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    outfile<<"rate      "<<vHH<<setw(10)<<vWH<<setw(10)<<tt<<setw(10)<<ttll<<setw(10)<<tauF<<setw(10)<<Zll<<setw(10)<<Ztau<<setw(10)<<sTop<<setw(10)<<VV<<"       Projected event rates"<<endl; 
    outfile<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    outfile<<endl;
    outfile<<" "<<1<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tau+MET trg"<<endl;
    outfile<<" "<<2<<"  lnN"<<setw(7)<<1+trHH/vHH<<setw(10)<<1+trWH/vWH<<setw(10)<<1+trPtt/tt<<setw(10)<<1+trPttll/ttll<<setw(10)<<1.00<<setw(10)<<1+trPZll/Zll<<setw(10)<<1+trPZtau/Ztau<<setw(10)<<1+trPsTop/sTop<<setw(10)<<1+trPVV/VV<<"    e+jet+MHT trig"<<endl;
    outfile<<" "<<3<<"  lnN"<<setw(7)<<1.06<<setw(10)<<1.06<<setw(10)<<1.06<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.06<<setw(10)<<1.06<<setw(10)<<1.06<<"    tau ID"<<endl;
    outfile<<" "<<4<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.15<<setw(10)<<1.00<<setw(10)<<1.15<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    jet-> tau miss ID"<<endl;
    outfile<<" "<<5<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sytauF/tauF<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tau-fakes syst"<<endl;
    outfile<<" "<<6<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErtauF/tauF<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tau-fakes stat"<<endl;
    outfile<<" "<<7<<"  lnN"<<setw(7)<<1-jesNHH/vHH<<"/"<<1+jesPHH/vHH<<setw(7)<<1-jesNWH/vWH<<"/"<<1+jesPWH/vWH<<setw(7)<<1-jesNtt/tt<<"/"<<1+jesPtt/tt<<setw(7)<<1-jesNttll/ttll<<"/"<<1+jesPttll/ttll<<setw(7)<<1.00<<setw(10)<<1-jesNZll/Zll<<"/"<<1+jesPZll/Zll<<setw(7)<<1-jesNZtau/Ztau<<"/"<<1+jesPZtau/Ztau<<setw(7)<<1-jesNsTop/sTop<<"/"<<1+jesPsTop/sTop<<setw(7)<<1-jesNVV/VV<<"/"<<1+jesPVV/VV<<"    JES/JER/MET scale"<<endl;
    outfile<<" "<<8<<"  lnN"<<setw(7)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.00<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<setw(10)<<1.02<<"    Lepton Efficiency"<<endl;
    outfile<<" "<<9<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Lepton Veto"<<endl;
    outfile<<10<<"  lnN"<<setw(7)<<1+buHH/vHH<<setw(10)<<1+buWH/vWH<<setw(10)<<1+buPtt/tt<<setw(10)<<1+buPttll/ttll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+buPsTop/sTop<<setw(10)<<1.00<<"    b-tagging"<<endl;
    outfile<<11<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+buPZll/Zll<<setw(10)<<1+buPZtau/Ztau<<setw(10)<<1.00<<setw(10)<<1+buPVV/VV<<"    b-mis-tagging"<<endl;
    outfile<<12<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Hadr: QCD stat."<<endl;
    outfile<<13<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Hadr: QCD syst."<<endl;
    outfile<<14<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Hadr: Non-QCD-stat"<<endl;
    outfile<<15<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Hadr: Non-QCD W->tau->mu"<<endl;
    outfile<<16<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    Hadr: Muon selection"<<endl;
    outfile<<17<<"  lnN"<<setw(7)<<1+sErHH/vHH<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    HH     MC stat"<<endl;
    outfile<<18<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1+sErWH/vWH<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    WH     MC stat"<<endl;
    outfile<<19<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    embedding MC stat"<<endl;
    outfile<<20<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErtt/tt<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tt->ltau MC stat "<<endl;
    outfile<<21<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErttll/ttll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tt->ll  MC stat"<<endl;
    outfile<<22<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    tt->jjtau MC stat"<<endl;
    outfile<<23<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErZll/Zll<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<"    Z->ee,mumu MC stat"<<endl;
    outfile<<24<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErZtau/Ztau<<setw(10)<<1.00<<setw(10)<<1.00<<"    Z->tautau  MC stat"<<endl;
    outfile<<25<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    W+jets MC stat"<<endl;
    outfile<<26<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErsTop/sTop<<setw(10)<<1.00<<"    single top MC stat"<<endl;
    outfile<<27<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1+sErVV/VV<<"    diboson MC stat"<<endl;
    outfile<<28<<"  lnN"<<setw(7)<<0.90<<"/"<<1.07<<setw(7)<<0.90<<"/"<<1.07<<setw(10)<<0.90<<"/"<<1.07<<setw(7)<<0.90<<"/"<<1.07<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    ttbar cross-section"<<endl;
    outfile<<29<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<"    W+jets cross section"<<endl;
    outfile<<30<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.08<<setw(10)<<1.00<<setw(10)<<"    singleTop cross-section"<<endl; 
    outfile<<31<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.04<<setw(10)<<1.04<<setw(10)<<1.00<<setw(10)<<1.00<<"    Z->ll cross section"<<endl;
    outfile<<32<<"  lnN"<<setw(7)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.00<<setw(10)<<1.04<<"    diboson cross-section"<<endl;
    outfile<<33<<"  lnN"<<setw(7)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.00<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<setw(10)<<1.045<<"    Luminosity Error"<<endl;
    outfile<<34<<"  lnN"<<setw(7)<<1+puHH/vHH<<setw(10)<<1+puWH/vWH<<setw(10)<<1+putt/tt<<setw(10)<<1+puttll/ttll<<setw(10)<<1.00<<setw(10)<<1+puZll/Zll<<setw(10)<<1+puZtau/Ztau<<setw(10)<<1+pusTop/sTop<<setw(10)<<1+puVV/VV<<"    pileup"<<endl;

*/

} // namespace tableutils 
