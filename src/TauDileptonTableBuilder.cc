#include "LIP/TopTaus/interface/TauDileptonTableBuilder.hh"

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "LIP/TopTaus/interface/CommonDefinitions.hh"
// System headers
#include <cmath>
#include <map>

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
  void TauDileptonTableBuilder::mcTable( int detailed, bool includeSoverB , bool printAllErrors, bool higgs, TString key, TString name, bool systset1, bool systset2, bool systset3, bool systset4){ 
    
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
    double fhh(1),fhw(1); // this is for the datacards (so that they are normalized to 165 (xsecttbar)
    
    
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
    bool includeSM,hh_contrib,hw_contrib(true);
  ////////////////////////////////////////////
    
    
    if(higgs){ hh_contrib = true;  hw_contrib = true; includeSM = false; }
  else     { hh_contrib = false; hw_contrib = false; includeSM = true; }
    
    
    vector<TString> hwEffHitoNames;
    hwEffHitoNames.push_back("WH80");  hwEffHitoNames.push_back("WH100"); hwEffHitoNames.push_back("WH120"); hwEffHitoNames.push_back("WH140"); 
    hwEffHitoNames.push_back("WH150"); hwEffHitoNames.push_back("WH155"); hwEffHitoNames.push_back("WH160"); 
    
    vector<TString> hhEffHitoNames;
    hhEffHitoNames.push_back("HH80");  hhEffHitoNames.push_back("HH100"); hhEffHitoNames.push_back("HH120"); hhEffHitoNames.push_back("HH140"); 
    hhEffHitoNames.push_back("HH150"); hhEffHitoNames.push_back("HH155");
    hhEffHitoNames.push_back("HH160"); 
    
    vector<TString> hhFiles; vector<TString> hwFiles;
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

  TString pFile(outputArea_+TString("nomt-NEWJ-V1-mc-EL-20GeV/out-ttbar_539.root"));       TFile * processedFile     = TFile::Open(pFile); 
  TString dFile(outputArea_+TString("nomt-NEWJ-V1-data-EL-20GeV/out-Electron_1.root"));   TFile * processedDataFile = TFile::Open(dFile); 

  vector<TFile *> processedHWFiles; vector<TFile *> processedHHFiles;
  for(uint i = 0; i<hhFiles.size(); i++ ){  processedHHFiles.push_back( TFile::Open( hhFiles[i] ) ); }
  for(uint i = 0; i<hwFiles.size(); i++ ){  processedHWFiles.push_back( TFile::Open( hwFiles[i] ) ); }

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
  TString hName_tesplus;     hName_tesplus  = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_tesplus_")  + leptauStringMC; 
  TString hName_tesminus;    hName_tesminus = key+TString("/YieldsMC")+TString("/cutflow_yields_mc_tesminus_") + leptauStringMC;
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
  TString hName_higgs_e_tau_tesplus;   hName_higgs_e_tau_tesplus  = key+TString("/Yields") + TString("/e_tau/cutflow_yields_tesplus_")  + etauString; 
  TString hName_higgs_e_tau_tesminus;  hName_higgs_e_tau_tesminus = key+TString("/Yields") + TString("/e_tau/cutflow_yields_tesminus_") + etauString;
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
  TString hName_higgs_m_tau_tesplus;   hName_higgs_m_tau_tesplus  = key+TString("/Yields") + TString("/m_tau/cutflow_yields_tesplus_")  + mtauString;  
  TString hName_higgs_m_tau_tesminus;  hName_higgs_m_tau_tesminus = key+TString("/Yields") + TString("/m_tau/cutflow_yields_tesminus_") + mtauString;
  TString hName_higgs_m_tau_bplus;     hName_higgs_m_tau_bplus    = key+TString("/Yields") + TString("/m_tau/cutflow_yields_bplus_")    + mtauString;  
  TString hName_higgs_m_tau_bminus;    hName_higgs_m_tau_bminus   = key+TString("/Yields") + TString("/m_tau/cutflow_yields_bminus_")   + mtauString;
  TString hName_higgs_m_tau_trigger;   hName_higgs_m_tau_trigger  = key+TString("/Yields") + TString("/m_tau/cutflow_triggErr_")        + mtauString;

  // MC split channel plots
  TH2D * h(0); TH2D * h_plus(0); TH2D * h_minus(0); 
  TH2D * h_uncplus(0); TH2D * h_uncminus(0); TH2D * h_jerplus(0); TH2D * h_jerminus(0); TH2D * h_tesplus(0); TH2D * h_tesminus(0); TH2D * h_bplus(0); TH2D * h_bminus(0);TH2D * h_trigger(0);

  // data plots
  TH1D * hData(0) ; TH1D * hData_e(0); TH1D * hData_m(0);
  // MC e/m channels
  TH1D * h_e(0); TH1D * h_m(0);

  vector< TH1D * > h_hh_e, h_hh_e_plus, h_hh_e_minus, h_hh_e_uncplus, h_hh_e_uncminus, h_hh_e_jerplus, h_hh_e_jerminus, h_hh_e_tesplus, h_hh_e_tesminus, h_hh_e_bplus, h_hh_e_bminus, h_hh_e_trigger;
  vector< TH1D * > h_hh_m, h_hh_m_plus, h_hh_m_minus, h_hh_m_uncplus, h_hh_m_uncminus, h_hh_m_jerplus, h_hh_m_jerminus, h_hh_m_tesplus, h_hh_m_tesminus, h_hh_m_bplus, h_hh_m_bminus, h_hh_m_trigger;
    
  vector< TH1D * > h_hw_e, h_hw_e_plus, h_hw_e_minus, h_hw_e_uncplus, h_hw_e_uncminus, h_hw_e_jerplus, h_hw_e_jerminus, h_hw_e_tesplus, h_hw_e_tesminus, h_hw_e_bplus, h_hw_e_bminus, h_hw_e_trigger; 
  vector< TH1D * > h_hw_m, h_hw_m_plus, h_hw_m_minus, h_hw_m_uncplus, h_hw_m_uncminus, h_hw_m_jerplus, h_hw_m_jerminus, h_hw_m_tesplus, h_hw_m_tesminus, h_hw_m_bplus, h_hw_m_bminus, h_hw_m_trigger;

  if(processedFile){ 
    h          = (TH2D*) processedFile->Get(hName);  
    h_plus     = (TH2D*) processedFile->Get(hName_plus); 
    h_minus    = (TH2D*) processedFile->Get(hName_minus); 
    h_uncplus  = (TH2D*) processedFile->Get(hName_uncplus);
    h_uncminus = (TH2D*) processedFile->Get(hName_uncminus);
    h_jerplus  = (TH2D*) processedFile->Get(hName_jerplus);
    h_jerminus = (TH2D*) processedFile->Get(hName_jerminus);
    h_tesplus  = (TH2D*) processedFile->Get(hName_tesplus);
    h_tesminus = (TH2D*) processedFile->Get(hName_tesminus);
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
      h_hh_e_tesplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_tesplus) ); h_hh_m_tesplus.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_tesplus) );   
      h_hh_e_tesminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_tesminus)); h_hh_m_tesminus.push_back((TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_tesminus));  
      h_hh_e_bplus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_bplus) );   h_hh_m_bplus.push_back(   (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_bplus)   );   
      h_hh_e_bminus.push_back(  (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_bminus));   h_hh_m_bminus.push_back(  (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_bminus)  );  
      if (triggerunc ){ 
        h_hh_e_trigger.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_e_tau_trigger));  h_hh_m_trigger.push_back( (TH1D *)(processedHHFiles[i])->Get(hName_higgs_m_tau_trigger) );  
      }else{ h_hh_e_trigger.push_back(0);  h_hh_m_trigger.push_back(0);  }
    }
    else { 
      h_hh_e.push_back(0);  h_hh_e_plus.push_back(0); h_hh_e_minus.push_back(0); h_hh_e_trigger.push_back(0); 
      h_hh_e_uncplus.push_back(0); h_hh_e_uncminus.push_back(0); h_hh_e_jerplus.push_back(0); h_hh_e_jerminus.push_back(0); h_hh_e_tesplus.push_back(0); h_hh_e_tesminus.push_back(0);
      h_hh_e_bplus.push_back(0);   h_hh_e_bminus.push_back(0);

      h_hh_m.push_back(0);  h_hh_m_plus.push_back(0); h_hh_m_minus.push_back(0); h_hh_m_trigger.push_back(0);
      h_hh_m_uncplus.push_back(0); h_hh_m_uncminus.push_back(0); h_hh_m_jerplus.push_back(0); h_hh_m_jerminus.push_back(0); h_hh_m_tesplus.push_back(0); h_hh_m_tesminus.push_back(0);
      h_hh_m_bplus.push_back(0); h_hh_m_bminus.push_back(0); 

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
      h_hw_e_tesplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_tesplus) ); h_hw_m_tesplus.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_tesplus) ); 
      h_hw_e_tesminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_tesminus)); h_hw_m_tesminus.push_back((TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_tesminus)); 
      h_hw_e_bplus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_bplus) );   h_hw_m_bplus.push_back(   (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_bplus)   ); 
      h_hw_e_bminus.push_back(  (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_bminus));   h_hw_m_bminus.push_back(  (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_bminus)  ); 
      if(triggerunc ){
        h_hw_e_trigger.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_e_tau_trigger));  h_hw_m_trigger.push_back( (TH1D *)(processedHWFiles[i])->Get(hName_higgs_m_tau_trigger));
      }else{ h_hw_e_trigger.push_back(0);  h_hw_m_trigger.push_back(0);}
    }
    else { 
      h_hw_e.push_back(0); h_hw_e_plus.push_back(0); h_hw_e_minus.push_back(0); h_hw_e_trigger.push_back(0);
      h_hw_e_uncplus.push_back(0); h_hw_e_uncminus.push_back(0); h_hw_e_jerplus.push_back(0); h_hw_e_jerminus.push_back(0); h_hw_e_tesplus.push_back(0); h_hw_e_tesminus.push_back(0);
      h_hw_e_bplus.push_back(0); h_hw_e_bminus.push_back(0);

      h_hw_m.push_back(0); h_hw_m_plus.push_back(0); h_hw_m_minus.push_back(0); h_hw_m_trigger.push_back(0);
      h_hw_m_uncplus.push_back(0); h_hw_m_uncminus.push_back(0); h_hw_m_jerplus.push_back(0); h_hw_m_jerminus.push_back(0); h_hw_m_tesplus.push_back(0); h_hw_m_tesminus.push_back(0);
      h_hw_m_bplus.push_back(0); h_hw_m_bminus.push_back(0);

    }
  }

  double dataTotal[14], dataTotal_e[14], dataTotal_m[14], syst_error_plus[14], syst_error_minus[14];
   
//   // copy all data to matrix ( x = 14 steps of selection ,y =16 = 13 ev types + 1 total + 2 for S/B and S/sqrt B )
//   double d[14][16], e[14][16], d_plus[14][16], d_minus[14][16], 
//          d_uncplus[14][16],    d_uncminus[14][16], 
//          d_jerplus[14][16],    d_jerminus[14][16], 
//          d_bplus[14][16],      d_bminus[14][16], 
//          d_trigger[14][16],
//          systErrPlus[14][16],  systErrMinus[14][16];

  // copy all data to matrix ( x = 14 steps of selection ,y =18 = 15 ev types + 1 total + 2 for S/B and S/sqrt B )
  double d[14][18], e[14][18], d_plus[14][18], d_minus[14][18], 
         d_uncplus[14][18],    d_uncminus[14][18], 
         d_jerplus[14][18],    d_jerminus[14][18], 
         d_bplus[14][18],      d_bminus[14][18], 
         d_tesplus[14][18],      d_tesminus[14][18], 
         d_trigger[14][18],
         systErrPlus[14][18],  systErrMinus[14][18];
         
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
      double data_hh_e_tesplus, data_hh_e_tesminus, data_hh_m_tesplus, data_hh_m_tesminus;
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

      if( h_hh_e_tesplus[samp]  !=0 ){ data_hh_e_tesplus  = h_hh_e_tesplus[samp]->GetBinContent(x+1);   } else { data_hh_e_tesplus = 0 ;       }
      if( h_hh_e_tesminus[samp] !=0 ){ data_hh_e_tesminus = h_hh_e_tesminus[samp]->GetBinContent(x+1) ; } else { data_hh_e_tesminus= 0 ;       } 
      if( h_hh_m_tesplus[samp]  !=0 ){ data_hh_m_tesplus  = h_hh_m_tesplus[samp]->GetBinContent(x+1);   } else { data_hh_m_tesplus = 0 ;       }
      if( h_hh_m_tesminus[samp] !=0 ){ data_hh_m_tesminus = h_hh_m_tesminus[samp]->GetBinContent(x+1);  } else { data_hh_m_tesminus= 0 ;       }
 
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
      double syst_tesplus_e  = fabs( data_hh_e[x]-data_hh_e_tesplus  );  double syst_tesplus_m  = fabs( data_hh_m[x]-data_hh_m_tesplus );
      double syst_tesminus_e = fabs( data_hh_e[x]-data_hh_e_tesminus );  double syst_tesminus_m = fabs( data_hh_m[x]-data_hh_m_tesminus );
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
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0; 
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset2){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;            
      }
      if(systset3){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0;  
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
      }
      if(systset4){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0; syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;                  
      }

      double temp_e = pow(syst_tauid_e,2) +  pow(syst_eff_e,2) + pow(syst_ttbar_cs_e,2) +pow(syst_lum_e,2);
      double temp_m = pow(syst_tauid_m,2) +  pow(syst_eff_m,2) + pow(syst_ttbar_cs_m,2) +pow(syst_lum_m,2);

      if(triggerunc){ temp_e += pow(syst_trigger_e,2 ); temp_m += pow(syst_trigger_m,2 ); }

      data_syst_e_plus.push_back(  sqrt( pow(syst_plus_e,2)   + pow(syst_uncplus_e,2)  + pow(syst_jerplus_e,2)  + pow(syst_bplus_e,2) + pow(syst_tesplus_e,2) + temp_e )); 
      data_syst_e_minus.push_back( sqrt( pow(syst_minus_e,2)  + pow(syst_uncminus_e,2) + pow(syst_jerminus_e,2) + pow(syst_bminus_e,2) + pow(syst_tesminus_e,2) + temp_e )); 
      data_syst_m_plus.push_back(  sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2) + pow(syst_tesplus_m,2) + temp_m )); 
      data_syst_m_minus.push_back( sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + pow(syst_tesminus_m,2) + temp_m ));

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
      double data_hw_e_tesplus, data_hw_e_tesminus, data_hw_m_tesplus, data_hw_m_tesminus;
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

      if( h_hw_e_tesplus[samp]  !=0 ){ data_hw_e_tesplus  = h_hw_e_tesplus[samp]->GetBinContent(x+1); } else { data_hw_e_tesplus = 0 ;      }
      if( h_hw_e_tesminus[samp] !=0 ){ data_hw_e_tesminus = h_hw_e_tesminus[samp]->GetBinContent(x+1);} else { data_hw_e_tesminus= 0 ;      } 
      if( h_hw_m_tesplus[samp]  !=0 ){ data_hw_m_tesplus  = h_hw_m_tesplus[samp]->GetBinContent(x+1); } else { data_hw_m_tesplus = 0 ;      }
      if( h_hw_m_tesminus[samp] !=0 ){ data_hw_m_tesminus = h_hw_m_tesminus[samp]->GetBinContent(x+1);} else { data_hw_m_tesminus= 0 ;      }

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

      double syst_tesplus_e  = fabs( data_hw_e[x]-data_hw_e_tesplus );  double syst_tesplus_m  = fabs( data_hw_m[x]-data_hw_m_tesplus );
      double syst_tesminus_e = fabs( data_hw_e[x]-data_hw_e_tesminus ); double syst_tesminus_m = fabs( data_hw_m[x]-data_hw_m_tesminus );

      double syst_bplus_e  = fabs( data_hw_e[x]-data_hw_e_bplus );      double syst_bplus_m  = fabs( data_hw_m[x]-data_hw_m_bplus );
      double syst_bminus_e = fabs( data_hw_e[x]-data_hw_e_bminus );     double syst_bminus_m = fabs( data_hw_m[x]-data_hw_m_bminus );

      double syst_lum_e      = LUM_ERR*data_hw_e[x];                    double syst_lum_m      = LUM_ERR*data_hw_m[x];
      double syst_eff_e      = LEP_EFF_ERR*data_hw_e[x];                double syst_eff_m      = LEP_EFF_ERR*data_hw_m[x];
      double syst_ttbar_cs_e = TTBAR_CS_ERR*data_hw_e[x];               double syst_ttbar_cs_m = TTBAR_CS_ERR*data_hw_m[x];  
      double syst_trigger_e  = data_hw_e_trigger;                       double syst_trigger_m  = data_hw_m_trigger;
      double syst_tauid_e(0);                                           double syst_tauid_m(0);

      if( x >= TAUSTEP){ syst_tauid_e = TAU_ID_ERR*data_hw_e[x]; syst_tauid_m = TAU_ID_ERR*data_hw_m[x]; }  

//       if(systset1){ 
//         syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
//         syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
//         syst_trigger_e = 0;  syst_trigger_m=0;
//       }
//       if(systset2){
//         syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
//         syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0;
//         syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
//         syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
//         syst_trigger_e = 0;  syst_trigger_m=0;
//       }
//       if(systset3){
//         syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
//         syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
//         syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
//         syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
//         syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
//       }
      if(systset1){ 
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0; 
        syst_trigger_e = 0;  syst_trigger_m=0;
      }
      if(systset2){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;            
      }
      if(systset3){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0; syst_tesplus_e = 0; syst_tesplus_m=0;  syst_tesminus_e = 0; syst_tesminus_m=0;  
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e  = 0; syst_minus_m =0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;               
      }
      if(systset4){
        syst_tauid_e =0;  syst_tauid_m=0; syst_lum_e =0; syst_lum_m=0; syst_eff_e=0; syst_eff_m=0; syst_ttbar_cs_e=0; syst_ttbar_cs_m=0; 
        syst_plus_e = 0;  syst_plus_m=0;  syst_minus_e = 0; syst_minus_m=0; syst_bplus_e = 0; syst_bplus_m=0; syst_bminus_e = 0; syst_bminus_m=0;
        syst_uncplus_e = 0;  syst_uncplus_m=0;  syst_uncminus_e = 0; syst_uncminus_m=0;
        syst_jerplus_e = 0;  syst_jerplus_m=0;  syst_jerminus_e = 0; syst_jerminus_m=0;   
        syst_trigger_e = 0;  syst_trigger_m=0;                  
      }

      double temp_e = pow(syst_tauid_e,2) +  pow(syst_eff_e,2) + pow(syst_ttbar_cs_e,2) +pow(syst_lum_e,2);
      double temp_m = pow(syst_tauid_m,2) +  pow(syst_eff_m,2) + pow(syst_ttbar_cs_m,2) +pow(syst_lum_m,2);

      if(triggerunc){ temp_e += pow(syst_trigger_e,2 ); temp_m += pow(syst_trigger_m,2 ); }


      data_syst_e_plus.push_back(  sqrt( pow(syst_plus_e,2)   + pow(syst_uncplus_e,2)  + pow(syst_jerplus_e,2)  + pow(syst_bplus_e,2)  + pow(syst_tesplus_e,2)  + temp_e )); 
      data_syst_e_minus.push_back( sqrt( pow(syst_minus_e,2)  + pow(syst_uncminus_e,2) + pow(syst_jerminus_e,2) + pow(syst_bminus_e,2) + pow(syst_tesminus_e,2)  + temp_e )); 
      data_syst_m_plus.push_back(  sqrt( pow(syst_plus_m,2)   + pow(syst_uncplus_m,2)  + pow(syst_jerplus_m,2)  + pow(syst_bplus_m,2)  + pow(syst_tesplus_m,2)  + temp_m )); 
      data_syst_m_minus.push_back( sqrt( pow(syst_minus_m,2)  + pow(syst_uncminus_m,2) + pow(syst_jerminus_m,2) + pow(syst_bminus_m,2) + pow(syst_tesminus_m,2)  + temp_m ));
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    }


    hw_e[samp] = data_hw_e; hw_m[samp] = data_hw_m; hw_e_err[samp] = data_hw_e_err; hw_m_err[samp] = data_hw_m_err;
    hw_e_syst_plus[samp] = data_syst_e_plus; hw_e_syst_minus[samp] = data_syst_e_minus; 
    hw_m_syst_plus[samp] = data_syst_m_plus; hw_m_syst_minus[samp] = data_syst_m_minus; 


  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  // reset matrix /////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int x=0;x<14;x++){ 
    //hh_e[x]= 0; hh_m[x] = 0; hw_e[x]= 0; hw_m[x] = 0;
    dataTotal[x]=0; syst_error_plus[x]=0; syst_error_minus[x]=0; 
    dataTotal_e[x]=0; dataTotal_m[x]=0;

    d_e[x]=0; d_m[x]=0; d_e_systErrPlus[x]=0; d_e_systErrMinus[x]=0; d_m_systErrPlus[x]=0; d_m_systErrMinus[x]=0;

    for(int y=0;y<17;y++){ 
      d[x][y]=0; e[x][y]=0; d_trigger[x][y]=0; d_plus[x][y]=0; d_minus[x][y]=0;
      d_uncplus[x][y]=0; d_uncminus[x][y]=0; d_jerplus[x][y]=0; d_jerminus[x][y]=0; d_tesplus[x][y]=0; d_tesminus[x][y]=0; d_bplus[x][y]=0; d_bminus[x][y]=0; 
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
    for(int y=0;y<15;y++){  
   
      // yields //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(h)          { d[x][y]          = ( h->GetBinContent(x+1,y+1)         ) ;  e[x][y] = ( h->GetBinError(x+1,y+1) );            } 
      if(h_plus)     { d_plus[x][y]     = ( h_plus->GetBinContent(x+1,y+1)    ) ;                                                    } 
      if(h_minus)    { d_minus[x][y]    = ( h_minus->GetBinContent(x+1,y+1)   ) ;                                                    }
      if(h_uncplus)  { d_uncplus[x][y]  = ( h_uncplus->GetBinContent(x+1,y+1) ) ;                                                    } 
      if(h_uncminus) { d_uncminus[x][y] = ( h_uncminus->GetBinContent(x+1,y+1)) ;                                                    }
      if(h_jerplus)  { d_jerplus[x][y]  = ( h_jerplus->GetBinContent(x+1,y+1) ) ;                                                    } 
      if(h_jerminus) { d_jerminus[x][y] = ( h_jerminus->GetBinContent(x+1,y+1)) ;                                                    }
      if(h_tesplus)  { d_tesplus[x][y]  = ( h_tesplus->GetBinContent(x+1,y+1) ) ;                                                    } 
      if(h_tesminus) { d_tesminus[x][y] = ( h_tesminus->GetBinContent(x+1,y+1)) ;                                                    }
      if(h_bplus)    { d_bplus[x][y]    = ( h_bplus->GetBinContent(x+1,y+1)   ) ;                                                    } 
      if(h_bminus)   { d_bminus[x][y]   = ( h_bminus->GetBinContent(x+1,y+1)  ) ;                                                    }
      if (triggerunc && h_trigger)  { d_trigger[x][y]  = ( h_trigger->GetBinContent(x+1,y+1) ) ;                                     }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      cout << "Name: " << hName << ", value: " << d[x][y] << ", stat error: " << e[x][y] <<endl;

//       if( XSECMEASUREMENT && ( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH)  ){  
//         double scalefactor( XSEC_EXP_ / 165. );  //note it should be 165
      if( XSECMEASUREMENT && ( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH //y==LEPLEP_TDCH
			       || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH)  ){  
        double scalefactor( XSEC_EXP_ / 234. );  //note it should be 165 // 234!!!
        d[x][y]          *= scalefactor ; 
        d_plus[x][y]     *= scalefactor ;
        d_minus[x][y]    *= scalefactor ;     
        d_uncplus[x][y]  *= scalefactor ; 
        d_uncminus[x][y] *= scalefactor ;
        d_jerplus[x][y]  *= scalefactor ;
        d_jerminus[x][y] *= scalefactor ;
	d_tesplus[x][y]  *= scalefactor ;
        d_tesminus[x][y] *= scalefactor ;
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
      double syst_tesplus  = fabs(d[x][y]-d_tesplus[x][y]);    
      double syst_tesminus = fabs(d[x][y]-d_tesminus[x][y]); 
      double syst_bplus    = fabs(d[x][y]-d_bplus[x][y]);    
      double syst_bminus   = fabs(d[x][y]-d_bminus[x][y]); 
      double syst_lum      = LUM_ERR*d[x][y];
      double syst_eff_lep  = LEP_EFF_ERR*d[x][y];
      double syst_trigger  = d_trigger[x][y];


      double sample_cs(0);

      if(      y==QCD_TDCH )                                                                                                               { sample_cs = QCD_CS_ERR ;      }
//       else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
      else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH //y==LEPLEP_TDCH
	       || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
      else if( y==SINGLETOP_TDCH )                                                                                                         { sample_cs = SINGLETOP_CS_ERR; }              
      else if( y==WJETS_TDCH     )                                                                                                         { sample_cs = WJETS_CS_ERR;     }              
      else                                                                                                                                 { sample_cs = OTHER_CS_ERR;     }

      double syst_cs = sample_cs * d[x][y];
     
      double syst_tauid(0);                                             
 
//       if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=LEPLEP_TDCH     ){ syst_tauid = TAU_ID_ERR*d[x][y];     }  
//       else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==LEPLEP_TDCH )  ){ syst_tauid = TAU_MISSID_ERR*d[x][y]; }
      if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=EE_TDCH && y!=EMU_TDCH && y!= MUMU_TDCH    ){ syst_tauid = TAU_ID_ERR*d[x][y];   }    //LEPLEP_TDCH     ){ syst_tauid = TAU_ID_ERR*d[x][y];     }  
      else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH  )  ){ syst_tauid = TAU_MISSID_ERR*d[x][y]; }    //LEPLEP_TDCH )  ){ syst_tauid = TAU_MISSID_ERR*d[x][y]; }

      // signal theory unc.
      double syst_theoryeff(0); if( XSECMEASUREMENT && ( y==ETAU_TDCH || y==MUTAU_TDCH )) syst_theoryeff = SIGNAL_EFF_ERR*d[x][y];

      if(systset1){ 
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ;
        syst_bplus = 0;  syst_bminus = 0; syst_theoryeff=0;
        syst_trigger=0;  syst_tesplus = 0;  syst_tesminus = 0;
      }
      if(systset2){
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ; 
        syst_uncplus = 0; syst_uncminus= 0;  syst_jerplus = 0;  syst_jerminus = 0; syst_plus=0; syst_minus=0; syst_theoryeff=0;
        syst_trigger=0; syst_tesplus = 0;  syst_tesminus = 0;
      }
      if(systset3){
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ;
        syst_bplus = 0;  syst_bminus = 0; syst_theoryeff=0; syst_tesplus = 0;  syst_tesminus = 0;
        syst_uncplus = 0; syst_uncminus= 0;  syst_jerplus = 0;  syst_jerminus = 0; syst_plus=0; syst_minus=0; syst_theoryeff=0;
        
      }
      if(systset4){
        syst_tauid = 0;  syst_lum = 0; syst_eff_lep=0; syst_cs=0; ; 
        syst_uncplus = 0; syst_uncminus= 0;  syst_jerplus = 0;  syst_jerminus = 0; syst_plus=0; syst_minus=0; syst_theoryeff=0;
        syst_trigger=0; syst_bplus = 0;  syst_bminus = 0;
      }
     
      // GOOD ERRORS
      double temp;
      temp = pow(syst_lum,2)+ pow(syst_tauid,2) + pow(syst_cs,2) + pow(syst_eff_lep,2)+pow(syst_theoryeff,2);
      if(triggerunc){ temp += pow(syst_trigger,2); }
      systErrPlus[x][y]  = sqrt( pow(syst_plus,2)  + pow(syst_uncplus,2)  + pow(syst_jerplus,2)  + pow(syst_bplus,2)  + pow(syst_tesplus,2)  + temp ); 
      systErrMinus[x][y] = sqrt( pow(syst_minus,2) + pow(syst_uncminus,2) + pow(syst_jerminus,2) + pow(syst_bminus,2) + pow(syst_tesminus,2)  + temp );
      
   
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
    double signal_tesplus(d_tesplus[x][0] + d_tesplus[x][1]), signal_tesminus( d_tesminus[x][0] + d_tesminus[x][1]);
    double signal_bplus(  d_bplus[x][0]   + d_bplus[x][1]  ), signal_bminus(   d_bminus[x][0]   + d_bminus[x][1]);
    double signal_trigger = pow( d_trigger[x][0],2) + pow(d_trigger[x][1],2) ;
    double err_signal2 = (e[x][0])*(e[x][0]) + (e[x][1])*(e[x][1]) ;

    
    // standard model background 
    double bkg(0),bkg_plus(0), bkg_minus(0), bkg_uncplus(0), bkg_uncminus(0), bkg_jerplus(0), bkg_jerminus(0), bkg_tesplus(0), bkg_tesminus(0), bkg_bplus(0), bkg_bminus(0), err_bkg2(0), bkg_trigger(0); 
    for(int y2=2;y2<15;y2++){ 
      bkg += d[x][y2]; err_bkg2 += (e[x][y2])*(e[x][y2]);
      bkg_plus     += d_plus[x][y2];     bkg_minus   += d_minus[x][y2];  
      bkg_uncminus += d_uncminus[x][y2]; bkg_uncplus += d_uncplus[x][y2];       
      bkg_jerminus += d_jerminus[x][y2]; bkg_jerplus += d_jerplus[x][y2];
      bkg_tesminus += d_tesminus[x][y2]; bkg_tesplus += d_tesplus[x][y2];
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

    d_plus[x][15]    = signal_plus  + bkg_plus;    d_uncplus[x][15]  = signal_uncplus  + bkg_uncplus;   d_jerplus[x][15]  = signal_jerplus  + bkg_jerplus; 
    d_bplus[x][15]   = signal_bplus + bkg_bplus; d_tesplus[x][15]  = signal_tesplus  + bkg_tesplus; 
    
    d_minus[x][15]   = signal_minus + bkg_minus;   d_uncminus[x][15] = signal_uncminus + bkg_uncminus;  d_jerminus[x][15] = signal_jerminus + bkg_jerminus;
    d_bminus[x][15]  = signal_bminus  + bkg_bminus; d_tesminus[x][15] = signal_tesminus + bkg_tesminus;

    d_trigger[x][15] = sqrt( signal_trigger +  bkg_trigger );


    d[x][15]=signal+bkg; e[x][15] = sqrt( err_signal2 + err_bkg2 );
    // S/B
    d[x][16] = sob;      e[x][16] = err_sob;
    // S/Sqrt(B)
    d[x][17] = sosqrtb;  e[x][17] = err_sosqrtb;


    double temp_syst_plus      = fabs(d_plus[x][15]    - d[x][15]);   
    double temp_syst_minus     = fabs(d_minus[x][15]   - d[x][15]);
    double temp_syst_uncplus   = fabs(d_uncplus[x][15] - d[x][15]);   
    double temp_syst_uncminus  = fabs(d_uncminus[x][15]- d[x][15]);
    double temp_syst_jerplus   = fabs(d_jerplus[x][15] - d[x][15]);   
    double temp_syst_jerminus  = fabs(d_jerminus[x][15]- d[x][15]);
    double temp_syst_tesplus   = fabs(d_tesplus[x][15] - d[x][15]);   
    double temp_syst_tesminus  = fabs(d_tesminus[x][15]- d[x][15]);
    double temp_syst_bplus     = fabs(d_bplus[x][15]   - d[x][15]);
    double temp_syst_bminus    = fabs(d_bminus[x][15]  - d[x][15]);

    //    double temp_syst_ttbar      = TTBAR_CS_ERR     * ( d[x][ETAU_TDCH]+d[x][MUTAU_TDCH]+d[x][LEPJETS_TDCH]+d[x][LEPLEP_TDCH]+d[x][TAUJETS_TDCH]+d[x][TAUTAU_TDCH]+d[x][ALLJETS_TDCH]);
    double temp_syst_ttbar      = TTBAR_CS_ERR     * ( d[x][ETAU_TDCH]+d[x][MUTAU_TDCH]+d[x][LEPJETS_TDCH]+d[x][EE_TDCH]+d[x][EMU_TDCH]+d[x][MUMU_TDCH]
						       +d[x][TAUJETS_TDCH]+d[x][TAUTAU_TDCH]+d[x][ALLJETS_TDCH]);
    double temp_syst_wjets      = WJETS_CS_ERR     * ( d[x][WJETS_TDCH] );
    double temp_syst_singletop  = SINGLETOP_CS_ERR * ( d[x][SINGLETOP_TDCH]);
    double temp_syst_diboson    = OTHER_CS_ERR     * ( d[x][EWKDI_TDCH] );
    double temp_syst_zjets      = OTHER_CS_ERR     * ( d[x][ZJETS_TDCH] + d[x][ZTAUTAU_TDCH]);
    double temp_syst_qcd        = QCD_CS_ERR       * ( d[x][QCD_TDCH] );

    double temp_syst_tauid      = 0; if( x>=TAUSTEP ){ temp_syst_tauid     = TAU_ID_ERR     *(d[x][ETAU_TDCH]+d[x][MUTAU_TDCH]+d[x][TAUJETS_TDCH]+d[x][TAUTAU_TDCH]);}
    //    double temp_syst_misstauid  = 0; if( x>=TAUSTEP ){ temp_syst_misstauid = TAU_MISSID_ERR *(d[x][LEPJETS_TDCH]+d[x][LEPLEP_TDCH]+d[x][ALLJETS_TDCH]);              }
    double temp_syst_misstauid  = 0; if( x>=TAUSTEP ){ temp_syst_misstauid = TAU_MISSID_ERR *(d[x][LEPJETS_TDCH]+d[x][EE_TDCH]+d[x][EMU_TDCH]+d[x][MUMU_TDCH]+d[x][ALLJETS_TDCH]);              }

    double temp_syst_lum      = LUM_ERR*d[x][15];
    double temp_syst_eff_lep  = LEP_EFF_ERR*d[x][15];

    double syst_theoryeff = 0;
    double temp_syst_trigger  = d_trigger[x][15] ;

    if (XSECMEASUREMENT) syst_theoryeff = SIGNAL_EFF_ERR*signal;
     
    if(systset1){ 
      temp_syst_misstauid = 0;
      temp_syst_eff_lep   = 0;
      temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
      temp_syst_bplus =0; temp_syst_bminus = 0;  temp_syst_tesplus = 0; temp_syst_tesminus = 0;
      syst_theoryeff=0;
      temp_syst_trigger=0;
    }
    if(systset2){
      temp_syst_misstauid = 0;
      temp_syst_eff_lep   = 0;
      temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
      temp_syst_plus = 0; temp_syst_minus = 0; temp_syst_tesplus = 0; temp_syst_tesminus = 0;
      temp_syst_uncplus = 0; temp_syst_uncminus = 0;
      temp_syst_jerplus = 0; temp_syst_jerminus = 0;
      syst_theoryeff=0;
      temp_syst_trigger=0;
    }
   if(systset3){    
     temp_syst_misstauid = 0;
     temp_syst_eff_lep   = 0;
     temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
     temp_syst_bplus =0; temp_syst_bminus = 0; temp_syst_tesplus = 0; temp_syst_tesminus = 0;
     temp_syst_plus = 0; temp_syst_minus = 0;
     temp_syst_uncplus = 0; temp_syst_uncminus = 0;
     temp_syst_jerplus = 0; temp_syst_jerminus = 0;
     syst_theoryeff=0;
   }
    if(systset4){
      temp_syst_misstauid = 0;
      temp_syst_eff_lep   = 0;
      temp_syst_tauid =0; temp_syst_qcd=0; temp_syst_ttbar =0; temp_syst_wjets=0; temp_syst_singletop=0; temp_syst_diboson=0; temp_syst_zjets=0; temp_syst_lum=0;
      temp_syst_plus = 0; temp_syst_minus = 0; temp_syst_bplus = 0; temp_syst_bminus = 0;
      temp_syst_uncplus = 0; temp_syst_uncminus = 0;
      temp_syst_jerplus = 0; temp_syst_jerminus = 0;
      syst_theoryeff=0;
      temp_syst_trigger=0;
    }


   // Good errors
   double temp(0);
   temp = pow(temp_syst_lum,2) + pow(temp_syst_tauid,2)+pow(temp_syst_misstauid,2)+pow(temp_syst_eff_lep,2)+ 
   pow(temp_syst_wjets,2)+pow(temp_syst_qcd,2)+pow(temp_syst_ttbar,2)+pow(temp_syst_singletop,2)+pow(temp_syst_diboson,2)+pow(temp_syst_zjets,2) + pow(syst_theoryeff,2);

   //cout<<"\n debug : temp should be zero : "<<temp<<endl;

   if( triggerunc ){ temp += pow(temp_syst_trigger,2);  }

   syst_error_plus[x]  = sqrt( pow(temp_syst_plus,2)  + pow(temp_syst_uncplus,2) + pow(temp_syst_jerplus,2) + pow(temp_syst_tesplus,2) + pow(temp_syst_bplus,2) + temp );
   syst_error_minus[x] = sqrt( pow(temp_syst_minus,2) + pow(temp_syst_uncminus,2) + pow(temp_syst_jerminus,2) + pow(temp_syst_tesminus,2) + pow(temp_syst_bminus,2) + temp );

   cout<<"\n debug : trigger unc = "<<temp_syst_trigger<<" syst plus = "<<syst_error_plus[x]<<endl;
  

  }
  /////////////////////////////////////////////////////////////////////////////////////////////

  TString fName = name+key+TString(".tex");
  TString vName("-v2"), systName("");


  if(systset1)             systName=TString("-JESMETJER");
  if(systset2)             systName=TString("-BTAG");
  if(systset3)             systName=TString("-TRIGGER");
  if(systset4)             systName=TString("-TAUSCALE");
  if( (systset1 && systset2) || (systset1 && systset3) || (systset2 && systset3)  ) systName=TString("-ERROR"); // not valid condition
  


  if(!hh_contrib && !hw_contrib){ fName = name+key+TString("-SM")+systName+vName+TString(".tex");}
  if(hh_contrib && hw_contrib ) { fName = name+key+TString("-HIGGS")+systName+vName+TString(".tex");}

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
    //    "lepton+jets", "dilepton", 
    "lepton+jets", "ee", "e\\mu", "\\mu\\mu",
    "tau+jets", "$\\tau\\tau$", "full had.", 
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

  fprintf(f,"\\hline \n");
  //fprintf(f,"\\hline \n");
 
  // size of extra samples
  int numbHWSamples = processedHWFiles.size();
  int numbHHSamples = processedHHFiles.size();

  int extraLines_hw = 2*numbHWSamples;
  int extraLines_hh = 2*numbHHSamples;
  int extraLines    = extraLines_hw + extraLines_hh;
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
   
    if(includeSM){
    // PRINT extra background sources /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( l ==1 ){
    
      if(hh_contrib || hw_contrib) fprintf(f,"\\hline \n\n"); 

      for(int k=0; k<15;k++){ // 15 sources of background

        //if(k== 3 || k ==2 || k==5 || k == 6 || k==7 || k==11 ) continue;

        TString d2;
        if     (k ==0 ) { d2=TString("tau dilepton from $t\\bar{t}$ & $e\\tau_h$ & ");}
        else if(k == 1) { d2=TString(" & $\\mu\\tau_h$ & ");}
        else if(k == 2) { fprintf(f,"\\hline ");  d2=TString(" & lepton+jets & ");}
//         else if(k == 3) { d2=TString(" & dilepton & ");}
//         else if(k == 4) { d2=TString(" other $t\\bar{t}$ & tau+jets & ");}
//         else if(k == 5) { d2=TString(" & $\\tau\\tau$ & ");}
//         else if(k == 6) { d2=TString(" & full had. & ");}
//         else if(k == 7) { fprintf(f,"\\hline");  d2=TString(" & W+jets & ");}
//         else if(k == 8) { d2=TString(" & Z(ll)+jets & ");}
//         else if(k == 9) { d2=TString("  other bkg & Z($\\tau\\tau$)+jets & ");}
//         else if(k == 10){ d2=TString(" & single top & ");}
//         else if(k == 11){ d2=TString(" & QCD & "); }
//         else if(k == 12){ d2=TString(" & WW,WZ,ZZ & ");}
        else if(k == 3) { d2=TString(" & ee & ");}
        else if(k == 4) { d2=TString(" & e$\\mu$ & ");}
        else if(k == 5) { d2=TString(" & $\\mu\\mu$ & ");}
        else if(k == 6) { d2=TString(" other $t\\bar{t}$ & tau+jets & ");}
        else if(k == 7) { d2=TString(" & $\\tau\\tau$ & ");}
        else if(k == 8) { d2=TString(" & full had. & ");}
        else if(k == 9) { fprintf(f,"\\hline");  d2=TString(" & W+jets & ");}
        else if(k == 10) { d2=TString(" & Z(ll)+jets & ");}
        else if(k == 11) { d2=TString("  other bkg & Z($\\tau\\tau$)+jets & ");}
        else if(k == 12){ d2=TString(" & single top & ");}
        else if(k == 13){ d2=TString(" & QCD & "); }
        else if(k == 14){ d2=TString(" & WW,WZ,ZZ & ");}

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
        fprintf(f, data4, d[COL2][15], d[COL3][15], d[COL4][15], d[COL5][15], d[COL6][15], d[COL7][15], d[COL8][15], d[COL9][15]);
	if(printAllErrors) {
        d5=TString(" & stat $\\pm$ syst & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ &  $\\pm$ %8.1f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ \\\\ \n"); 
        line.Append(d5);  data = line.Data();
	  fprintf(f, data, e[COL2][15], syst_error_plus[COL2], syst_error_minus[COL2],
		  e[COL3][15], syst_error_plus[COL3], syst_error_minus[COL3],    
		  e[COL4][15], syst_error_plus[COL4], syst_error_minus[COL4], 
		  e[COL5][15], syst_error_plus[COL5], syst_error_minus[COL5],
		  e[COL6][15], syst_error_plus[COL6], syst_error_minus[COL6],
		  e[COL7][15], syst_error_plus[COL7], syst_error_minus[COL7],
		  e[COL8][15], syst_error_plus[COL8], syst_error_minus[COL8],
		  e[COL9][15], syst_error_plus[COL9], syst_error_minus[COL9]);
	}
      } 


      if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
        fprintf(f, data_e, d_e[COL2], d_e[COL3], d_e[COL4], d_e[COL5], d_e[COL6], d_e[COL7], d_e[COL8]);
        fprintf(f, data_m, d_m[COL2], d_m[COL3], d_m[COL4], d_m[COL5], d_m[COL6], d_m[COL7], d_m[COL8]);
        d4= TString(" MC(SM)  &  total & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f\\\\ \n");
        line4.Append(d4);  const char * data4 = line4.Data();
        fprintf(f, data4, d[COL2][15], d[COL3][15], d[COL4][15], d[COL5][15], d[COL6][15], d[COL7][15], d[COL8][15]);
	if(printAllErrors) {
        d5=TString(" & stat $\\pm$ syst & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ &  $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ & $\\pm$ %8.0f $^{+ %8.0f }_{- %8.0f }$ \\\\ \n"); 
        line.Append(d5);  data = line.Data();
        fprintf(f, data, e[COL2][15], syst_error_plus[COL2], syst_error_minus[COL2],
          e[COL3][15], syst_error_plus[COL3], syst_error_minus[COL3],    
          e[COL4][15], syst_error_plus[COL4], syst_error_minus[COL4], 
          e[COL5][15], syst_error_plus[COL5], syst_error_minus[COL5],
          e[COL6][15], syst_error_plus[COL6], syst_error_minus[COL6],
          e[COL7][15], syst_error_plus[COL7], syst_error_minus[COL7],
          e[COL8][15], syst_error_plus[COL8], syst_error_minus[COL8]
        );
	}
      } 

    }
    if(l==1 && includeSoverB ){
      fprintf(f,"\n \\hline"); fprintf(f,"\n \\hline"); 
      TString d1(""),d2(""),line1(""),line2(""); const char * data1; const char * data2;
      
      if( STARTINGPOINT == STARTING_AT_LJETS_ ){
        d1=TString(" & S/B &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n");
        line1.Append(d1); data1 = line1.Data();
        fprintf(f, data1, d[COL2][16], d[COL3][16], d[COL4][16], d[COL5][16], d[COL6][16], d[COL7][16] , d[COL8][16], d[COL9][16]);

        d2=TString(" & S/$\\sqrt{B}$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");
        line2.Append(d2); data2 = line2.Data();
        fprintf(f, data2, d[COL2][17], d[COL3][17], d[COL4][17], d[COL5][17],  d[COL6][17],  d[COL7][17],  d[COL8][17],  d[COL9][17]);
      }
      if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){
        d1=TString(" & S/B &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n");
        line1.Append(d1); data1 = line1.Data();
        fprintf(f, data1, d[COL2][16], d[COL3][16], d[COL4][16], d[COL5][16], d[COL6][16], d[COL7][16], d[COL8][16]);

        d2=TString(" & S/$\\sqrt{B}$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n");
        line2.Append(d2); data2 = line2.Data();
        fprintf(f, data2, d[COL2][17], d[COL3][17], d[COL4][17], d[COL5][17],  d[COL6][17], d[COL7][17],  d[COL8][17]);
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
void TauDileptonTableBuilder::summaryTable( bool detailed, bool higgs, bool systset1, bool systset2, bool systset3){ 


  //  double ttbar_init(20416081); //SIMPLE evt processed
  double ttbar_init(8228517); // evt processed from debug.txt
  // double ttbar_init(4939452); // evt in minitrees
  
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
  double lum(18072.2); double lumFactor(1); double scalefactor(1.); double scale(0.51393073); //double scale(0.02184544); // NEW VALUES

  // HIGGS ===>
  //double tau_fakes(54.9);   double tau_fakes_stat2(6.7*6.7); double tau_fakes_syst2(8.3*8.3); 
  // XSEC  ===> 
  //double tau_fakes(54.8);     double tau_fakes_stat2(6.5*6.5); double tau_fakes_syst2(8.2*8.2);
  //double tau_fakes(54.3);     double tau_fakes_stat2(6.4*6.4); double tau_fakes_syst2(8.1*8.1);    
  //  double tau_fakes(54.3);     double tau_fakes_stat2(0); double tau_fakes_syst2(8.1*8.1+6.4*6.4);    //no mt, new single top xsec (all syst included)
  //  double tau_fakes(111.9);     double tau_fakes_stat2(0); double tau_fakes_syst2(8.1*8.1+6.4*6.4);    //no mt, new single top xsec (all syst included)
  //  double tau_fakes(165.4*0.69);     double tau_fakes_stat2( ( (0.02*0.02)/(0.69*0.69) + (23.7*23.7)/(165.4*165.4) ) * 165.4*165.4*0.69*0.69); double tau_fakes_syst2(8.1*8.1+6.4*6.4);    // faux from prelim likelihood fit
  //  double tau_fakes(108.8);     double tau_fakes_stat2( 7.3 *7.3); double tau_fakes_syst2(8.1*8.1+6.4*6.4);    // faux from prelim likelihood fit
  //  sqrt( (eA.B)^2 + (A.eB)^2 )
  //  double tau_fakes(206.92*0.69);     double tau_fakes_stat2( 22.5*22.5*0.69*0.69 + 0.02*0.02*206.92*206.92 ); double tau_fakes_syst2(8.1*8.1+6.4*6.4);    // faux from prelim likelihood fit in 4625/pb
  //fn  double tau_fakes(1679.77);     double tau_fakes_stat2( 0); double tau_fakes_syst2(144);    // MT > 40 GeV applied
  double tau_fakes(1931);     double tau_fakes_stat2( 0); double tau_fakes_syst2(170*170);    // from likelihood fit in 4903/pb
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

  int TAU_NOMT, OS_NOMT;

  if (APPLY_MT_CUT_) { TAU_NOMT = TAU_STEP2; OS_NOMT = OS_STEP2; }
  else { TAU_NOMT = TAU_STEP2 - 1; OS_NOMT = OS_STEP2 - 1; }
  
  TAUSTEP = TAU_NOMT; BTAG1STEP = BTAG1_STEP2; /*STEP=BTAG1_STEP2; */ STEP = OS_NOMT;
// nomt //  if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;    COL4=BTAG1_STEP2; COL5=TAU_NOMT; COL6=OS_NOMT;    COL7=BTAG2_STEP2; COL8=R_STEP1; }
// nomt //  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;   COL3=BTAG1_STEP2 ;  COL4=TAU_NOMT;   COL5=OS_NOMT;  COL6=BTAG2_STEP2; COL7=R_STEP2; }                 

  if(     STARTINGPOINT == STARTING_AT_LJETS_         ){ COL2 = JET3_STEP2; COL3=MET_STEP2;   COL4=MT_STEP2;  COL5=BTAG1_STEP2; COL6=TAU_NOMT; COL7=OS_NOMT;    COL8=BTAG2_STEP2; COL9=R_STEP1; }
  else if( STARTINGPOINT == STARTING_AT_LJETSPLUSMET_ ){ COL2 = MET_STEP2;  COL3=MT_STEP2; COL4=BTAG1_STEP2 ;  COL5=TAU_NOMT;   COL6=OS_NOMT;  COL7=BTAG2_STEP2; COL8=R_STEP2; }                 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // table options /////////////////////////////////////////////////////////////////
  bool incDocument(true);
  bool showAllHiggs(true), processHH(true),processHW(true);
  if(!higgs){ showAllHiggs=false; processHH=false; processHW=false;}
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
  
  TFile * processedMCFile = TFile::Open(outputArea_+TString("nomt-NOBT-V1-mc-MU-20GeV/out-mc.root"));
  TFile * processedDataFile = TFile::Open(outputArea_+TString("nomt-NEWJ-V1-data-MU-20GeV/out-data.root")); 

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
  vector<TString> hhFiles; vector<TString> hwFiles;
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
 
  vector<TFile *> processedHWFiles; vector<TFile *> processedHHFiles;
  for(uint i = 0; i<hhFiles.size(); i++ ){  processedHHFiles.push_back( TFile::Open( hhFiles[i] ) ); }
  for(uint i = 0; i<hwFiles.size(); i++ ){  processedHWFiles.push_back( TFile::Open( hwFiles[i] ) ); }
  cout << "Open hh and wh files if any" << endl;
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
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





  // HH / WH data/error holders
  map< pair< TString, int > , vector<double> > hh, hw, hh_err, hw_err;
  map< pair< TString, int > , vector<double> > hh_syst_plus, hh_syst_minus, hw_syst_plus, hw_syst_minus;

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
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << "Processed hw and hh histos, if any" << endl;



    // fill data points////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int x=0;x<12;x++){
    
      // fill data /////////////////////////////////////////////////////////////////////////////////
      pair<TString,int> dk(algo,x); if( hData[algo] ) dataTotal[dk]= hData[algo]->GetBinContent(x+1);
      //////////////////////////////////////////////////////////////////////////////////////////////

      // fill values with histogram data /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int y=0;y<17;y++){  

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

	//        if( XSECMEASUREMENT && y==LEPLEP_TDCH ){   
        if( XSECMEASUREMENT && ( y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH ) ){   
          // Note : we re-scale ttbar lepton+lepton ( we get it from MC)
          double scalefactor( XSEC_EXP_/234.); //note it should be 234 @ 8 TeV //note it should be 165 @ 7 TeV
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
	//        else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
        else if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ sample_cs = TTBAR_CS_ERR;     }
        else if( y == SINGLETOP_TDCH )                                                                                                       { sample_cs = SINGLETOP_CS_ERR; }              
        else if( y == WJETS_TDCH     )                                                                                                       { sample_cs = WJETS_CS_ERR;     }              
        else                                                                                                                                 { sample_cs = OTHER_CS_ERR;     }
        

	//        if(XSECMEASUREMENT && y==LEPLEP_TDCH ){ sample_cs = XSEC_EXP_ERR_; }
        if(XSECMEASUREMENT && ( y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH ) ){ sample_cs = XSEC_EXP_ERR_; }

        syst_cs = sample_cs * dMC[k];
 
        double syst_theory(0);
        if( XSECMEASUREMENT && (y==ETAU_TDCH || y==MUTAU_TDCH) ){ syst_theory=SIGNAL_EFF_ERR*dMC[k];}


//         if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=LEPLEP_TDCH   ){ syst_tauid = TAU_ID_ERR*dMC[k];     }  
//         else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==LEPLEP_TDCH )){ syst_tauid = TAU_MISSID_ERR*dMC[k]; }
	if     ( x >= TAUSTEP && y!=LEPJETS_TDCH && y!=ZJETS_TDCH  && y!=EE_TDCH && y!=EMU_TDCH && y!= MUMU_TDCH    ){ syst_tauid = TAU_ID_ERR*dMC[k];     }
	else if( x >= TAUSTEP && ( y==LEPJETS_TDCH || y==ZJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH  )  ){ syst_tauid = TAU_MISSID_ERR*dMC[k]; }

 
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
    //    vector<int> codes; codes.push_back(LEPLEP_TDCH); codes.push_back(ZJETS_TDCH); codes.push_back(ZTAUTAU_TDCH); codes.push_back(SINGLETOP_TDCH); codes.push_back(EWKDI_TDCH); 
    vector<int> codes; codes.push_back(EE_TDCH); codes.push_back(EMU_TDCH); codes.push_back(MUMU_TDCH); //codes.push_back(LEPLEP_TDCH);
    codes.push_back(ZJETS_TDCH); codes.push_back(ZTAUTAU_TDCH); codes.push_back(SINGLETOP_TDCH); codes.push_back(EWKDI_TDCH); 
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

     
    for(int y=0;y<17;y++){

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

//       if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==LEPLEP_TDCH || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ 
      if( y==ETAU_TDCH || y==MUTAU_TDCH || y==LEPJETS_TDCH || y==EE_TDCH || y==EMU_TDCH || y==MUMU_TDCH
	  || y==TAUJETS_TDCH || y==TAUTAU_TDCH || y==ALLJETS_TDCH ){ 
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

  processedMCFile->Close();
  processedDataFile->Close();

 

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


} // namespace tableutils 
