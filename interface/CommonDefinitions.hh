#ifndef _CommonDefinitions_hh
#define _CommonDefinitions_hh
// FIXME: honestly, most of this stuff should be moved to cfg files (numerical parameters)
// Sample enums and such can (and should) remain here
// xsecs and their errors should remain here. Or in SampleProcessor.

#include <TString.h>

namespace commondefinitions{

  enum tanc_fr_Disc     { TANCFR_ONE   = 24, TANCFR_HALF   = 25, TANCFR_QUARTER = 26, TANCFR_TENTH = 27 };
  
  enum tau_hpsIso_Disc  { HPS_LOOSEISO = 25, HPS_MEDIUMISO = 26, HPS_TIGHTISO   = 27 };
  
  enum TauType          { PFTAU, CALOTAU, HPSTAU, PFLOWTAU, HPSTANCTAU, UNKNOWNTAU };

  
  //extern double LUM_ERR            = 0.045; // check it 
  extern double LUM_ERR            ;
  
  extern double TAU_ID_ERR         ;  
  
  extern double QCD_CS_ERR         ;
  
  extern double TTBAR_CS_ERR       ;   
  
  extern double SINGLETOP_CS_ERR   ;
  extern double OTHER_CS_ERR       ;
  
  extern double LEP_EFF_ERR        ;        // 3 electrons, 2 muons
  extern double TAU_MISSID_ERR     ;        // BULSHIT //0.15;
  
  extern double WJETS_CS_ERR       ;
  
  //  extern double BTAG_EFF_ERR       = 0.10;
  extern double BTAG_EFF_ERR       ;
  
  extern double SIGNAL_EFF_ERR     ;  
  
  extern double BTAG_eff_R_ ;
  extern double BTAG_eff_F_ ;
  
  
  extern double XSEC_EXP_   ;//164.4;       //  was 158, 164.4  2.8 (stat.)  11.9 (syst.)  7.4
  extern double XSEC_EXP_ERR_;          // 14.29/164.4 => exp error, theoretical -> 165 +4/-9 +/-7
  
  // SFb_error taken from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt (2011 prescription)
  
  ///*
  //  float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  //  float ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  //
  //   Tagger: TCHEL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  //   SFb = 0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)));
  //*/
  //   extern double SFb_error[] = {
  //     0.0244956,
  //     0.0237293,
  //     0.0180131,
  //     0.0182411,
  //     0.0184592,
  //     0.0106444,
  //     0.011073,
  //     0.0106296,
  //     0.0175259,
  //     0.0161566,
  //     0.0158973,
  //     0.0186782,
  //     0.0371113,
  //     0.0289788 
  //   };
  
  /*
    Tagger: CSVL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
    SFb = 1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x)));
  */
  extern double SFb_error[];  
  
  
  /*  
      extern double BTAG_eff_R_ = 0.812272;
      extern double BTAG_eff_F_ = 0.17211;
  */
  
  
  /*
  // 30 GeV
  Real Bs : numerator is : 8320.26 denominator is : 10257.8 Computed efficiency is : 0.811112 +/- 0.00278575
  Fake 2 Bs : numerator is : 7633.33 denominator is : 44280 Computed efficiency is : 0.172388 +/- 0.00187882
  
  // 20 GeV
  Real Bs : numerator is : 9635.2 denominator is : 11898.7 Computed efficiency is : 0.809769 +/- 0.00280644
  Fake 2 Bs : numerator is : 10421.4 denominator is : 63299.4 Computed efficiency is : 0.164637 +/- 0.00147377
  
  */
  
  extern std::vector<double> brHtaunu_; // BR(H->taunu)
  extern std::vector<double> brHtb_; // BR(H->tb)
  
  enum puconditions{ NOPU, PUPLUS, PUMINUS };
  
  
  
  enum other{
    JETPGID_OLD=1, // get parton jet and associated pgid 
    JETPGID_NEW=32 // get jet flavour directly from jet
  };
  
  
  
  enum inputChannels{
    TTBAR=0, ALLJETS, EJETS, MUJETS, TAUJETS, EE, EMU, ETAU, MUMU, MUTAU, TAUTAU ,
    ZJETS=11, ZTAUTAU=12, WJETS=13, TOPS=14, TOPT=15, TOPW=16,QCD=17,OTHER=18, DATA=20, HIGGS=31, WW=80, WZ=81, ZZ=82,
    WQQ = 100, WENU=101, WMUNU=102, WTAUNU=103, WUNKNOWN =104 
  };
  
  
  enum processLike{ IGNORE=0, ETAU_, MUTAU_, MUMU_, EMU_, EE_, TTBAR_DDBKG_ ,TTBAR_MCBKG_}; 
  
  
  //  enum taudileptonSteps{ 
  //   RECO_STEP1 =0, TRIGGER_STEP1,  MINITREE_STEP1, LEP_STEP1, JET1_STEP1 , JET2_STEP1, JET3_STEP1, MET_STEP1, TAU_STEP1, OS_STEP1, BTAG1_STEP1, BTAG2_STEP1,
  //   RECO_STEP2 =0, TRIGGER_STEP2,  MINITREE_STEP2, LEP_STEP2, JET1_STEP2 , JET2_STEP2, JET3_STEP2, MET_STEP2, BTAG1_STEP2, TAU_STEP2, OS_STEP2, BTAG2_STEP2 
  //  };
  
  enum taudileptonSteps{ 
    // nomt //    RECO_STEP1 =0, TRIGGER_STEP1,  MINITREE_STEP1, LEP_STEP1, JET1_STEP1 , JET2_STEP1, JET3_STEP1, MET_STEP1, MT_STEP1, TAU_STEP1, OS_STEP1, BTAG1_STEP1, R_STEP1, BTAG2_STEP1,
    // nomt //   RECO_STEP2 =0, TRIGGER_STEP2,  MINITREE_STEP2, LEP_STEP2, JET1_STEP2 , JET2_STEP2, JET3_STEP2, MET_STEP2, MT_STEP2, BTAG1_STEP2, TAU_STEP2, OS_STEP2, R_STEP2, BTAG2_STEP2 
    RECO_STEP1 =0, TRIGGER_STEP1,  MINITREE_STEP1, LEP_STEP1, JET1_STEP1 , JET2_STEP1, JET3_STEP1, MET_STEP1, MT_STEP1, TAU_STEP1, OS_STEP1, BTAG1_STEP1, R_STEP1, BTAG2_STEP1,
    RECO_STEP2 =0, TRIGGER_STEP2,  MINITREE_STEP2, LEP_STEP2, JET1_STEP2 , JET2_STEP2, JET3_STEP2, MET_STEP2, MT_STEP2, BTAG1_STEP2, TAU_STEP2, OS_STEP2, R_STEP2, BTAG2_STEP2 
    
  };

  enum dileptonSteps{
    DIL_RECO_STEP1 =0, DIL_TRIGGER_STEP1,  DIL_MINITREE_STEP1, DIL_LEP_STEP1, DIL_JET1_STEP1 , DIL_MT_STEP1, DIL_JET2_STEP1, DIL_JET3_STEP1, DIL_MET_STEP1, DIL_TAU_STEP1, DIL_OS_STEP1, DIL_BTAG1_STEP1, DIL_R_STEP1, DIL_BTAG2_STEP1,
    DIL_RECO_STEP2 =0, DIL_TRIGGER_STEP2,  DIL_MINITREE_STEP2, DIL_LEP_STEP2, DIL_JET1_STEP2 , DIL_MT_STEP2, DIL_JET2_STEP2, DIL_JET3_STEP2, DIL_MET_STEP2, DIL_BTAG1_STEP2, DIL_TAU_STEP2, DIL_OS_STEP2, DIL_R_STEP2, DIL_BTAG2_STEP2 
  };


  
  enum triggermode { STARTING_AT_LJETS_, STARTING_AT_LJETSPLUSMET_};
  
  

  enum channels{ 
    
    // TTBAR Channels /////////////////////////////////////
    ETAU_CH, MUTAU_CH, EE_CH, EMU_CH, MUMU_CH,
    EJETS_CH, MUJETS_CH, TAUJETS_CH, TAUTAU_CH, ALLJETS_CH,  
    ///////////////////////////////////////////////////////
    
    WJETS_CH,
    
    ZJETS_CH, 
    ZTAUTAU_CH,                                               
    
    SINGLETOP_CH, 
    QCD_CH, 
    
    EWKDI_CH,   
    
    WQQ_CH,WENU_CH, WMUNU_CH, WTAUNU_CH, 
    
    CTRLWW_CH,
    HIGGS_CH,
    TBH_CH
    
  };
  
  
  enum taudilCh{ 
    ETAU_TDCH,  MUTAU_TDCH, LEPJETS_TDCH, LEPLEP_TDCH, TAUJETS_TDCH, TAUTAU_TDCH, ALLJETS_TDCH, 
    WJETS_TDCH,
    ZJETS_TDCH, 
    ZTAUTAU_TDCH,                                               
    SINGLETOP_TDCH, 
    QCD_TDCH, 
    EWKDI_TDCH,   
    WQQ_TDCH,WENU_TDCH, WMUNU_TDCH, WTAUNU_TDCH, 
    CTRLWW_TDCH,
    HIGGS_TDCH
  };
  
  
  enum pgid {
    PGID_UNKNOWN  = 0,
    PGID_D        = 1,
    PGID_U        = 2,
    PGID_S        = 3,
    PGID_C        = 4,
    PGID_B        = 5,
    PGID_G        = 21,
    PGID_ELECTRON = 11
  };
  
  
  enum urlCodes { 
    
   DATA_URL = 0, 
   
   TTBAR_URL, 
   
   TT0J_URL, TT1J_URL, TT2J_URL, TT3J_URL, TT4J_URL,
   
   W_ENU_URL, W_MUNU_URL, W_TAUNU_URL,
   
   // DY_10LL_URL, DY_EE_URL, DY_MM_URL, DY_TAUTAU_URL,
   DY_10TO50_URL, DY_FROM50_URL,
   
   S_URL, T_URL, W_URL, A_S_URL, A_T_URL, A_W_URL,
   
   QCD_EM20TO30_URL, 
   QCD_EM30TO80_URL, 
   QCD_EM80TO170_URL,
   
   QCD_BCTOE20TO30_URL,  
   QCD_BCTOE30TO80_URL,  
   QCD_BCTOE80TO170_URL, 
   
   QCD_PHOTON30TO50_URL, 
   QCD_PHOTON50TO80_URL, 
   QCD_PHOTON80TO120_URL, 
   QCD_PHOTON120TO170_URL,
   QCD_MU_URL, 
   
   W0J_URL,
   W1JA_URL,
   W1JB_URL,    
   W2JA_URL,
   W2JB_URL, 
   W3JA_URL,
   W3JB_URL,
   W4JA_URL,      
   W4JB_URL,
   W5JA_URL,
   W5JB_URL,

   WJMADGRAPH_URL,

   WW_URL, WZ_URL, ZZ_URL,
   
   CTRLWW_URL,
   
   HH80_URL,   WH80_URL,
   WH100_URL,  HH100_URL,
   WH120_URL,  HH120_URL,
   WH140_URL,  HH140_URL,
   WH150_URL,  HH150_URL,
   WH155_URL,  HH155_URL,
   WH160_URL,  HH160_URL,
   
   TBH180_URL,  
   TBH190_URL,  
   TBH200_URL,  
   TBH220_URL, 
   TBH250_URL,
   TBH300_URL,
   
   HTB180_URL,  
   HTB190_URL,  
   HTB200_URL,  
   HTB220_URL, 
   HTB250_URL,
   HTB300_URL
  };

  
  
  //RUNNING CONDITIONS ////////////////////////////////////////////////
  //extern bool run2012_         = true;
  extern bool run2012_         ;
  extern bool APPLY_MT_CUT_    ;
  //extern bool APPLY_MT_CUT_    = true;
  extern bool applybtagweight_ ;
  extern bool MODE_ ; //MODE_ = STARTING_AT_LJETSPLUSMET_;
  /*extern bool  eChONmuChOFF_ = true; */ 
  extern bool eChONmuChOFF_ ;
  //extern double LUM_ = 2217.089; //RUN2011A
  // extern double LUM_ =  3929; // /pb  //RUN2011A+RUN2011B
  //  extern double LUM_ = 4064;
  // extern double LUM_ = 4625.21;
  extern double LUM_;
  extern int BTAGIND_; extern double BTAG_CUT_;
  extern TString puFileName_; // Pileup file
  /////////////////////////////////////////////////////////////////////
  
  // I/O areas
  extern TString inputArea_;
  extern TString outputArea_;
  //
  
 // namespace commondefinitions
}  
#endif // _CommonDefinitions_hh
