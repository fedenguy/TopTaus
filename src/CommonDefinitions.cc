#include "LIP/TopTaus/interface/CommonDefinitions.hh" // Include it for the initialization of "MODE", below. FIXME: check whether MODE can be initialized in the physicsAnalyzer or such


// FIXME: honestly, most of this stuff should be moved to cfg files (numerical parameters)
// Sample enums and such can (and should) remain here
// xsecs and their errors should remain here. Or in SampleProcessor.

namespace commondefinitions{
  
  double LUM_ERR            = 0.044; // ICHEP 2012 https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
  
  double TAU_ID_ERR         = 0.06;  
  
  double QCD_CS_ERR         = 0.50;
  
  double TTBAR_CS_ERR       = 0.10;   
  
  double SINGLETOP_CS_ERR   = 0.08;
  double OTHER_CS_ERR       = 0.04;
  
  double LEP_EFF_ERR        = 0.02;        // 3 electrons, 2 muons
  double TAU_MISSID_ERR     = 0.15;        // BULSHIT //0.15;
  
  double WJETS_CS_ERR       = 0.05;
  
  //  double BTAG_EFF_ERR       = 0.10;
  double BTAG_EFF_ERR       = 0.;
  
  double SIGNAL_EFF_ERR     = 0.04;  
  
  double BTAG_eff_R_ = 0.809769;
  double BTAG_eff_F_ = 0.164637;
  
  
  double XSEC_EXP_    = 234;//164.4;       //  was 158, 164.4  2.8 (stat.)  11.9 (syst.)  7.4
  double XSEC_EXP_ERR_= 16./234;  // FIXME: update error        // 14.29/164.4 => exp error, theoretical -> 165 +4/-9 +/-7
  
  // SFb_error taken from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt (2011 prescription)
  
  ///*
  //  float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  //  float ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  //
  //   Tagger: TCHEL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  //   SFb = 0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)));
  //*/
  //   double SFb_error[] = {
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
  double SFb_error[] = {
    0.0188743,
    0.0161816,
    0.0139824,
    0.0152644,
    0.0161226,
    0.0157396,
    0.0161619,
    0.0168747,
    0.0257175,
    0.026424,
    0.0264928,
    0.0315127,
    0.030734,
    0.0438259 };
  
  
  
  /*  
      double BTAG_eff_R_ = 0.812272;
      double BTAG_eff_F_ = 0.17211;
  */
  
  
  /*
  // 30 GeV
  Real Bs : numerator is : 8320.26 denominator is : 10257.8 Computed efficiency is : 0.811112 +/- 0.00278575
  Fake 2 Bs : numerator is : 7633.33 denominator is : 44280 Computed efficiency is : 0.172388 +/- 0.00187882
  
  // 20 GeV
  Real Bs : numerator is : 9635.2 denominator is : 11898.7 Computed efficiency is : 0.809769 +/- 0.00280644
  Fake 2 Bs : numerator is : 10421.4 denominator is : 63299.4 Computed efficiency is : 0.164637 +/- 0.00147377
  
  */
  

  // BR(H->tauNu) for the heavy tbh samples (from 180 to 300 GeV/c^2)
  std::vector<double> brHtaunu_; 
  std::vector<double> brHtb_;
  
//  double ftbh[] ={ // Selection  mu==200&&Xt==2000&&m2==200&&tanb==5
//    0.831852,
//    0.407966,
//    0.228564,
//    0.106446,
//    0.04129783,
//    0.02514242
//  };
    
  
  
  //RUNNING CONDITIONS ////////////////////////////////////////////////
  bool run2012_         = true;
  //  bool run2012_         = false;
  bool APPLY_MT_CUT_    = false;
  //bool APPLY_MT_CUT_    = true;
  bool applybtagweight_ = false;
  bool MODE_ = STARTING_AT_LJETS_; //MODE_ = STARTING_AT_LJETSPLUSMET_;
  /*bool  eChONmuChOFF_ = true; */ 
  bool eChONmuChOFF_ = false;
  double LUM_ = 18072.17; // 2012 ABCD final ntuples
  //  double LUM_ = 11514.17; // 2012 ABC final ntuples
  //  double LUM_ = 5040.32; // 2013 AB final ntuples

  //  double LUM_ = 803.38; // 2012 A final ntuples
  //  double LUM_ = 82.52; // 2012 Ar final ntuples
  //  double LUM_ = 4154.; // 2012 B final ntuples
  //  double LUM_ = 482.26; // 2012 C1 final ntuples
  //  double LUM_ = 5992.; // 2012 C2 final ntuples
  //  double LUM_ = 6558.; // 2012 D final ntuples
  
  int BTAGIND_; double BTAG_CUT_;
  TString puFileName_;
  /////////////////////////////////////////////////////////////////////
  
  // I/O
  TString inputArea_;
  TString outputArea_;
  //////
 // namespace commondefinitions
  
}
