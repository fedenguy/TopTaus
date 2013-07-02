#include "LIP/TopTaus/interface/SampleProcessor.hh"

// System headers
#include <iostream>
#include <iomanip>

// ROOT headers



SampleProcessor::SampleProcessor(double tauPtCut, TString inputArea, TString outputArea):
  sampleTauPtCut_(tauPtCut)
{
  
  commondefinitions::inputArea_ = inputArea;
  commondefinitions::outputArea_ = outputArea;
    
    outFile_=0;
    
    checkleptonsFromTauDecay_ = false;
    //checkleptonsFromTauDecay_ = true;  // used to compute tables with leptons from tau decays (tt->WbWb->tauhtaul)

    //pdfweights_=true;
    pdfweights_=false; // this should allways be false, it will be overwritten by the process_ttbar

    jetpgid_=JETPGID_OLD;

    fast_=true;

    //trigEffStudy_= true;
    trigEffStudy_= false;
   
    pu_ = NOPU;

    //pu_ = PUPLUS;
    //pu_ = PUPLUS;

//    keys_.push_back("TC");
//    keys_.push_back("PF");  
//    keys_.push_back("TaNC0");  // tanc FR    1  %
//    keys_.push_back("TaNC1");  // tanc FR  0.5  %
//    keys_.push_back("TaNC2");  // tanc FR  0.25 %
//    keys_.push_back("TaNC3");  // tanc FR  0.1  %
//    keys_.push_back("HPS0");   // HPS loose iso
//    keys_.push_back("HPS1");   // HPS medium iso
//    keys_.push_back("HPS2");   // HPS high iso
      keys_.push_back("PFlow");  // PF HPS loose iso
 
    char data[10]; sprintf(data,"%dGeV",(int)sampleTauPtCut_);
    
    // output folder and input folder
    
	// Folders at LIP
	
    //oFolder_ = TString("/soft/lip-sw/cmssw/users/nalmeida/new/nalmeida/sw424/Physics/UserCode/nalmeida/TopTauDileptons2012/");  
    //iFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3-B/");  
    //    oFolder_ = TString("/lustre/data3/cmslocal/vischia/tau_dilepton/outputFiles444_3/");
    //    iFolder_ = TString("/lustre/data3/cmslocal/samples/CMSSW_4_4_4/mc/minitrees-muon-Fall11/");

    // Temp 2012
    //    oFolder_ = TString("/lustre/data3/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/");
    oFolder_ = TString(commondefinitions::outputArea_);
    iFolder_ = TString(commondefinitions::inputArea_);
    //    iFolder_ = TString("/lustre/data3/cmslocal/samples/CMSSW_4_4_4/mc/minitrees-muon-Fall11/");
    //    iFolder_ = TString("/lustre/data3/cmslocal/vischia/tau_dilepton/first2011/");

    // Folders at CERN
    //oFolder_ = TString("/afs/cern.ch/work/n/nalmeida/private/CMSSW_5_2_5/src/lipcms/Physics/TopTauDileptons2012/test/");   
	
    //iFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3-B/");
    // if( run2012_ ) iFolder_ = TString("/lustre/lip.pt/data/cmslocal/fnguyen/");   
    // Temp 2012
   //    if( run2012_ ) iFolder_ = TString("/lustre/data3/cmslocal/samples/CMSSW_5_3_7_patch4/");   

    if( pu_ == PUPLUS  ) oFolder_ += TString("puplus-"); else if ( pu_ == PUMINUS ) oFolder_ += TString("puminus-");

    TString mtStr("nomt-");    if ( APPLY_MT_CUT_)    { mtStr  = TString("mt-");}
    TString lepStr("EL-");     if (!eChONmuChOFF_)    { lepStr = TString("MU-");}
    TString bwStr("");         if (!bwStr)            { bwStr  = TString("bw-");}
    TString vers("2011-V1-");

   
    if( run2012_ ) vers = TString("2012-V1-");


    oFolder_ += mtStr+vers+TString("mc-")+lepStr+bwStr+TString(data)+TString("/");
    

    if( eChONmuChOFF_ ) iDFolder_ = TString(commondefinitions::inputArea_);  
    //else              iDFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/data/mTrees-v3/");  
    //   else                iDFolder_ = TString("/lustre/data3/cmslocal/samples/CMSSW_4_4_4/data/minitrees-muon-2011/");  
    // Temp 2012
    else                iDFolder_ = TString(commondefinitions::inputArea_);  


    // LIP
    //if( run2012_ ) iDFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_8_patch4b/data/minitrees-muon-2011-newTau/");  
    //oDFolder_  = TString("/afs/cern.ch/work/n/nalmeida/private/CMSSW_5_2_5/src/lipcms/Physics/TopTauDileptons2012/test");

    // Temp 2012
    if( run2012_ ) iDFolder_ = TString(commondefinitions::inputArea_);  
    //if( run2012_ ) iDFolder_ = TString("/lip-sw/cmssw/users/nalmeida/new/CMSSW_5_2_5/src/lipcms/Physics/TopTauDileptons2012/mt-2012-V1-data-MU-20GeV/");  
    //    oDFolder_  = TString("/lustre/data3/cmslocal/vischia/tau_dilepton/outputFiles444_3/");
    // Temp 2012
    //   oDFolder_  = TString("/lustre/data3/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/");
    oDFolder_  = TString(commondefinitions::outputArea_);


// Fall 2010 samples   : https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionFall2010
// NLO cross sections  : https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries#crosssections
// NNLO cross sections : https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf


// Summer11
//http://cms.cern.ch/iCMS/prep/requestmanagement?status=Done&pwg=HIG&campid=Summer11_R1&dsn=
//http://cms.cern.ch/iCMS/prep/requestmanagement?pwg=TOP&campid=Summer11_R1
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections




    oDFolder_ += mtStr+vers+TString("data-")+lepStr+TString(data)+TString("/");
    
    
    for (int i=0; i<200; i++ ){ defaultXSections_.push_back(0); }
    
    if( run2012_ ) {
      defaultXSections_[DATA_URL]      = 0.;         
      defaultXSections_[TTBAR_URL]     = 234.;  //Kidonakis - approx. NNLO - Inclusive - 234+10-)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))) - Top mass: m(top)=173GeV, Default scales: muF=muR=mt, PDF: MSTW2008 NNLO 
      
      defaultXSections_[W_ENU_URL]     = 12503.0; //NNLO 37509.0  36257.2/3 = 12085.7 - NNLO - W->lv, l=e,m,t - Inclusive W production, BR(W->lv) included, l=e,m,t, PDF error also includes alphas 
      defaultXSections_[W_MUNU_URL]    = 12503.0;
      defaultXSections_[W_TAUNU_URL]   = 12503.0; 

      defaultXSections_[WJMADGRAPH_URL]  = 37509.0; // 36257.2;// 30400; //36257.2;
      
      // for the _filter sample:     xsec: 11050.0 eff: 0.069 
      defaultXSections_[DY_10TO50_URL]          = 11050.0; // from http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=DYJetsToLL_M-10To50filter_8TeV-madgraph&campid=Summer12_DR53X
      defaultXSections_[DY_FROM50_URL]          = 2950.; // from http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball&campid=Summer12_DR53X

      defaultXSections_[S_URL]                  = 3.79;     
      defaultXSections_[T_URL]                  = 56.4;    
      defaultXSections_[W_URL]                  = 11.1;  
      
      
      defaultXSections_[A_S_URL]                = 1.76;    
      defaultXSections_[A_T_URL]                = 30.7;    
      defaultXSections_[A_W_URL]                = 11.1;

      //fn defaultXSections_[QCD_EM20TO30_URL]       = 2454400.; not considered!!!
      defaultXSections_[QCD_EM30TO80_URL]       = 4615893.; // 74330000. eff 0.0621 http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6&campid=Summer12
      defaultXSections_[QCD_EM80TO170_URL]      = 183294.9; // 1191000. eff 0.1539 http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6&campid=Summer12
      
      //fn defaultXSections_[QCD_BCTOE20TO30_URL]    = 132160; not considered!!! 
      defaultXSections_[QCD_BCTOE30TO80_URL]    = 167040.; // 74240000. eff 0.00225 http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=QCD_Pt_30_80_BCtoE_TuneZ2star_8TeV_pythia6&campid=Summer12
      defaultXSections_[QCD_BCTOE80TO170_URL]   = 12981.9; //1191000. 0.0109 http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=QCD_Pt_80_170*&campid=Summer12         
      
      defaultXSections_[QCD_MU_URL]             = 134680.; //364000000. eff 3.7E-4 http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6&campid=Summer12
      
      
      defaultXSections_[QCD_PHOTON30TO50_URL]   = 19931.62;        
      defaultXSections_[QCD_PHOTON50TO80_URL]   = 3322.309;         
      defaultXSections_[QCD_PHOTON80TO120_URL]  = 558.2865;        
      defaultXSections_[QCD_PHOTON120TO170_URL] = 108.0068;
      
      //fn xsections from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
      defaultXSections_[WW_URL]                 = 54.838;
      defaultXSections_[WZ_URL]                 = 33.72;
      defaultXSections_[ZZ_URL]                 = 17.627;
      
      defaultXSections_[CTRLWW_URL]             = 234; 
      defaultXSections_[WH80_URL]               = 234;            
      defaultXSections_[HH80_URL]               = 234;                                           
      defaultXSections_[WH100_URL]              = 234;                                          
      defaultXSections_[HH100_URL]              = 234;                                        
      defaultXSections_[WH120_URL]              = 234;                                          
      defaultXSections_[HH120_URL]              = 234;                                          
      defaultXSections_[WH140_URL]              = 234;                                          
      defaultXSections_[HH140_URL]              = 234; 
      defaultXSections_[WH150_URL]              = 234; 
      defaultXSections_[HH150_URL]              = 234; 
      defaultXSections_[WH155_URL]              = 234; 
      defaultXSections_[HH155_URL]              = 234; 
      defaultXSections_[WH160_URL]              = 234;                                           
      defaultXSections_[HH160_URL]              = 234;                                                                                 
      
      // 1.1 is for plots and tables: for plots it's already implemented in the data/plotter/samples.xml. For tables must be implemented via switchable parameter in cfg file. 
      // 234 is for datacards.
      defaultXSections_[TBH180_URL]             = 234. * 0.831852  ; //1.1; //234;
      defaultXSections_[TBH190_URL]             = 234. * 0.407966  ; //1.1; //234;  
      defaultXSections_[TBH200_URL]             = 234. * 0.228564  ; //1.1; //234;  
      defaultXSections_[TBH220_URL]             = 234. * 0.106446  ; //1.1; //234; 
      defaultXSections_[TBH250_URL]             = 234. * 0.04129783; //1.1; //234;
      defaultXSections_[TBH300_URL]             = 234. * 0.02514242; //1.1; //234;
    
      defaultXSections_[HTB180_URL]             = 234. * (1 - 0.831852  ); //1.1; //234;
      defaultXSections_[HTB190_URL]             = 234. * (1 - 0.407966  ); //1.1; //234;  
      defaultXSections_[HTB200_URL]             = 234. * (1 - 0.228564  ); //1.1; //234;  
      defaultXSections_[HTB220_URL]             = 234. * (1 - 0.106446  ); //1.1; //234; 
      defaultXSections_[HTB250_URL]             = 234. * (1 - 0.04129783); //1.1; //234;
      defaultXSections_[HTB300_URL]             = 234. * (1 - 0.02514242); //1.1; //234;
    
	
    } else {
      
      defaultXSections_[DATA_URL]      = 0.;         
      defaultXSections_[TTBAR_URL]     = 165;  //NNLO aproximation +- 10 old was 157.5 
      
      defaultXSections_[W_ENU_URL]     = 10438; //NNLO 31314/3 =10438 //NLO 28000/3 = 9333 // xsec 7899 WToENu_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1                 
      defaultXSections_[W_MUNU_URL]    = 10438;
      defaultXSections_[W_TAUNU_URL]   = 10438; 
      
      
      // ALPGEN SAMPLES
      /*
	W0Jets           /W0Jets_TuneD6T_7TeV-alpgen-tauola                                                    2.003e+04
	W1Jets_0to100    /W1Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO 	    3.693e+03
	W1Jets_100to300  /W1Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO     7.197e+01
	W2Jets_0to100    /W2Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO 	    9.434e+02
	W2Jets_100to300  /W2Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO     6.718e+01
	W3Jets_0to100    /W3Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO 	    2.087e+02
	W3Jets_100to300  /W3Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO     3.243e+01
	W4Jets_0to100    /W4Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO 	    4.446e+01
	W4Jets_100to300  /W4Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO     1.138e+01
	W5Jets_0to100    /W5Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO 	    1.111e+01
	W5Jets_100to300  /W5Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO     3.789e+00 
      */
      
      // total 
      // w0 ->      2.024e+04 //Z2
      // W0 ->D6    2.003e+04
      
      defaultXSections_[W0J_URL]         = 20030.*1.3;
      defaultXSections_[W1JA_URL]        = 3693.*1.3;
      defaultXSections_[W1JB_URL]        = 71.97*1.3;
      defaultXSections_[W2JA_URL]        = 943.4*1.3;
      defaultXSections_[W2JB_URL]        = 67.18*1.3;
      defaultXSections_[W3JA_URL]        = 208.7*1.3;
      defaultXSections_[W3JB_URL]        = 32.43*1.3;
      defaultXSections_[W4JA_URL]        = 44.46*1.3;
      defaultXSections_[W4JB_URL]        = 11.38*1.3;
      defaultXSections_[W5JA_URL]        = 11.11*1.3;
      defaultXSections_[W5JB_URL]        = 3.789*1.3;
      
      defaultXSections_[WJMADGRAPH_URL]  = 31314;
      
      defaultXSections_[TT0J_URL]        = 49.17*2.1768;
      defaultXSections_[TT1J_URL]        = 17.97*2.1768 ;
      defaultXSections_[TT2J_URL]        = 4.279*2.1768 ;
      defaultXSections_[TT3J_URL]        = 0.7951*2.1768 ;
      defaultXSections_[TT4J_URL]        = 0.1404*2.1768 ;
      
      
      /*
	defaultXSections_[DY_HIGHLL_URL]     = 3048;
	defaultXSections_[DY_EE_URL]         = 650;
	defaultXSections_[DY_MM_URL]         = 650;
	defaultXSections_[DY_TAUTAU_URL]     = 650;
      */
      
      //defaultXSections_[DY_10TO50_URL]          = 310;  //NLO 310 // xsec = 310  eff=0.18 DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola (on madgraph we do not apply eff)

      // the number is LO cross section (3x310 pb) with a k-factor to correct for NLO (1.33). The actual source of this - I don't know.
      defaultXSections_[DY_10TO50_URL]          = 12369;
      defaultXSections_[DY_FROM50_URL]          = 3048; //NNLO 3048 //NLO 2800// xsec= 2321  eff=0.44 DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola 
      
      
      /* Fall 2010 samples
	 defaultXSections_[S_URL]                  = 1.49;     
	 defaultXSections_[T_URL]                  = 20.93 ;    
	 defaultXSections_[W_URL]                  = 10.6;   
      */
      
      
      

      /*
	
      t+ (s-channel) 	MCFM 	-- 	NLO 	Inclusive 	2.72 	+0.07-0.06 (±0.08), Total = +0.11 - 0.10 	Inclusive top decay, PDF errors also includes alphas
      t- (s-channel) 	MCFM 	-- 	NLO 	Inclusive 	1.49 	+0.04-0.03 (±0.08), Total = +0.09-0.08 	Inclusive top decay, PDF errors also includes alphas 
      
      t+ (t-channel) 	MCFM 	-- 	NLO 	Inclusive 	42.6 	+0.9-0.8 (±2.2), Total = +2.4-2.3 	Inclusive top decay, PDF errors also includes alphas
      t- (t-channel) 	MCFM 	-- 	NLO 	Inclusive 	22.0 	+0.6-0.3 (±0.8), Total = +0.10-0.8 	Inclusive top decay, PDF errors also includes alphas 
      
      W+tbar   	MCFM 	-- 	NLO 	Inclusive 	5.3 	±0.6 	--
      W-t 	        MCFM 	-- 	NLO 	Inclusive 	5.3 	±0.6 	-- 
      */
      
      
      
      
      defaultXSections_[S_URL]                  = 2.72;     
      defaultXSections_[T_URL]                  = 42.6;    
      //defaultXSections_[W_URL]                = 5.3;   //WARNING :: This is old
      defaultXSections_[W_URL]                  = 7.87;  
      
      
      defaultXSections_[A_S_URL]                = 1.49;    
      defaultXSections_[A_T_URL]                = 22;    
      //defaultXSections_[A_W_URL]              = 5.3;   //WARNING :: use the new value
      defaultXSections_[A_W_URL]                = 7.87;  
      

      
      
      
      defaultXSections_[QCD_EM20TO30_URL]       = 2454400.;     //xsec =236000000 eff=0.0104 QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6
      defaultXSections_[QCD_EM30TO80_URL]       = 3866200.;     //xsec =59480000  eff=0.065  QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia6 
      defaultXSections_[QCD_EM80TO170_URL]      = 139500.;      //xsec =900000    eff=0.155  QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6 
      
      defaultXSections_[QCD_BCTOE20TO30_URL]    = 132160;       //xsec = 236000000 eff=0.00056 QCD_Pt-20to30_BCtoE_TuneZ2_7TeV-pythia6 
      defaultXSections_[QCD_BCTOE30TO80_URL]    = 136804;       //xsec = 59480000  eff=0.00230 QCD_Pt-30to80_BCtoE_TuneZ2_7TeV-pythia6 
      defaultXSections_[QCD_BCTOE80TO170_URL]   = 9360;         //xsec = 900000    eff=0.0104  QCD_Pt-80to170_BCtoE_TuneZ2_7TeV-pythia6   
      
      defaultXSections_[QCD_PHOTON30TO50_URL]   = 16690;        //xsec = 1.669e+04   G_Pt_30to50_TuneZ2_7TeV_pythia6                 
      defaultXSections_[QCD_PHOTON50TO80_URL]   = 2722;         //xsec = 2.722e+03   G_Pt_50to80_TuneZ2_7TeV_pythia6                         
      defaultXSections_[QCD_PHOTON80TO120_URL]  = 447.2;        //xsec = 4.472e+02   G_Pt_80to120_TuneZ2_7TeV_pythia6
      defaultXSections_[QCD_PHOTON120TO170_URL] = 84.1;	        //xsec = 8.417e+01   G_Pt_120to170_TuneZ2_7TeV_pythia6
      
      //defaultXSections_[QCD_MU_URL]           = 349988;       //xsec=296600000 	eff=0.00118  QCD_Pt-20_MuEnrichedPt-10_TuneZ2_7TeV-pythia6 //79688.;                                    
      defaultXSections_[QCD_MU_URL]             = 84679;        //xsec=296600000 	eff=0.002855 QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6      
      

      // xsections taken from CMS-AN-10-337
      defaultXSections_[WW_URL]                 = 43;
      defaultXSections_[WZ_URL]                 = 18.2;
      defaultXSections_[ZZ_URL]                 = 5.9;
      
      defaultXSections_[CTRLWW_URL]             = 165; 
      defaultXSections_[WH80_URL]               = 165;            
      defaultXSections_[HH80_URL]               = 165;                                           
      defaultXSections_[WH100_URL]              = 165;                                          
      defaultXSections_[HH100_URL]              = 165;                                        
      defaultXSections_[WH120_URL]              = 165;                                          
      defaultXSections_[HH120_URL]              = 165;                                          
      defaultXSections_[WH140_URL]              = 165;                                          
      defaultXSections_[HH140_URL]              = 165; 
      defaultXSections_[WH150_URL]              = 165; 
      defaultXSections_[HH150_URL]              = 165; 
      defaultXSections_[WH155_URL]              = 165; 
      defaultXSections_[HH155_URL]              = 165; 
      defaultXSections_[WH160_URL]              = 165;                                           
      defaultXSections_[HH160_URL]              = 165;                                                                                 
      
      defaultXSections_[TBH250_URL]                = 165;           
    }
  }




void SampleProcessor::init(vector<TString> listOfurls, vector<double> listOfXSections, bool fullStats, double lum, uint n_data){

  
  /*
    
  // QCD MVA
  qcdMVA_   = new MVAComputer("/u/home/cms/vischia/public/tau_dilepton/mvaWeights_428TrainedDiJetDataPFlow_pt20.mva");
  
  // W plus jets MVA
  wjetsMVA_ =  new MVACompute("/u/home/cms/vischia/public/tau_dilepton/mvaWeights_428TrainedWMuDataPFlow_pt20.mva");
  
  */
  
  nToProcess_ = n_data;
 
  if(listOfurls.size()!= listOfXSections.size() ){ 
    cout<<endl<<"\n\n Number of files does not match number of xsections... \n\n\n"<<endl; 
    exit(0);
  }
  
  // Reset urls, xsections,files,histos, scales,reades,events ////////
  listOfurls_.clear();
  listOfXSections_.clear();
  listOfFiles_.clear();
  listOfHistos_.clear();
  listOfScales_.clear();
  listOfReaders_.clear();
  listOfReadersMC_.clear();
  //listOfSpies_.clear();

  listOfEventsToProcess_.clear();
  listOfEventsToProcessMC_.clear();

  listOfEvents_.clear();
 
  
  URLS_=listOfurls.size();

  for(unsigned int i=0; i< URLS_;i++){ 

    listOfurls_.push_back(listOfurls[i]); 
    listOfXSections_.push_back(listOfXSections[i]); 
    listOfFiles_.push_back(0); 
    listOfHistos_.push_back(0); 
    listOfScales_.push_back(0); 
    listOfReaders_.push_back(0);
    listOfReadersMC_.push_back(0);  
    //listOfSpies_.push_back(0);
    listOfEventsToProcess_.push_back(0); 
    listOfEventsToProcessMC_.push_back(0); 
    listOfEvents_.push_back(0); 
  }
  ///////////////////////////////////////////////////////////////////

  fullStats_ = fullStats;
  lum_       = lum;

  init();
  
}


void SampleProcessor::initProcessedFiles(){

  TString dataFolder  ("xxx"); 	 	 
  TString mcFolder    ("xxx");   

  listOfurls_.clear(); 
  listOfurls_.push_back(dataFolder+TString("out-data.root"));
  listOfurls_.push_back(mcFolder+TString("out-ttbar.root"));
  listOfurls_.push_back(mcFolder+TString("out-singletop.root"));  
  listOfurls_.push_back(mcFolder+TString("out-wjets.root"));  
  listOfurls_.push_back(mcFolder+TString("out-zjets.root"));  
  listOfurls_.push_back(mcFolder+TString("out-dibosons.root"));  
  listOfurls_.push_back(mcFolder+TString("out-qcd.root"));  

  listOfProcessedFiles_.clear();

  // open files ///////////////////////////////////////////////////////////////////////////////////////////////
  for(unsigned int i=0;i <listOfurls_.size();i++){
    listOfProcessedFiles_.push_back( TFile::Open(listOfurls_[i]) );
    if(! listOfProcessedFiles_[i] ){ cout <<endl<<"\n Can't find file : "<<(listOfurls_[i])<<endl; continue;}  
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

}



void SampleProcessor::init(){


  // data location in the tree ///////////////////////
  TString cutflow("myEventSelector/Selection/cutflow");
  TString tag("data"); 
  TString id("myEventSelector");
  ////////////////////////////////////////////////////


  // open files ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(unsigned int i=0;i <URLS_;i++){
    listOfFiles_[i] = TFile::Open(listOfurls_[i]);  
    if(! listOfFiles_[i] ){ cout <<endl<<"\n Can't find file : "<<(listOfurls_[i])<<endl; continue;}

    listOfHistos_[i] = (TH1D *) (listOfFiles_[i])->Get(cutflow);
    listOfEvents_[i] = (listOfHistos_[i])->GetBinContent(1);
    listOfScales_[i] = (lum_*listOfXSections_[i])/listOfEvents_[i];   

    if( listOfXSections_[i] == 0){ listOfScales_[i]=1;}  // Data

    listOfReaders_[i]     = new event::Reader(); 
    listOfReadersMC_[i]   = new event::Reader();

    listOfEventsToProcess_[i]   = (listOfReaders_[i])->AssignMiniEventTreeFrom(listOfFiles_[i],tag,id);
    listOfEventsToProcessMC_[i] = (listOfReadersMC_[i])->AssignMiniEventTreeFrom(listOfFiles_[i],"mc",id);


    //SPY ///////////////////////////////////////////////////////////
    TString extra;
    if(eChONmuChOFF_) extra=TString("_spy_eltau.root");
    else extra=TString("_spy_mutau.root");
    TString spyName = outFileName_ + extra;
    //TTree * spyTree = listOfReaders_[i]->PrepareToSpyEvents(spyName);
    //cout<<endl<<"Spy Events ... on file ... "<<spyName<<endl;
    //////////////////////////////////////////////////////////////////
    
    if(fullStats_){
      cout<<endl<<"\n processing file  : "<<(listOfurls_[i])
          <<"\n Tau pt cut was         : "<<setw(12)<<setiosflags(ios::fixed) << setprecision(8)<<sampleTauPtCut_
          <<"\n sample -> scale        : "<<setw(12)<<setiosflags(ios::fixed) << setprecision(8)<<(listOfScales_[i])
          <<"\n evt processed          : "<<setw(8)<<setprecision(1)<<(listOfEvents_[i])
          <<"\n evt in minitree        : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcess_[i])
          <<"\n MC evt in minitree     : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcessMC_[i])
          <<"\n Xsection employed      : "<<(listOfXSections_[i])
          <<"\n Normalized to lumi (pb): "<<LUM_<<endl;
      
      infoFile_<<endl<<"\n processing file : "<<(listOfurls_[i])
	       <<"\n Tau pt cut was      : "<<setw(12)<<setiosflags(ios::fixed) << setprecision(8)<<sampleTauPtCut_
	       <<"\n sample -> scale     : "<<setw(12)<<setiosflags(ios::fixed) << setprecision(8)<<(listOfScales_[i])
	       <<"\n evt processed       : "<<setw(8)<<setprecision(1)<<(listOfEvents_[i])
	       <<"\n evt in minitree     : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcess_[i])
	       <<"\n MC evt in minitree  : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcessMC_[i])
	       <<"\n Xsection employed : "<<(listOfXSections_[i])
	       <<"\n Normalized to lumi (pb): "<<LUM_<<endl;
      
      if(trigEffStudy_)   {  infoFile_<<" Trigger efficiency is ON  "; cout<<" Trigger efficiency is ON ";  }
      else                {  infoFile_<<" Trigger efficiency is OFF "; cout<<" Trigger efficiency is OFF "; }
      if(applybtagweight_){  infoFile_<<" BTAG weights are ON  ";      cout<<" BTAG weights are ON ";       }
      else                {  infoFile_<<" BTAG weights are OFF ";      cout<<" BTAG weights are OFF";       }



    }else{ 

      cout<<endl<<"\n processing file : "<<(listOfurls_[i])<<"\n evt in minitree     : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcess_[i]);
      infoFile_<<"\n processing file : "<<(listOfurls_[i])<<"\n evt in minitree     : "<<setw(8)<<setprecision(1)<<(listOfEventsToProcess_[i]);
 
      if( nToProcess_ > listOfEventsToProcess_[i] ){ 
        cout<<endl<<" WARNING :: required number of events too large : "<<nToProcess_<<" sample will be fully processed..."<<endl;
        infoFile_<<"  WARNING :: required number of events too large : "<<nToProcess_<<" sample will be fully processed..."<<endl;
   
      }else {      
        cout<<" but only "<<nToProcess_<<" will be processed... "<<endl;
        infoFile_<<" but only "<<nToProcess_<<" will be processed... "<<endl;

        listOfScales_[i] = (lum_*listOfXSections_[i])/nToProcess_;   
        listOfEventsToProcess_[i] = nToProcess_;
      }

      if(trigEffStudy_)   {  infoFile_<<" Trigger efficiency is ON ";  cout<<" Trigger efficiency is ON ";  }
      else                {  infoFile_<<" Trigger efficiency is OFF "; cout<<" Trigger efficiency is OFF "; } 
      if(commondefinitions::applybtagweight_){  infoFile_<<" BTAG weights are ON  ";      cout<<" BTAG weights are ON ";       }
      else                {  infoFile_<<" BTAG weights are OFF ";      cout<<" BTAG weights are OFF";       }

    }
  }// finish processing urls /////////////////////////////////////////////////////////////////////////////////////////////////////////
  

}



int SampleProcessor::codeChannel(int i, int urlCode){

  int ret(666);
  if( 
    urlCode == WH80_URL  || urlCode == HH80_URL  ||
    urlCode == WH100_URL || urlCode == HH100_URL ||
    urlCode == WH120_URL || urlCode == HH120_URL || 
    urlCode == WH140_URL || urlCode == HH140_URL || 
    urlCode == WH150_URL || urlCode == HH150_URL || 
    urlCode == WH155_URL || urlCode == HH155_URL ||
    urlCode == WH160_URL || urlCode == HH160_URL 
  )                                          {  ret = HIGGS_CH;     }
  else if( 
	  urlCode == TBH180_URL || urlCode == HTB180_URL || 
	  urlCode == TBH190_URL || urlCode == HTB190_URL || 
	  urlCode == TBH200_URL || urlCode == HTB200_URL || 
	  urlCode == TBH220_URL || urlCode == HTB220_URL || 
	  urlCode == TBH250_URL || urlCode == HTB250_URL || 
	  urlCode == TBH300_URL	|| urlCode == HTB300_URL
	  ){  ret = TBH_CH;       }                   
  else if( i == ETAU                        ){  ret = ETAU_CH;      } 
  else if( i == MUTAU                       ){  ret = MUTAU_CH;     } 
  else if( i == EJETS                       ){  ret = EJETS_CH;     }
  else if( i == MUJETS                      ){  ret = MUJETS_CH;    }
  else if( i == EE                          ){  ret = EE_CH;        }
  else if( i == MUMU                        ){  ret = MUMU_CH;      }
  else if( i == EMU                         ){  ret = EMU_CH;       }
  else if( i == TAUJETS                     ){  ret = TAUJETS_CH;   }
  else if( i == TAUTAU                      ){  ret = TAUTAU_CH;    }
  else if( i == ALLJETS                     ){  ret = ALLJETS_CH;   }
  else if( i == WW || i == WZ || i== ZZ     ){  ret = EWKDI_CH;     }
  else if( i == WQQ                         ){  ret = WQQ_CH;       }
  else if( i == WENU                        ){  ret = WENU_CH;      }
  else if( i == WMUNU                       ){  ret = WMUNU_CH;     }
  else if( i == WTAUNU                      ){  ret = WTAUNU_CH;    }  
  else if( i == ZJETS                       ){  ret = ZJETS_CH;     } 
  else if( i == ZTAUTAU                     ){  ret = ZTAUTAU_CH;   } 
  else if( i == TOPS                        ){  ret = SINGLETOP_CH; } 
  else if( i == TOPT                        ){  ret = SINGLETOP_CH; } 
  else if( i == TOPW                        ){  ret = SINGLETOP_CH; } 
  else if( i == QCD                         ){  ret = QCD_CH;       } 

  return ret;

}



int SampleProcessor::tdChannel(int i){
  
  int ret(666);
  if     ( i == ETAU_CH                              ){  ret = ETAU_TDCH;      } 
  else if( i == MUTAU_CH                             ){  ret = MUTAU_TDCH;     } 
  else if( i == EJETS_CH || i == MUJETS_CH           ){  ret = LEPJETS_TDCH;   }
  else if( i == EE_CH || i == EMU_CH || i == MUMU_CH ){  ret = LEPLEP_TDCH;    }
  else if( i == TAUJETS_CH                           ){  ret = TAUJETS_TDCH;   }
  else if( i == TAUTAU_CH                            ){  ret = TAUTAU_TDCH;    }
  else if( i == ALLJETS_CH                           ){  ret = ALLJETS_TDCH;   }
  else if( i == EWKDI_CH                             ){  ret = EWKDI_TDCH;     }
  else if( i == WQQ_CH                               ){  ret = WQQ_TDCH;       }
  else if( i == WENU_CH                              ){  ret = WENU_TDCH;      }
  else if( i == WMUNU_CH                             ){  ret = WMUNU_TDCH;     }
  else if( i == WTAUNU_CH                            ){  ret = WTAUNU_TDCH;    }  
  else if( i == WJETS_CH                             ){  ret = WJETS_TDCH;     }
  else if( i == ZJETS_CH                             ){  ret = ZJETS_TDCH;     } 
  else if( i == ZTAUTAU_CH                           ){  ret = ZTAUTAU_TDCH;   } 
  else if( i == SINGLETOP_CH                         ){  ret = SINGLETOP_TDCH; } 
  else if( i == QCD_CH                               ){  ret = QCD_TDCH;       } 

  return ret;

}





void SampleProcessor::process_ttbar(){

  //APLGEN TTBAR
  //process(false, TT0J_URL, iFolder_ + TString("TT0Jets_Alpgen.root"), oFolder_ + TString("out-TT0_Alpgen.root",keys_); 
  //process(false, TT1J_URL, iFolder_ + TString("TT1Jets_Alpgen.root"), oFolder_ + TString("out-TT1_Alpgen.root",keys_); 
  //process(false, TT2J_URL, iFolder_ + TString("TT2Jets_Alpgen.root"), oFolder_ + TString("out-TT2_Alpgen.root",keys_); 
  //process(false, TT3J_URL, iFolder_ + TString("TT3Jets_Alpgen.root"), oFolder_ + TString("out-TT3_Alpgen.root",keys_); 
  //process(false, TT4J_URL, iFolder_ + TString("TT4Jets_Alpgen.root"), oFolder_ + TString("out-TT4_Alpgen.root",keys_); 
  //MADGRAPH
  //process(false, TTBAR_URL, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar-madgraph.root"),keys_);                          
  //process(false, TTBAR_URL, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar-signal-madgraph.root"),keys_, TTBAR_SIGNAL_ );   
  //process(false, TTBAR_URL, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar-other-madgraph.root"),keys_, TTBAR_OTHER_ );     

  url_= TTBAR_URL; 


  //NORMAL SAMPLE
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1.root"),keys_);
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2.root"),keys_);
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar.root"),keys_);

  //process(false, url_, iFolder_ + TString("AllTTJets_TuneZ2star_8TeV-madgraph-tauola.root"), oFolder_+TString("out-ttbar.root"),keys_);
  
  
  

  //SAMPLE WITH WEIGHTS
  //pdfweights_=true; process(false, url_, iFolder_ + TString("ttbar-pdf.root"), oFolder_+TString("out-ttbar-pdf.root"),keys_); pdfweights_=false;
  //pdfweights_=true; process(false, url_, iFolder_ + TString("ttbar-pdf-nofilter-B.root"), oFolder_+TString("out-ttbar-pdf.root"),keys_); pdfweights_=false;
  ///lustre/data3/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3-B/ttbar-pdf-nofilter-B.root


  // Trigger efficiencies study ////////////////////////////////////////////////////////////////////////
  //iFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3/");
  //process(false, url_, iFolder_ + TString("trig-ttbar.root"), oFolder_+TString("out-ttbar.root"),keys_);
  //////////////////////////////////////////////////////////////////////////////////////////////////////    





}




void SampleProcessor::process_ttbar_unc(){

  url_= TTBAR_URL;
/*
  process(false, url_, iFolder_ + TString("ttbar-scaleup_424.root"),   oFolder_+TString("out-ttbar-scaleup.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-scaledown_424.root"), oFolder_+TString("out-ttbar-scaledown.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-matchup_424.root"),   oFolder_+TString("out-ttbar-pdfup.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-matchdown_424.root"), oFolder_+TString("out-ttbar-pdfdown.root"),keys_);   

  process(false, url_, iFolder_ + TString("ttbar-m161p5_424.root"),   oFolder_+TString("out-ttbar-m161p5.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-m166p5_424.root"),   oFolder_+TString("out-ttbar-m166p5.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-m169p5_424.root"),   oFolder_+TString("out-ttbar-m169p5.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-m175p5_424.root"),   oFolder_+TString("out-ttbar-m175p5.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-m178p5_424.root"),   oFolder_+TString("out-ttbar-m178p5.root"),keys_);    
  process(false, url_, iFolder_ + TString("ttbar-m184p5_424.root"),   oFolder_+TString("out-ttbar-m184p5.root"),keys_);  

*/

  // this is negleted for now 
/*
  process(false, url_, iFolder_ + TString("ttbar-m161p5_424.root"),   oFolder_+TString("out-ttbar-m161p5_etau.root"),keys_,ETAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m166p5_424.root"),   oFolder_+TString("out-ttbar-m166p5_etau.root"),keys_,ETAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m169p5_424.root"),   oFolder_+TString("out-ttbar-m169p5_etau.root"),keys_,ETAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m175p5_424.root"),   oFolder_+TString("out-ttbar-m175p5_etau.root"),keys_,ETAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m178p5_424.root"),   oFolder_+TString("out-ttbar-m178p5_etau.root"),keys_,ETAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m184p5_424.root"),   oFolder_+TString("out-ttbar-m184p5_etau.root"),keys_,ETAU_);  
  
  process(false, url_, iFolder_ + TString("ttbar-m161p5_424.root"),   oFolder_+TString("out-ttbar-m161p5_mutau.root"),keys_,MUTAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m166p5_424.root"),   oFolder_+TString("out-ttbar-m166p5_mutau.root"),keys_,MUTAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m169p5_424.root"),   oFolder_+TString("out-ttbar-m169p5_mutau.root"),keys_,MUTAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m175p5_424.root"),   oFolder_+TString("out-ttbar-m175p5_mutau.root"),keys_,MUTAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m178p5_424.root"),   oFolder_+TString("out-ttbar-m178p5_mutau.root"),keys_,MUTAU_);    
  process(false, url_, iFolder_ + TString("ttbar-m184p5_424.root"),   oFolder_+TString("out-ttbar-m184p5_mutau.root"),keys_,MUTAU_);    

  process(false, url_, iFolder_ + TString("ttbar-m161p5_424.root"),   oFolder_+TString("out-ttbar-m161p5_other.root"),keys_,TTBAR_OTHER_);    
  process(false, url_, iFolder_ + TString("ttbar-m166p5_424.root"),   oFolder_+TString("out-ttbar-m166p5_other.root"),keys_,TTBAR_OTHER_);    
  process(false, url_, iFolder_ + TString("ttbar-m169p5_424.root"),   oFolder_+TString("out-ttbar-m169p5_other.root"),keys_,TTBAR_OTHER_);    
  process(false, url_, iFolder_ + TString("ttbar-m175p5_424.root"),   oFolder_+TString("out-ttbar-m175p5_other.root"),keys_,TTBAR_OTHER_);    
  process(false, url_, iFolder_ + TString("ttbar-m178p5_424.root"),   oFolder_+TString("out-ttbar-m178p5_other.root"),keys_,TTBAR_OTHER_);    
  process(false, url_, iFolder_ + TString("ttbar-m184p5_424.root"),   oFolder_+TString("out-ttbar-m184p5_other.root"),keys_,TTBAR_OTHER_);    

*/
}



void SampleProcessor::process_ttbar_mutau(){

  url_= TTBAR_URL;


  //  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_etau.root"),keys_,  ETAU_  );    
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1_mutau.root"),keys_, MUTAU_ ); 
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2_mutau.root"),keys_, MUTAU_ ); 
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_mutau.root"),keys_, MUTAU_ ); 
}

void SampleProcessor::process_ttbar_mumu(){

  url_= TTBAR_URL;


  //  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_etau.root"),keys_,  ETAU_  );    
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1_mutau.root"),keys_, MUTAU_ ); 
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2_mutau.root"),keys_, MUTAU_ ); 
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_mumu.root"),keys_, MUMU_ ); 
}

void SampleProcessor::process_ttbar_emu(){

  url_= TTBAR_URL;


  //  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_etau.root"),keys_,  ETAU_  );    
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1_mutau.root"),keys_, MUTAU_ ); 
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2_mutau.root"),keys_, MUTAU_ ); 
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_emu.root"),keys_, EMU_ ); 
}


void SampleProcessor::process_ttbar_ddbkg(){

  url_= TTBAR_URL;
   
  // dd bkg
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1_ddbkg.root"),keys_, TTBAR_DDBKG_);   
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2_ddbkg.root"),keys_, TTBAR_DDBKG_);   
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_ddbkg.root"),keys_, TTBAR_DDBKG_);   


}

void SampleProcessor::process_ttbar_mcbkg(){

  url_= TTBAR_URL;

  // mc bkg
  //  process(false, url_, iFolder_ + TString("ttbar_v1.root"), oFolder_+TString("out-ttbar_v1_mcbkg.root"),keys_, TTBAR_MCBKG_ );
  //  process(false, url_, iFolder_ + TString("ttbar_v2.root"), oFolder_+TString("out-ttbar_v2_mcbkg.root"),keys_, TTBAR_MCBKG_ );
  process(false, url_, iFolder_ + TString("ttbar.root"), oFolder_+TString("out-ttbar_mcbkg.root"),keys_, TTBAR_MCBKG_ );

}


void SampleProcessor::process_singletop(){

  
  //SINGLE TOP
  url_= S_URL; process(false, url_, iFolder_ + TString("stop_s.root"),     oFolder_ + TString("out-stop_s.root"),keys_);   
  //  url_= T_URL; process(false, url_, iFolder_ + TString("stop_t_v1.root"),     oFolder_ + TString("out-stop_t_v1.root"),keys_);  
  //  url_= T_URL; process(false, url_, iFolder_ + TString("stop_t_v2.root"),     oFolder_ + TString("out-stop_t_v2.root"),keys_);  
  url_= T_URL; process(false, url_, iFolder_ + TString("stop_t.root"),     oFolder_ + TString("out-stop_t.root"),keys_);  
  url_= W_URL; process(false, url_, iFolder_ + TString("stop_tW-DR.root"), oFolder_ + TString("out-stop-DR_tW.root"),keys_); 

  //ANTI SINGLE TOP
  url_= A_S_URL; process(false, url_, iFolder_ + TString("stopbar_s.root"),     oFolder_ + TString("out-stopbar_s.root"),keys_);   
  url_= A_T_URL; process(false, url_, iFolder_ + TString("stopbar_t.root"),     oFolder_ + TString("out-stopbar_t.root"),keys_);  
  url_= A_W_URL; process(false, url_, iFolder_ + TString("stopbar_tW-DR.root"), oFolder_ + TString("out-stopbar-DR_tW.root"),keys_);   


}


void SampleProcessor::process_wjets(){

 

  // WARNING W+jets with new pgid info
  url_= WJMADGRAPH_URL;
  //  process(false, url_, iFolder_ + TString("WJetsToLNu_v1.root"), oFolder_+TString("out-wjets_v1.root"),keys_);
  //  process(false, url_, iFolder_ + TString("WJetsToLNu_v2.root"), oFolder_+TString("out-wjets_v2.root"),keys_);
  process(false, url_, iFolder_ + TString("WJetsToLNu.root"), oFolder_+TString("out-wjets.root"),keys_);


  // Trigger efficiencies study //////////////////////////////////////////////////////////////////////////
  // iFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3/");
  // process(false, url_, iFolder_ + TString("trig-WJetsToLNu.root"), oFolder_+TString("out-wjets.root"),keys_);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////



}

void SampleProcessor::process_zjets_from50(){

  // MADGRAPH Z+JETS

  url_ = DY_FROM50_URL; 
  process(false, url_, iFolder_ + TString("dy_from50.root"), oFolder_+TString("out-dy_from50.root"),keys_); 

}

void SampleProcessor::process_zjets_10to50(){
  
  // MADGRAPH Z+JETS
  url_ = DY_10TO50_URL;
  process(false, url_, iFolder_ + TString("dy_10_50.root"), oFolder_+TString("out-dy_10to50.root"),keys_);    
}

void SampleProcessor::process_hh_higgs(){

  url_ = HH80_URL;  process(false, url_, iFolder_ + TString("hh-pythia-m80.root"),  oFolder_+TString("out-hh-pythia-m80.root"),keys_);  
  url_ = HH100_URL; process(false, url_, iFolder_ + TString("hh-pythia-m100.root"), oFolder_+TString("out-hh-pythia-m100.root"),keys_); 
  url_ = HH120_URL; process(false, url_, iFolder_ + TString("hh-pythia-m120.root"), oFolder_+TString("out-hh-pythia-m120.root"),keys_); 
  url_ = HH140_URL; process(false, url_, iFolder_ + TString("hh-pythia-m140.root"), oFolder_+TString("out-hh-pythia-m140.root"),keys_); 
  url_ = HH150_URL; process(false, url_, iFolder_ + TString("hh-pythia-m150.root"), oFolder_+TString("out-hh-pythia-m150.root"),keys_); 
  //  url_ = HH155_URL; process(false, url_, iFolder_ + TString("hh-pythia-m155.root"), oFolder_+TString("out-hh-pythia-m155.root"),keys_); 
  url_ = HH160_URL; process(false, url_, iFolder_ + TString("hh-pythia-m160.root"), oFolder_+TString("out-hh-pythia-m160.root"),keys_); 


}

void SampleProcessor::process_hh_higgs_bychannel(){

//  url_ = HH80_URL;  process(false, url_, iFolder_ + TString("hh-pythia-m80.root"),  oFolder_+TString("out-hh-pythia-m80_etau.root"), keys_, ETAU_);  
//  url_ = HH100_URL; process(false, url_, iFolder_ + TString("hh-pythia-m100.root"), oFolder_+TString("out-hh-pythia-m100_etau.root"),keys_, ETAU_); 
//  url_ = HH120_URL; process(false, url_, iFolder_ + TString("hh-pythia-m120.root"), oFolder_+TString("out-hh-pythia-m120_etau.root"),keys_, ETAU_); 
//  url_ = HH140_URL; process(false, url_, iFolder_ + TString("hh-pythia-m140.root"), oFolder_+TString("out-hh-pythia-m140_etau.root"),keys_, ETAU_); 
//  url_ = HH150_URL; process(false, url_, iFolder_ + TString("hh-pythia-m150.root"), oFolder_+TString("out-hh-pythia-m150_etau.root"),keys_, ETAU_); 
//  url_ = HH155_URL; process(false, url_, iFolder_ + TString("hh-pythia-m155.root"), oFolder_+TString("out-hh-pythia-m155_etau.root"),keys_, ETAU_); 
//  url_ = HH160_URL; process(false, url_, iFolder_ + TString("hh-pythia-m160.root"), oFolder_+TString("out-hh-pythia-m160_etau.root"),keys_, ETAU_); 

  
  url_ = HH80_URL;  process(false, url_, iFolder_ + TString("hh-pythia-m80.root"),  oFolder_+TString("out-hh-pythia-m80_mutau.root"), keys_, MUTAU_);  
  url_ = HH100_URL; process(false, url_, iFolder_ + TString("hh-pythia-m100.root"), oFolder_+TString("out-hh-pythia-m100_mutau.root"),keys_, MUTAU_); 
  url_ = HH120_URL; process(false, url_, iFolder_ + TString("hh-pythia-m120.root"), oFolder_+TString("out-hh-pythia-m120_mutau.root"),keys_, MUTAU_); 
  url_ = HH140_URL; process(false, url_, iFolder_ + TString("hh-pythia-m140.root"), oFolder_+TString("out-hh-pythia-m140_mutau.root"),keys_, MUTAU_); 
  url_ = HH150_URL; process(false, url_, iFolder_ + TString("hh-pythia-m150.root"), oFolder_+TString("out-hh-pythia-m150_mutau.root"),keys_, MUTAU_); 
  // url_ = HH155_URL; process(false, url_, iFolder_ + TString("hh-pythia-m155.root"), oFolder_+TString("out-hh-pythia-m155_mutau.root"),keys_, MUTAU_); 
  url_ = HH160_URL; process(false, url_, iFolder_ + TString("hh-pythia-m160.root"), oFolder_+TString("out-hh-pythia-m160_mutau.root"),keys_, MUTAU_); 

}


void SampleProcessor::process_wh_higgs(){

  url_ = WH80_URL;  process(false, url_, iFolder_ + TString("wh-pythia-m80.root"), oFolder_+TString("out-wh-pythia-m80.root"),keys_  ); 
  // For 90 must add url keys. Still, useless w/out hh, which does not exist.
  //url_ = WH90_URL;  process(false, url_, iFolder_ + TString("wh-pythia-m90.root"), oFolder_+TString("out-wh-pythia-m90.root"),keys_  ); 
  url_ = WH100_URL; process(false, url_, iFolder_ + TString("wh-pythia-m100.root"), oFolder_+TString("out-wh-pythia-m100.root"),keys_); 
  url_ = WH120_URL; process(false, url_, iFolder_ + TString("wh-pythia-m120.root"), oFolder_+TString("out-wh-pythia-m120.root"),keys_); 
  url_ = WH140_URL; process(false, url_, iFolder_ + TString("wh-pythia-m140.root"), oFolder_+TString("out-wh-pythia-m140.root"),keys_); 
  url_ = WH150_URL; process(false, url_, iFolder_ + TString("wh-pythia-m150.root"), oFolder_+TString("out-wh-pythia-m150.root"),keys_); 
  url_ = WH155_URL; process(false, url_, iFolder_ + TString("wh-pythia-m155.root"), oFolder_+TString("out-wh-pythia-m155.root"),keys_); 
  url_ = WH160_URL; process(false, url_, iFolder_ + TString("wh-pythia-m160.root"), oFolder_+TString("out-wh-pythia-m160.root"),keys_); 


//    -rw-r--r-- 1 vischia cms  1.9M Aug 28 15:22 wh-pythia-m100.root
//    -rw-r--r-- 1 vischia cms  144M Aug 28 15:22 wh-pythia-m140.root
//    -rw-r--r-- 1 vischia cms  136M Aug 28 15:22 wh-pythia-m155.root
//    -rw-r--r-- 1 vischia cms  135M Aug 28 15:23 wh-pythia-m160.root
//    -rw-r--r-- 1 vischia cms  105M Aug 28 15:23 hh-pythia-m140.root



}


void SampleProcessor::process_wh_higgs_bychannel(){

//  url_ = WH80_URL;  process(false, url_, iFolder_ + TString("wh-pythia-m80.root"),  oFolder_+TString("out-wh-pythia-m80_etau.root"),keys_,  ETAU_); 
//  url_ = WH100_URL; process(false, url_, iFolder_ + TString("wh-pythia-m100.root"), oFolder_+TString("out-wh-pythia-m100_etau.root"),keys_, ETAU_); 
//  url_ = WH120_URL; process(false, url_, iFolder_ + TString("wh-pythia-m120.root"), oFolder_+TString("out-wh-pythia-m120_etau.root"),keys_, ETAU_); 
//  url_ = WH140_URL; process(false, url_, iFolder_ + TString("wh-pythia-m140.root"), oFolder_+TString("out-wh-pythia-m140_etau.root"),keys_, ETAU_); 
//
//  url_ = WH150_URL; process(false, url_, iFolder_ + TString("wh-pythia-m150.root"), oFolder_+TString("out-wh-pythia-m150_etau.root"),keys_, ETAU_); 
//  url_ = WH155_URL; process(false, url_, iFolder_ + TString("wh-pythia-m155.root"), oFolder_+TString("out-wh-pythia-m155_etau.root"),keys_, ETAU_); 
//  url_ = WH160_URL; process(false, url_, iFolder_ + TString("wh-pythia-m160.root"), oFolder_+TString("out-wh-pythia-m160_etau.root"),keys_, ETAU_); 

  url_ = WH80_URL;  process(false, url_, iFolder_ + TString("wh-pythia-m80.root"),  oFolder_+TString("out-wh-pythia-m80_mutau.root"),keys_,  MUTAU_); 
  url_ = WH100_URL; process(false, url_, iFolder_ + TString("wh-pythia-m100.root"), oFolder_+TString("out-wh-pythia-m100_mutau.root"),keys_, MUTAU_); 
  url_ = WH120_URL; process(false, url_, iFolder_ + TString("wh-pythia-m120.root"), oFolder_+TString("out-wh-pythia-m120_mutau.root"),keys_, MUTAU_); 
  url_ = WH140_URL; process(false, url_, iFolder_ + TString("wh-pythia-m140.root"), oFolder_+TString("out-wh-pythia-m140_mutau.root"),keys_, MUTAU_); 
  url_ = WH150_URL; process(false, url_, iFolder_ + TString("wh-pythia-m150.root"), oFolder_+TString("out-wh-pythia-m150_mutau.root"),keys_, MUTAU_); 
  url_ = WH155_URL; process(false, url_, iFolder_ + TString("wh-pythia-m155.root"), oFolder_+TString("out-wh-pythia-m155_mutau.root"),keys_, MUTAU_); 
  url_ = WH160_URL; process(false, url_, iFolder_ + TString("wh-pythia-m160.root"), oFolder_+TString("out-wh-pythia-m160_mutau.root"),keys_, MUTAU_); 
  

}


void SampleProcessor::process_tbh_higgs_bychannel(){
  url_ = TBH180_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m180.root"), oFolder_+TString("out-tbh-pythia-m180_mutau.root"),keys_, MUTAU_); 
  url_ = TBH190_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m190.root"), oFolder_+TString("out-tbh-pythia-m190_mutau.root"),keys_, MUTAU_); 
  url_ = TBH200_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m200.root"), oFolder_+TString("out-tbh-pythia-m200_mutau.root"),keys_, MUTAU_); 
  url_ = TBH220_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m220.root"), oFolder_+TString("out-tbh-pythia-m220_mutau.root"),keys_, MUTAU_); 
  url_ = TBH250_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m250.root"), oFolder_+TString("out-tbh-pythia-m250_mutau.root"),keys_, MUTAU_); 
  url_ = TBH300_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m300.root"), oFolder_+TString("out-tbh-pythia-m300_mutau.root"),keys_, MUTAU_); 
  
}


void SampleProcessor::process_tbh_higgs(){
  url_ = TBH180_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m180.root"), oFolder_+TString("out-tbh-pythia-m180.root"),keys_); 
  url_ = TBH190_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m190.root"), oFolder_+TString("out-tbh-pythia-m190.root"),keys_); 
  url_ = TBH200_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m200.root"), oFolder_+TString("out-tbh-pythia-m200.root"),keys_); 
  url_ = TBH220_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m220.root"), oFolder_+TString("out-tbh-pythia-m220.root"),keys_); 
  url_ = TBH250_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m250.root"), oFolder_+TString("out-tbh-pythia-m250.root"),keys_); 
  url_ = TBH300_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("tbh-pythia-m300.root"), oFolder_+TString("out-tbh-pythia-m300.root"),keys_); 
}


void SampleProcessor::process_htb_higgs_bychannel(){
  //  url_ = HTB180_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m180.root"), oFolder_+TString("out-htb-pythia-m180_mutau.root"),keys_, MUTAU_); 
  //  url_ = HTB190_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m190.root"), oFolder_+TString("out-htb-pythia-m190_mutau.root"),keys_, MUTAU_); 
  //  url_ = HTB200_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m200.root"), oFolder_+TString("out-htb-pythia-m200_mutau.root"),keys_, MUTAU_); 
  //  url_ = HTB220_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m220.root"), oFolder_+TString("out-htb-pythia-m220_mutau.root"),keys_, MUTAU_); 
  url_ = HTB250_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m250.root"), oFolder_+TString("out-htb-pythia-m250_mutau.root"),keys_, MUTAU_); 
  //  url_ = HTB300_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m300.root"), oFolder_+TString("out-htb-pythia-m300_mutau.root"),keys_, MUTAU_); 
  
}


void SampleProcessor::process_htb_higgs(){
  //  url_ = HTB180_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m180.root"), oFolder_+TString("out-htb-pythia-m180.root"),keys_); 
  //  url_ = HTB190_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m190.root"), oFolder_+TString("out-htb-pythia-m190.root"),keys_); 
  //  url_ = HTB200_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m200.root"), oFolder_+TString("out-htb-pythia-m200.root"),keys_); 
  //  url_ = HTB220_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m220.root"), oFolder_+TString("out-htb-pythia-m220.root"),keys_); 
  url_ = HTB250_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m250.root"), oFolder_+TString("out-htb-pythia-m250.root"),keys_); 
  //  url_ = HTB300_URL; process(false, url_, TString("/lustre/ncg.ingrid.pt/cmslocal/vischia/ttbarDileptonAnalysis/topMassAnalysis/tDilbh/") + TString("htb-pythia-m300.root"), oFolder_+TString("out-htb-pythia-m300.root"),keys_); 


}


void SampleProcessor::process_qcd(){
   
  url_= QCD_EM30TO80_URL;     process(false,url_, iFolder_ + TString("qcd_EM_Pt30to80.root"),     oFolder_ + TString("out-qcd_EM_Pt30to80.root"),keys_); 
  url_= QCD_EM80TO170_URL;    process(false,url_, iFolder_ + TString("qcd_EM_Pt80to170.root"),    oFolder_ + TString("out-qcd_EM_Pt80to170.root"),keys_);    
  url_= QCD_BCTOE30TO80_URL;  process(false,url_, iFolder_ + TString("qcd_30to80_BCtoE.root"),    oFolder_ + TString("out-qcd_30to80_BCtoE.root"),keys_);  
  url_= QCD_BCTOE80TO170_URL; process(false,url_, iFolder_ + TString("qcd_80to170_BCtoE.root"),   oFolder_ + TString("out-qcd_80to170_BCtoE.root"),keys_);  
 
  url_= QCD_PHOTON30TO50_URL;   process(false,url_, iFolder_ + TString("PhotonJets_30to50.root"),   oFolder_ + TString("out-PhotonJets_30to50.root"),keys_);  
  url_= QCD_PHOTON50TO80_URL;   process(false,url_, iFolder_ + TString("PhotonJets_50to80.root"),   oFolder_ + TString("out-PhotonJets_50to80.root"),keys_);  
  url_= QCD_PHOTON80TO120_URL;  process(false,url_, iFolder_ + TString("PhotonJets_80to120.root"),  oFolder_ + TString("out-PhotonJets_80to120.root"),keys_); 
  url_= QCD_PHOTON120TO170_URL; process(false,url_, iFolder_ + TString("PhotonJets_120to170.root"), oFolder_ + TString("out-PhotonJets_120to170.root"),keys_); 

  // mu + jets specific samples
  url_= QCD_MU_URL; process(false, url_, iFolder_ + TString("qcdmu15_20toinf.root"), oFolder_ + TString("out-qcdmu15_20toinf.root"),keys_);   

 
}

void SampleProcessor::process_dibosons(){

  //iFolder_ = TString("/lustre/lip.pt/data/cmslocal/samples/CMSSW_4_2_X/mc/mTrees-v3/");  

  url_= WW_URL; process(false, url_, iFolder_ + TString("WW.root"), oFolder_+TString("out-ww.root"),keys_); 
  url_= WZ_URL; process(false, url_, iFolder_ + TString("WZ.root"), oFolder_+TString("out-wz.root"),keys_);
  url_= ZZ_URL; process(false, url_, iFolder_ + TString("ZZ.root"), oFolder_+TString("out-zz.root"),keys_);

}




void SampleProcessor::process_data_RunA(){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA.root"),oDFolder_+TString("out-Muon_RunA.root"),keys_);
  }
/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}



void SampleProcessor::process_data_RunAr(){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_06AugRecover.root"),oDFolder_+TString("out-Muon_RunA_06AugRecover.root"),keys_);
  }
/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}



void SampleProcessor::process_data_RunB(int i){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){
    switch(i){
    case 1:
      process(true,url_,iDFolder_+TString("Muon_RunB_1.root"),oDFolder_+TString("out-Muon_RunB_1.root"),keys_);
      break;
    case 2:
      process(true,url_,iDFolder_+TString("Muon_RunB_2.root"),oDFolder_+TString("out-Muon_RunB_2.root"),keys_);
      break;
    case 3:
      process(true,url_,iDFolder_+TString("Muon_RunB_3.root"),oDFolder_+TString("out-Muon_RunB_3.root"),keys_);
      break;
    case 4:
      process(true,url_,iDFolder_+TString("Muon_RunB_4.root"),oDFolder_+TString("out-Muon_RunB_4.root"),keys_);
      break;
    case 5:
      process(true,url_,iDFolder_+TString("Muon_RunB_5.root"),oDFolder_+TString("out-Muon_RunB_5.root"),keys_);
      break;
    case 6:
      process(true,url_,iDFolder_+TString("Muon_RunB_6.root"),oDFolder_+TString("out-Muon_RunB_6.root"),keys_);
      break;
    case 7:
      process(true,url_,iDFolder_+TString("Muon_RunB_7.root"),oDFolder_+TString("out-Muon_RunB_7.root"),keys_);
      break;
    case 8:
      process(true,url_,iDFolder_+TString("Muon_RunB_8.root"),oDFolder_+TString("out-Muon_RunB_8.root"),keys_);
      break;
    case 9:
      process(true,url_,iDFolder_+TString("Muon_RunB_9.root"),oDFolder_+TString("out-Muon_RunB_9.root"),keys_);
      break;
    case 10:
      process(true,url_,iDFolder_+TString("Muon_RunB_10.root"),oDFolder_+TString("out-Muon_RunB_10.root"),keys_);
      break;
    case 11:
      process(true,url_,iDFolder_+TString("Muon_RunB_11.root"),oDFolder_+TString("out-Muon_RunB_11.root"),keys_);
      break;
    case 12:
      process(true,url_,iDFolder_+TString("Muon_RunB_12.root"),oDFolder_+TString("out-Muon_RunB_12.root"),keys_);
      break;
    case 13:
      process(true,url_,iDFolder_+TString("Muon_RunB_13.root"),oDFolder_+TString("out-Muon_RunB_13.root"),keys_);
      break;
    default:
      cout << "ERROR CODE - must be in the range [1,13]" << endl;
      break;
    }
    
  }


/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}



void SampleProcessor::process_data_RunC1(){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunC1.root"),oDFolder_+TString("out-Muon_RunC1.root"),keys_);
  }
/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}



void SampleProcessor::process_data_RunC2(int i){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){ 
    //    process(true,url_,iDFolder_+TString("Muon_RunC2.root"),oDFolder_+TString("out-Muon_RunC2.root"),keys_);
    switch(i){
    case 1:
      process(true,url_,iDFolder_+TString("Muon_RunC2_1.root"),oDFolder_+TString("out-Muon_RunC2_1.root"),keys_);
      break;
    case 2:
      process(true,url_,iDFolder_+TString("Muon_RunC2_2.root"),oDFolder_+TString("out-Muon_RunC2_2.root"),keys_);
      break;
    case 3:
      process(true,url_,iDFolder_+TString("Muon_RunC2_3.root"),oDFolder_+TString("out-Muon_RunC2_3.root"),keys_);
      break;
    case 4:
      process(true,url_,iDFolder_+TString("Muon_RunC2_4.root"),oDFolder_+TString("out-Muon_RunC2_4.root"),keys_);
      break;
    case 5:
      process(true,url_,iDFolder_+TString("Muon_RunC2_5.root"),oDFolder_+TString("out-Muon_RunC2_5.root"),keys_);
      break;
    case 6:
      process(true,url_,iDFolder_+TString("Muon_RunC2_6.root"),oDFolder_+TString("out-Muon_RunC2_6.root"),keys_);
      break;
    case 7:
      process(true,url_,iDFolder_+TString("Muon_RunC2_7.root"),oDFolder_+TString("out-Muon_RunC2_7.root"),keys_);
      break;
    case 8:
      process(true,url_,iDFolder_+TString("Muon_RunC2_8.root"),oDFolder_+TString("out-Muon_RunC2_8.root"),keys_);
      break;
    case 9:
      process(true,url_,iDFolder_+TString("Muon_RunC2_9.root"),oDFolder_+TString("out-Muon_RunC2_9.root"),keys_);
      break;
    case 10:
      process(true,url_,iDFolder_+TString("Muon_RunC2_10.root"),oDFolder_+TString("out-Muon_RunC2_10.root"),keys_);
      break;
    case 11:
      process(true,url_,iDFolder_+TString("Muon_RunC2_11.root"),oDFolder_+TString("out-Muon_RunC2_11.root"),keys_);
      break;
    case 12:
      process(true,url_,iDFolder_+TString("Muon_RunC2_12.root"),oDFolder_+TString("out-Muon_RunC2_12.root"),keys_);
      break;
    case 13:
      process(true,url_,iDFolder_+TString("Muon_RunC2_13.root"),oDFolder_+TString("out-Muon_RunC2_13.root"),keys_);
      break;
    default:
      cout << "ERROR CODE - must be in the range [1,13]" << endl;
      break;
    }

  }
/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}



void SampleProcessor::process_data_RunD(int i){

  url_ = DATA_URL;

  // DON't FORGET TO SET MUTAU ON
  if( ! eChONmuChOFF_ && ! run2012_){ 
    process(true,url_,iDFolder_+TString("Muon_RunA_ReReco.root"),oDFolder_+TString("out-Muon_RunA_ReReco.root"),keys_);
  }

  if( ! eChONmuChOFF_ && run2012_){ 
    //    process(true,url_,iDFolder_+TString("Muon_RunD.root"),oDFolder_+TString("out-Muon_RunD.root"),keys_);
    switch(i){
    case 1:
      process(true,url_,iDFolder_+TString("Muon_RunD_1.root"),oDFolder_+TString("out-Muon_RunD_1.root"),keys_);
      break;
    case 2:
      process(true,url_,iDFolder_+TString("Muon_RunD_2.root"),oDFolder_+TString("out-Muon_RunD_2.root"),keys_);
      break;
    case 3:
      process(true,url_,iDFolder_+TString("Muon_RunD_3.root"),oDFolder_+TString("out-Muon_RunD_3.root"),keys_);
      break;
    case 4:
      process(true,url_,iDFolder_+TString("Muon_RunD_4.root"),oDFolder_+TString("out-Muon_RunD_4.root"),keys_);
      break;
    case 5:
      process(true,url_,iDFolder_+TString("Muon_RunD_5.root"),oDFolder_+TString("out-Muon_RunD_5.root"),keys_);
      break;
    case 6:
      process(true,url_,iDFolder_+TString("Muon_RunD_6.root"),oDFolder_+TString("out-Muon_RunD_6.root"),keys_);
      break;
    case 7:
      process(true,url_,iDFolder_+TString("Muon_RunD_7.root"),oDFolder_+TString("out-Muon_RunD_7.root"),keys_);
      break;
    case 8:
      process(true,url_,iDFolder_+TString("Muon_RunD_8.root"),oDFolder_+TString("out-Muon_RunD_8.root"),keys_);
      break;
    case 9:
      process(true,url_,iDFolder_+TString("Muon_RunD_9.root"),oDFolder_+TString("out-Muon_RunD_9.root"),keys_);
      break;
    case 10:
      process(true,url_,iDFolder_+TString("Muon_RunD_10.root"),oDFolder_+TString("out-Muon_RunD_10.root"),keys_);
      break;
    case 11:
      process(true,url_,iDFolder_+TString("Muon_RunD_11.root"),oDFolder_+TString("out-Muon_RunD_11.root"),keys_);
      break;
    case 12:
      process(true,url_,iDFolder_+TString("Muon_RunD_12.root"),oDFolder_+TString("out-Muon_RunD_12.root"),keys_);
      break;
    case 13:
      process(true,url_,iDFolder_+TString("Muon_RunD_13.root"),oDFolder_+TString("out-Muon_RunD_13.root"),keys_);
      break;
    default:
      cout << "ERROR CODE - must be in the range [1,13]" << endl;
      break;
    }

  }
/*
  else{ 
    // inclusive electron trigger
    process(true,url_,iDFolder_+TString("Electron_1.root"), oDFolder_+TString("out-Electron_1.root"),keys_);
  }
*/
/*
  else if(! eChONmuChOFF_ && ! run2012_ ) {
  electron+2jets+mht samples
  process(true,url_,iDFolder_+TString("Ele_Met_eJTrigMatch_4.root"), oDFolder_+TString("out-eJTrigMatch_4.root"),keys_);
  }
*/  
}






























void SampleProcessor::process_triggersamples_mht(){

  url_ = DATA_URL;

  iDFolder_=TString("/lustre/data3/cmslocal/samples/CMSSW_4_2_X/data/mTrees-v3/");


  process(true,url_,iDFolder_+TString("Ele_MetTrig_4.root"),oDFolder_+TString("out-Ele_MetTrig_4.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_MetTrig_5.root"),oDFolder_+TString("out-Ele_MetTrig_5.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_MetTrig_6.root"),oDFolder_+TString("out-Ele_MetTrig_6.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_MetTrig_7.root"),oDFolder_+TString("out-Ele_MetTrig_7.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_MetTrig_8.root"),oDFolder_+TString("out-Ele_MetTrig_8.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_MetTrig_9.root"),oDFolder_+TString("out-Ele_MetTrig_9.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_MetTrig_10.root"),oDFolder_+TString("out-Ele_MetTrig_10.root"),keys_);




}


void SampleProcessor::process_triggersamples_mhtjets(){

  url_ = DATA_URL;

  iDFolder_=TString("/lustre/data3/cmslocal/samples/CMSSW_4_2_X/data/mTrees-v3-B/");
/*
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_4.root"),oDFolder_+TString("out-Ele_JetMetTrig_4.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_5.root"),oDFolder_+TString("out-Ele_JetMetTrig_5.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_6.root"),oDFolder_+TString("out-Ele_JetMetTrig_6.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_7.root"),oDFolder_+TString("out-Ele_JetMetTrig_7.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_8.root"),oDFolder_+TString("out-Ele_JetMetTrig_8.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_9.root"),oDFolder_+TString("out-Ele_JetMetTrig_9.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_10.root"),oDFolder_+TString("out-Ele_JetMetTrig_10.root"),keys_);
*/

  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_4.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_4.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_5.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_5.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_6.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_6.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_7.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_7.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_8.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_8.root"),keys_); 
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_9.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_9.root"),keys_);
  process(true,url_,iDFolder_+TString("Ele_JetMetTrig_matched_10.root"),oDFolder_+TString("out-Ele_JetMetTrig_matched_10.root"),keys_);

}


