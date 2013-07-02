import FWCore.ParameterSet.Config as cms

#########
# Parameter file to give as input to runTauDileptonPDFBuilderFitter and produceLandSShapes
#########

LandSShapesProducerParSet = cms.PSet(
    #    outFolder        = cms.string("shapesForDatacard/"),
    outFolder        = cms.string("shapes/"),
    #    outFolder        = cms.string("testMultiDimensionalShapes/"),
    outputFileName  = cms.string("shapes"),
#    massPointName   = cms.vstring("m180", "m190", "m200", "m220", "m250", "m300"),
    massPointName   = cms.vstring("m250"),

    
    #    baseMCDir        = cms.string("/lustre/data3/cmslocal/vischia/tau_dilepton/outputFiles444_3/mt-2011-V1-mc-MU-20GeV/"),
    #    baseDataDir      = cms.string("/lustre/data3/cmslocal/vischia/tau_dilepton/outputFiles444_3/mt-2011-V1-data-MU-20GeV/"),
    
    
    baseDir = cms.vstring("/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-data-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-data-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/",
                          "/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/"
                          ),
    
    
    # obsolete format # for fig7 #    mcBkgFileName      = cms.vstring(    "out-ttbar_mutau.root", "out-ttbar_other.root", "out-singletop.root",  "out-wjets.root"  , "out-dibosons.root", "out-zjets.root"    ),
    #    inputFileName      = cms.vstring( "out-data.root", "out-wh-pythia-", "out-hh-pythia-", "out-data_rescaled.root", "out-ttbar_mutau.root", "out-ttbar_mcbkg.root", "out-singletop.root", "out-dibosons.root" ),
    inputFileName      = cms.vstring( "out-data.root", "out-heavyHiggs-pythia-", "out-data_rescaled.root", "out-ttbar_mutau.root", "out-ttbar_mcbkg.root", "out-singletop.root", "out-dibosons.root" ),
        # is it *really* needed that "3" that was told it is necessary in the lands framework?
    #    sampleName      = cms.vstring( "data_obs", "WH3", "HH3", "tau_fake3", "tt_ltau3", "tt_ll3", "singleTop3","di_boson3", "Z_tautau3",  "Z_eemumu3" ),
    sampleName      = cms.vstring( "data_obs", "TBH", "tau_fake", "tt_ltau", "tt_ll", "singleTop","di_boson", "Z_tautau",  "Z_eemumu" ),
    
    fancySampleName      = cms.vstring("data", "H^{+}[m_{H^{#pm}} = 250 GeV]", "misidentified #tau_{h}", "t#bar{t} #rightarrow #mu#tau_{h}", "other t#bar{t}", "Single t", "Diboson","Z #rightarrow #tau#tau", "Z #rightarrow ee,#mu#mu"),
    
    # obsolete format # for fig 7 #    mcBkgSampleColor    = cms.vint32( 614, 8, 824, 809,  596, 831, 831, 835 ),
    sampleColour    = cms.vint32( 1,  1, 799,  614,    8,   824,  596, 831, 809),
    #sampleColour    = cms.vint32( 1,  616, 799, 614, 8, 824,  596, 831, 809),
    sampleFillStyle = cms.vint32( 1001, 0, 1001, 1001, 1001,  1001, 1001, 1001, 1001),
    
    doMultiDimensionalShapes = cms.bool(False),
    unsplitUncertainties = cms.bool(False),
    
    minitree = cms.vstring("m_tau_Selected",
                           "m_tau_Selected",
                           "m_tau_DataDriven",
                           "m_tau_Selected",
                           "m_tau_Selected",
                           "m_tau_Selected",
                           "m_tau_Selected"
                           ),
    
    isFromData  = cms.vint32( 1, 0, 1, 0, 0, 0, 0, 0, 0),
    isDDbkg     = cms.vint32( 0, 0, 1, 0, 0, 0, 0, 0, 0),
    isSignal    = cms.vint32( 0, 1, 0, 0, 0, 0, 0, 0, 0),
    
    # jes +/-, met +/-, jer +/-
    systComponents = cms.vstring(
    "_plus"    , "_minus"    ,
    "_uncplus" , "_uncminus" ,
    "_jerplus" , "_jerminus" 
    ),
    
    systFancyComponents = cms.vstring(
    "_jesUp"    , "_jesDown"    ,
    "_metUp" , "_metDown" ,
    "_jerUp" , "_jerDown" 
    ),
    
    #osCutEff               = cms.double(0.69*220.92/242.915), #    osCutEff               = cms.double(0.577749097), #223.82/387.4),
#    osCutEff               = cms.double(2054.99/4232.500009), #FIXME: subtract residual shape 0.699), # // FIXME: fixme 220.92/242.915), #    osCutEff               = cms.double(0.577749097), #223.82/387.4),
    osCutEff = cms.double( 2810.89/5736.58),
    #fitType     = cms.vstring( HIGGS2BKG, HIGGS3BKG, SM2BKG, SM3BKG)
    #    fitType     = cms.vint32( 0),
    #, 2, 3, 0),

#    mins        = cms.vdouble(  0.0),
#    maxs        = cms.vdouble(  1.001),
#    bins        = cms.vint32(     10),




    vars        = cms.vstring( "rc_t"    ),
    mins        = cms.vdouble(   -0.099  ),
    maxs        = cms.vdouble(    1.001  ),
    bins        = cms.vint32(      11     ),
    hmin        = cms.vdouble(      0    ),
    hmax        = cms.vdouble(    200    ),
    unbinned    = cms.vint32(       0    ),
    smoothOrder = cms.vint32(       3    ),   


# vars        = cms.vstring(    "pt_t"  ),
# mins        = cms.vdouble(         0  ),
# maxs        = cms.vdouble(       500  ),
# bins        = cms.vint32(         20  ),
# hmin        = cms.vdouble(         0  ),
# hmax        = cms.vdouble(       200  ),
# unbinned    = cms.vint32(          0  ),
# smoothOrder = cms.vint32(          3  ),



#    vars        = cms.vstring( "rc_t"  ,"Dphi_tau_met"  ),
#    mins        = cms.vdouble(   -0.099,     0           ),
#    maxs        = cms.vdouble(    1.001,   4.8          ),
#    bins        = cms.vint32(      6  ,    5          ),
#    hmin        = cms.vdouble(      0  ,     0          ),
#    hmax        = cms.vdouble(    200  ,   200          ),
#    unbinned    = cms.vint32(       0  ,     0          ),
#    smoothOrder = cms.vint32(       3  ,     3          ),   


#    vars        = cms.vstring( "rc_t"    , "pt_t"  ),
#    mins        = cms.vdouble(   -0.099  ,      0   ),
#    maxs        = cms.vdouble(    1.001  ,    200  ),
#    bins        = cms.vint32(      6    ,     5  ),
#    hmin        = cms.vdouble(      0    ,      0  ),
#    hmax        = cms.vdouble(    200    ,    100  ),
#    unbinned    = cms.vint32(       0    ,      0  ),
#    smoothOrder = cms.vint32(       3    ,      5  ),   


# for multiD    vars        = cms.vstring("rc_t"  , "Dphi_tau_met"   ),
# for multiD    mins        = cms.vdouble(  -0.099,      0           ),
# for multiD    maxs        = cms.vdouble(   1.001,    4.8           ),
# for multiD    bins        = cms.vint32(     11  ,     10           ),
# for multiD    hmin        = cms.vdouble(     0  ,      0           ),
# for multiD    hmax        = cms.vdouble(   200  ,    200           ),
# for multiD    unbinned    = cms.vint32(      0  ,      0           ),
# for multiD    smoothOrder = cms.vint32(      3  ,      3           ),   

#    "pt_t", ,      "pt_t" 
#         0, ,          -20
#       200, ,         200 
#        20, ,          11 
#         0, ,           0 
#       100, ,         100 
#         0, ,           0 
#         5, ,           5 
    
    
    #    vars        = cms.vstring("rc_t",   "pt_l", "pt_met", "multiplicity_j", "btagmultiplicity_j",    ),
    #    mins        = cms.vdouble(     0,        0,      0,        2,                0,                  ),
    #    maxs        = cms.vdouble(   1.152,    200,    300,      7,                4,                    ),
    #    bins        = cms.vint32(     12,       20,     30,       5,                4,                   ),
    #    hmin        = cms.vdouble(     0,        0,      0,        0,                0,                  ),
    #    hmax        = cms.vdouble(   200,      100,    100,      400,              600,                  ),
    #    unbinned    = cms.vint32(      0,        0,      0,        0,                0,                  ),
    #    smoothOrder = cms.vint32(      3,        5,      5,        0,                0,                  )
    
    #    vars        = cms.vstring( "pt_l", "pt_met", "multiplicity_j", "btagmultiplicity_j","rc_t",   "Dphi_tau_met"    ),
    #    mins        = cms.vdouble(      0,      0,        2,                0,                   0,        0            ),
    #    maxs        = cms.vdouble(    200,    300,      7,                4,                   1.152,    4.8            ),
    #    bins        = cms.vint32(      20,     30,       5,                4,                   12,       10            ),
    #    hmin        = cms.vdouble(      0,      0,        0,                0,                   0,        0            ),
    #    hmax        = cms.vdouble(    100,    100,      400,              600,                 200,      200            ),
    #    unbinned    = cms.vint32(       0,      0,        0,                0,                   0,        0            ),
    #    smoothOrder = cms.vint32(       5,      5,        0,                0,                   3,        3            )
    
    )

TauDileptonPDFBuilderFitterParSet = cms.PSet(
    outFolder        = cms.string("likelihoodFit/"),
    resultsFileName  = cms.string("likelihoodResults"),
    
    baseMCDir        = cms.string("/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-mc-MU-20GeV/"),
    baseDataDir      = cms.string("/lustre/ncg.ingrid.pt/cmslocal/samples/CMSSW_5_3_7_patch4/chiggs/nomt-2012-V1-data-MU-20GeV/"),

    signalFileNameWH   = cms.string("out-wh-pythia-m120.root"),
    signalFileNameHH   = cms.string("out-hh-pythia-m120.root"),
    
    dataFileName       = cms.string("out-data.root"),
    ddBkgFileName      = cms.string("out-data_rescaled.root"),
    ttbarmcBkgFileName = cms.string("out-ttbar_mutau.root"),
    mcBkgFileName      = cms.string("out-mcbkg.root"),
    
    minitreeSelected   = cms.string("m_tau_Selected"),
    minitreeDataDriven = cms.string("m_tau_DataDriven"),

    #fitType     = cms.vstring( HIGGS2BKG = 2, HIGGS3BKG = 3, SM2BKG = 0, SM3BKG = 1)
    fitType     = cms.vint32( 1),
    #, 2, 3, 0),
    useOS       = cms.bool(False),

    #    osCutEff               = cms.double(0.69 * 220.92/242.915),
    osCutEff               = cms.double(0.69 * 220.92/242.915),
    cHiggsBR               = cms.double(0.01),
    
    vars        = cms.vstring( "pt_l",  "multiplicity_j", "btagmultiplicity_j",  "Dphi_tau_met"    ),
    mins        = cms.vdouble(      0,       2,                0,                     0            ),
    maxs        = cms.vdouble(    200,     7,                4,                     4.8            ),
    bins        = cms.vint32(      20,      5,                4,                     10            ),
    hmin        = cms.vdouble(      0,       0,                0,                     0            ),
    hmax        = cms.vdouble(    100,     400,              600,                   200            ),
    unbinned    = cms.vint32(       0,       0,                0,                     0            ),
    smoothOrder = cms.vint32(       5,       0,                0,                     2            )

#"pt_met",
#     0,  
#   300,  
#    30,  
#     0,  
#   100,  
#     0,  
#     5,  

#"rc_t",
#     0,
#   1.2,
#    12,
#     0,
#   200,
#     0,
#     3,

    )


