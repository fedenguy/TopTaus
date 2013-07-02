-----------------------------------
Repository for tauDilepton analysis
-----------------------------------
* This is the manual. Read the fucking manual +BEFORE+ running the code.

* For help: pietro.vischia atPlzNoSpam gmail.com / pietro.vischia atPlzNoSpam cern.ch

* Will work for ttbar cross section measurement as well as charged higgs search/tau polarization

* This started as a porting to a compiled binary approach like LIP/Top, in order
to remove levels of messiness in the original code to be compiled every time together with
its libraries - whoah running time is already 1/2 than before, LoL

* MOVED TO GITHUB.COM ON APRIL 30th, 2013 because of CERN CVS service shutdown

------------
Installation
------------

# 537
setenv SCRAM_ARCH slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_7_patch4
cd CMSSW_5_3_7_patch4/src/
cmsenv
cvs co -p UserCode/LIP/TopTaus/TAGS_2012.txt | sh
scram b -j8

# Production version for Higgs Combination Tool
setenv SCRAM_ARCH slc5_amd64_gcc472 
cmsrel CMSSW_6_1_1 ### must be >= 6.1.1, as older versions have bugs
cd CMSSW_6_1_1/src 
cmsenv
addpkg HiggsAnalysis/CombinedLimit V03-01-08
scramv1 b

# 444
setenv SCRAM_ARCH slc5_amd64_gcc434
scramv1 project CMSSW CMSSW_4_4_4
cd CMSSW_4_4_4/src
cvs co -p UserCode/LIP/TopTaus/TAGS_2011.txt | sh

---------------------------------
Instructions for running analysis
---------------------------------

# Local pattuple test run (local data file at LIP)
# obsolete, must update: cmsRun LIP/TopTaus/test/createDataPattuple_cfg.py /lustre/data3/cmslocal/samples/CMSSW_5_2_5/test/Run2012B_SingleMu_AOD_PromptReco-v1_000_193_998_0C7DCC80-4E9D-E111-B22A-001D09F25267.root pattuple.root inclusive_mu


Finding HLT trigger on a per-run-range basis, together with the prescale
------------------------------------------------------------------------
ssh lxplus.cern.ch # needed for access to database
cvs co HLTrigger/Tools/python
scram b
cd HLTrigger/Tools/python
source $CMSSW_BASE/LIP/TopTaus/scripts/getTriggers.sh
This will fetch the list of triggers, and depending on what you comment will print or not the prescales table. (see TODO list)


Multiple likelihood fit for improving tau fakes uncertainty
-----------------------------------------------------------
runTauDileptonPDFBuilderFitter tauDileptonAnalysisParSets_cfg.py

Relevant parameters: straightforward cfg file. Must write doc here, though

Produce shapes:
---------------
produceLandSShapes tauDileptonAnalysisParSets_cfg.py --produceOnly (true|false [false])

Relevant parameters: straightforward cfg file.
Command options: 
- produceOnly: (necessary because in plots HH and WH are rescaled, in rootfile not)
               in TBH won't be necessary,  so fuck it.
	       -True: just produces rootfile with all the shapes, no plot.
	       -False: [default] just produces plots, no rootfile.

Run physicsAnalysis - this corresponds to the obsolete PhysicsAnalyzer code
------------------------------------------------------------------
- interactively:
  physicsAnalysis test/physicsAnalysisParSets_cfg.py sample
  (where "sample" is ttbar, wjets, etc. Strings can be found in bin/physicsAnalysis.cc)
- on the LIP batch:
  cd scripts/lip-batch/
  sh combineResults.sh path/to/dir clean --> this cleans the output directory in which you want to store the outputs	
  EDIT submit-jobs.sh in order to change the script path to your local installation.
  EDIT the *output* paths in test/physicsAnalysisParSets_cfg.py . The input are fixed (53X production)
  sh submit-jobs.sh
  qstat -u username --> this checks whether your jobs are running and if they did already finish
  sh combineResults.sh path/to/dir hadd [AB|ABC|ABCD] --> this hadds the relevant output files
  RUN physicsAnalysis test/physicsAnalysisParSets_cfg.py doPlots   ---> this produces plots
  RUN physicsAnalysis test/physicsAnalysisParSets_cfg.py doTables  ---> this produces tables and datacards (if "higgs" is turned on)
- have fun

Produce plots
-------------
cd LIP/TopTaus  # this is necessary, because at the moment bin/physicsAnalysis.cc is configured
                # in such a way as to fetch plots from "data/plotter/". Must change to configurable paths via cfg file
mkdir plots
physicsAnalysis test/physicsAnalysisParSets_cfg.py doPlots

Produce tables
--------------
physicsAnalysis test/physicsAnalysisParSets_cfg.py doTables --doDatacards [true|false]
		--doDatacards: true means datacards are produced, and tbh yields are rescaled to the ttbar production cross-section
			       false means datacards are not produced, and tbh yields are rescaled to the production cross section of the tbh samples (1.1pb)

Produce fakerate
----------------
for a list of all the options:
    doTauFakesStudy test/physicsAnalysisParSets_cfg.py --help 
for running everything:
    doTauFakesStudy test/physicsAnalysisParSets_cfg.py --do all

Datacards
---------
Install repository (CMS members only):
svn co svn+ssh://vischia@svn.cern.ch/reps/lipcms/Physics/datacards lipcms/Physics/datacards
chiggs datacards are in lipcms/Physics/datacards/chiggs/
xsec datacards are in lipcms/Physics/datacards/xsec/
See lipcms/Physics/datacards/Readme.txt for additional info


CHANGELOG for major updates:
2013-05-08: Added H+->tb yields to summary table, and fixed in a consistent way its xsec*BR
	    H+->taunu and H+->tb are now stacked in plots
	    Added H+->tb decay to datacards for mH+ = 250GeV/c^2
	    Committed revision 718 for chiggs datacards
2013-05-07: Reworked table builder for not having to run twice (tables, datacards)
	    Almost implemented decoupling of event analysis and application of a selection
	    Switched to text names for nuisance parameters (since february combine likes no numbers - for event counting this is not a problem, but for consistency I switched)
	    Removed unused nuisance parameters (full-hadronic ones)
2013-05-06: Fixed jer smearing/met propagation (condensed in single loop)
	    Fixed normalization for datacards and shapes
	    committed revision 716 for chiggs datacards - supersedes all previous revisions
2013-05-03: Enabled TDRStyle for kNN monitoring plots
	    Not show anymore ratio in R plot
2013-05-02: Enabled analsis of the new H+->tb (noH+->\tau\nu) samples
2013-05-01: Fixed bug in jer smearing propagation to met
2013-04-30: REPOSITORY MOVED TO https://github.com/vischia/TopTaus.git. No more updates will be committed to CERN CVS. After a week of testing at github, LIP/TopTaus/ files will be zeroed with the exception of the Readme.txt
	    plot fixes
	    jer smearing to tf
2013-04-29: committed revision 715 for chiggs datacards
	    cleanup and fix of shapes producer
	    jer and propagation to met fix
	    updated fakes values in summary table
2013-04-27: new kinematic variables added
	    plots ranges for eta tuned
	    some prefixes for optimization
	    metscaling testing
	    cleanup in plots
	    cleanup in todolist
2013-04-24: 2.41AM: added stat+syst bands in ratio plots
	    wJetsAnalysis was not using jerc (bugfix)
	    fixed btag bug for wJetsAnalysis
	    updated btag scale factors
2013-04-22: fixed bug in variables for OS shapes (it was fixed in 2011 and was never committed to the 2012 branch
	    added automatic calculation and table of fakes
2013-04-19: added index.html automatically put inside plots directory for easy and fancy propagation of the plots
	    configured LandsShapesProducer class for heavyChiggs
	    configured TauDileptonTableBuilder for producing datacards with shapes
	    committed revision 714 for chiggs datacards
2013-04-18: added mjj vs mjjb study (discrimination sig/ttbarbkg)
	    added residual mc shape subtraction from dd shape K-S test
	    added combine installation instructions and datacard automatic maker	
	    committed revision 702 for chiggs datacards
	    added doTauFakesStudy (complete tau fakes study)
2013-04-15: added chiggs BRs to TauDileptonTableBuilder, almost finished tau fakes software
2013-04-08: configured combineResults.sh script for different data hadding schemes 
2013-03-27: fixed jet smearing factors
	    fixed index of muon trigger eff muon, just for completeness
	    (it was already correct since the collection is already pt-ordered)
2013-03-22: muon trigger efficiencies scale factors added
	    pileup files splitted per run (data/pileup)
	    lumi per run added to CommonDefinitions.cc
2013-03-21: lepton-jet masses plots added
2013-03-20: fixed pileup reweighting
	    added table with the codes for antiLepton discriminators for taus
	    chosen for chiggs: antiMuonTight for highPt, againstElectronMVA3Medium
	    added class and executable for training distributions for fakes
2013-03-19: fixed tables for SM ttbar
	    added table for heavy charged higgs samples
2013-03-18: new data ntuples with improved splitting.
	    Modified combineResults.sh in order to support directory cleaning
	    New pileup files, and chosen value: 70300
	    data/plotter/samples.xml now has just final name path. Full path is read from cfg file (outputArea parameter)
2013-03-05: ported to NCG
	    improved lip-batch scripting system: now the user path is set once only in submit-jobs.sh
	    I/O and PU filenames moved to cfg files
2013-02-24: fixes: plotter perfectly working
2013-02-22: plotter and table added and integrated into physicsAnalysis.
2013-02-13: updated with 2012 cross sections - must check sample-specific ones
	    jet pt, btag, mlj added for leadingjet/nlj/nnlj
2013-02-11: updated JEC uncertainty sources
	    code is running on samples
	    lip-batch submission implemented
2013-02-08: ended base porting of the code. Now must move stuff to config file.
2013-02-07: added scripts/getTriggers.sh for easy trigger and prescales fetching
	    (will be crucial in particular for etau channel, later on)
2013-01-23: Likelihood fitter for tau fakes has RooDataset which will take into account correlations
            when performing the combined fit. Must just switch likelihoods to use it and to
	    product of likelihoods
2012-12-17: 2011 tagged (V12-12-17), start 2012 heavy developing and code porting
2012-09-23: memory problem solved. Now both lands and mlf are usable with multiple variables
2012-09-21: shapes producer for lands and multiple likelihood fit for fake rates are ready to use.

FIXED TAGS:
- V13-03-07: full working heavy charged higgs tag
- V12-12-21 / V12-12-17: 2011 code.

CHANGELOG for datacards:
rev 718 (2013-05-08): cHiggs: H+->tb decay added to mH+=250GeV/c^2 datacard
rev 716 (2013-05-06): cHiggs: fixed yields and shapes. This version supersedes all the previous.
rev 715 (2012-04-29): cHiggs: moved to use pt_tau shapes. Fixed datacards with values after met correction.
rev 714 (2012-04-19): fixed first shapes files and datacards with support for shapes (R_tau) for heavy charged higgs		
rev 713 (2012-04-19): first shapes files and datacards with support for shapes (R_tau) for heavy charged higgs	
rev 702 (2012-04-18): first counting-only chiggs datacards for heavy charged higgs


TODO:
- fix uncertainty band plotting in case the plot is normalized to 1 (Rtau)
- update tau ID uncertainties
- update uncertainties in HistogramPlotter
- move to txt file the cout of the fakes computation in TauFakesHelper.cc (latex table and everything)
- reweight MC distros con q/g fraction evaluated from MC
- From my 2013-04-22 presentation at TauPOG: recalculate fakes feeding the wmu sample the qcd mva and vice versa, try training with more kinematical variables
- move all 234.s in the code (TauDileptonTableBuilder in particular) to the value from SampleProcessor (or move the xsecs to CommonDefinitions, better) 
- vstrings in physicsAnalysisParSets_cfg.py for pileup syst
- move path of plots output directory to command line or cfg file.
- split JEC uncertainty sources into the 16 components.
- scripts/getTriggers.hs: implement running through sh getTriggers.sh PARS, in order to allow to choose
  at runtime whether to print the prescales table or just the triggers

- Make the unbinned fit work, caralho

- Apply single RooDataset per sample to LandSShapesProducer too

- Rework of LandSShapesProducer and PDFBuilderFitter
  --> New structure: -> single class for building datasets and models
                     -> one inheriting class for shapes production, storing and plotting (LandSShapesProducer)
		     -> one inheriting class for fitting, saving and plotting (PDFBuilderFitter)
		     -> perhaps one single executable with switch via ParSet in config file

- Pattuple producer 
  --> Acquire and improve version from 2013-01-23 - must upload from lipcms svn - in order to use
      multicrab in the future

- Ntuples producer --> Convert plugin

- Analysis code:
                 --> move hardcoded parameters to configuration file via edm::ParameterSet (minimize the need for recompiling)
		 --> import batch submission and modify process in order to run on a given subset of events.
		     (add s.th like void processEvents(uint firstEv, uint lastEv, string dataset){ if(string...) CutflowAnalyzer::process_ttbar(.. process(firstEv, lastEv) ..)  }
		     Doing that needs to redefine the total number of events in the event reweighting according to xsec and N_mc_evts
- Tau energy scale for main analysis code. 

- Split fake contribution into e and mu
  --> Evan/Monica 2013-01-23: "there are different scale factors for e and mu fakes to be applied,
      so you better split the two components in the cutflow yields table"

- Fix crash after natural end when running interactively (low priority)      
