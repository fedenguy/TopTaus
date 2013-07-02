#!/bin/bash

# usage: sh combineResults.sh outputArea [hadd|clean] [AB|ABC|ABCD]

if [ "${2}" = "hadd" ]; then

    DATAPLACE="${1}/nomt-2012-V1-mc-MU-20GeV/"
    
    cd $DATAPLACE

    
    # QCD
    hadd -f $DATAPLACE/out-qcd.root $DATAPLACE/out-qcd_*root $DATAPLACE/out-qcdmu*.root $DATAPLACE/out-PhotonJets*root
    
    # zjets
    hadd -f $DATAPLACE/out-zjets.root $DATAPLACE/out-dy*root

    # single top
    hadd -f $DATAPLACE/out-singletop.root $DATAPLACE/out-stop*.root
    
    # diboson
    hadd -f $DATAPLACE/out-dibosons.root $DATAPLACE/out-ww.root $DATAPLACE/out-wz.root $DATAPLACE/out-zz.root

    # mc total
    hadd -f $DATAPLACE/out-mc.root $DATAPLACE/out-qcd.root $DATAPLACE/out-wjets.root $DATAPLACE/out-zjets.root $DATAPLACE/out-dibosons.root $DATAPLACE/out-ttbar.root $DATAPLACE/out-singletop.root &

    # data driven bkg, using mc
    hadd -f $DATAPLACE/out-ddbkg.root $DATAPLACE/out-qcd.root $DATAPLACE/out-wjets.root $DATAPLACE/out-ttbar_ddbkg.root

    # mc bkg (w.r.t. top xsec measurement: ttbar is signal for this)
    hadd -f $DATAPLACE/out-mcbkg.root $DATAPLACE/out-zjets.root $DATAPLACE/out-singletop.root $DATAPLACE/out-dibosons.root $DATAPLACE/out-ttbar_mcbkg.root

    # other cont
    hadd -f $DATAPLACE/out-ttbar_other.root  $DATAPLACE/out-ttbar_ddbkg.root $DATAPLACE/out-ttbar_mcbkg.root &

    # allmcmkg (w.r.t. charged higgs search: ttbar is includes as irreducible background)
    hadd -f $DATAPLACE/out-allmcbkg.root $DATAPLACE/out-ttbar.root $DATAPLACE/out-mcbkg.root &


    # higgs  mass points
    #hadd -f $DATAPLACE/out-higgs-pythia-m80.root $DATAPLACE/out-hh-pythia-m80.root $DATAPLACE/out-wh-pythia-m80.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m80_mutau.root $DATAPLACE/out-hh-pythia-m80_mutau.root $DATAPLACE/out-wh-pythia-m80_mutau.root
    
    #hadd -f $DATAPLACE/out-higgs-pythia-m100.root $DATAPLACE/out-hh-pythia-m100.root $DATAPLACE/out-wh-pythia-m100.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m100_mutau.root $DATAPLACE/out-hh-pythia-m100_mutau.root $DATAPLACE/out-wh-pythia-m100_mutau.root
    
    # hadd -f $DATAPLACE/out-higgs-pythia-m120.root $DATAPLACE/out-hh-pythia-m120.root $DATAPLACE/out-wh-pythia-m120.root &
    # hadd -f $DATAPLACE/out-higgs-pythia-m120_mutau.root $DATAPLACE/out-hh-pythia-m120_mutau.root $DATAPLACE/out-wh-pythia-m120_mutau.root &
    
    #hadd -f $DATAPLACE/out-higgs-pythia-m140.root $DATAPLACE/out-hh-pythia-m140.root $DATAPLACE/out-wh-pythia-m140.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m140_mutau.root $DATAPLACE/out-hh-pythia-m140_mutau.root $DATAPLACE/out-wh-pythia-m140_mutau.root
    
    #hadd -f $DATAPLACE/out-higgs-pythia-m150.root $DATAPLACE/out-hh-pythia-m150.root $DATAPLACE/out-wh-pythia-m150.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m150_mutau.root $DATAPLACE/out-hh-pythia-m150_mutau.root $DATAPLACE/out-wh-pythia-m150_mutau.root
    
    #hadd -f $DATAPLACE/out-higgs-pythia-m155.root $DATAPLACE/out-hh-pythia-m155.root $DATAPLACE/out-wh-pythia-m155.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m155_mutau.root $DATAPLACE/out-hh-pythia-m155_mutau.root $DATAPLACE/out-wh-pythia-m155_mutau.root
    
    #hadd -f $DATAPLACE/out-higgs-pythia-m160.root $DATAPLACE/out-hh-pythia-m160.root $DATAPLACE/out-wh-pythia-m160.root
    #hadd -f $DATAPLACE/out-higgs-pythia-m160_mutau.root $DATAPLACE/out-hh-pythia-m160_mutau.root $DATAPLACE/out-wh-pythia-m160_mutau.root

    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m180.root $DATAPLACE/out-tbh-pythia-m180.root $DATAPLACE/out-htb-pythia-m180.root
    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m190.root $DATAPLACE/out-tbh-pythia-m190.root $DATAPLACE/out-htb-pythia-m190.root
    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m200.root $DATAPLACE/out-tbh-pythia-m200.root $DATAPLACE/out-htb-pythia-m200.root
    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m220.root $DATAPLACE/out-tbh-pythia-m220.root $DATAPLACE/out-htb-pythia-m220.root
    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m250.root $DATAPLACE/out-tbh-pythia-m250.root $DATAPLACE/out-htb-pythia-m250.root
    #hadd -f $DATAPLACE/out-heavyHiggs-pythia-m300.root $DATAPLACE/out-tbh-pythia-m300.root $DATAPLACE/out-htb-pythia-m300.root

    cp  $DATAPLACE/out-tbh-pythia-m180.root $DATAPLACE/out-heavyHiggs-pythia-m180.root
    cp  $DATAPLACE/out-tbh-pythia-m190.root $DATAPLACE/out-heavyHiggs-pythia-m190.root
    cp  $DATAPLACE/out-tbh-pythia-m200.root $DATAPLACE/out-heavyHiggs-pythia-m200.root
    cp  $DATAPLACE/out-tbh-pythia-m220.root $DATAPLACE/out-heavyHiggs-pythia-m220.root
    hadd -f $DATAPLACE/out-heavyHiggs-pythia-m250.root $DATAPLACE/out-tbh-pythia-m250.root $DATAPLACE/out-htb-pythia-m250.root
    cp  $DATAPLACE/out-tbh-pythia-m300.root $DATAPLACE/out-heavyHiggs-pythia-m300.root 
    
    DATAPLACE="${1}/nomt-2012-V1-data-MU-20GeV/"
    
     #DATAPLACE="/lustre/data3/cmslocal/vischia/tau_dilepton/outputFiles444_3/mt-2011-V1-data-MU-20GeV/"

    cd $DATAPLACE

#    if [ "${3}" = "AB" ]; then
#    	hadd -f out-data.root out-Muon_RunA.root out-Muon_RunA_06AugRecover.root out-Muon_RunB_*.root
#    elif ["${3}" = "ABC" ]; then
#        hadd -f out-data.root out-Muon_RunA.root out-Muon_RunA_06AugRecover.root out-Muon_RunB_*.root out-Muon_RunC1.root out-Muon_RunC2_*.root
#    elif ["${3}" = "ABCD" ]; then
    hadd -f out-data.root out-Muon_RunA.root out-Muon_RunA_06AugRecover.root out-Muon_RunB_*.root out-Muon_RunC1.root out-Muon_RunC2_*.root out-Muon_RunD_*.root 
#    fi

elif [ "${2}" = "clean" ]; then
  # Clean processed files
    DATAPLACE="${1}/nomt-2012-V1-mc-MU-20GeV/"
    rm $DATAPLACE/out-*.root
    rm $DATAPLACE/out-*.txt
    DATAPLACE="${1}/nomt-2012-V1-data-MU-20GeV/"
    rm $DATAPLACE/out-*.root
    rm $DATAPLACE/out-*.txt
    #rm -rf $DATAPLACE/out-qcd.root
    #rm -rf $DATAPLACE/out-zjets.root
    #rm -rf $DATAPLACE/out-singletop.root
    #rm -rf $DATAPLACE/out-dibosons.root
    #rm -rf $DATAPLACE/out-mc.root
    #DATAPLACE="/lip-sw/cmssw/users/nalmeida/new/CMSSW_5_2_5/src/lipcms/Physics/TopTauDileptons2012/mt-2011-V1-data-MU-20GeV/"
    #rm -rf $DATAPLACE/out-data.root 

fi