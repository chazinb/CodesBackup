#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "  "
    echo "  ./merge.sh ../rootfiles/<systematic>/<analysis>"
    echo "  ./merge.sh ../minitrees/<systematic>/<analysis>"
    echo "  "
    exit -1
fi

FOLDER="$1"

pushd $FOLDER

hadd -f -k 01_Data.root       *03Feb2017*

hadd -f -k 02_WZTo3LNu.root  WZTo3LNu.root

hadd -f -k 03_VZ.root        ZZTo4L__*.root ZZTo2L2Q__part*.root ZZTo2L2Nu__part*.root WZTo2L2Q__part*.root

hadd -f -k 04_TTTo2L2Nu.root TTTo2L2Nu__part*.root

hadd -f -k 05_ST.root        ST_tW_antitop.root ST_tW_top.root

hadd -f -k 06_WW.root        WWTo2L2Nu.root GluGluWWTo2L2Nu_MCFM.root

hadd -f -k 07_ZJets.root     DYJetsToLL_M-10to50.root DYJetsToLL_M-50__part*.root

hadd -f -k 09_TTW.root       TTWJetsToLNu.root TTWJetsToLNu_ext2.root TTWJetsToQQ.root

hadd -f -k 10_TTZ.root       TTZToQQ.root TTZToLLNuNu_M-10.root 

hadd -f -k 11_HWW.root       GluGluHToWWTo2L2NuAMCNLO_M125.root VBFHToWWTo2L2Nu_M125.root GluGluHToTauTau_M125.root VBFHToTauTau_M125.root HWminusJ_HToWW_M125.root HWplusJ_HToWW_M125.root

hadd -f -k 13_VVV.root       WWW.root WWZ.root WZZ.root WWG.root 

popd
