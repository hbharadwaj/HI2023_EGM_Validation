# HI2023_EGM_Validation

Branch containing the gen configs for Official 2023 PbPb MC requests submitted here: https://twiki.cern.ch/twiki/bin/view/CMS/MCFor2023PbPb5p36TeV
Official Requests were submitted for GEN in CMSSW_13_0_18_HeavyIon.

The instructions to fill the tables were taken from here: https://paper.dropbox.com/doc/Instruction-to-MC-Production--CQbLdImV5TvLNs989XB29p55Ag-wK2N6gYsHbMyWiYIZuaof

To get these values you need to run the output configuration file of cmsDriver command:
```
    cmsRun -e -j log.xml config.py 
```
After all the events are processed, you will see the printout like this:

 ```
   ------------------------------------
    GenXsecAnalyzer:
    ------------------------------------
    Before Filter: total cross section = 1.551e+09 +- 8.971e+07 pb
    Filter efficiency (taking into account weights)= (16) / (100) = 1.600e-01 +- 3.666e-02
    Filter efficiency (event-level)= (16) / (100) = 1.600e-01 +- 3.666e-02    [TO BE USED IN MCM]
    
    After filter: final cross section = 2.482e+08 +- 5.866e+07 pb
    After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
    After filter: final equivalent lumi for 1M events (1/fb) = 4.029e-06 +- 9.520e-07
```
- The filter efficiency is the value of `Filter efficiency (event-level)` . The number of events for the test is normally determined to get the efficiency relative error < 10%.
- The cross-section is the value of `After filter: final cross section`
- The size/event is the output .root file size / number of output events
- The time/event  is 1/EventThroughput, and the value of “EventThroughput” can be found in `log.xml` → Note: this is the time per output event instead of the time per input event.

---

### Private MC studies

Private MC was produced with GEN in CMSSW_13_0_16_HeavyIon. Instructions and samples listed here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiEgamma2024

CMSSW_13_0_16_HeavyIon was not available in slc7 machines, only in el8 for some reason. GEN-SIM was made in lxplus8 with CMSSW_13_0_16_HeavyIon

Place gen fragment file under the path 'Configuration/genproduction/python/' ahead of these steps. Gen Fragments available for QCDPhoton30, Ze10e10 and EmEnrichedDijet30.

STEP 1: GEN-SIM
```
cmsrel CMSSW_13_0_16_HeavyIon
cd CMSSW_13_0_16_HeavyIon/src
mkdir -p Configuration/genproduction/python/
cmsenv

# Build
scram b -j 8

cmsDriver.py NameOfFragment --mc --eventcontent RAWSIM --pileup HiMixGEN --datatier GEN-SIM --conditions 130X_mcRun3_2023_realistic_HI_v18 --beamspot MatchHI --step GEN,SIM --scenario HeavyIons --geometry DB:Extended --era Run3_pp_on_PbPb --fileout file:step1.root --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --nThreads 8 --no_exec -n 10
```

STEP 2: RAW-DIGI

Further steps were made in slc7 with CMSSW_13_2_10
```
cmsrel CMSSW_13_2_10
cd CMSSW_13_2_10/src
cmsenv

cmsDriver.py step2 --mc --eventcontent RAWSIM --pileup HiMix --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --conditions 132X_mcRun3_2023_realistic_HI_v9 --step DIGI:pdigi_hi_nogen,L1,DIGI2RAW,HLT:HIon --geometry DB:Extended --era Run3_pp_on_PbPb_2023 --filein file:step1.root --fileout file:step2.root --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --nThreads 8 --no_exec -n 10
```

STEP 3: AOD/RECO

Could not directly make MiniAOD from DIGI
```
cmsDriver.py step3 --mc --eventcontent AODSIM --datatier AODSIM --conditions 132X_mcRun3_2023_realistic_HI_v9 --customise_commands "process.hltSiStripRawToDigi.ProductLabel='rawDataCollector';process.hltScalersRawToDigi.scalersInputTag='rawD ataCollector'" --step REPACK:DigiToApproxClusterRaw,RAW2DIGI,L1Reco,RECO --era Run3_pp_on_PbPb_approxSiStripClusters_2023 --filein file:step2.root --fileout file:step3.root --nThreads 8 --no_exec -n 10
```

STEP 3: MiniAOD

```
cmsDriver.py step4 --mc --eventcontent MINIAODSIM --datatier MINIAODSIM --conditions 132X_mcRun3_2023_realistic_HI_v9 --step PAT --geometry DB:Extended --era Run3_pp_on_PbPb_2023 --filein file:step3.root --fileout file:step4.root --nThreads 8 --no_exec -n 10
```
