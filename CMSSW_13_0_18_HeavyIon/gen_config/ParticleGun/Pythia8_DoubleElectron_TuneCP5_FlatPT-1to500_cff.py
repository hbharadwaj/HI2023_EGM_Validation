import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8PtGun",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(500.0),
        MinPt = cms.double(1.0),
        ParticleID = cms.vint32(11),
        AddAntiParticle = cms.bool(True),
        MaxEta = cms.double(3.1),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-3.1),
        MinPhi = cms.double(-3.14159265359) ## in radians
    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
    psethack = cms.string('double electron pt 1.0 to 500'),
    firstRun = cms.untracked.uint32(1),
    PythiaParameters = cms.PSet(parameterSets = cms.vstring())
)

# with command lines under CMSSW_13_0_18_HeavyIon:
# cmsDriver.py Configuration/GenProduction/python/ParticleGun/Pythia8_DoubleElectron_TuneCP5_FlatPT-1to500_cff.py --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 130X_mcRun3_2023_realistic_HI_v18 --beamspot MatchHI --step GEN,SIM --scenario HeavyIons --geometry DB:Extended --era Run3_pp_on_PbPb --pileup HiMixGEN --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --nThreads 4 --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 100


# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 0.000e+00 +- 0.000e+00 pb
# Filter efficiency (taking into account weights)= (100) / (100) = 1.000e+00 +- 0.000e+00
# Filter efficiency (event-level)= (100) / (100) = 1.000e+00 +- 0.000e+00    [TO BE USED IN MCM]

# After filter: final cross section = 0.000e+00 +- 0.000e+00 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 0.000e+00 +- 0.000e+00

# 3.655 sec/output event, 714 kB/output event
