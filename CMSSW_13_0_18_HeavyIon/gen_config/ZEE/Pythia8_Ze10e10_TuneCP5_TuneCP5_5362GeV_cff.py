import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.0),
                         #crossSection = cms.untracked.double(425.6),
                         comEnergy = cms.double(5362.0),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'WeakSingleBoson:ffbar2gmZ = on',
            '23:onMode = off',
            '23:onIfAny = 11',
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters',
                                    )
        )
                         )

eegenfilter = cms.EDFilter("MCParticlePairFilter",
    Status         = cms.untracked.vint32(1, 1),
    MinPt          = cms.untracked.vdouble(10, 10),
    MaxEta         = cms.untracked.vdouble(2.5, 2.5),
    MinEta         = cms.untracked.vdouble(-2.5, -2.5),
    MinInvMass     = cms.untracked.double(60.0),
    MaxInvMass     = cms.untracked.double(120.0),
    ParticleCharge = cms.untracked.int32(-1),
    ParticleID1    = cms.untracked.vint32(11),
    ParticleID2    = cms.untracked.vint32(11)
)

ProductionFilterSequence = cms.Sequence(generator*eegenfilter)

# with command lines under CMSSW_13_0_18_HeavyIon:
# cmsDriver.py Configuration/GenProduction/python/Pythia8_Ze10e10_TuneCP5_TuneCP5_5362GeV_cff.py --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 130X_mcRun3_2023_realistic_HI_v18 --beamspot MatchHI --step GEN,SIM --scenario HeavyIons --geometry DB:Extended --era Run3_pp_on_PbPb --pileup HiMixGEN --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --nThreads 4 --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 10000

# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 3.349e+03 +- 1.771e+01 pb
# Filter efficiency (taking into account weights)= (1101) / (10000) = 1.101e-01 +- 3.130e-03
# Filter efficiency (event-level)= (1101) / (10000) = 1.101e-01 +- 3.130e-03    [TO BE USED IN MCM]

# After filter: final cross section = 3.688e+02 +- 1.066e+01 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 2.712e+00 +- 7.847e-02


# 0.1602 sec/output event, 96 kB/output event
