import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        parameterSets = cms.vstring('pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters'),
        processParameters = cms.vstring('HardQCD:all = on',
                                        'PromptPhoton:all = on',
                                        'PhaseSpace:pTHatMin = 15.',
                                        'PhaseSpace:pTHatMax = 9999.',
                                        'PhaseSpace:bias2Selection = on',
                                        'PhaseSpace:bias2SelectionPow = 2.0',
                                        'PhaseSpace:bias2SelectionRef = 15'),
    ),
    comEnergy = cms.double(5362.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('PYTHIA 8, Tune CP5, (unquenched) photons in NN (pt-hat > 15 GeV) at sqrt(s) = 5.36 TeV')
    )

photonFilter = cms.EDFilter("PythiaFilterMultiMother",
                            Status = cms.untracked.int32(1),
                            MinPt = cms.untracked.double(15.0),
                            ParticleID = cms.untracked.int32(22),
                            MotherIDs = cms.untracked.vint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1)
)


ProductionFilterSequence = cms.Sequence(generator*photonFilter)


# with command lines under CMSSW_13_0_18_HeavyIon:
# cmsDriver.py Configuration/GenProduction/python/QCDPhoton/Pythia8_AllQCDPhoton15_TuneCP5_5362GeV_bias_cff.py --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 130X_mcRun3_2023_realistic_HI_v18 --beamspot MatchHI --step GEN,SIM --scenario HeavyIons --geometry DB:Extended --era Run3_pp_on_PbPb --pileup HiMixGEN --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --nThreads 4 --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 50000

# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 4.407e+08 +- 1.181e+06 pb
# Filter efficiency (taking into account weights)= (17.0058) / (24055.5) = 7.069e-04 +- 1.299e-04
# Filter efficiency (event-level)= (60) / (50000) = 1.200e-03 +- 1.548e-04    [TO BE USED IN MCM]

# After filter: final cross section = 3.116e+05 +- 5.727e+04 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 3.210e-03 +- 5.900e-04

# 0.01133 sec/output event, 1.2783 kB/output event
