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
                                        'PhaseSpace:pTHatMin = 120.',
                                        'PhaseSpace:pTHatMax = 9999.'),
    ),
    comEnergy = cms.double(5362.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('PYTHIA 8, Tune CP5, (unquenched) photons in NN (pt-hat > 120 GeV) at sqrt(s) = 5.36 TeV')
    )

photonFilter = cms.EDFilter("PythiaFilterMultiMother",
                            Status = cms.untracked.int32(1),
                            MinPt = cms.untracked.double(15.0),
                            ParticleID = cms.untracked.int32(22),
                            MotherIDs = cms.untracked.vint32(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1)
)


ProductionFilterSequence = cms.Sequence(generator*photonFilter)


# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 7.140e+04 +- 1.616e+02 pb
# Filter efficiency (taking into account weights)= (608) / (50000) = 1.216e-02 +- 4.901e-04
# Filter efficiency (event-level)= (608) / (50000) = 1.216e-02 +- 4.901e-04    [TO BE USED IN MCM]

# After filter: final cross section = 8.682e+02 +- 3.505e+01 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 1.152e+00 +- 4.651e-02

# 0.05642 sec/output event, 13.7407 kB/output event