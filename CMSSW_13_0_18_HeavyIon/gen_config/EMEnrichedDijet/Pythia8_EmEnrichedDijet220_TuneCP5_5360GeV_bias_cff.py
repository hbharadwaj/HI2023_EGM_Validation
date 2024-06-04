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
                                        'PhaseSpace:pTHatMin = 220.',
                                        'PhaseSpace:pTHatMax = 9999.'),
    ),
    comEnergy = cms.double(5362.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('PYTHIA 8, Tune CP5, (unquenched) EM-enriched Dijets in NN (pt-hat > 220 GeV) at sqrt(s) = 5.36 TeV')
    )

partonFilter = cms.EDFilter("MCSingleParticleFilter",
                            ParticleID = cms.untracked.vint32(1, 2, 3, 4, 5, 6, #quarks
                                                              21, 22), #gluon, photon
                            Status = cms.untracked.vint32(2, 2, 2, 2, 2, 2,
                                                          2, 1)
                            )
neutralMesonFilter = cms.EDFilter("MCSingleParticleFilter",
                                  ParticleID = cms.untracked.vint32(111, #pi0
                                                                    221, #eta
                                                                    331, #eta'
                                                                    223), #omega
                                  Status = cms.untracked.vint32(2, #pi0
                                                                2, #eta
                                                                2, #eta'
                                                                2), #omega
                                  MinEta = cms.untracked.vdouble(-3,
                                                                 -3,
                                                                 -3,
                                                                 -3),
                                  MaxEta = cms.untracked.vdouble(3,
                                                                 3,
                                                                 3,
                                                                 3),
                                  MinPt = cms.untracked.vdouble(35,
                                                                35,
                                                                35,
                                                                35)

                                  )


ProductionFilterSequence = cms.Sequence(generator*partonFilter*neutralMesonFilter)


# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 3.038e+03 +- 4.725e+01 pb
# Filter efficiency (taking into account weights)= (695) / (1000) = 6.950e-01 +- 1.456e-02
# Filter efficiency (event-level)= (695) / (1000) = 6.950e-01 +- 1.456e-02    [TO BE USED IN MCM]

# After filter: final cross section = 2.111e+03 +- 5.509e+01 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 4.737e-01 +- 1.237e-02


# 8.8295 sec/output event, 815.7684 kB/output event