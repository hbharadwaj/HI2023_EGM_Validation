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
                                        'PhaseSpace:pTHatMin = 50.',
                                        'PhaseSpace:pTHatMax = 9999.'),
    ),
    comEnergy = cms.double(5362.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('PYTHIA 8, Tune CP5, (unquenched) EM-enriched Dijets in NN (pt-hat > 50 GeV) at sqrt(s) = 5.36 TeV')
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
# Before Filter: total cross section = 3.938e+06 +- 2.082e+04 pb
# Filter efficiency (taking into account weights)= (394) / (10000) = 3.940e-02 +- 1.945e-03
# Filter efficiency (event-level)= (394) / (10000) = 3.940e-02 +- 1.945e-03    [TO BE USED IN MCM]

# After filter: final cross section = 1.552e+05 +- 7.706e+03 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 6.444e-03 +- 3.201e-04


# 0.94869 sec/output event, 40.9621 kB/output event