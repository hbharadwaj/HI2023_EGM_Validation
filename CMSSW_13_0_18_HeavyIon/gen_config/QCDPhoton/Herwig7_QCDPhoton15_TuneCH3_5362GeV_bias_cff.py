import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Herwig7Settings.Herwig7CH3TuneSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7PSWeightsSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7StableParticlesForDetector_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7_7p1SettingsFor7p2_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7LHECommonSettings_cfi import *

#refer to: https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=QCD_PtGT15_TuneCH3_5p02TeV_herwig7

generator = cms.EDFilter("Herwig7GeneratorFilter",
#generator = cms.EDFilter("ThePEGGeneratorFilter",
                         herwig7CH3SettingsBlock,
                         herwig7PSWeightsSettingsBlock,
                         herwig7StableParticlesForDetectorBlock,                         
                         herwig7p1SettingsFor7p2Block,
                         configFiles = cms.vstring(),
                         crossSection = cms.untracked.double(-1),
                         dataLocation = cms.string('${HERWIGPATH:-6}'),
                         eventHandlers = cms.string('/Herwig/EventHandlers'),
                         filterEfficiency = cms.untracked.double(1.0),
                         generatorModule = cms.string('/Herwig/Generators/EventGenerator'),
                         #configFiles = cms.vstring('LHC-Matchbox.in'),
                         hw_user_settings = cms.vstring(
                             'read snippets/PPCollider.in', 
                             #'cd /Herwig/EventHandlers',
                             #'set EventHandler:LuminosityFunction:Energy 13000*GeV',
                             #'set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 5020.0',
                             #'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.0*GeV',
                             'cd /Herwig/Generators',
                             'set EventGenerator:EventHandler:LuminosityFunction:Energy 5362.0',
                             'set /Herwig/Shower/ShowerHandler:IntrinsicPtGaussian 2.0*GeV',
                             'cd /',
                             #'read snippets/PPCollider.in',
                             'mkdir /Herwig/Weights',
                             'cd /Herwig/Weights',
                             'create ThePEG::ReweightMinPT reweightMinPT ReweightMinPT.so',
                             'cd /Herwig/MatrixElements/',
                             'insert SubProcess:MatrixElements[0] MEGammaJet',
                             'insert SubProcess:Preweights[0] /Herwig/Weights/reweightMinPT',
                             'cd /',
                            'set /Herwig/Weights/reweightMinPT:Power 4.5',
                             'set /Herwig/Weights/reweightMinPT:Scale 15*GeV',
                            'set /Herwig/Cuts/JetKtCut:MinKT 15*GeV',
                             'set /Herwig/Cuts/JetKtCut:MaxKT 6000*GeV',
                             'set /Herwig/Cuts/Cuts:MHatMin 0.0*GeV',
                             'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',
                             #'set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 5020.0',
                             #'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.0*GeV',
                         ),
                         parameterSets = cms.vstring(
                             'herwig7CH3PDF',
                             'herwig7CH3AlphaS', 
                             'herwig7CH3MPISettings',
                             'hw_PSWeights_settings',
                             'herwig7StableParticlesForDetector',
                             
                             'hw_user_settings',
                             'hw_7p1SettingsFor7p2',
                         ),
                         repository = cms.string('${HERWIGPATH}/HerwigDefaults.rpo'),
                         run = cms.string('InterfaceMatchboxTest'),
                         runModeList = cms.untracked.string('read,run'),
)

#from GeneratorInterface.Core.ExternalGeneratorFilter import ExternalGeneratorFilter
#generator = ExternalGeneratorFilter(_generator)
ProductionFilterSequence = cms.Sequence(generator)

# ------------------------------------
# GenXsecAnalyzer:
# ------------------------------------
# Before Filter: total cross section = 3.253e+04 +- 5.683e+03 pb
# Filter efficiency (taking into account weights)= (3.08419) / (3.08419) = 1.000e+00 +- 0.000e+00
# Filter efficiency (event-level)= (200) / (200) = 1.000e+00 +- 0.000e+00    [TO BE USED IN MCM]
                                                                                                                                                            
# After filter: final cross section = 3.253e+04 +- 5.683e+03 pb
# After filter: final fraction of events with negative weights = 0.000e+00 +- 0.000e+00
# After filter: final equivalent lumi for 1M events (1/fb) = 3.074e-02 +- 5.371e-03

# 1.23625 sec/output event, 1009.86 kB/output event
