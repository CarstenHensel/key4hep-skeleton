from Gaudi.Configuration import *
import os
import sys


from Configurables import k4DataSvc
from Configurables import MarlinProcessorWrapper, EDM4hep2LcioTool, Lcio2EDM4hepTool
from Configurables import EventDataSvc
from k4FWCore import IOSvc
from k4MarlinWrapper.io_helpers import IOHandlerHelper
from k4FWCore.parseArgs import parser

inputFiles = [] #TODO: add input files 

# setting up the input
alg_list = []
evt_svc = EventDataSvc("EventDataSvc")
evt_svc.OutputLevel = INFO
svc_list = [evt_svc]
io_svc = IOSvc()

io_handler = IOHandlerHelper(alg_list, io_svc)
io_handler.add_reader(reco_args.inputFiles)


from Configurables import SelectEvents
myalg = SelectEvents()
myalg.cross_section = 199.21827
myalg.n_events_generated = 43200
myalg.processName = '2f_z_eehiq'
myalg.processID = 15780
myalg.targetLumi = 1000.0
myalg.root_output_file = myalg.root
myalg.RecoParticleColl = 'PandoraPFOs'
myalg.IsolatedLeptonsColl = 'IsolatedLeptons'
myalg.EventHeaderColl = 'EventHeader'
myalg.MCParticleColl = 'MCParticlesSkimmed'
myalg.JetFinderColl = 'MyJets'
myalg.OutputLevel = INFO


# adding my jet finder here:
myJetFinder = MarlinProcessorWrapper("MyJetFinder")
myJetFinder.ProcessorType = "FastJetProcessor"
myJetFinder.Parameters = {
    "algorithm": ["ValenciaPlugin", "1.2", "1.0", "0.7"],
    "clusteringMode": ["ExclusiveNJets", "2"],
    "jetOut": ["MyJets"],
    "recParticleIn": ["PandoraPFOs"],
    "recParticleOut": ["PFOsFromJets"],
    "recombinationScheme": ["E_scheme"],
    "storeParticlesInJets": ["true"],
}
myJetFinder.OutputLevel = INFO


e2lConv = EDM4hep2LcioTool("EDM4hep2Lcio")
l2eConv = Lcio2EDM4hepTool("Lcio2EDM4hep")
e2lConv.convertAll = False
e2lConv.collNameMapping = {
        'PrimaryVertex':'PrimaryVertex',
        'PandoraPFOs':'PandoraPFOs',
        "PandoraClusters": "PandoraClusters",
        "MarlinTrkTracks": "MarlinTrkTracks",
        "EventHeader": "EventHeader",
        "MCParticlesSkimmed": "MCParticlesSkimmed"
        }


l2eConv.convertAll = False
l2eConv.collNameMapping = {
     'PandoraPFOs': 'PandoraPFOs',
     'MyJets': 'MyJets',
     'PFOsFromJets': 'PFOsFromJets',
     "PandoraClusters": "PandoraClusters",
     "MarlinTrkTracks": "MarlinTrkTracks",
     "EventHeader": "EventHeader",
     "MCParticlesSkimmed": "MCParticlesSkimmed"
}

e2lConv.OutputLevel = INFO
l2eConv.OutputLevel = INFO

myJetFinder.EDM4hep2LcioTool = e2lConv
myJetFinder.Lcio2EDM4hepTool = l2eConv


# Isolated Lepton Processor
myIsolatedLeptonTaggingProcessor = MarlinProcessorWrapper("MyIsolatedLeptonTaggingProcessor")
myIsolatedLeptonTaggingProcessor.OutputLevel = INFO
myIsolatedLeptonTaggingProcessor.ProcessorType = "IsolatedLeptonTaggingProcessor"
myIsolatedLeptonTaggingProcessor.Parameters = {
                                               "CosConeLarge": ["0.95"],
                                               "CosConeSmall": ["0.98"],
                                               "CutOnTheISOElectronMVA": ["2.0"],
                                               "CutOnTheISOMuonMVA": ["0.7"],
                                               "DirOfISOElectronWeights": ["/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinReco/v01-32/Analysis/IsolatedLeptonTagging/example/isolated_electron_weights"],
                                               "DirOfISOMuonWeights": ["/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinReco/v01-32/Analysis/IsolatedLeptonTagging/example/isolated_muon_weights_woYoke"],
                                               "InputPandoraPFOsCollection": ["PandoraPFOs"],
                                               "InputPrimaryVertexCollection": ["PrimaryVertex"],
                                               "IsSelectingOneIsoLep": ["false"],
                                               "MaxD0SigForElectron": ["50"],
                                               "MaxD0SigForMuon": ["20"],
                                               "MaxEOverPForElectron": ["1.3"],
                                               "MaxEOverPForMuon": ["0.3"],
                                               "MaxZ0SigForElectron": ["50"],
                                               "MaxZ0SigForMuon": ["20"],
                                               "MinEOverPForElectron": ["0.5"],
                                               "MinEecalOverTotEForElectron": ["0.9"],
                                               "MinEyokeForMuon": ["1.2"],
                                               "MinPForElectron": ["5"],
                                               "MinPForMuon": ["5"],
                                               "OutputIsoLeptonsCollection": ["IsolatedLeptons"],
                                               "OutputPFOsWithoutIsoLepCollection": ["PandoraPFOsWithoutIsoLep"],
                                               "UseYokeForMuonID": ["false"]
                                               }


edm4hep2LcioConv = EDM4hep2LcioTool("EDM4hep2Lcio")
lcio2edm4hepConv = Lcio2EDM4hepTool("Lcio2EDM4hep")
edm4hep2LcioConv.convertAll = False
edm4hep2LcioConv.collNameMapping = {
        'PrimaryVertex':'PrimaryVertex',
        'PandoraPFOs':'PandoraPFOs',
        "PandoraClusters": "PandoraClusters",
        "MarlinTrkTracks": "MarlinTrkTracks",
        "EventHeader": "EventHeader"
      }


lcio2edm4hepConv.convertAll = False
lcio2edm4hepConv.collNameMapping = {
     'PandoraPFOs': 'PandoraPFOs',
     'IsolatedLeptons': 'IsolatedLeptons',
     'PandoraPFOsWithoutIsoLep': 'PandoraPFOsWithoutIsoLep',
     "PandoraClusters": "PandoraClusters",
     "MarlinTrkTracks": "MarlinTrkTracks",
     "EventHeader": "EventHeader",
     "MCParticlesSkimmed": "MCParticlesSkimmed",
     'MyJets': 'MyJets'
     }

#MyIsolatedLeptonTaggingProcessor.EDM4hep2LcioTool = edm4hep2LcioConv
MyIsolatedLeptonTaggingProcessor.Lcio2EDM4hepTool = lcio2edm4hepConv




# setting up a dummy processor to avoid global variable issue
monitor = MarlinProcessorWrapper("EventNumber")
monitor.ProcessorType = "Statusmonitor"
monitor.Parameters = {"HowOften": ["1"], "Verbosity": ["MESSAGE"]}


alg_list.extend([myIsolatedLeptonTaggingProcessor, MyJetFinder, myalg])

io_handler.finalize_converters()

from k4FWCore import ApplicationMgr
ApplicationMgr( TopAlg = alg_list,
                EvtSel="NONE",
                EvtMax=3,
                ExtSvc=[evt_svc],
                OutputLevel=INFO,
               )
