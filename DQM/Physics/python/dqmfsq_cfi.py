import FWCore.ParameterSet.Config as cms

DQMFSQ= cms.EDAnalyzer('DQMFSQ',
                     triggerResultsCollection = cms.InputTag("TriggerResults", "", "HLT"),
                     HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
                     LabelPFJet       = cms.string("ak4PFJetsCHS"),
                     LabelCastorJet   = cms.string("ak5CastorJets"),
                     LabelTrack       = cms.string("generalTracks"),
                     LabelBeamSpot    = cms.string("offlineBeamSpot"),
                     LabelVertex      = cms.string("offlinePrimaryVertices"),
                     
)
