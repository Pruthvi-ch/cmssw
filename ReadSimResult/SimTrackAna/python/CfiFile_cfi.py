import FWCore.ParameterSet.Config as cms

gen = cms.EDAnalyzer('TrackAnalyzer'
     ,tracks = cms.untracked.InputTag('generalTracks')
)
