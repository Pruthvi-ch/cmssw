import FWCore.ParameterSet.Config as cms

process = cms.Process("Track")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        #'file:/home/idas/t3store3/root_files/TTJets_8TeV_53X.root'
        'file:/home/idas/test/cmssw/CMSSW_11_3_3/src/BPH-RunII/BPH-RunIIAutumn18DR_step2.root'
        #'file:/home/idas/test/cmssw/CMSSW_11_3_3/src/root_files/SingleMuFlatPt2To100_test3.root'
        #'file:/home/idas/test/cmssw/CMSSW_11_3_3/src/CloseByParticleGun/23293.0_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root'
    )
#     inputCommands=cms.untracked.vstring(
# #        'keep *',
#         'drop recoCaloTauDiscriminator_*__RECO',
#         'drop recoCaloTauTagInfos_*__RECO',
#         'drop recoCaloTaus_*__RECO'
#     )
)

# tracks = cms.PSet( 
#     s = cms.string("generalTracks"),
#     u = cms.string(""),
#     t = cms.string("HLT"),
#     #v = cms.string( 'thing one', "thing two")
# )

process.gen = cms.EDAnalyzer('TrackAnalyzer',
                             tracks1 = cms.untracked.InputTag("generalTracks")
                         )

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('histo.root')
 )

process.p = cms.Path(process.gen)
