import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
process = cms.Process('PROD',Phase2C11)

#process.load("SimGeneral.HepPDTESSource.pdt_cfi")
#process.load("Configuration.Geometry.GeometryExtended2026D76Reco_cff")

#from Configuration.Eras.Modifier_phase2_hgcalV12_cff import phase2_hgcalV12
#process = cms.Process('PROD',phase2_hgcalV12)

# from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
# process = cms.Process('PROD',Phase2C9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D83Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/home/idas/t3store3/root_files/SingleMuFlatPt2To100_cfi_py_GEN_SIM.root'
        #'file:/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-09-08-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default.root'
        #'file:/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-09-08-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D83.root'
        #'file:/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-09-08-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D86.root'
        #'file:/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-09-08-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D86_higheta.root'
        #'file:/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-09-08-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D83_higheta.root'
        'file:/home/idas/t3store3/root_files/SingleMuFlatPt2To100_cfi_py_GEN_geo_default.root'
        #'file:/home/idas/t3store3/root_files/SingleMuFlatPt2To100_cfi_py_GEN_geo_mod.root'
        #'file:/home/idas/test/cmssw/CMSSW_12_1_X_2021-09-13-2300/src/SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D86_higheta.root'
    )
#     inputCommands=cms.untracked.vstring(
# #        'keep *',
#         'drop recoCaloTauDiscriminator_*__RECO',
#         'drop recoCaloTauTagInfos_*__RECO',
#         'drop recoCaloTaus_*__RECO'
#     )
)


process.prodEE = cms.EDAnalyzer('GeantRead',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             Detector   = cms.string("HGCalEESensitive"),
                         )


process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
)

process.prodHEB = process.prodHEF.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
    Detector   = cms.string("HGCalHEScintillatorSensitive"),
)

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('geantoutput.root')
 )

process.p = cms.Path(process.prodEE*process.prodHEF*process.prodHEB)
#process.p = cms.Path(process.prodEE)
