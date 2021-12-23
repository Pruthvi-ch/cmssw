import FWCore.ParameterSet.Config as cms


# from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
# process = cms.Process('PROD',Phase2C11)

#process.load("SimGeneral.HepPDTESSource.pdt_cfi")
#process.load("Configuration.Geometry.GeometryExtended2026D76Reco_cff")

#from Configuration.Eras.Modifier_phase2_hgcalV12_cff import phase2_hgcalV12
#process = cms.Process('PROD',phase2_hgcalV12)

# from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
# process = cms.Process('PROD',Phase2C9)

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D83Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D83_higheta.root')
    #fileNames = cms.untracked.vstring('file:SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D86_higheta.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/step1_D86.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/MuFlatPt/D83/SingleMuFlatPt2To100_D83_step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/MuFlatPt/D86/SingleMuFlatPt2To100_D86_step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/CloseByPhoton_Official/38693.D83/step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/CloseByPhoton_Official/38693.D86/step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/MuDeltaPt/D83/step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/MuDeltaPt/D86/step1.root')
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1.root')
                                      #'file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1_0.root'
                                      #'file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1_1.root',
                                      #'file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1_2.root',
                                      #'file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1_3.root',
                                      #'file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D83/step1_4.root'
    #fileNames = cms.untracked.vstring('file:/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D86/step1.root')
    #fileNames = cms.untracked.vstring('file:/eos/user/p/psuryade/SimOut/DeltaPt/step1.root')
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/p/psuryade/private/CMSSW/CMSSW_12_2_X_2021-11-22-1100/src/step1.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


process.prodEE = cms.EDAnalyzer('ml_ntuple',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             Detector   = cms.string("HGCalEESensitive"),
                         )


process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
)


#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('ml_ntuple.root')
 )

process.p = cms.Path(process.prodEE*process.prodHEF)
