##################################################################################
# Nanoanalyzer configuration file for all years                                  #
# Use for HT Condor and VM                                                       #
# Uncomment & comment relevant lines before you run it                           #
##################################################################################

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Nano")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff") # for 2011
#process.load("Configuration.StandardSequences.Geometry_cff") # for 2010
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

################################### VARPASSING ###################################
# Use VarParsing to specify your input directory and output file
# Comment Varparsing if you want to submit job using CRAB
# because CRAB does not support input directory as we put input dataset directly
#options = VarParsing('analysis')
#options.register(
#    "inputDir",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Input directory with inputs"
#    )
#options.register(
#    "outputName",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Output name"
#    )
#options.parseArguments()

#if (options.inputDir == ""):
#    sys.exit("Directory to find input file where???")
#else:
    # 2010 VM
    #InputDir = "/home/cms-opendata/CMSSW_4_2_8/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2011 NAF/HTC
    # InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_5_3_32/src/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2015 NAF/HTC
    #InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/" + options.inputDir

##################################################################################

###################################### GLOBAL TAG ################################
# Change the global tag accordingly
# ParkingBPH 2021 UL
#process.GlobalTag.globaltag = '106X_dataRun2_v35'
# 2021 data
#process.GlobalTag.globaltag = '120X_dataRun3_Prompt_v2'
#process.GlobalTag.globaltag = '121X_dataRun3_v13'
# 2022 data
#process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v8'
#process.GlobalTag.globaltag = '123X_dataRun3_Express_v10'
process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v12'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Nano')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#####################################  JSON FILE #################################
# Change the directory and JSON file accordingly
# Only uncomment if you run in Data
# ParkingBPH 2021 UL
#goodJSON = '/nfs/dust/cms/user/yangq2/goodJson/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
# pilot 2021
#goodJSON = './BT21GOOD_withTKon.json'
# 2022A, data taken up to June 3
#goodJSON = './Collisions22AGOOD_withALLON.json'


##################################################################################

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

##################################################################################

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(*fileinPut)
#)

# To submit batch job using CRAB or test locally (2nd option), use this:


inputFiles= [

    '/store/data/Run2022C/ParkingDoubleElectronLowMass1/MINIAOD/PromptReco-v1/000/356/170/00000/07cf47be-67de-4b52-8956-261221ac18a9.root',
    '/store/data/Run2022C/ParkingDoubleElectronLowMass2/MINIAOD/PromptReco-v1/000/356/170/00000/57130a2d-1e3e-4013-9236-e38cdfd81181.root',
    '/store/data/Run2022C/ParkingDoubleElectronLowMass3/MINIAOD/PromptReco-v1/000/356/170/00000/2f693ac8-1454-4889-954c-9a77e07c82a8.root',
    '/store/data/Run2022C/ParkingDoubleElectronLowMass4/MINIAOD/PromptReco-v1/000/356/170/00000/8e4588b3-e391-4e2b-8507-a72dad99a244.root',
    '/store/data/Run2022C/ParkingDoubleElectronLowMass5/MINIAOD/PromptReco-v1/000/356/170/00000/e247937b-47c5-4088-b39d-e0e631882072.root',
    '/store/data/Run2022C/ParkingDoubleElectronLowMass0/MINIAOD/PromptReco-v1/000/356/170/00000/45c0f2ed-eb5b-4292-abc8-3117424d9432.root'
#    '/store/data/Run2022C/HLTPhysics/MINIAOD/PromptReco-v1/000/355/872/00000/725217a6-902f-48e4-84ff-ec18ec794c66.root',
#    '/store/data/Run2022C/HLTPhysics/MINIAOD/PromptReco-v1/000/355/872/00000/97c45763-b1e9-490d-8f99-f5133b889fa6.root',
#    '/store/data/Run2022C/HLTPhysics/MINIAOD/PromptReco-v1/000/355/872/00000/cc0ee657-ec41-4470-a3c6-40d386088c1e.root'
]

#'/store/data/Run2022C/HLTPhysics/MINIAOD/PromptReco-v1/000/356/003/00000/194e1847-a63a-44e0-805a-30c56ae0f48a.root'


process.source = cms.Source("PoolSource",
# for crab
                            fileNames = cms.untracked.vstring(inputFiles)
)


##################################################################################

# Process the lumi
# Only uncomment if you run in Data
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('test.root')
                                   )


# Process the analyzer
process.nano = cms.EDAnalyzer('NanoAnalyzer',
                              electrons = cms.InputTag("slimmedElectrons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              packedpfcandidates = cms.InputTag('packedPFCandidates'),
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              # Change this:
                              # If HTC/VM:
                              #outFile = cms.string(options.outputName),
                              # If interactive:
                              #outFile = cms.string('test.root'), 
                              # If CRAB:
                              # make sure the name is same as the crab config
#                              outFile = cms.string('22ZeroBias136.root'),

                              # Change this:
                              # If MC:
                              #isData = cms.bool(False)
                              # If Data:
#                              isData = cms.bool(True)
)
process.p = cms.Path(process.nano)
