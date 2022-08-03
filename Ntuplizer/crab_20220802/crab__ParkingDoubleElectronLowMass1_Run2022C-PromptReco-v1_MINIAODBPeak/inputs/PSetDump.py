import FWCore.ParameterSet.Config as cms

process = cms.Process("Nano")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2022C/ParkingDoubleElectronLowMass1/MINIAOD/PromptReco-v1/000/356/381/00000/d0ea2b7f-5b63-4ebb-aec8-95d4adc464e6.root')
)
process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
        0.000157, -3e-06
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.maxLuminosityBlocks = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.MEtoEDMConverter = cms.EDProducer("MEtoEDMConverter",
    Frequency = cms.untracked.int32(50),
    MEPathToSave = cms.untracked.string(''),
    Name = cms.untracked.string('MEtoEDMConverter'),
    Verbosity = cms.untracked.int32(0)
)


process.nano = cms.EDAnalyzer("NanoAnalyzer",
    HLT = cms.InputTag("TriggerResults","","HLT"),
    electrons = cms.InputTag("slimmedElectrons"),
    packedpfcandidates = cms.InputTag("packedPFCandidates"),
    triggerobjects = cms.InputTag("slimmedPatTrigger"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.DQMStore = cms.Service("DQMStore")


process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        enable = cms.untracked.bool(True),
        enableStatistics = cms.untracked.bool(True),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.untracked.bool(False),
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.untracked.string('WARNING'),
        threshold = cms.untracked.string('INFO'),
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(False),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.optional.untracked.bool,
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.optional.untracked.string,
        threshold = cms.optional.untracked.string,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    debugModules = cms.untracked.vstring(),
    default = cms.untracked.PSet(
        limit = cms.optional.untracked.int32,
        lineLength = cms.untracked.int32(80),
        noLineBreaks = cms.untracked.bool(False),
        noTimeStamps = cms.untracked.bool(False),
        reportEvery = cms.untracked.int32(1),
        statisticsThreshold = cms.untracked.string('INFO'),
        threshold = cms.untracked.string('INFO'),
        timespan = cms.optional.untracked.int32,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    files = cms.untracked.PSet(
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            enableStatistics = cms.untracked.bool(False),
            extension = cms.optional.untracked.string,
            filename = cms.optional.untracked.string,
            lineLength = cms.optional.untracked.int32,
            noLineBreaks = cms.optional.untracked.bool,
            noTimeStamps = cms.optional.untracked.bool,
            output = cms.optional.untracked.string,
            resetStatistics = cms.untracked.bool(False),
            statisticsThreshold = cms.optional.untracked.string,
            threshold = cms.optional.untracked.string,
            allowAnyLabel_=cms.optional.untracked.PSetTemplate(
                limit = cms.optional.untracked.int32,
                reportEvery = cms.untracked.int32(1),
                timespan = cms.optional.untracked.int32
            )
        )
    ),
    suppressDebug = cms.untracked.vstring(),
    suppressFwkInfo = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    allowAnyLabel_=cms.optional.untracked.PSetTemplate(
        limit = cms.optional.untracked.int32,
        reportEvery = cms.untracked.int32(1),
        timespan = cms.optional.untracked.int32
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('singleeletrigmatch.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL',
        'ZDC',
        'CASTOR',
        'EcalBarrel',
        'EcalEndcap',
        'EcalPreshower',
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService",
    maxExtrapolationTimeInSec = cms.uint32(0)
)


process.EcalLaserCorrectionServiceMC = cms.ESProducer("EcalLaserCorrectionServiceMC",
    appendToDataLabel = cms.string('')
)


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    fromDD4hep = cms.untracked.bool(False),
    fromDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerAdditionalParametersPerDet = cms.ESProducer("TrackerAdditionalParametersPerDetESModule",
    appendToDataLabel = cms.string('')
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
    usePhase2Stacks = cms.bool(False)
)


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.ak10PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL1Fastjet',
        'ak10PFCHSL2Relative',
        'ak10PFCHSL3Absolute',
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL1Offset',
        'ak10PFCHSL2Relative',
        'ak10PFCHSL3Absolute',
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL2Relative',
        'ak10PFCHSL3Absolute'
    )
)


process.ak10PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL2Relative',
        'ak10PFCHSL3Absolute',
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2Relative')
)


process.ak10PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L3Absolute')
)


process.ak10PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak10PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL1Fastjet',
        'ak10PFL2Relative',
        'ak10PFL3Absolute',
        'ak10PFResidual'
    )
)


process.ak10PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL1Offset',
        'ak10PFL2Relative',
        'ak10PFL3Absolute',
        'ak10PFResidual'
    )
)


process.ak10PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL2Relative',
        'ak10PFL3Absolute'
    )
)


process.ak10PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL2Relative',
        'ak10PFL3Absolute',
        'ak10PFResidual'
    )
)


process.ak10PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2Relative')
)


process.ak10PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L3Absolute')
)


process.ak10PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2L3Residual')
)


process.ak1PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL1Fastjet',
        'ak1PFCHSL2Relative',
        'ak1PFCHSL3Absolute',
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL1Offset',
        'ak1PFCHSL2Relative',
        'ak1PFCHSL3Absolute',
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL2Relative',
        'ak1PFCHSL3Absolute'
    )
)


process.ak1PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL2Relative',
        'ak1PFCHSL3Absolute',
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2Relative')
)


process.ak1PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L3Absolute')
)


process.ak1PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak1PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL1Fastjet',
        'ak1PFL2Relative',
        'ak1PFL3Absolute',
        'ak1PFResidual'
    )
)


process.ak1PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL1Offset',
        'ak1PFL2Relative',
        'ak1PFL3Absolute',
        'ak1PFResidual'
    )
)


process.ak1PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL2Relative',
        'ak1PFL3Absolute'
    )
)


process.ak1PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL2Relative',
        'ak1PFL3Absolute',
        'ak1PFResidual'
    )
)


process.ak1PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2Relative')
)


process.ak1PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L3Absolute')
)


process.ak1PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2L3Residual')
)


process.ak2PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL1Fastjet',
        'ak2PFCHSL2Relative',
        'ak2PFCHSL3Absolute',
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL1Offset',
        'ak2PFCHSL2Relative',
        'ak2PFCHSL3Absolute',
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL2Relative',
        'ak2PFCHSL3Absolute'
    )
)


process.ak2PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL2Relative',
        'ak2PFCHSL3Absolute',
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2Relative')
)


process.ak2PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L3Absolute')
)


process.ak2PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak2PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL1Fastjet',
        'ak2PFL2Relative',
        'ak2PFL3Absolute',
        'ak2PFResidual'
    )
)


process.ak2PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL1Offset',
        'ak2PFL2Relative',
        'ak2PFL3Absolute',
        'ak2PFResidual'
    )
)


process.ak2PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL2Relative',
        'ak2PFL3Absolute'
    )
)


process.ak2PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL2Relative',
        'ak2PFL3Absolute',
        'ak2PFResidual'
    )
)


process.ak2PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2Relative')
)


process.ak2PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L3Absolute')
)


process.ak2PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2L3Residual')
)


process.ak3PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL1Fastjet',
        'ak3PFCHSL2Relative',
        'ak3PFCHSL3Absolute',
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL1Offset',
        'ak3PFCHSL2Relative',
        'ak3PFCHSL3Absolute',
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL2Relative',
        'ak3PFCHSL3Absolute'
    )
)


process.ak3PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL2Relative',
        'ak3PFCHSL3Absolute',
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2Relative')
)


process.ak3PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L3Absolute')
)


process.ak3PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak3PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL1Fastjet',
        'ak3PFL2Relative',
        'ak3PFL3Absolute',
        'ak3PFResidual'
    )
)


process.ak3PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL1Offset',
        'ak3PFL2Relative',
        'ak3PFL3Absolute',
        'ak3PFResidual'
    )
)


process.ak3PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL2Relative',
        'ak3PFL3Absolute'
    )
)


process.ak3PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL2Relative',
        'ak3PFL3Absolute',
        'ak3PFResidual'
    )
)


process.ak3PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2Relative')
)


process.ak3PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L3Absolute')
)


process.ak3PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2L3Residual')
)


process.ak4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute',
        'ak4CaloL6SLB'
    )
)


process.ak4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute',
        'ak4CaloResidual'
    )
)


process.ak4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Offset',
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Offset',
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute',
        'ak4CaloResidual'
    )
)


process.ak4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute',
        'ak4CaloL6SLB'
    )
)


process.ak4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative',
        'ak4CaloL3Absolute',
        'ak4CaloResidual'
    )
)


process.ak4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTFastjet',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTFastjet',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute',
        'ak4JPTResidual'
    )
)


process.ak4JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute',
        'ak4JPTResidual'
    )
)


process.ak4JPTL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset',
        'ak4JPTL2Relative',
        'ak4JPTL3Absolute',
        'ak4JPTResidual'
    )
)


process.ak4JPTL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Fastjet')
)


process.ak4L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Offset')
)


process.ak4PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet',
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet',
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute',
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Offset',
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Offset',
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute',
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL2Relative',
        'ak4PFCHSL3Absolute',
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ak4PFL2Relative',
        'ak4PFL3Absolute'
    )
)


process.ak4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ak4PFL2Relative',
        'ak4PFL3Absolute',
        'ak4PFL6SLB'
    )
)


process.ak4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ak4PFL2Relative',
        'ak4PFL3Absolute',
        'ak4PFResidual'
    )
)


process.ak4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Offset',
        'ak4PFL2Relative',
        'ak4PFL3Absolute'
    )
)


process.ak4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Offset',
        'ak4PFL2Relative',
        'ak4PFL3Absolute',
        'ak4PFResidual'
    )
)


process.ak4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative',
        'ak4PFL3Absolute'
    )
)


process.ak4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative',
        'ak4PFL3Absolute',
        'ak4PFL6SLB'
    )
)


process.ak4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative',
        'ak4PFL3Absolute',
        'ak4PFResidual'
    )
)


process.ak4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ak4TrackL2Relative',
        'ak4TrackL3Absolute'
    )
)


process.ak4TrackL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4TrackL2Relative',
        'ak4TrackL3Absolute'
    )
)


process.ak4TrackL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL1Fastjet',
        'ak5PFCHSL2Relative',
        'ak5PFCHSL3Absolute',
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL1Offset',
        'ak5PFCHSL2Relative',
        'ak5PFCHSL3Absolute',
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL2Relative',
        'ak5PFCHSL3Absolute'
    )
)


process.ak5PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL2Relative',
        'ak5PFCHSL3Absolute',
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2Relative')
)


process.ak5PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL1Fastjet',
        'ak5PFL2Relative',
        'ak5PFL3Absolute',
        'ak5PFResidual'
    )
)


process.ak5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL1Offset',
        'ak5PFL2Relative',
        'ak5PFL3Absolute',
        'ak5PFResidual'
    )
)


process.ak5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL2Relative',
        'ak5PFL3Absolute'
    )
)


process.ak5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL2Relative',
        'ak5PFL3Absolute',
        'ak5PFResidual'
    )
)


process.ak5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)


process.ak5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)


process.ak5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.ak6PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL1Fastjet',
        'ak6PFCHSL2Relative',
        'ak6PFCHSL3Absolute',
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL1Offset',
        'ak6PFCHSL2Relative',
        'ak6PFCHSL3Absolute',
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL2Relative',
        'ak6PFCHSL3Absolute'
    )
)


process.ak6PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL2Relative',
        'ak6PFCHSL3Absolute',
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2Relative')
)


process.ak6PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L3Absolute')
)


process.ak6PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL1Fastjet',
        'ak6PFL2Relative',
        'ak6PFL3Absolute',
        'ak6PFResidual'
    )
)


process.ak6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL1Offset',
        'ak6PFL2Relative',
        'ak6PFL3Absolute',
        'ak6PFResidual'
    )
)


process.ak6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL2Relative',
        'ak6PFL3Absolute'
    )
)


process.ak6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL2Relative',
        'ak6PFL3Absolute',
        'ak6PFResidual'
    )
)


process.ak6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2Relative')
)


process.ak6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L3Absolute')
)


process.ak6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2L3Residual')
)


process.ak7CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset',
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute',
        'ak7CaloL6SLB'
    )
)


process.ak7CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Fastjet',
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute',
        'ak7CaloResidual'
    )
)


process.ak7CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak7CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset',
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset',
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute',
        'ak7CaloResidual'
    )
)


process.ak7CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute',
        'ak7CaloL6SLB'
    )
)


process.ak7CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative',
        'ak7CaloL3Absolute',
        'ak7CaloResidual'
    )
)


process.ak7CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2Relative')
)


process.ak7CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7CaloJetsSoftMuonTagInfos")
)


process.ak7CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2L3Residual')
)


process.ak7JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTFastjet',
        'ak7JPTL2Relative',
        'ak7JPTL3Absolute'
    )
)


process.ak7JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTFastjet',
        'ak7JPTL2Relative',
        'ak7JPTL3Absolute',
        'ak7JPTResidual'
    )
)


process.ak7JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset',
        'ak7JPTL2Relative',
        'ak7JPTL3Absolute'
    )
)


process.ak7JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset',
        'ak7JPTL2Relative',
        'ak7JPTL3Absolute',
        'ak7JPTResidual'
    )
)


process.ak7JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset',
        'ak7JPTL2Relative',
        'ak7JPTL3Absolute'
    )
)


process.ak7L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Fastjet')
)


process.ak7L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Offset')
)


process.ak7PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet',
        'ak7PFCHSL2Relative',
        'ak7PFCHSL3Absolute'
    )
)


process.ak7PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL1Fastjet',
        'ak7PFCHSL2Relative',
        'ak7PFCHSL3Absolute',
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL1Offset',
        'ak7PFCHSL2Relative',
        'ak7PFCHSL3Absolute',
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL2Relative',
        'ak7PFCHSL3Absolute'
    )
)


process.ak7PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL2Relative',
        'ak7PFCHSL3Absolute',
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2Relative')
)


process.ak7PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L3Absolute')
)


process.ak7PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak7PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ak7PFL2Relative',
        'ak7PFL3Absolute'
    )
)


process.ak7PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ak7PFL2Relative',
        'ak7PFL3Absolute',
        'ak7PFL6SLB'
    )
)


process.ak7PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Fastjet',
        'ak7PFL2Relative',
        'ak7PFL3Absolute',
        'ak7PFResidual'
    )
)


process.ak7PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Offset',
        'ak7PFL2Relative',
        'ak7PFL3Absolute'
    )
)


process.ak7PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Offset',
        'ak7PFL2Relative',
        'ak7PFL3Absolute',
        'ak7PFResidual'
    )
)


process.ak7PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative',
        'ak7PFL3Absolute'
    )
)


process.ak7PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative',
        'ak7PFL3Absolute',
        'ak7PFL6SLB'
    )
)


process.ak7PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative',
        'ak7PFL3Absolute',
        'ak7PFResidual'
    )
)


process.ak7PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2Relative')
)


process.ak7PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L3Absolute')
)


process.ak7PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7PFJetsSoftMuonTagInfos")
)


process.ak7PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2L3Residual')
)


process.ak8PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL1Fastjet',
        'ak8PFCHSL2Relative',
        'ak8PFCHSL3Absolute',
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL1Offset',
        'ak8PFCHSL2Relative',
        'ak8PFCHSL3Absolute',
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL2Relative',
        'ak8PFCHSL3Absolute'
    )
)


process.ak8PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL2Relative',
        'ak8PFCHSL3Absolute',
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2Relative')
)


process.ak8PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L3Absolute')
)


process.ak8PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak8PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL1Fastjet',
        'ak8PFL2Relative',
        'ak8PFL3Absolute',
        'ak8PFResidual'
    )
)


process.ak8PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL1Offset',
        'ak8PFL2Relative',
        'ak8PFL3Absolute',
        'ak8PFResidual'
    )
)


process.ak8PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL2Relative',
        'ak8PFL3Absolute'
    )
)


process.ak8PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL2Relative',
        'ak8PFL3Absolute',
        'ak8PFResidual'
    )
)


process.ak8PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2Relative')
)


process.ak8PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L3Absolute')
)


process.ak8PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2L3Residual')
)


process.ak9PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL1Fastjet',
        'ak9PFCHSL2Relative',
        'ak9PFCHSL3Absolute',
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL1Offset',
        'ak9PFCHSL2Relative',
        'ak9PFCHSL3Absolute',
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL2Relative',
        'ak9PFCHSL3Absolute'
    )
)


process.ak9PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL2Relative',
        'ak9PFCHSL3Absolute',
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2Relative')
)


process.ak9PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L3Absolute')
)


process.ak9PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak9PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL1Fastjet',
        'ak9PFL2Relative',
        'ak9PFL3Absolute',
        'ak9PFResidual'
    )
)


process.ak9PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL1Offset',
        'ak9PFL2Relative',
        'ak9PFL3Absolute',
        'ak9PFResidual'
    )
)


process.ak9PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL2Relative',
        'ak9PFL3Absolute'
    )
)


process.ak9PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL2Relative',
        'ak9PFL3Absolute',
        'ak9PFResidual'
    )
)


process.ak9PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2Relative')
)


process.ak9PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L3Absolute')
)


process.ak9PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2L3Residual')
)


process.ctppsBeamParametersFromLHCInfoESSource = cms.ESProducer("CTPPSBeamParametersFromLHCInfoESSource",
    appendToDataLabel = cms.string(''),
    beamDivX45 = cms.double(0.1),
    beamDivX56 = cms.double(0.1),
    beamDivY45 = cms.double(0.1),
    beamDivY56 = cms.double(0.1),
    lhcInfoLabel = cms.string(''),
    vtxOffsetX45 = cms.double(0.01),
    vtxOffsetX56 = cms.double(0.01),
    vtxOffsetY45 = cms.double(0.01),
    vtxOffsetY56 = cms.double(0.01),
    vtxOffsetZ45 = cms.double(0.01),
    vtxOffsetZ56 = cms.double(0.01),
    vtxStddevX = cms.double(0.02),
    vtxStddevY = cms.double(0.02),
    vtxStddevZ = cms.double(0.02)
)


process.ctppsInterpolatedOpticalFunctionsESSource = cms.ESProducer("CTPPSInterpolatedOpticalFunctionsESSource",
    appendToDataLabel = cms.string(''),
    lhcInfoLabel = cms.string(''),
    opticsLabel = cms.string('')
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.ic5CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset',
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute',
        'ic5CaloL6SLB'
    )
)


process.ic5CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Fastjet',
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute',
        'ic5CaloResidual'
    )
)


process.ic5CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ic5CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset',
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset',
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute',
        'ic5CaloResidual'
    )
)


process.ic5CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute',
        'ic5CaloL6SLB'
    )
)


process.ic5CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative',
        'ic5CaloL3Absolute',
        'ic5CaloResidual'
    )
)


process.ic5CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2Relative')
)


process.ic5CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L3Absolute')
)


process.ic5CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5CaloJetsSoftMuonTagInfos")
)


process.ic5CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2L3Residual')
)


process.ic5PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ic5PFL2Relative',
        'ic5PFL3Absolute'
    )
)


process.ic5PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'ic5PFL2Relative',
        'ic5PFL3Absolute',
        'ic5PFL6SLB'
    )
)


process.ic5PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Fastjet',
        'ic5PFL2Relative',
        'ic5PFL3Absolute',
        'ic5PFResidual'
    )
)


process.ic5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ic5PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Offset',
        'ic5PFL2Relative',
        'ic5PFL3Absolute'
    )
)


process.ic5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Offset',
        'ic5PFL2Relative',
        'ic5PFL3Absolute',
        'ic5PFResidual'
    )
)


process.ic5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative',
        'ic5PFL3Absolute'
    )
)


process.ic5PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative',
        'ic5PFL3Absolute',
        'ic5PFL6SLB'
    )
)


process.ic5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative',
        'ic5PFL3Absolute',
        'ic5PFResidual'
    )
)


process.ic5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2Relative')
)


process.ic5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L3Absolute')
)


process.ic5PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5PFJetsSoftMuonTagInfos")
)


process.ic5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2L3Residual')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.kt4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset',
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute',
        'kt4CaloL6SLB'
    )
)


process.kt4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Fastjet',
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute',
        'kt4CaloResidual'
    )
)


process.kt4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset',
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset',
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute',
        'kt4CaloResidual'
    )
)


process.kt4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute',
        'kt4CaloL6SLB'
    )
)


process.kt4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative',
        'kt4CaloL3Absolute',
        'kt4CaloResidual'
    )
)


process.kt4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2Relative')
)


process.kt4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L3Absolute')
)


process.kt4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4CaloJetsSoftMuonTagInfos")
)


process.kt4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2L3Residual')
)


process.kt4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'kt4PFL2Relative',
        'kt4PFL3Absolute'
    )
)


process.kt4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'kt4PFL2Relative',
        'kt4PFL3Absolute',
        'kt4PFL6SLB'
    )
)


process.kt4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Fastjet',
        'kt4PFL2Relative',
        'kt4PFL3Absolute',
        'kt4PFResidual'
    )
)


process.kt4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Offset',
        'kt4PFL2Relative',
        'kt4PFL3Absolute'
    )
)


process.kt4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Offset',
        'kt4PFL2Relative',
        'kt4PFL3Absolute',
        'kt4PFResidual'
    )
)


process.kt4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative',
        'kt4PFL3Absolute'
    )
)


process.kt4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative',
        'kt4PFL3Absolute',
        'kt4PFL6SLB'
    )
)


process.kt4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative',
        'kt4PFL3Absolute',
        'kt4PFResidual'
    )
)


process.kt4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2Relative')
)


process.kt4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L3Absolute')
)


process.kt4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4PFJetsSoftMuonTagInfos")
)


process.kt4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2L3Residual')
)


process.kt6CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet',
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset',
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute',
        'kt6CaloL6SLB'
    )
)


process.kt6CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Fastjet',
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute',
        'kt6CaloResidual'
    )
)


process.kt6CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt6CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset',
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset',
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute',
        'kt6CaloResidual'
    )
)


process.kt6CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute',
        'kt6CaloL6SLB'
    )
)


process.kt6CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative',
        'kt6CaloL3Absolute',
        'kt6CaloResidual'
    )
)


process.kt6CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2Relative')
)


process.kt6CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L3Absolute')
)


process.kt6CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6CaloJetsSoftMuonTagInfos")
)


process.kt6CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2L3Residual')
)


process.kt6PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'kt6PFL2Relative',
        'kt6PFL3Absolute'
    )
)


process.kt6PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet',
        'kt6PFL2Relative',
        'kt6PFL3Absolute',
        'kt6PFL6SLB'
    )
)


process.kt6PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Fastjet',
        'kt6PFL2Relative',
        'kt6PFL3Absolute',
        'kt6PFResidual'
    )
)


process.kt6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt6PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Offset',
        'kt6PFL2Relative',
        'kt6PFL3Absolute'
    )
)


process.kt6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Offset',
        'kt6PFL2Relative',
        'kt6PFL3Absolute',
        'kt6PFResidual'
    )
)


process.kt6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative',
        'kt6PFL3Absolute'
    )
)


process.kt6PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative',
        'kt6PFL3Absolute',
        'kt6PFL6SLB'
    )
)


process.kt6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative',
        'kt6PFL3Absolute',
        'kt6PFResidual'
    )
)


process.kt6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2Relative')
)


process.kt6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L3Absolute')
)


process.kt6PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6PFJetsSoftMuonTagInfos")
)


process.kt6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2L3Residual')
)


process.multipleScatteringParametrisationMakerESProducer = cms.ESProducer("MultipleScatteringParametrisationMakerESProducer",
    appendToDataLabel = cms.string('')
)


process.muonGeometryConstants = cms.ESProducer("MuonGeometryConstantsESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(False)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    ),
    siPixelQualityLabel = cms.string('')
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ),
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('123X_dataRun3_Prompt_v12'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ),
        cms.PSet(
            slope = cms.double(-1.5610227),
            tmax = cms.double(10.0),
            tzero = cms.double(11.977461)
        ),
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(208),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(208),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.prefer("es_hardcode")

process.endOfProcess = cms.Sequence(process.MEtoEDMConverter)


process.p = cms.Path(process.nano)


