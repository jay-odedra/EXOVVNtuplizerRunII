from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_SingleMu.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.totalUnits = 400

#config.Data.outLFNDirBase = '/store/user/jodedra/Run3/' #% (getUsernameFromSiteDB())
#config.Data.publication = True
#config.Data.outputDatasetTag = 'winter21'
config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.ignoreGlobalBlacklist = True
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_*']
#config.Site.blacklist = ['T2_US_Purdue']


config.General.workArea = 'crab_20220802/ParkingSingleMu'

#name = 'EphemeralZeroBias8'
#name = 'ZeroBias10'

name='/ParkingSingleMuon2/Run2022C-PromptReco-v1/MINIAODBPeak'
config.Data.outLFNDirBase = '/store/user/jodedra/Run3/ParkingSingleMuon/'
config.Data.publication = True
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/ParkingSingleMuon2/Run2022C-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'SUMMER22'


#config.Data.inputDataset = '/' + name + '/ytakahas-winter21-30e4efebaefe52740b9ca928c2409cd7/USER'


#for a in [1,2,3,4,5,6,7]:
    
#
#'/EphemeralZeroBias1/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias2/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias3/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias4/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias5/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias6/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias7/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias8/Run2018D-PromptReco-v2/MINIAOD'
#
#'/EphemeralZeroBias1/Run2018D-v1/RAW'
#'/EphemeralZeroBias2/Run2018D-v1/RAW'
#'/EphemeralZeroBias3/Run2018D-v1/RAW'
#'/EphemeralZeroBias4/Run2018D-v1/RAW'
#'/EphemeralZeroBias5/Run2018D-v1/RAW'
#'/EphemeralZeroBias6/Run2018D-v1/RAW'
#'/EphemeralZeroBias7/Run2018D-v1/RAW'
#'/EphemeralZeroBias8/Run2018D-v1/RAW'




#config.General.requestName = 'EphemeralZeroBias3'
#config.Data.inputDataset = '/EphemeralZeroBias3/Run2018D-PromptReco-v2/MINIAOD'
#config.Data.secondaryInputDataset = '/EphemeralZeroBias3/Run2018D-v1/RAW'



