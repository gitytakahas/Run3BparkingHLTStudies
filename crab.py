from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'purity.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 4

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits=10

config.Data.outLFNDirBase = '/store/user/ytakahas/purityStudy/' #% (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'winter21'
config.Site.storageSite = 'T3_CH_PSI'
#config.Site.ignoreGlobalBlacklist = True
config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_US_Caltech', 'T2_US_MIT', 'T2_US_Florida', 'T2_US_Vanderbilt', 'T2_US_Nebraska', 'T2_US_Wisconsin']
config.Site.blacklist = ['T2_US_Purdue']


config.General.workArea = 'crab_Nov2021_massprod'

config.General.requestName = 'EphemeralZeroBias8'
config.Data.inputDataset = '/EphemeralZeroBias8/Run2018D-PromptReco-v2/MINIAOD'
config.Data.secondaryInputDataset = '/EphemeralZeroBias8/Run2018D-v1/RAW'

#config.General.requestName = 'EphemeralZeroBias3'
#config.Data.inputDataset = '/EphemeralZeroBias3/Run2018D-PromptReco-v2/MINIAOD'
#config.Data.secondaryInputDataset = '/EphemeralZeroBias3/Run2018D-v1/RAW'



