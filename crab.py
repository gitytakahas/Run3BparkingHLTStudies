from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt_devpath2.py'
config.JobType.maxMemoryMB = 2800
config.JobType.numCores = 4

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 8
#config.Data.totalUnits=1000

config.Data.outLFNDirBase = '' #'/store/user/ytakahas/Trigger/Ephermeral_data'
config.Data.publication = False
config.Data.outputDatasetTag = 'winter21'
config.Site.storageSite = 'T2_CH_CSCS'
config.Site.ignoreGlobalBlacklist = True

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    
    config.General.workArea = 'crab_folder_June2021'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)
 

    #config.General.requestName = 'DYToLLM50'
    #config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
    #submit(config)

    #config.General.requestName = 'QCDPt300To470'
    #config.Data.inputDataset = '/QCD_Pt-300To470_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
    #submit(config)

    #config.General.requestName = 'ZprimeToEE6TeV'
    #config.Data.inputDataset = '/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU0to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
    #submit(config)

    for ii in range(1, 2):

        dname = 'EphemeralZeroBias' + str(ii)

        config.General.requestName = dname
        config.Data.inputDataset = '/' + dname + '/Run2018D-v1/RAW'
        config.Data.outLFNDirBase = '/store/user/ytakahas/Trigger/' + dname
        submit(config)


