pnfs="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/"

outdir="job"


# creating HLTAnalysis trees 
#python getDataset.py --chunk 4 --jdir ${outdir} --odir ${pnfs} --name Winter21_official_mc_v2 --filelist fileList/Run3Winter21DR.txt --type mc
#python getDataset.py --chunk 4 --jdir ${outdir} --odir ${pnfs} --name Winter21_official_data_all --filelist fileList/datalist_all --type data


# create EgHLT trees ... 
#python getDataset_noedm_v2.py --chunk 3 --jdir ${outdir} --odir ${pnfs} --name Winter21_official_data_all --filelist $pnfs/../purityStudy_CMSSW12_1_0 --type data



# rate estimate based on HLTAnalysis trees
python getDataset_noedm.py --chunk 3 --jdir ${outdir} --odir ${pnfs} --name HLT_data_rate_all --filelist ${pnfs}/job/Winter21_official_data_all --type rate

# eff. estimate based on HLTAnalysis trees
#python getDataset_noedm.py --chunk 10 --jdir ${outdir} --odir ${pnfs} --name HLT_mc_Eff --filelist ${pnfs}/job/Winter21_official_mc --type eff_hlt

# data dist.
#python getDataset_noedm.py --chunk 10 --jdir ${outdir} --odir ${pnfs} --name HLT_data_dist --filelist ${pnfs}/job/Winter21_official_data --type runTauDisplay_BcJpsiTauNu_data

# mc dist.
#python getDataset_noedm.py --chunk 20 --jdir ${outdir} --odir ${pnfs} --name HLT_mc_dist --filelist ${pnfs}/job/Winter21_official_mc --type runTauDisplay_BcJpsiTauNu

