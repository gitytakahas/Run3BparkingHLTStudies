pnfs="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/"

outdir="job"

#python getDataset.py --chunk 10 --jdir ${outdir} --odir ${pnfs} --name mc_4gev_22 --filelist fileList/mc.list --type mc
#python getDataset.py --chunk 1 --jdir ${outdir} --odir ${pnfs} --name data_4gev_22 --filelist fileList/datalist8 --type data


# rate estimate only!
python getDataset_noedm.py --chunk 20 --jdir ${outdir} --odir ${pnfs} --name data_4gev_22_rate --filelist ${pnfs}/job/data_4gev_22

