#!/bin/bash
# 
#SBATCH -p standard
#SBATCH --account=t3
#SBATCH --time 11:59:00
#SBATCH -e cn-test.err  # replace default slurm-SLURM_JOB_ID.err
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out

echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#
#if [ -r CMSSW_12_0_0_pre4/src ] ; then
#    echo "release CMSSW_12_0_0_pre4 already exists"
#else
#    scram p CMSSW CMSSW_12_0_0_pre4
#fi
#
#cd CMSSW_12_0_0_pre4/src
#eval `scram runtime -sh`
#cd -

echo 'cmssw release = ' $CMSSW_BASE

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}


#cmsRun RTYPE.py INPUT maxEvents=-1 outputFile=$TMPDIR/TMPOUTPUT
python3 /work/ytakahas/work/HLT/CMSSW_12_0_0_pre4/src/Analysis/HLTAnalyserPy/test/makeRun3Ntup.py INPUT -o $TMPDIR/FINALOUTPUT


#xrdcp -f $TMPDIR/TMPOUTPUT root://t3dcachedb03.psi.ch/EDMOUTFILE
xrdcp -f $TMPDIR/FINALOUTPUT root://t3dcachedb03.psi.ch/OUTFILE


rm -rf /scratch/$USER/${SLURM_JOB_ID}

date
