import subprocess, os, sys
from optparse import OptionParser, OptionValueError

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


cdir = os.getcwd()

import datetime

d_today = datetime.datetime.now()

d_today = str(d_today).split('.')[0]
d_today = d_today.replace(' ', '-').replace(':','')

print d_today



usage = "usage: python getDataset.py"
parser = OptionParser(usage)

parser.add_option("-d", "--filelist", default="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/data/", type="string", help="filelist", dest="filelist")
parser.add_option("-c", "--chunk", default=1, type="int", help="chunk", dest="chunk")
parser.add_option("-j", "--jdir", default="job", type="string", help="job output dir", dest="jdir")
parser.add_option("-o", "--odir", default="job", type="string", help="output dir", dest="odir")
parser.add_option("-n", "--name", default="None", type="string", help="name", dest="name")

(options, args) = parser.parse_args()


#listoffiles = []
#
#for line in open(options.filelist):
#
#    if line.find('.root')==-1: continue
#
#    line = line.rstrip()
#
#    listoffiles.append(line)

#print listoffiles

#print options.filelist

import glob

#listoffiles = glob.glob(options.filelist + '/EDM*.root')
listoffiles = glob.glob(options.filelist + '/Myroot*.root')

#for ii in listoffiles:
#    print ii


print len(listoffiles), 'files detected'

#import pdb; pdb.set_trace()

listoffiles = list(chunks(listoffiles,options.chunk))

print len(listoffiles), 'chunks are created'

jobdir = cdir + '/' + options.jdir + '/' + options.name #+ '_' + options.analysis + '_' + options.type
outdir = options.odir + '/' + options.jdir + '/' + options.name

ensureDir(outdir)

if os.path.isdir(jobdir):
    ans = raw_input("Directory " + jobdir + " already exists. Delete?\n")
    
    if ans in ['Y', 'yes', 'Yes', 'y']:
        print 'deleted'
        os.system('rm -rf ' + jobdir)
    else:
        print 'quit...'
        sys.exit(1)

ensureDir(jobdir)


for ijob, filename in enumerate(listoffiles):


#    if ijob > 2: break

    
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
    outfile = outdir + '/Myroot_' + str(ijob) + '_rate.root'
#    edmoutfile = outdir + '/EDM_' + str(ijob) + '.root'

    os.system("cp job_template_noedm.sh " + jobscript)
    
    with open(jobscript) as f:
        data_lines = f.read()
        
    inputcommand = ' '.join(filename)
    
#    for ff in filename:
#        inputcommand += 'inputFiles=' + ff + ' '


    data_lines = data_lines.replace('INPUT', inputcommand).replace('FINALOUTPUT', 'Myroot_' + str(ijob) + '_rate.root').replace('OUTFILE', outfile)
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)

    command = 'sbatch -p short --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
#    command = 'sbatch -p standard --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
    os.system(command)
