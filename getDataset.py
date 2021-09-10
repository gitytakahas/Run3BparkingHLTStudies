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

test=[
1100, 
1101, 
1102, 
1103, 
1104, 
1105, 
1106, 
1107, 
1108, 
1109, 
110, 
1110, 
1111, 
1112, 
1113, 
1114, 
1115, 
1116, 
1117, 
1118, 
1119, 
111, 
1120, 
1121, 
1122, 
1123, 
1124, 
1125, 
1126, 
1127, 
1128, 
1129, 
112, 
1130, 
1131, 
1132, 
1133, 
1134, 
1135, 
1136, 
1137, 
1138, 
1139, 
113, 
1140, 
1141, 
1142, 
1143, 
1144, 
1145, 
1146, 
1147, 
1148, 
1149, 
114, 
1150, 
1151, 
1152, 
1153, 
1154, 
1155, 
1156, 
1157, 
1158, 
1159, 
115, 
1160, 
1161, 
1162, 
1163, 
1164, 
1165, 
1166, 
1167, 
1168, 
1169, 
116, 
1170, 
1171, 
1172, 
1173, 
1174, 
1175, 
1176, 
1177, 
1178, 
1179, 
117, 
1180, 
1181, 
1182, 
1183, 
1184, 
1185, 
1186, 
1187, 
1188, 
1189, 
118, 
1190, 
1191, 
1192, 
1193, 
1194, 
1195, 
1196, 
1197, 
1198, 
1199, 
119, 
11
]


usage = "usage: python getDataset.py"
parser = OptionParser(usage)

parser.add_option("-d", "--filelist", default="fileList/mc.list", type="string", help="filelist", dest="filelist")
parser.add_option("-c", "--chunk", default=1, type="int", help="chunk", dest="chunk")
parser.add_option("-j", "--jdir", default="job", type="string", help="job output dir", dest="jdir")
parser.add_option("-o", "--odir", default="job", type="string", help="output dir", dest="odir")
parser.add_option("-n", "--name", default="None", type="string", help="name", dest="name")
parser.add_option("-t", "--type", default="mc", type="string", help="type", dest="type")

(options, args) = parser.parse_args()


listoffiles = []

for line in open(options.filelist):

    if line.find('.root')==-1: continue

    line = line.rstrip()

    listoffiles.append(line)

#print listoffiles

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

    if options.type=='data':
        if ijob in test: continue

#    if ijob > 1000: break

    
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
    outfile = outdir + '/Myroot_' + str(ijob) + '.root'
    edmoutfile = outdir + '/EDM_' + str(ijob) + '.root'

    os.system("cp job_template.sh " + jobscript)
    
    with open(jobscript) as f:
        data_lines = f.read()
        
    inputcommand = ''

    for ff in filename:
        inputcommand += 'inputFiles=' + ff + ' '


    data_lines = data_lines.replace('INPUT', inputcommand).replace('TMPOUTPUT', 'Myroot_' + str(ijob) + '_tmp.root').replace('FINALOUTPUT', 'Myroot_' + str(ijob) + '.root').replace('EDMOUTFILE', edmoutfile).replace('OUTFILE', outfile).replace('RTYPE', options.type)
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)

#    command = 'sbatch -p short --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
    command = 'sbatch -p standard --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
    os.system(command)
