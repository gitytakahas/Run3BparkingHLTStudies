from ROOT import TDatime, TGraph, TFile, TH1F

from optparse import OptionParser, OptionValueError
usage = "usage: python estimate.py"
parser = OptionParser(usage)

parser.add_option('-c', '--create', action="store_true", default=False, dest='create')
parser.add_option('-p', '--pseudo', action="store_true", default=False, dest='pseudo')
parser.add_option("-t", "--l1tol", default=80000, type="float", help="max. l1 tolerance", dest="l1tol")

(options, args) = parser.parse_args()

print(options)


import numpy as np

lumi_level = 20160

l1_ptrange = np.arange(5, 11., 0.5).tolist() 
hlt_ptrange = np.arange(4, 11., 0.5).tolist() 
lumi_range = np.arange(0, lumi_level, 20).tolist() 



drdict = {
    4.0:1.0,
    4.5:0.9,
    5.0:0.9,
    5.5:0.8,
    6.0:0.8,
    6.5:0.8,
    7.0:0.8,
    7.5:0.7,
    8.0:0.7,
    8.5:0.7,
    9.0:0.7,
    9.5:0.6,
    10.0:0.6,
    10.5:0.6,
    11.0:0.6,
    11.5:0.5,
    12.0:0.5,
    12.5:0.5,
    13.0:0.5,
    13.5:0.4,
    14.0:0.4,
}   


# http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
# truely inclusive cross-section
#Sigma_B = 1.5250e+08 # pb
fB = 0.4
Sigma_B = 4.6940e+08 # fB*pb
Br_kee = 4.7e-7

print('Cross-section', Sigma_B)
print('Branching ratio=', Br_kee)

#lumi=10 # [nb-1 s-1]
#eff=2e-4

#rate=lumi*Sigma_B*pow(10.,-3)*Br_kee*eff

#print(lumi,Sigma_B,pow(10.,-3),Br_kee,eff)

#print('rate=',rate)


if options.create:

    ##https://cmsoms.cern.ch/cms/runs/report?cms_run=324980&cms_run_sequence=GLOBAL-RUN
    ##"324980": [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]],
    ##golden = [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]
    golden = [[53, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]

    times = []
    lumis = []

    for line in open('/work/ytakahas/work/Trigger/CMSSW_11_1_0/src/Analysis/HLTAnalyserPy/test/LumiData_2018_20200401.csv', 'r'):

        if line.find('324980:7321')==-1: continue

        line = line.rstrip().split(',')
        ls = line[1].split(':')[0]


        flag = False

        for lrange in golden:
            if int(ls) >= min(lrange) and int(ls) <= max(lrange):
                flag = True
    

        if not flag: continue

        time = line[2].split(' ')
    
        instL = float(line[5])*0.0001

        time_ = TDatime(int(time[0].split('/')[2])+2000, int(time[0].split('/')[0]), int(time[0].split('/')[1]), int(time[1].split(':')[0]), int(time[1].split(':')[1]), int(time[1].split(':')[2]))

        times.append(time_.Convert())
        lumis.append(instL)




    min_times = min(times)

    times = [number - min_times for number in times]
    

    graph = TGraph()
    graph.SetName('lumiprofile')
    graph.SetTitle('lumiprofile')

#    graph_inv = TGraph()
#    graph_inv.SetName('lumiprofile_inv')
#    graph_inv.SetTitle('lumiprofile_inv')

    idx = 0
    for time, lumi in zip(times, lumis):
        graph.SetPoint(idx, time, lumi)
#        graph_inv.SetPoint(idx, lumi, time)

        idx += 1



    file = TFile('estimate.root', 'recreate')
    graph.Write()
#    graph_inv.Write()
    file.Write()
    file.Close()







if options.pseudo:

    ##https://cmsoms.cern.ch/cms/runs/report?cms_run=324980&cms_run_sequence=GLOBAL-RUN
    ##"324980": [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]],
    ##golden = [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]
    golden = [[53, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]

    times = []
    lumis = []


    for line in open('/work/ytakahas/work/Trigger/CMSSW_11_1_0/src/Analysis/HLTAnalyserPy/test/LumiData_2018_20200401.csv', 'r'):

        if line.find('324980:7321')==-1: continue

        line = line.rstrip().split(',')
        ls = line[1].split(':')[0]


        flag = False

        for lrange in golden:
            if int(ls) >= min(lrange) and int(ls) <= max(lrange):
                flag = True
    

        if not flag: continue

        time = line[2].split(' ')
    
        instL = float(line[5])*0.0001

        time_ = TDatime(int(time[0].split('/')[2])+2000, int(time[0].split('/')[0]), int(time[0].split('/')[1]), int(time[1].split(':')[0]), int(time[1].split(':')[1]), int(time[1].split(':')[2]))

        times.append(time_.Convert())
        lumis.append(instL)



    graph = TGraph()
    graph.SetName('lumiprofile')
    graph.SetTitle('lumiprofile')


    min_times = min(times)
    idx = 0

    for ii in lumi_range:
        graph.SetPoint(idx, ii, 2.)

        idx+=1

    times = [number - min_times + lumi_level  for number in times]
    


#    graph_inv = TGraph()
#    graph_inv.SetName('lumiprofile_inv')
#    graph_inv.SetTitle('lumiprofile_inv')

    for time, lumi in zip(times, lumis):
        graph.SetPoint(idx, time, lumi + 0.2)

        if time > 39600:
            break


#        graph_inv.SetPoint(idx, lumi, time)

        idx += 1



    file = TFile('estimate_pseudo.root', 'recreate')
    graph.Write()
#    graph_inv.Write()
    file.Write()
    file.Close()






file2read = 'estimate.root'

if options.pseudo:
    file2read = 'estimate_pseudo.root'
    
file = TFile(file2read)


graph = file.Get('lumiprofile')
#graph_inv = file.Get('lumiprofile_inv')

l1_file = TFile('/work/ytakahas/work/Trigger/CMSSW_11_1_0/src/Run3BparkingTriggerStudies/bandwidth.root')
l1rate = l1_file.Get('otherrate')

hlt_file = TFile('rob/Run3BparkingHLTStudies/roc.root')

total = 0
total_prompt = 0


 


#2E+34,1.7E+34,1.5E+34,1.3E+34,1.1E+34,9E+33,6E+33

switch_lumi = [(2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.)]

#menudict = {}

#for sl in switch_lumi:
#    menudict[sl] = graph_inv.Eval(sl)
#}

max_bw = options.l1tol
print('max_l1 tolerance=', max_bw)

max_bw_hlt = {
    2.0:1515,
    1.7:1740,
    1.5:1929,
    1.3:2163.5,
    1.1:2463,
    0.9:2929
}
    
    

#print(menudict)

print('n, xmin, xmax = ', graph.GetN(), graph.GetPointX(0), graph.GetPointX(graph.GetN()-1))

#h_integral = TH1F('n', 'n', graph.GetN(), graph.GetPointX(0), graph.GetPointX(graph.GetN()-1))
#h_prompt_integral = TH1F('n_prompt', 'n_prompt', graph.GetN(), graph.GetPointX(0), graph.GetPointX(graph.GetN()-1))

h_integral = TGraph()
h_prompt_integral = TGraph()

h_integral.SetName('n')
h_integral.SetTitle('n')

h_prompt_integral.SetName('n_prompt')
h_prompt_integral.SetTitle('n_prompt')

total_lumi = 0
idx_total = 0
idx_total_prompt = 0

for ii in range(graph.GetN()):

    if ii==0: continue

    instL = graph.GetPointY(ii)
    time_duration = graph.GetPointX(ii) - graph.GetPointX(ii-1)

#    print(instL, time_duration)

    total_lumi += instL*10*time_duration
    
    which_lumi = -1

    for sl in switch_lumi:

#        print(sl[0], sl[1])
        if sl[1] <= instL and instL <= sl[0]:
            which_lumi = sl[0]
    
#    print(which_lumi)

    # determine the L1 threshold ... the maximum one that can accommodated with L1 rate ... 

    if which_lumi==-1: 
        print('WARNING!!: no corresponding lumis!!!')
        continue

#    import pdb; pdb.set_trace()

    parking_bw = max_bw - l1rate.Eval(which_lumi)

    # scan which trigger will give us the best shot ... 


    roc_prompt = hlt_file.Get('inv_pt10p0')
    eff_prompt = roc_prompt.Eval(100)

    print(ii, 'eff_prompt, instL', eff_prompt, instL)

    n_prompt = instL*10.*fB*Sigma_B*pow(10.,-3)*Br_kee*eff_prompt*time_duration
    total_prompt += n_prompt
#    h_prompt_integral.SetBinContent(ii+1, total_prompt)
    h_prompt_integral.SetPoint(idx_total_prompt,  graph.GetPointX(ii), total_prompt)
    idx_total_prompt += 1

    which_pt_trigger = -1

    diff = 1000000

    for pt in l1_ptrange:
        rate_ = l1_file.Get('doubleE' + str(pt).replace('.','p') + '_dR' + str(drdict[pt]).replace('.','p') + '0_fit')

#        import pdb; pdb.set_trace()
#        print(pt, 'doubleE' + str(pt).replace('.','p') + '_dR' + str(drdict[pt]).replace('.','p') + '0_fit')
        if rate_.Eval(which_lumi) < parking_bw:
#        if abs(rate_.Eval(which_lumi) - parking_bw) < diff:
#            diff = abs(rate_.Eval(which_lumi) - parking_bw)
            which_pt_trigger = pt
            break


    if which_pt_trigger==-1: 
        print('L1 not available ... !')
        continue

        

    # now check the HLT curve ... 
    roc = hlt_file.Get('inv_pt' + str(which_pt_trigger).replace('.','p'))
    
    # scan the point until it fits the budget ... 

    eff = -1


#    wheretolook = min(which_pt_trigger - 6, 4)
#    wheretolook /= 0.5
#    print(which_pt_trigger, 6, wheretolook)

#    hltrate = roc.GetPointX(int(wheretolook))

#    if hltrate > max_bw_hlt[which_lumi]:
#        print('This HLT rate is way too large !!!!!!!!!!!!!!!!!! (lumi, hltrate, bw)=', which_lumi, hltrate, max_bw_hlt[which_lumi])


    for ip in range(roc.GetN()):
        hltrate = roc.GetPointX(ip)
        if hltrate < max_bw_hlt[which_lumi]:
            eff = roc.GetPointY(ip)
            break
   
    
#    print('instL, lumi, l1 rate, parking_bw, max_bw, max_bw_hlt, eff, which_pt_trigger = ', instL, which_lumi, l1rate.Eval(which_lumi), parking_bw, max_bw, max_bw_hlt[which_lumi], eff, which_pt_trigger)

   
    n = instL*10.*fB*Sigma_B*pow(10.,-3)*Br_kee*eff*time_duration
   
    total += n

    h_integral.SetPoint(idx_total, graph.GetPointX(ii), total)
    idx_total += 1
     



#import pdb; pdb.set_trace()
#(5e6)/53312 = 93.79
# The run 324980 has 

factor = 40./0.528

print('# of Total Kee events =', total, 'with time duration', graph.GetPointX(graph.GetN()-1), 'sec')
print('Assuming 40/fb of data, we will get ', total*factor, 'events')

print('# of prompt Kee events =', total_prompt, 'with time duration', graph.GetPointX(graph.GetN()-1), 'sec')
print('Assuming 40/fb of data, we will get ', total_prompt*factor, 'events')

print('total lumi = ', total_lumi, 'nb')


file2out = 'n_' + str(options.l1tol) + '.root'

if options.pseudo:
    file2out = 'n_' + str(options.l1tol) + '_pseudo.root'

nfile = TFile(file2out, 'recreate')

h_integral.Write()
h_prompt_integral.Write()
graph.Write()
nfile.Write()
nfile.Close()


