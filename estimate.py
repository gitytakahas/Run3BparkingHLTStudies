from ROOT import TDatime, TGraph, TFile, TH1F, Double, TCanvas, TLegend, gROOT, gStyle, TH2F
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
import numpy as np
import os, sys, copy

#gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)


######### input parameters

def returnPoint(graph, ii):
    
    xx = Double(0.)
    yy = Double(0.)
    
    graph.GetPoint(ii, xx, yy)

    return xx, yy

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


ensureDir('plots/')

lumi_level = 20160 # 6h*3600s/h
# http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
# truely inclusive cross-section
fB = 0.4
Sigma_B = 4.6940e+08 # pb
Br_kee = 4.7e-7

#rate=lumi*Sigma_B*pow(10.,-3)*Br_kee*eff

prompt_hlt_bw = 100.
integrated_lumi = 25.
#############################


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



# create Luminosity profile ... 

lumiprofile2output = 'root/lumiprof.root'

if not os.path.exists(lumiprofile2output):

    ##https://cmsoms.cern.ch/cms/runs/report?cms_run=324980&cms_run_sequence=GLOBAL-RUN
    ##"324980": [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]],
    golden = [[53, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]

    times_ = []
    lumis = []

    for line in open('LumiData/LumiData_2018_20200401.csv', 'r'):

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

        times_.append(time_.Convert())
        lumis.append(instL)




    min_times = min(times_)

    times = [number - min_times for number in times_]
    

    graph = TGraph()
    graph.SetName('lumiprofile')
    graph.SetTitle('lumiprofile')

    idx = 0
    for time, lumi in zip(times, lumis):
        graph.SetPoint(idx, time, lumi)

        idx += 1




    graph_ll = TGraph()
    graph_ll.SetName('lumiprofile_ll')
    graph_ll.SetTitle('lumiprofile_ll')


    idx = 0

    for ii in lumi_range:
        graph_ll.SetPoint(idx, ii, 2.)

        idx+=1

    times = [number - min_times + lumi_level  for number in times_]


    for time, lumi in zip(times, lumis):
        graph_ll.SetPoint(idx, time, lumi + 0.2)

        if time > 39600:
            break

        idx += 1



    file = TFile(lumiprofile2output, 'recreate')
    graph.Write()
    graph_ll.Write()
    file.Write()
    file.Close()



    
file = TFile(lumiprofile2output)

graph_norm = file.Get('lumiprofile')
graph_ll = file.Get('lumiprofile_ll')


# This will be replaced once Sebastian derive realistic L1/HLT rates
l1_file = TFile('root/bandwidth.root')
l1rate = l1_file.Get('otherrate')

#hlt_file = TFile('rob/Run3BparkingHLTStudies/roc.root')

total = 0
total_prompt = 0
 
# PS column
# Index,Name,Emergency,2.2E+34 ,2E+34, 1.7E+34, 1.5E+34, 1.3E+34, 1.1E+34, 9E+33, 6E+33

# conversion
# nPU = (instL+0.0011904)/0.0357388
# instL = nPU*0.0357338 - 0.0011904

switch_lumi = [(2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.)]
switch_npu = [48, 42, 36, 30, 25, 17]

max_bw_hlt = {
    2.0:1515,
    1.7:1740,
    1.5:1929,
    1.3:2163.5,
    1.1:2463,
    0.9:2929
}
    


#    for graph, name in zip([graph_norm, graph_ll], ['norm', 'll']):

for graph, name in zip([graph_norm], ['norm']):

    hists = []

    #for l1tol in [80, 90]:
    for l1tol in [80000]:


        h_integral = TGraph()
        h_integral.SetName('n_l1tol' + str(l1tol) + '_' +  name)
        h_integral.SetTitle('n_l1tol' + str(l1tol) + '_' + name)

        h_prompt_integral = TGraph()
        h_prompt_integral.SetName('n_prompt_l1tol' + str(l1tol) + name)
        h_prompt_integral.SetTitle('n_prompt_l1tol' + str(l1tol) + name)


        total_lumi = 0
        idx_total = 0
        idx_total_prompt = 0

        for ii in range(graph.GetN()):

            if ii==0: continue

#            import pdb; pdb.set_trace()

            t1, instL = returnPoint(graph, ii)
            t2, instL_ = returnPoint(graph, ii-1)

#            instL = graph.GetPointY(ii)
            time_duration = t1 - t2 #graph.GetPointX(ii) - graph.GetPointX(ii-1)

            total_lumi += instL*10*time_duration
            
            which_lumi = -1
            which_npu = -1

            for index, sl in enumerate(switch_lumi):

                if sl[1] <= instL and instL <= sl[0]:
                    which_lumi = sl[1]
                    which_npu = switch_npu[index]
    

            if which_lumi==-1: 
                print('WARNING!!: no corresponding lumis!!!')
                continue


            parking_bw = l1tol - l1rate.Eval(which_lumi)

            print(instL, l1tol, which_lumi, which_npu, l1rate.Eval(which_lumi), parking_bw)

            hlt_file = TFile('root/roc_hlt_pu' + str(which_npu) + '.root')


            which_pt_trigger = 10.5
            flag_park = False

            for pt in l1_ptrange:

                rate_ = l1_file.Get('doubleE' + str(pt) + ', dR < ' + str(drdict[pt]) + '0')
                
#                print('check=', 'doubleE' + str(pt).replace('.','p') + '_dR' + str(drdict[pt]).replace('.','p') + '0_fit', rate_)

                if rate_.Eval(which_lumi) <= parking_bw:
                    which_pt_trigger = pt
                    flag_park = True
                    break



#            print(which_lumi, which_npu, which_pt_trigger, parking_bw)





#            rochlt = hlt_file.Get('rep_l1pt' + str(which_pt_trigger).replace('.','p') + '_hltpt' + str(which_pt_trigger - 1).replace('.','p'))
                




            roc_prompt = hlt_file.Get('inv_pt' + str(which_pt_trigger).replace('.','p'))
            eff_prompt = -1.

            for ip in range(roc_prompt.GetN()):
                hltrate, eff = returnPoint(roc_prompt, ip)
                if hltrate < prompt_hlt_bw:
                    eff_prompt = eff
                    break




###            diff = 10000. 
###            eff_prompt = -1.
###            rate_prompt = -1.
###
###            for ip in range(roc_prompt.GetN()):
###                hltrate, dummy = returnPoint(roc_prompt, ip)
###                
###                diff_ = abs(hltrate - prompt_hlt_bw)
###
###                if diff_ < diff:
###                    diff = diff_
###                    dummy, eff_prompt = returnPoint(roc_prompt, ip)
###                    rate_prompt = hltrate
####                    break



            if eff_prompt==-1:
                print('WARNING!!: this cannot happen!')



            n_prompt = instL*10.*fB*Sigma_B*pow(10.,-3)*Br_kee*eff_prompt*time_duration
            total_prompt += n_prompt
            h_prompt_integral.SetPoint(idx_total_prompt,  t1, total_prompt)
#            print(idx_total_prompt,  t1, total_prompt)


            idx_total_prompt += 1




            if not flag_park: continue



#                print('check = rep_l1pt' + str(which_pt_trigger).replace('.','p') + '_hltpt' + str(which_pt_trigger - 1).replace('.','p'))

            rochlt = hlt_file.Get('rep_l1pt' + str(which_pt_trigger).replace('.','p') + '_hltpt' + str(which_pt_trigger - 1).replace('.','p'))
            eff_hlt, dummy = returnPoint(rochlt, 0)    

         
#            for ip in range(roc.GetN()):
#                hltrate = roc.GetPointX(ip)
#                if hltrate < max_bw_hlt[which_lumi]:
#                    eff_hlt = roc.GetPointY(ip)
#                    break
  
   
            n = instL*10.*fB*Sigma_B*pow(10.,-3)*Br_kee*eff_hlt*time_duration
   
            total += n

            h_integral.SetPoint(idx_total, t1, total)
            idx_total += 1

#                print(idx_total, t1, total)


        print(name, '--> l1 total rate = ', l1tol)
 
        factor = integrated_lumi/(total_lumi*pow(10.,-6))

        totalt, dummy = returnPoint(graph, graph.GetN()-1)

        print('# of Total Kee events =', total, 'with time duration', totalt, 'sec')
        print('Assuming ' + str(integrated_lumi) + '/fb of data, we will get', total*factor, 'events')
        print('factor = ', integrated_lumi, total_lumi*pow(10.,-6), factor)
        
        print('# of prompt Kee events =', total_prompt, 'with time duration', totalt, 'sec')
        print('Assuming ' + str(integrated_lumi) + '/fb of data, we will get ', total_prompt*factor, 'events')
        
        print('total lumi = ', total_lumi, 'nb^-1')
        
        hists.append(copy.deepcopy(h_prompt_integral))
        hists.append(copy.deepcopy(h_integral))
        

        ##### drawing ... 


#    import pdb; pdb.set_trace()
    canvas = TCanvas('can_' +  name)
    canvas.SetGridx()
    canvas.SetGridy()
    
    leg = TLegend(0.2, 0.68,0.5,0.86)
    
    applyLegendSettings(leg)
    leg.SetTextSize(0.03)
    
    ymax = max([hist.GetMaximum() for hist in hists])

    frame = TH2F('frame_' + name, 'frame_' + name, graph.GetN(), 0, totalt, 100,0,100)
    frame.GetXaxis().SetTitle('Time (s)')
    frame.GetYaxis().SetTitle('Instantaneous Luminosity (E34)')
    frame.Draw()


    for ii, hist in enumerate(hists):

        hist.SetMaximum(4)
        hist.SetLineColor(ii+2)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(ii+1)
        hist.SetMarkerSize(1)

        hist.Draw('plsame')

        leg.AddEntry(hist, hist.GetName().replace('bw_', 'L1 max = '), 'lep')


    graph.SetLineWidth(3)
    graph.Draw('same')

#    prompt.SetLineColor(4)
#    prompt.SetLineWidth(2)
#    prompt.Draw('same')

    leg.AddEntry(graph, 'Fill7321 (0.55/fb)', 'lep')
#    leg.AddEntry(prompt, 'Prompt (5kHz at L1, 100Hz at HLT)', 'lep')

    leg.Draw()

    canvas.RedrawAxis();

#    filename = 'plots/estimate.gif'

#    if options.pseudo:
#        filename = 'plots/estimate_pseudo.gif'

    filename = 'plots/' + name + '.gif'

    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")
    
    canvas.SaveAs(filename)
    canvas.SaveAs(filename.replace('.gif', '.pdf'))








#        file2out = 'n_' + str(options.l1tol) + '.root'
#
#        if options.pseudo:
#            file2out = 'n_' + str(options.l1tol) + '_pseudo.root'
#
#        nfile = TFile(file2out, 'recreate')
#
#        h_integral.Write()
#        h_prompt_integral.Write()
#        graph.Write()
#        nfile.Write()
#        nfile.Close()
