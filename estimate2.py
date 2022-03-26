from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, kGray
from ctypes import c_double
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
import numpy as np
import os, sys, copy
from common import path
from common import createLumiProfiles

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

################################################################################
# input parameters

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
ensureDir('plots/')

Sigma_B = 4.6940e+11 # fb (http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html)
fB = 0.4
Br_kee = 2*4.5e-7
integrated_lumi = 25.

l1_ptrange = np.arange(4, 11, 0.5).tolist() 

drdict = {
    4.0:0.9,
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

# Level-1 total rate estimate for full CMS menu
l1_file = TFile(path+'ee/l1_bandwidth.root')
l1rate = l1_file.Get('otherrate')

# Di-electron trigger rate estimates from data
l1_file_official = TFile(path+'ee/l1_bandwidth_official.root')

# conversion
# nPU = (Linst+0.0011904)/0.0357388
# Linst = nPU*0.0357338 - 0.0011904
switch_lumi = [(2.2, 2.0), (2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), 
               (0.6, 0.45), (0.45, 0.3), (0.3, 0.15), (0.15, 0.)]
switch_npu = [56, 56, 48, 42, 36, 30, 25,
              17, 17, 17, 17] # NEED TO UPDATE FOR LOWER LINST

# HLT bandwidth from Sara's presentation
# https://indico.cern.ch/event/1032638/contributions/4336416/
max_bw_hlt = {
    2.2:1515,
    2.0:1515,
    1.7:1740,
    1.5:1929,
    1.3:2163.5,
    1.1:2463,
    0.9:2929,
    0.6:3791,
    0.45:3791,0.3:3791,0.15:3791 # NEED TO UPDATE FOR LOWER LINST
}

################################################################################
# Read (same) luminosity profiles (again) ...

output='root/lumiprof.root'
createLumiProfiles(output)    
file = TFile(output)

################################################################################
# 
  
l1_total = 90000.
dimuon = 5000.
l1_total_adj = l1_total - dimuon

summary = []

profiles = [
    file.Get('falling_from_1p8e34_original'),
    file.Get('falling_from_1p8e34'),
    file.Get('falling_from_0p9e34'),
    file.Get('levelled_at_2p0e34'),
    file.Get('levelled_at_1p7e34'),
    file.Get('levelled_at_1p5e34'),
    file.Get('levelled_at_1p3e34'),
    file.Get('levelled_at_1p1e34'),
    file.Get('levelled_at_0p9e34'),
    file.Get('levelled_at_0p6e34'),
    file.Get('levelled_at_4p5e33'),
    file.Get('levelled_at_3p0e33'),
    file.Get('levelled_at_1p5e33')
]

for idx,profile in enumerate(profiles):
    name = profile.GetName()
    print()
    print("name:",name)

    graphs = []
    for colour,allocation in enumerate([5000,10000,20000]):

        graph = TGraph()
        graph.SetName('name'+'_allocation_' + str(allocation))
        graph.SetTitle('L1 allocation: ' + str(int(allocation/100.)/10.) + ' kHz')
        graph.SetLineColor(colour+1)
        graph.SetMarkerSize(0)
        graph.SetLineWidth(3)

        ls_cntr = 0
        total_lumi = 0
        total_count = 0

        old_lumi = -1
        time_elapsed = 0.
        for ii in range(profile.GetN()):

            if ii==0: continue
            Linst = max(profile.GetPointY(ii), 0.1)
            time_duration = profile.GetPointX(ii) - profile.GetPointX(ii-1)
            time_elapsed += time_duration
            total_lumi += Linst*1.e-5 * time_duration

            # Correct down by trigger rate overspend (currently ~7 kHz)
            l1rate_correction = max(0., l1rate.Eval(switch_lumi[0][0]) - l1_total)
            
            # check which luminosity it is ... 
            which_lumi = -1
            which_npu = -1
            for index, sl in enumerate(switch_lumi):

                if sl[1] < Linst and Linst <= sl[0]:
                    #switch = True
                    which_lumi = sl[0]
                    which_npu = switch_npu[index]
                    break

            if which_lumi==-1: 
                print('   WARNING!!: no corresponding lumis!!!',"Linst",Linst)
                continue

            switch = True if which_lumi != old_lumi else False
            #@@print(ii,time_elapsed,Linst,sl[0],sl[1],which_lumi,which_npu,switch,old_lumi)
            old_lumi = which_lumi

            spare = max(0.,l1_total - l1rate.Eval(which_lumi) + l1rate_correction)
            spare -= dimuon # Account for allocation to di-muon trigger

#            print('L=', Linst, 
#                  ', which lumi=', which_lumi, 
#                  '(npu = ', which_npu, 
#                  '), l1 rate =', l1rate.Eval(which_lumi), 
#                  ', l1 b/w =', spare)

            # Determine L1 threshold ... 
            which_l1pt = l1_ptrange[-1] 
            flag_park = False
            l1_ee_rate = None
            for pt in l1_ptrange:
                histo_l1_rate = l1_file_official.Get('L1_DoubleEG'+
                                                     str(pt).replace('.','p').replace('p0','')+
                                                     'er1p22_dR_' + str(drdict[pt]).replace('.','p'))
                l1_rate = histo_l1_rate.Eval(which_npu)*1000
                # vvv HACK HACK HACK vvv
                if   which_lumi == 0.45 : l1_rate *= 0.45/0.6
                elif which_lumi == 0.30 : l1_rate *= 0.30/0.6
                elif which_lumi == 0.15 : l1_rate *= 0.15/0.6
                # ^^^ HACK HACK HACK ^^^
                l1_ok = l1_rate < allocation or (l1_rate-allocation) < spare
                if l1_ok:
                    which_l1pt = pt
                    flag_park = True
                    l1_ee_rate = l1_rate
                    break
            if not flag_park: 
                print('\t L1 not available')
                continue

            #@@print('\t which_l1pt=', which_l1pt, 'with l1 ee rate = ', l1_ee_rate)

            # Determine HLT threshold ...
            hlt_file = TFile('root/roc_hlt_pu' + str(which_npu) + '.root')
            roc = hlt_file.Get('inv_pt' + str(which_l1pt).replace('.','p'))
            eff_hlt = -1
            for ip in range(roc.GetN()):
                if roc.GetPointX(ip) < max_bw_hlt[which_lumi]:
                    eff_hlt = roc.GetPointY(ip)
                    break
            if eff_hlt==-1:
                print('!!!!!! This cannot happen !!!!')
            #@@print('\t parking: HLT eff. =', eff_hlt, 'max b/w=', max_bw_hlt[which_lumi])

            count = Linst * 1.e-5 * fB * Sigma_B * Br_kee * eff_hlt * time_duration
            total_count += count
            graph.SetPoint(ls_cntr, profile.GetPointX(ii), total_count)

            ls_cntr += 1
            if switch:
                print(
                    " ".join(["ii:",str("{:5.0f}".format(ii)),
                              "dt:",str("{:4.1f}".format(time_duration)),
                              "time:",str("{:5.0f}".format(time_elapsed)),
                              "Peak:",str("{:4.2f}".format(which_lumi)),
                              "Linst:",str("{:.3f}".format(Linst)),
                              "Lint:",str("{:.3f}".format(total_lumi)),
                              "spare:",str("{:5.0f}".format(spare)),
                              "pT:",str(which_l1pt),
                              "rate:",str("{:5.0f}".format(l1_ee_rate)),
                              "Eff:",str("{:.6f}".format(eff_hlt)),
                              #"dN:",str("{:.4f}".format(count)),
                              "N:",str("{:6.3f}".format(total_count))
                          ]))
 
        graphs.append(copy.deepcopy(graph))

        print("SUMMARY:") 
        print('Lumi profile name:       ', name) 
        print('L1 total rate:           ', l1_total)
        print('Di-ele allocation:       ', allocation)
        print('Dimuon allocation:       ', dimuon)
        total_time = profile.GetPointX(profile.GetN()-1)
        print('Time duration            ', total_time, 's')
        print('Total lumi:              ', total_lumi,"fb^-1")
        print('Total # or Kee events:   ', total_count)
        factor = integrated_lumi/(total_lumi)#*pow(10.,-6))
        print('Target lumi:             ', integrated_lumi)
        print('Luminosity factor:       ', factor)
        print('Expected # of Kee events:', total_count*factor)
        summary.append((name,allocation,total_lumi,total_count))

    ##########
    # Create canvas
    canvas = TCanvas('canvas_' +  name)
    canvas.SetGridx()
    canvas.SetGridy()

    for only_profile in [True,False]:

        # Maximum for y-axis
        ymax = 2.5 if only_profile else 12.#max([graph.GetPointY(graph.GetN()-1) for graph in graphs])*1.2

        # Create fram for lumi profile
        frame = TH2F('frame_'+name, 'frame_'+name, 
                     profile.GetN(), 0, total_time, # time axis
                     100, 0, ymax ) # counts axis
        frame.GetXaxis().SetTitle('Time (s)')
        if only_profile : frame.GetYaxis().SetTitle('L_{inst} [E34 Hz/cm^{2}]')
        else : frame.GetYaxis().SetTitle('Cumulative Kee event count')
        frame.Draw()

        # Draw lumi profile
        if not only_profile: profile.SetLineColor(kGray)
        profile.SetLineWidth(4)
        profile.SetMarkerSize(0)
        profile.Draw('same')

        # Create legend
        lower = 0.80 if only_profile else 0.80-len(graphs)*0.05
        leg = TLegend(0.2,lower,0.5,0.85)
        applyLegendSettings(leg)
        leg.SetTextSize(0.04)
        leg.AddEntry(profile, 'L_{int}'+' = {:.2f}/fb'.format(total_lumi), 'l')

        # Draw graphs and legend
        if not only_profile:
            for ii, graph in enumerate(graphs[::-1]): # in reverse order
                graph.Draw('lsame')
                leg.AddEntry(graph, graph.GetTitle(), 'l')
        leg.Draw()
        canvas.RedrawAxis();

        # Add CMS labels
        l2=add_Preliminary()
        l2.Draw("same")
        l3=add_CMS()
        l3.Draw("same")

        # Save canvas
        filename = 'plots/'+name+'.pdf'
        if not only_profile: filename = filename.replace("plots/","plots/estimates_for_")
        canvas.SaveAs(filename)

# Print summary
for name,allocation,lumi,count in summary:
    print('Name:',"{:30s}".format(name),
          'L1 alloc [Hz]:',"{:6.0f}".format(allocation),
          'Lint [/fb] per fill:',"{:6.3f}".format(lumi),
          'Counts, per fill:',"{:5.1f}".format(count),
          ', per /fb:',"{:5.1f}".format(count/lumi),
          ', per 25/fb:',"{:5.1f}".format(integrated_lumi*count/lumi)
    )
