from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, kGray
from ctypes import c_double
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
import numpy as np
import os, sys, copy
from common import path

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

lumi_level = 20160 # 6h*3600s/h

# http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
# truely inclusive cross-section
Sigma_B = 4.6940e+11 # fb
fB = 0.4
Br_kee = 2*4.5e-7

integrated_lumi = 25.

lumi_range = np.arange(0, lumi_level, 20).tolist() 

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
    
################################################################################
# Determine various metrics and construct menu

# Level-1 total rate estimate for full CMS menu
l1_file = TFile(path+'ee/l1_bandwidth.root')
l1rate = l1_file.Get('otherrate')

# Di-electron trigger rate estimates from data
l1_file_official = TFile(path+'ee/l1_bandwidth_official.root')

# conversion
# nPU = (instL+0.0011904)/0.0357388
# instL = nPU*0.0357338 - 0.0011904
switch_lumi = [(2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), (0.6, 0.1)]
switch_npu = [56, 48, 42, 36, 30, 25, 17]

# HLT bandwidth from Sara's presentation
# https://indico.cern.ch/event/1032638/contributions/4336416/
max_bw_hlt = {
    2.0:1515,
    1.7:1740,
    1.5:1929,
    1.3:2163.5,
    1.1:2463,
    0.9:2929,
    0.6:3791
}

################################################################################
# Initialise
# For each L1 pT threshold: 
#   For each Peak Linst:
#     Calc spare capacity
#     Get L1 rate, Get HLT rate
#     Check if within capacity

l1_pts = np.arange(4.,11.,0.5).tolist()
hlt_pts = np.arange(4.,11.,0.5).tolist()
print(l1_pts)
print(hlt_pts)
print

def print_table(debug,value,l1_total=90000.,dimuon=5000.,allocation=5000.):
    print("Table:",value)
    dct = {}
    hlt_files = [ TFile('root/roc_hlt_pu'+str(npu)+'.root') for npu in switch_npu ] 
    for ii,l1_pt in enumerate(l1_pts) :
        histo_l1_rate = l1_file_official.Get('L1_DoubleEG'+
                                             str(l1_pt).replace('.','p').replace('p0','')+
                                             'er1p22_dR_' + str(drdict[l1_pt]).replace('.','p'))
        if debug : print("L1 LOOP:",ii,l1_pt,histo_l1_rate) # DEBUG
        if not debug and value is not "menu" : print("{:4.1f}".format(l1_pt),end="")
        l1rate_correction = max(0., l1rate.Eval(switch_lumi[0][0]) - l1_total) # (Corrected down by ~7 kHz)
        for jj,(lumi,npu,hlt_file) in enumerate(reversed(zip(switch_lumi,switch_npu,hlt_files))) :
            peak_lumi = lumi[0]
            spare = max(0.,l1_total - l1rate.Eval(peak_lumi) + l1rate_correction)
            spare -= dimuon # Account for allocation to di-muon trigger
            l1_rate = histo_l1_rate.Eval(npu)*1000
            l1_ok = l1_rate < allocation or (l1_rate-allocation) < spare
            hlt_roc = hlt_file.Get('inv_pt' + str(l1_pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()
            hlt_max = max_bw_hlt[peak_lumi]
            hlt_pt = -1
            hlt_rate = -1
            hlt_eff = -1
            hlt_ok = True
            for kk in range(hlt_n):
                ip = kk#hlt_n - kk - 1
                hlt_pt = hlt_pts[ip]
                hlt_rate = hlt_roc.GetPointX(ip)
                if hlt_rate < hlt_max:
                    hlt_eff = hlt_roc.GetPointY(ip)
                    if debug : print("    HLT SCAN:",hlt_n,ip,hlt_pt,hlt_max,hlt_rate,hlt_eff)
                    break
            if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')

            if peak_lumi not in dct and l1_ok and hlt_ok : 
                dct[peak_lumi] = (l1_pt,hlt_pt,l1_rate,hlt_rate,hlt_eff)

            if not debug : 
                #l1_ok = True; hlt_ok = True #@@
                if value == "l1_rate" : 
                    string = " & {:5.1f}".format(l1_rate/1000.)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "spare" : 
                    string = " & {:5.1f}".format(spare/1000.)
                    print(string,end="")
                elif value == "l1_pt" : 
                    string = " & {:5.1f}".format(l1_pt)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_pt" : 
                    string = " & {:5.1f}".format(hlt_pt)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_rate" : 
                    string = " & {:5.2f}".format(hlt_rate/1000.)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "eff" : 
                    string = " & {:5.2f}".format(hlt_eff*1.e4)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "lint" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12. # Assume 12-hour fill
                    string = " & {:5.2f}".format(Lint)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fill" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12.
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fb" : 
                    Lint = integrated_lumi # defined above (25.?)
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if l1_rate > l1_total : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "menu" : 
                    pass
                else : 
                    print("???")
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(l1_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG
        if not debug and value is not "menu" : print("\\\\")
    if value == "menu" :
        Lint = integrated_lumi # defined above
        print("Lint:", Lint)
        for jj,lumi in enumerate(switch_lumi) :
            peak_lumi = lumi[0]
            if peak_lumi not in dct : 
                print(peak_lumi," \\\\")
                continue
            (l1_pt,hlt_pt,l1_rate,hlt_rate,hlt_eff) = dct[peak_lumi]
            count = Lint * fB * Sigma_B * Br_kee * hlt_eff
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12.
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff
            print(peak_lumi," & ",\
                "{:5.1f}".format(l1_pt)," & ",\
                "{:5.1f}".format(hlt_pt)," & & ",\
                "{:5.1f}".format(l1_rate/1000.)," & ",\
                "{:5.2f}".format(hlt_rate/1000.)," & ",\
                "{:5.2f}".format(hlt_eff/1.e-4)," & ",\
                "{:5.2f}".format(Lint_per_fill)," & ",\
                "{:5.2f}".format(count_per_fill)," & ",\
                "{:5.1f}".format(count)," \\\\ ")
    print()

for value in ["l1_rate",
              "spare",
              "l1_pt",
              "hlt_pt",
              "hlt_rate",
              "eff",
              "lint",
              "counts_per_fill",
              "counts_per_fb",
]: print_table(debug=False,value=value,l1_total=90000.,dimuon=5000.,allocation=5000.)

l1_total=90000.
for allocation in [0.,5000.,10000.,20000.]: 
    print("allocation:",allocation)
    print_table(debug=False,value="menu",l1_total=l1_total,dimuon=5000.,allocation=allocation)
