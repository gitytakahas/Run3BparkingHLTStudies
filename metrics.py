from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, kGray
from ctypes import c_double
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
import numpy as np
import os, sys, copy
from common import *
import json

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-c', '--corrected', action="store_true", default=False, dest='corrected', help="apply trigger rates correction factor")
parser.add_option('-l', '--limit', action="store_true", default=False, dest='limit', help="limit HLT rate to 300 Hz")
(options, args) = parser.parse_args()

################################################################################
# input parameters

ensureDir('plots/')

# Level-1 total rate estimate for full CMS menu
l1_file = TFile(path+'ee/l1_bandwidth.root')
l1_his = l1_file.Get('otherrate')

# Di-electron trigger rate estimates from data
ee_file = TFile(path+'ee/l1_bandwidth_official.root')

integrated_lumi = 25.

l1_pts = np.arange(4.,11.,0.5).tolist()
hlt_pts = np.arange(4.,11.,0.5).tolist()

# conversion
# nPU = (Linst+0.0011904)/0.0357388
# Linst = nPU*0.0357338 - 0.0011904
switch_lumi = [(2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), (0.6, 0.)]
switch_npu = [56, 48, 42, 36, 30, 25, 17]

l1_max = 95000.
dimuon = 0.

idx_2e34 = 0
l1_rate_corr = l1_his.Eval(switch_lumi[idx_2e34][0]) - l1_max

corrs_dict = {}
if options.corrected:
    for npu in [56, 48, 42, 36, 30, 25, 17]:
        filename = path+'rates/corrections_' + str(npu) + '.json'
        infile = open(filename,'r')
        dct = json.load(infile)
        corrs_dict[npu] = dct
        infile.close()

################################################################################
# Initialise
# For each L1 pT threshold: 
#   For each Peak Linst:
#     Calc spare capacity
#     Get L1 rate, Get HLT rate
#     Check if within capacity

def print_table(debug,value,l1_max=95000.,dimuon=0.,allocation=5000.):
    print("Table:",value)
    dct = {}
    hlt_files = [ TFile(path+'ee/roc_hlt_pu'+str(npu)+'.root') for npu in switch_npu ] 
    for ii,l1_pt in enumerate(l1_pts) :
        histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                    str(l1_pt).replace('.','p').replace('p0','')+
                                    'er1p22_dR_' + str(drdict[l1_pt]).replace('.','p'))
        if debug : print("L1 LOOP:",ii,l1_pt,histo_ee_rate) # DEBUG
        if not debug and value != "menu" : print("{:4.1f}".format(l1_pt),end="")
        for jj,(lumi,npu,hlt_file) in enumerate(reversed(list(zip(switch_lumi,switch_npu,hlt_files)))) :
            peak_lumi = lumi[0]
            l1_rate = l1_his.Eval(peak_lumi)
            spare = l1_max+l1_rate_corr - l1_rate
            #spare -= dimuon # Account for allocation to di-muon trigger
            ee_rate = histo_ee_rate.Eval(npu)*1000
            l1_ok = ee_rate < allocation or ee_rate < (spare+allocation)

            hlt_roc = hlt_file.Get('inv_pt' + str(l1_pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

            # Extract corrections to rates from JSON
            if options.limit:
                hlt_max = 300. #max_bw_hlt[peak_lumi] # <-- Limit to 300 Hz for prompt reco
                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                hlt_pt = -1
                for kk in range(hlt_n):
                    ip = kk#hlt_n - kk - 1
                    # print("test",kk,ip,hlt_pts[ip],hlt_roc.GetPointX(ip),hlt_max,hlt_roc.GetPointY(ip))
                    hlt_rate = hlt_roc.GetPointX(ip)
                    hlt_eff = hlt_roc.GetPointY(ip)
                    hlt_pt = hlt_pts[ip]
                    if options.corrected:
                        l1ptstr = str(l1_pt)
                        hltptstr = str(hlt_pt)
                        if npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[npu].keys() and \
                                hltptstr in corrs_dict[npu][l1ptstr].keys():
                            corr = corrs_dict[npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",npu,l1ptstr,hltptstr)
                    if hlt_rate < hlt_max:
                        if debug : print("    HLT SCAN:",hlt_n,ip,hlt_pt,hlt_max,hlt_rate,hlt_eff)
                        break
                if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')
            else:
                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                hlt_pt = hlt_threshold_dict.get(l1_pt,4.0)
                hltpt_list = np.arange(4, 11, 0.5).tolist()
                index = hltpt_list.index(hlt_pt)
                if index < hlt_roc.GetN():
                    hlt_rate = hlt_roc.GetPointX(index)
                    hlt_eff = hlt_roc.GetPointY(index)
                    if options.corrected:
                        l1ptstr = str(l1_pt)
                        hltptstr = str(hlt_pt)
                        if npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[npu].keys() and \
                                hltptstr in  corrs_dict[npu][l1ptstr].keys():
                            corr = corrs_dict[npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",npu,l1ptstr,hltptstr)
                else: print("cannot find hltpt")

            #print("test",ii,jj,npu,l1_pt,hlt_pt,hlt_eff,hlt_rate)

            if peak_lumi not in dct and l1_ok and hlt_ok : 
                dct[peak_lumi] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff)

            if not debug : 
                #l1_ok = True; hlt_ok = True #@@
                if value == "ee_rate" : 
                    string = " & {:5.1f}".format(ee_rate/1000.)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "spare" : 
                    string = " & {:5.1f}".format(spare/1000.)
                    print(string,end="")
                elif value == "l1_pt" : 
                    string = " & {:5.1f}".format(l1_pt)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_pt" : 
                    string = " & {:5.1f}".format(hlt_pt)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_rate" : 
                    string = " & {:5.2f}".format(hlt_rate/1000.)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "eff" : 
                    string = " & {:5.2f}".format(hlt_eff*1.e4)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "lint" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12. # Assume 12-hour fill
                    string = " & {:5.2f}".format(Lint)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fill" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12.
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fb" : 
                    Lint = integrated_lumi # defined above (25.?)
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "menu" : 
                    pass
                else : 
                    print("???")
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(ee_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG
        if not debug and value != "menu" : print("\\\\")
    if value == "menu" :
        Lint = integrated_lumi # defined above
        print("Lint:", Lint)
        for jj,lumi in enumerate(switch_lumi) :
            peak_lumi = lumi[0]
            if peak_lumi not in dct : 
                print(peak_lumi," \\\\")
                continue
            (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff) = dct[peak_lumi]
            count = Lint * fB * Sigma_B * Br_kee * hlt_eff
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12.
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff
            print(peak_lumi," & ",\
                "{:5.1f}".format(l1_pt)," & ",\
                "{:5.1f}".format(hlt_pt)," & & ",\
                "{:5.1f}".format(ee_rate/1000.)," & ",\
                "{:5.2f}".format(hlt_rate/1000.)," & ",\
                "{:5.2f}".format(hlt_eff/1.e-4)," & ",\
                "{:5.2f}".format(Lint_per_fill)," & ",\
                "{:5.2f}".format(count_per_fill)," & ",\
                "{:5.1f}".format(count)," \\\\ ")
    print()

for value in ["ee_rate",
              "spare",
              "l1_pt",
              "hlt_pt",
              "hlt_rate",
              "eff",
              "lint",
              "counts_per_fill",
              "counts_per_fb",
]: print_table(debug=False,value=value,l1_max=l1_max,dimuon=dimuon,allocation=5000.)

for allocation in [0.,5000.,10000.,20000.]: 
#for allocation in [0.,3600.,6100.,18900.]: 
    print("allocation:",allocation)
    print_table(debug=False,value="menu",l1_max=l1_max,dimuon=dimuon,allocation=allocation)
