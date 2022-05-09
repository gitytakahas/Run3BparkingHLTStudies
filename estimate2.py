from __future__ import print_function
from common import *
import copy, os, sys
from ctypes import c_double
from DisplayManager import DisplayManager, add_Preliminary, add_Private, add_CMS, add_lumi, applyLegendSettings
import numpy as np
from officialStyle import officialStyle
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, TGaxis, gPad, TLine, TPaveText, TLatex, kBlue, kRed, kBlack
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

l1_ptrange = np.arange(4, 11, 0.5).tolist() 

# Conversion:
# nPU = (Linst+0.0011904)/0.0357388
# Linst = nPU*0.0357338 - 0.0011904
switch_lumi = [(2.2, 2.0),
               (2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.7),
               (0.7, 0.47), (0.47, 0.24), (0.24, 0.06), (0.06, 0.)]
switch_npu = [56, # NEED TO UPDATE FOR HIGHEST LINST?
              56, 48, 42, 36, 30, 25,
              17, 17, 17, 17] # NEED TO UPDATE FOR LOWER LINST
  
l1_max = 95000.
#dimuon = 5000.
#l1_max -= dimuon

# Normalise rate to l1_max for 2E34 (i.e. remove any over/underspend)
idx_2e34 = 1
l1_rate_corr = l1_his.Eval(switch_lumi[idx_2e34][0]) - l1_max

summary = []

# Read luminosity profiles ...
output='root/lumiprof.root'
createLumiProfiles(output,synthetic=True)
file_profile = TFile(output)

profiles = [
# FALLING
#    file_profile.Get('falling_from_1p8e34_original'),
#    file_profile.Get('falling_from_1p8e34'),
#    file_profile.Get('falling_from_0p9e34'),
# ORIG
#    file_profile.Get('levelled_at_2p0e34'),
#    file_profile.Get('levelled_at_1p7e34'),
#    file_profile.Get('levelled_at_1p5e34'),
#    file_profile.Get('levelled_at_1p3e34'),
#    file_profile.Get('levelled_at_1p1e34'),
#    file_profile.Get('levelled_at_0p9e34'),
#    file_profile.Get('levelled_at_0p6e34'),
#    file_profile.Get('levelled_at_4p5e33'),
#    file_profile.Get('levelled_at_3p0e33'),
#    file_profile.Get('levelled_at_1p5e33')
# FILIP
#    file_profile.Get('levelled_at_2p1e34'),
#    file_profile.Get('levelled_at_1p8e34'),
#    file_profile.Get('levelled_at_1p4e34'),
#    file_profile.Get('levelled_at_0p9e34'),
#    file_profile.Get('levelled_at_0p7e34'),
#    file_profile.Get('levelled_at_0p4e34'),
#    file_profile.Get('levelled_at_0p2e34'),
# HYBRID
    file_profile.Get('levelled_at_2p0e34'),
    file_profile.Get('levelled_at_1p7e34'),
    file_profile.Get('levelled_at_1p5e34'),
    file_profile.Get('levelled_at_1p3e34'),
    file_profile.Get('levelled_at_1p1e34'),
    file_profile.Get('levelled_at_0p9e34'),
    file_profile.Get('levelled_at_0p7e34'),
    file_profile.Get('levelled_at_0p4e34'),
    file_profile.Get('levelled_at_0p2e34'),
]

corrs_dict = {}
if options.corrected:
    for npu in [56, 48, 42, 36, 30, 25, 17]:
        filename = path+'rates/corrections_' + str(npu) + '.json'
        infile = open(filename,'r')
        dct = json.load(infile)
        corrs_dict[npu] = dct
        infile.close()

################################################################################
# Process ...

for idx,profile in enumerate(profiles):
    name = profile.GetName()
    print()
    print("name:",name)

    # Spare bandwidth
    gspare = TGraph()
    gspare.SetName('spare')
    gspare.SetTitle('L1 spare')
    gspare.SetMarkerSize(0)
    gspare.SetMarkerColor(kRed+1)
    gspare.SetLineStyle(1)
    gspare.SetLineColor(kRed+1)
    gspare.SetLineWidth(3)
    switches=[]

    graphs = []
    allocations = [0,5000,10000,20000] # [3600,6100,18900]
    for style,allocation in enumerate(allocations):

        graph = TGraph()
        graph.SetName('name'+'_allocation_' + str(allocation))
        graph.SetTitle(str(int(allocation/100.)/10.))
        graph.SetMarkerSize(0)
        graph.SetMarkerColor(kBlue+2)
        #graph.SetLineStyle(style+1)
        graph.SetLineStyle(1)
        graph.SetLineColor(kBlue+len(allocations)-style)
        #graph.SetLineWidth(4)
        graph.SetLineWidth(style+1)

        ls_cntr = 0
        total_lumi = 0
        total_count = 0

        old_lumi = -1
        time_elapsed = 0.
        for ii in range(profile.GetN()):

            if ii==0: continue
            Linst = max(profile.GetPointY(ii), 1.e-9)
            time_duration = profile.GetPointX(ii) - profile.GetPointX(ii-1)
            time_elapsed += time_duration
            total_lumi += Linst*1.e-5 * time_duration
            
            # check which luminosity it is ... 
            which_lumi = -1
            which_npu = -1
            for index, sl in enumerate(switch_lumi):

                if sl[1] < Linst and Linst <= sl[0]:
                    which_lumi = sl[0]
                    which_npu = switch_npu[index]
                    break

            if which_lumi==-1: 
                print('   WARNING!!: no corresponding lumis!!!',"Linst",Linst)
                continue

            switch = True if which_lumi != old_lumi else False
            #@@print(ii,time_elapsed,Linst,sl[0],sl[1],which_lumi,which_npu,switch,old_lumi)
            old_lumi = which_lumi

            l1_rate = l1_his.Eval(which_lumi)
            l1_rate = hackRate(l1_rate,which_lumi) # <<< HACK HACK HACK
            spare = l1_max+l1_rate_corr - l1_rate
            if style==0: gspare.SetPoint(ls_cntr, profile.GetPointX(ii), spare)

#            print('L=', Linst, 
#                  ', which lumi=', which_lumi, 
#                  '(npu = ', which_npu, 
#                  '), l1 rate =', l1_his.Eval(which_lumi), 
#                  ', l1 b/w =', spare)

            # Determine L1 threshold ... 
            which_l1pt = l1_ptrange[-1] 
            flag_park = False
            l1_ee_rate = None
            for pt in l1_ptrange:
                histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                            str(pt).replace('.','p').replace('p0','')+
                                            'er1p22_dR_' + str(drdict[pt]).replace('.','p'))
                ee_rate = histo_ee_rate.Eval(which_npu)*1000
                ee_rate = hackRate(ee_rate,which_lumi) # <<< HACK HACK HACK
                l1_ok = ee_rate < allocation or ee_rate < (spare+allocation)
                if l1_ok:
                    which_l1pt = pt
                    flag_park = True
                    l1_ee_rate = ee_rate
                    break
            if not flag_park: 
                #print('\t L1 not available')
                continue

            #@@print('\t which_l1pt=', which_l1pt, 'with l1 ee rate = ', l1_ee_rate)

            # Determine HLT threshold ...
            hlt_file = TFile(path+'ee/roc_hlt_pu' + str(which_npu) + '.root')
            hlt_roc = hlt_file.Get('inv_pt' + str(which_l1pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

#            hlt_eff = -1
#            for ip in range(hlt_roc.GetN()):
#                if hlt_roc.GetPointX(ip) < max_bw_hlt[which_lumi]:
#                    hlt_eff = hlt_roc.GetPointY(ip)
#                    break
#            if hlt_eff==-1:
#                print('!!!!!! This cannot happen !!!!')
#            #@@print('\t parking: HLT eff. =', hlt_eff, 'max b/w=', max_bw_hlt[which_lumi])

            # Extract corrections to rates from JSON
            if options.limit:
                hlt_max = 300. #max_bw_hlt[peak_lumi] # <-- Limit to 300 Hz for prompt reco
                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                which_hltpt = -1
                for kk in range(hlt_n):
                    ip = kk#hlt_n - kk - 1
                    # print("test",kk,ip,hlt_pts[ip],hlt_roc.GetPointX(ip),hlt_max,hlt_roc.GetPointY(ip))
                    hlt_rate = hlt_roc.GetPointX(ip)
                    hlt_eff = hlt_roc.GetPointY(ip)
                    which_hltpt = np.arange(4, 11, 0.5).tolist()[ip]
                    if options.corrected:
                        l1ptstr = str(which_l1pt)
                        hltptstr = str(which_hltpt)
                        if which_npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[which_npu].keys() and \
                                hltptstr in corrs_dict[which_npu][l1ptstr].keys():
                            corr = corrs_dict[which_npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",which_npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",which_npu,l1ptstr,hltptstr)
                    if hlt_rate < hlt_max:
                        break
                if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')
            else:
                hlt_rate = -1
                hlt_eff = -1
                which_hltpt = hlt_threshold_dict.get(which_l1pt,4.0)
                hltpt_list = np.arange(4, 11, 0.5).tolist()
                index = hltpt_list.index(which_hltpt)
                if index < hlt_roc.GetN():
                    hlt_rate = hlt_roc.GetPointX(index)
                    hlt_eff = hlt_roc.GetPointY(index)
                    if options.corrected:
                        l1ptstr = str(which_l1pt)
                        hltptstr = str(which_hltpt)
                        if which_npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[which_npu].keys() and \
                                hltptstr in  corrs_dict[which_npu][l1ptstr].keys():
                            corr = corrs_dict[which_npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",which_npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",which_npu,l1ptstr,hltptstr)
                    else: print("cannot find hltpt")

            #print("test",which_npu,which_l1pt,which_hltpt,hlt_eff,hlt_rate)
            #quit()

            count = Linst * 1.e-5 * fB * Sigma_B * Br_kee * hlt_eff * time_duration
            total_count += count
            graph.SetPoint(ls_cntr, profile.GetPointX(ii), total_count)

            ls_cntr += 1
            if switch:
                if idx==0: switches.append((time_elapsed,which_l1pt,l1_ee_rate,spare))
                print(
                    " ".join(["ii:",str("{:5.0f}".format(ii)),
                              "dt:",str("{:4.1f}".format(time_duration)),
                              "time:",str("{:5.0f}".format(time_elapsed)),
                              "Peak:",str("{:4.2f}".format(which_lumi)),
                              "Linst:",str("{:.3f}".format(Linst)),
                              "Lint:",str("{:.3f}".format(total_lumi)),
                              "spare:",str("{:5.0f}".format(spare)),
                              "pT:",str("{:4.1f}".format(which_l1pt)),
                              "rate:",str("{:5.0f}".format(l1_ee_rate)),
                              "Eff:",str("{:.6f}".format(hlt_eff)),
                              #"dN:",str("{:.4f}".format(count)),
                              "N:",str("{:6.3f}".format(total_count))
                          ]))
 
        graphs.append(copy.deepcopy(graph))

        print("SUMMARY:") 
        print('Lumi profile name:       ', name) 
        print('L1 total rate:           ', l1_max)
        print('Di-ele allocation:       ', allocation)
        #print('Dimuon allocation:       ', dimuon)
        total_time = profile.GetPointX(profile.GetN()-1)
        print('Time duration            ', total_time, 's')
        print('Total lumi:              ', total_lumi,"fb^-1")
        print('Total # or Kee events:   ', total_count)
        factor = integrated_lumi/(total_lumi)#*pow(10.,-6))
        print('Target lumi:             ', integrated_lumi)
        print('Luminosity factor:       ', factor)
        print('Expected # of Kee events:', total_count*factor)
        summary.append((name,allocation,total_lumi,total_count))

    # First, lumi profile (only), then add estimates
    for only_profile in [True,False]:

        # Margins
        gStyle.SetPadTopMargin(0.08)
        gStyle.SetPadLeftMargin(0.14)
        if only_profile: gStyle.SetPadRightMargin(0.03)
        else: gStyle.SetPadRightMargin(0.26) # Extra space for multiple axes

        # Create canvas
        canvas = TCanvas('canvas_' +  name)
        canvas.SetGridx()
        canvas.SetGridy()
        gPad.SetTicks(1,0) # No tick marks on RHS

        # Create frame for lumi profile
        ymax = 2.5 # Maximum for y-axis (Linst)
        frame = TH2F('frame_'+name, 'frame_'+name, 
                     profile.GetN(), 0, total_time, # time axis
                     100, 0, ymax ) # counts axis
        frame.GetXaxis().SetTitle('Time [s]')
        frame.GetYaxis().SetTitle('L_{inst} [10^{34} Hz/cm^{2}]')
        frame.SetTitleOffset(1.2)
        frame.Draw()

        # Draw lumi profile
        profile.SetLineColor(kBlack)
        profile.SetLineWidth(3)
        profile.SetMarkerSize(0)
        profile.Draw('same')
        canvas.Update()

        # Create legend
        t = TLatex()
        t.SetTextAlign(11)
        t.SetTextSize(0.035)
        if not only_profile: t.DrawLatexNDC(0.47,0.87,"Allocation [kHz]")
        upper = 0.86
        lower = upper-len(graphs)*0.04
        leg = TLegend(0.45,lower,0.7,upper)
        applyLegendSettings(leg)
        leg.SetTextSize(0.035)
                
        # Draw graphs and right-hand axes
        if not only_profile:

            # Spare capacity
            right_max = l1_max * ymax / 2.0
            scale = gPad.GetUymax()/right_max
            scaleGraph(gspare,scale)
            gspare.Draw("same")
            xshift = (gPad.GetUxmax()-gPad.GetUxmin())*0.25
            axis2 = TGaxis(gPad.GetUxmax()+xshift,
                           gPad.GetUymin(),
                           gPad.GetUxmax()+xshift,
                           gPad.GetUymax(),
                           0,right_max,510,"+L")
            axis2.SetLineColor(kRed+1)
            axis2.SetLabelColor(kRed+1)
            axis2.SetTitleOffset(1.4)
            axis2.SetTitle('Spare L1 trigger capacity [Hz]')
            axis2.SetTitleColor(kRed+1)
            axis2.Draw()

            # Cumu candidates
            right_max = 12.5
            scale = gPad.GetUymax()/right_max
            for ii, graph in enumerate(graphs[::-1]): # in reverse order
                scaleGraph(graph,scale)
                graph.Draw('lsame')
                leg.AddEntry(graph, graph.GetTitle(), 'l')
            axis1 = TGaxis(gPad.GetUxmax(),
                           gPad.GetUymin(),
                           gPad.GetUxmax(),
                           gPad.GetUymax(),
                           0,right_max,510,"+L")
            axis1.SetLineColor(kBlue+2)
            axis1.SetLabelColor(kBlue+2)
            axis1.SetTitleOffset(1.2)
            axis1.SetTitle('Cumulative number of Kee candidates')
            axis1.SetTitleColor(kBlue+2)
            axis1.Draw()

        # Draw legend
        if not only_profile: leg.Draw()

#        lines = []
#        texts = []
#        last = 0.
#        for time,pt,rate,spare in switches:
#            print("line:",time,pt,rate,spare)
#            if time < 100 or time - last > 1000. :
#                lines.append(TLine(time,0.,time,10.))
#                l=lines[-1]
#                l.SetLineColor(kRed)
#                l.SetLineStyle(2)
#                l.SetLineWidth(1)
#                l.Draw()
#                #texts.append(TLatex(time+100.,1.,"Rate={:6.0f} Spare={:6.0f}".format(rate,spare)))
#                texts.append(TLatex(time+100.,1.,"Spare={:6.0f}".format(rate,spare)))
#                t = texts[-1]
#                t.SetTextAngle(90)
#                t.SetTextAlign(13)
#                t.SetTextColor(kRed)
#                t.SetTextSize(0.02)
#                t.Draw()
#            last = time

        # Add CMS labels
        l2=add_Private(text="#it{Private Work}",x=0.27)
        l2.Draw("same")
        l3=add_CMS(x=0.14)
        l3.Draw("same")
        l4=add_lumi(total_lumi,x=0.83 if only_profile else 0.6)
        l4.Draw("same")

        # Save canvas
        #canvas.RedrawAxis()
        canvas.Update()
        filename = 'plots/'+name+'.pdf'
        if not only_profile: filename = filename.replace("plots/","plots/estimates_for_")
        canvas.SaveAs(filename)
        del canvas

# Print summary
for name,allocation,lumi,count in summary:
    print('Name:',"{:30s}".format(name),
          'L1 alloc [Hz]:',"{:6.0f}".format(allocation),
          'Lint [/fb] per fill:',"{:6.3f}".format(lumi),
          'Counts, per fill:',"{:5.1f}".format(count),
          ', per /fb:',"{:5.1f}".format(count/lumi),
          ', per 25/fb:',"{:5.1f}".format(integrated_lumi*count/lumi)
    )

# Tables for estimates

dct = {}
for name,allocation,lumi,count in summary:
    if name not in dct.keys() : dct[name] = {}
    dct[name][allocation] = (lumi,count)

for name,counts in dct.items():
    print(name," ",end="")
    print(counts[0][0]," ",end="")
    for alloc in [0,5000,10000,20000]: print("{:5.1f}".format(counts[alloc][1])," ",end="")
    print()
print()

nfills = {}
nfills['levelled_at_2p0e34'] = 14
nfills['levelled_at_1p7e34'] = 14
nfills['levelled_at_1p5e34'] = 14
nfills['levelled_at_1p3e34'] = 14
nfills['levelled_at_1p1e34'] = 14
nfills['levelled_at_0p9e34'] = 14
nfills['levelled_at_0p7e34'] = 3
nfills['levelled_at_0p4e34'] = 3
nfills['levelled_at_0p2e34'] = 3

for name,counts in dct.items():
    print(name," ",end="")
    print("{:0.2f}".format(counts[0][0])," ",end="")
    print(nfills[name]," ",end="")
    for alloc in [0,5000,10000,20000]: print("{:5.1f}".format(counts[alloc][1]*nfills[name])," ",end="")
    print()
print()
