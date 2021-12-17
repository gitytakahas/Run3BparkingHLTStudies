import copy, math, os
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TColor, kLightTemperature, TGraphErrors, kRed, TLegend, kGray
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
from array import array
import numpy as np
import itertools

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
#gStyle.SetTitleOffset(1.0,"X")
#gStyle.SetTitleOffset(1.0,"Y")

from optparse import OptionParser, OptionValueError
usage = "usage: python estimate.py"
parser = OptionParser(usage)

parser.add_option('-p', '--pseudo', action="store_true", default=False, dest='pseudo')

(options, args) = parser.parse_args()

print(options)



def applyHistStyle(h, i):
#    print h, i
    h.SetLineColor(colours[i])
    h.SetMarkerColor(colours[i])
    h.SetMarkerSize(0)
    h.SetLineStyle(styles[i])
    h.SetLineWidth(3)
    h.SetStats(False)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isEff = False, isRatio=False, isLegend=False, x=None, y=None, label="test"):

    display = DisplayManager(pname, isLog, isRatio, 0.6, 0.7, x, y, label)
    display.draw_legend = isLegend
    display.isEff = isEff
    
    display.Draw(hists, titles)



ensureDir('plots/')


lumi = None
prompt = None
hists = []

for bw in ['80', '90']:

    file2read = 'n_' + bw + '000.0.root'
    
    if options.pseudo:
        file2read = 'n_' + bw + '000.0_pseudo.root'

    file = TFile(file2read)
    
    if lumi==None:
        lumi = copy.deepcopy(file.Get('lumiprofile'))
        prompt = copy.deepcopy(file.Get('n_prompt'))

    hist = file.Get('n')
    hist.SetName('bw_' + bw + ' kHz')
    hist.SetTitle('bw_' + bw + ' kHz')
    
    hists.append(copy.deepcopy(hist))


canvas = TCanvas()
canvas.SetGridx()
canvas.SetGridy()

leg = TLegend(0.2, 0.68,0.5,0.86)

applyLegendSettings(leg)
leg.SetTextSize(0.03)

ymax = max([hist.GetMaximum() for hist in hists])

frame = TH2F('frame', 'frame', lumi.GetN(), lumi.GetPointX(0), lumi.GetPointX(lumi.GetN()-1), 100,0,4)
frame.GetXaxis().SetTitle('Time (s)')
#frame.GetYaxis().SetTitle('# of Kee events')
frame.GetYaxis().SetTitle('Instantaneous Luminosity (E34)')
frame.Draw()


for ii, hist in enumerate(hists):

    hist.SetMaximum(4)
    
#    if ii==0:
#        hist.Draw('apl')
#    else:
#    hist.Draw('plsame')

    hist.SetLineColor(ii+1)
    hist.SetLineWidth(3)
    hist.SetMarkerColor(ii+1)
    hist.SetMarkerSize(1)


    leg.AddEntry(hist, hist.GetName().replace('bw_', 'L1 max = '), 'lep')

#lumi.SetLineColor(kGray)
lumi.SetLineWidth(3)
lumi.Draw('same')

prompt.SetLineColor(4)
prompt.SetLineWidth(2)
#prompt.Draw('same')

#leg.AddEntry(lumi, 'Fill7321 (0.55/fb)', 'lep')
#leg.AddEntry(prompt, 'Prompt (5kHz at L1, 100Hz at HLT)', 'lep')

#leg.Draw()

canvas.RedrawAxis();

filename = 'plots/estimate.gif'

if options.pseudo:
    filename = 'plots/estimate_pseudo.gif'

l2=add_Preliminary()
l2.Draw("same")
l3=add_CMS()
l3.Draw("same")


canvas.SaveAs(filename)
canvas.SaveAs(filename.replace('.gif', '.pdf'))


