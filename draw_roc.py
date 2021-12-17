import copy, math, os, sys
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TLegend, TGraph
from officialStyle import officialStyle
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label
import numpy as np

l1_ptrange = np.arange(5, 11.5, 0.5).tolist() 
hlt_ptrange = np.arange(4, 11.5, 0.5).tolist() 

print('l1', l1_ptrange)
print('hlt', hlt_ptrange)

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-w', '--weight', action="store_true", default=True, dest='weight')
parser.add_option("-p", "--pu", default=1.0, type="float", help="target PU", dest="pu")
(options, args) = parser.parse_args()


effrefs = {
    'Mu12_IP6':0.00011,
    'Mu9_IP6':0.00025,
    'Mu9_IP5':0.00028,
    'Mu8_IP5':0.00040,
    'Mu7_IP4':0.00066
    }

def applyLegendSettings2(leg):     
    leg.SetBorderSize(0)           
    leg.SetFillColor(10)           
    leg.SetLineColor(0)            
    leg.SetFillStyle(0)            
    leg.SetTextSize(0.05)          
    leg.SetTextFont(42)     


def returnGraph(name, rates, effs):
    graph = TGraph()
    graph.SetName('pt' + str(name).replace('.0', ''))
    graph.SetTitle('L1 di-e p_{T} > ' + str(name).replace('.0','') + 'GeV')


    idx = 0
    for rate, eff in zip(rates, effs):

        graph.SetPoint(idx, eff, rate)

        idx += 1

    return copy.deepcopy(graph)


colours = [1, 2, 4, 6, 8, 13, 15]

# 12, 10, 9, 8, 7
colours_l1 = [1, 12,  9, 7, 5]
styles = [1, 2, 4, 3, 5, 1, 1,1]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def applyLegendSettings(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.05)
#    leg.SetTextFont(42)



def createPdf(rootfile, pname, xt, yt, prec):
    canvas = TCanvas(pname)

    canvas.cd()
    histo = rootfile.Get(pname)

    histo.Draw("TEXT COL")
    histo.GetXaxis().SetTitle(xt) #"sub-leading gen. electron p_{T} (GeV)")
    histo.GetYaxis().SetTitle(yt) #"leading gen. electron p_{T} (GeV)")
    gStyle.SetPaintTextFormat('.' + str(prec) + 'f');


    histo.SetStats(0)
    histo.SetMarkerSize(0.7)

    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")

    canvas.SaveAs("plots/hlt_" + pname + ".pdf")
    canvas.SaveAs("plots/hlt_" + pname + ".gif")

    canvas.SetLogz()
    canvas.SaveAs("plots/hlt_" + pname + "_log.pdf")
    canvas.SaveAs("plots/hlt_" + pname + "_log.gif")

    return canvas



def createROCPdf(effmap, file_rate, file_ref, name):

    graphs = []
    graphs_inv = []

    for l1pt in l1_ptrange:

        rates = []
        effs = []
    
        for hltpt in hlt_ptrange:


#            print(l1pt, hltpt, type(l1pt), type(hltpt))
            xbin = effmap.GetXaxis().FindBin(l1pt)
            ybin = effmap.GetYaxis().FindBin(hltpt)
            
#            rate = ratemap.GetBinContent(xbin, ybin)
#            print('getting ... l1_' + str(l1pt).replace('.','p') + '_hlt' + str(hltpt).replace('.','p') + '_fit')
            rate = file_rate.Get('l1_' + str(l1pt).replace('.','p') + '_hlt_' + str(hltpt).replace('.','p') + '_fit').Eval(options.pu)
            eff = effmap.GetBinContent(xbin, ybin)
            
            rates.append(rate)
            effs.append(eff)

        graph = returnGraph(l1pt, rates, effs)
        graph.SetName('pt' + str(l1pt).replace('.','p'))
        graph.SetTitle('pt' + str(l1pt).replace('.','p'))

        
        graph_inv = returnGraph(l1pt, effs, rates)
        graph_inv.SetName('inv_pt' + str(l1pt).replace('.','p'))
        graph_inv.SetTitle('inv_pt' + str(l1pt).replace('.','p'))

        graphs.append(copy.deepcopy(graph))
        graphs_inv.append(copy.deepcopy(graph_inv))


    graphs_ref = []


    for iref, refname in enumerate(['Mu12_IP6', 'Mu9_IP6', 'Mu9_IP5', 'Mu8_IP5', 'Mu7_IP4']):

        graph_ref = TGraph()
        graph_ref.SetName('ref_' + refname)
        graph_ref.SetTitle('ref_' + refname)

        rate_ref = file_ref.Get(refname)

#        print(refname, 0, effrefs[refname], rate_ref.Eval(options.pu)*2544.)
        graph_ref.SetPoint(0, effrefs[refname], rate_ref.Eval(options.pu)*2544.)
        if iref==0:
            graph_ref.SetMarkerStyle(30)
        else:
            graph_ref.SetMarkerStyle(29)

        graph_ref.SetMarkerColor(colours_l1[iref])
        graph_ref.SetMarkerSize(3)
        graphs_ref.append(graph_ref)

    canvas = TCanvas(name)
    canvas.SetLogy()
    canvas.SetLogx()
    canvas.SetGridx()
    canvas.SetGridy()

    frame_roc = TH2F(name, name, 100, 0.000003, 0.01, 1000, 1., 100000)

    frame_roc.GetXaxis().SetTitle('L1 x HLT Trigger efficiency')
    frame_roc.GetYaxis().SetTitle('HLT Trigger Rate (Hz)')
    frame_roc.Draw()

    leg = TLegend(0.2, 0.25,0.5,0.81)

    applyLegendSettings(leg)
    leg.SetTextSize(0.03)


    for idx, graph in enumerate(graphs):
#        print(idx)

        col = idx + 1
        if col>=10:
            col += 1

        graph.SetLineColor(col)
        graph.SetMarkerColor(col)
        graph.SetMarkerSize(1)
        graph.Draw('plsame')
        leg.AddEntry(graph, 'Level-1 p_{T} > ' + graph.GetTitle().replace('pt','').replace('p','.'), 'l')



    for graph_ in graphs_ref:
        graph_.Draw('psame')
        leg.AddEntry(graph_, graph_.GetName().replace('ref_',''), 'p')

    leg.Draw()
    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")


    leg2 = TLegend(0.19, 0.81,0.3,0.85)
    applyLegendSettings2(leg2)
    leg2.SetTextSize(0.032)
    leg2.SetTextFont(62)
    leg2.AddEntry(frame_roc, 'L = {0:.1f}'.format(float(options.pu*0.0357338 - 0.0011904)) + 'E34 (PU ~ ' + str(options.pu) + ')', '')
    leg2.Draw()



    canvas.SaveAs('plots/' + name + '_' + str(options.pu).replace('.','p') + '.gif')
    canvas.SaveAs('plots/' + name + '_' + str(options.pu).replace('.','p') + '.pdf')

    file = TFile('roc_' + name + '.root', 'recreate')
    for idx, graph in enumerate(graphs):
        graph.Write()
    for idx, graph in enumerate(graphs_inv):
        graph.Write()
    file.Write()
    file.Close()



ensureDir('plots')

file_rate = TFile('root/ratemap.root')
#ratemap = file_rate.Get('rate')
#ratemap_mass = file_rate.Get('rate_mass')

createPdf(file_rate, 'rate', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)
createPdf(file_rate, 'rate_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)

file_eff = TFile('root/effmap.root')
effmap = file_eff.Get('gall')

file_ref = TFile('single-mu/dict.root')
#ref = file_ref.Get('Mu9_IP6')

createPdf(file_eff, 'e1', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
createPdf(file_eff, 'e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
createPdf(file_eff, 'e1e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
createPdf(file_eff, 'all', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
createPdf(file_eff, 'gall', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )
createPdf(file_eff, 'gall_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )

#sys.exit(1)

createROCPdf(effmap, file_rate, file_ref, 'roc_hlt')
#createROCPdf(effmap, file_rate, file_ref, 'roc_mass_hlt')


