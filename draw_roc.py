import copy, math, os, sys
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TLegend, TGraph
from officialStyle import officialStyle
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings, applyLegendSettings2
import numpy as np

l1_ptrange = np.arange(5, 12, 1.0).tolist() 
hlt_ptrange = np.arange(4, 12, 1.0).tolist() 

print('l1', l1_ptrange)
print('hlt', hlt_ptrange)

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-w', '--weight', action="store_true", default=False, dest='weight')
(options, args) = parser.parse_args()


effrefs = {
    #'Mu12_IP6':0.00011,
    #'Mu9_IP6':0.00025,
    'Mu9_IP6':0.00000329
    #'Mu9_IP5':0.00028,
    #'Mu8_IP5':0.00040,
    #'Mu7_IP4':0.00066
    }

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



def createROCPdf(effmap, file_rate, file_ref, npu, name):

    graphs = []
    graph_envelope = TGraph()
    graph_envelope.SetName('envelope')
    graph_envelope.SetTitle('envelope')
    graph_envelope_inv = TGraph()
    graph_envelope_inv.SetName('envelope_inv')
    graph_envelope_inv.SetTitle('envelope_inv')

    ii = 0
    
    for l1pt in l1_ptrange:

        rates = []
        effs = []
    
        for hltpt in hlt_ptrange:


#            print(l1pt, hltpt, type(l1pt), type(hltpt))
            xbin = effmap.GetXaxis().FindBin(l1pt)
            ybin = effmap.GetYaxis().FindBin(hltpt)
            
#            rate = ratemap.GetBinContent(xbin, ybin)
#            print('getting ... l1_' + str(l1pt).replace('.','p') + '_hlt' + str(hltpt).replace('.','p') + '_fit')
            rate = file_rate.Get('l1_' + str(l1pt).replace('.','p') + '_hlt_' + str(hltpt).replace('.','p') + '_fit').Eval(npu)
            eff = effmap.GetBinContent(xbin, ybin)
            
            rates.append(rate)
            effs.append(eff)

            ##### 
            if hltpt == l1pt - 1.:

#                print('(l1, hlt) = ', l1pt, hltpt)

                graph_rep = TGraph()
                graph_rep.SetName('rep_l1pt' + str(l1pt).replace('.','p') + '_hltpt' + str(hltpt).replace('.','p'))
                graph_rep.SetTitle('rep_l1pt' + str(l1pt).replace('.','p') + '_hltpt' + str(hltpt).replace('.','p'))
                
                graph_rep.SetPoint(0, eff, rate)

                graph_rep.SetMarkerColor(1)
                graph_rep.SetMarkerSize(4)
                graph_rep.SetMarkerStyle(42)

                graphs.append(copy.deepcopy(graph_rep))

                graph_envelope.SetPoint(ii, eff, rate)
                graph_envelope_inv.SetPoint(ii, rate, eff)
                ii += 1


        graph = returnGraph(l1pt, rates, effs)
        graph.SetMarkerSize(1)
        graph.SetName('pt' + str(l1pt).replace('.','p'))
        graph.SetTitle('pt' + str(l1pt).replace('.','p'))

        
        graph_inv = returnGraph(l1pt, effs, rates)
        graph_inv.SetName('inv_pt' + str(l1pt).replace('.','p'))
        graph_inv.SetTitle('inv_pt' + str(l1pt).replace('.','p'))

        graphs.append(copy.deepcopy(graph))
        graphs.append(copy.deepcopy(graph_inv))

    graphs.append(copy.deepcopy(graph_envelope))
    graphs.append(copy.deepcopy(graph_envelope_inv))

    graphs_ref = []


    for iref, refname in enumerate(['Mu9_IP6']):#['Mu12_IP6', 'Mu9_IP6', 'Mu9_IP5', 'Mu8_IP5', 'Mu7_IP4']):

        graph_ref = TGraph()
        graph_ref.SetName('ref_' + refname)
        graph_ref.SetTitle('ref_' + refname)

        rate_ref = file_ref.Get(refname)

#        print(refname, 0, effrefs[refname], rate_ref.Eval(options.pu)*2544.)
        graph_ref.SetPoint(0, effrefs[refname], rate_ref.Eval(npu)*2544.)
        if iref==0:
            graph_ref.SetMarkerStyle(30)
        else:
            graph_ref.SetMarkerStyle(29)

        graph_ref.SetMarkerColor(colours_l1[iref])
        graph_ref.SetMarkerSize(3)
        graphs_ref.append(graph_ref)


    # no envelope
    makeCanvas(name, options.weight, False, graphs, graphs_ref, npu)

    # with envelope
    makeCanvas(name, options.weight, True, graphs, graphs_ref, npu)



def makeCanvas(name, weight, envelope, graphs, graphs_ref, npu):

    canvas = TCanvas(name)
    canvas.SetLogy()
    canvas.SetLogx()
    canvas.SetGridx()
    canvas.SetGridy()

    if not weight: 
        frame_roc = TH2F(name, name, 100, 0.000003, 0.01, 1000, 0.1, 100000)
        frame_roc.GetXaxis().SetTitle('L1 x HLT Trigger efficiency')
    else : 
        frame_roc = TH2F(name, name, 100, 0.0000003, 0.0003, 1000, 0.1, 100000)
        frame_roc.GetXaxis().SetTitle('L1 x HLT x analysis eff.')


    frame_roc.GetYaxis().SetTitle('HLT Trigger Rate (Hz)')
    frame_roc.Draw()

    
    leg = TLegend(0.2, 0.25,0.5,0.65)
        
    applyLegendSettings(leg)
    leg.SetTextSize(0.03)

    
    if not envelope:

        col_ = 1

        for idx, graph in enumerate(graphs):
            if graph.GetName().find('rep')==-1 and graph.GetName().find('envelope')==-1 and graph.GetName().find('inv')==-1:

                col = col_
                if col >=10:
                    col += 1

                    
                graph.SetLineColor(col)
                graph.SetMarkerColor(col)
                leg.AddEntry(graph, 'Level-1 p_{T} > ' + graph.GetTitle().replace('pt','').replace('p','.'), 'l')

                col_ += 1

            graph.Draw('plsame')

#        leg.Draw()

    else:

        for idx, graph in enumerate(graphs):
            if graph.GetName()!='envelope': continue
            
#            print(idx, graph.GetName())

            graph.SetLineStyle(2)
            graph.SetLineWidth(3)
            graph.SetMarkerSize(2)
            graph.SetMarkerStyle(21)
            graph.Draw('plsame')
            leg.AddEntry(graph, 'L1 p_{T} = [5, 11, 0.5], HLT = L1 - 1 GeV', 'l')

    leg.Draw()

    leg2 = TLegend(0.19, 0.81,0.3,0.85)
    applyLegendSettings2(leg2)
    leg2.SetTextSize(0.032)
    leg2.SetTextFont(62)
    leg2.AddEntry(frame_roc, 'L = {0:.1f}'.format(float(npu*0.0357338 - 0.0011904)) + 'E34 (PU ~ ' + str(npu) + ')', '')
    leg2.Draw()


    for graph_ in graphs_ref:
        graph_.Draw('psame')
        leg.AddEntry(graph_, graph_.GetName().replace('ref_',''), 'p')


    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")


    suffix1 = '_envelope' if envelope else ''
    suffix2 = '_weighted' if weight else ''

    canvas.SaveAs('plots/' + name + '_' + str(npu).replace('.','p') +suffix1+suffix2+ '.gif')
    canvas.SaveAs('plots/' + name + '_' + str(npu).replace('.','p') +suffix1+suffix2+ '.pdf')
        

    if not envelope:
        file = TFile('root/' + name + '_pu' + str(int(npu)) + '.root', 'recreate')

        for idx, graph in enumerate(graphs):
            graph.Write()

        file.Write()
        file.Close()

#    for idx, graph in enumerate(graphs_inv):
#        graph.Write()

        
#    graph_envelope.Write()
#    graph_envelope_inv.Write()





ensureDir('plots')
ensureDir('root')

file_rate = TFile('/eos/cms/store/group/phys_bphys/bpark/RootFiles4Run3Parking/ee/ratemap4roc.root')

#createPdf(file_rate, 'rate', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)
#createPdf(file_rate, 'rate_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)

file_eff = None

if not options.weight: file_eff = TFile('/eos/cms/store/group/phys_bphys/bpark/RootFiles4Run3Parking/ee/effmap4roc.root')
else: file_eff = TFile('/eos/cms/store/group/phys_bphys/bpark/RootFiles4Run3Parking/ee/effmap4roc_weighted.root')

effmap = file_eff.Get('gall')

file_ref = TFile('/eos/cms/store/group/phys_bphys/bpark/RootFiles4Run3Parking/single-mu/obs_rate_summary_fit.root')
#ref = file_ref.Get('Mu9_IP6')

#createPdf(file_eff, 'e1', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'e1e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'all', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'gall', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )
#createPdf(file_eff, 'gall_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )

#sys.exit(1)

for npu in [56, 48, 42, 36, 30, 25, 17]:

    createROCPdf(effmap, file_rate, file_ref, npu, 'roc_hlt')

#createROCPdf(effmap, file_rate, file_ref, 'roc_mass_hlt')


