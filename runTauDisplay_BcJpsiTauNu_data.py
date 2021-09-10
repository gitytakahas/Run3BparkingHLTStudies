#import os, math, sys
#from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num

from TreeProducerBcJpsiTauNu import *
from DeltaR import deltaR, deltaPhi
import copy
import random
import numpy as np
import itertools

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def hlt_criteria(chain, idx):
    
    flag = False
    
    if chain.eg_et[idx] > 7 \
       and chain.eg_pms2[idx] < 10000 \
       and chain.eg_invEInvP[idx] < 0.2 \
       and chain.eg_trkDEtaSeed[idx] < 0.01 \
       and chain.eg_trkDPhi[idx] < 0.2 \
       and chain.eg_trkChi2[idx] < 50 \
       and chain.eg_trkValidHits[idx] >= 10 \
       and chain.eg_trkNrLayerIT[idx] >= 2:

        flag = True

    return flag 
            
    


from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)

parser.add_option("-o", "--out", default='out_data_4gev_22.root', type="string", help="output filename", dest="out")

parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/data_4gev_22/Myroot.root', type="string", help="file", dest="file")

(options, args) = parser.parse_args()

print(options)


out = TreeProducerBcJpsiTauNu(options.out)

chain = ROOT.TChain('egHLTRun3Tree', 'tree')
chain.AddFile(options.file)

Nevt = chain.GetEntries()

chain.SetBranchStatus('L1Upgrade*', 0)
chain.SetBranchStatus('path*', 0)
chain.SetBranchStatus('eg_gen*', 0)


print('Total Number of events = ', Nevt)
evtid = 0

drdict = {
    3:1,
    4:1,
    5:0.9,
    6:0.8,
    7:0.75,
    8:0.7,
    9:0.65,
    10:0.6,
    11:0.55,
    12:0.5,
    13:0.45,
    14:0.4,
}   


for evt in range(Nevt):
    chain.GetEntry(evt)

    if evt%10000==0: print('{0:.2f}'.format(float(evt)/float(Nevt)*100.), '% processed')

#    if evt==10000:
#        break

    l1eles = []

    for l1index in range(len(chain.eg_l1eg_et)):

        if abs(chain.eg_l1eg_eta[l1index]) > 1.2: continue
#        if chain.eg_l1eg_et[l1index] < 5.: continue

        l1eles.append(l1index)


    for pt in [5,6,7,8]:

        flag = False

        for ele1,ele2 in itertools.combinations(l1eles,2):

            if chain.eg_l1eg_et[ele1] < pt: continue
            if chain.eg_l1eg_et[ele2] < pt: continue

            if deltaR(chain.eg_l1eg_eta[ele1], chain.eg_l1eg_phi[ele1], chain.eg_l1eg_eta[ele2], chain.eg_l1eg_phi[ele2]) > float(drdict[pt]): continue

            flag = True
            break

        getattr(out, 'l1_doubleE' + str(pt))[0] = flag



###    l1eles = []
###    for l1index in range(len(chain.eg_l1eg_et)):
###
###        if abs(chain.eg_l1eg_eta[l1index]) > 1.2: continue
###        if chain.eg_l1eg_et[l1index] < 5.: continue
###
###        l1eles.append(l1index)
###
###
###    evtFlag = False
###    for ele1,ele2 in itertools.combinations(l1eles,2):
###
###        if deltaR(chain.eg_l1eg_eta[ele1], chain.eg_l1eg_phi[ele1],
###                  chain.eg_l1eg_eta[ele2], chain.eg_l1eg_phi[ele2]) < 0.9: 
###            
###            evtFlag = True        
###
###    if not evtFlag:
###        continue

#
#        if dr < 0.4 and max_dr > dr:
#                
#            max_dr = dr
#            l1index = jj



#    for ii in range(len(chain.eg_gen_pt)):
#        if abs(chain.eg_gen_eta[ii]) > 1.2: continue
#
#
#        out.gen_pt[0] = chain.eg_gen_pt[ii]
#        out.gen_eta[0] = chain.eg_gen_eta[ii]
#
#
#        l1index = -1
#        max_dr = 99
#        
#        for jj in range(len(chain.eg_l1eg_et)):
#
#            dr = deltaR(chain.eg_gen_eta[ii], chain.eg_gen_phi[ii],
#                      chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])
#
#            if dr < 0.4 and max_dr > dr:
#                
#                max_dr = dr
##                l1o = (chain.eg_l1eg_et[jj], chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])
#                l1index = jj
#                
#
#
#        if l1index!=-1:
#            out.l1_pt[0] = chain.eg_l1eg_pt[l1index]
#            out.l1_eta[0] = chain.eg_l1eg_eta[l1index]
#            out.l1_phi[0] = chain.eg_l1eg_phi[l1index]
#            out.l1_dr[0] = max_dr
#        else:
#            out.l1_pt[0] = -1
#            out.l1_eta[0] = -1
#            out.l1_phi[0] = -1
#            out.l1_dr[0] = -1
#
#
#        hltindex = -1
#        max_dr = 99
#
#        for jj in range(len(chain.eg_et)):
#
#            dr = deltaR(chain.eg_gen_eta[ii], chain.eg_gen_phi[ii],
#                        chain.eg_eta[jj], chain.eg_phi[jj]) 
#
#            if dr < 0.4 and max_dr > dr:
#                
#                max_dr = dr
##                hlto = (chain.eg_et[jj], chain.eg_eta[jj], chain.eg_phi[jj])
#                hltindex = jj
#
#
#

#    eles = []
    for hltindex in range(len(chain.eg_et)):

        if abs(chain.eg_eta[hltindex]) > 1.2: continue

        out.hlt_pt[0] = chain.eg_et[hltindex]
        out.hlt_eta[0] = chain.eg_eta[hltindex]
        out.hlt_phi[0] = chain.eg_phi[hltindex]
        out.hlt_energy[0] = chain.eg_energy[hltindex]
        out.hlt_rawEnergy[0] = chain.eg_rawEnergy[hltindex]
        out.hlt_phiWidth[0] = chain.eg_phiWidth[hltindex]
        out.hlt_nrClus[0] = chain.eg_nrClus[hltindex]
        out.hlt_seedId[0] = chain.eg_seedId[hltindex]
        out.hlt_seedDet[0] = chain.eg_seedDet[hltindex]
        out.hlt_seednrCrystals[0] = chain.eg_seednrCrystals[hltindex]
        out.hlt_sigmaIEtaIEta[0] = chain.eg_sigmaIEtaIEta[hltindex]
        out.hlt_sigmaIEtaIEtaNoise[0] = chain.eg_sigmaIEtaIEtaNoise[hltindex]
        out.hlt_ecalPFIsol[0] = chain.eg_ecalPFIsol[hltindex]
        out.hlt_hcalPFIsol[0] = chain.eg_hcalPFIsol[hltindex]
        out.hlt_trkIsol[0] = chain.eg_trkIsol[hltindex]
        out.hlt_trkChi2[0] = chain.eg_trkChi2[hltindex]
        out.hlt_trkMissHits[0] = chain.eg_trkMissHits[hltindex]
        out.hlt_trkNrLayerIT[0] = chain.eg_trkNrLayerIT[hltindex]
        out.hlt_trkValidHits[0] = chain.eg_trkValidHits[hltindex]
        out.hlt_hcalHForHoverE[0] = chain.eg_hcalHForHoverE[hltindex]
        out.hlt_invESeedInvP[0] = chain.eg_invESeedInvP[hltindex]
        out.hlt_invEInvP[0] = chain.eg_invEInvP[hltindex]
        out.hlt_trkDEta[0] = chain.eg_trkDEta[hltindex]
        out.hlt_trkDEtaSeed[0] = chain.eg_trkDEtaSeed[hltindex]
        out.hlt_trkDPhi[0] = chain.eg_trkDPhi[hltindex]
        out.hlt_pms2[0] = chain.eg_pms2[hltindex]
        out.hlt_dr[0] = -1
        out.hlt_dphi[0] = -1
        out.hlt_deta[0] = -1
        out.hlt_dxy[0] = chain.eg_gsfdxy[hltindex]        

        out.isgjson[0] = chain.isgjson
        out.instL[0] = chain.instL
        out.npu[0] = chain.npu
        
        out.tree.Fill()


        
#        eles.append(hltindex)

    
#    flag = False
#
#    for ele1,ele2 in itertools.combinations(eles,2):
#        if not hlt_criteria(chain, ele1): continue
#        if not hlt_criteria(chain, ele2): continue
#        
#        flag = True


    


##    for pt in ptrange:
##
##        flag = False
##
##        for ii, e1 in enumerate(hlt_electrons):
##            
##            if flag: break
##
##            for jj, e2 in enumerate(hlt_electrons):
##
##                if jj <= ii: continue
##
##
##                if e1[0] >= pt and e2[0] >= pt and deltaR(e1[1], e1[2], e2[1], e2[2]) < float(drdict[pt]):
##                        
##                    flag = True
##                    break
##                        
##
##        getattr(out, 'doubleE' + str(pt))[0] = flag





#    for ii in range(len(chain.eg_et)):
#        if chain.eg_et[ii] > 5.


#    if options.priority=='pt':
#            
#        for itau in range(len(chain.JpsiTau_tau_pt)):
#            if chain.JpsiTau_tau_vprob[itau] < 0.1: continue
#            if chain.JpsiTau_tau_fls3d[itau] < 3.: continue
#            if chain.JpsiTau_tau_mass[itau] > 1.7: continue
#            if bool(chain.JpsiTau_tau_pi1_trigMatch[itau])==False and bool(chain.JpsiTau_tau_pi2_trigMatch[itau])==False and bool(chain.JpsiTau_tau_pi3_trigMatch[itau])==False: 
##                print 'trigger matching was not satisifed ...'
#                continue
#            
#            # you can add tau mass cuts here
#
#            tindex_ = itau
#            break


#        tindex_ = 0
            

#    out.puweight[0] = ROOT.Double(putool.getWeight(chain.nPuVtxTrue[0]))


    evtid += 1
    

print(Nevt, 'evt processed.', evtid, 'evt has matching')

out.endJob()
