#import os, math, sys
#from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num

from TreeProducerBcJpsiTauNu import *
from DeltaR import deltaR, deltaPhi
import copy
import random
import numpy as np
import itertools

def hlt_criteria(chain, idx, pt):
    
    flag = False


    
    if chain.eg_et[idx] > float(pt) \
       and chain.eg_pms2[idx] < 10000 \
       and chain.eg_invEInvP[idx] < 0.2 \
       and chain.eg_trkDEtaSeed[idx] < 0.01 \
       and chain.eg_trkDPhi[idx] < 0.2 \
       and chain.eg_trkChi2[idx] < 40 \
       and chain.eg_trkValidHits[idx] >= 5 \
       and chain.eg_trkNrLayerIT[idx] >= 2:

        flag = True


    return flag




def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)

parser.add_option("-o", "--out", default='out_mc_gabagaba.root', type="string", help="output filename", dest="out")
#parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/mc_widedPhiMaxLowEtGrad/Myroot.root', type="string", help="file", dest="file")
parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/mc_gabagaba/Myroot.root', type="string", help="file", dest="file")

(options, args) = parser.parse_args()

print(options)


out = TreeProducerBcJpsiTauNu(options.out)

chain = ROOT.TChain('egHLTRun3Tree', 'tree')
chain.AddFile(options.file)

Nevt = chain.GetEntries()

print('Total Number of events = ', Nevt)
evtid = 0
#evtid_7 = 0

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


#    print('-'*80)

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
####        dr = deltaR(chain.eg_l1eg_eta[ele1], chain.eg_l1eg_phi[ele1],
####                    chain.eg_l1eg_eta[ele2], chain.eg_l1eg_phi[ele2])
###
####        print(ele1, ele2, 'dr=', dr)
###
###
###        if deltaR(chain.eg_l1eg_eta[ele1], chain.eg_l1eg_phi[ele1],
###                  chain.eg_l1eg_eta[ele2], chain.eg_l1eg_phi[ele2]) < 0.9: 
###
###            evtFlag = True
###
###    if not evtFlag:
###        continue



    for ii in range(len(chain.eg_gen_pt)):
        if abs(chain.eg_gen_eta[ii]) > 1.2: continue


#        print('gen', ii, chain.eg_gen_pt[ii], chain.eg_gen_eta[ii], chain.eg_gen_phi[ii])

        out.gen_pt[0] = chain.eg_gen_pt[ii]
        out.gen_eta[0] = chain.eg_gen_eta[ii]
        out.gen_phi[0] = chain.eg_gen_phi[ii]
        out.gen_energy[0] = chain.eg_gen_energy[ii]


        l1index = -1
        maxl1_dr = 99
        maxl1_dphi = 99
        maxl1_deta = 99
        
        for jj in range(len(chain.eg_l1eg_et)):


            if abs(chain.eg_l1eg_eta[jj]) > 1.2: continue

            dr = deltaR(chain.eg_gen_eta[ii], chain.eg_gen_phi[ii],
                        chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])

#            print('\t l1eg:', chain.eg_l1eg_et[jj], chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj], '----->', dr, maxl1_dr)
            if maxl1_dr > dr:

#                print('This is good!')
                
                maxl1_dphi = deltaPhi(chain.eg_l1eg_phi[jj], chain.eg_gen_phi[ii])
                maxl1_deta = chain.eg_l1eg_eta[jj] - chain.eg_gen_eta[ii]
                maxl1_dr = dr
#                l1o = (chain.eg_l1eg_et[jj], chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])


                l1index = jj
                


        if l1index!=-1:
            out.l1_pt[0] = chain.eg_l1eg_et[l1index]
            out.l1_eta[0] = chain.eg_l1eg_eta[l1index]
            out.l1_phi[0] = chain.eg_l1eg_phi[l1index]
            out.l1_dr[0] = maxl1_dr
            out.l1_dphi[0] = maxl1_dphi
            out.l1_deta[0] = maxl1_deta

        else:
            out.l1_pt[0] = -1
            out.l1_eta[0] = -1
            out.l1_phi[0] = -1
            out.l1_dr[0] = -1
            out.l1_dphi[0] = -1
            out.l1_deta[0] = -1


        hltindex = -1
        maxhlt_dr = 99
        maxhlt_dphi = 99
        maxhlt_deta = 99

        for jj in range(len(chain.eg_et)):

            dr = deltaR(chain.eg_gen_eta[ii], chain.eg_gen_phi[ii],
                        chain.eg_eta[jj], chain.eg_phi[jj]) 

            if abs(chain.eg_eta[jj]) > 1.2: continue
            
            if maxhlt_dr > dr:
                
                maxhlt_dr = dr
                maxhlt_dphi = deltaPhi(chain.eg_phi[jj], chain.eg_gen_phi[ii])
                maxhlt_deta = chain.eg_eta[jj] - chain.eg_gen_eta[ii]
#                hlto = (chain.eg_et[jj], chain.eg_eta[jj], chain.eg_phi[jj])
                hltindex = jj



        if hltindex!=-1:
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
            out.hlt_dr[0] = maxhlt_dr
            out.hlt_dphi[0] = maxhlt_dphi
            out.hlt_deta[0] = maxhlt_deta
            out.hlt_dxy[0] = chain.eg_gsfdxy[hltindex]

        else:
            out.hlt_pt[0] = -1
            out.hlt_eta[0] = -1
            out.hlt_phi[0] = -1
            out.hlt_energy[0] = -1
            out.hlt_rawEnergy[0] = -1
            out.hlt_phiWidth[0] = -1
            out.hlt_nrClus[0] = -1
            out.hlt_seedId[0] = -1
            out.hlt_seedDet[0] = -1
            out.hlt_seednrCrystals[0] = -1
            out.hlt_sigmaIEtaIEta[0] = -1
            out.hlt_sigmaIEtaIEtaNoise[0] = -1
            out.hlt_ecalPFIsol[0] = -1
            out.hlt_hcalPFIsol[0] = -1
            out.hlt_trkIsol[0] = -1
            out.hlt_trkChi2[0] = -1
            out.hlt_trkMissHits[0] = -1
            out.hlt_trkNrLayerIT[0] = -1
            out.hlt_trkValidHits[0] = -1
            out.hlt_hcalHForHoverE[0] = -1
            out.hlt_invESeedInvP[0] = -1
            out.hlt_invEInvP[0] = -1
            out.hlt_trkDEta[0] = -1
            out.hlt_trkDEtaSeed[0] = -1
            out.hlt_trkDPhi[0] = -1
            out.hlt_pms2[0] = -1
            out.hlt_dr[0] = -1
            out.hlt_dphi[0] = -1
            out.hlt_deta[0] = -1
            out.hlt_dxy[0] = -1
        
        out.tree.Fill()



##    eles = []
##    for hltindex in range(len(chain.eg_et)):
##
##        if abs(chain.eg_eta[hltindex]) > 1.2: continue
##
##        eles.append(hltindex)
##
##
###    for pt in [7]:
##
##    flag_7 = False
##
##    for ele1,ele2 in itertools.combinations(eles,2):
##
##        if not hlt_criteria(chain, ele1, 7): continue
##        if not hlt_criteria(chain, ele2, 7): continue
##        
##        if deltaR(chain.eg_eta[ele1], chain.eg_phi[ele1],
##                  chain.eg_eta[ele2], chain.eg_phi[ele2]) < 0.9: 
##
##            flag_7 = True
        

#        flag_7 = True
#        break

#        getattr(out, 'doubleE' + str(pt) + '_nomatch')[0] = flag




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
#
#    if flag_7:
#        evtid_7 += 1
    

#print(Nevt, 'evt processed.', evtid, evtid_7, 'evt has matching')
print(Nevt, 'evt processed.', evtid, 'evt has matching')

out.endJob()
