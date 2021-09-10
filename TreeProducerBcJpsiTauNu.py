import ROOT
import math 
from TreeProducerCommon import *


class TreeProducerBcJpsiTauNu(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, **kwargs):
        print('TreeProducerBsTauTau is called for', name)
        super(TreeProducerBcJpsiTauNu, self).__init__(name,**kwargs)

        self.addBranch('gen_pt', 'f')
        self.addBranch('gen_eta', 'f')
        self.addBranch('gen_phi', 'f')
        self.addBranch('gen_energy', 'f')

        self.addBranch('l1_pt', 'f')
        self.addBranch('l1_eta', 'f')
        self.addBranch('l1_phi', 'f')
        self.addBranch('l1_dr', 'f')
        self.addBranch('l1_dphi', 'f')
        self.addBranch('l1_deta', 'f')

        self.addBranch('hlt_pt', 'f')
        self.addBranch('hlt_eta', 'f')
        self.addBranch('hlt_phi', 'f')
        self.addBranch('hlt_energy', 'f')
        self.addBranch('hlt_rawEnergy', 'f')
        self.addBranch('hlt_phiWidth', 'f')
        self.addBranch('hlt_nrClus', 'i')
        self.addBranch('hlt_seedId', 'f')
        self.addBranch('hlt_seedDet', 'i')
        self.addBranch('hlt_seednrCrystals', 'i')
        self.addBranch('hlt_sigmaIEtaIEta', 'f')
        self.addBranch('hlt_sigmaIEtaIEtaNoise', 'f')
        self.addBranch('hlt_ecalPFIsol', 'f')
        self.addBranch('hlt_hcalPFIsol', 'f')
        self.addBranch('hlt_trkIsol', 'f')
        self.addBranch('hlt_trkChi2', 'f')
        self.addBranch('hlt_trkMissHits', 'f')
        self.addBranch('hlt_trkNrLayerIT', 'i')
        self.addBranch('hlt_trkValidHits', 'i')
        self.addBranch('hlt_hcalHForHoverE', 'f')
        self.addBranch('hlt_invESeedInvP', 'f')
        self.addBranch('hlt_invEInvP', 'f')
        self.addBranch('hlt_trkDEta', 'f')
        self.addBranch('hlt_trkDEtaSeed', 'f')
        self.addBranch('hlt_trkDPhi', 'f')
        self.addBranch('hlt_pms2', 'f')
        self.addBranch('hlt_dr', 'f')
        self.addBranch('hlt_dphi', 'f')
        self.addBranch('hlt_deta', 'f')
        self.addBranch('hlt_dxy', 'f')

        self.addBranch('isgjson', '?')
        self.addBranch('instL', 'f')
        self.addBranch('npu', 'i')


        for pt in [5,6,7,8]:
#            self.addBranch('doubleE' + str(pt),                  '?')
            self.addBranch('l1_doubleE' + str(pt),                  '?')
#            self.addBranch('doubleE' + str(pt) + '_nomatch',                  '?')



#        for pt in ptrange:
#            self.addBranch('doubleE' + str(pt),                  '?')
