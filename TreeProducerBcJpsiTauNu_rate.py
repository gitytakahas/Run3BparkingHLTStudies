import ROOT
import math 
from TreeProducerCommon import *


class TreeProducerBcJpsiTauNu_rate(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, **kwargs):
        print('TreeProducerBsTauTau_rate is called for', name)
        super(TreeProducerBcJpsiTauNu_rate, self).__init__(name,**kwargs)

        for pt in l1_ptrange:
            self.addBranch('l1_doubleE' + str(pt),                  '?')
#            self.addBranch('doubleE' + str(pt) + '_nomatch',                  '?')

        for pt in hlt_ptrange:
            self.addBranch('doubleE' + str(pt),                  '?')

        self.addBranch('isgjson', '?')
        self.addBranch('instL', 'f')
        self.addBranch('npu', 'i')
