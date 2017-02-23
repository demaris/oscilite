from numpy import *
from scipy.stats import entropy
class AnalysisCML:

    def __init__(self, initLattice,doSpins=1,doEntropy=1,binSpec=16):
        """
        take cml lattisce initial state as parameter, and controls for what stats to gather
        """
        #global matrix, numCells
        self.doEntropy=doEntropy
        self.doSpins=doSpins
        self.binSpec=binSpec
        self.last=initLattice
        self.spin=0
        # scalar counts of spin transitions in current and prev iteration, for control
        self.spinTrans=0
        self.lastSpinTrans=0
        self.spinTrend=0
        # value of entropy
        self.entropy=0
        self.cells = size(initLattice,0) * size(initLattice,1)
        self.bins=0
        self.edges=0
        self.cumBins=0
        self.count=0


    def update(self,lattice,iter,histrange='KKfull'):
        """
        Update stats for CML matrix
        """
        # compute all the stats
        # bins histogram
        # use temp vars for unscaled
        if histrange=='KKfull':
            bins,self.edges=histogram(lattice,bins=self.binSpec,range=(-1.0,1.0))
            self.cumBins=self.cumBins+self.bins
                # don't need normalization with density=True

        self.bins= bins/float(self.cells)

        # shannon entropy over binsh
        if self.doEntropy:
            self.entropy=entropy(self.bins)

        # spins and transitions
        if self.doSpins:
            self.spin=lattice-self.last
            self.spin[where(self.spin>=0)]=1.0
            self.spin[where(self.spin<0)]=-1.0

        # the next block saves current state needed to compute spin and spin transitions
        if self.count>=1:
            if self.doSpins:
                self.last=lattice

                if self.count>1:
                    self.spinTrans=len(where(self.spin!=self.lastSpin)[0].tolist())
                    self.spinTrend=self.spinTrans-self.lastSpinTrans
                self.lastSpin=self.spin
                self.lastSpinTrans=self.spinTrans

        self.count = self.count + 1



