from numpy import *
from scipy.stats import entropy
class AnalysisCML:

    def __init__(self, init_lattice, do_spins=1, do_entropy=1, bin_spec=16):
        """
        take cml lattice initial state as parameter, and controls for what stats to gather
        """
        #global matrix, num_cells
        self.do_entropy=do_entropy
        self.do_spins=do_spins
        self.bin_spec=bin_spec
        self.last=init_lattice
        self.spin=0
        # scalar counts of spin transitions in current and prev iteration, for control
        self.last_spin=0
        self.spin_trans=0
        self.last_spin_trans=0
        self.spin_trend=0
        # value of entropy
        self.entropy=0
        self.cells = size(init_lattice, 0) * size(init_lattice, 1)
        self.bins=0
        self.edges=0
        self.cum_bins=0
        self.iter_count=0

    def update(self,lattice,histrange='plusminus'):
        """
        Update stats for CML matrix
        histrange plusminus indicates that the Kaneko logistic map form is used with range -1:1
        One might provide another map with a different range and adjust bin range accordingly
        """

        if histrange=='plusminus':
            bins,self.edges=histogram(lattice,bins=self.bin_spec,range=(-1.0,1.0))
            self.cum_bins= self.cum_bins + self.bins

        self.bins= bins/float(self.cells)

        # shannon entropy over binsh
        if self.do_entropy:
            self.entropy=entropy(self.bins)

        # spins and transitions
        if self.do_spins:
            self.spin=lattice-self.last
            # if trending up or reversing direction to up
            self.spin[where(self.spin>=0)]=1.0
            # if trending down or reversing direction to down
            self.spin[where(self.spin<0)]=-1.0

        # the next block saves current state needed to compute spin and spin transitions
        if self.iter_count>=1:   # must be
            if self.do_spins:
                self.last=lattice
                if self.iter_count>1:
                    self.spin_trans=len(where(self.spin != self.last_spin)[0].tolist())
                    self.spin_trend= self.spin_trans - self.last_spin_trans
                self.last_spin=self.spin
                self.last_spin_trans=self.spin_trans
        self.iter_count = self.iter_count + 1



