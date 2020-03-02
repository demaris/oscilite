__author__ = 'david'

"""
The competetive CML class:  responsible for initializing defaults, iterating and managing the lattice (matrix)
"""

from numpy import *
from scipy.signal import convolve2d

class CompetitiveCML:

    def __init__(self, lattice,l=5.0,a=0.26):
        """
        lattice: initial state of lattice
        l:
        a:
        """
        #global matrix, num_cells
        self.matrix=lattice
        self.computeMod=1
        self.numCells = size(lattice,0) * size(lattice,1)
        self.l=l
        # replicate an array to conform to size of lattice
        
        self.a=a
        self.k=1.0
        # make a coupling kernel to sum the surround
        self.ckern=array([(0,self.k,0),(self.k,0,self.k),(0,self.k,0)])
        self.iter=0


    def iterate(self):
        """
        see arxiv.org/pdf/1204.2463.pdf
        """

        self.iter += 1
        # ricker
        # x(t+1)=l x(t)exp(-(x(t))+a*sumSurround(x(t)))
        # ccml
        # x(t+1) x(t)*ricker[
        #convolution with 1.0 symmetric kernel sums map
        sum_surround=convolve2d(self.matrix, self.ckern, mode='same', boundary='wrap')
        #
        temp=self.matrix+self.a*sum_surround
        #self.matrix = self.l* (1+power(intermediate,self.b))

        self.matrix = self.l*self.matrix*exp(-temp)
