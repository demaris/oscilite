__author__ = 'demaris'

"""
The diffusive CML class:  responsible for initializing defaults, iterating and managing the matrix
if the global coupling gg is set to zero, the dynamics are the coupled map lattice.
The following parameters result in different pattern formation regimes or phases:
frozen chaos: a<1.5
The sites are divided into the clusters with various sizes. The divided pattern can be regarded as an attractor.
With the different initial condition, the pattern of the clusters will be changed, thus it is known that the many attractors co-exist in the system.
The number of attractors increases (at least) exponentially with the increase of N.

pattern selection: e.g., a=1.71, gl=0.4 gg=0
The sites are divided into the clusters with almost the same sizes. The size of the clusters depends on the parameters.

chaotic Brownian motion of defect: e.g., a=1.85, gl=0.1
In the case of the above pattern selection, some phase defects sometimes remain in the system.
These defects fluctuates chaotically like Brownian motion.

defect turbulence: e.g., a=1.895, gl=0.1, gg=0
The many defects are generated and collide each other like the turbulence.

spatiotemporal intermittency: e.g., a=1.75, gl=0.3
The each site transits between the coherent state and the chaotic state intermittently.

fully developed spatiotemporal chaos: e.g., a=2.0, gl=0.3
Almost all the sites oscillate chaotically independently.

traveling wave: e.g., a=1.47, gl=0.5, gg=0
The wave of clusters travels at the very low speed.
"""


from numpy import *
from numpy.random import rand
from scipy.signal import convolve2d


class DiffusiveCML:

    def __init__(self, lattice,kern='symm4',dynamicsParms='',gl=0.4,gg=0.0,a=1.74,name='name',wait=0,localIter=1):
        """
        Parameters:
        lattice             initial state of matrix
        kern                coupling kernel from known sets
        dynamicsParms       dynamical evolution parms from known sets
        gl                  local coupling, range 0-1  (typically .5 or less)
        gg                  global coupling range 0-1 (typically .3 or less)
        a                   alpha nonlinearity parameter; low is periodic, high (<2.0) chaotic
        name                pass in a name to save parameters in a file when you have something interesting
        wait                cycles to count before computation to slow down CML graphics for small arrays

        """

        self.computeMod=20
        self.matrix = lattice
        self.numCells = size(lattice,0) * size(lattice,1)
        self.kernType=kern
        # wait enables to run a wait loop before iterating to slow down for human visualization
        self.wait=wait
        # localIter can run multiple diffusion cycles before nonlinear map
        self.localIter=localIter
        # good parms gg=.1,gl=.4,a=1.7
        # gg.05, same
        # gg 0.05, gl 0.5
        if dynamicsParms=='':
            self.gl=0.4
            self.gg=0.1
            self.a=1.74
        if dynamicsParms=='travWave':
            self.a=1.47
            self.gl=0.5
            self.gg=0
        if dynamicsParms=='chaoticBrown':
            self.a=1.85
            self.gl=0.1
            self.gg=0
        # defaults tested here should match argument defaults
        if gl!=0.4: self.gl=gl
        if gg>=0.0: self.gg=gg
        if a!=1.4: self.a=a

        if self.kernType == 'symm4':
            self.cc=self.gl/4
            self.dkern=array([(0,self.cc,0),(self.cc,0,self.cc),(0,self.cc,0)])
        elif self.kernType == 'symm8':
            self.cc=self.gl/8
            self.dkern=array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])
        elif self.kernType =='asymm':
            self.cc=self.gl/5
            self.dkern=array([(0,self.cc,0),(self.cc,0,self.cc),(0,self.cc,self.cc)])
        elif self.kernType == 'magic11':
            self.gl=0.404040404
            self.dkern=[[ 0.08080808,  0.01010101,  0.06060606],
                         [ 0.03030303,  0.0,  0.07070707],
                         [ 0.04040404,  0.09090909,  0.02020202]]
        self.iter=0

    def kernelUpdate(self):
        if self.kernType == 'symm4':
            self.cc=self.gl/4
            self.dkern=array([(0,self.cc,0.0),(self.cc,0,self.cc),(0,self.cc,0)])
        elif self.kernType == 'symm8':
            self.cc=self.gl/8
            self.dkern=array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])
        elif self.kernType =='asymm':
            self.cc=self.gl/5
            self.dkern=array([(0,self.cc,0.0),(self.cc,0,self.cc),(0,self.cc,self.cc)])
        elif self.kernType == 'magic11':
            self.gl=0.404040
            self.dkern=[[ 0.08080808,  0.01010101,  0.06060606],
                         [ 0.03030303,  0.0,  0.07070707],
                         [ 0.04040404,  0.09090909,  0.02020202]]
        elif self.kernType == 'zero':
            self.cc=0.0
            self.dkern=array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])


    def iterate(self):
        """
        Iterate / convolve the matrix
        """
        self.iter += 1
        if self.wait>0:
            count=0
            while count<self.wait:
                count+=1

        # diffusion
        #react=(1-self.a*self.matrix**2)

        for i in range(self.localIter):
            diff=convolve2d(self.matrix, self.dkern, mode='same', boundary='wrap')
            # scale before adding to keep value in <-1,+1> bounds
            diffScaled=((1-self.gl) * self.matrix + diff)

        # apply the map after scaling
        self.matrix = (1-self.gg) * (1- (self.a * (diffScaled**2))) + (self.gg/self.numCells) * sum(self.matrix)
        #self.matrix=(1-self.gg) * (1-self.gl*react)+diff + (self.gg/self.numCells) * sum(self.matrix)