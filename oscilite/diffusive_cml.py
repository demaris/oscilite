__author__ = 'demaris'

"""
The diffusive CML class:  responsible for initializing defaults, iterating and managing the matrix
if the global coupling gg is set to zero, the dynamics are strictly the classic coupled map lattice.
If gl is set to zero, the dynamics are the globally coupled map (similar to logistic map on small world network
To give some guidance of interesting parameters, a java CML simulation suggests some parameters for 1D
http://brain.cc.kogakuin.ac.jp/~kanamaru/Chaos/e/CMLGCM
inspired by Kunihiko Kaneko and Ichiro Tsuda "Complex Systems: Chaos and Beyond, a Constructive Approach with
Applications in the Life Sciences", Springer Verlag, 2000
The following parameters result in different pattern formation regimes or phases, however they are developed for a 1-D model,
with a symmetric kernel. They may still be useful as guidelines to search for phase regime transitions, as bounds or
as starting points for some parameter search.

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
import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage.morphology import generate_binary_structure

class DiffusiveCML:

    def __init__(self, lattice, kern='symm4', gl=0.5, gg=0.0, a=1.47, a_spread=0.0, local_iter=1, \
                 use_pinned=False, pinned_mask='', use_activation=False):
        """

        :param lattice: initial state of matrix
        :param kern: coupling kernel drawn from predeclared values in this function, typically with orientation sensitivity
        :param gl: local coupling (gamma in Kaneko papers), range 0-1  (typically .5 or less)
        :param gg: global coupling range 0-1 (typically .12 or less)
        :param a: alpha nonlinearity control parameter; low (<1.5) is period doubling cascade, high (<2.0) chaotic with windows
        :param a_spread: a gradient over alpha is defined in the vertical direction
        :param local_iter: number of diffusion cycles before each reaction.  Classical cml is 1 diffusion cycle then one nonlinearity cycle.
        :param use_pinned:  a mask is applied after the nonlinear map, allowing sites to be pinned to a fixed value serving as control
        :param pinned_mask: the mask to be used for pinning
        :param use_activation: a spreading activation mask with its own dynamical evolution is used after the lattice and only active sites are updated
        """
        self.matrix = lattice

        self.num_cells = np.size(lattice, 0) * np.size(lattice, 1)

        self.kern_type=kern
        # activation array is used to control which sites are transmitted to next generation, with spreading activation
        # from any seed sites
        self.activation = np.zeros_like(lattice)
        self.use_activation=use_activation
        # non-zero values in mask array are used to pin a lattice site to a high value, controlling the oscillations

        self.pinned_mask=pinned_mask
        self.use_pinned=use_pinned

        # local_iter allows multiple diffusion cycles before nonlinear map is applied
        self.local_iter=local_iter
        # defaults tested here should match argument defaults
        self.gl=gl
        self.gg=gg
        self.a=a
        self.a_spread=a_spread
        if self.kern_type == 'symm4':
            self.cc=self.gl/4
            self.dkern=np.array([(0,self.cc,0),(self.cc,0,self.cc),(0,self.cc,0)])
        elif self.kern_type == 'symm4n':
                self.cc = self.gl / 4
                self.dkern = np.array([(0, self.cc, 0), (self.cc, -1.0, self.cc), (0, self.cc, 0)])
        elif self.kern_type == 'symm8':
            self.cc=self.gl/8
            self.dkern=np.array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])
        elif self.kern_type == 'asymm':
            self.cc=self.gl/5
            self.dkern=np.array([(0,self.cc,0),(self.cc,0,self.cc),(0,self.cc,self.cc)])
        elif self.kern_type == 'magic11':
            self.gl=0.404040404
            self.dkern=[[ 0.08080808,  0.01010101,  0.06060606],
                         [ 0.03030303,  0.0,  0.07070707],
                         [ 0.04040404,  0.09090909,  0.02020202]]
        elif self.kern_type == 'pinwheel':
            self.cc=self.gl/20
            self.dkern=[[0, self.cc, self.cc, 0, 0, self.cc, 0],
                        [self.cc, 0, self.cc, 0, self.cc, 0, self.cc],
                        [0, self.cc, 0, 0, 0, self.cc, self.cc],
                        [0, 0, 0, 0, 0, 0, 0,],
                        [self.cc, self.cc, 0, 0, 0, self.cc, 0],
                        [self.cc, 0, self.cc, 0, self.cc, 0, self.cc],
                        [0, self.cc, 0, 0, self.cc, self.cc, 0]]
        # this is the structuring kernel for dilation, which implements a spreading activation from active sites
        # second arg is connnectivity, with 2 giving a 3x3 all ones structring element
        self.activation_struct=generate_binary_structure(2, 2)
        self.iter=0
        if self.pinned_mask=='':
            self.pinned_mask=np.zeros_like(lattice)

    def alpha_update(self):
        # check for reasonable bounds on a_spread; non-negative, < stable range in map
        if self.a_spread < 0.0 or self.a_spread > 2.0:
            print ('a_spread value outside range 0:2')
            return
        if self.a_spread> 0.0:
            # generate a as matrix rather than scalar
            if np.isscalar(self.a):
                # it's scalar
                low=self.a
            else:
                # was already a vector alpha
                low=self.a[0,0]

            high=low+self.a_spread
            if high > 2.0:
                # print('alpha spread high ',high,' saturated, clamping max to 2.0')
                high=2.0
            a_vals =  np.linspace(low,high,np.size(self.matrix[0]))
            a_rep = np.tile(a_vals,(np.size(a_vals),1))
            self.a=np.rot90(a_rep,1)
        else:
            # reset back to scalar from corner value of array
            if np.isscalar(self.a):
                self.a=self.a
            else:
                self.a=self.a[0,0]
        # print self.a

    def kernel_update(self):
        """"
        Convolution kernel for diffusive coupling can be updated dynamically
        """
        """
        Convolution kernel for diffusive coupling can be updated dynamically
        """
        if self.kern_type == 'symm4':
            # kernel has some edge sensitivity
            self.cc=self.gl/4
            self.dkern=np.array([(0,self.cc,0.0),(self.cc,0,self.cc),(0,self.cc,0)])

        elif self.kern_type == 'symm8':
            # symmetric in all directions
            self.cc=self.gl/8
            self.dkern=np.array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])
        elif self.kern_type == 'asymm':
            self.cc=self.gl/5
            self.dkern=np.array([(0,self.cc,0.0),(self.cc,0,self.cc),(0,self.cc,self.cc)])
        elif self.kern_type == 'magic11':
            self.gl=0.404040
            self.dkern=[[ 0.08080808,  0.01010101,  0.06060606],
                         [ 0.03030303,  0.0,  0.07070707],
                         [ 0.04040404,  0.09090909,  0.02020202]]
        elif self.kern_type == 'zero':
            self.cc=0.0
            self.dkern=np.array([(self.cc,self.cc,self.cc),(self.cc,0,self.cc),(self.cc,self.cc,self.cc)])

        elif self.kern_type == 'pinwheel':
            self.cc = self.gl / 20
            self.dkern = [[0, self.cc, self.cc, 0, 0, self.cc, 0],
                      [self.cc, 0, self.cc, 0, self.cc, 0, self.cc],
                      [0, self.cc, 0, 0, 0, self.cc, self.cc],
                      [0, 0, 0, 0, 0, 0, 0, ],
                      [self.cc, self.cc, 0, 0, 0, self.cc, 0],
                      [self.cc, 0, self.cc, 0, self.cc, 0, self.cc],
                      [0, self.cc, 0, 0, self.cc, self.cc, 0]]

    def iterate(self):
        """
        Iterate / convolve the matrix using image toolkit, then apply nonlinear map with mask
        """
        self.iter += 1

        # update activation map if in use
        if self.use_activation:
            self.activation = binary_dilation(self.activation, structure=self.activation_struct).astype(
                self.activation.dtype)

            # scale lattice by activation map before nonlinear map and re-apply pinMask if active; otherwise
            # pin sites from last step would be lost by activation scaling in next step
            self.matrix = self.matrix * self.activation

            if self.use_pinned:
                # apply mask values where non-zero
                pin_points = np.where(self.pinned_mask != 0)
                self.matrix[pin_points] = self.pinned_mask[pin_points]
        # diffusion
        for i in range(self.local_iter):
            diff=convolve2d(self.matrix, self.dkern, mode='same', boundary='wrap')
            # scale before adding to keep value in <-1,+1> bounds
            diff_scaled=((1-self.gl) * self.matrix + diff)

        # apply the map after scaling and accounting for any global coupling
        if self.gg > 0:
            self.matrix = (1-self.gg) * (1- (self.a * (diff_scaled**2))) + (self.gg / self.num_cells) * sum(self.matrix)
        else:
            self.matrix = (1- (self.a * (diff_scaled**2)))
        # mask by activation map
        # apply pin mask if enabled
        if self.use_pinned:
            # apply mask values where non-zero
            pin_points=np.where(self.pinned_mask != 0)
            self.matrix[pin_points]=self.pinned_mask[pin_points]