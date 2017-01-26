# -*- coding: utf-8 -*-
"""
Hybrid CML displayed via a RawImageWidget in QT
For glitch free sound, you want a long buffer in pyo and order 1 on numpy zoom
"""

import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.ndimage import zoom
from diffusiveCML import DiffusiveCML
from initCML import *
from analysisCML import *
# enable animation
plt.ion()
sidelen=100
# scale will scale the matrix in graphics window
scale=4
# establish the number of bins (partitions in the state space
binCount=32
# create and show window
# drawmod draws only on modulo drawmod lattice iterations. This is useful to limit framerate or find a cycle avoiding flicker
drawmod=2
# Various initial lattice styles
cmlInit=''
#initLattice=imageCML('./shapes150.jpg')
#initLattice=primesSquare(sidelen)
#initLattice=magicSquare(sidelen)
#initLattice=np.rot90(initLattice,2)
#initLattice=imageCML('./sri_mandala.jpg');
#initLattice=randomPing(sidelen,sidelen,scaleFactor=0.0)
initLattice=randomCML(sidelen,sidelen)
#initLattice=magicSquare(sidelen)
#initLattice=primesSquare(sidelen)
#initLattice=randbin(sidelen,sidelen)
#print initLattice
# wait variable can slow things down by running a counter inside
#why do we have a, gl, and gg in here as well as initCML?
#cml = DiffusiveCML(initLattice,kern='magic11',a=1.8,gl=0.03,gg=0.01,wait=0)
cml = DiffusiveCML(initLattice,kern='symm4',a=1.5,gl=0.7,gg=0.00,wait=0);
#cml = DiffusiveCML(initLattice)
stats = AnalysisCML(initLattice,binSpec=binCount)
#cml = DiffusiveCML(initLattice,kern='magic11')
#cml = CompetitiveCML(initLattice)


last_render_time = 0
# show initial state and create image
if cmlInit == 'image':
    llshow = (cml.matrix + 1) * 128
else:
    llshow = zoom(((cml.matrix) + 1) * 128, scale, order=3)
im=plt.imshow(llshow)

while True:

    # perform diffusion and reaction steps
    cml.iterate()
    print cml.iter
    # do some measurements on the whole lattice
    stats.update(cml.matrix,cml.iter)
    # control example by spins (ref.  J.C. Perez, The New Way of Artificial Intelligence
    # experiment with spin control - number of spin transitions > threshold, or else decrease alpha
    # if a is chaotic, it will search and find a more stable (but probably still chaotic) value reducing spin transitions
    # This can be considered a form of unsupervised learning and essentially a similarity preserving hash into the partition space
    if stats.spinTrend>100:
        print "reducing alpha"
        cml.a=cml.a-.05
    # stop diffusion
    if cml.iter==1000:
        cml.kernType='zero'
        cml.kernelUpdate()
        #cml.a=1.7
        print "stopping diffusion and quenching a= ",cml.a
    # keep track of computation and render performance
    if (cml.iter>1 and cml.iter % drawmod==0):
        #calculating fps with goal of 30fps
        current_time = time.time()
        render_time = current_time - last_render_time
        last_render_time = current_time
        fps = round(1/render_time)
        if cml.iter==10:
            print(fps,'fps')
        # if an image is big, don't do the scaling with scipy zoom but rather use it direct.
        if cmlInit=='image':
            llshow=(cml.matrix+1)*128
        else:
            #print 'scaling'
            llshow=zoom(((cml.matrix)+1)*128, scale, order=3)
            #llshow=zoom(((stats.spin)+1)*128, scale, order=3)
        ## Display the data
        # zero out background in image to prevent flashing
        llshow[np.where(llshow==llshow[0,0])]=0
        # following should be faster but I don't see anything
        #im.set_data(llshow)
        #plt.draw()
        plt.imshow(llshow)
        # draw histogram



