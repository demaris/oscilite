import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import time
from scipy.ndimage import zoom
from diffusiveCML import DiffusiveCML
from initCML import *
from analysisCML import *

class ControllerCML:

    def __init__(self,cml,stats,im=None,histax=None,scale=4,drawmod=1,drawtype='state',wait=0,localIter=1,last_render_time=0):
       self.cml=cml
       self.stats=stats
       self.im=im
       self.histax=histax
       self.scale=scale
       self.drawmod=drawmod
       self.drawtype=drawtype
       self.wait=wait
       self.localIter=localIter
       self.last_render_time=last_render_time

# create window frame
fig, ax = plt.subplots()
ax.patch.set_facecolor('black')
plt.subplots_adjust(left=0.25, bottom=0.45)

sidelen=100
# scale will scale the matrix in graphics window
scale=4
# establish the number of bins (partitions in the state space
binCount=64
# drawmod draws only on modulo drawmod lattice iterations. This is useful to limit framerate or find a cycle avoiding flicker

# Various initial lattice styles
cmlInit=''
#initLattice=imageCML('./shapes150.jpg')
# primes will be slow for a big lattice
#initLattice=primesSquare(sidelen)
#initLattice=np.rot90(initLattice,2)
#initLattice=randomPing(sidelen,sidelen,scaleFactor=0.0)
#initLattice=randomCML(sidelen,sidelen)
# following needs magic library installed
#initLattice=magicSquare(sidelen)
# create controller class to pass to update
initLattice=randbin(sidelen,sidelen)

# predefined kernels:  'symm4','symm8','asymm', 'magic11

cml = DiffusiveCML(initLattice,kern='symm4',a=1.75,gl=0.1,gg=0.05,wait=0);

stats = AnalysisCML(initLattice,binSpec=binCount)

# show initial state and create image
if cmlInit == 'image':
    # translate float in interval -1,1 to colormap
    llshow = (cml.matrix + 1) * 128
else:
    llshow = zoom(((cml.matrix) + 1) * 128, scale, order=2)
im=plt.imshow(llshow)
# gist_heat, flag, gist_ncar,
plt.set_cmap('Spectral')
cml_ax=plt.axis([0, sidelen, sidelen, .001])

# create control sliders for alpha, local and global coupling
axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axalpha = plt.axes([0.2, 0.2, 0.65, 0.03],axisbg=axcolor)
axlocal = plt.axes([0.2, 0.25, 0.65, 0.03],axisbg=axcolor)
axglob = plt.axes([0.2, 0.3, 0.65, 0.03],axisbg=axcolor)
histax = plt.axes([0.2, 0.04, 0.65, 0.1])
histax.set_ylim(0.0,0.5)
histax.set_xlim(-1.0,1.0)
# creat parameter controll sliders
alpha = Slider(axalpha, 'alpha', 0.1, 2.0, valinit=cml.a)
local = Slider(axlocal, 'local coupling', 0.01, 0.5, valinit=cml.gl)
glob = Slider(axglob, 'global coupling', 0.01, 0.2, valinit=cml.gg)

# callback for when any slider changes
def update_parms(val):
    cml.a = alpha.val
    cml.gl = local.val
    cml.gg = glob.val
    fig.canvas.draw_idle()

# callback for parameter sliders
alpha.on_changed(update_parms)
local.on_changed(update_parms)
glob.on_changed(update_parms)

#rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
rax = plt.axes([0.025, 0.5, 0.15, 0.15])
radio = RadioButtons(rax, ('state', 'spin'), active=0)

def visual_func(label):
    global drawtype
    drawtype=label
    fig.canvas.draw_idle()
radio.on_clicked(visual_func)

# create controller to encapsulate various rendering parameters
controller=ControllerCML(cml,stats,im=im,histax=histax,scale=4,drawmod=2,drawtype='state',wait=0,localIter=1,last_render_time=0)

def update_cml(cntl):
    global cml, stats, alpha
    global im, histax, drawmod, last_render_time

    # perform diffusion and reaction steps
    cntl.cml.iterate()
    # print inter to console (should make a label perhaps
    if cntl.cml.iter % 500==0: print cntl.cml.iter
    # compute measurements on the whole lattice
    cntl.stats.update(cml.matrix,cml.iter)
    # control example by spins (ref.  J.C. Perez, The New Way of Artificial Intelligence
    # experiment with spin control - number of spin transitions > threshold, or else decrease alpha
    # if a is chaotic, it will search and find a more stable (but probably still chaotic) value reducing spin transitions
    # This can be considered a form of unsupervised learning and essentially a similarity preserving hash into the partition space
    """
    if stats.spinTrend>100:
        print "reducing alpha"
        cml.a=cml.a-.05
        alpha.val=cml.a
    # stop diffusion
    if cml.iter==1000:
        cml.kernType='zero'
        cml.kernelUpdate()
        #cml.a=1.7
        print "stopping diffusion and quenching a= ",cml.a
    """
    # keep track of computation and render performance
    if (cntl.cml.iter>1 and cntl.cml.iter % cntl.drawmod==0):
        #calculating fps with goal of 30fps
        current_time = time.time()
        render_time = current_time - cntl.last_render_time
        cntl.last_render_time = current_time
        fps = round(1/render_time)
        if cntl.cml.iter==10:
            print(fps,'fps')
        # if an image is big, don't do the scaling with scipy zoom but rather use it direct.
        if cmlInit=='image':
            llshow=(cntl.cml.matrix+1)*128
        else:
            #print 'scaling'
            if cntl.drawtype=='state': llshow=zoom(((cntl.cml.matrix)+1)*128, cntl.scale, order=2)
            if cntl.drawtype=='spin': llshow=zoom(((cntl.stats.spin)+1)*128, cntl.scale, order=2)
        ## Display the data
        # zero out background in image to prevent flashing
        llshow[np.where(llshow==llshow[0,0])]=0
        cntl.im.set_data(llshow)

        width=2.0 / len(cntl.stats.edges)
        # set barplot centers; next line only works with uniform bins
        center= (cntl.stats.edges[:-1]+cntl.stats.edges[1:]) / 2
        cntl.histax.cla()
        cntl.histax.autoscale(enable=False)
        cntl.histax.bar(center,stats.bins,width=width,align='center')
        cntl.im.figure.canvas.draw()
        cntl.histax.figure.canvas.draw()


# Create a new timer object. Set the interval to number of milliseconds
# (1000 is default) and tell the timer what function should be called.
timer = fig.canvas.new_timer(interval=1)
timer.add_callback(update_cml,controller)
# pass controller instance to avoid global vars in update_cml
#timer.add_callback(update_cml,cml_controller)
timer.start()
plt.show()





