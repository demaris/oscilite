import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import time
from scipy.ndimage import zoom
from competitive_cml import CompetitiveCML
import init_cml
from analysis_cml import AnalysisCML

class ViewCML:

    def __init__(self, cml, stats, im=None, axhist=None, scale=4, draw_mod=1, draw_type='state', wait=0, local_iter=1, last_render_time=0):
       self.cml=cml
       self.stats=stats
       self.im=im
       self.axhist=axhist
       self.scale=scale
       self.draw_mod=draw_mod
       self.draw_type=draw_type
       self.wait=wait
       self.local_iter=local_iter
       self.last_render_time=last_render_time

    def visual_func(self,label):
        view.draw_type = label
        view.im.figure.canvas.draw_idle()


    def update(self):
        # perform diffusion and reaction steps
        view.cml.iterate()
        # print inter to console (should make a label perhaps
        if view.cml.iter % 10 == 0: print view.cml.iter
        # compute measurements on the whole lattice
        view.stats.update(cml.matrix)

        # keep track of computation and render performance
        if (view.cml.iter > 1 and view.cml.iter % view.draw_mod == 0):
            # calculating fps with goal of 30fps
            current_time = time.time()
            render_time = current_time - view.last_render_time
            view.last_render_time = current_time
            fps = round(1 / render_time)
            if view.cml.iter == 10:
                print(fps, 'fps')
            # if an image is big, don't do the scaling with scipy zoom but rather use it direct.
            if zoom_style == 'image':
                llshow = (view.cml.matrix + 1) * 128
                # zero out background in image to prevent flashing
                # llshow[np.where(llshow == llshow[0, 0])] = 0
            else:
                # print 'scaling'
                if view.draw_type == 'state':
                    # competitive values always positive so no shift required
                    drawoffset = 0
                    # compute normalization factor to put data in range 0-127
                    norm_val = np.max(cml.matrix)
                    #print norm_val
                    llshow = zoom(((cml.matrix) / norm_val + drawoffset) * 128, scale, order=2)
                    #print llshow
                elif view.draw_type == 'spin': llshow = zoom(((view.stats.spin) + 1) * 128, view.scale, order=0)
            ## Display the data
            view.im.set_data(llshow)

            width = 2.0 / len(view.stats.edges)
            # set barplot centers; next line only works with uniform bins
            center = (view.stats.edges[:-1] + view.stats.edges[1:]) / 2
            view.axhist.cla()
            view.axhist.autoscale(enable=False)
            view.axhist.bar(center, stats.bins, width=width, align='center')
            view.im.figure.canvas.draw()
            view.axhist.figure.canvas.draw()


# create window frame
fig, ax = plt.subplots()
ax.patch.set_facecolor('black')
plt.subplots_adjust(left=0.25, bottom=0.45)
plt.suptitle('Excitation Level, random re-seed on parameter changes')
fig.canvas.set_window_title("Competitive Coupled Map Lattice")

sidelen=100

# establish the number of bins (partitions in the state space
bin_count=64
# draw_mod draws only on modulo draw_mod lattice iterations. This is useful to limit framerate or find a cycle avoiding flicker
scale=4
# Various initial lattice styles
zoom_style= ''
#initLattice=imageCML('./shapes150.jpg')
# primes will be slow for a big lattice
#initLattice=primesSquare(sidelen)
#initLattice=np.rot90(initLattice,2)
#initLattice=randomPing(sidelen,sidelen,scaleFactor=0.0)
#initLattice=randomCML(sidelen,sidelen)
# following needs magic library installed
#initLattice=magicSquare(sidelen)

init_lattice=init_cml.random_cml_pos(sidelen,sidelen)
zoom_style= 'rand'
# predefined kernels:  'symm4','symm8','asymm', 'magic11

cml = CompetitiveCML(init_lattice,l=5,a=0.26);
stats = AnalysisCML(init_lattice, bin_spec=bin_count)

# show initial state and create image
if zoom_style == 'image':
    # translate float in interval -1,1 to colormap
    llshow = (cml.matrix + 1) * 128
else:
    # competitive values always positive so no shift required
    drawoffset=0
    # compute normalization factor to put data in range 0-127
    norm_val=np.max(cml.matrix)
    llshow = zoom(((cml.matrix) / norm_val + drawoffset) * 128, scale, order=2)
im=plt.imshow(llshow)
# gist_heat, flag, gist_ncar,
plt.set_cmap('gist_heat')

cml_ax=plt.axis([0, sidelen, sidelen, .001])

# create control sliders for alpha, local and global coupling
axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axalpha = plt.axes([0.2, 0.2, 0.65, 0.03],axisbg=axcolor)
axlocal = plt.axes([0.2, 0.25, 0.65, 0.03],axisbg=axcolor)
axhist = plt.axes([0.2, 0.04, 0.65, 0.1])
axhist.set_ylim(0.0,1.0)
axhist.set_xlim(0,1)
# creat parameter contrll sliders
lmda = Slider(axalpha, 'lambda', 0.1, 5.0, valinit=cml.l)
exp = Slider(axlocal, 'exponent a', 0.01,1.0, valinit=cml.a)

# callback for when any slider changes
def update_parms(val):
    cml.l = lmda.val
    cml.a = exp.val
    # also reinit, since Competitive will have reached a stable state
    cml.matrix=init_cml.random_cml_pos(sidelen,sidelen)
    fig.canvas.draw_idle()

# callback for parameter sliders
lmda.on_changed(update_parms)
exp.on_changed(update_parms)

#rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
rax = plt.axes([0.025, 0.5, 0.15, 0.15])
radio = RadioButtons(rax, ('state', 'spin'), active=0)

# create controller to encapsulate various rendering parameters
view=ViewCML(cml, stats, im=im, axhist=axhist, scale=4, draw_mod=1, draw_type='state', wait=0, local_iter=1, last_render_time=0)

radio.on_clicked(view.visual_func)

# Create a new timer object. Set the interval to number of milliseconds
# (1000 is default) and tell the timer what function should be called.
timer = fig.canvas.new_timer(interval=1)
timer.add_callback(view.update)
# pass controller instance to avoid global vars in update_cml
#timer.add_callback(update_cml,cml_controller)
timer.start()
plt.show()





