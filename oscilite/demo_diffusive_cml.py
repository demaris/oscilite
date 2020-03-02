#import PyQt5
import matplotlib
#matplotlib.use("TkAgg")
#matplotlib.use("Qt5Agg");  matplotlib.rcParams['backend.qt5']='PySide2'
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

import time
from scipy.ndimage import zoom
from diffusive_cml import DiffusiveCML
import init_cml
from analysis_cml import AnalysisCML

"""
demo_diffusive is a demo instanatiating a ViewCML is a class to maintain various parameters for a demonstration and matplotlib controls.
Several features of oscillite not in classical CML literature are illustrated, including spins, pinning or 'crystal defect' sites, and activation.
"""
class ViewCML:
    """

    """
    def __init__(self,cml,stats,im=None,axhist=None,scale=4,drawmod=1,drawtype='state'):
        """

        :param cml: the previously created lattice to be displayed
        :param stats: the stats object created to monitor the lattice
        :param im: the image object to be updated after some number of steps (default 1)
        :param axhist: an axis to draw a histogram associated with the lattice
        :param scale: a scale value to increase the size of the lattice view by pixel replication and smoothing
        :param drawmod: a modulus which can be used to skip drawing every site in case of flickering dynamics
        :param drawtype: visualization type, one of 'state','spin','pinned',activation,
        """
        self.cml=cml
        self.stats=stats
        self.im=im
        self.axhist=axhist
        self.scale=scale
        # draw_mod draws only on modulo draw_mod lattice iterations, to limit framerate for small lattices,
        # or to avoid visualization flicker in a demo when you're in a state of high synchronization
        self.drawmod=drawmod
        # what variable is drawn in the image (state, spin)
        self.drawtype=drawtype
        # if spin_control True alpha is updated based on spin 'frustration', i.e. variance in local neighborhood
        self.spin_control=False
        self.lfo=False
        self.a_stack=[cml.a]
        self.gl_stack=[cml.gl]
        self.gg_stack=[cml.gg]
        self.cycle=5
        self.last_render_time=0


    def display_func(self,label):
        # set draw_type to match selected radio button
        view.drawtype = label
        view.im.figure.canvas.draw_idle()

    def actions_func(self, label):
        # control lattice evolution options
        # free and pinned are mutually exclusive, others can run with either free or pinned
        if label == 'free':
            # default, set others False if selected again
            self.cml.use_pinned=False
            self.spin_control=False
            self.lfo=False
            self.use_activation=False

        elif label =='pinned':
            self.cml.use_pinned=True
        elif label =='spin_control':
            self.spin_control=True
        elif label == 'LFO':
             self.lfo=True
        elif label == 'activation':
              self.cml.use_activation=True

        view.im.figure.canvas.draw_idle()

    # callback handling model parm slider changes
    def update_parms(self,val):

        self.cml.a = a_slider.val
        # insert only if different so that slider updates in lfo mode don't push same value
        if self.cml.a!=self.a_stack[0]:
            self.a_stack.insert(0,cml.a)
            # print "alpha stack ",self.a_stack,"cycled every ",self.cycle
        self.cml.gl = gl_slider.val
        if self.cml.gl!=self.gl_stack[0]:
            self.gl_stack.insert(0,cml.gl)
            # print "local coupling stack ", self.gl_stack, "cycled every ", self.cycle
        self.cml.gg = gg_slider.val
        if self.cml.gg!=self.gg_stack[0]:
            self.gg_stack.insert(0,cml.gg)
            # print "global coupling stack ", self.gg_stack, "cycled every ", self.cycle
        self.cml.a_spread=a_spread_slider.val
        self.cml.alpha_update()
        fig.canvas.draw_idle()

    def update(self):
        # perform diffusion and reaction steps
        view.cml.iterate()
        # print inter to console
        if view.cml.iter % 50 == 0: print(view.cml.iter)
        # compute measurements on the whole lattice
        view.stats.update(cml.matrix)
        # control example by spins (ref.  J.C. Perez, The New Way of Artificial Intelligence
        # experiment with spin control - number of spin transitions > threshold, or else decrease alpha
        # if a is chaotic, it will search and find a more stable (but probably still chaotic) value reducing spin transitions
        # This can be considered a form of unsupervised learning; essentially a similarity preserving hash into the partition space

        if view.spin_control:
            print ("spin_trend ",view.stats.spin_trend)
            if view.stats.spin_trend>100:

                view.cml.a=cml.a-.05
                a_slider.set_val(view.cml.a)
                print("reducing alpha", view.cml.a)
                fig.canvas.draw_idle()
            # stop diffusion
            if cml.iter==1000:
                cml.kern_type= 'zero'
                cml.kernel_update()
                #cml.a=1.7
                print("stopping diffusion and quenching a= ",cml.a)

        # if modulo cycle, switch to other of last two values slider values cl, alpha
        if view.lfo:
            # every 'cycle' iterations, swap top of stack values, set current alpha, update slider
            if view.cml.iter % view.cycle == 0:
                view.a_stack.insert(1,view.a_stack.pop(0))
                #print "cyclic reset cml.a",view.a_stack
                view.cml.a=view.a_stack[0]
                a_slider.set_val(view.cml.a)

                view.gl_stack.insert(1, view.gl_stack.pop(0))
                #print "cyclic reset cml.gl", view.gl_stack
                view.cml.gl = view.gl_stack[0]
                gl_slider.set_val(view.cml.gl)

                view.gg_stack.insert(1, view.gg_stack.pop(0))
                #print "cyclic reset cml.gg", view.gg_stack
                view.cml.gg = view.gg_stack[0]
                gg_slider.set_val(view.cml.gg)

                # since alpha may now be a spatial field rather than scalar, call the function that syncs that with
                # the switching alpha value
                view.cml.alpha_update()

        # keep track of computation and render performance
        if (view.cml.iter > 1 and view.cml.iter % view.drawmod == 0):
            # calculating fps
            current_time = time.time()
            render_time = current_time - view.last_render_time
            view.last_render_time = current_time
            fps = round(1 / render_time)
            if view.cml.iter % 50 == 0:
                print(fps, 'fps')
            # if an image is big, don't do the scaling with scipy zoom but rather use it direct.
            if zoom_style == 'image':
                llshow = (view.cml.matrix + 1) * 128
            else:
                # print 'scaling'
                if view.drawtype == 'state': llshow = zoom(((view.cml.matrix) + 1) * 128, view.scale, order=2)
                if view.drawtype == 'spin': llshow = zoom(((view.stats.spin) + 1) * 128, view.scale, order=2)
                if view.drawtype == 'pinned': llshow = zoom(((view.cml.pinned_mask)) * 128, view.scale, order=2)
                if view.drawtype == 'activation': llshow = zoom(((view.cml.activation)) * 128, view.scale, order=2)
            ## Display the data

            view.im.set_data(llshow)
            width = 2.0 / len(view.stats.edges)
            # set barplot centers; next line only works with uniform bins
            center = (view.stats.edges[:-1] + view.stats.edges[1:]) / 2
            view.axhist.cla()
            view.axhist.autoscale(enable=False)
            view.axhist.bar(center, stats.bins, width=width, align='center')
            ax.draw_artist(im)
            view.im.figure.canvas.draw_idle()
            view.im.figure.canvas.flush_events()
            #view.im.figure.canvas.draw()
            view.axhist.figure.canvas.draw()


# end view constructor and event handlers; create GUI
# create window frame
fig, ax = plt.subplots()
plt.suptitle('State of each lattice site')
fig.canvas.set_window_title("Hybrid Diffusive-Global Coupled Map Lattice")
ax.patch.set_facecolor('black')
plt.subplots_adjust(left=0.25, bottom=0.45)
# set the array size, square by default
# some display and setup code probably needs adjustment for non-square lattices
sidelen=100
# establish the number of bins (partitions in the state space) for the histogram display
bin_count=64
# controls scaling of graphics i.e how many pixels per site; rendering order on image draw controls blurring
scale=4
# var zoom_style controls zooming (pixel replication to make a bigger display), default is to zoom
zoom_style=''
# Various initial lattice styles.
#init_lattice=init_cml.image_cml('./MNISTdigits_crop.png')
# uncomment if initializing with an image
#zoom_style='image'
# primes will be slow for a big lattice
#init_lattice=init_cml.primes_square(sidelen)
#init_lattice=np.rot90(initLattice,2)
#init_lattice=init_cml.random_ping(sidelen,sidelen,scaleFactor=0.0)
#init_lattice=init_cml.random_cml(sidelen,sidelen)
# following needs magic library installed
#init_lattice=init_cml.magicSquare(sidelen)


init_lattice=init_cml.random_cml(sidelen,sidelen)
# setup pinned sites for test of pinning
pinmap=init_cml.random_binary(sidelen,sidelen,sparsity=0.2)

# predefined kernels:  'symm4','symm8','asymm', 'magic11','pinwheel'
cml = DiffusiveCML(init_lattice,kern='pinwheel',a=1.73,gl=.07,gg=0.0,use_pinned=False,pinned_mask=pinmap,use_activation=False)
stats = AnalysisCML(init_lattice, bin_spec=bin_count)

# setup activation array for demo in case it is selected
cml.activation=init_cml.center_ping_binary(sidelen,sidelen)
#cml.activation=initCML.randbin(sidelen,sidelen,sparsity=0.01)

# show initial state and create image
# translate float in interval -1,1 to colormap
if zoom_style == 'image':
    llshow = (cml.matrix + 1) * 128
else:
    llshow = zoom(((cml.matrix) + 1) * 128, scale, order=2)
im=plt.imshow(llshow)
# gist_heat, flag, gist_ncar,
plt.set_cmap('spring')
cml_ax=plt.axis([0, cml.matrix.shape[1], cml.matrix.shape[0], .001])
# create control sliders for alpha, local and global coupling
axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
# axes no longer takes axisbg?
axa_spread = plt.axes([0.2, 0.20, 0.65, 0.03])

axalpha = plt.axes([0.2, 0.25, 0.65, 0.03])
axlocal = plt.axes([0.2, 0.30, 0.65, 0.03])
axglob = plt.axes([0.2, 0.35, 0.65, 0.03])
axhist = plt.axes([0.2, 0.04, 0.65, 0.1])
axhist.set_ylim(0.0,0.4)
axhist.set_xlim(-1.0,1.0)
# creat parameter control sliders
a_slider = Slider(axalpha, 'alpha', 0.1, 2.0, valinit=cml.a)
gl_slider = Slider(axlocal, 'local coupling', 0.01, 0.5, valinit=cml.gl)
gg_slider = Slider(axglob, 'global coupling', 0.01, 0.2, valinit=cml.gg)
a_spread_slider = Slider(axa_spread, 'alpha spread', 0.0,2.0, valinit=cml.a_spread)

# display selection buttions
display_ax = plt.axes([0.025, 0.5, 0.20, 0.15])
display_ctl = RadioButtons(display_ax, ('state', 'spin', 'pinned', 'activation'), active=0)
display_ax.set_title('Lattice Display', size=12)
# evolution control selection buttons
actions_ax = plt.axes([0.025, 0.8, 0.20, 0.15])
actions_ctl = RadioButtons(actions_ax, ('free', 'pinned', 'spin_control', 'LFO', 'activation'), active=0)
actions_ax.set_title('Evolution Control', size=12)

# create instance of view and controller objects to encapsulate various rendering parameters and interpret GUI clicks
view=ViewCML(cml,stats,im=im,axhist=axhist,scale=scale,drawmod=1,drawtype='state')
display_ctl.on_clicked(view.display_func)
actions_ctl.on_clicked(view.actions_func)
# callback for parameter sliders
a_slider.on_changed(view.update_parms)
a_spread_slider.on_changed(view.update_parms)
gl_slider.on_changed(view.update_parms)
gg_slider.on_changed(view.update_parms)


# Create a new timer object. Set the interval to number of milliseconds
# (1000 is default) and tell the timer what function should be called.
timer = fig.canvas.new_timer(interval=1)
timer.add_callback(view.update)
timer.start()
plt.show()

