
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
import matplotlib.rcsetup as rcsetup
from matplotlib.widgets import Slider, Button, RadioButtons
import time
from scipy.ndimage import zoom
from diffusive_cml import DiffusiveCML
import init_cml
from analysis_cml import AnalysisCML
import mido
from numpy import where, size
class ViewCML:

    def __init__(self,cml,stats,im=None,axhist=None,scale=4,drawmod=1,drawtype='state'):
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
       self.send_allowed=True

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

    def lock_manager(self):

        print "setting allow video"
        self.send_allowed=True
        # sprinkle some noise

        #lock_timer.stop()

    def update(self):
        # check midi signals for parameter changes

        for msg in inport.iter_pending():
            print msg
            # assuming we get only cc changes
            tmp = "%s" % (msg)
            tmp = tmp.split(' ')
            # control number or note third
            cc = tmp[2].split('=')[1]
            print cc
            val = tmp[3].split('=')[1]
            last_midi_event[cc] = val
            # if cc > 16 it must be the note 48 reset button
            if int(cc) > 16:
                print 'reset cml to random'
                self.cml.matrix=init_cml.random_cml(sidelenx,sideleny)


        # sync values in gui at some point
        # scale sliders to allowable range and set parms
        # chaos level a
        a_min=0.0; a_max=2.0;
        a_range=a_max-a_min;
        self.cml.a = a_min + int(last_midi_event['11'])/128.0 * a_range
        #print 'chaos level', self.cml.a
        #a_slider.set_val(view.cml.a)
        gl_min = 0;
        gl_max = 0.2;
        gl_range = gl_max - gl_min;
        self.cml.gl = gl_min + int(last_midi_event['12']) / 128.0 * gl_range
        #print 'local coupling', self.cml.gl
        #gl_slider.set_val(self.cml.gl)
        gg_min = 0;
        gg_max = 0.1;
        gg_range = gg_max - gg_min;
        self.cml.gg = gg_min + int(last_midi_event['13']) / 128.0 * gg_range
        #print 'global coupling', self.cml.gg
        #gg_slider.set_val(self.cml.gg)
        a_range_min = 0;
        a_range_max = 2.0;
        a_range_range = a_range_max - a_range_min;
        self.cml.a_spread = a_range_min + (int(last_midi_event['14']) / 128.0) * a_range_range
        #print 'chaos gradient', self.cml.a_spread
        #a_spread_slider.set_val(self.cml.a_spread)
        self.cml.alpha_update()
        #fig.canvas.draw_idle()
        # perform diffusion and reaction steps

        view.cml.iterate()
        # print inter to console
        if view.cml.iter % 50 == 0: print view.cml.iter
        # compute measurements on the whole lattice
        view.stats.update(self.cml.matrix)
        # control example by spins (ref.  J.C. Perez, The New Way of Artificial Intelligence
        # experiment with spin control - number of spin transitions > threshold, or else decrease alpha
        # if a is chaotic, it will search and find a more stable (but probably still chaotic) value reducing spin transitions
        # This can be considered a form of unsupervised learning; essentially a similarity preserving hash into the partition space
        #print 'x ',size(view.cml.matrix, 0),'y ',size(view.cml.matrix, 1)
        if view.spin_control:
            print "spin_trend ",view.stats.spinTrend
            if view.stats.spinTrend>100:

                view.cml.a=cml.a-.05
                a_slider.set_val(view.cml.a)
                print "reducing alpha", view.cml.a
                fig.canvas.draw_idle()
            # stop diffusion
            if cml.iter==1000:
                cml.kern_type= 'zero'
                cml.kernel_update()
                #cml.a=1.7
                print "stopping diffusion and quenching a= ",cml.a

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
            if view.cml.iter == 10:
                print(fps, 'fps')
            # if an image is big, don't do the scaling with scipy zoom but rather use it direct.
            if zoom_style == 'image':
                llshow = (view.cml.matrix + 1) * 128
            else:
                # print 'scaling'
                tmp=True
                #llshow = zoom(((view.cml.matrix) + 1) * 128, view.scale, order=2)
                #if view.draw_type == 'spin': llshow = zoom(((view.stats.spin) + 1) * 128, view.scale, order=1)
                #if view.draw_type == 'pinned': llshow = zoom(((view.cml.pinnedMask)) * 128, view.scale, order=1)
                #if view.draw_type == 'activation': llshow = zoom(((view.cml.activation)) * 128, view.scale, order=1)
            ## Display the data

            view.im.set_data(zoom(((view.cml.matrix) + 1) * 128, view.scale, order=2))
            width = 2.0 / len(view.stats.edges)
            # set barplot centers; next line only works with uniform bins
            center = (view.stats.edges[:-1] + view.stats.edges[1:]) / 2
            view.axhist.cla()
            view.axhist.autoscale(enable=False)
            view.axhist.bar(center, stats.bins, width=width, align='center')
            #view.im.figure.canvas.draw()
            ax.draw_artist(im)
            view.im.figure.canvas.update()
            view.im.figure.canvas.flush_events()
            #view.axhist.figure.canvas.draw()
            # check for stats.bins exceeding threshold and send midi event if send_allowed
            # choose first (lowest bin over threshold if more than one
            # skip last two bins which tend to grab a lot of state
            # in fact just grab from center of 64 bins
            max_bin_val = max(stats.bins[15:47])
            #print 'max bin val', max_bin_val

            if max_bin_val > midi_threshold:
                max_bin = where(stats.bins == max_bin_val)[0].tolist()
                #print 'max bin index', max_bin[0]-15
                if view.send_allowed:
                # construct message
                    #print 'sending note number ',max_bin[0]-15
                    # shift note down by bin offset

                    outport.send(mido.Message('note_on',channel=0,note=max_bin[0]-15))
                    view.send_allowed = False
                    # restart timer
                    #
                    #print "new clips blocked"

                else:
                    temp=True

                    #print 'clip blocked'

# end view constructor and event handlers;
# setup midi ports


outport=mido.open_output('IAC Driver IAC Bus 1')
last_midi_event={'11':'64','12':'64','13':'0','14':'0'}
inport= mido.open_input('Network Session 1')
midi_threshold=0.03

# check backend avail
print 'current backend ',plt.get_backend()
print 'all ', rcsetup.all_backends
# create GUI


# create window frame
fig, ax = plt.subplots()
plt.suptitle('State of each lattice site')
fig.canvas.set_window_title("Hybrid Diffusive-Global Coupled Map Lattice")
ax.patch.set_facecolor('black')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
plt.subplots_adjust(left=0.25, bottom=0.45)

sidelenx=120
sideleny=120

# establish the number of bins (partitions in the state space) for the histogram display
bin_count=64
# controls scaling of graphics i.e how many pixels per site; rendering order on image draw controls blurring
scale=1
# var zoom_style controls zooming (pixel replication to make a bigger display), default is to zoom
zoom_style=''
# Various initial lattice styles.
#init_latticeattice=init_cml.image_cml('./MNISTdigits.png')
# uncomment if initializing with an image
#zoom_style='image'
# primes will be slow for a big lattice
#init_lattice=init_cml.primes_square(sidelen)
#init_lattice=np.rot90(initLattice,2)
#init_lattice=init_cml.random_ping(sidelen,sidelen,scaleFactor=0.0)
#init_lattice=init_cml.random_cml(sidelen,sidelen)
# following needs magic library installed
#init_lattice=init_cml.magicSquare(sidelen)


init_lattice=init_cml.random_cml(sideleny,sidelenx)
# setup pinned sites forx
pinmap=init_cml.random_binary(sideleny,sidelenx,sparsity=0.2)

# predefined kernels:  'symm4','symm8','asymm', 'magic11','pinwheel'
cml = DiffusiveCML(init_lattice,kern='symm4',a=1.73,gl=.07,gg=0.0,use_pinned=False,pinned_mask=pinmap,use_activation=False)
stats = AnalysisCML(init_lattice, bin_spec=bin_count,do_spins=1)

# setup activation array for demo in case it is selected
cml.activation=init_cml.center_ping_binary(sideleny,sidelenx)
#cml.activation=initCML.randbin(sidelen,sidelen,sparsity=0.01)

# show initial state and create image
# translate float in interval -1,1 to colormap
if zoom_style == 'image':
    llshow = (cml.matrix + 1) * 128
else:
    llshow = zoom(((cml.matrix) + 1) * 128, scale, order=0)
im=plt.imshow(llshow)
# gist_heat, flag, gist_ncar,
plt.set_cmap('gist_ncar')
cml_ax=plt.axis([0, sidelenx, sideleny, .001])
# create control sliders for alpha, local and global coupling
axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axa_spread = plt.axes([0.2, 0.20, 0.65, 0.03],axisbg=axcolor)
axalpha = plt.axes([0.2, 0.25, 0.65, 0.03],axisbg=axcolor)
axlocal = plt.axes([0.2, 0.30, 0.65, 0.03],axisbg=axcolor)
axglob = plt.axes([0.2, 0.35, 0.65, 0.03],axisbg=axcolor)
axhist = plt.axes([0.2, 0.04, 0.65, 0.1])
axhist.set_ylim(0.0,0.4)
axhist.set_xlim(-1.0,1.0)
# creat parameter control sliders
a_slider = Slider(axalpha, 'alpha', 0.1, 2.0, valinit=cml.a)
gl_slider = Slider(axlocal, 'local coupling', 0.01, 0.5, valinit=cml.gl)
gg_slider = Slider(axglob, 'global coupling', 0.01, 0.2, valinit=cml.gg)
a_spread_slider = Slider(axa_spread, 'alpha spread', 0.0,2.0, valinit=cml.a_spread)

# display selection buttions
display_ax = plt.axes([0.025, 0.5, 0.15, 0.15])
display_ctl = RadioButtons(display_ax, ('state', 'spin', 'pinned', 'activation'), active=0)
display_ax.set_title('Lattice Display', size=12)
# evolution control selection buttons
actions_ax = plt.axes([0.025, 0.8, 0.15, 0.15])
actions_ctl = RadioButtons(actions_ax, ('free', 'pinned', 'spin_control', 'LFO', 'activation'), active=0)
actions_ax.set_title('Evolution Control', size=12)

# create instance of view and controller objects to encapsulate various rendering parameters and interpret GUI clicks
view=ViewCML(cml,stats,im=im,axhist=axhist,scale=scale,drawmod=1,drawtype='state')
"""
display_ctl.on_clicked(view.display_func)
actions_ctl.on_clicked(view.actions_func)
# callback for parameter sliders
a_slider.on_changed(view.update_parms)
a_spread_slider.on_changed(view.update_parms)
gl_slider.on_changed(view.update_parms)
gg_slider.on_changed(view.update_parms)

"""
# Create a new timer object. Set the interval to number of milliseconds
# (1000 is default) and tell the timer what function should be called.
timer = fig.canvas.new_timer(interval=1)
timer.add_callback(view.update)
timer.start()
# create a lock disallowing too frequent notes triggering new overlay videos

lock_timer=fig.canvas.new_timer(interval=12000)
lock_timer.add_callback(view.lock_manager)
lock_timer.start()

plt.show()

