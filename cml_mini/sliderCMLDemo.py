
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import time
from scipy.ndimage import zoom
from diffusiveCML import DiffusiveCML
import initCML
from analysisCML import AnalysisCML

class ViewCML:

    def __init__(self,cml,stats,im=None,histax=None,scale=4,drawmod=1,drawtype='state'):
       self.cml=cml
       self.stats=stats
       self.im=im
       self.histax=histax
       self.scale=scale
       # drawmod draws only on modulo drawmod lattice iterations, to limit framerate for small lattices,
       # or to avoid visualization flicker in a demo
       # when you're in a state of high synchronization
       self.drawmod=drawmod
       # what variable is drawn in the image (state, spin)
       self.drawtype=drawtype
       # if spinControl True alpha is updated based on spin frustration
       self.spinControl=False
       self.lfo=False
       self.a_stack=[cml.a]
       self.gl_stack=[cml.gl]
       self.gg_stack=[cml.gg]
       self.cycle=5
       self.last_render_time=0

    def display_func(self,label):
        # set drawtype to match selected radio button
        view.drawtype = label
        view.im.figure.canvas.draw_idle()

    def actions_func(self, label):
        # control evolution options
        # set value true,
        if label == 'free':
            # default, set others False if selected again
            self.cml.usePinned=False
            self.spinControl=False
            self.lfo=False
            self.useActivation=False
        elif label =='pinned':
            self.cml.usePinned=True
        elif label =='spinControl':
            self.spinControl=True

        elif label == 'LFO':
             self.lfo=True
        elif label == 'activation':
              self.cml.useActivation=True
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
        self.cml.alphaUpdate()
        fig.canvas.draw_idle()

    def update(self):
        # perform diffusion and reaction steps
        view.cml.iterate()
        # print inter to console (todo- label in GUI)
        if view.cml.iter % 50 == 0: print view.cml.iter
        # compute measurements on the whole lattice
        view.stats.update(cml.matrix, cml.iter)
        # control example by spins (ref.  J.C. Perez, The New Way of Artificial Intelligence
        # experiment with spin control - number of spin transitions > threshold, or else decrease alpha
        # if a is chaotic, it will search and find a more stable (but probably still chaotic) value reducing spin transitions
        # This can be considered a form of unsupervised learning; essentially a similarity preserving hash into the partition space

        if view.spinControl:
            print "spinTrend ",view.stats.spinTrend
            if view.stats.spinTrend>100:

                view.cml.a=cml.a-.05
                a_slider.set_val(view.cml.a)
                print "reducing alpha", view.cml.a
                fig.canvas.draw_idle()
            # stop diffusion
            if cml.iter==1000:
                cml.kernType='zero'
                cml.kernelUpdate()
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
                view.cml.alphaUpdate()

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
            if cmlInit == 'image':
                llshow = (view.cml.matrix + 1) * 128
            else:
                # print 'scaling'
                if view.drawtype == 'state': llshow = zoom(((view.cml.matrix) + 1) * 128, view.scale, order=2)
                if view.drawtype == 'spin': llshow = zoom(((view.stats.spin) + 1) * 128, view.scale, order=2)
                if view.drawtype == 'pinned': llshow = zoom(((view.cml.pinnedMask)) * 128, view.scale, order=2)
                if view.drawtype == 'activation': llshow = zoom(((view.cml.activation)) * 128, view.scale, order=2)
            ## Display the data


            # zero out background in image to prevent flashing
            #llshow[np.where(llshow == llshow[0, 0])] = 0
            view.im.set_data(llshow)

            width = 2.0 / len(view.stats.edges)
            # set barplot centers; next line only works with uniform bins
            center = (view.stats.edges[:-1] + view.stats.edges[1:]) / 2
            view.histax.cla()
            view.histax.autoscale(enable=False)
            view.histax.bar(center, stats.bins, width=width, align='center')
            view.im.figure.canvas.draw()
            view.histax.figure.canvas.draw()


# end view constructor and event handlers; create GUI
# create window frame
fig, ax = plt.subplots()
plt.suptitle('State of each lattice site')
fig.canvas.set_window_title("Hybrid Diffusive-Global Coupled Map Lattice")
ax.patch.set_facecolor('black')
plt.subplots_adjust(left=0.25, bottom=0.45)

sidelen=100

# establish the number of bins (partitions in the state space) for the histogram display
binCount=64
# controls scaling of graphics i.e how many pixels per site; rendering order on image draw controls blurring
scale=1
# Various initial lattice styles.  var cmlInit controls zooming (pixel replication to make a bigger display
#cmlInit='image'
cmlInit=''
#initLattice=initCML.imageCML('./0683_01.png')
# primes will be slow for a big lattice
#initLattice=primesSquare(sidelen)
#initLattice=np.rot90(initLattice,2)
#initLattice=randomPing(sidelen,sidelen,scaleFactor=0.0)
#initLattice=randomCML(sidelen,sidelen)
# following needs magic library installed
#initLattice=magicSquare(sidelen)


initLattice=initCML.randomCML(sidelen,sidelen)
initPinnedMask=initCML.randbin(sidelen,sidelen,sparsity=0.2)
# predefined kernels:  'symm4','symm8','asymm', 'magic11'

cml = DiffusiveCML(initLattice,kern='pinwheel',a=1.73,gl=.07,gg=0.0,usePinned=False,pinnedMask=initPinnedMask,useActivation=True)
# setup activation array for demo. scaleFactor zero gives 1.0 in array center, 0 elsewhere because random values in range are multipled by 0
cml.activation=initCML.centerPingBin(sidelen,sidelen)
#cml.activation=initCML.randbin(sidelen,sidelen,sparsity=0.01)


stats = AnalysisCML(initLattice,binSpec=binCount)

# show initial state and create image
# translate float in interval -1,1 to colormap
if cmlInit == 'image':

    llshow = (cml.matrix + 1) * 128
else:
    llshow = zoom(((cml.matrix) + 1) * 128, scale, order=0)
im=plt.imshow(llshow)
# gist_heat, flag, gist_ncar,
plt.set_cmap('gist_ncar')
cml_ax=plt.axis([0, sidelen, sidelen, .001])

# create control sliders for alpha, local and global coupling
axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axa_spread = plt.axes([0.2, 0.20, 0.65, 0.03],axisbg=axcolor)
axalpha = plt.axes([0.2, 0.25, 0.65, 0.03],axisbg=axcolor)
axlocal = plt.axes([0.2, 0.30, 0.65, 0.03],axisbg=axcolor)
axglob = plt.axes([0.2, 0.35, 0.65, 0.03],axisbg=axcolor)
histax = plt.axes([0.2, 0.04, 0.65, 0.1])
histax.set_ylim(0.0,0.4)
histax.set_xlim(-1.0,1.0)
# creat parameter control sliders
a_slider = Slider(axalpha, 'alpha', 0.1, 2.0, valinit=cml.a)
gl_slider = Slider(axlocal, 'local coupling', 0.01, 0.5, valinit=cml.gl)
gg_slider = Slider(axglob, 'global coupling', 0.01, 0.2, valinit=cml.gg)
a_spread_slider = Slider(axa_spread, 'alpha spread', 0.0,2.0, valinit=cml.a_spread)

# display selection buttions
displayAx = plt.axes([0.025, 0.5, 0.15, 0.15])
displayCtl = RadioButtons(displayAx, ('state', 'spin','pinned','activation'), active=0)
displayAx.set_title('Lattice Display',size=12)
# evolution control selection buttons
actionsAx = plt.axes([0.025, 0.8, 0.15, 0.15])
actionsCtl = RadioButtons(actionsAx, ('free','pinned','spinControl','LFO','activation'), active=0)
actionsAx.set_title('Evolution Control',size=12)


# create instance of view and controller objects to encapsulate various rendering parameters and interpret GUI clicks
view=ViewCML(cml,stats,im=im,histax=histax,scale=scale,drawmod=1,drawtype='state')
displayCtl.on_clicked(view.display_func)
actionsCtl.on_clicked(view.actions_func)
# callback for parameter sliders
a_slider.on_changed(view.update_parms)
a_spread_slider.on_changed(view.update_parms)
gl_slider.on_changed(view.update_parms)
gg_slider.on_changed(view.update_parms)

# Create a new timer object. Set the interval to number of milliseconds
# (1000 is default) and tell the timer what function should be called.
timer = fig.canvas.new_timer(interval=1)
timer.add_callback(view.update)
# pass controller instance to avoid global vars in update_cml
#timer.add_callback(update_cml,cml_controller)
timer.start()
plt.show()

