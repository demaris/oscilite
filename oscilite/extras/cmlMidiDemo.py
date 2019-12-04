from diffusive_cml import DiffusiveCML
from analysis_cml import AnalysisCML
import init_cml
import mido
import signal
import threading
import time
import numpy as np
import math
# defines some scales, frome lawriecapb.co.uk/theblog/indexphdp/archives/881
Blues =[41,42,43,46,48,51,53,54,55,58,60,63,65,66,67,70,72,75,77,78,79,82,84,87,89,90,91,94,96,99,101,102]
BluesDiminished =[48,49,51,52,54,55,56,58,60,61,63,64,66,67,68,70,72,73,75,76,78,79,80,82,84,85,87,88,90,91,92,94]
Dorian =[25,27,30,32,34,37,39,42,44,46,49,51,54,56,58,61,63,66,68,70,73,75,78,80,82,85,87,90,92,94,97,99]
FullMino =[51,53,55,56,57,58,59,60,62,63,65,67,68,69,70,71,72,74,75,77,79,80,81,82,83,84,86,87,89,91,92,93]
HarmonicMajor =[44,47,48,50,52,53,55,56,59,60,62,64,65,67,68,71,72,74,76,77,79,80,83,84,86,88,89,91,92,95,96,98]
class InterpretLattice:
    def __init__(self,sequence,fun,args,kwargs): # take out references here since they came in via sequence set
        self.sequence=sequence
        self.fun=fun
        self.args=args
        self.kwargs=kwargs
        self.pause=0.0
        self.iter=0
        print self.kwargs
        print "warmup ", self.sequence.warmup
        self.iter=0
        while self.iter<warmup:
            self.sequence.lattice.iterate()
            self.sequence.analysis.update(self.sequence.lattice.matrix,self.iter)
            self.iter = self.iter + 1

    def advance(self):
        # adjust params based on schedule
        # ignoring schedule for now, do random walk in param space from initial values
        #self.sequence.lattice.a=self.sequence.lattice.a+
        #check bounds and set to extrema if exceed min or max
        #self.sequence.lattice.gl=self.sequence.lattice.gl+rand

        # update convolution kernel after resetting gl

        # update lattice and stats
        self.sequence.lattice.iterate()
        self.sequence.analysis.update(self.sequence.lattice.matrix,self.iter)
        self.iter = self.iter + 1
        self.pause=self.fun(self.args,self.kwargs)
        time.sleep(self.pause)
        self.advance()

# assume multiple sequence generators
class Sequencer:
    def __init__(self,lattice,analysis,out_port,schedule,warmup=10,rate=1.0,loop=False,loop_count=0,interpreter_fun=''):
        self.lattice= lattice
        self.analysis= analysis
        self.out = out_port
        self.schedule=schedule #list of one or more lists
        self.warmup=warmup
        self.loop=loop #True or False to execute schedule for loop_count
        self.loop_count=loop_count # negative value means run forever until killed
        self.fun=None
        self.args=[]
        self.kwargs=[]
        self.rate=rate
        self.note_queue=[]

    def set_interpreter(self,fun,*args,**kwargs):
        self.fun=fun
        self.args=args
        self.kwargs=kwargs

    def start(self):
        if self.fun is None:
            print "Sequence aborted, no interpreter function defined. Use set_interpreter(fun,args,kwargs)"

        else:
            interpreter = InterpretLattice(self, self.fun, self.args, self.kwargs)
            t = threading.Timer(self.rate, interpreter.advance())
            t.start()


# example of an interpreter function which can access lattice and stats via parent arg
# interpreters map stats to note, velocity, duration, pause
# interpeter functions start notes, post notes (if more than one) to parent, kill old timer and create timers
# for note off and advance
# interpreter functions should return a pause (rest, breathing value) before next advance
def melody(parent,kwdict):
    queue=[]
    seq=parent[0]
    dur=kwdict['dur']
    threshold=kwdict['threshold']
    tension_max=1
    #print seq.analysis
    # choose note as highest valued element, vol from pop, duration from delta to second highest
    sorted=np.argsort(seq.analysis.bins)
    print seq.analysis.bins
    print sorted
    # run through 4 highest values
    for i in range(4):
        n=i+1 # we're going to index from end of sorted with -1 but range would give -0
        seq.out.send(mido.Message('note_on',note=Blues[sorted[-n]%32],velocity=100))
        # make a queue in case we do more than one note
        queue.append(Blues[sorted[-n]%32])
        # compute duration from some stats
        # look at delta between highest value and next highest
        tension=seq.analysis.bins[sorted[-n]]/seq.analysis.bins[sorted[-(n+1)]]
        # let's not get crazy with tension, limit to 4
        tension=min([tension,4])
        # snap to nearest integer or not
        #tension=math.floor(tension+0.5)
        tension_max=max([tension_max,tension])
        tension=max([tension,1])  # don't want tension to be zero for sleep timer below
        print "tension ",tension
        # create timer based on duration callback should kill previous timer, create new timer to make note off and advance interpeter
        time.sleep(dur * tension)
        # turn off note turned on
        seq.out.send(mido.Message('note_off',note=queue.pop()))
    # return value is used to pause with sleep() after a note
    return 1.0 * tension_max

if __name__ == '__main__':
    with mido.open_output('IAC Driver IAC Bus 1') as out:
        # create CML and analysis, schedule lists
        initLattice = init_cml.randomCML(100, 100)
        # predefined kernels:  'symm4','symm8','asymm', 'magic11','pinwheel'
        cml = DiffusiveCML(initLattice, kern='pinwheel', a=1.7, gl=.2, gg=0.08, usePinned=False,useActivation=False)
        stats = AnalysisCML(initLattice, bin_spec=64)
        # create a schedule
        # sequence of lists of form a,gl,gg,steps
        # a is alpha, nonlinearity parameter of maps
        # gl is local coupling constant
        # gg is global mean field coupling constant
        warmup=20
        schedule=[[1.73,0.2,0.05,5],[1.5,0.05,0.05,5]]
        # create sequencer for melody
        sequencer=Sequencer(cml,stats,out,schedule,warmup,rate=3.0,interpreter_fun=melody)
        # all interpreters take a parent sequence as first parameter, then any parameters specific to the algorithm
        sequencer.set_interpreter(melody,sequencer,threshold=0.2,dur=.125)
        sequencer.start()
        # maybe create another sequencer for counter melody or chords

