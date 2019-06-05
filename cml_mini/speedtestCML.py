from diffusiveCML import DiffusiveCML
import initCML
import time
# just initialize CML and run at full speed to measure performance, without pinning, activation, or analysis as in the demo
initLattice=initCML.randomCML(100,100)
cml = DiffusiveCML(initLattice,kern='pinwheel',a=1.73,gl=.07,gg=0.0,usePinned=False,useActivation=False)
cycles=1000
i=0
print "Starting iterations"
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print "Speed was ",cps," lattice updates per second, large pinwheel kernel"

cml = DiffusiveCML(initLattice,kern='symm4',a=1.73,gl=.07,gg=0.0,usePinned=False,useActivation=False)

i=0
print "Starting iterations"
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print "Speed was ",cps," lattice updates per second, symm4 kernel, no global coupling"


cml = DiffusiveCML(initLattice,kern='symm4',a=1.73,gl=.07,gg=0.1,usePinned=False,useActivation=False)
i=0
print "Starting iterations"
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print "Speed was ",cps," lattice updates per second, symm4 kernel, with global coupling"


