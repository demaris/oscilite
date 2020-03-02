from diffusive_cml import DiffusiveCML
import init_cml
import time
# just initialize CML and run at full speed to measure performance, without pinning, activation, or analysis as in the demo
init_lattice=init_cml.random_cml(100, 100)
cml = DiffusiveCML(init_lattice,kern='pinwheel',a=1.73,gl=.07,gg=0.0,use_pinned=False,use_activation=False)
cycles=1000
i=0
print('Starting iterations')
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print('Speed was ',cps,' lattice updates per second, large pinwheel kernel')

cml = DiffusiveCML(init_lattice,kern='symm4',a=1.73,gl=.07,gg=0.0,use_pinned=False,use_activation=False)

i=0
print('Starting iterations')
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print('Speed was ',cps,' lattice updates per second, symm4 kernel, no global coupling')


cml = DiffusiveCML(init_lattice,kern='symm4',a=1.73,gl=.07,gg=0.1,use_pinned=False,use_activation=False)
i=0
print('Starting iterations')
start=time.clock()
while i<cycles:
    cml.iterate()
    i+=1
end=time.clock()
secs=end-start
cps=cycles/secs
print('Speed was ',cps,' lattice updates per second, symm4 kernel, with global coupling')


