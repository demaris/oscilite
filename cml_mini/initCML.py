import numpy as np
# These functions create rectangular values useful as initial conditions for research

from PIL import Image

def imageCML(imagePath,scaleFactor=1.0):
    img=Image.open(imagePath)
    img=img.convert('L')
    ll=np.array(img.getdata(),float).reshape(img.size[1], img.size[0])
    # scale into range 0,1
    ll=ll/np.max(ll) * scaleFactor
    return ll

def randomCML(xside,yside,cmlType='KK',scaleFactor=1.0):
    ll=np.random.rand(xside,yside)
    if cmlType == 'KK':
        ll=((ll*1.99)-.999)*scaleFactor
    else:
        ll=ll*scaleFactor
    return ll

# create random matrix of positive values between zero and one
def randomCMLpos(xside,yside,scaleFactor=1.0):
    ll=np.random.rand(xside,yside)*scaleFactor

    return ll

# random field centered around zero with a single point in center
# produces a mandala-like structure with symmetric kernels
def randomPing(xside,yside,cmlType='KK',scaleFactor=.000000000001):
    ll=np.random.rand(xside,yside)*scaleFactor
    # KK (Kaneko) type logistic map has values -1:1, rather than more common logistic 0:1
    if cmlType == 'KK':
        ll=((ll*1.999)-.999*scaleFactor)
        ll[xside/2,yside/2]=1.0
    else:
        # in case we add a domain 0 to 1 map
        ll[xside/2,yside/2]=.99
    return ll

def centerPingBin(xside,yside):
    ll=np.zeros((xside,yside))
    ll[xside/2,yside/2]=1
    print ll[xside/2-1:xside/2+1,yside/2-1:yside/2+1]
    return ll

# create a sparse pattern of ones in a zero background
def randbin(xside,yside,sparsity=0.5,scaleFactor=1.0):
    ll=np.random.rand(xside, yside)
    ll[np.where(ll>=1.0-sparsity)]=1.0
    ll[np.where(ll<1.0-sparsity)]=0.0
    return ll

# The rest fun but possibly irrelevant for research.
# However interesting things happen

# requires magic_square package: pip magic_square
#from magic_square import
"""
def magicSquare(n):
    ll=magic(n)/(n*n*1.0)
    return ll
"""


# this is for lulz; it's pretty slow in pure python
def primesSquare(n):
    N=n*n
    primes  = []
    chkthis = 2
    while len(primes) < N:
        ptest    = [chkthis for i in primes if chkthis%i == 0]
        primes  += [] if ptest else [chkthis]
        chkthis += 1

    ll=np.reshape(primes,(n,n))/(primes[N-1]*1.0) # mult by one to get floats otherwise you get all zero
    return ll
