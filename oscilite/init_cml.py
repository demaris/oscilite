import numpy as np
# These functions create rectangular values useful as initial conditions for research

from PIL import Image

def image_cml(image_path,scale_factor=1.0):
    img=Image.open(image_path)
    img=img.convert('L')
    ll=np.array(img.getdata(),float).reshape(img.size[1], img.size[0])
    # scale into range 0,1
    ll=ll/np.max(ll) * scale_factor
    return ll

def random_cml(xside,yside,cmlType='KK',scale_factor=1.0):
    ll=np.random.rand(xside,yside)
    if cmlType == 'KK':
        ll=((ll*1.99)-.999)*scale_factor
    else:
        ll=ll*scaleFactor
    return ll

# create random matrix of positive values between zero and one
def random_cml_pos(xside,yside,scale_factor=1.0):
    ll=np.random.rand(xside,yside)*scale_factor

    return ll

# random field centered around zero with a single point in center
# produces a mandala-like structure with symmetric kernels
def random_ping(xside,yside,cmlType='KK',scale_factor=.000000000001):
    ll=np.random.rand(xside,yside)*scale_factor
    # KK (Kaneko) type logistic map has values -1:1, rather than more common logistic 0:1
    if cmlType == 'KK':
        ll=((ll*1.999)-.999*scale_factor)
        ll[xside/2,yside/2]=1.0
    else:
        # in case we add a domain 0 to 1 map
        ll[xside/2,yside/2]=.99
    return ll

def center_ping_binary(xside,yside):
    ll=np.zeros((xside,yside))
    ll[xside/2,yside/2]=1
    print ll[xside/2-1:xside/2+1,yside/2-1:yside/2+1]
    return ll

# create a sparse pattern of ones in a zero background
def random_binary(xside,yside,sparsity=0.5):
    ll=np.random.rand(xside, yside)
    ll[np.where(ll>=1.0-sparsity)]=1.0
    ll[np.where(ll<1.0-sparsity)]=0.0
    return ll

# The rest is fun but possibly irrelevant for research.
# However interesting things happen!

# requires magic_square package: pip magic_square
"""
from magic_square import magic

def magic_square(n):
    ll=magic(n)/(n*n*1.0)
    return ll
"""

# this is for lulz; it's pretty slow in pure python
# create an initial condition based on a field of prime values scaled by the max prime, mod side length
def primes_square(n):
    N=n*n
    primes  = []
    chkthis = 2
    while len(primes) < N:
        ptest    = [chkthis for i in primes if chkthis%i == 0]
        primes  += [] if ptest else [chkthis]
        chkthis += 1
    ll=np.reshape(primes,(n,n))/(primes[N-1]*1.0) # mult by one to get floats otherwise you get all zero
    return ll
