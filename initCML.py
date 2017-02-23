import numpy as np
from numpy.random import rand

# some initialization types need library support not included in this version; try pip magic_square
#from magic_square import
"""
def magicSquare(n):
    ll=magic(n)/(n*n*1.0)
    return ll
"""

#from PIL import Image
"""
def imageCML(imagePath,scaleFactor=1.0):
    img=Image.open(imagePath)
    #img=Image.open('./sri_mandala.jpg')
    img=img.convert('L')

    ll=np.array(img.getdata(),float).reshape(img.size[1], img.size[0])
    xlen=img.size[1]
    ylen=img.size[0]
    # scale into range 0,1
    ll=ll/np.max(ll) * scaleFactor
    return ll
"""
def randomCML(xside,yside,cmlType='KK',scaleFactor=1.0):
    ll=np.random.rand(xside,yside)
    if cmlType == 'KK':
        ll=((ll*1.99)-.999)*scaleFactor
    else:
        ll=ll*scaleFactor
    return ll

def randomCMLpos(xside,yside,cmlType='KK',scaleFactor=1.0):
    ll=np.random.rand(xside,yside)
    if cmlType == 'KK':
        ll=((ll*1.999)-.999)*scaleFactor
    else:
        ll=ll*scaleFactor
    return ll


def randomPing(xside,yside,cmlType='KK',scaleFactor=.000000000001):
    ll=np.random.rand(xside,yside)*scaleFactor
    if cmlType == 'KK':
        ll=((ll*1.999)-.999*scaleFactor)
        ll[xside/2,yside/2]=.99
    else:
        # in case we add a domain 0 to 1 map
        ll[xside/2,yside/2]=.99
    return ll

def randbin(xside,yside,scaleFactor=1.0):
    ll=np.random.rand(xside, yside)
    ll[np.where(ll>=.5)]=1.0
    ll[np.where(ll<.5)]=0.0
    return ll

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
