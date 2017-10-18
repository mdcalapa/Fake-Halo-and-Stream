#Making a fake data set with velocity dispersion at each point and a power law density fall off of points to test emcee for differentiating the stream from the smooth halo

import numpy as np
import math
import matplotlib.pyplot as plt
import pdb

#fake data set of stars - pick enough so there's at least 100 points between radius 100 and 110
#give each a random MC'd radius
#then another array of equal size of gauss velocities

for w in range(2,7):
        x=np.random.rand(10**w)*150.0
        x=x[x>=1.0]
        #prob3d=x**(-3.0)*np.log10(150.0)*3.0
        prob3d=x**(-1.0)*np.log10(150.0)*1.0
        y=np.random.rand(len(x))*np.max(prob3d)*1.1
        keep=np.where(y < prob3d)
        discard = np.where(y >= prob3d)
        r_t = x[keep]
        if len(r_t[np.where((r_t>=100) & (r_t<=110))])>= 100:
                break

z=np.random.uniform(low=0,high=1,size=len(r_t)) #random number for cos(gb)

#r_t=x[keep]
g_t=z

g_t=np.array(g_t)
r_t=np.array(r_t)
v1=r_t[(7.0<=r_t)]
meanrho=v1[(9.0>=v1)]

#vel time
mu, sigma = 0.0, 100.0 # mean and standard deviation
q=len(r_t)
v = np.random.normal(mu, sigma, q)

#see if stars right density

#hist,bin_edges=np.histogram(x[keep],bins=100,range=(0,150))
#dens=np.zeros(len(hist))
#densr=np.zeros(len(hist))

#for i in range(0,len(hist)):
#        dens[i]=hist[i]/((4./3.)*math.pi*(bin_edges[i+1]**3.0-bin_edges[i]**3.0))
#        densr[i]=(bin_edges[i+1]+bin_edges[i])/2.0


rmin=np.min(r_t)
rmax=np.max(r_t)

alph = np.arange(-3.1,8.9,0.4)
sig= np.arange(20.0,920.0,30.0)

def lbit(meanrho,gb,v,r,sig,alph,rmin,rmax):
    xdos=np.zeros(len(r))
    for l in range(0,len(r)):
        const=(1./(3-alph))*(rmax**(3-alph) - rmin**(3-alph))
        rho=(r[l]**(2.0-alph))/const
        vp=np.exp(-(((v[l])**2.0)/(2.0*(sig**2.0))))
        xdos[l]=vp*rho/(np.sqrt(2*math.pi*sig**2.0))
    p = np.sum(np.log10(xdos))
    #pdb.set_trace()
    return p

L=np.zeros((len(sig),len(alph)))
for k in range(len(sig)):
    print(k)
    for j in range(0,len(alph)):
        L[k,j]=lbit(meanrho,g_t,v,np.sort(x[keep]),sig[k],alph[j],rmin,rmax)
#pdb.set_trace()

#plt.plot(sig,L[:,0])

#COUNTOR ROUTINE
CS=plt.contour(sig, alph, L*(1e-5), colors='black')
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel('Sigma')
plt.ylabel('Alpha')
plt.show()

plt.imshow(L, interpolation='none', origin='lower',extent=[min(alph),max(alph),min(sig),max(sig)],aspect='auto')
plt.show()
