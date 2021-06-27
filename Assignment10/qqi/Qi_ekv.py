#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.special import lambertw as lw


# In[2]:


tox=1e-7
tsi=500e-7
q=1.6e-19
ϵ0=8.84e-14
ϵs= 11.4*ϵ0
ϵox= 3.9*ϵ0
Cox=ϵox/tox
Na=1e+16
ni=1e+10
n0=(ni**2/Na)
μ=1000 
r=1
W=50e-7
L=r*W
ϕt=0.026
ϕf=ϕt*np.log(Na/ni)
ϕms=4.26 - 4.05
Qox=0
β=μ*Cox*(W/L)

Vfb = ϕms - Qox/Cox

γ=np.sqrt(2*q*ϵs*Na)/(Cox)


# In[3]:


def ψp(vg):
    
    y = vg - Vfb - (γ**2)*( np.sqrt( (vg - Vfb)/(γ**2) + (1/4) ) - (1/2) )
    
    n = 1 + γ/(2*np.sqrt(y))
    
    q=-2*n*ϕt*Cox
    
    y1=y/ϕt 
    ϕf1=ϕf/ϕt
    γ1=γ/np.sqrt(ϕt)
    
    
    
    vp = y1 - 2*ϕf1 - np.log( (4*n/γ1)*np.sqrt(y1))
    
    return  q , vp, n


# In[4]:


def Qi(Vcb=0,i=ϕf,f=4,m=50 ):
    vcb = Vcb
    Vg = np.linspace(i,f,m)

    VVP=[]

    Qekvi=[]
    
    SOL=[]

    for vg in Vg:
    
        w = 2*np.exp( ψp(vg)[1] - vcb )
    
        sol = lw(w)
    
        vv_p = ψp(vg)[1] - vcb
    
        VVP.append(vv_p)
    
    
        sol = abs(sol)
    
        qi= (sol/2)*( ψp(vg)[0] )
        
        SOL.append(sol/2)
    
        Qekvi.append(qi)
    
    Qekvi = np.asarray(Qekvi)
    n = ψp(vg)[2]
    SOL=np.asarray(SOL)
    
    return SOL, Vg, n


# In[ ]:




