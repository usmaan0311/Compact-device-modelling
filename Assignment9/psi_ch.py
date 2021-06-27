#!/usr/bin/env python
# coding: utf-8

# In[1]:


#%matplotlib widget
import scipy.integrate as sci, matplotlib.pyplot as plt, numpy as np
from MOSFET import *
from PNR import *
from NR import *
import time


# # charges:
# 
# $Q_{B}^{'} = -\sqrt{2q\epsilon_{s}N_{A}}\sqrt{\psi_{s}(x)}$
# 
# $Q_{B} = W\int_{0}^{L}Q_{B}^{'}\,dx$
# 
# $Q_{G}^{'} = C_{ox}^{'}( V_{GB} - V_{FB} - \psi_{s}(x) ) - Q_{o}^{'}$
# 
# $Q_{G} = W\int_{0}^{L}Q_{G}^{'}\,dx$
# 
# $Q_{I}^{'} = -\sqrt{2q\epsilon_{s}N_{A}}[ \sqrt{\psi_{s}(x) + \phi_{t}e^{\frac{\psi_{s}(x)\ -\ (2\phi_{f} + V_{cb})}{\phi_{t}}}    } - \sqrt{\psi_{s}(x)}]$
# 
# - From $Q_{I}^{'}$,  $Q_{S}$ and $Q_{D}$ can be calculated as follows:
# 
# $Q_{D} = W\int_{0}^{L}\frac{x}{L}Q_{I}^{'}\,dx$
# 
# $Q_{S} = W\int_{0}^{L}(1 - \frac{x}{L})Q_{I}^{'}\,dx$

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
μ=1000 
r=1
W=50e-7
L=r*W
ϕt=0.026
ϕf=ϕt*np.log(Na/ni)
ϕms=4.26 - 4.05
Qox=0

Vfb = ϕms - Qox/Cox

γ=np.sqrt(2*q*ϵs*Na)/(Cox)

def fun(Ψs,Vgb):
    
    f=μ*Cox*(( Vgb - Vfb + ϕt  )*Ψs - 0.5*(Ψs**2) - (2/3)*γ*(Ψs**(3/2)) + ϕt*γ*(Ψs**(1/2)))
    
    return f

L=500e-7
vg=4
Vg=np.linspace(1,2,6)
Vds=np.linspace(0,5,6)


def Ψch(Ecb,Vg):
    
    for vg in Vg:

        x0,y0=poissonMOS(Vgs1=vg,ΔEcb=0)
        xL,yL=poissonMOS(Vgs1=vg,ΔEcb=Ecb)


        ψ0=y0[0][0]
        f0=fun(ψ0,vg)

        ψL=yL[0][0]
        fL=fun(ψL,vg)

    

        X=np.linspace(0,L,100)
    
        ψx=[]
        
    
        
   
        for x in X:
            g = lambda ψ : (x/L) - ( fun(ψ,vg) - f0 )/( fL - f0 )
            k=Newton(g,f0)
            ψx.append(k[0])
        return ψx, X


# In[ ]:




