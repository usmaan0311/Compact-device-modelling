#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from NR import *                                                       


# In[2]:


def PoissonNR(v_in, v_f, Ec0=0, Ecf=5, n=1,m=100):
    q=1.6e-19
    ϵ0=8.84e-14
    ϵs = 11.4*ϵ0
    ϵox = 3.9*ϵ0
    Na=1e+16
    ni=1e+10
    tox=1e-7
    Cox = ϵox/tox
    ϕt = 0.026
    ϕms= 4.26 - 4.05
    Qox=0
    ϕf = ϕt*np.log(Na/ni)
    ΔECB=np.linspace(Ec0,Ecf,n)
    Vgb = np.linspace(v_in,v_f,m)
    Ψg = []
    for Ecb in ΔECB:
        Gate=[]
        ψ_y=[]
        dis=[]
        Ψ=[]
        
        
        γ = np.sqrt(2*q*ϵs*Na)/(Cox)
        VFB = ϕms - (Qox/Cox)

        
        iter = []
        G=[]
    
            
        for vgb in Vgb:
            g=0
    
            if (vgb < 0):
                g1 = ( (-γ/2) + np.sqrt( abs((-γ/2)**2 - (VFB - abs(vgb)) )) )**2
        
                g2 = -2*ϕt*np.log((abs(vgb - VFB))/(γ*np.sqrt(ϕt)))
        
                g= min(g1,g2)
        
            else:
                g1 = ( (-γ/2) + np.sqrt( abs((-γ/2)**2 - (VFB - abs(vgb)) )) )**2
        
                g3 = Ecb + 2*( ϕf + ϕt*np.log((abs(abs(vgb) - VFB))/(γ*np.sqrt(ϕt))) )
        
                g=min(g1,g3)
        
        
            
        
    
            
            print(f"guess is {g}\n")
    
            f = lambda ψs : VFB + ψs + γ*np.sqrt( abs(ψs) + ϕt*np.exp((ψs - 2*ϕf - Ecb)/ϕt) -ϕt + ϕt*np.exp(-ψs/ϕt) - np.exp(-2*ϕf/ϕt)*( ϕt*np.exp(-Ecb/ϕt) + abs(ψs) )) - abs(vgb)
            y=Newton(f,g)
            Ψ.append(y[0])
            iter.append(y[1])
            G.append(g)
        
        Ψg.append(Ψ)
        
        return Ψ, Vgb


# In[3]:


def PoissonNRD(v_in, v_f, Ec0=0, Ecf=5, n=1,m=100):
    q=1.6e-19
    ϵ0=8.84e-14
    ϵs = 11.4*ϵ0
    ϵox = 3.9*ϵ0
    #Ec0=0
    #Ecf=10
    Na=1e+16
    ni=1e+10
    tox=1e-7
    Cox = ϵox/tox
    ϕt = 0.026
    ϕms= 4.26 - 4.05
    Qox=0
    ϕf = ϕt*np.log(Na/ni)
    ΔECB=np.linspace(Ec0,Ecf,m)
    #v_in=-3
    #v_f=7
    Vgb = np.linspace(v_in,v_f,n)
    Ψg = []
    for vgb in Vgb:
        Gate=[]
        ψ_y=[]
        dis=[]
        Ψ=[]
        
        
        γ = np.sqrt(2*q*ϵs*Na)/(Cox)
        VFB = ϕms - (Qox/Cox)

        
        iter = []
        G=[]
    
            
        for Ecb in ΔECB:
            g=0
    
            if (vgb < 0):
                g1 = ( (-γ/2) + np.sqrt( abs((-γ/2)**2 - (VFB - abs(vgb)) )) )**2
        
                g2 = -2*ϕt*np.log((abs(vgb - VFB))/(γ*np.sqrt(ϕt)))
        
                g= min(g1,g2)
        
            else:
                g1 = ( (-γ/2) + np.sqrt( abs((-γ/2)**2 - (VFB - abs(vgb)) )) )**2
        
                g3 = Ecb + 2*( ϕf + ϕt*np.log((abs(abs(vgb) - VFB))/(γ*np.sqrt(ϕt))) )
        
                g=min(g1,g3)
        
        
            
        
    
        
            print(f"guess is {g}\n")
    
            f = lambda ψs : VFB + ψs + γ*np.sqrt( abs(ψs) + ϕt*np.exp((ψs - 2*ϕf - Ecb)/ϕt) -ϕt + ϕt*np.exp(-ψs/ϕt) - np.exp(-2*ϕf/ϕt)*( ϕt*np.exp(-Ecb/ϕt) + abs(ψs) )) - abs(vgb)
            y=Newton(f,g)
            Ψ.append(y[0])
            iter.append(y[1])
            G.append(g)
             
        
        Ψg.append(Ψ)
    
        return Ψ, ΔECB

