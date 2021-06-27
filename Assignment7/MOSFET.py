#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.integrate as sci, matplotlib.pyplot as plt, numpy as np


# ## Poisson's Equation for MOSFET:
# 
# - $\frac{d^{2}\psi}{dy^{2}} = - \frac{qN_{a}}{\epsilon_{s}}[( e^{-\frac{\psi (y)}{\phi_{t}}} - 1) - e^{-\frac{2\phi_{f}}{\phi_{t}}}( e^{\frac{\psi (y)\ -\ \Delta E_{CB}}{\phi_{t}}} - 1 )  ] $

# In[2]:


def poissonMOS(Na=1e+16,Vgs1=1, tox=1e-7,Vds=5,ΔEcb=1 ):
    q=1.6e-19
    eps_0=8.84e-14
    eps_s = 11.4*eps_0
    eps_ox = 3.9*eps_0
    #Na = 1e+16
    ni = 1e+10
    tsi = 100e-7
    phi_ms = 4.26 - 4.05
    #Vgs1 = np.linspace(-1,1.5,15)
    Vgs = Vgs1 - phi_ms
    #ΔEcb = linspace(0,Vds,10)
    #tox = 1e-7
    Cox = eps_ox/tox
    phi_t = 0.026
    phi_f = phi_t*np.log(Na/ni)
    tsidep=np.sqrt((2*eps_s*2*phi_f/(q*Na)))
    t = np.linspace(0.0 , 500e-7, 1000)
    tspan=(0,tsidep)
    
    def poisson(x,y):
        return np.vstack((y[1], -(q*Na/eps_s)*( np.exp(-y[0]/phi_t) - 1 - np.exp(-2*(phi_f/phi_t))*( np.exp((y[0] - ΔEcb)/phi_t) - 1 ) )))
    
    def bc(ya,yb):
        return [ya[1] + (Cox/(eps_s))*( Vgs - ya[0] ), yb[0]]
    
    y_a0=0
    y_b0=0
    
    y=[-y_a0*np.ones(t.size), y_b0*np.ones(t.size)]
    
    
    sol = sci.solve_bvp(poisson, bc,t,y,tol=1e-4,max_nodes=1000000, verbose=2 )
    
    return sol.x, sol.y


# In[4]:


if __name__=="__main__":
    Gate=[]
    ψ_y=[]
    dis=[]
    ϕt=0.026
    Na=1e+16
    ni=1e+10
    ϕf = ϕt*np.log(Na/ni)
    v_in=-4
    
    
    vgs=np.linspace(v_in,4,100)
    plt.figure()
    for Vgs1 in vgs:
        y=poissonMOS(Vgs1=Vgs1)
        Gate.append(y[1][0][0])
        ψ_y.append(y[1][0])
        dis.append(y[0])
        plt.plot(y[0],y[1][0],label=f"Vgs={Vgs1}")
        
    plt.ylabel("Surface potential")
    plt.ylabel("Depth (um)")
    plt.legend()
    plt.show()

    plt.plot(vgs,Gate, label=f"ΔEcb = 1 V")
    plt.legend()


# In[ ]:




