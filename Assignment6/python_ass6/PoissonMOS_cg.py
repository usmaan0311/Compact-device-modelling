#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.integrate as sci, matplotlib.pyplot as plt, numpy as np


# In[2]:


def poissonMOS(Na=1e+16,Vgs1=1, tox=1e-7 ):
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
    #tox = 1e-7
    Cox = eps_ox/tox
    phi_t = 0.026
    phi_f = phi_t*np.log(Na/ni)
    tsidep=np.sqrt((2*eps_s*2*phi_f/(q*Na)))
    t = np.linspace(0.0 , 500e-7, 1000)
    tspan=(0,tsidep)
    
    def poisson(x,y):
        return np.vstack((y[1], -(q*Na/eps_s)*( np.exp(-y[0]/phi_t) - 1 - np.exp(-2*(phi_f/phi_t))*( np.exp(y[0]/phi_t) - 1 ) )))
    
    def bc(ya,yb):
        return [ya[1] + (Cox/(eps_s))*( Vgs - ya[0] ), yb[0]]
    
    y_a0=0
    y_b0=0
    
    y=[-y_a0*np.ones(t.size), y_b0*np.ones(t.size)]
    
    
    sol = sci.solve_bvp(poisson, bc,t,y,tol=1e-4,max_nodes=1000000, verbose=2 )
    
    return sol.x, sol.y


# In[30]:


if __name__!="__main__":
    Gate=[]
    
    vgs=np.linspace(-3,4,100)
    plt.figure()
    for Vgs1 in vgs:
        y=poissonMOS(Vgs1=Vgs1)
        Gate.append(y[1][0][0])
        plt.plot(y[0],y[1][0],label=f"{Vgs1}")
        
    plt.ylabel("Surface potential")
    plt.ylabel("Depth (um)")
    plt.legend()
    plt.show()


# In[31]:

vgs=np.linspace(-3,4,100)
plt.plot(vgs,Gate)


# In[32]:


gate=np.asarray(Gate)
ψs=np.average(gate[-6:-1])


# In[33]:


print(f"ψs = {np.round(ψs,3)} V")


# ## Low Frequency Capacitance $C_{gg}$
# $C_{gg} = C_{ox}( 1 - \frac{\partial\psi}{\partial V_{gs}})$

# In[38]:


def LowC(tox=1e-7):
    dψ=np.diff(gate)
    dVgs=np.diff(vgs)
    dψ_dVgs=dψ/dVgs
    eps_0=8.84e-14
    eps_ox = 3.9*eps_0
    Cox = eps_ox/tox
    Cgg = Cox*( 1 - dψ_dVgs)
    return Cgg
    



# In[39]:


Cgg=LowC()


# In[40]:


plt.plot(vgs[0:-1],Cgg,'g',linewidth=5, label="Low frequency")
plt.ylabel('$C_{gg}$ $ (\dfrac{F}{cm^{2}})$')
plt.xlabel('$V_{g}$ $ (V)$')
plt.legend()
plt.savefig('Cgg.jpeg', dpi=1200, bbox_inches='tight')


# In[ ]:




