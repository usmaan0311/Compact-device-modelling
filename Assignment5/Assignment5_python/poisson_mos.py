#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.integrate as sci, matplotlib.pyplot as plt, numpy as np


# In[2]:


q=1.6e-19
eps_0=8.84e-14
eps_s = 11.4*eps_0
eps_ox = 3.9*eps_0
Na = 1e+16
ni = 1e+10
tsi = 100e-7
phi_ms = 4.26 - 4.05
Vgs1 = np.linspace(-1,1.5,15)
Vgs =1 #Vgs1 - phi_ms
tox = 1e-7
Cox = eps_ox/tox
phi_t = 0.026
phi_f = phi_t*np.log(Na/ni)
tsidep=np.sqrt((2*eps_s*2*phi_f/(q*Na)))
t = np.linspace(0.0 , 100e-7, 1000000)
tspan=(0,tsidep)


# ## Poisson's Equation
# 
# $ \frac{d^{2}\psi}{dy^{2}} = -\frac{qN_{A}}{\epsilon_{s}}[ e^{-\psi (y)/\phi_{t}} - 1 - e^{-2\phi_{F}/\phi_{t}}( e^{\psi (y)/\phi_{t}} - 1 ) ]$
# 
# 
# $\frac{d\psi}{dy} = y1$
# 
# $y1^{'} = -\frac{qN_{A}}{\epsilon_{s}}[ e^{-\psi (y)/\phi_{t}} - 1 - e^{-2\phi_{F}/\phi_{t}}( e^{\psi (y)/\phi_{t}} - 1 ) ]$
# 
# 
# 
# $ \begin{equation}
# S (y) =
# \begin{bmatrix}
# \psi (y)  \\
# \psi^{'} (y)  
# \end{bmatrix}
# \end{equation}$
# 
# 
# 
# $ \begin{equation}
# \frac{dS (y)}{dy} =
# \begin{bmatrix}
# \psi^{'} (y)  \\
# \psi^{''} (y)  
# \end{bmatrix}
# \end{equation} $  $ \begin{equation}
# =
# \begin{bmatrix}
# \psi^{'} (y)  \\
# -\frac{qN_{A}}{\epsilon_{s}}[ e^{-\psi (y)/\phi_{t}} - 1 - e^{-2\phi_{F}/\phi_{t}}( e^{\psi (y)/\phi_{t}} - 1 ) ]  
# \end{bmatrix}
# \end{equation} $
# 
# $ \begin{equation}
# \frac{dS (y)}{dy} =
# \begin{bmatrix}
# S_{2} (y)  \\
# -\frac{qN_{A}}{\epsilon_{s}}[ e^{-S_{1} (y)/\phi_{t}} - 1 - e^{-2\phi_{F}/\phi_{y}}( e^{S_{1} (t)/\phi_{t}} - 1 ) ]  
# \end{bmatrix}
# \end{equation} $
# 
# $ \begin{equation}
# \frac{dS (y)}{dy} =
# \begin{bmatrix}
# y[1]  \\
# -\frac{qN_{A}}{\epsilon_{s}}[ e^{-y[0]/\phi_{t}} - 1 - e^{-2\phi_{F}/\phi_{y}}( e^{y[0]/\phi_{t}} - 1 ) ]  
# \end{bmatrix}
# \end{equation} $
# 
# 
# 
# 
# 
 

 

# In[3]:


def poisson(x,y):
    return np.vstack((y[1], -(q*Na/eps_s)*( np.exp(-y[0]/phi_t) - 1 - np.exp(-2*(phi_f/phi_t))*( np.exp(y[0]/phi_t) - 1 ) )))


# In[4]:


#def poisson(x,y):
 #   dy= y1
  #  
   # dy1=-(q*Na/eps_s)*( np.exp(-y[/phi_t) - 1 - np.exp(-2*(phi_f/phi_t))*( np.exp(y[0]/phi_t) - 1 ) )


# In[5]:


def bc(ya,yb):
    return [ya[-1]-0, yb[0] + (Cox/(eps_s))*( Vgs - ya[0] )]


# In[6]:


#def bc(ya,yb):
 #   return ya[1], yb[0] + (Cox/(eps_s))*( vgs - ya[0] )


# In[7]:


y_a0=0
y_b0=0


# In[8]:


#y_a = y_a0*np.ones((2,t.size))
#y_b = y_b0*np.ones((2,t.size))
#y0=[y_a,y_b]
y=[-y_a0*np.ones(t.size), y_b0*np.ones(t.size)]
y[0][0]=-0.5
y[1][0]=0



# In[9]:


sol = sci.solve_bvp(poisson, bc,t,y,tol=1e-10, max_nodes=100000,verbose=2)


# In[10]:


print(sol.message)


# In[11]:


sol.success


# In[12]:


y_plot = sol.sol(t)[0]


# In[13]:


plt.plot(t,y_plot)


# In[ ]:




