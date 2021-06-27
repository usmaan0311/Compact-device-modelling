#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import trapz
from MOSFET import *


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
W=500e-7
L=r*W
β = μ*Cox*(W/L)
ϕt=0.026
ϕf=ϕt*np.log(Na/ni)
ϕms=4.26 - 4.05
Qox=0

Vfb = ϕms - Qox/Cox

γ=np.sqrt(2*q*ϵs*Na)/(Cox)


# # $I_{D}$ :
# 
# 
# * $ I_{D} = -\mu(\frac{W}{L})C_{ox}(q_{im} - \alpha_{m}\phi_{t})\phi$
# 
# * $\phi = \phi_{sd} - \phi_{ss}$
# * $ \alpha_{m} = 1 + \frac{\gamma}{2\sqrt{\phi_{m} - \phi_{t}}}$
# 
# * $ q_{im} = -\frac{\gamma\phi_{t}\Delta(\phi_{m},\zeta_{m})}{\sqrt{\phi_{m}\ -\ \phi_{t}\ +\ \phi_{t}\Delta(\phi_{m},\zeta_{m})}\ +\ \sqrt{(\phi_{m}\ -\ \phi_{t})}}$
# 
# * $\phi_{m} = \frac{(\phi_{ss}\ +\ \phi_{sd})}{2}$
# 
# * $\Delta(\phi_{m},\zeta_{m}) = \frac{1}{2}[\Delta(\phi_{ss},V_{sb})\ +\ \Delta(\phi_{sd},V_{sb}\ +\ V_{ds} )]\ -\ \frac{\phi^{2}}{4\gamma^{2}\phi_{t}} $
# 
# * $\Delta(\phi_{s},\zeta) = e^{\frac{(\phi_{s}\ -\ 2\phi_{f}\ -\ \zeta)}{\phi_{t}}}$
# 
# # $Q_{i}$, $\phi_{s}$ :
# 
# * $\phi_{s} = \phi_{m}\ +\ H[1\ -\ \sqrt{1\ -\ (\frac{2\phi}{HL})(y\ -\ y_{m}) }]$
# 
# * $H = \phi_{t}\ -\ \frac{q_{im}}{\alpha_{m}}$
# 
# * $y_{m} = (\frac{L}{2})(1\ +\ \frac{\phi}{4H})$
# 
# * $q_{i} = q_{im}\ +\ \alpha_{m}(\phi_{s}\ -\ \phi_{m})$
# 

# In[3]:


def Φ(vg,vd, vs=0):
    #x,y=poissonMOS(Vgs1=vg, ΔEcb=vd )
    
    ϕss = poissonMOS(Vgs1=vg, ΔEcb=0 )[1][0][0]
    
    ϕsd = poissonMOS(Vgs1=vg, ΔEcb=vd )[1][0][0]
    
    ϕ = ϕsd - ϕss
    
    ϕm = (ϕss + ϕsd)/2
    
    αm = 1 + γ/( 2*np.sqrt(abs(ϕm - ϕt)) )
    
    Δm = 0.5*( Δ(ϕss, vs) + Δ(ϕsd, vs + vd) ) - (ϕ**2)/(4*(γ**2)*ϕt)
    
    qim = - ( γ*ϕt*Δm )/( np.sqrt( abs(ϕm - ϕt + ϕt*Δm) ) + np.sqrt( abs(ϕm - ϕt ) )  )
    
    
    H = ϕt - (qim/αm)
    
    ym = (L/2)*(1 + (ϕ/(4*H)) )
    
    Qg = -qim + γ*np.sqrt(ϕm - ϕt) + ((ϕ**2)/(12*H))
    Qd = 0.5*( qim + (αm*ϕ/6)*(1 - ( ϕ/(2*H) ) - (ϕ**2)/(20*(H**2) ) )  )
    Qb = -γ*np.sqrt( ϕm - ϕt ) - ( (1 - αm)*(ϕ**2) )/(12*H)
    Qs = -Qg - Qd - Qb
    
    return ϕ, ϕm, αm, Δm, qim, H, ym, Qg, Qd, Qb, Qs



def Δ(ϕs,ζ):
    
    δ = np.exp( (ϕs - 2*ϕf - ζ)/ϕt )
    return δ

 


# In[4]:


Vgs= np.linspace(Vfb + 0.2, 5, 6)

Vds = np.linspace(0,5,40)

for vg in Vgs:
    
    ID=[]
    for vd in Vds:
        
        ϕ, αm, qim = Φ(vg=vg, vd=vd)[0], Φ(vg=vg, vd=vd)[2], Φ(vg=vg, vd=vd)[4]
        
        Id = -β*(qim - αm*ϕt)*ϕ
        
        ID.append(Id)
        
    plt.plot(Vds,ID, label=f"Vgs = {round(vg,2)}")
    
plt.legend()
plt.show()
        
        
        
        
        


# In[5]:


Vgs= np.linspace(Vfb + 0.2, 5, 40)

Vds = np.linspace(0,5,6)

fig,ax = plt.subplots()
ax1=ax.twinx()

for vd in Vds:
    
    ID=[]
    for vg in Vgs:
        
        ϕ, αm, qim = Φ(vg=vg, vd=vd)[0], Φ(vg=vg, vd=vd)[2], Φ(vg=vg, vd=vd)[4]
        
        Id = -β*(qim - αm*ϕt)*ϕ
        
        ID.append(Id)
        
    ax.semilogy(Vgs,ID, label=f"Vds = {round(vd,2)}")
    ax1.plot(Vgs,ID, label=f"Vds = {round(vd,2)}")
    
plt.legend()
plt.show()


# In[6]:


y = np.linspace(0,L,100)
Vgs= np.linspace(Vfb + 0.2, 5, 6)
vds=4

for vg in Vgs:
    ϕ, ϕm, αm, qim, H, ym=Φ(vg,vds, vs=0)[0], Φ(vg,vds, vs=0)[1], Φ(vg,vds, vs=0)[2], Φ(vg,vds, vs=0)[4], Φ(vg,vds, vs=0)[5], Φ(vg,vds, vs=0)[6] 
    
    ψs = ϕm + H*( 1 - np.sqrt( 1 - (2*ϕ/(H*L))*(y - ym)  ) )
    
    qi = qim 
    
    plt.plot(y,ψs, label=f"Vgs = {round(vg,2)}")
    
plt.legend()
plt.show()


# In[7]:


y = np.linspace(0,L,100)
Vgs= np.linspace(Vfb + 0.2, 5, 6)
vds=4

for vg in Vgs:
    ϕ, ϕm, qim, H, ym=Φ(vg,vds, vs=0)[0], Φ(vg,vds, vs=0)[1], Φ(vg,vds, vs=0)[4], Φ(vg,vds, vs=0)[5], Φ(vg,vds, vs=0)[6] 
    
    ψs = ϕm + H*( 1 - np.sqrt( 1 - (2*ϕ/(H*L))*(y - ym)  ) )
    
    qi = qim + αm*(ψs - ϕm)
    
    plt.plot(y,abs(qi), label=f"Vgs = {round(vg,2)}")
    
plt.legend()
plt.show()


# # Transcapacitances:
# 
# * $Q_{G} = -q_{im}\ +\ \gamma\sqrt{(\phi_{m}\ -\ \phi_{t})}\ +\ (\frac{\phi^{2}}{12H})$
# 
# * $Q_{D} = \frac{1}{2}[ q_{im}\ +\ \frac{\alpha_{m}\phi}{6}(1\ -\ \frac{\phi}{2H}\ -\ \frac{\phi^{2}}{20H^{2}}) ] $
# 
# * $Q_{B} = -\gamma\sqrt{\phi_{m}\ -\ \phi_{t}}\ -\ \frac{(1\ -\ \alpha_{m})\phi^{2}}{12H}$
# 
# * $Q_{S} = -Q_{G}\ -\ Q_{D}\ -\ Q_{B}$
# 
# * $C_{ij} = (2\delta_{ij}\ -\ 1)(\frac{\partial Q_{i}}{\partial V_{j}})\ ,\ where\ i,\ j\ =\ G,\ S,\ D,\ B $

# In[ ]:


Vgs= np.linspace(Vfb + 0.2, 5, 40)

vds = 2.0 # np.linspace(0,5,6)
QG=[]
QD=[]
QB=[]
QS=[]
for vg in Vgs:
    y=Φ(vg,vds)
    
    Qg = y[7] 
    Qd = y[8] 
    Qb = y[9] 
    Qs = y[10]
    
    QG.append(Qg)
    QD.append(Qd)
    QB.append(Qb)
    #QS.append(Qs)

QG = np.asarray(QG)
QD = np.asarray(QD)
QB = np.asarray(QB)
#QS = np.asarray(QS)



dVg = np.diff(Vgs)
dQg = np.diff(QG)
dQd = np.diff(QD)
dQb = np.diff(QB)
#dQs = np.diff(QS)

Cgg = dQg/dVg
Cgd = -dQd/dVg
Cgb = -dQb/dVg
#Cgs = -dQs/dVg



plt.plot(Vgs[:-1], Cgg, label="$C_{gg}$")
plt.plot(Vgs[:-1], Cgd, label="$C_{gd}$")
plt.plot(Vgs[:-1], Cgb, label="$C_{gb}$")
#plt.plot(Vgs[:-1], Cgs, label="$C_{gs}$")

plt.legend()
plt.show()
    


# In[17]:


#Vgs= np.linspace(Vfb + 0.2, 5, 40)

#vds = 2.0 # np.linspace(0,5,6)
#QB=[]
#for vg in Vgs:
#    y=Φ(vg,vds)[9]
#    
#    QB.append(y)
    
#QB = np.asarray(QB)
#dQb = np.diff(QB)
#Cgb = -dQb/dVg
#
#plt.plot(Vgs[:-1], Cgb, label="$C_{gb}$")


# In[16]:


#plt.plot(Vgs,QB)


# In[ ]:




