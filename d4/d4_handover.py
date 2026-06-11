#!/usr/bin/env python3
"""
d4_handover.py - the tangential contact-handover law of a {111} glide wall
==========================================================================

Tests the handover law derived from the lossless-rolling principle:

  (1) RATE-BASED: rolling is a no-slip condition on relative surface
      VELOCITY, so the tangential spring stores the path-integral of the
      slip that rolling cannot remove (route-dependent).
  (2) LOSSLESS FADE: a contact breaking while still carrying tangential
      force would radiate a dissipative snap, forbidden by the same
      principle that justifies rolling; so k_t must fade continuously to
      zero as the contact opens.

With the fade the misfit profile gamma(f) becomes single-valued and free
of the spurious contact events a hard cut-off produces, and the faulted
endpoint carries no stored shear. The bracket test below shows the
CHARACTER is robust (clean, monotonic gamma for every fade) while the
VALUE still spans ~0.99-1.12: the fade RATE - the tangential bond's
load-dependence, the shear-sector analogue of the radial xi = -7 - is
the one remaining microscopic input, and it is not fixed by the
gravitational (central-force) sector. The minimal single-junction choice
k_t propto k_n returns the Morse-tracking edge (~0.99, near Frenkel) and
makes the equilibrium ratio strain-independent (eta_t = 7).

Author: Mitchell A. Cox   Date: June 2026
"""

import numpy as np
from scipy.optimize import minimize

# ============================================================
# Testing the handover law derived from the lossless principle.
#
# CLAIM (to test): the rolling constraint is a RATE (non-holonomic)
# constraint -> the tangential spring stores the path-integral of
# FRUSTRATED slip, and lossless transmission FORCES k_t to fade
# continuously to zero as a contact opens (else breaking = a
# dissipative snap). The decisive prediction: with a smooth fade,
# the misfit surface gamma(f) becomes BRACKET-INDEPENDENT (no
# arbitrary s_break), curing the pathology the frozen/snap laws had.
#
# Test = compute gamma(f) for several fade choices; if gamma_max and
# gamma_SF agree across them, the law is robust (cracked).
# ============================================================
s=1/np.sqrt(2); w_hat=np.array([1,1,1,0])/np.sqrt(3); t_hat=np.array([1,1,-2,0])/np.sqrt(6)
dhop=1/np.sqrt(3); AM=7/3; R=0.5
Gt=np.outer(t_hat,w_hat)-np.outer(w_hat,t_hat)
Dm=1/(2*AM*AM)
def Vmorse(r): return Dm*(1-np.exp(-AM*(r-1)))**2
def kn_eff(r):  # normal stiffness V''(r), the natural engagement measure
    e=np.exp(-AM*(r-1)); return 2*Dm*AM*AM*(2*e*e-e)

# fade functions for k_t(r) = r_ratio * fade(r) * kn0 ; all -> 0 at opening, =1 at r=1
def fade_morse(r):   return np.clip(kn_eff(r)/kn_eff(1.0),0,None)        # tracks k_n
def fade_overlap(r): return np.clip(1-(r-1)/0.30,0,1)                    # linear, range 0.30
def fade_overlap2(r):return np.clip(1-(r-1)/0.45,0,1)                    # linear, range 0.45 (bracket test)
def fade_exp(r):     return np.clip(np.exp(-2*AM*(r-1)),0,1)*(r<1.6)     # exponential

def gen(box=5):
    out=[]
    for a in range(-box,box+1):
        for b in range(-box,box+1):
            for c in range(-box,box+1):
                m=a+b+c
                if m<-3 or m>3: continue
                for e in range(-box,box+1):
                    if (a+b+c+e)%2: continue
                    out.append((m,s*np.array([a,b,c,e])))
    return out
nodes=gen()
lowers=[(m,p) for m,p in nodes if m<=0]
reps={m:min((p for mm,p in nodes if mm==m), key=lambda p:np.linalg.norm(p-(p@w_hat)*w_hat)) for m in (1,2,3)}
MU=[];ML=[];V0=[]
for mu,rep in reps.items():
    for ml,p in lowers:
        v0=rep-p
        if np.linalg.norm(v0)<1.9: MU.append(mu);ML.append(ml);V0.append(v0)
MU=np.array(MU);ML=np.array(ML);V0=np.array(V0); Npair=len(V0)
ROT=[-2,-1,0,1,2]; li={m:k for k,m in enumerate(ROT)}
IU=np.array([li.get(m,-1) for m in MU]); IL=np.array([li.get(m,-1) for m in ML])

def run(r_ratio, fade, nf=41):
    fs=np.linspace(0,dhop,nf)
    XI=np.zeros((Npair,4)); th_p=np.zeros(len(ROT)); dh_p=0.0; f_p=0.0
    g=np.zeros(nf); E0=None; x=np.zeros(1+len(ROT))
    for k,f in enumerate(fs):
        XIc=XI.copy()
        def energy(x):
            dh=x[0]; th=x[1:]
            v=V0+(-f*t_hat+dh*w_hat); r=np.linalg.norm(v,axis=1)
            E=Vmorse(r).sum()
            n=v/r[:,None]
            kt=r_ratio*fade(r)*kn_eff(1.0)
            dU=(-(f-f_p))*t_hat+(dh-dh_p)*w_hat
            dUp=dU-(n@dU)[:,None]*n
            tu=np.where(IU>=0, th[np.maximum(IU,0)]-th_p[np.maximum(IU,0)],0.0)
            tl=np.where(IL>=0, th[np.maximum(IL,0)]-th_p[np.maximum(IL,0)],0.0)
            xin=XIc+dUp-R*(tu+tl)[:,None]*(n@Gt.T)
            E+=0.5*np.sum(kt*np.einsum('ij,ij->i',xin,xin))
            return E
        res=minimize(energy,x,method='Powell',options={'xtol':1e-9,'ftol':1e-11,'maxiter':3000})
        x=res.x
        if E0 is None: E0=res.fun
        g[k]=res.fun-E0
        dh=x[0]; th=x[1:]
        v=V0+(-f*t_hat+dh*w_hat); r=np.linalg.norm(v,axis=1); n=v/r[:,None]
        dU=(-(f-f_p))*t_hat+(dh-dh_p)*w_hat; dUp=dU-(n@dU)[:,None]*n
        tu=np.where(IU>=0, th[np.maximum(IU,0)]-th_p[np.maximum(IU,0)],0.0)
        tl=np.where(IL>=0, th[np.maximum(IL,0)]-th_p[np.maximum(IL,0)],0.0)
        XI=XIc+dUp-R*(tu+tl)[:,None]*(n@Gt.T)
        # reset where contact has opened (fade ~ 0)
        opened=fade(r)<1e-3; XI[opened]=0
        th_p=th.copy(); dh_p=dh; f_p=f
    return fs,g

r0=1/(3*np.pi-5)
print("BRACKET-INDEPENDENCE TEST (handover = rate-based + smooth fade):")
print(f"{'fade law':>22s}  {'gamma_max':>10s}  {'gamma_SF':>10s}  {'f_max/d':>8s}")
for name,fade in [("Morse-tracking",fade_morse),("linear range .30",fade_overlap),
                  ("linear range .45",fade_overlap2),("exponential",fade_exp)]:
    fs,g=run(r0,fade)
    gm=g.max(); gd=g[-1]; fm=fs[g.argmax()]/dhop
    print(f"{name:>22s}  {gm:10.4f}  {gd:10.4f}  {fm:8.2f}")
