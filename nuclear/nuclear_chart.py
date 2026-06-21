"""
The nuclear chart from the assembled mass formula.

Cosserat Supersolid Lattice framework. Plots the binding-per-nucleon curve (and
its iron peak) and the valley of stability from the four terms of the
semi-empirical mass formula. The Coulomb coefficient is the framework's
(3/5) * alpha * hbar c / r0; the volume, surface, and asymmetry coefficients use
standard measured values here, since a_V and a_S are the pieces the framework has
yet to derive from the saturation density. Produces nuclear_chart.pdf.
"""
import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
mpl.rcParams.update({"font.size":10.5,"axes.linewidth":0.8,"figure.dpi":140,"savefig.dpi":140})
alpha=1/137.035999; hbarc=197.3269; r0=1.20
a_V,a_S,a_C,a_A,a_P=15.75,17.8,0.6*alpha*hbarc/r0,23.7,11.2
def B(Z,A):
    N=A-Z; d=0.0
    if Z%2==0 and N%2==0: d=a_P/np.sqrt(A)
    elif Z%2==1 and N%2==1: d=-a_P/np.sqrt(A)
    return a_V*A-a_S*A**(2/3)-a_C*Z*(Z-1)/A**(1/3)-a_A*(A-2*Z)**2/A+d
def Zstar(A):
    zs=np.arange(1,A); return zs[np.argmax([B(z,A) for z in zs])]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10.6,4.4))

# (a) binding curve
As=np.arange(4,252)
BA=np.array([B(Zstar(A),A)/A for A in As])
ax1.plot(As,BA,'-',color="#ee6c4d",lw=1.8,label="lattice mass formula",zorder=2)
meas={4:7.074,12:7.680,16:7.976,40:8.551,56:8.790,62:8.795,90:8.710,
      120:8.505,150:8.262,180:8.022,208:7.867,238:7.570}
ax1.plot(list(meas),list(meas.values()),'o',ms=5,color="#1b4965",
         label="measured (AME2020)",zorder=3)
ax1.axvline(58,color='0.6',ls=':',lw=1)
ax1.text(64,7.0,"iron-nickel\npeak",fontsize=9,color='0.35')
ax1.set_xlabel("mass number  $A$"); ax1.set_ylabel(r"binding per nucleon  $B/A$  (MeV)")
ax1.set_title("(a)  Binding curve and the iron peak")
ax1.legend(loc="lower right",fontsize=9,frameon=False); ax1.grid(alpha=0.25,lw=0.5)
ax1.set_ylim(6.8,9.1)

# (b) valley of stability (Segre chart)
Av=np.arange(4,252,1)
Zv=np.array([Zstar(A) for A in Av]); Nv=Av-Zv
ax2.plot(Nv,Zv,'-',color="#ee6c4d",lw=2,label=r"lattice $Z^\star(A)$",zorder=2)
ax2.plot([0,150],[0,150],'--',color='0.6',lw=1,label="$N=Z$")
stable={"O":(8,8),"Ca":(20,20),"Fe":(30,26),"Zr":(50,40),"Sn":(70,50),
        "Ce":(82,58),"Hf":(106,72),"Pb":(126,82),"U":(146,92)}
for nm,(N,Z) in stable.items():
    ax2.plot(N,Z,'o',ms=6,color="#1b4965",zorder=3)
    ax2.annotate(nm,(N,Z),textcoords="offset points",xytext=(5,-9),fontsize=8.5)
ax2.set_xlabel("neutron number  $N$"); ax2.set_ylabel("proton number  $Z$")
ax2.set_title("(b)  The valley of stability")
ax2.legend(loc="upper left",fontsize=9,frameon=False); ax2.grid(alpha=0.25,lw=0.5)
ax2.set_xlim(0,150); ax2.set_ylim(0,100)

fig.tight_layout(); fig.savefig("nuclear_chart.pdf",bbox_inches="tight")
print("saved nuclear_chart.pdf")
