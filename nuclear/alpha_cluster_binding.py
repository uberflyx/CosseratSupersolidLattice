"""Figure: larger nuclei as close-packed alpha clusters on FCC.
(a) Hafstad-Teller law with FCC bond counts; (b) the binding curve."""
import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
mpl.rcParams.update({"font.size":10.5,"axes.linewidth":0.8,"figure.dpi":140,"savefig.dpi":140})

# n, name, A, Z, B_meas (AME2020, MeV), FCC bond count p
data = [
 (1,"$^4$He",4,2,28.296,0),(2,"$^8$Be",8,4,56.500,1),(3,"$^{12}$C",12,6,92.162,3),
 (4,"$^{16}$O",16,8,127.619,6),(5,"$^{20}$Ne",20,10,160.645,8),(6,"$^{24}$Mg",24,12,198.257,11),
 (7,"$^{28}$Si",28,14,236.537,15),(8,"$^{32}$S",32,16,271.781,18),
 (9,"$^{36}$Ar",36,18,306.716,21),(10,"$^{40}$Ca",40,20,342.052,24)]
Ba=28.296
n=np.array([d[0] for d in data]); A=np.array([d[2] for d in data])
B=np.array([d[4] for d in data]); p=np.array([d[5] for d in data])
extra=B-n*Ba
m=p>0
slope=(p[m]@extra[m])/(p[m]@p[m])

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10.4,4.3))

# (a) Hafstad-Teller line
ax1.plot(p[m],extra[m],'o',ms=6,color="#1b4965",zorder=3)
pp=np.linspace(0,p.max()*1.05,50)
ax1.plot(pp,slope*pp,'-',color="#ee6c4d",lw=1.8,
         label=f"one dock per bond\n$B_{{\\rm bond}}={slope:.2f}$ MeV")
ax1.axhline(0,color='0.7',lw=0.6)
for d in data:
    if d[5]>0: ax1.annotate(d[1],(d[5],d[4]-d[0]*Ba),textcoords="offset points",
                            xytext=(6,-3),fontsize=8.5)
ax1.set_xlabel("FCC alpha-alpha bond count  $p$")
ax1.set_ylabel(r"binding beyond free alphas  $B-nB_\alpha$  (MeV)")
ax1.set_title("(a)  Bond-counting law from FCC packing")
ax1.legend(loc="upper left",fontsize=9,frameon=False)
ax1.grid(alpha=0.25,lw=0.5)

# (b) binding curve B/A vs A
B_pred=n*Ba+p*slope
ax2.plot(A,B/A,'o',ms=6,color="#1b4965",label="measured (AME2020)",zorder=3)
ax2.plot(A,B_pred/A,'s',ms=5.5,mfc='none',color="#ee6c4d",mew=1.4,
         label="framework: $nB_\\alpha+pB_{\\rm bond}$")
ax2.set_xlabel("mass number  $A$")
ax2.set_ylabel(r"binding per nucleon  $B/A$  (MeV)")
ax2.set_title("(b)  The rising binding curve, reproduced")
ax2.legend(loc="lower right",fontsize=9,frameon=False)
ax2.grid(alpha=0.25,lw=0.5)
ax2.set_ylim(6.6,9.0)

fig.tight_layout()
fig.savefig("alpha_cluster_binding.pdf",bbox_inches="tight")
# report deviations
print("nucleus   B_meas   B_pred   dev(MeV)   dev(%)")
for d,bp in zip(data,B_pred):
    print(f"{d[1]:<9} {d[4]:>7.2f} {bp:>7.2f}  {bp-d[4]:>+7.2f}   {(bp-d[4])/d[4]*100:>+5.2f}")
print("\nsaved alpha_cluster_binding.pdf")
