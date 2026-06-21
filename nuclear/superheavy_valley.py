"""
The beta-stability valley extended into the superheavies.

Cosserat Supersolid Lattice framework. Finds the most beta-stable isotope of each
superheavy element by locating the mass number A whose most-bound Z (over the
assembled mass formula) equals the target element, then compares the lattice valley
to the standard liquid-drop SEMF, to the synthesised nuclei, and to the shell-magic
island at N=184. The Coulomb coefficient is the framework's (3/5) alpha hbar c / r0;
a_V, a_S, a_A use measured values, the pieces the framework has yet to derive.
"""
import numpy as np
alpha=1/137.035999; hbarc=197.3269; r0=1.20
aC=0.6*alpha*hbarc/r0
aV,aS,aA,aP=15.75,17.8,23.7,11.2
def B(Z,A,c=(aV,aS,aA,aC,aP)):
    aV,aS,aA,aC,aP=c; N=A-Z; d=0.0
    if Z%2==0 and N%2==0: d=aP/np.sqrt(A)
    elif Z%2==1 and N%2==1: d=-aP/np.sqrt(A)
    return aV*A-aS*A**(2/3)-aC*Z*(Z-1)/A**(1/3)-aA*(A-2*Z)**2/A+d
def Zstar(A,c=(aV,aS,aA,aC,aP)):
    zs=np.arange(1,A); return int(zs[np.argmax([B(z,A,c) for z in zs])])
def valleyA(Ztarget,c=(aV,aS,aA,aC,aP)):
    """A whose beta-stable Z equals Ztarget (the valley isotope of that element)"""
    best,bestA=999,None
    for A in range(2*Ztarget, 3*Ztarget+60):
        if abs(Zstar(A,c)-Ztarget)<best:
            best=abs(Zstar(A,c)-Ztarget); bestA=A
    return bestA
semf=(15.75,17.8,23.7,0.711,11.2)
synth={104:267,106:269,108:270,110:281,112:285,114:289,115:290,116:293,117:294,118:294}
elem={104:"Rf",106:"Sg",108:"Hs",110:"Ds",112:"Cn",114:"Fl",115:"Mc",116:"Lv",117:"Ts",118:"Og",120:"Ubn",126:"Uhh"}
print("Beta-stability valley into the superheavies (corrected)")
print(f"{'Z':>4}{'el':>5}{'A* lattice':>11}{'N* lat':>7}{'A* SEMF':>8}"
      f"{'synth A':>9}{'synth N':>8}{'island N=184 A':>15}")
for Z in [92,104,106,108,110,112,114,115,116,117,118,120,126]:
    Al=valleyA(Z); As=valleyA(Z,semf); s=synth.get(Z)
    sa=f"{s}" if s else "-"; sn=f"{s-Z}" if s else "-"
    print(f"{Z:>4}{elem.get(Z,'U'):>5}{Al:>11}{Al-Z:>7}{As:>8}{sa:>9}{sn:>8}{184+Z:>15}")
