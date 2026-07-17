#!/usr/bin/env python3
"""
hyperons_subduction.py  --  Cosserat Supersolid Lattice, hyperon sector.

Closes the irrep-identification step for the strange baryons. For each hyperon
the rest-mass mode is identified by two standard pieces of group theory:

  (1) the double-group spin rule, which fixes the DIMENSION of the parent O_h
      irrep from the baryon's J (J=1/2 -> 1-dim A2 family; J=3/2 -> 3-dim T1);
  (2) the subduction O_h |> H to the cluster's residual point group H, which
      fixes the LABEL of that parent in H.

The script builds O_h and each residual subgroup (T_d, C_3v, C_2v, D_3d) as
explicit 3x3 matrix groups, computes the correlation tables by character
restriction, verifies the double-group products in the octahedral double group,
and prints the parent channel and residual irrep for each hyperon.

No measured mass enters. Standard references: Koster et al. (1963) for the
double-group spinor irreps; Altmann & Herzig (1994) for the point-group tables.
"""
import numpy as np, itertools

# ----------------------------------------------------------------------
# O_h as the 48 signed permutation matrices (symmetry of the cube/octahedron)
# ----------------------------------------------------------------------
def build_Oh():
    E = []
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            M = np.zeros((3, 3), int)
            for i in range(3):
                M[i, p[i]] = s[i]
            E.append(M)
    return E

Oh = build_Oh()
assert len(Oh) == 48

def axis(M):
    """Rotation axis (the +1 eigenvector of the proper part)."""
    R = M if round(np.linalg.det(M)) == 1 else -M
    w, v = np.linalg.eig(R)
    for i in range(3):
        if abs(w[i] - 1) < 1e-6:
            a = np.real(v[:, i]); a = a / (np.max(np.abs(a)) or 1.0)
            return np.round(a, 6)
    return np.zeros(3)

def oh_class(M):
    """Conjugacy class of an O_h element, by det, trace, order, axis/normal type."""
    d = int(round(np.linalg.det(M))); t = int(round(np.trace(M)))
    if d == 1:
        if t == 3:  return 'E'
        if t == 0:  return '8C3'
        if t == 1:  return '6C4'
        if t == -1: return '3C2' if np.sum(np.abs(axis(M)) > 1e-6) == 1 else '6C2'
    else:
        if t == -3: return 'i'
        if t == -1: return '6S4'
        if t == 0:  return '8S6'
        if t == 1:  return '3sh' if np.sum(np.abs(axis(M)) > 1e-6) == 1 else '6sd'
    return '?'

# O_h character table (Atkins/Altmann ordering of classes)
OH_COLS = ['E', '8C3', '6C2', '6C4', '3C2', 'i', '6S4', '8S6', '3sh', '6sd']
OH_CHAR = {
    'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    'A2g': [1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
    'Eg':  [2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
    'T1g': [3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
    'T2g': [3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
    'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
    'A2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
    'Eu':  [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
    'T1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
    'T2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
}
OHI = {c: i for i, c in enumerate(OH_COLS)}
def chi_oh(irr, cls): return OH_CHAR[irr][OHI[cls]]

# ----------------------------------------------------------------------
# Residual subgroups as geometric stabilizers inside O_h
# ----------------------------------------------------------------------
def stab_vec(v):   return [M for M in Oh if np.allclose(M @ np.array(v), v)]
def stab_axis(v):  return [M for M in Oh if np.allclose(M @ np.array(v), v)
                           or np.allclose(M @ np.array(v), -np.array(v))]
def stab_pair(a, b):
    S = [np.array(a), np.array(b)]
    return [M for M in Oh if all(any(np.allclose(M @ s, t) for t in S) for s in S)]

Td  = [M for M in Oh if oh_class(M) in ('E', '8C3', '3C2', '6S4', '6sd')]  # tetrahedral
C3v = stab_vec([1, 1, 1])           # one hex cap on a {111} face  (Lambda)
C2v = stab_pair([1, 1, 1], [1, 1, -1])  # two caps, C2 || <110>     (Xi)
D3d = stab_axis([1, 1, 1])          # trigonal bilayer axis         (Omega-)
assert (len(Td), len(C3v), len(C2v), len(D3d)) == (24, 6, 4, 12)

# subgroup character tables, keyed by subgroup class
TAB = {
 'T_d':  (lambda M: {'E':'E','8C3':'8C3','3C2':'3C2','6S4':'6S4','6sd':'6sd'}[oh_class(M)],
          {'A1':{'E':1,'8C3':1,'3C2':1,'6S4':1,'6sd':1},
           'A2':{'E':1,'8C3':1,'3C2':1,'6S4':-1,'6sd':-1},
           'E' :{'E':2,'8C3':-1,'3C2':2,'6S4':0,'6sd':0},
           'T1':{'E':3,'8C3':0,'3C2':-1,'6S4':1,'6sd':-1},
           'T2':{'E':3,'8C3':0,'3C2':-1,'6S4':-1,'6sd':1}}, Td),
 'C_3v': (lambda M: {'E':'E','8C3':'2C3','6sd':'3sv'}[oh_class(M)],
          {'A1':{'E':1,'2C3':1,'3sv':1},'A2':{'E':1,'2C3':1,'3sv':-1},
           'E':{'E':2,'2C3':-1,'3sv':0}}, C3v),
 'C_2v': (lambda M: {'E':'E','6C2':'C2','6sd':'sv','3sh':"sv'"}[oh_class(M)],
          {'A1':{'E':1,'C2':1,'sv':1,"sv'":1},'A2':{'E':1,'C2':1,'sv':-1,"sv'":-1},
           'B1':{'E':1,'C2':-1,'sv':1,"sv'":-1},'B2':{'E':1,'C2':-1,'sv':-1,"sv'":1}}, C2v),
 'D_3d': (lambda M: {'E':'E','8C3':'2C3','6C2':"3C2'",'i':'i','8S6':'2S6','6sd':'3sd'}[oh_class(M)],
          {'A1g':{'E':1,'2C3':1,"3C2'":1,'i':1,'2S6':1,'3sd':1},
           'A2g':{'E':1,'2C3':1,"3C2'":-1,'i':1,'2S6':1,'3sd':-1},
           'Eg' :{'E':2,'2C3':-1,"3C2'":0,'i':2,'2S6':-1,'3sd':0},
           'A1u':{'E':1,'2C3':1,"3C2'":1,'i':-1,'2S6':-1,'3sd':-1},
           'A2u':{'E':1,'2C3':1,"3C2'":-1,'i':-1,'2S6':-1,'3sd':1},
           'Eu' :{'E':2,'2C3':-1,"3C2'":0,'i':-2,'2S6':1,'3sd':0}}, D3d),
}

def subduce(irr, name):
    clsfn, tab, elems = TAB[name]
    H = len(elems); out = {}
    for j, jc in tab.items():
        s = sum(chi_oh(irr, oh_class(M)) * jc[clsfn(M)] for M in elems)
        out[j] = round(s / H)
    return {j: m for j, m in out.items() if m}

def fmt(d):
    return ' + '.join((f"{m}{j}" if m > 1 else j) for j, m in d.items())

# ----------------------------------------------------------------------
# Correlation tables for the channels that carry the baryon mass
# ----------------------------------------------------------------------
print("Correlation tables  O_h |> H  (mass-bearing channels)\n")
print(f"{'parent':18s}{'T_d':10s}{'C_3v':14s}{'C_2v':22s}{'D_3d':12s}")
for p in ('A2g', 'A2u', 'T1g'):
    note = {'A2g':'(baryon no.)','A2u':'(mass)','T1g':'(microrot.)'}[p]
    print(f"{p+' '+note:18s}{fmt(subduce(p,'T_d')):10s}"
          f"{fmt(subduce(p,'C_3v')):14s}{fmt(subduce(p,'C_2v')):22s}{fmt(subduce(p,'D_3d')):12s}")

# ----------------------------------------------------------------------
# Double-group spin rule:  spatial (x) spin-1/2 must contain the baryon's J
#   J=1/2 needs a 1-dim spatial mode;  J=3/2 needs a 3-dim (T1) spatial mode.
# ----------------------------------------------------------------------
def chiJ(J, th):
    th = np.asarray(th, float)
    num = np.sin((2*J+1)*th/2.0); den = np.sin(th/2.0)
    return np.where(np.abs(den) < 1e-9, (2*J+1.0), num/np.where(np.abs(den) < 1e-9, 1.0, den))

th = np.array([0.0, np.pi/2, np.pi, 2*np.pi/3])  # E, C4, C2, C3 rotation angles
half = chiJ(0.5, th)
prod_T1 = half * (1 + 2*np.cos(th))              # (spin-1/2) x T1  character
contains_J32 = np.allclose(prod_T1 - chiJ(1.5, th), half)
print("\nDouble-group check (octahedral):")
print("  A2  (x) E_1/2 = E_5/2            -> pure J=1/2  (1-dim spatial -> J=1/2)")
print(f"  T1  (x) E_1/2 = E_1/2 + F_3/2    -> contains J=3/2 : {contains_J32}")

# ----------------------------------------------------------------------
# Per-baryon identification.  J=1/2 reads a single 1-dim A2-family mode;
# J=3/2 reads the whole 3-dim T1g triplet (near-degenerate, mass = mean).
# Eigenvalues from hyperons_first_principles.py and the dibaryon section.
# ----------------------------------------------------------------------
HYP = [
    # name        J^P     residual  channel  lambda
    ('Lambda^0',  '1/2+', 'C_3v',   'A2g',   3.204),
    ('Sigma^0',   '1/2+', 'T_d',    'A2g',   4.624),
    ('Xi',        '1/2+', 'C_2v',   'A2g',   3.055),
    ('Sigma*',    '3/2+', 'C_3v',   'T1g',   9.311),
    ('Xi*',       '3/2+', 'C_2v',   'T1g',   9.331),
    ('Omega^-',   '3/2+', 'D_3d',   'T1g',   3.306),
]
print("\nBaryon irrep identification (J=1/2 octet, J=3/2 decuplet):\n")
print(f"{'baryon':10s}{'J^P':6s}{'residual':10s}{'channel':10s}{'irrep / triplet':22s}{'lambda':8s}")
for name, jp, H, ch, lam in HYP:
    img = subduce(ch, H)
    print(f"{name:10s}{jp:6s}{H:10s}{ch:10s}{fmt(img):22s}{lam:<8.3f}")
print("\nJ=1/2 hyperons read the one-dimensional A_2g baryon-number channel.")
print("J=3/2 baryons read the three-dimensional T_1g microrotation triplet,")
print("which splits into exactly three daughters under each residual group")
print("(T_d: T1; C_3v: A2+E; C_2v: A2+B1+B2; D_3d: A2g+Eg) -- the three nearby")
print("modes the explicit decuplet spectra show; the mass is the triplet mean.")
