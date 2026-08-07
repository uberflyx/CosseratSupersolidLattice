import numpy as np, itertools

# ---- O_h as 48 signed permutation matrices (FCC point group) ----
def build_Oh():
    E=[]
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1,-1),repeat=3):
            M=np.zeros((3,3),int)
            for i in range(3): M[i,p[i]]=s[i]
            E.append(M)
    return E
Oh=build_Oh()
key=lambda M: M.tobytes()
idx={key(M):i for i,M in enumerate(Oh)}
def mul(A,B): return A@B

def stab_vec(v):
    v=np.array(v); return [M for M in Oh if np.allclose(M@v,v)]

# cuboctahedron feature stabilizers (the building-block attachment sites)
feat = {
 'C2v (vertex, <110>)' : stab_vec([1,1,0]),   # 12 vertices  -> O_h/C2v = cuboctahedron
 'C3v (tri-face, <111>)': stab_vec([1,1,1]),   # 8 triangular faces (hex-cap / void site)
 'C4v (sq-face, <100>)' : stab_vec([0,0,1]),   # 6 square faces
}
for name,H in feat.items():
    print(f"{name:24s} order {len(H):2d}   index in O_h = {48//len(H)}")

def double_cosets(H,K):
    Hset=[h for h in H]; Kset=[k for k in K]
    seen=set(); reps=[]; sizes=[]
    for g in Oh:
        if key(g) in seen: continue
        dc=set()
        for h in Hset:
            hg=mul(h,g)
            for k in Kset:
                dc.add(key(mul(hg,k)))
        reps.append(g); sizes.append(len(dc)); seen|=dc
    return reps,sizes

print("\nDouble cosets  H \\ O_h / K   (inequivalent fusions of two clusters along the named features):\n")
names=list(feat); 
print(f"{'feature A':24s}{'feature B':24s}{'# double cosets':16s}{'orbit sizes (degeneracies)'}")
results={}
for i,a in enumerate(names):
    for j,b in enumerate(names):
        if j<i: continue
        reps,sizes=double_cosets(feat[a],feat[b])
        results[(a,b)]=(len(reps),sorted(sizes,reverse=True))
        print(f"{a:24s}{b:24s}{len(reps):<16d}{sorted(sizes,reverse=True)}  (sum {sum(sizes)})")

# sanity: |H\\G/K| equals number of H-orbits on G/K (Burnside cross-check for one pair)
def orbits_H_on_GmodK(H,K):
    # cosets gK
    cosets=[]; seen=set()
    for g in Oh:
        c=frozenset(key(mul(g,k)) for k in K)
        if next(iter(c)) not in seen:
            cosets.append(c); seen|=c
    # H acts on cosets by left mult; count orbits
    cset=[c for c in cosets]
    rep_of={}
    for ci,c in enumerate(cset):
        for e in c: rep_of[e]=ci
    seen=set(); norb=0
    for ci in range(len(cset)):
        if ci in seen: continue
        norb+=1; stack=[ci]
        while stack:
            x=stack.pop()
            if x in seen: continue
            seen.add(x)
            # apply each h
            g0=Oh[idx[next(iter(cset[x]))]]
            for h in H:
                img=mul(h,g0)
                stack.append(rep_of[key(img)])
    return norb
a,b='C3v (tri-face, <111>)','C3v (tri-face, <111>)'
print(f"\ncross-check  |C3v\\O_h/C3v| via double cosets = {results[(a,b)][0]}, "
      f"via H-orbits on O_h/C3v = {orbits_H_on_GmodK(feat[a],feat[b])}")

print("\n\n=== junction symmetry of each double coset (residual = H ∩ gKg^{-1}) ===")
def residual(H,K,g):
    ginv=np.linalg.inv(g).round().astype(int)
    gKg=[ (g@k@ginv).round().astype(int) for k in K]
    gKkeys=set(key(m) for m in gKg)
    return sum(1 for h in H if key(h) in gKkeys)
def rep_type(g):
    d=int(round(np.linalg.det(g))); t=int(round(np.trace(g)))
    if d==1:
        return {3:'E (identity)',0:'C3',1:'C4',-1:'C2'}[t]
    else:
        return {-3:'i (inversion)',-1:'S4',0:'S6',1:'mirror σ'}[t]
for (a,b) in [('C3v (tri-face, <111>)','C3v (tri-face, <111>)'),
              ('C2v (vertex, <110>)','C2v (vertex, <110>)'),
              ('C4v (sq-face, <100>)','C4v (sq-face, <100>)')]:
    H,K=feat[a],feat[b]; reps,sizes=double_cosets(H,K)
    print(f"\n{a}  ⋈  {b}:")
    for g,s in sorted(zip(reps,sizes),key=lambda z:z[1]):
        print(f"   orbit {s:2d}   junction symmetry order {residual(H,K,g):2d}   rep ~ {rep_type(g)}")


# ======================================================================
# Does the mirrored grain actually leave the lattice?
# ======================================================================
# Counting orientations is not the same as counting NEW junctions.  A
# reoriented copy is a second grain only if it sits off the FCC lattice; if
# the mirror in a feature plane is already a lattice symmetry composed with a
# translation, it maps the crystal to itself and no boundary forms, so the
# "molecular" junction is really a coherent overlap under another name.
# This is the test that decides which family each feature contributes to.

NN_SHELL = [np.array(v) for v in itertools.product((-1,0,1), repeat=3)
            if sum(abs(x) for x in v) == 2]
SHELL = [np.array((0,0,0))] + NN_SHELL

def mirror_in_feature_plane(axis):
    """Reflect the shell in the plane of its outermost feature on `axis`."""
    a = np.array(axis, float); a /= np.linalg.norm(a)
    h = max(float(v @ a) for v in SHELL)
    return [v - 2*(float(v @ a) - h)*a for v in SHELL]

def on_lattice(pts):
    """Every node an FCC site: integer coordinates with even sum."""
    for v in pts:
        r = np.round(v)
        if not np.allclose(v, r, atol=1e-9) or int(round(sum(r))) % 2:
            return False
    return True

print("\n\n=== does the mirrored grain leave the lattice? ===")
print("(only a copy that leaves it is a genuine second grain)\n")
print(f"{'feature':22s}{'shared':8s}{'nodes':7s}{'on-lattice':12s}{'verdict'}")
Akeys = {tuple(int(round(x)) for x in v) for v in SHELL}
for label, axis in [('tri-face <111>', (1,1,1)),
                    ('sq-face <100>',  (1,0,0)),
                    ('vertex <110>',   (1,1,0))]:
    B = mirror_in_feature_plane(axis)
    lat = on_lattice(B)
    Bk = {tuple(int(round(x)) for x in v) for v in B} if lat else set()
    shared = (len(Akeys & Bk) if lat
              else sum(1 for p in SHELL for q in B
                       if np.allclose(p, q, atol=1e-9)))
    nodes = 26 - shared
    verdict = ("TRUE grain boundary" if not lat else
               "collapses to a coherent overlap")
    print(f"{label:22s}{shared:<8d}{nodes:<7d}{str(lat):12s}{verdict}")

print("\nThe {111} mirror is not an element of O_h, so it alone produces a")
print("copy off the lattice: three shared face atoms across 23 positioned")
print("nodes, the antiprismatic ABA twin.  The {100} and <110> mirrors are")
print("lattice symmetries, and their 'molecular' junctions are the <200> and")
print("<220> coherent overlaps of docking_enumerate.py.  So the lattice")
print("admits exactly ONE genuine two-shell grain boundary, and it is the")
print("{111} twin: the stacking fault that also carries strangeness and")
print("colour.  That is why the triangular face carries the physics.")
