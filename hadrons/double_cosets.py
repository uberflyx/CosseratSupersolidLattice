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
