import numpy as np, itertools
np.set_printoptions(suppress=True)

# ---- FCC cuboctahedron coordination shell (NN distance = sqrt2) ----
NN = [np.array(v) for v in set(itertools.permutations((1,1,0)))
      | set(itertools.permutations((1,-1,0)))
      | set(itertools.permutations((-1,-1,0)))]
NN = [v for v in NN if np.dot(v,v)==2]
assert len(NN)==12
shellA = [np.array((0,0,0))] + NN          # 13 atoms: centre + 12 vertices
def key(v): return tuple(int(x) for x in v)
Aset = set(key(a) for a in shellA)

# FCC sites = integer points with even coordinate sum
def fcc_sites(rmax2):
    out=[]
    R=int(np.ceil(np.sqrt(rmax2)))+1
    for x in range(-R,R+1):
        for y in range(-R,R+1):
            for z in range(-R,R+1):
                if (x+y+z)%2==0 and 0<x*x+y*y+z*z<=rmax2:
                    out.append((x,y,z))
    return out

# O_h orbit rep of a vector (to group equivalent separations)
def oh_orbit_key(v):
    vs=[]
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1,-1),repeat=3):
            vs.append(tuple(s[i]*v[p[i]] for i in range(3)))
    return min(vs)

# classify the shape of a shared-atom set
def shape(pts):
    n=len(pts)
    if n==0: return "none"
    if n==1: return "1 node (vertex)"
    P=np.array(pts)
    # pairwise squared distances
    d2=sorted(set(int(round(np.dot(P[i]-P[j],P[i]-P[j]))) for i in range(n) for j in range(n) if i<j))
    if n==3 and d2==[2]: return "3 atoms: triangular face (all NN)"
    if n==4 and d2==[2,4]: return "4 atoms: square face"
    if n==2 and d2==[2]: return "2 atoms (edge)"
    return f"{n} atoms (dists^2 {d2})"

print("LATTICE-COHERENT dockings  (cluster B = FCC translate of A; two grains of one crystal)\n")
print(f"{'separation d':16s}{'|d|^2':6s}{'shared':8s}{'shape of shared set':34s}{'named mode'}")
seen=set()
known={1:'Mode V (dibaryon d*)',4:'Mode P (deuteron)'}
for d in fcc_sites(8):           # overlap possible only if |d| <= 2*sqrt2 -> |d|^2<=8
    ok=oh_orbit_key(d)
    if ok in seen: continue
    seen.add(ok)
    dv=np.array(d)
    Bset=set(key(a+dv) for a in shellA)
    shared=[np.array(s) for s in (Aset & Bset)]
    sh=shape([tuple(s) for s in shared])
    mode=known.get(len(shared),'')
    # midpoint atom present? (the V bridge)
    mid = key(dv/2) if all(c%2==0 for c in d) else None
    note = '  (shares midpoint bridge)' if (mid and mid in (Aset&Bset)) else ''
    print(f"{str(d):16s}{int(np.dot(d,d)):<6d}{len(shared):<8d}{sh:34s}{mode}{note}")

print("\nMOLECULAR dockings  (cluster B reoriented, face-to-face; two separate grains bonded at interface)")
print("  = the double cosets H \\ O_h / K of the feature stabilisers (computed separately):")
print("    triangular-face <-> triangular-face : 4 orientations, 2 lattice-coherent (eclipsed=FCC, antiprism=twin)  -> Mode F / Mode 2F family")
print("    square-face     <-> square-face     : 3 orientations, 2 coherent")
print("    vertex          <-> vertex          : 5 orientations, 2 coherent")

# count the distinct lattice-coherent translational interfaces (by shared shape, excluding pure point-touch)
print("\nDISTINCT shared-interface types found among lattice-coherent overlaps:")
seen=set(); types={}
for d in fcc_sites(8):
    dv=np.array(d); Bset=set(key(a+dv) for a in shellA)
    shared=tuple(sorted(Aset & Bset))
    s=shape([tuple(x) for x in shared])
    types[s]=types.get(s,0)
for s in sorted(types): print("   -", s)
