"""
Larger nuclei in the framework = close-packed alpha clusters.

The A=4 wall (FCC has no 5 mutually-touching sites) makes the alpha the densest
single compound. Beyond A=4 the nucleus is a packing of alphas. Alphas, like any
hard units, pack densest in the FCC arrangement, so a compact cluster of n alphas
has a definite number p(n) of nearest-neighbour alpha-alpha contacts, fixed by
geometry. If each contact is one inter-shell dock of energy B_bond, the nuclear
binding is
        B(n alpha) = n * B_alpha + p(n) * B_bond,
the Hafstad-Teller alpha-model law (1938), here with p from FCC geometry and
B_bond a single dock rather than a fit.

We (1) compute p(n) from FCC packing, (2) compare to the bond counts the data
demands, (3) predict the binding curve, with an alpha-alpha Coulomb term.
"""
import numpy as np
from itertools import combinations

# ---------- FCC packing geometry: p(n) = max alpha-alpha contacts ----------
def fcc_sites(R):
    out, r = [], int(np.ceil(R))
    for x in range(-r, r+1):
        for y in range(-r, r+1):
            for z in range(-r, r+1):
                if (x+y+z) % 2 == 0 and x*x+y*y+z*z <= R*R+1e-9:
                    out.append((x, y, z))
    return np.array(out)

NN = {(a,b,0) for a in(-1,1) for b in(-1,1)} | \
     {(a,0,b) for a in(-1,1) for b in(-1,1)} | \
     {(0,a,b) for a in(-1,1) for b in(-1,1)}

def bonds(subset, idx):
    """count nearest-neighbour contacts inside a chosen subset of FCC sites"""
    s = set(subset)
    c = 0
    for site in subset:
        for d in NN:
            nb = (site[0]+d[0], site[1]+d[1], site[2]+d[2])
            if nb in s and nb > site:   # count each pair once
                c += 1
    return c

def compact_cluster(n, pool):
    """greedy: grow an n-site cluster maximising internal contacts"""
    pool = [tuple(p) for p in pool]
    # seed at the most central site
    cen = min(pool, key=lambda p: p[0]**2+p[1]**2+p[2]**2)
    chosen = [cen]
    while len(chosen) < n:
        best, best_gain = None, -1
        cset = set(chosen)
        for cand in pool:
            if cand in cset:
                continue
            gain = sum(1 for d in NN
                       if (cand[0]+d[0],cand[1]+d[1],cand[2]+d[2]) in cset)
            # tie-break toward the centroid to stay compact
            if gain > best_gain or (gain == best_gain and best is not None and
               (cand[0]**2+cand[1]**2+cand[2]**2) < (best[0]**2+best[1]**2+best[2]**2)):
                best, best_gain = cand, gain
        chosen.append(best)
    return chosen

pool = fcc_sites(4.0)
p_fcc = {}
for n in range(1, 11):
    cl = compact_cluster(n, pool)
    p_fcc[n] = bonds(cl, None)

# ---------- the data: alpha-conjugate nuclei (N=Z), AME2020 binding energies ----------
# B in MeV (to be verified against AME2020); Z = 2n
nuclei = {
    1:  ("He-4",  28.296),
    2:  ("Be-8",  56.500),
    3:  ("C-12",  92.162),
    4:  ("O-16",  127.619),
    5:  ("Ne-20", 160.645),
    6:  ("Mg-24", 198.257),
    7:  ("Si-28", 236.537),
    8:  ("S-32",  271.781),
    9:  ("Ar-36", 306.716),
    10: ("Ca-40", 342.052),
}

B_alpha = 28.296          # MeV, the alpha's own binding (= n=1 entry)
# bond energy the data demands, pair by pair (subtract alpha self-binding, divide by p)
print("n   nucleus   p_FCC   B_meas    n*B_a    (B-n*B_a)/p = B_bond")
for n,(name,B) in nuclei.items():
    p = p_fcc[n]
    extra = B - n*B_alpha
    bb = extra/p if p>0 else float('nan')
    print(f"{n:>2}  {name:<7}  {p:>4}   {B:>7.2f}   {n*B_alpha:>7.2f}    {bb:>6.3f}")
