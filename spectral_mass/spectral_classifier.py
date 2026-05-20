"""
spectral_classifier.py
======================

Companion script for the chapter "Spectral classification of the building
blocks" in the FCC defect-graph monograph.  Reproduces every numerical entry
of that chapter from FCC geometry alone.

Method
------
1.  For each cluster (hex cap, cuboctahedron, coordination shell, Born
    cluster), build the FCC-embedded vertex coordinates.
2.  Build the 3n x 3n dynamical matrix with central nearest-neighbour springs
    of unit stiffness, plus second-NN springs for the Born cluster.
3.  Identify the subgroup of O_h that maps the cluster to itself.
4.  Decompose each dynamical-matrix eigenspace into irreps of that subgroup
    (O_h for the O_h-symmetric clusters, D_3d for the hex cap) using the
    standard character formula.
5.  Print the spectrum: eigenvalue, multiplicity, irrep content.

Authors: M. Cox, with Claude (Anthropic).
License: MIT.  Reference: monograph Sec. sec:cluster_spectral.
"""

import numpy as np
from scipy.linalg import eigh
from itertools import product


ELL = 1.0
A_LAT = np.sqrt(2.0) * ELL


# ---------------------------------------------------------------- clusters

def fcc_nn_vectors():
    a = A_LAT
    nns = []
    for sx, sy in [(1, 1), (-1, -1), (1, -1), (-1, 1)]:
        nns.append((sx * a / 2, sy * a / 2, 0))
        nns.append((sx * a / 2, 0, sy * a / 2))
        nns.append((0, sx * a / 2, sy * a / 2))
    return np.array(nns)


def cluster_hex_cap():
    e1 = np.array([1., -1., 0.]) / np.sqrt(2)
    e2 = np.array([1., 1., -2.]) / np.sqrt(6)
    coords = [np.array([0., 0., 0.])]
    for i in range(6):
        a = i * np.pi / 3
        coords.append(ELL * (np.cos(a) * e1 + np.sin(a) * e2))
    return np.array(coords), "Hex cap (7 nodes)"


def cluster_cuboctahedron():
    return fcc_nn_vectors(), "Cuboctahedron (12 nodes)"


def cluster_coord_shell():
    return np.vstack([np.array([[0., 0., 0.]]), fcc_nn_vectors()]), \
           "Coordination shell (13 nodes)"


# Lambda baryon geometry: the four inactive {111} directions, the hex-cap
# extension that builds on one of them, and the Lambda ground-state cluster.
INACTIVE_DIRS = [
    (-1, -1, -1),
    (+1, +1, -1),
    (+1, -1, +1),
    (-1, +1, +1),
]


def hex_cap_extension_on_inactive_dir(direction):
    """Three FCC sites on the {111} plane projecting at -2 a/sqrt(3) along
    the inactive direction."""
    a_fcc = np.sqrt(2.0)
    sublat = [np.array([0, 0, 0]), np.array([a_fcc/2, a_fcc/2, 0]),
              np.array([a_fcc/2, 0, a_fcc/2]), np.array([0, a_fcc/2, a_fcc/2])]
    fcc_sites = set()
    for ni, nj, nk in product(range(-3, 4), repeat=3):
        for s in sublat:
            p = np.array([ni*a_fcc, nj*a_fcc, nk*a_fcc]) + s
            fcc_sites.add(tuple(np.round(p, 8)))
    fcc_arr = np.array(sorted(fcc_sites,
                              key=lambda p: (p[0]**2 + p[1]**2 + p[2]**2,
                                             p[0], p[1], p[2])))
    d = np.array(direction, dtype=float) / np.sqrt(3)
    target_proj = -2.0 * a_fcc / np.sqrt(3)
    projs = fcc_arr @ d
    on_plane = fcc_arr[np.abs(projs - target_proj) < 1e-6]
    perp = np.linalg.norm(on_plane - np.outer(on_plane @ d, d), axis=1)
    return on_plane[np.argsort(perp)[:3]]


def build_lambda_cluster():
    """Lambda = proton coordination shell + one hex-cap extension (no voids).

    The octet baryon's defect is the bare screw dislocation (no T_d voids
    needed for spin-1/2). Strange flavour adds one hex extension on an
    inactive (-1,-1,-1) {111} face. Total: 13 + 3 = 16 atoms.
    """
    shell, _ = cluster_coord_shell()
    ext = hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[0])
    cluster = np.vstack([shell, ext])
    unique = []
    for p in cluster:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in unique):
            unique.append(p)
    return np.array(unique)


def cluster_born():
    centre = np.array([[0., 0., 0.]])
    co = fcc_nn_vectors()
    second = []
    for s in [+1, -1]:
        for ax in range(3):
            v = np.zeros(3)
            v[ax] = s * A_LAT
            second.append(v)
    return np.vstack([centre, co, np.array(second)]), "Born cluster (19 nodes)"


# -------------------------------------------------------- dynamical matrix

def build_dynamical_matrix(coords, K_nn=1.0, K_2nn=0.0, ell=ELL, tol=1e-6):
    """3n x 3n Hessian for central-force springs.

    Phi[i*3+a, j*3+b] = -K rhat_ij^a rhat_ij^b   for bonded i != j
    Phi[i*3+a, i*3+b] = K sum_{j bonded} rhat_ij^a rhat_ij^b
    """
    n = len(coords)
    Phi = np.zeros((3 * n, 3 * n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r = coords[j] - coords[i]
            d = np.linalg.norm(r)
            K = None
            if abs(d - ell) < tol:
                K = K_nn
            elif abs(d - ell * np.sqrt(2)) < tol:
                K = K_2nn
            if K is not None and K > 0:
                rhat = r / d
                outer = np.outer(rhat, rhat)
                Phi[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= K * outer
                Phi[3 * i:3 * i + 3, 3 * i:3 * i + 3] += K * outer
    return Phi


# ---------------------------------------------------------------- O_h group

def generate_Oh():
    elements = []
    seen = set()
    for perm in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        for signs in product([1, -1], repeat=3):
            R = np.zeros((3, 3))
            for i, p in enumerate(perm):
                R[i, p] = signs[i]
            key = tuple(R.flatten())
            if key not in seen:
                seen.add(key)
                elements.append(R)
    assert len(elements) == 48
    return elements


OH_CLASSES = ['E', 'C_3', 'C_2(diag)', 'C_4', 'C_2(cube)',
              'i', 'S_6', 'sigma_d', 'S_4', 'sigma_h']
OH_CLASS_SIZES = {'E': 1, 'C_3': 8, 'C_2(diag)': 6, 'C_4': 6, 'C_2(cube)': 3,
                  'i': 1, 'S_6': 8, 'sigma_d': 6, 'S_4': 6, 'sigma_h': 3}
OH_CHARACTERS = {
    'A_1g': {'E':1,'C_3':1,'C_2(diag)':1,'C_4':1,'C_2(cube)':1,'i':1,'S_6':1,'sigma_d':1,'S_4':1,'sigma_h':1},
    'A_2g': {'E':1,'C_3':1,'C_2(diag)':-1,'C_4':-1,'C_2(cube)':1,'i':1,'S_6':1,'sigma_d':-1,'S_4':-1,'sigma_h':1},
    'E_g':  {'E':2,'C_3':-1,'C_2(diag)':0,'C_4':0,'C_2(cube)':2,'i':2,'S_6':-1,'sigma_d':0,'S_4':0,'sigma_h':2},
    'T_1g': {'E':3,'C_3':0,'C_2(diag)':-1,'C_4':1,'C_2(cube)':-1,'i':3,'S_6':0,'sigma_d':-1,'S_4':1,'sigma_h':-1},
    'T_2g': {'E':3,'C_3':0,'C_2(diag)':1,'C_4':-1,'C_2(cube)':-1,'i':3,'S_6':0,'sigma_d':1,'S_4':-1,'sigma_h':-1},
    'A_1u': {'E':1,'C_3':1,'C_2(diag)':1,'C_4':1,'C_2(cube)':1,'i':-1,'S_6':-1,'sigma_d':-1,'S_4':-1,'sigma_h':-1},
    'A_2u': {'E':1,'C_3':1,'C_2(diag)':-1,'C_4':-1,'C_2(cube)':1,'i':-1,'S_6':-1,'sigma_d':1,'S_4':1,'sigma_h':-1},
    'E_u':  {'E':2,'C_3':-1,'C_2(diag)':0,'C_4':0,'C_2(cube)':2,'i':-2,'S_6':1,'sigma_d':0,'S_4':0,'sigma_h':-2},
    'T_1u': {'E':3,'C_3':0,'C_2(diag)':-1,'C_4':1,'C_2(cube)':-1,'i':-3,'S_6':0,'sigma_d':1,'S_4':-1,'sigma_h':1},
    'T_2u': {'E':3,'C_3':0,'C_2(diag)':1,'C_4':-1,'C_2(cube)':-1,'i':-3,'S_6':0,'sigma_d':-1,'S_4':1,'sigma_h':1},
}

D3D_CLASSES = ['E', 'C_3', 'C_2(diag)', 'i', 'S_6', 'sigma_d']
D3D_CLASS_SIZES = {'E': 1, 'C_3': 2, 'C_2(diag)': 3, 'i': 1, 'S_6': 2, 'sigma_d': 3}
D3D_CHARACTERS = {
    'A_1g': {'E':1,'C_3':1,'C_2(diag)':1,'i':1,'S_6':1,'sigma_d':1},
    'A_2g': {'E':1,'C_3':1,'C_2(diag)':-1,'i':1,'S_6':1,'sigma_d':-1},
    'E_g':  {'E':2,'C_3':-1,'C_2(diag)':0,'i':2,'S_6':-1,'sigma_d':0},
    'A_1u': {'E':1,'C_3':1,'C_2(diag)':1,'i':-1,'S_6':-1,'sigma_d':-1},
    'A_2u': {'E':1,'C_3':1,'C_2(diag)':-1,'i':-1,'S_6':-1,'sigma_d':1},
    'E_u':  {'E':2,'C_3':-1,'C_2(diag)':0,'i':-2,'S_6':1,'sigma_d':0},
}


def classify_operation(R, tol=1e-4):
    det = np.linalg.det(R)
    trace = np.trace(R)

    def axis_180(R):
        ev, vc = np.linalg.eig(R)
        for k in range(3):
            if abs(ev[k].real - 1) < 0.1 and abs(ev[k].imag) < 0.1:
                a = vc[:, k].real
                return a / np.linalg.norm(a)
        return None

    def axis_mirror(R):
        ev, vc = np.linalg.eig(R)
        for k in range(3):
            if abs(ev[k].real + 1) < 0.1 and abs(ev[k].imag) < 0.1:
                a = vc[:, k].real
                return a / np.linalg.norm(a)
        return None

    def is_cube_axis(a):
        s = np.sort(np.abs(a))[::-1]
        return s[0] > 0.99 and s[1] < 0.1

    if det > 0:
        if abs(trace - 3) < tol: return 'E'
        if abs(trace - 0) < tol: return 'C_3'
        if abs(trace - 1) < tol: return 'C_4'
        if abs(trace + 1) < tol:
            a = axis_180(R)
            return 'C_2(cube)' if a is not None and is_cube_axis(a) else 'C_2(diag)'
    else:
        if abs(trace + 3) < tol: return 'i'
        if abs(trace + 0) < tol: return 'S_6'
        if abs(trace + 1) < tol: return 'S_4'
        if abs(trace - 1) < tol:
            a = axis_mirror(R)
            return 'sigma_h' if a is not None and is_cube_axis(a) else 'sigma_d'
    return 'unknown'


# -------------------- displacement representation and eigenspace decomp

def vertex_perm(R, coords, tol=1e-6):
    n = len(coords)
    perm = []
    for i in range(n):
        rot = R @ coords[i]
        m = None
        for j in range(n):
            if np.linalg.norm(rot - coords[j]) < tol:
                m = j
                break
        if m is None:
            return None
        perm.append(m)
    return perm if len(set(perm)) == n else None


def disp_rep(R, coords, tol=1e-6):
    perm = vertex_perm(R, coords, tol)
    if perm is None:
        return None
    n = len(coords)
    M = np.zeros((3 * n, 3 * n))
    for i in range(n):
        j = perm[i]
        M[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
    return M


def cluster_subgroup(coords, oh_elements):
    return [R for R in oh_elements if vertex_perm(R, coords) is not None]


def group_eigvalues(eigvals, tol=1e-4):
    idx_sorted = np.argsort(eigvals)
    groups, cur, cur_i = [], None, []
    for idx in idx_sorted:
        v = eigvals[idx]
        if cur is None or abs(v - cur) > tol:
            if cur_i:
                groups.append((cur, cur_i))
            cur = v
            cur_i = [idx]
        else:
            cur_i.append(idx)
    if cur_i:
        groups.append((cur, cur_i))
    return groups


def decompose(eigvecs, indices, coords, subgroup, characters, classes, sizes):
    V = eigvecs[:, indices]
    chi = {c: 0.0 for c in classes}
    for R in subgroup:
        M = disp_rep(R, coords)
        if M is None:
            continue
        chi_g = np.trace(V.T @ M @ V)
        cls = classify_operation(R)
        if cls in chi:
            chi[cls] += chi_g
    H = sum(sizes.values())
    out = {}
    for ir, ch in characters.items():
        n = sum(ch[c] * chi[c] for c in classes) / H
        if abs(n) > 0.01:
            out[ir] = round(n)
    return out


# ---------------------------------------------------------------- driver

def report(builder, K_2nn, oh):
    coords, name = builder()
    sub = cluster_subgroup(coords, oh)
    H = len(sub)
    Phi = build_dynamical_matrix(coords, K_nn=1.0, K_2nn=K_2nn)
    ev, ec = eigh(Phi)
    groups = group_eigvalues(ev)

    if H == 48:
        chars, classes, sizes, sym = OH_CHARACTERS, OH_CLASSES, OH_CLASS_SIZES, 'O_h'
    elif H == 12:
        chars, classes, sizes, sym = D3D_CHARACTERS, D3D_CLASSES, D3D_CLASS_SIZES, 'D_3d'
    else:
        chars = None
        sym = f'|H|={H}'

    print(f'\n{name}: n={len(coords)}, |H|={H}, symmetry={sym}')
    print(f'  {"lambda/K":>10s}  {"mult":>4s}  irrep content')
    print(f'  {"-"*10}  {"-"*4}  {"-"*55}')
    for val, idx in groups:
        if chars is not None:
            d = decompose(ec, idx, coords, sub, chars, classes, sizes)
            dstr = ' + '.join([f'{n}*{ir}' if n > 1 else ir for ir, n in d.items()])
        else:
            dstr = f'(low-symmetry, |H|={H}; not classified)'
        marker = ' [ZERO]' if abs(val) < 1e-5 else ''
        print(f'  {val:8.4f}{marker:7s}  {len(idx):4d}  {dstr}')


def main():
    print('=' * 78)
    print('Spectral classification of FCC building blocks')
    print('Dynamical-matrix decomposition under cluster point group')
    print('=' * 78)

    oh = generate_Oh()
    print(f'\nGenerated O_h group: |O_h| = {len(oh)} elements')

    for builder, K_2 in [
        (cluster_cuboctahedron, 0.0),
        (cluster_coord_shell, 0.0),
        (cluster_born, 1.0),
        (cluster_hex_cap, 0.0),
    ]:
        report(builder, K_2, oh)

    print('\nReproduces Tables ref:cubocta_dynmat, ref:shell_dynmat,')
    print('and ref:hex_cap_dynmat of chapter sec:cluster_spectral.')


if __name__ == '__main__':
    main()
