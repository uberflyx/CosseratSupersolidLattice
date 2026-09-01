"""Exact two-orbit radial blocks for the accommodated baryons, checked against the
numerical Cosserat matrices.

The radial microrotation ("hedgehog") phi_i ~ r_i with u = 0 is curl-free, so it is an
exact eigenvector of the coupling and its Rayleigh quotient is a bond count:

    H_shell = 1 + (24 + n_out)/12 = 3 + n_out/12,
    H_off   = -(1/sqrt(12 m)) sum_<ij> rhat_i . rhat_j,

with n_out the bonds leaving the twelve cuboctahedral vertices and m the extension's
node count.  Only two cosines occur in the family: sqrt(2/3) on a void bond and
sqrt(3)/2 on a cap bond, so every block entry is exact.

The Sigma row is a closure (the hedgehog lands on its root with weight 0.90).  The cap
hosts scatter the hedgehog (weights 0.17 and 0.18), so their rows are reported for
completeness and are NOT closures; see the monograph's table caption.  The Born row is
the check: its integer block reproduces the independently derived A_1u pair.
"""
import sys, numpy as np, sympy as sp
sys.path.insert(0,'/home/claude/CosseratSupersolidLattice/spectral_mass')
from delta_first_principles import cluster_delta, build_cosserat_matrix_two_d as build
from hyperons_first_principles import cluster_octet
M_E=sp.Rational(51099895,100000000); M_0=M_E*sp.Rational(1370359992,10000000)
def spines(host,idx):
    nh=len(host); v=np.zeros(6*nh)
    for i in idx:
        r=host[i]; n=np.linalg.norm(r)
        if n>1e-9: v[3*nh+3*i:3*nh+3*i+3]=r/n
    return v/np.linalg.norm(v)
shell=cluster_delta()[:13]; delta=cluster_delta()
born=np.vstack([shell,np.sqrt(2)*np.vstack([np.eye(3),-np.eye(3)])])
cases=[(delta,17,"Sigma",[5,-sp.sqrt(2),9],1193.15),
       (cluster_octet(1)[0],16,"Lambda",[sp.Rational(9,2),-sp.sqrt(3)/2,sp.Rational(10,3)],1115.683),
       (cluster_octet(2)[0],19,"Xi",[5,-sp.sqrt(6)/2,sp.Rational(10,3)],1318.29),
       (born,19,"Born",[6,-2,5],None)]
for host,N,lab,(a,b,c),pdg in cases:
    nh=len(host); M=build(host); s1=spines(host,range(13)); s2=spines(host,range(13,nh))
    num=np.array([[s1@M@s1,s1@M@s2],[s2@M@s1,s2@M@s2]])
    exact=sp.Matrix([[a,b],[b,c]]); err=abs(num-np.array(exact.evalf(),dtype=float)).max()
    root=min(exact.eigenvals(),key=lambda e: sp.re(sp.N(e))); root=sp.nsimplify(sp.simplify(root))
    m=sp.N(N*M_0-N*(4-root)*M_E,9)
    r=f"  PDG {pdg}  {float(100*(m-pdg)/pdg):+.3f}%" if pdg else ""
    print(f"{lab:7s} exact block matches numeric to {err:.2e};  root {root} = {sp.N(root,8)};  m = {m}{r}")
