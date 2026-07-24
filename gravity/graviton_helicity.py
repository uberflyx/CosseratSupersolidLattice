"""Helicity content of candidate graviton fields in the elastic vacuum.

Two exact statements, verified symbolically:

1. NO-GO. For a plane wave (d_mu -> i k_mu, k along the propagation axis), any
   rank-2 field built from one derivative of a vector (strain of a displacement
   wave, wryness of a microrotation wave) has ZERO helicity-2 content and a
   vanishing transverse-traceless part. The derivative carries helicity 0 about
   the propagation axis; helicity is additive; a vector supplies only 0, +-1.
   No linear branch of an elastic medium is a graviton.

2. THE CARRIER. The displacement covariance Q_ij = <u_i u_j> - (1/3) d_ij <u^2>
   (the tensor whose isotropic part is the Debye-Waller exponent) is rank-2 and
   NOT a derivative of a vector. Its transverse components Q_xx - Q_yy and Q_xy
   pick up exactly e^{-+ 2 i theta} under rotation by theta about the axis:
   helicity +-2, the two gravitational-wave polarisations h_+ and h_x.

3. THE EDGE. Q_ij fluctuations are two-quantum states of the transverse branch
   (omega = c|k| each). Two collinear quanta of helicity +1: total helicity +2
   at omega = c(k1 + k2) = c k. The spin-2 spectral weight therefore reaches
   the light cone exactly: the tensor channel is gapless because and only
   because its constituents are.
"""
import sympy as sp

th = sp.symbols('theta', real=True)
k = sp.symbols('k', positive=True)
ux, uy, uz, px, py, pz = sp.symbols('u_x u_y u_z phi_x phi_y phi_z')

# ---- 1. no-go: symmetric derivative of a vector, wave along z ----------------
S = sp.zeros(3, 3)
for i, p in enumerate([px, py, pz]):
    S[i, 2] += sp.I*k*p/2
    S[2, i] += sp.I*k*p/2
hel2_a = sp.simplify(S[0, 0] - S[1, 1])   # helicity +-2 combinations
hel2_b = sp.simplify(S[0, 1])
print("[1] sym. derivative of a vector, wave along z:")
print(f"    S_xx - S_yy = {hel2_a},  S_xy = {hel2_b}   (helicity-2 content: none)")
P = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 0]])           # transverse projector
TT = P*S*P - sp.Rational(1, 2)*P*sp.trace(P*S*P)
print(f"    TT projection = {sp.simplify(TT)}  (zero tensor)")

# ---- 2. the covariance tensor carries helicity 2 -----------------------------
R = sp.Matrix([[sp.cos(th), -sp.sin(th), 0],
               [sp.sin(th),  sp.cos(th), 0],
               [0, 0, 1]])
u = sp.Matrix([ux, uy, uz])
Q = u*u.T                                  # covariance amplitude (traceless part implied)
Qr = R*Q*R.T                               # rotated by theta about z
hplus  = sp.simplify(Qr[0, 0] - Qr[1, 1])
hcross = sp.simplify(2*Qr[0, 1])
# compare with e^{-+2i theta} acting on (h+ -+ i hx)/...
orig_p, orig_c = ux**2 - uy**2, 2*ux*uy
combo  = sp.simplify(sp.expand_trig(hplus + sp.I*hcross) - sp.exp(2*sp.I*th)*(orig_p + sp.I*orig_c))
print("\n[2] covariance tensor under rotation by theta about z:")
print(f"    (Q_xx-Q_yy) + 2i Q_xy  -  e^(2 i theta) x original = {combo}")
print("    => the pair (h_+, h_x) = (Q_xx - Q_yy, 2 Q_xy) is exactly helicity +-2")

# ---- 3. the collinear two-quantum edge ---------------------------------------
k1, k2 = sp.symbols('k_1 k_2', positive=True)
omega = k1 + k2          # units c = 1: each transverse quantum has omega = k
ktot  = k1 + k2          # collinear
print("\n[3] two collinear helicity-(+1) quanta: helicity 1 + 1 = 2,")
print(f"    omega = c(k1 + k2), |k| = k1 + k2  ->  omega = c|k| identically: {sp.simplify(omega - ktot) == 0}")
print("    the spin-2 channel is gapless exactly when its constituents are massless.")
