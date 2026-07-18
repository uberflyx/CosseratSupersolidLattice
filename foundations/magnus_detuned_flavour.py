"""Detuned two-channel Magnus problem, subspace-projected.

Channel a (glide leg): width 1, edge 1. Channel b (stacking leg): width
sqrt(2), edge 1/sqrt(2) = 0.7071 (Pythagoras on the interlayer bond).
Time-harmonic coupled system, psi ~ e^{-i om t}:
    (H_a - om^2) psi_a - i gamma om B psi_b = 0
    (H_b - om^2) psi_b + i gamma om B psi_a = 0
For fixed om this is the Hermitian pencil A(om) psi = om^2 psi with
    A(om) = [[H_a, -i gamma om B], [ i gamma om B, H_b]],
solved self-consistently. Projection onto the lowest M modes of each
channel keeps every iteration cheap; convergence in M is checked.
"""
import numpy as np

N, L, M = 768, 150.0, 320
x = (np.arange(N)-N/2)*(2*L/N); dx = x[1]-x[0]
k = np.fft.fftfreq(N, d=dx)*2*np.pi
F = np.fft.fft(np.eye(N), axis=0)
Lam = (np.fft.ifft(np.abs(k)[:,None]*F, axis=0)).real
Lam = 0.5*(Lam+Lam.T)

def Hw(w):
    return Lam + np.diag((1.0/w)*(x**2-w**2)/(x**2+w**2))

wa, wb = 1.0, np.sqrt(2.0)
ea, va = np.linalg.eigh(Hw(wa))
eb, vb = np.linalg.eigh(Hw(wb))
th_a, th_b = va[:,1], vb[:,1]          # odd threshold (width) modes
print("edges: a=1.0000, b=%.4f | threshold E: a=%.5f, b=%.5f"
      % (1/wb, ea[1], eb[1]))

Ua, Ub = va[:,:M], vb[:,:M]            # low subspaces
Da, Db = np.diag(ea[:M]), np.diag(eb[:M])
Bfull = np.diag(1.0/(x**2+1.0))        # circulation profile, healing ~ w_par
Bab = Ua.T @ Bfull @ Ub                # M x M cross block

tha = np.zeros(2*M); tha[:M] = Ua.T @ th_a          # threshold-a in basis
thb = np.zeros(2*M); thb[M:] = Ub.T @ th_b          # threshold-b in basis

def solve(gamma, E0, target, iters=80):
    """Track the mode of max overlap with `target`, iterating om^2 = lam."""
    E = E0
    for _ in range(iters):
        om = np.sqrt(max(E, 1e-12))
        A = np.block([[Da, -1j*gamma*om*Bab],
                      [1j*gamma*om*Bab.T, Db]])
        lam, U = np.linalg.eigh(A)
        ov = np.abs(U.conj().T @ target)**2
        j = np.argmax(ov)
        En = lam[j]
        if abs(En - E) < 1e-11:
            return En, ov[j]
        E = 0.5*(E + En)
    return E, ov[j]

print("\ntracking the LOWER threshold (stacking, edge 0.7071):")
print(" gamma   E_track   below b-edge?  ov")
for g in [0.0, 0.2, 0.4, 0.68, 0.81, 1.2]:
    E, ov = solve(g, 0.5, thb + 0j)
    print(f" {g:4.2f}   {E:+.4f}   {'YES' if E < 1/wb - 1e-4 else 'no '}"
          f"          {ov:.3f}")

print("\ntracking the UPPER threshold (glide, edge 1.0):")
print(" gamma   E_track   below b-edge?  ov")
for g in [0.0, 0.2, 0.4, 0.68, 0.81, 1.2]:
    E, ov = solve(g, 1.0, tha + 0j)
    print(f" {g:4.2f}   {E:+.4f}   {'YES' if E < 1/wb - 1e-4 else 'no '}"
          f"          {ov:.3f}")

# convergence check in M at one working point
for Mtest in [200, 320]:
    Ua_, Ub_ = va[:,:Mtest], vb[:,:Mtest]
    Da_, Db_ = np.diag(ea[:Mtest]), np.diag(eb[:Mtest])
    Bab_ = Ua_.T @ Bfull @ Ub_
    t_ = np.zeros(2*Mtest); t_[Mtest:] = Ub_.T @ th_b
    E = 0.5
    for _ in range(80):
        om = np.sqrt(max(E,1e-12))
        A = np.block([[Da_, -1j*0.68*om*Bab_],[1j*0.68*om*Bab_.T, Db_]])
        lam, U = np.linalg.eigh(A)
        j = np.argmax(np.abs(U.conj().T @ (t_+0j))**2)
        En = lam[j]
        if abs(En-E) < 1e-11: break
        E = 0.5*(E+En)
    print(f"M={Mtest}: E_bound(gamma=0.68) = {E:+.5f}")
