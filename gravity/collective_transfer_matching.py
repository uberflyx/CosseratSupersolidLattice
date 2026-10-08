#!/usr/bin/env python3
"""Check candidate collective transfer operators and their source normalisation.

This is an algebraic diagnostic, not a derivation of Newton's constant.
It distinguishes a coherent linear field from independent product amplitudes,
and checks which intermediate channel kernel gives a single rank-one trace.
It also tests metric inertia and source overlap on the actual FCC cluster.
Finally it shows that the endpoint identity is the inertial (f-sum) limit of
the channel contraction, whose static limit is one, and that requiring
wavelength-independent static links selects the vertex's slide-spin covector.
The electromagnetic amplitude is an input; no gravitational datum is fitted.
"""

import itertools
import json
from collections import deque
import numpy as np


def fcc_cluster_graph():
    """Centre, twelve nearest neighbours, six axial next-nearest neighbours."""
    sites = [(0, 0, 0)]
    sites += [p for p in itertools.product((-1, 0, 1), repeat=3)
              if sum(x*x for x in p) == 2]
    sites += [p for p in itertools.product((-2, 0, 2), repeat=3)
              if sum(x*x for x in p) == 4]
    sites = np.array(sites)
    edges = [(i, j) for i, j in itertools.combinations(range(len(sites)), 2)
             if np.sum((sites[i]-sites[j])**2) == 2]
    return sites, edges


def microscopic_normalisation_checks(alpha):
    """Necessary conditions in two restricted models, not a derived G.

    First grant one physical TT polarisation, [h,p]=i, and the Hamiltonian
    U/2 sum p**2 + K/2 sum_bonds (h_i-h_j)**2, with each bond counted once.
    h is physical metric strain, with polarisation norm e_ij e_ij=2.
    Required coefficients use the proposed alpha_G only as a target.

    Next give nineteen continuum fields identical kinetic coefficients and
    positive nearest-neighbour mixing (h_i-q_ij*h_j)**2. Test the smallest
    possible local source residue of an exactly massless mode, given a
    bound on each link ratio. This is not a tunnelling-amplitude bound.
    """
    prefactor = (1+1/np.pi)*(1-17*alpha/18)
    target = prefactor*alpha**19
    sites, edges = fcc_cluster_graph()
    directions = sites[1:13]/np.sqrt(2)
    moment = directions.T @ directions
    assert np.allclose(moment, 4*np.eye(3), atol=1e-14)

    # Independent discrete dispersion checks the continuum spatial factor.
    # a=1 here; eigenvalue of the FCC graph Laplacian is sum(1-cos(k.delta)).
    dispersion = []
    for direction in ((0, 0, 1), (1, 1, 0), (1, 1, 1)):
        direction = np.array(direction, dtype=float)
        direction /= np.linalg.norm(direction)
        ka = 1e-3
        laplacian = np.sum(2*np.sin(ka*(directions @ direction)/2)**2)
        ratio = float(laplacian/(ka*ka))
        assert abs(ratio-2) < ka*ka
        dispersion.append(ratio)

    u_required = 16*np.sqrt(2)*np.pi*target
    k_required = .5/u_required
    # Node units: hbar=c=a=E0=1; primitive volume=1/sqrt(2).
    volume = 1/np.sqrt(2)
    inertia = 1/(2*u_required*volume)
    stiffness = k_required/volume
    einstein = 1/(32*np.pi*target)
    assert np.isclose(inertia/einstein, 1, rtol=1e-14)
    assert np.isclose(stiffness/einstein, 1, rtol=1e-14)
    assert np.isclose(2*u_required*k_required, 1, rtol=1e-14)

    adjacency = [[] for _ in sites]
    for i, j in edges:
        adjacency[i].append(j)
        adjacency[j].append(i)
    distances = []
    for source in range(len(sites)):
        d = [-1]*len(sites)
        d[source] = 0
        queue = deque([source])
        while queue:
            i = queue.popleft()
            for j in adjacency[i]:
                if d[j] < 0:
                    d[j] = d[i]+1
                    queue.append(j)
        assert min(d) == 0
        distances.append(d)
    distances = np.array(distances)
    assert len(edges) == 60
    assert distances.max() == 4
    assert np.array_equal(np.bincount(distances[0]), [1, 12, 6])
    outer = int(np.flatnonzero(np.sum(sites*sites, axis=1) == 4)[0])
    assert np.array_equal(np.bincount(distances[outer]), [1, 4, 9, 4, 1])
    rows = []
    for exponent in (.5, 1.):
        q = alpha**(-exponent)
        bounds = 1/np.sum(q**(2*distances), axis=1)
        # Compare the full graph enumeration with the distance polynomials.
        central = 1/(1+12*q*q+6*q**4)
        outer_bound = 1/(1+4*q*q+9*q**4+4*q**6+q**8)
        assert np.isclose(bounds[0]/central, 1, rtol=1e-14)
        assert np.isclose(bounds.min()/outer_bound, 1, rtol=1e-14)
        rows.append({'max_link_ratio_alpha_exponent':exponent,
                     'central_residue_bound':float(central),
                     'minimum_local_residue_bound':float(outer_bound),
                     'bound_over_target_alpha_G':float(outer_bound/target)})

    # Exhibit a zero mode attaining the bound at moderate q, where numerical
    # diagonalisation is well conditioned. Every FCC edge remains present.
    q = 2.
    w = q**distances[outer]
    incidence = np.zeros((len(edges), len(sites)))
    for row, (i, j) in enumerate(edges):
        incidence[row, i] = 1
        incidence[row, j] = -w[i]/w[j]
    normalised = w/np.linalg.norm(w)
    mass_matrix = incidence.T @ incidence
    eigenvalues, eigenvectors = np.linalg.eigh(mass_matrix)
    residual = float(np.linalg.norm(mass_matrix @ normalised))
    assert residual < 1e-13
    assert eigenvalues[1] > .01  # a single null mode on this connected graph
    assert abs(abs(eigenvectors[:, 0] @ normalised)-1) < 1e-13
    assert abs(normalised[outer]**2-1/np.sum(q**(2*distances[outer]))) < 1e-15

    # All orientations of equal non-unit link magnitudes fail a triangle's
    # cycle-product condition: three signs cannot sum to zero.
    for signs in itertools.product((-1, 1), repeat=3):
        assert sum(signs) != 0
        triangle = np.array([[1, -q**signs[0], 0],
                             [0, 1, -q**signs[1]],
                             [-q**signs[2], 0, 1]])
        assert abs(np.linalg.det(triangle)) > .1

    return {'target_alpha_G':target,
            'fcc_bond_count':len(edges),
            'graph_diameter':int(distances.max()),
            'fcc_laplacian_over_k2':dispersion,
            'required_U_over_E0':float(u_required),
            'required_K_over_E0':float(k_required),
            'speed_over_c_with_ordinary_K':float(np.sqrt(2*u_required)),
            'local_source_mixing_bounds':rows,
            'attained_null_mode_residual':residual}


def physical_response_checks(node_count=19):
    """Test a conditional embedding of the vertex in kinetic coordinates.

    Q=(sqrt(rho)*u,sqrt(J)*phi), x=k*ell, z=nu*ell/c.
    The nonchiral continuum kernel uses J=gamma/c**2. Assigning the PN
    channel vector to these coordinates with real relative phase is an
    additional assumption, not a microscopic derivation of that vertex.
    """
    a = np.pi / (np.pi - 1.0)
    b = 4.0 / (np.pi - 2.0)
    t = np.sqrt(b * (a - 1.0))
    n = 1.0 / np.sqrt(np.pi)
    v = np.array([1.0, n], dtype=complex)
    S = np.outer(v, v.conj())  # alpha factored out
    max_error = 0.0
    rows = []
    for x in (0.001, 0.03, 0.1, 0.3, 1.0):
        for z in (0.0, 0.01, 0.3, 1.0):
            M = np.array([[z*z + a*x*x, 1j*t*x],
                          [-1j*t*x, z*z + x*x + b]])
            # Strip the scalar acoustic denominator for channel comparison.
            # The physical Green function still contains that denominator.
            R = (z*z + x*x) * np.linalg.inv(M)
            contraction = float((v.conj() @ R @ v).real)
            exact = 1.0 + n*n*z*z / (z*z + a*x*x + b)
            max_error = max(max_error, abs(contraction - exact))
            chain = S.copy()
            for _ in range(node_count - 1):
                chain = chain @ R @ S
            endpoint = float(np.trace(chain).real)
            closed = float(np.trace(np.linalg.matrix_power(R @ S, node_count)).real)
            assert np.isclose(endpoint, (1+n*n)*exact**(node_count-1), rtol=3e-14)
            assert np.isclose(closed, exact**node_count, rtol=3e-14)
            if x == 0.1:
                rows.append({'k_ell':x, 'nu_ell_over_c':z,
                             'contraction':contraction,
                             'endpoint_trace_prefactor':endpoint,
                             'fully_propagated_loop_prefactor':closed})
    assert max_error < 1e-14

    # Changing only the vertex phase is a different physical source.
    x = 0.1
    M = np.array([[a*x*x, 1j*t*x], [-1j*t*x, x*x+b]])
    D = np.linalg.inv(M)
    R = x*x*D
    phase_rows = []
    for phase in (0.0, np.pi/2, -np.pi/2):
        vp = np.array([1.0, n*np.exp(1j*phase)])
        value = float((vp.conj() @ R @ vp).real)
        exact = 1 + 2*n*t*x*np.sin(phase)/(b+a*x*x)
        assert abs(value-exact) < 1e-14
        phase_rows.append({'relative_phase':phase, 'contraction':value})

    # Changing field coordinates AND source coordinates preserves the result.
    L = np.array([[2.0, 0.1j], [0.0, 0.7j]])
    M_new = L.conj().T @ M @ L
    v_new = L.conj().T @ v
    original = v.conj() @ np.linalg.solve(M, v)
    transformed = v_new.conj() @ np.linalg.solve(M_new, v_new)
    covariance_error = float(abs(transformed/original-1))
    assert covariance_error < 1e-14

    # A local added translational potential keeps a translational amputated
    # scattering vertex. Rotation appears on the external response lines.
    eu = np.array([1.0, 0.0])
    f = 0.07
    P = np.outer(eu, eu)
    T = f/(1+f*D[0,0])*P
    perturbed = np.linalg.inv(M+f*P)
    resolvent_error = float(np.linalg.norm(perturbed-(D-D@T@D))/np.linalg.norm(perturbed))
    assert resolvent_error < 1e-13
    assert abs(D[1,0]/D[0,0]-1j*t*x/(x*x+b)) < 1e-14

    # A core changing only the slip potential has the same rotational
    # determinant as the vacuum. Test with multiple coupled site coordinates.
    A = np.array([[4.,1.,0.], [1.,5.,.2], [0.,.2,3.]])
    C = np.array([[3.,.2], [.2,2.]])
    B = np.array([[.4,.1], [-.3,.2], [.1,.5]])
    V0 = np.diag([.3,.2,.4])
    Vcore = np.diag([-.1,.6,.5])
    effective = A-B@np.linalg.solve(C,B.T)
    H0 = np.block([[A+V0,B], [B.T,C]])
    Hcore = np.block([[A+Vcore,B], [B.T,C]])
    full_ratio = np.linalg.det(Hcore)/np.linalg.det(H0)
    reduced_ratio = np.linalg.det(effective+Vcore)/np.linalg.det(effective+V0)
    determinant_error = float(abs(full_ratio/reduced_ratio-1))
    assert determinant_error < 1e-14
    return {'max_contraction_error':max_error, 'samples':rows,
            'phase_samples':phase_rows,
            'coordinate_covariance_relative_error':covariance_error,
            'translational_potential_resolvent_relative_error':resolvent_error,
            'core_vacuum_determinant_relative_error':determinant_error}


def sum_rule_checks(node_count=19):
    """The endpoint identity as the inertial limit of the channel contraction.

    In kinetic coordinates Q=(sqrt(rho)*u, sqrt(J)*phi) the contraction
    p(z,x) = v^T R(z,x) v = 1 + N^2 z^2/(z^2 + a x^2 + b) rises from 1 (static)
    to 1+N^2 (inertial).  The inertial limit is the f-sum (Thomas-Reiche-Kuhn)
    rule for the vertex coordinate w = v.Q: sum_n (E_n-E_0)|<n|w|0>|^2 =
    (hbar^2/2) v.v whatever the potential.  The rule is checked on a
    two-dimensional quantum model, with and without anharmonic and Peierls
    terms, by finite differences that converge as the grid is refined.
    """
    import scipy.sparse as sparse
    import scipy.sparse.linalg as splinalg
    a = np.pi / (np.pi - 1.0)
    b = 4.0 / (np.pi - 2.0)
    t = np.sqrt(b * (a - 1.0))
    n2 = 1.0 / np.pi
    v = np.array([1.0, np.sqrt(n2)])

    # The coupling number is the fraction of the clamped-rotation shear
    # stiffness mu+kappa_c lost when the rotation is free: 1 - mu_bar/(mu+kappa_c).
    mu = 1.0
    kappa = 2.0 * mu / (np.pi - 2.0)
    mubar = mu + kappa / 2.0
    assert abs(kappa / (2.0 * (mu + kappa)) - n2) < 1e-15
    assert abs(1.0 - mubar / (mu + kappa) - n2) < 1e-15

    def contraction(z, x):
        M = np.array([[z*z + a*x*x, 1j*t*x], [-1j*t*x, z*z + x*x + b]])
        return float((v @ ((z*z + x*x) * np.linalg.inv(M)) @ v).real)

    static = [contraction(0.0, x) for x in (0.01, 1.0, 2*np.pi*np.sqrt(3))]
    inertial = [contraction(1e6, x) for x in (0.01, 1.0, 2*np.pi*np.sqrt(3))]
    assert max(abs(s - 1.0) for s in static) < 1e-13
    assert max(abs(s - (1.0 + n2)) for s in inertial) < 1e-9
    # Static links and an inertial readout give alpha^19 (1+N^2) exactly.
    sigma = np.outer(v, v)
    R0 = (1.0) * np.linalg.inv(np.array([[a, 1j*t], [-1j*t, 1.0 + b]]))
    chain = sigma.astype(complex)
    for _ in range(node_count - 1):
        chain = chain @ R0 @ sigma
    endpoint = float(np.trace(chain).real)
    assert abs(endpoint - (1.0 + n2)) < 1e-12
    # Readout at one radian per light crossing (reading A's condensate rate)
    # and at (v_p/c)^2 = 1/alpha_G (reading B's).  The continuum value at the
    # latter is only indicative: the node's internal structure is omitted.
    at_compton = contraction(1.0, 0.0)
    assert abs(at_compton - (1.0 + n2/(1.0 + b))) < 1e-14
    # The same rate read as a real frequency, z**2 -> -1 (below the optical gap b).
    at_compton_real = 1.0 + n2*(-1.0)/(-1.0 + b)
    assert abs(at_compton_real - (1.0 - n2/(b - 1.0))) < 1e-15
    z_b = 3.0e40
    reading_b_shortfall = n2 * (a * (2*np.pi*np.sqrt(3))**2 + b) / z_b**2

    # f-sum rule on a two-dimensional model in kinetic coordinates (hbar = 1).
    def ratio(m, anharmonic, half_width=6.0):
        xs = np.linspace(-half_width, half_width, m)
        h = xs[1] - xs[0]
        X, Y = np.meshgrid(xs, xs, indexing='ij')
        V = 0.5*(1.3*X**2 + 2.1*Y**2) + 0.6*X*Y
        if anharmonic:
            V = V + 0.08*X**4 + 0.4*np.cos(2.0*(X + np.sqrt(n2)*Y))
        o = np.ones(m)
        D2 = sparse.diags([-o[:-2]/12, 4*o[:-1]/3, -2.5*o, 4*o[:-1]/3, -o[:-2]/12],
                          [-2, -1, 0, 1, 2]) / h**2
        eye = sparse.identity(m)
        H = (-0.5*(sparse.kron(D2, eye) + sparse.kron(eye, D2))
             + sparse.diags(V.ravel())).tocsc()
        E, psi = splinalg.eigsh(H, k=1, sigma=V.min() - 1.0, which='LM')
        psi = psi[:, 0]

        def first_moment(op):
            phi = op.ravel() * psi
            return phi @ (H @ phi) - E[0] * (phi @ phi)
        return first_moment(X + np.sqrt(n2)*Y) / first_moment(X)

    fsum = {m: (ratio(m, False), ratio(m, True)) for m in (81, 161)}
    for m, (harm, anh) in fsum.items():
        assert abs(harm - anh) < 2e-5
    assert abs(fsum[161][1] - (1.0 + n2)) < 2e-6
    return {'static_contraction':static, 'inertial_contraction':inertial,
            'endpoint_with_static_links_over_alpha19':endpoint,
            'contraction_at_one_radian_per_crossing':at_compton,
            'readout_at_reading_A_over_trace':at_compton/(1.0 + n2),
            'contraction_at_one_radian_per_crossing_real_frequency':at_compton_real,
            'real_frequency_readout_over_trace':at_compton_real/(1.0 + n2),
            'reading_B_shortfall_bound':reading_b_shortfall,
            'f_sum_ratio_by_grid':{str(m):list(r) for m, r in fsum.items()},
            'one_plus_N2':1.0 + n2}


def vertex_selection_checks():
    """Which slide-spin covector has wavelength-independent static links.

    Crystal variables (u, phi), static transverse stiffness
    K(k) = [[mu_tot k^2, i kappa k], [-i kappa k, gamma k^2 + 2 kappa]].
    The compliance along w = vu*u + vphi*phi equals vu^2/(mubar k^2) at every
    k only for vphi/vu = +-N sqrt(gamma/mubar); no microinertia enters.
    The inertial (f-sum) weight of that covector is 1 + N^2 gamma/(J c^2).
    This names the vertex the endpoint form assumes; it is not a derivation of
    it, and it fails once the slide is pinned by a uniform restoring force.
    Units: mu = ell = 1, gamma = mu ell^2, rho = 1 so c^2 = mubar.
    """
    mu, ell = 1.0, 1.0
    kappa = 2.0 * mu / (np.pi - 2.0)
    gamma = mu * ell**2
    mubar, mutot = mu + kappa/2.0, mu + kappa
    n2 = kappa / (2.0 * mutot)
    lever = np.sqrt(n2 * gamma / mubar)

    def weight(ratio, phase, k):
        K = np.array([[mutot*k*k, 1j*kappa*k], [-1j*kappa*k, gamma*k*k + 2*kappa]])
        v = np.array([1.0, ratio*np.exp(1j*phase)])
        return float((v.conj() @ np.linalg.solve(K, v)).real * mubar * k * k)

    ks = (0.01, 0.3, 1.0, 3.0, 2*np.pi*np.sqrt(3))
    selected = [weight(lever, 0.0, k) for k in ks]
    assert max(abs(w - 1.0) for w in selected) < 1e-12
    assert max(abs(weight(-lever, 0.0, k) - 1.0) for k in ks) < 1e-12
    slip_only = [weight(0.0, 0.0, k) for k in ks]
    assert abs(weight(0.0, 0.0, 1e4) - (1.0 - n2)) < 1e-6
    longer = [weight(1.2*lever, 0.0, k) for k in ks]
    phased = [weight(lever, 0.3, k) for k in ks]
    # A real covector with unit weight at one k has it at every k: check one.
    assert abs(weight(lever, 0.0, 0.5) - 1.0) < 1e-13
    dense = np.linspace(0.01, 12.0, 4000)
    phase_shift = [weight(lever, 0.3, k) - 1.0 for k in dense]
    worst = int(np.argmax(phase_shift))
    assert max(slip_only) - min(slip_only) > 0.3
    assert max(longer) - min(longer) > 0.05 and min(longer) >= 1.0 - 1e-12
    assert max(phased) - min(phased) > 0.05
    rho = 1.0
    c2 = mubar / rho
    j_light = gamma / c2
    j_alt = rho * ell**2
    inertial = {name: 1.0 + lever**2 * rho / J
                for name, J in (('J=gamma/c^2', j_light), ('J=rho*ell^2', j_alt))}
    assert abs(inertial['J=gamma/c^2'] - (1.0 + 1.0/np.pi)) < 1e-14
    assert abs(inertial['J=rho*ell^2'] - (1.0 + n2*mu/mubar)) < 1e-14
    return {'lever_arm_over_ell':lever,
            'stiffness_length_sqrt_gamma_over_mubar':float(np.sqrt(gamma/mubar)),
            'selected_static_weights':selected,
            'slip_only_static_weights':slip_only,
            'longer_lever_static_weights':longer,
            'phased_static_weights':phased,
            'phase_0p3_max_shift':float(phase_shift[worst]),
            'phase_0p3_max_shift_at_k_ell':float(dense[worst]),
            'contact_radius_covector_inertial_weight':1.0 + (0.5*ell)**2*rho/j_light,
            'inertial_weight':inertial}


def calculate(node_count=19, alpha=1.0 / 137.035999177):
    v = np.array([1.0, 1.0 / np.sqrt(np.pi)])
    channel_norm = float(v @ v)
    vertex = alpha * np.outer(v, v)
    s = np.ones(node_count) / np.sqrt(node_count)
    coherent = np.outer(s, s)
    pu = np.diag([1.0, 0.0])

    # A coherent superposition in a direct sum remains a one-field amplitude.
    transfer = np.diag([alpha, 0.14530929])
    field = np.kron(s, np.array([1.0, 0.0]))
    direct = np.kron(np.eye(node_count), transfer)
    direct_residual = float(np.linalg.norm(direct @ field - alpha * field))

    # Minimal graph-coupled scalar extension, M U''=(V I+k L)U.
    # Internal difference springs leave the uniform mode unchanged.
    graph_residual = None
    if node_count == 19:
        sites, edges = fcc_cluster_graph()
        laplacian = np.zeros((19, 19))
        for a, b in edges:
            laplacian[a, a] += 1
            laplacian[b, b] += 1
            laplacian[a, b] -= 1
            laplacian[b, a] -= 1
        embedding = np.zeros((38, 2))
        embedding[:19, 0] = s
        embedding[19:, 1] = s
        residuals = []
        for x in (0., .13, .5, .91):
            potential = 3.7*np.sin(np.pi*x)**2
            one = np.array([[0., 1.], [potential, 0.]])
            many = np.block([[np.zeros((19,19)), np.eye(19)],
                             [potential*np.eye(19)+2.3*laplacian,
                              np.zeros((19,19))]])
            residuals.append(np.linalg.norm(many@embedding-embedding@one))
        graph_residual = float(max(residuals))
        assert graph_residual < 1e-13

    # Serial rank-one vertices and a translation-projected chain are different
    # operators. Compute both before comparing with their closed forms.
    serial = np.linalg.matrix_power(vertex, node_count)
    projected = vertex.copy()
    for _ in range(node_count - 1):
        projected = projected @ pu @ vertex
    proposed = alpha**node_count * np.kron(coherent, np.outer(v, v))
    serial_weight = float(np.trace(serial) / alpha**node_count)
    projected_weight = float(np.trace(projected) / alpha**node_count)
    proposed_weight = float(np.trace(proposed) / alpha**node_count)
    assert direct_residual < 1e-14
    assert np.isclose(serial_weight, channel_norm**node_count, rtol=2e-14)
    assert np.isclose(projected_weight, channel_norm, rtol=2e-14)
    assert np.isclose(proposed_weight, channel_norm, rtol=2e-14)

    # Explicit small tensor products verify the independent-product identity;
    # the N=19 result follows without allocating a 2**19 square matrix.
    tensor_errors = []
    for n in (2, 3, 4):
        product = vertex.copy()
        for _ in range(n - 1):
            product = np.kron(product, vertex)
        tensor_errors.append(float(abs(np.trace(product)/(alpha*channel_norm)**n - 1)))
    assert max(tensor_errors) < 1e-13

    # R_eta interpolates between translation projection and identity.
    # Sigma R Sigma = alpha*(v.T R v)*Sigma; each intermediate R matters.
    kernels = []
    for eta in (0.0, 0.01, 0.1, 1.0):
        R = np.diag([1.0, eta])
        chain = vertex.copy()
        for _ in range(node_count - 1):
            chain = chain @ R @ vertex
        weight = float(np.trace(chain) / alpha**node_count)
        exact = channel_norm * (1.0 + eta/np.pi)**(node_count-1)
        assert np.isclose(weight, exact, rtol=2e-14)
        kernels.append({'rotational_kernel_weight':eta,'trace_prefactor':weight})

    # A homogeneous equation A X = 0 is unchanged by A -> z*A, whereas the
    # sourced response A X = J changes by 1/z. Its pole residue is extra data.
    A = np.array([[3., -1.], [-1., 2.]])
    source = np.array([1., 0.])
    z = 7.
    response = np.linalg.solve(A, source)
    scaled_response = np.linalg.solve(z*A, source)
    assert np.allclose(scaled_response, response/z, rtol=1e-14)

    return {
        'node_count':node_count,
        'single_amplitude':alpha,
        'independent_product_amplitude':alpha**node_count,
        'coherent_linear_amplitude':alpha,
        'coherent_linear_residual':direct_residual,
        'difference_spring_coherent_reduction_residual':graph_residual,
        'single_trace':channel_norm,
        'serial_or_independent_vertex_trace_prefactor':serial_weight,
        'translation_projected_trace_prefactor':projected_weight,
        'proposed_collective_trace_prefactor':proposed_weight,
        'tensor_product_relative_errors':tensor_errors,
        'intermediate_kernels':kernels,
        'source_response_rescaling':float(scaled_response[0]/response[0]),
        'conditional_physical_response':physical_response_checks(node_count),
        'sum_rule_readout':sum_rule_checks(node_count),
        'vertex_selection':vertex_selection_checks(),
        'microscopic_normalisation':microscopic_normalisation_checks(alpha),
    }


if __name__ == '__main__':
    print(json.dumps(calculate(), indent=2))
