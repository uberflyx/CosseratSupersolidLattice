"""Every constant entering the light vector sector, classified by what kind of thing it is.

WHY CLASSIFY BEFORE MECHANISING.  If the lattice is ontological then every genuine
degree of freedom means something, and finding the mechanism behind each is the
programme.  Not every number in a calculation is a degree of freedom, though, and three
kinds are traps:

  A SCHEME CONVENTION is a choice of bookkeeping, not a fact about the world.  The
  Breit-Wigner mass and the pole mass of the same resonance differ by 12 MeV, three
  residue conventions for c_0 span 0.087 to 0.112, and Gamma_ee differs by 3% according
  to whether vacuum polarisation is removed.  Attaching a mechanism to a scheme choice
  invents physics for a convention, and it will not survive a referee who uses the other
  convention.

  A COMPOSITE is a number already determined by others in the table.  Mechanising it
  separately double-counts, which is how a set of agreements comes to look stronger
  than the evidence behind it.  Five such were found in this sector in a single audit.

  A COMPARISON VALUE belongs to the experiment, not the framework.  It needs verifying,
  not explaining.

What is left after those are removed is the real target list, and it is short.

STATUS CODES
  derived   fixed by alpha and m_e through a stated construction; already mechanised
  input     a measurement the framework imports and does not claim to produce
  fitted    adjusted to data; a genuine target, since a fitted number is a missing
            mechanism wearing a value
  stated    asserted with a reason but not derived; a target, and a weaker claim than
            the corpus sometimes makes for it
  scheme    a convention; NOT a target
  composite determined by other entries; NOT a separate target
  compare   an experimental value used for comparison; NOT a target
"""

# (name, value, status, what it physically is / why it is not a target)
INVENTORY = [
    # ---- the two inputs -------------------------------------------------
    ("alpha", "1/137.035999177", "derived",
     "Peierls-Nabarro tunnelling amplitude of a screw dislocation; the published result"),
    ("m_e", "0.51099895069 MeV", "input",
     "the single scale-setting input; sets SI, cancels from every ratio"),

    # ---- scales, all pure numbers times m_0 ------------------------------
    ("m_0", "m_e/alpha = 70.025 MeV", "derived", "node mass"),
    ("sigma", "(2 pi m_0)^2", "derived",
     "ribbon energy per unit length: fault energy gamma folded with the equilibrium"
     " width w of the partial pair, sigma = gamma w"),
    ("Lambda_QCD", "pi m_0", "derived", "half-wave probe on the ribbon"),
    ("E_BZ", "pi sqrt2 m_0", "derived", "D4 zone face, shortest reciprocal vector / 2"),
    ("L*", "2 m_q / sigma = 0.642 fm", "derived",
     "ribbon length at which a fresh partial-antipartial pair pays for itself"),

    # ---- graph and lattice integers -------------------------------------
    ("N_c", "3", "derived", "the three-fold stacking degeneracy ABC; colour"),
    ("Z_co", "12", "derived", "FCC coordination number"),
    ("11", "node count", "derived", "nodes in the crossed-fault cluster carrying the rho"),
    ("4.891", "cluster eigenvalue", "derived",
     "B_3u mode of the 11-node crossed-fault graph; an output of diagonalisation"),
    ("6.241", "cluster eigenvalue", "derived", "B_1u mode of the same graph; the omega"),
    ("m_a1 = 1227.37", "MeV", "derived",
     "18-node Born cluster at its soft rigid-microrotation eigenvalue"),

    # ---- the genuine targets --------------------------------------------
    ("k", "0.957", "fitted",
     "TARGET. One normalisation common to rho, omega and phi leptonic widths. The rho"
     " alone asks 0.966 and the isoscalar pair alone 0.941, agreeing to 2.5% with no"
     " chance to negotiate, so it is one factor and not three. Four candidates excluded:"
     " the chiral-limit decay constant, the N_c exponent in f_pi, per-channel chiral"
     " partners, and tower saturation of the sum rules."),
    ("delta", "3.09 degrees", "fitted",
     "TARGET. omega-phi mixing away from ideal. Needs the off-diagonal element of the"
     " isoscalar mass matrix, which the framework does not supply. Check it already"
     " passes: the nonet angle lands at 38.4 degrees against the standard 39."),
    ("alpha_s(IR)", "0.50", "stated",
     "TARGET. Light-sector coupling in the Cornell solve. Affects the tower, which"
     " carries about 7% of a_mu, so the corpus tolerates it; it is still a number the"
     " lattice should produce."),
    ("m_q", "N_c^2 m_0 / 2 = 315 MeV", "stated",
     "TARGET, partially. The N_c^2 counting is structural, the factor 1/2 is not"
     " derived in the corpus. It now carries more weight than it did, since L* depends"
     " on it directly."),
    ("c_2", "|c_2| <= 0.112", "stated",
     "TARGET. Second recurrence coupling to pions. Bounded, not derived; the bound is"
     " what makes the truncation systematic one-sided."),
    ("B_pipi(rho')", "0.15", "stated",
     "TARGET. Two-pion branching of the first recurrence, asserted from the observed"
     " dominance of 4 pi."),

    # ---- imported measurements ------------------------------------------
    ("m_omega, m_phi", "782.66, 1019.46 MeV", "input",
     "used as pole positions in the dispersion integral"),
    ("Gamma_ee(omega)", "0.641 keV", "input",
     "measured; derived only up to delta, so it is an input until delta is"),
    ("Gamma_rho'", "400 MeV", "input", "recurrence width; sensitivity quoted at +-1 unit"),
    ("quarkonia", "PDG masses and widths", "input", "anchor the Cornell tower ratios"),

    # ---- schemes: not targets -------------------------------------------
    ("BW vs pole mass", "775.26 vs 763.7 MeV", "scheme",
     "the same rho, 12 MeV apart. Which one a formula predicts is a convention"),
    ("c_0 residue convention", "0.087 / 0.099 / 0.112", "scheme",
     "three readings of f_n ~ sqrt(Gamma_ee m^p); the loosest is carried"),
    ("VP-removed vs dressed", "3% on Gamma_ee", "scheme",
     "whether the leptonic width includes photon vacuum polarisation"),
    ("GS vs Breit-Wigner", "line shape", "scheme",
     "checked rather than assumed: the truncation systematic survives the change"),

    # ---- composites: already counted ------------------------------------
    ("Regge slope alpha'", "0.822 GeV^-2", "composite",
     "NOT an independent measure of sigma. A two-parameter fit to a bending trajectory"
     " spreads the bend across slope and intercept; at fixed Adler intercept the first"
     " rung returns the derived tension to 0.23%"),
    ("rho(1450) prediction", "1464 MeV", "composite",
     "the same test as the doublet-centre check, rearranged; the measured centre is"
     " built from the measured rho' and a1"),
    ("half-splitting", "0.318 GeV^2", "composite", "same test again"),
    ("a2 offset", "34 MeV", "composite",
     "one slice of a rung spread over 236 MeV, not an anomaly of its own"),
    ("rung spreads", "17.6%, 3.8%", "composite",
     "spin splitting. The ladder is spin-blind, so it predicts nothing here and this"
     " must not be scored as though it predicted zero"),

    # ---- comparison values ----------------------------------------------
    ("Delta-alpha_had (data)", "0.027609(112)", "compare", "KNT dispersive evaluation"),
    ("a_mu (data / lattice)", "692.78(2.42) / 713(6)", "compare",
     "3.1 sigma apart from each other, so no prediction passes both"),
    ("Gamma_ee(rho) (PDG)", "7.04(6) keV", "compare",
     "constrained fit; the six published determinations agree at chi2/dof = 0.76"),
]

ORDER = ["derived", "fitted", "stated", "input", "scheme", "composite", "compare"]
HEAD = {"derived": "MECHANISED: fixed by alpha and m_e",
        "fitted": "TARGETS, fitted: a value where a mechanism should be",
        "stated": "TARGETS, stated: asserted with a reason, not derived",
        "input": "IMPORTED measurements: honest inputs, not claims",
        "scheme": "SCHEME CONVENTIONS: not targets, and mechanising one would be an error",
        "composite": "COMPOSITES: already determined by other entries; not separate targets",
        "compare": "COMPARISON VALUES: belong to the experiments"}

if __name__ == "__main__":
    for st in ORDER:
        rows = [r for r in INVENTORY if r[2] == st]
        print("\n%s  (%d)" % (HEAD[st], len(rows)))
        print("-" * 78)
        for name, val, _, note in rows:
            print("  %-22s %-26s %s" % (name, val, note if len(note) < 90 else ""))
            if len(note) >= 90:
                words, line = note.split(), "     "
                for w in words:
                    if len(line) + len(w) > 76:
                        print(line)
                        line = "     "
                    line += w + " "
                print(line)
    n_t = sum(1 for r in INVENTORY if r[2] in ("fitted", "stated"))
    print("\n" + "=" * 78)
    print("%d entries.  %d are genuine mechanisation targets." % (len(INVENTORY), n_t))
    print("%d are already mechanised, %d are honest imports, and %d are schemes,"
          % (sum(1 for r in INVENTORY if r[2] == "derived"),
             sum(1 for r in INVENTORY if r[2] == "input"),
             sum(1 for r in INVENTORY if r[2] in ("scheme", "composite", "compare"))))
    print("composites or comparison values, which must NOT be mechanised.")
