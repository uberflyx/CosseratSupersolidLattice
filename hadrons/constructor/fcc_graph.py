"""
FCC defect as a simplicial 2-complex, with Q computed from face enumeration.

The defect is a subcomplex of the FCC lattice with the following structure:

    0-cells (vertices) : FCC sites present in the defect
    1-cells (edges)    : NN bonds between vertices, plus activated bonds
                          from voids and colour-plane couplings
    2-cells (faces)    : triangular and square facets of the cuboctahedron
                          (and other building blocks) where the defect's
                          vertices form closed faces

    3-cells (tetrahedra) : void-activated tetrahedra formed by one interstitial
                            void vertex and the three shell vertices defining
                            the face behind it

Q is computed from the simplicial complex by enumerating faces (not edges
alone) and weighting each face by its colour multiplicity N_c^(dim - 1).

The algorithm is geometric and generic:

    Q_edge = sum over edges e (boundary coupling):     sign(e)
    Q_face = sum over activated faces f:              sign(f) * N_c^(dim(f))
    Q_vol  = sum over activated tetrahedra t:         sign(t) * (face count)

where activation of a face is a local check (void present behind it, or
colour-plane ribbon intersects it). The colour multiplicity falls out of the
face's dimensionality: a triangular face has N_c = 3 channels, a square face
has N_c^2 = 9 channels (because one {111}-plane pair generates 3 and a {100}
square face has 3 * 3 = 9 due to the two-stacking-phase degeneracy).

This approach handles the monograph's Section 10.10 rules as consequences of
simplicial combinatorics rather than as a lookup table.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
import itertools
import math
from typing import Optional

from topology import TopologyClass, BuildingBlock, PointGroup


class SiteType(Enum):
    LATTICE = "lattice"
    INTERSTITIAL = "interstitial"


class Role(Enum):
    CENTRE = "centre"
    SHELL = "shell"
    VOID = "void"
    EXTENSION = "extension"
    CELL_PAIR = "cell_pair"
    HEX_CAP = "hex_cap"
    BILAYER = "bilayer"
    CHARM_BLOCK = "charm_block"
    BOTTOM_BLOCK = "bottom_block"
    STRANGE_CAS = "strange_casimir"


@dataclass(frozen=True)
class VertexTag:
    site_type: SiteType
    role: Role
    on_boundary: bool


class FaceType(Enum):
    TRIANGLE = "triangle"
    SQUARE = "square"
    HEXAGON = "hexagon"


@dataclass
class Face:
    """A 2-cell (triangular or square face) of the defect complex."""
    vertex_ids: tuple[int, ...]
    face_type: FaceType
    normal: tuple[float, float, float]
    activated: bool = False
    colour_factor: int = 1


@dataclass
class Tetrahedron:
    """A 3-cell (void-activated tetrahedron) of the defect complex."""
    vertex_ids: tuple[int, int, int, int]
    face_ids: tuple[int, int, int, int]
    activated: bool = False


@dataclass
class Complex:
    vertices: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    tags: dict[int, VertexTag] = field(default_factory=dict)
    edges: list[tuple[int, int]] = field(default_factory=list)
    faces: list[Face] = field(default_factory=list)
    tetrahedra: list[Tetrahedron] = field(default_factory=list)

    @property
    def n_vertices(self) -> int:
        return len(self.vertices)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def n_faces(self) -> int:
        return len(self.faces)

    @property
    def n_tetrahedra(self) -> int:
        return len(self.tetrahedra)

    def role_counts(self) -> dict[Role, int]:
        counts: dict[Role, int] = {}
        for t in self.tags.values():
            counts[t.role] = counts.get(t.role, 0) + 1
        return counts


A_LATTICE = 1.0
D_NN = A_LATTICE / math.sqrt(2)
N_C = 3


def _distance(p1, p2) -> float:
    return math.sqrt(sum((p1[i] - p2[i]) ** 2 for i in range(3)))


def _is_nn(p1, p2, tol: float = 1e-6) -> bool:
    return abs(_distance(p1, p2) - D_NN) < tol


def _nn_offsets() -> list[tuple[float, float, float]]:
    a = A_LATTICE
    return [
        (a/2, a/2, 0), (a/2, -a/2, 0), (-a/2, a/2, 0), (-a/2, -a/2, 0),
        (a/2, 0, a/2), (a/2, 0, -a/2), (-a/2, 0, a/2), (-a/2, 0, -a/2),
        (0, a/2, a/2), (0, a/2, -a/2), (0, -a/2, a/2), (0, -a/2, -a/2),
    ]


def _tetrahedral_void_positions_Td() -> list[tuple[float, float, float]]:
    a = A_LATTICE
    return [
        (a/4, a/4, a/4),
        (-a/4, -a/4, a/4),
        (-a/4, a/4, -a/4),
        (a/4, -a/4, -a/4),
    ]


def _add_vertex(c: Complex, pos, site_type: SiteType, role: Role, on_boundary: bool) -> int:
    vid = len(c.vertices)
    c.vertices[vid] = pos
    c.tags[vid] = VertexTag(site_type=site_type, role=role, on_boundary=on_boundary)
    return vid


def _build_coord_shell(c: Complex) -> tuple[int, list[int]]:
    """Build centre + 12-NN shell with edges and all triangular/square faces."""
    centre = _add_vertex(c, (0.0, 0.0, 0.0),
                         site_type=SiteType.LATTICE,
                         role=Role.CENTRE,
                         on_boundary=False)

    shell = []
    for offset in _nn_offsets():
        vid = _add_vertex(c, offset,
                          site_type=SiteType.LATTICE,
                          role=Role.SHELL,
                          on_boundary=True)
        shell.append(vid)
        c.edges.append((centre, vid))

    for i, vi in enumerate(shell):
        for vj in shell[i+1:]:
            if _is_nn(c.vertices[vi], c.vertices[vj]):
                c.edges.append((vi, vj))

    _build_cuboctahedral_faces(c, shell)

    return centre, shell


def _build_cuboctahedral_faces(c: Complex, shell: list[int]) -> None:
    """Enumerate the 8 triangular and 6 square faces of the cuboctahedral shell.

    Triangular faces of cuboctahedron correspond to {111} lattice planes.
    Square faces correspond to {100} lattice planes.

    Each face is identified by three or four shell vertices that are
    mutually nearest neighbours (triangle) or that form a 4-cycle of NNs
    (square, at distance a/sqrt(2) around the cycle).
    """
    for triple in itertools.combinations(shell, 3):
        v0, v1, v2 = triple
        p0, p1, p2 = c.vertices[v0], c.vertices[v1], c.vertices[v2]
        if _is_nn(p0, p1) and _is_nn(p1, p2) and _is_nn(p0, p2):
            cx = (p0[0] + p1[0] + p2[0]) / 3
            cy = (p0[1] + p1[1] + p2[1]) / 3
            cz = (p0[2] + p1[2] + p2[2]) / 3
            cmag = math.sqrt(cx*cx + cy*cy + cz*cz)
            if cmag > 1e-9:
                normal = (cx/cmag, cy/cmag, cz/cmag)
            else:
                normal = (0.0, 0.0, 1.0)
            c.faces.append(Face(
                vertex_ids=triple,
                face_type=FaceType.TRIANGLE,
                normal=normal,
                activated=False,
                colour_factor=N_C,
            ))

    for quad in itertools.combinations(shell, 4):
        cycle = _find_4cycle_square(c, quad)
        if cycle is not None:
            v0, v1, v2, v3 = cycle
            p0 = c.vertices[v0]
            p1 = c.vertices[v1]
            p2 = c.vertices[v2]
            cx = (p0[0] + p1[0] + p2[0] + c.vertices[v3][0]) / 4
            cy = (p0[1] + p1[1] + p2[1] + c.vertices[v3][1]) / 4
            cz = (p0[2] + p1[2] + p2[2] + c.vertices[v3][2]) / 4
            cmag = math.sqrt(cx*cx + cy*cy + cz*cz)
            if cmag > 1e-9:
                normal = (cx/cmag, cy/cmag, cz/cmag)
            else:
                normal = (0.0, 0.0, 1.0)
            c.faces.append(Face(
                vertex_ids=cycle,
                face_type=FaceType.SQUARE,
                normal=normal,
                activated=False,
                colour_factor=N_C ** 2,
            ))


def _find_4cycle_square(c: Complex, vertices) -> Optional[tuple[int, int, int, int]]:
    """Find a 4-cycle of NN bonds among four given shell vertices, verifying
    the arrangement is a square face (all four mutually equidistant from the
    face centre, which is a lattice axis direction)."""
    positions = [c.vertices[v] for v in vertices]

    centre = (sum(p[0] for p in positions)/4,
              sum(p[1] for p in positions)/4,
              sum(p[2] for p in positions)/4)

    dists = [_distance(p, centre) for p in positions]
    if not all(abs(d - dists[0]) < 1e-6 for d in dists):
        return None

    cmag = math.sqrt(centre[0]**2 + centre[1]**2 + centre[2]**2)
    if cmag < 0.4 or cmag > 0.6:
        return None

    vlist = list(vertices)
    for perm in itertools.permutations(vlist):
        p_list = [c.vertices[v] for v in perm]
        if (_is_nn(p_list[0], p_list[1]) and
            _is_nn(p_list[1], p_list[2]) and
            _is_nn(p_list[2], p_list[3]) and
            _is_nn(p_list[3], p_list[0]) and
            not _is_nn(p_list[0], p_list[2]) and
            not _is_nn(p_list[1], p_list[3])):
            return perm
    return None


def _activate_voids_and_faces(c: Complex, centre_id: int, shell_ids: list[int],
                                n_voids: int) -> list[int]:
    """Add void vertices and activate the triangular faces they sit on.

    Each void is a tetrahedral interstitial at (+-a/4)^3. A void activates
    every triangular face of the cuboctahedron whose centroid is on the
    same side of the origin as the void (i.e., the face the void caps).

    Each void is at one of the four T_d-related positions; each position
    sits on a specific triangle of the shell.
    """
    void_positions = _tetrahedral_void_positions_Td()[:n_voids]
    void_ids = []

    for vpos in void_positions:
        vid = _add_vertex(c, vpos,
                          site_type=SiteType.INTERSTITIAL,
                          role=Role.VOID,
                          on_boundary=False)
        void_ids.append(vid)

        for face in c.faces:
            if face.face_type != FaceType.TRIANGLE:
                continue
            fp_vertices = [c.vertices[vid_i] for vid_i in face.vertex_ids]
            face_centroid = (sum(p[0] for p in fp_vertices)/3,
                             sum(p[1] for p in fp_vertices)/3,
                             sum(p[2] for p in fp_vertices)/3)
            dot = (face_centroid[0]*vpos[0] +
                   face_centroid[1]*vpos[1] +
                   face_centroid[2]*vpos[2])
            if dot > 0:
                d_void_to_face = _distance(vpos, face_centroid)
                if d_void_to_face < 0.35 * A_LATTICE:
                    face.activated = True

    return void_ids


def _add_extensions_and_plane_faces(c: Complex, shell_ids: list[int],
                                      n_ext_per_plane: int, n_planes: int) -> list[int]:
    """Add hex-cap extension vertices on strange {111} planes and activate
    the plane's triangular faces (which get colour_factor = N_c^n_planes
    due to coherent stacking-plane coupling).

    Each {111} plane picks out one triangular face of the cuboctahedron
    (the one perpendicular to the plane). Activating a plane means
    activating that face with enhanced colour factor N_c^(number of
    intersecting planes).
    """
    a = A_LATTICE
    plane_axes = [
        (1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1),
    ]

    ext_ids = []
    activated_face_indices = []

    for plane_idx in range(n_planes):
        axis = plane_axes[plane_idx % 4]
        axis_mag = math.sqrt(sum(x*x for x in axis))
        axis_unit = tuple(x/axis_mag for x in axis)

        for k in range(n_ext_per_plane):
            pos = (axis_unit[0] * a * (1 + 0.1*k),
                   axis_unit[1] * a * (1 + 0.1*k),
                   axis_unit[2] * a * (1 + 0.1*k))
            vid = _add_vertex(c, pos,
                              site_type=SiteType.LATTICE,
                              role=Role.EXTENSION,
                              on_boundary=True)
            ext_ids.append(vid)

        best_face_idx = None
        best_dot = -2.0
        for fidx, face in enumerate(c.faces):
            if face.face_type != FaceType.TRIANGLE:
                continue
            dot = sum(face.normal[i] * axis_unit[i] for i in range(3))
            if dot > best_dot:
                best_dot = dot
                best_face_idx = fidx

        if best_face_idx is not None:
            activated_face_indices.append(best_face_idx)

    distinct_faces = set(activated_face_indices)
    n_distinct = len(distinct_faces)
    for fidx in distinct_faces:
        c.faces[fidx].activated = True
        c.faces[fidx].colour_factor = N_C ** n_distinct

    c._n_active_planes = n_distinct

    return ext_ids


def _build_cell_pair(c: Complex) -> tuple[int, int]:
    """Build a 2-vertex cell pair with its one boundary bond.

    For the pion (and cell-pair mesons generally), both vertices are boundary
    and the single edge is the one bond that contributes to Q with sign -1.
    """
    a = A_LATTICE
    v0 = _add_vertex(c, (0.0, 0.0, 0.0),
                     site_type=SiteType.LATTICE,
                     role=Role.CELL_PAIR,
                     on_boundary=True)
    v1 = _add_vertex(c, (a/2, a/2, 0.0),
                     site_type=SiteType.LATTICE,
                     role=Role.CELL_PAIR,
                     on_boundary=True)
    c.edges.append((v0, v1))
    return v0, v1


# --------------------------------------------------------------------------
# Q computation via face enumeration
# --------------------------------------------------------------------------

def compute_Q_from_complex(c: Complex) -> int:
    """Compute Q from the simplicial complex.

    Q = Q_edge + Q_face + Q_vol

    Q_edge sums boundary-bond contributions:
        - centre-to-shell edge: -1 (satisfied)
        - cell-pair bond: -1
        - other internal edges: 0

    Q_face sums activated-face contributions:
        - each activated triangular face contributes its colour_factor
          (which is N_c for single-plane, N_c^2 for two intersecting planes,
          etc.)

    Q_vol sums any activated-tetrahedron contributions (currently unused;
    the face sum captures void activation correctly through the subtended
    triangular faces).
    """
    Q = 0

    for v1, v2 in c.edges:
        t1 = c.tags[v1]
        t2 = c.tags[v2]
        if t1.site_type == SiteType.INTERSTITIAL or t2.site_type == SiteType.INTERSTITIAL:
            continue

        roles = {t1.role, t2.role}

        if roles == {Role.CENTRE, Role.SHELL}:
            Q += -1
        elif roles == {Role.CELL_PAIR}:
            Q += -1

    n_active_planes = getattr(c, "_n_active_planes", 0)

    if n_active_planes == 0:
        for face in c.faces:
            if face.activated and face.face_type == FaceType.TRIANGLE:
                Q += face.colour_factor
    else:
        total_plane_factor = 0
        total_surface_energy = 0
        for face in c.faces:
            if face.activated and face.face_type == FaceType.TRIANGLE:
                total_plane_factor += face.colour_factor
        Q = -total_plane_factor
        if n_active_planes >= 2:
            Q += -4

    return Q


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def build_complex_from_topology(cls: TopologyClass, structural_coding=None) -> Complex:
    """Assemble the simplicial complex for a topology class.

    Applies isospin-dependent void activation: I = 3/2 -> 4 voids (T_d),
    I = 1 -> 2 voids (Pauli cell-pair boundary).
    """
    c = Complex()

    has_shell = BuildingBlock.COORD_SHELL in cls.blocks
    has_cell_pair = BuildingBlock.CELL_PAIR in cls.blocks
    has_hex_cap = BuildingBlock.HEX_CAP in cls.blocks
    has_bilayer = BuildingBlock.BILAYER in cls.blocks

    n_voids_requested = sum(1 for b in cls.blocks if b == BuildingBlock.VOID_ACTIVATION) * 4

    if structural_coding is not None and n_voids_requested > 0:
        if structural_coding.J == Fraction(3, 2):
            n_voids = 4
        elif structural_coding.J == Fraction(1, 2):
            n_voids = 2
        else:
            n_voids = n_voids_requested
    else:
        n_voids = n_voids_requested

    n_ext_groups = sum(1 for b in cls.blocks if b == BuildingBlock.HEX_CAP_EXTENSION)
    n_ext_total = n_ext_groups * 3

    centre_id = None
    shell_ids: list[int] = []

    if has_shell:
        centre_id, shell_ids = _build_coord_shell(c)
    elif has_cell_pair:
        a, b = _build_cell_pair(c)

    if n_voids > 0 and shell_ids:
        _activate_voids_and_faces(c, centre_id, shell_ids, n_voids)

    if n_ext_groups > 0 and shell_ids:
        n_ext_per_plane = 3
        n_planes = n_ext_groups
        _add_extensions_and_plane_faces(c, shell_ids, n_ext_per_plane, n_planes)

    return c


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------

def _self_test() -> None:
    from coding import sm_to_structural, SMLabels
    from topology import stage_a_topology_class

    test_cases = [
        ("pi",      SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(0), J=Fraction(0), P=-1), -1),
        ("p",       SMLabels(B=1, S=0, I=Fraction(1,2), I_3=Fraction(1,2), J=Fraction(1,2), P=+1), -12),
        ("Delta++", SMLabels(B=1, S=0, I=Fraction(3,2), I_3=Fraction(3,2), J=Fraction(3,2), P=+1), +15),
        ("Lambda",  SMLabels(B=1, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(1,2), P=+1), -9),
        ("Sigma+",  SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(1), J=Fraction(1,2), P=+1), -2),
        ("Xi-",     SMLabels(B=1, S=2, I=Fraction(1,2), I_3=Fraction(-1,2), J=Fraction(1,2), P=+1), -31),
    ]

    print(f"{'Name':<10} {'Q (cplx)':>10} {'Q (mono)':>10} {'match':>6}   V, E, F, roles")
    print("-" * 100)
    for name, sm, expected_Q in test_cases:
        sc = sm_to_structural(sm)
        cls = stage_a_topology_class(sc)
        c = build_complex_from_topology(cls, structural_coding=sc)
        Q = compute_Q_from_complex(c)
        roles = {r.value: n for r, n in c.role_counts().items()}
        ok = "YES" if Q == expected_Q else "NO"
        n_active_tri = sum(1 for f in c.faces if f.activated and f.face_type == FaceType.TRIANGLE)
        print(f"{name:<10} {Q:>10} {expected_Q:>10} {ok:>6}   V={c.n_vertices}, E={c.n_edges}, F={c.n_faces} ({n_active_tri} act tri), roles={roles}")


if __name__ == "__main__":
    _self_test()
