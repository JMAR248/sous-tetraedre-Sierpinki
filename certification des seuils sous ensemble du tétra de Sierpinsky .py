from fractions import Fraction
import math
import mpmath as mp

mp.mp.dps = 60
alpha = mp.log(4) / mp.log(8)  # Exactement 2/3
M = 4**4  # 256 tetraedres d'ordre 4

# Sommets du tetraedre regulier T0 d'arete 1 dans R^3.
# Convention : (x, y, z) represente le point (x, y*sqrt(3), z*sqrt(6)).
T0 = (
    (Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(1, 2), Fraction(1, 2), Fraction(0)),
    (Fraction(1, 2), Fraction(1, 6), Fraction(1, 3)),
)


def g1(p):
    return (p[0] / 8, p[1] / 8, p[2] / 8)


def g2(p):
    return (p[0] / 8 + Fraction(7, 8), p[1] / 8, p[2] / 8)


def g3(p):
    return (p[0] / 8 + Fraction(7, 16), p[1] / 8 + Fraction(7, 16), p[2] / 8)


def g4(p):
    return (
        p[0] / 8 + Fraction(7, 16),
        p[1] / 8 + Fraction(7, 48),
        p[2] / 8 + Fraction(7, 24),
    )


G = (None, g1, g2, g3, g4)


# ------------------------------------------------------------
# Tetraedres d'ordre r.
# Convention : T_{i1...ir} = g_{i1} o ... o g_{ir}(T).
# ------------------------------------------------------------
def tetraedres(r):
    L = [{"v": T0, "a": ""}]
    for _ in range(r):
        L = [
            {"v": tuple(G[i](p) for p in T["v"]), "a": str(i) + T["a"]}
            for T in L
            for i in (1, 2, 3, 4)
        ]
    return L


# ------------------------------------------------------------
# Geometrie 3D exacte (fractions rationnelles uniquement).
# ------------------------------------------------------------
def vec(a, b):
    return (b[0] - a[0], b[1] - a[1], b[2] - a[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(c, u):
    return (c * u[0], c * u[1], c * u[2])


def dot(u, v):
    return u[0] * v[0] + 3 * u[1] * v[1] + 6 * u[2] * v[2]


def d2p(a, b):
    u = vec(a, b)
    return dot(u, u)


def point_segment_d2(p, a, b):
    ab = vec(a, b)
    ap = vec(a, p)
    q = dot(ab, ab)
    if q == 0:
        return d2p(p, a)
    t = dot(ap, ab)
    if t <= 0:
        return d2p(p, a)
    if t >= q:
        return d2p(p, b)
    return dot(ap, ap) - (t * t) / q


def segment_segment_d2(p1, q1, p2, q2):
    """Distance au carre exacte entre deux segments 3D [p1, q1] et [p2, q2]."""
    d1 = vec(p1, q1)
    d2 = vec(p2, q2)
    r = vec(p2, p1)
    a = dot(d1, d1)
    e = dot(d2, d2)
    f = dot(d2, r)

    c = dot(d1, r)
    b = dot(d1, d2)
    denom = a * e - b * b

    if denom != 0:
        s = (b * f - c * e) / denom
        if s < 0:
            s = Fraction(0)
        elif s > 1:
            s = Fraction(1)
    else:
        s = Fraction(0)

    t = (b * s + f) / e
    if t < 0:
        t = Fraction(0)
        s = -c / a
        if s < 0:
            s = Fraction(0)
        elif s > 1:
            s = Fraction(1)
    elif t > 1:
        t = Fraction(1)
        s = (b - c) / a
        if s < 0:
            s = Fraction(0)
        elif s > 1:
            s = Fraction(1)

    c1 = add(p1, scale(s, d1))
    c2 = add(p2, scale(t, d2))
    return d2p(c1, c2)


def point_triangle_d2(p, a, b, c):
    """Distance au carre exacte entre un point p et un triangle 3D (a, b, c)."""
    ab = vec(a, b)
    ac = vec(a, c)
    ap = vec(a, p)

    d1 = dot(ab, ap)
    d2 = dot(ac, ap)
    if d1 <= 0 and d2 <= 0:
        return d2p(p, a)

    bp = vec(b, p)
    d3 = dot(ab, bp)
    d4 = dot(ac, bp)
    if d3 >= 0 and d4 <= d3:
        return d2p(p, b)

    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        proj = add(a, scale(v, ab))
        return d2p(p, proj)

    cp = vec(c, p)
    d5 = dot(ab, cp)
    d6 = dot(ac, cp)
    if d6 >= 0 and d5 <= d6:
        return d2p(p, c)

    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        proj = add(a, scale(w, ac))
        return d2p(p, proj)

    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        proj = add(b, scale(w, vec(b, c)))
        return d2p(p, proj)

    denom = va + vb + vc
    v = vb / denom
    w = vc / denom
    proj = add(a, add(scale(v, ab), scale(w, ac)))
    return d2p(p, proj)


def edges(V):
    return (
        (V[0], V[1]),
        (V[0], V[2]),
        (V[0], V[3]),
        (V[1], V[2]),
        (V[1], V[3]),
        (V[2], V[3]),
    )


def faces(V):
    return (
        (V[0], V[1], V[2]),
        (V[0], V[1], V[3]),
        (V[0], V[2], V[3]),
        (V[1], V[2], V[3]),
    )


def tetra_d2(A, B):
    """Distance au carre exacte entre deux tetraedres reguliers disjoints A et B."""
    VA, VB = A["v"], B["v"]
    EA, EB = edges(VA), edges(VB)
    FA, FB = faces(VA), faces(VB)

    min_d2 = min(
        segment_segment_d2(a1, a2, b1, b2) for a1, a2 in EA for b1, b2 in EB
    )
    for p in VA:
        for f0, f1, f2 in FB:
            d = point_triangle_d2(p, f0, f1, f2)
            if d < min_d2:
                min_d2 = d
    for p in VB:
        for f0, f1, f2 in FA:
            d = point_triangle_d2(p, f0, f1, f2)
            if d < min_d2:
                min_d2 = d
    return min_d2


# ------------------------------------------------------------
# Graphe et couverture de sommets (identique au code 2D).
# ------------------------------------------------------------
def bits(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask ^= b


def edge(adj, active):
    for v in bits(active):
        N = adj[v] & active
        if N:
            b = N & -N
            return v, b.bit_length() - 1
    return None


def degree(adj, v, active):
    return (adj[v] & active).bit_count()


def matching_lb(adj, active, k):
    n = 0
    while True:
        e = edge(adj, active)
        if e is None:
            return n
        u, v = e
        active &= ~(1 << u)
        active &= ~(1 << v)
        n += 1
        if n > k:
            return n


def vertex_cover(adj, k):
    active0 = (1 << len(adj)) - 1
    memo = set()

    def solve(active, k):
        while True:
            e = edge(adj, active)
            if e is None:
                return True
            if k == 0:
                return False
            forced = None
            for v in bits(active):
                if degree(adj, v, active) > k:
                    forced = v
                    break
            if forced is None:
                break
            active &= ~(1 << forced)
            k -= 1
        if matching_lb(adj, active, k) > k:
            return False
        state = (active, k)
        if state in memo:
            return False
        memo.add(state)
        u, v = edge(adj, active)
        if degree(adj, v, active) > degree(adj, u, active):
            u, v = v, u
        return solve(active & ~(1 << u), k - 1) or solve(
            active & ~(1 << v), k - 1
        )

    return solve(active0, k)


# ------------------------------------------------------------
# Distances exactes, decouverte et passage s -> N(P).
# ------------------------------------------------------------
def all_distances(T):
    return [
        (tetra_d2(T[i], T[j]), i, j)
        for i in range(len(T))
        for j in range(i + 1, len(T))
    ]


def graph_ge(n, D, s2):
    A = [0] * n
    for d2, i, j in D:
        if d2 >= s2:
            A[i] |= 1 << j
            A[j] |= 1 << i
    return A


def discover(T, D, m):
    n = len(T)
    k = n - m
    values = sorted({d for d, _, _ in D})
    lo, hi = 0, len(values) - 1
    best = None
    while lo <= hi:
        q = (lo + hi) // 2
        s2 = values[q]
        if not vertex_cover(graph_ge(n, D, s2), k):
            best = s2
            lo = q + 1
        else:
            hi = q - 1
    return best


def mpq(q):
    return mp.mpf(q.numerator) / q.denominator


def Nmin(s2):
    s = mp.sqrt(mpq(s2))
    x = M * s**alpha
    return int(mp.floor(x)) + 1, x


def cells(N, r):
    c = 4 ** (4 - r)
    return (N + c - 1) // c


# ============================================================
# CHAINE COMPLETE (ORDRE 1 A 4)
# ============================================================
print("CHAINE EXACTE DE DECOUVERTE ET DE CERTIFICATION (TETRAEDRE 3D)\n")
print("N(P) est le nombre de tetraedres d'ordre 4 constituant P\n")

s2 = Fraction(9, 16)  # s1 = 3/4 -> s1^2 = 9/16
N, x = Nmin(s2)
print("s1 = 3/4")
print("N(P) >=", N)

for r in range(2, 4):
    m = cells(N, r)
    T = tetraedres(r)
    D = all_distances(T)
    s2 = discover(T, D, m)
    k = len(T) - m
    certified = not vertex_cover(graph_ge(len(T), D, s2), k)
    if not certified:
        raise RuntimeError(f"Echec de certification pour s{r}")
    N, x = Nmin(s2)
    print()
    print(f"ordre {r}")
    print(f"m{r} =", m)
    print(f"s{r}^2 =", s2)
    print(f"s{r} ~= {math.sqrt(float(s2)):.15f}")
    print("borne inferieure certifiee")
    print("N(P) >=", N)

# ------------------------------------------------------------
# Etape finale (ordre 4).
# ------------------------------------------------------------
print()
print("ordre 4")
if N < 256:
    raise RuntimeError("La chaine n'atteint pas la contradiction.")
print("N(P) >=", N)
print("Les 256 tetraedres d'ordre 4 sont tous presents (P = S_4).")
print("Deux tetraedres extremes au moins sont presents.")
print("D(P) = 1.")
print("Q(P) < 1 impliquerait W(P) > 1 : contradiction.")
