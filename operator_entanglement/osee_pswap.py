"""
Operator-space entanglement entropy (OSEE) of the charge-1 PSWAP operator
O(l) = sum_a c_a P_a, across a spatial bipartition of the N=30 ring.

Referee 2, point 4 (conventional chaos diagnostics). We view O(l) as a vector
in operator (Pauli) space, split each Pauli string into its content on the left
arc L and right arc R, reshape the coefficients into c_{a_L,a_R}, and take the
entanglement entropy of the singular values (OSEE). The initial operator sigma_z
is at site 0; the ring is cut into two arcs L={0..N/2-1}, R={N/2..N-1}, so one
boundary sits immediately next to the initial site (site N-1 | site 0) and the
other is antipodal. (On a ring a two-arc cut has TWO boundaries, ~2x the
open-chain single-boundary value.)

Key expectation: the charge-1 orbit has only 2N^2 Pauli strings, so the Schmidt
rank across any cut is <= 2N^2 and OSEE <= log2(2N^2) ~ 10.8 bits at N=30 --
i.e. at most logarithmic in N (sub-volume-law).

"""
import numpy as np
import scipy.sparse as sp

N = 30
THETA = np.pi/300
N_STEPS = 2000
INIT_SITE = 0
SAMPLE_EVERY = 20                       # compute OSEE every this many compressed layers

# ---- run-length samplers for the canonical rules ----
def sampler(rule):
    def k2(pg, prr, rng):
        if rng.random() < pg: return 1
        r = 2
        while rng.random() < prr: r += 1
        return r
    def f(rng):
        if rule == "uniform_ER80":
            return int(rng.integers(1, 160)), int(rng.integers(1, 160))   # mean 80
        if rule == "k0_ER4":
            return int(rng.geometric(0.25)), int(rng.geometric(0.75))     # E[RA]=4, E[RB]=4/3
        if rule == "k1_ER4":
            return int(rng.geometric(0.25)), int(rng.geometric(0.25))     # both mean 4
        if rule == "k2_ER4":
            return k2(0.5, 5/6, rng), k2(0.5, 5/6, rng)                   # mean 4
        raise ValueError(rule)
    return f

# ---- charge-1 PSWAP machinery (same gate as reproduce_selfaverage.py) ----
def pair_T(phi):
    c, s = np.cos(phi), np.sin(phi)
    U = np.array([[1,0,0,0],[0,c,1j*s,0],[0,1j*s,c,0],[0,0,0,1]], dtype=complex)
    pa = [np.eye(2), np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]])]
    T = {}
    for p0 in range(4):
        for p1 in range(4):
            P = np.kron(pa[p0], pa[p1]); Pe = U.conj().T @ P @ U; out = []
            for q0 in range(4):
                for q1 in range(4):
                    Q = np.kron(pa[q0], pa[q1]); co = np.trace(Q.conj().T @ Pe)/4
                    if abs(co) > 1e-12: out.append((q0*4+q1, co))
            T[p0*4+p1] = out
    return T

def apply_pair(state, a, b, T):
    res = {}
    for s, c0 in state.items():
        for o, c2 in T[s[a]*4+s[b]]:
            q0, q1 = divmod(o, 4); ls = list(s); ls[a] = q0; ls[b] = q1
            t = tuple(ls); res[t] = res.get(t, 0) + c0*c2
    return {k: v for k, v in res.items() if abs(v) > 1e-12}

def full_layer(state, offset, T, N):
    cur = state
    for j in range(N//2):
        cur = apply_pair(cur, (offset+2*j) % N, (offset+2*j+1) % N, T)
    return cur

def build_basis(N):
    Tref = pair_T(np.pi/4)
    z0 = tuple([3 if i == INIT_SITE else 0 for i in range(N)])
    orbit = {z0}; changed = True
    while changed:
        changed = False
        for s in list(orbit):
            for off in (0, 1):
                for st2 in full_layer({s: 1.0+0j}, off, Tref, N):
                    if st2 not in orbit: orbit.add(st2); changed = True
    basis = sorted(orbit); return basis, {s: i for i, s in enumerate(basis)}

def build_layer_matrix(basis, idx, offset, phi, N):
    T = pair_T(phi); K = len(basis); rows, cols, data = [], [], []
    for j, s in enumerate(basis):
        for st2, c in full_layer({s: 1.0+0j}, offset, T, N).items():
            rows.append(idx[st2]); cols.append(j); data.append(c)
    return sp.csr_matrix((data, (rows, cols)), shape=(K, K), dtype=complex)

# ---- OSEE across the bipartition ----
def bipartition_indices(basis, Lsites):
    Ls = sorted(Lsites); Rs = sorted(set(range(N)) - set(Lsites))
    lkey, rkey = {}, {}
    row = np.empty(len(basis), int); col = np.empty(len(basis), int)
    for i, s in enumerate(basis):
        a = tuple(s[j] for j in Ls); b = tuple(s[j] for j in Rs)
        if a not in lkey: lkey[a] = len(lkey)
        if b not in rkey: rkey[b] = len(rkey)
        row[i] = lkey[a]; col[i] = rkey[b]
    return row, col, len(lkey), len(rkey)

def osee(v, row, col, nrow, ncol):
    v = v / np.sqrt((np.abs(v)**2).sum())          # ensure normalized
    M = np.zeros((nrow, ncol), dtype=complex)
    M[row, col] = v
    sv = np.linalg.svd(M, compute_uv=False)
    p = sv**2; p = p[p > 1e-14]
    return float(-(p*np.log2(p)).sum())

def hamweight(basis):
    return np.array([sum(1 for x in s if x != 0) for s in basis], float)

if __name__ == "__main__":
    import time
    t0 = time.time()
    basis, idx = build_basis(N)
    K = len(basis)
    Lsites = set(range(0, N//2))                    # {0..14}; boundary N-1|0 next to init, and 14|15
    row, col, nrow, ncol = bipartition_indices(basis, Lsites)
    hamw = hamweight(basis)
    bound = np.log2(2*N*N)
    print(f"orbit K={K}  cut sizes(L,R)=({N//2},{N-N//2})  "
          f"distinct halves(L,R)=({nrow},{ncol})  OSEE bound=log2(2N^2)={bound:.2f} bits")
    print(f"build {time.time()-t0:.1f}s\n")

    rules = ["uniform_ER80", "k0_ER4", "k1_ER4", "k2_ER4"]
    cache = {}
    results = {}
    for rule in rules:
        samp = sampler(rule)
        rng = np.random.default_rng(12345)
        v = np.zeros(K, complex); v[idx[tuple([3 if i==INIT_SITE else 0 for i in range(N)])]] = 1.0
        layers, S, W = [0], [osee(v, row, col, nrow, ncol)], [hamw @ np.abs(v)**2]
        for t in range(1, N_STEPS+1):
            RA, RB = samp(rng)
            for (off, R) in ((0, RA), (1, RB)):
                key = (off, R)
                if key not in cache: cache[key] = build_layer_matrix(basis, idx, off, R*THETA, N)
                v = cache[key] @ v
            if t % SAMPLE_EVERY == 0:
                layers.append(t); S.append(osee(v, row, col, nrow, ncol)); W.append(hamw @ np.abs(v)**2)
        results[rule] = (np.array(layers), np.array(S), np.array(W))
        Sarr = np.array(S); Warr = np.array(W)
        Splateau = Sarr[len(Sarr)//2:].mean()
        print(f"{rule:14s}: OSEE max={Sarr.max():.2f} bits, late-plateau~{Splateau:.2f} bits "
              f"(bound {bound:.2f}); w/N late~{(Warr[len(Warr)//2:].mean()/N):.3f}; {time.time()-t0:.0f}s")
    np.savez("operator_entanglement/osee_results.npz",
             **{f"{r}_layers": results[r][0] for r in rules},
             **{f"{r}_osee": results[r][1] for r in rules},
             **{f"{r}_w": results[r][2] for r in rules}, N=N, bound=bound)
    print(f"\nsaved operator_entanglement/osee_results.npz  total {time.time()-t0:.0f}s")
