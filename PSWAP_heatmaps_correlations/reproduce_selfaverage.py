"""
Independent reconstruction of the PSWAP late-time space-time correlation
C(dx,dt) to cross-check the archived multi-seed data and quantify
self-averaging (Referee 2, point 1).

Model: charge-1 PSWAP brickwork on an N=30 ring, theta=pi/300, initial
operator Z on site N/2. Operator tracked as coefficients over the 2N^2 charge-1
orbit. Occupation d_i(t) = sum_{alpha nonidentity at i} |c_alpha|^2 (matches
core_nonumba.compute_observables site_density). C(dx,dt) matches
examples/_spacetime_corr_v3_latetime.py:compute_spacetime_corr (centered,
normalized so C(0,0)=1), on the late-time window t > 2*t_w.

Rules:
  k0: iid P(A)=0.75  -> R_A ~ Geom(0.25) [E=4], R_B ~ Geom(0.75) [E=4/3]
  k2: symmetric, pg=0.5, p_rr=5/6 -> R_A,R_B ~ _k2_runlength [E=4]  (both parities)

Usage: python reproduce_selfaverage.py k0   (or k2, the default).
"""
import numpy as np
import scipy.sparse as sp
import time, os, sys

CASE = sys.argv[1] if len(sys.argv) > 1 else "k2"
N = 30
THETA = np.pi/300
N_STEPS = 2000
INIT_SITE = N//2
MAX_DX = 15
MAX_DT = 80
HERE = os.path.dirname(__file__)
DATA = {"k0": "k0_ER4_multiseed_data.npz", "k2": "k2_ER4_multiseed_data.npz"}[CASE]
DATA = os.path.join(HERE, DATA)

# ---- run-length samplers ----
def sample_runs(rng):
    if CASE == "k0":
        return int(rng.geometric(0.25)), int(rng.geometric(0.75))
    # k2 symmetric: pg=0.5, p_rr=5/6 for both parities
    def k2(pg, prr):
        if rng.random() < pg: return 1
        r = 2
        while rng.random() < prr: r += 1
        return r
    return k2(0.5, 5/6), k2(0.5, 5/6)

# ---- single-pair PSWAP(phi) Pauli transition (charge-1 symmetric gate) ----
def pair_T(phi):
    c, s = np.cos(phi), np.sin(phi)
    U = np.array([[1,0,0,0],[0,c,1j*s,0],[0,1j*s,c,0],[0,0,0,1]], dtype=complex)
    pa = [np.eye(2), np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]),
          np.array([[1,0],[0,-1]])]
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

def build_nonI(basis, N):
    M = np.zeros((len(basis), N), dtype=bool)
    for j, s in enumerate(basis):
        for i in range(N):
            if s[i] != 0: M[j, i] = True
    return M

def evolve(layers, basis, idx, nonI, hamw, cache):
    """layers: iterable of (R_A, R_B). Returns site-density D and weight W per step."""
    K = len(basis)
    v = np.zeros(K, dtype=complex)
    v[idx[tuple([3 if i==INIT_SITE else 0 for i in range(N)])]] = 1.0
    n = len(layers)
    D = np.zeros((n+1, N)); W = np.zeros(n+1)
    p = np.abs(v)**2; D[0] = nonI.T @ p; W[0] = hamw @ p
    for t, (RA, RB) in enumerate(layers, 1):
        for (off, R) in ((0, int(RA)), (1, int(RB))):
            key = (off, R)
            if key not in cache: cache[key] = build_layer_matrix(basis, idx, off, R*THETA, N)
            v = cache[key] @ v
        p = np.abs(v)**2; D[t] = nonI.T @ p; W[t] = hamw @ p
    return D, W

def compute_spacetime_corr(D, t_start, t_end):
    hm = D[t_start:t_end, :].astype(float)
    mu = hm.mean(); var = ((hm - mu)**2).mean(); hmc = hm - mu; T = hmc.shape[0]
    C = np.zeros((2*MAX_DT+1, 2*MAX_DX+1))
    for it, dt in enumerate(range(-MAX_DT, MAX_DT+1)):
        if dt >= 0: t1 = hmc[:T-dt, :]; t2 = hmc[dt:, :]
        else:       t1 = hmc[-dt:, :];  t2 = hmc[:T+dt, :]
        for ix, dx in enumerate(range(-MAX_DX, MAX_DX+1)):
            C[it, ix] = (t1 * np.roll(t2, -dx, axis=1)).mean()
    return C/var if var > 1e-12 else C

def correlation(D, W):
    tw = int(np.argmax(W/N >= 0.5))
    t_start = min(2*tw, D.shape[0]-1 - MAX_DT - 100)
    return compute_spacetime_corr(D, int(t_start), D.shape[0]), tw

if __name__ == "__main__":
    t0 = time.time()
    basis, idx = build_basis(N)
    nonI = build_nonI(basis, N); hamw = nonI.sum(axis=1).astype(float)
    cache = {}
    print(f"CASE={CASE}  orbit K={len(basis)}  build {time.time()-t0:.1f}s")
    d = np.load(DATA, allow_pickle=True); ref = d["corr_avg"]

    # (1) TIGHT per-seed validation: reproduce archived corr_seed from archived run lengths
    print("\n(1) per-seed reproduction from archived run-length sequences:")
    for s in [42, 43, 44, 45, 46]:
        layers = d[f"layers_seed{s}"]
        D, W = evolve(layers, basis, idx, nonI, hamw, cache)
        C, tw = correlation(D, W)
        r = np.corrcoef(C.ravel(), d[f"corr_seed{s}"].ravel())[0, 1]
        print(f"    seed {s}: r(mine, archived corr_seed{s}) = {r:.3f}   peak={np.abs(C).max():.3f}")

    # (2) FRESH seeds: self-averaging + match to corr_avg
    print("\n(2) fresh independent seeds:")
    n_seeds = 20; Cs = []
    rngmaster = np.random.default_rng(12345)
    for k in range(n_seeds):
        rng = np.random.default_rng(2000+k)
        layers = [sample_runs(rng) for _ in range(N_STEPS)]
        D, W = evolve(layers, basis, idx, nonI, hamw, cache)
        C, tw = correlation(D, W); Cs.append(C)
    Cs = np.array(Cs); mine = Cs.mean(axis=0)
    r = np.corrcoef(ref.ravel(), mine.ravel())[0, 1]
    print(f"    {n_seeds}-seed mean vs archived corr_avg: r={r:.3f}  "
          f"peak mine={np.abs(mine).max():.3f} archived={np.abs(ref).max():.3f}")
    peak = np.abs(mine).max()
    sub = np.array([Cs[rngmaster.choice(n_seeds,5,replace=False)].mean(axis=0) for _ in range(200)])
    rms = np.sqrt(((sub-mine)**2).mean(axis=(1,2))).mean()
    pair = np.corrcoef(Cs[0].ravel(), Cs[1].ravel())[0,1]
    print(f"    self-averaging: 5-seed-mean RMS spread = {100*rms/peak:.1f}% of peak; "
          f"single-seed vs single-seed r = {pair:.3f}")
    print(f"total {time.time()-t0:.1f}s")
