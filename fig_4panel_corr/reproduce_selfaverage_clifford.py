"""
Self-averaging check for the CLIFFORD (SWAP/WC composite) late-time space-time
correlation at N=512 (Referee 2, point 1), using the paper's simulation code
(standalone copy in clifford_swap_wc.py, bit-identical to the original).

Runs 20 independent realizations of the i.i.d. composite rule (fA=0.83) at
N=512, computes the late-time C(dx,dt) (t>2*t_scr), and quantifies how tightly
a 5-realization mean is determined. Compares the 20-seed mean to the archived
comp_C_iid.

"""
import numpy as np, sys, os, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clifford_swap_wc import simulate_swap_wc   # composite (double-layer) sim

N = 512
FRAC_A = 0.83
MAX_LETTERS = 3 * N          # = 1536, matches archived comp_max
MAX_DX = 40; MAX_DT = 40     # archived comp_C_iid is 81x81
N_SEEDS = 20
HERE = os.path.dirname(__file__)

def make_iid_seq(frac_A, length, rng):
    return (rng.random(length) >= frac_A).astype(int)   # 0=A (SWAP), 1=B (WC)

def compute_spacetime_corr(heatmap, max_dx, max_dt, t_start, t_end):
    # identical to examples/_spacetime_corr_v3_latetime.py
    hm = heatmap[t_start:t_end, :].astype(float); T = hm.shape[0]
    mu = np.mean(hm); var = np.mean((hm-mu)**2); hmc = hm - mu
    C = np.zeros((2*max_dt+1, 2*max_dx+1))
    for it, dt in enumerate(range(-max_dt, max_dt+1)):
        if dt >= 0: t1 = hmc[:T-dt, :]; t2 = hmc[dt:, :]
        else:       t1 = hmc[-dt:, :];  t2 = hmc[:T+dt, :]
        for ix, dx in enumerate(range(-max_dx, max_dx+1)):
            C[it, ix] = np.mean(t1 * np.roll(t2, -dx, axis=1))
    return C/var if var > 1e-10 else C

if __name__ == "__main__":
    t0 = time.time()
    Cs = []
    for r in range(N_SEEDS):
        seq_rng = np.random.default_rng(7000 + r*100)
        gate_rng = np.random.default_rng(7000 + r*100 + 1)
        seq = make_iid_seq(FRAC_A, MAX_LETTERS, seq_rng)
        res = simulate_swap_wc(N, MAX_LETTERS, seq, gate_rng,
                               threshold=0.5, record_heatmap=True)
        hm = res['heatmap']
        t_scr = res['scr_time'] if res['hit'] else MAX_LETTERS
        t_start = min(2*t_scr, MAX_LETTERS - MAX_DT - 100)
        C = compute_spacetime_corr(hm, MAX_DX, MAX_DT, int(t_start), MAX_LETTERS)
        Cs.append(C)
        if r < 3 or r == N_SEEDS-1:
            print(f"  seed {r}: t_scr={t_scr}, {time.time()-t0:.1f}s")
    Cs = np.array(Cs); mine = Cs.mean(axis=0)
    print(f"{N_SEEDS} seeds done, total {time.time()-t0:.1f}s")

    # compare to archived comp_C_iid
    d = np.load(os.path.join(HERE, "fig_4panel_corr_data.npz"), allow_pickle=True)
    ref = d["comp_C_iid"]
    if ref.shape == mine.shape:
        r = np.corrcoef(ref.ravel(), mine.ravel())[0,1]
        print(f"\n20-seed mean vs archived comp_C_iid: r={r:.3f}  "
              f"peak mine={np.abs(mine).max():.3f} archived={np.abs(ref).max():.3f}")
    else:
        print(f"\nshape mismatch: mine {mine.shape} vs archived {ref.shape}")

    # self-averaging: spread of 5-seed subsample means
    peak = np.abs(mine).max()
    mr = np.random.default_rng(0)
    sub = np.array([Cs[mr.choice(N_SEEDS,5,replace=False)].mean(0) for _ in range(300)])
    rms = np.sqrt(((sub-mine)**2).mean(axis=(1,2))).mean()
    pair = np.corrcoef(Cs[0].ravel(), Cs[1].ravel())[0,1]
    print(f"self-averaging: 5-seed-mean RMS spread = {100*rms/peak:.1f}% of peak; "
          f"single-seed vs single-seed r = {pair:.3f}")
