"""fig_scrambling_N512.png: SWAP/WC scrambling time t_scr/N vs f_A at N=512.

Plotting values copied from quantum_batching/examples/_fig_scrambling_standalone.py.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps

ps.apply()

d512 = np.load(os.path.join(ps.DATA, "fig_scrambling_N512", "iid_sweep_N512.npz"))
N512 = int(d512['N'])
iid_fA, iid_t = [], []
for fA in d512['frac_A_values']:
    tag = f'fA{fA:.2f}'.replace('.', '')
    iid_fA.append(fA)
    iid_t.append(np.mean(d512[f'iid_{tag}_scr_times']) / N512)

# k=1 at N=512: (fA, lambda, t_scr/N) from control grid + asymmetric sweep
k1_grid_512 = [
    (0.50, -0.30, 0.207), (0.50, -0.15, 0.208), (0.50, 0.00, 0.207),
    (0.50, 0.20, 0.209), (0.50, 0.40, 0.207), (0.50, 0.60, 0.213),
    (0.67, -0.30, 0.198), (0.67, -0.15, 0.200), (0.67, 0.00, 0.199),
    (0.67, 0.20, 0.200), (0.67, 0.40, 0.205), (0.67, 0.60, 0.212),
    (0.75, -0.30, 0.197), (0.75, -0.15, 0.198), (0.75, 0.00, 0.200),
    (0.75, 0.20, 0.202), (0.75, 0.40, 0.207), (0.75, 0.60, 0.218),
    (0.83, -0.15, 0.201), (0.83, 0.00, 0.203), (0.83, 0.20, 0.213),
    (0.83, 0.40, 0.217), (0.83, 0.60, 0.238),
    (0.90, 0.00, 0.225), (0.90, 0.20, 0.245),
    (0.90, 0.40, 0.261), (0.90, 0.60, 0.335),
    (0.95, 0.00, 0.338), (0.95, 0.20, 0.391),
    (0.95, 0.40, 0.485), (0.95, 0.60, 0.638),
]
k1_extra_512 = [
    (0.50, -0.80, 0.208), (0.50, 0.50, 0.215), (0.50, 0.90, 0.280),
    (0.833, -0.15, 0.203), (0.833, 0.50, 0.232), (0.833, 0.80, 0.347),
    (0.950, 0.00, 0.339), (0.950, 0.50, 0.547), (0.950, 0.80, 1.040),
]
k1_by_fA = {}
for fA, lam, t in k1_grid_512 + k1_extra_512:
    k1_by_fA.setdefault(round(fA, 2), []).append(t)
k1_fA = sorted(k1_by_fA)
k1_t = [min(k1_by_fA[f]) for f in k1_fA]

k2_best_512 = {0.833: 0.202, 0.950: 0.333}

fig, ax = plt.subplots(figsize=(0.4 * ps.TEXTWIDTH, 2.05))

cliff_baseline = 0.27
pure_B = 0.248
ax.axhline(cliff_baseline, color='gray', ls=':', lw=1.0, zorder=0)
ax.axhline(pure_B, color='gray', ls='--', lw=0.9, zorder=0, alpha=0.6)

ax.plot(iid_fA, iid_t, 'o-', color='#2060c0', ms=3.8, lw=1.3,
        label='i.i.d. (k=0)', zorder=5)
ax.plot(k1_fA, k1_t, 's-', color='#d04020', ms=3.2, lw=1.1,
        label='best k=1', zorder=4)
k2_fA = sorted(k2_best_512)
ax.plot(k2_fA, [k2_best_512[f] for f in k2_fA], 'D', color='#20a040', ms=4.5,
        label='best k=2', zorder=6)

ax.text(0.30, cliff_baseline + 0.003, 'full Clifford', fontsize=7.5,
        color='dimgray', va='bottom')
ax.text(0.30, pure_B - 0.003, 'pure WC', fontsize=7.5, color='dimgray', va='top')

ax.set_xlabel('$f_A$ (fraction of SWAP letters)')
ax.set_ylabel('$t_{scr}$ / N')
ax.set_xlim(-0.02, 1.0)
ax.set_ylim(0.17, 0.40)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_title('SWAP/WC scrambling time (N = 512)', fontsize=8.5)
ax.legend(loc='upper right', fontsize=7.5)

fig.tight_layout()
ps.save(fig, "fig_scrambling_N512.png")
