"""cover_time_over_N.png: mean cover time / N vs N for k-Markov SWAP walks.

Shaded bands are mean +/- 1 sigma over 500 trials. Open symbols (with arrow)
mark sizes where fewer than 90% of runs covered within 200N steps.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps

ps.apply()

CT = os.path.join(ps.DATA, "cover_time_over_N")
K2S = CT


def load(*paths):
    """Merge several runs; later files take precedence at shared N."""
    rows = {}
    for p in paths:
        d = np.load(p)
        for i, n in enumerate(d['system_sizes']):
            rows[int(n)] = (d['means'][i] / n, d['stds'][i] / n, d['covered_fracs'][i])
    N = np.array(sorted(rows))
    m, s, c = (np.array([rows[n][j] for n in N]) for j in range(3))
    return N, m, s, c


series = [
    ("k=1, D = 9.5", '#d62728', 'o', 6.0,
     load(os.path.join(CT, "k1_p95_20260207_120350.npz"),
          os.path.join(CT, "k1_p95_20260207_111247.npz"))),
    ("k=2, D = 23", '#ff7f0e', 'D', 5.5,
     load(os.path.join(K2S, "k2_g97_r43_20260209_142810.npz"),
          os.path.join(K2S, "k2_g97_r43_extended.npz"))),
    ("k=2, D = 89", '#2ca02c', 's', 5.0,
     load(os.path.join(CT, "k2_g99_r80_20260207_162653.npz"),
          os.path.join(CT, "k2_g99_r80_20260206_182837.npz"))),
    ("k=3, D = 533", '#1f77b4', '^', 5.5,
     load(os.path.join(CT, "k3_g99_p95_e95_20260207_162653.npz"),
          os.path.join(CT, "k3_g99_p95_e95_20260207_070758.npz"))),
]

fig, ax = plt.subplots(figsize=(0.4 * ps.TEXTWIDTH, 2.15))

# k=0 uniform walk on a ring: t_C = N(N-1)/2
Nu = np.array([20, 12000])
ax.plot(Nu, (Nu - 1) / 2, color='gray', lw=1.1, label='k=0 uniform, D = 0.5')

for label, color, marker, ms, (N, m, s, cov) in series:
    ax.fill_between(N, np.maximum(m - s, 1e-3), m + s, color=color, alpha=0.15, lw=0)
    ax.plot(N, m, '-', color=color, lw=1.0)
    ok = cov >= 0.9
    ax.plot(N[ok], m[ok], marker, color=color, ms=ms * 0.62, ls='none', label=label)
    ax.plot(N[~ok], m[~ok], marker, mfc='white', mec=color, mew=1.0,
            ms=ms * 0.62, ls='none')
    for n, y in zip(N[~ok], m[~ok]):
        ax.annotate('', xy=(n, y * 2.2), xytext=(n, y * 1.18),
                    arrowprops=dict(arrowstyle='->', color=color, lw=0.9,
                                    mutation_scale=6, shrinkA=0, shrinkB=0))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(15, 1.2e4)
ax.set_ylim(0.8, 5e3)
ax.set_xlabel('System size N')
ax.set_ylabel('Cover time / N')
ax.grid(True, which='both', alpha=0.25, lw=0.5)
ax.legend(loc='upper left', fontsize=7, framealpha=0.9, handlelength=1.3,
          borderpad=0.3, labelspacing=0.2)

fig.tight_layout()
ps.save(fig, "cover_time_over_N.png")
