"""waterfall_dw_sector_Vpaper.png: domain-wall sector weight vs compressed layer
for six canonical cases (theta=pi/300, N=30), offset along the depth axis, with
a dashed line at each curve's late-time mean.

Reconstructed from the lost /tmp/waterfall_dw_v10.py.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps

ps.apply()

d = np.load(os.path.join(ps.DATA, "waterfall_dw_sector", "dw_sector_weights_canonical_data.npz"))

T = 1800
# front (depth 0) to back
cases = [
    ("uniform_ER80", "uniform $E[R]$=80", '#1b9e77'),
    ("k1_ER4", "$k$=1 sym, $E[R_A]$=$E[R_B]$=4", '#d95f02'),
    ("k2_ER32", "$k$=2 sym, $E[R_A]$=$E[R_B]$=32", '#7570b3'),
    ("k2_ER4", "$k$=2 sym, $E[R_A]$=$E[R_B]$=4", '#e7298a'),
    ("k0_ER4", "$k$=0 iid, $E[R_A]$=4,  $E[R_B]$=4/3", '#66a61e'),
    ("k2_tightloose", "$k$=2 tight/loose, $E[R_A]$=1.2,  $E[R_B]$=9.0", '#e6ab02'),
]

fig = plt.figure(figsize=(0.5 * ps.TEXTWIDTH, 2.8))
ax = fig.add_axes([0.1, -0.02, 0.9, 1.04], projection='3d')

t = np.arange(T + 1)
for i, (key, label, color) in enumerate(cases):
    w = d[f"{key}_dw_weight"][:T + 1]
    late = w[1000:].mean()
    ax.plot([100, T], [i, i], [late, late], color='gray', ls='--', lw=0.6)
    ax.plot(t, np.full_like(t, i, dtype=float), w, color=color, lw=0.8, label=label)

ax.set_xlim(0, T)
ax.set_ylim(-0.5, len(cases) - 0.5)
ax.set_zlim(0, 1.1)
ax.set_xticks([0, 500, 1000, 1500])
ax.set_yticks(range(len(cases)))
ax.set_yticklabels([])
ax.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.tick_params(axis='both', labelsize=7.5, pad=0)
ax.set_xlabel('Compressed layer number', labelpad=-3, fontsize=8.5)
# the 3D z-label lands on the far side at this azimuth; draw it by hand
fig.text(0.03, 0.5, 'DW sector weight', rotation=90, ha='center', va='center',
         fontsize=8.5)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor((0.96, 0.96, 0.96, 1.0))
    a._axinfo['grid']['linewidth'] = 0.4
ax.view_init(elev=13, azim=-125)
ax.set_box_aspect((1.6, 1.0, 0.9), zoom=1.08)

ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.3), fontsize=6,
          handlelength=1.4, framealpha=0.9, borderpad=0.3, labelspacing=0.2)

ps.save(fig, "waterfall_dw_sector_Vpaper.png")
