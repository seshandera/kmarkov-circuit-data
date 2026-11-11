"""PSWAP_heatmaps_correlations_Vpaper.png: 4 site-occupation heatmaps (seed 42)
over 4 seed-averaged late-time correlations C(dx, dt), PSWAP theta=pi/300, N=30.

Reconstructed from the lost /tmp/pswap_8panel_v3.py. Correlation panels use one
shared symmetric scale (vmax=0.71, read back from the published figure).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps

ps.apply()

CORR_VMAX = 0.71

cases = [
    # (multiseed file, title line(s), x-range of heatmap)
    ("uniform_ER80", "Uniform $E[R]$=80", 200),
    ("k0_ER4", "$k$=0 iid, $P(A)$=0.75\n$E[R_A]$=4,  $E[R_B]$=4/3", 800),
    ("k1_ER4", "$k$=1 sym\n$E[R_A]$=$E[R_B]$=4", 800),
    ("k2_tightloose", "$k$=2 tight/loose\n$E[R_A]$=1.2,  $E[R_B]$=9.0", 800),
]

fig, axes = plt.subplots(2, 4, figsize=(ps.TEXTWIDTH, 3.25),
                         gridspec_kw=dict(height_ratios=[1, 0.85]))

for j, (key, title, xmax) in enumerate(cases):
    ms = np.load(os.path.join(ps.PSWAP, f"{key}_multiseed_data.npz"))
    if key == "uniform_ER80":
        # the multiseed file for this case does not store per-seed densities
        single = np.load(os.path.join(ps.PSWAP, "uniform_ER80_data.npz"))
        dens = single['site_density']
    else:
        dens = ms['site_density_seed42']
    tws = ms['tws']
    tw, tw_bar = int(tws[0]), int(round(tws.mean()))
    Nq = int(ms['N'])
    max_dx, max_dt = int(ms['max_dx']), int(ms['max_dt'])

    ax = axes[0][j]
    him = ax.imshow(dens.T, aspect='auto', origin='lower', cmap='inferno',
                    interpolation='nearest', vmin=0, vmax=1,
                    extent=[0, dens.shape[0], 0, Nq])
    ax.axvline(tw, color='#7fb3bf', ls='--', lw=0.8)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, Nq)
    ax.set_yticks([0, 10, 20, 30])
    ax.set_xticks([0, 100, 200] if xmax == 200 else [0, 400, 800])
    ax.set_xticklabels(['0', '100', '200'] if xmax == 200 else ['0', '400', '800'])
    for lab, ha in zip(ax.get_xticklabels(), ['left', 'center', 'right']):
        lab.set_ha(ha)
    ax.set_xlabel('Compressed layer number', labelpad=1, fontsize=7.5)
    ax.set_title(f"{title}\n$t_w={tw}$,  $\\bar t_w={tw_bar}$", fontsize=7.5,
                 linespacing=1.15)
    if j == 0:
        ax.set_ylabel('Site')
    else:
        ax.set_yticklabels([])

    ax = axes[1][j]
    C = ms['corr_avg'][max_dt:]   # dt >= 0
    cim = ax.imshow(C, aspect='auto', origin='lower', cmap='RdBu_r',
                    interpolation='nearest', vmin=-CORR_VMAX, vmax=CORR_VMAX,
                    extent=[-max_dx - 0.5, max_dx + 0.5, -0.5, max_dt + 0.5])
    ax.set_xlim(-max_dx, max_dx)
    ax.set_ylim(0, max_dt)
    ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    ax.set_xticklabels(['\u221215', '', '\u22125', '0', '5', '', '15'])
    ax.get_xticklabels()[0].set_ha('left')
    ax.get_xticklabels()[-1].set_ha('right')
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_xlabel(r'$\Delta x$', labelpad=1)
    if j == 0:
        ax.set_ylabel(r'$\Delta t$')
    else:
        ax.set_yticklabels([])
    ax.text(0.5, -0.42, f'({"abcd"[j]})', transform=ax.transAxes,
            ha='center', va='top', fontsize=9, fontweight='bold')

fig.subplots_adjust(left=0.055, right=0.91, top=0.86, bottom=0.14,
                    wspace=0.14, hspace=0.6)
for row, im, ticks in [(0, him, [0, 0.2, 0.4, 0.6, 0.8]),
                       (1, cim, [-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])]:
    b = axes[row][3].get_position()
    cax = fig.add_axes([0.925, b.y0, 0.012, b.height])
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    cb.ax.tick_params(labelsize=7.5)

ps.save(fig, "PSWAP_heatmaps_correlations_Vpaper.png")
