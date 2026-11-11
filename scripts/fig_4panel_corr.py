"""fig_4panel_corr.png: heatmap + C baseline + Delta C composite + Delta C single-layer.

Plot-only version of quantum_batching/examples/_fig_4panel_corr.py, reading the
saved simulation output (fig_4panel_corr_data.npz) instead of re-simulating.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps

ps.apply()

d = np.load(os.path.join(ps.DATA, "fig_4panel_corr", "fig_4panel_corr_data.npz"))
N = int(d['N'])
fA = float(d['fA'])
lam = float(d['lam'])
max_dx = int(d['max_dx'])
max_dt = int(d['max_dt'])
comp_heatmap = d['comp_heatmap']
comp_heatmap_tscr = int(d['comp_heatmap_tscr'])

fig, axes = plt.subplots(2, 2, figsize=(0.5 * ps.TEXTWIDTH, 3.35))


def corr_panel(ax, data, vmin_floor, title, xlabel):
    nodc = data.copy()
    nodc[max_dt, max_dx] = 0
    vmax = max(np.percentile(np.abs(nodc), 98), vmin_floor)
    im = ax.imshow(data.T, aspect='auto', origin='lower', cmap='RdBu_r',
                   interpolation='nearest',
                   extent=[-max_dt, max_dt, -max_dx, max_dx],
                   vmin=-vmax, vmax=vmax)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03, fraction=0.08)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.ax.yaxis.get_offset_text().set_fontsize(6.5)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\Delta x$ (sites)', labelpad=1)
    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.set_xticklabels(['\u221240', '', '0', '', '40'])
    ax.set_yticks([-40, -20, 0, 20, 40])
    ax.axhline(0, color='gray', ls='-', lw=0.3)
    ax.axvline(0, color='gray', ls='-', lw=0.3)
    ax.plot([-max_dt, max_dt], [-max_dt, max_dt], 'k--', lw=0.4, alpha=0.3)
    return cb


# (a) operator heatmap
ax = axes[0][0]
max_t_show = int(1.5 * N)
hm_shifted = np.roll(comp_heatmap[:max_t_show], N // 2, axis=1)
ax.imshow(hm_shifted.T, aspect='auto', origin='lower', cmap='hot',
          interpolation='nearest', extent=[0, max_t_show / N, -N // 2, N // 2])
l1 = ax.axvline(0.269, color='#0066cc', ls=':', lw=1.6, alpha=0.9)
l2 = ax.axvline(comp_heatmap_tscr / N, color='#228B22', ls='--', lw=1.6, alpha=0.9)
ax.legend([l1, l2], ['full Clifford', f'i.i.d. $f_A$={fA}'], loc='upper right',
          fontsize=5.5, framealpha=0.85, handlelength=1.5, borderpad=0.25,
          labelspacing=0.2, handletextpad=0.4)
ax.set_title(f'(a) Heatmap: k=0 double\nlayer, $f_A$={fA}')
ax.set_xlabel('t / N (letters)')
ax.set_ylabel('site (centered)', labelpad=1)
ax.set_xticks([0, 0.5, 1.0, 1.5])
ax.set_yticks([-200, 0, 200])
ax.set_ylim(-N // 2, N // 2)

corr_panel(axes[0][1], d['comp_C_iid'], 0.01,
           '(b) $C(\\Delta x, \\Delta t)$: k=0 double\n'
           f'layer, $f_A$={fA}, late-time', r'$\Delta t$ (letters)')
cb_c = corr_panel(axes[1][0], d['comp_dC'], 0.005,
           '(c) $\\Delta C$: k=1 double layer\n'
           f'$f_A$={fA}, $\\lambda$={lam}, late-time', r'$\Delta t$ (letters)')
corr_panel(axes[1][1], d['sl_dC'], 0.005,
           '(d) $\\Delta C$: k=1 single layer\n'
           f'$f_A$={fA}, $\\lambda$={lam}, late-time', r'$\Delta t$ (sub-layers)')

for a in axes.flat:
    a.title.set_fontsize(7.5)
    a.xaxis.label.set_fontsize(8)
    a.yaxis.label.set_fontsize(8)
    a.tick_params(labelsize=7)

fig.tight_layout(h_pad=0.6, w_pad=0.3)
# widen panel (a) to span panel (c) plus its colorbar and tick labels
fig.canvas.draw()
pa, pc = axes[0][0].get_position(), axes[1][0].get_position()
to_fig = fig.transFigure.inverted()
cb_right = cb_c.ax.get_tightbbox().transformed(to_fig).x1
b_label_left = axes[0][1].yaxis.label.get_window_extent().transformed(to_fig).x0
right = min(cb_right, b_label_left - 0.07)  # leave clear space before (b)'s y-label
axes[0][0].set_position([pc.x0, pa.y0, right - pc.x0, pa.height])
ps.save(fig, "fig_4panel_corr.png")
