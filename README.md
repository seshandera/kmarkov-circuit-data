# Data for "Non-Markovian Gates and Operator Spreading"

Raw simulation data (NumPy `.npz` files) for the figures in the paper. Each subdirectory corresponds to one figure.

## Repository structure

### `cover_time_over_N/`
**Figure: cover_time_over_N.png** — Cover time / N vs system size for k=1, k=2, k=3 Markov sequences.

- `k1_p95_*.npz` — k=1 with P(switch)=0.95, D_eff=9.5. Multiple runs.
- `k2_g99_r80_*.npz` — k=2 with p_good=0.99, p_rr=0.80, D_eff=89.9.
- `k3_g99_p95_e95_*.npz` — k=3 with p_good=0.99, p_pair=0.95, p_exit=0.95, D_eff=533.
- `summary_*.npz` — Metadata (config names, system sizes).

Each file contains: `system_sizes`, `means`, `stds`, `medians`, `covered_fracs`, and per-N arrays `cover_times_n{N}` with 500 trials each. System sizes: N=20 to N=1600.

### `fig_4panel_corr/`
**Figure: fig_4panel_corr.png** — 4-panel figure: operator heatmap, spacetime correlation baseline C(Δx,Δt), and correlation differences ΔC for composite and single-layer SWAP/WC alphabets.

- `fig_4panel_corr_data.npz` — N=512, fA=0.83, λ=−0.15, 5 realizations, seed=42.

Keys: `comp_heatmap` (panel a), `comp_C_iid` (panel b), `comp_dC` (panel c), `sl_dC` (panel d), plus parameters `N`, `fA`, `lam`, `max_dx`, `max_dt`, `comp_heatmap_tscr`.

### `fig_scrambling_N512/`
**Figure: fig_scrambling_N512.png** — Scrambling time t_scr/N vs fA for SWAP/WC composite alphabet at N=512.

- `iid_sweep_N512.npz` — i.i.d. (k=0) scrambling times across 8 values of fA.

Keys: `N`, `n_trials`, `frac_A_values`, `iid_fA{XX}_scr_times`, `det_fA{XX}_scr_times` (scrambling time arrays for each fA value). The k=1 and k=2 best-case data points are hardcoded in the plotting script.

### `PSWAP_heatmaps_correlations/`
**Figure: PSWAP_heatmaps_correlations_Vpaper.png** — 8-panel figure with operator density heatmaps and spacetime correlations for k=0, k=1, k=2 Markov sequences driving PSWAP (θ=π/300) circuits on N=30 qubits.

- `k0_ER4_multiseed_data.npz` — k=0 (i.i.d.) with E[R]=4, 5 seeds.
- `k1_ER4_multiseed_data.npz` — k=1 Markov with E[R]=4, 5 seeds.
- `k2_ER4_multiseed_data.npz` — k=2 Markov with E[R]=4, 5 seeds.

Each file contains per-seed arrays: `site_density_seed{i}`, `weight_norm_seed{i}`, `corr_seed{i}`, `layers_seed{i}`, plus the seed-averaged `corr_avg`.

### `waterfall_dw_sector/`
**Figure: waterfall_dw_sector_Vpaper.png** — Domain-wall sector weight waterfall plot showing frequency content of weight oscillations under different batching schemes.

- `dw_sector_weights_canonical_data.npz` — Z-sector and DW-sector weights vs time for all canonical cases (θ=π/300, N=30).

## Loading data

```python
import numpy as np
data = np.load("PSWAP_heatmaps_correlations/k0_ER4_multiseed_data.npz")
print(list(data.keys()))
```

## Paper

S. Shandera et al., "Non-Markovian Gates and Operator Spreading," submitted to Physical Review Letters (2026).
