# Data and figure scripts for "Operator dynamics in k-Markov random circuits"

Raw simulation data (NumPy `.npz` files) and the plotting scripts for every figure in

> U. Akhouri, P.-J. Huang, E. Rose, S. Shandera, "Operator dynamics in k-Markov random circuits,"
> *Phys. Rev. Research* **8**, 043016 (2026). DOI: [10.1103/p2cj-1kkl](https://doi.org/10.1103/p2cj-1kkl) ·
> arXiv: [2603.18217](https://arxiv.org/abs/2603.18217)

The paper was first submitted under the title "Non-Markovian Gates and Operator Spreading." This repository
now matches the published version.

## Reproducing the figures

```bash
pip install numpy matplotlib
cd scripts
python fig_cover_time.py        # cover_time_over_N.png
python fig_scrambling.py        # fig_scrambling_N512.png
python fig_4panel_corr.py       # fig_4panel_corr.png
python fig_pswap_8panel.py      # PSWAP_heatmaps_correlations_Vpaper.png
python fig_waterfall_dw.py      # waterfall_dw_sector_Vpaper.png
```

Figures are written to `figures/` (set `FIG_OUT` to write elsewhere). `scripts/paperstyle.py` sets the
printed font sizes and resolves all data paths relative to the repository. The scripts read only the saved
data here; nothing is re-simulated. Different matplotlib versions can change tick-label formatting and
cropping slightly, but not the plotted data.

## Repository structure

### `cover_time_over_N/`
**Figure: `cover_time_over_N.png`.** Cover time / N against system size N for k=1, k=2 and k=3 Markov
sequences on the SWAP ring, with the k=0 uniform walk as a reference line. Every file stores 500 trials per N.

| Curve | Files (merged; the later file wins at a shared N) | N range |
|---|---|---|
| k=1, P(switch)=0.95, D_eff=9.5 | `k1_p95_20260207_120350.npz`, `k1_p95_20260207_111247.npz` | 20–8000 |
| k=2, p_good=0.97, p_rr=0.43, D_eff=23 | `k2_g97_r43_20260209_142810.npz`, `k2_g97_r43_extended.npz` | 20–8000 |
| k=2, p_good=0.99, p_rr=0.80, D_eff=89 | `k2_g99_r80_20260207_162653.npz`, `k2_g99_r80_20260206_182837.npz` | 20–8000 |
| k=3, p_good=0.99, p_pair=0.95, p_exit=0.95, D_eff=533 | `k3_g99_p95_e95_20260207_162653.npz`, `k3_g99_p95_e95_20260207_070758.npz` | 20–8000 |

Each file contains `system_sizes`, `means`, `stds`, `medians`, `covered_fracs`, and per-N arrays
`cover_times_n{N}`. Open markers in the figure are sizes where fewer than 90% of trials covered the ring,
so the plotted mean there is a lower bound. `k1_p95_20260207_144216.npz`, `k1_p95_20260207_162653.npz`
and `summary_20260207_162653.npz` are additional runs and metadata not used in the published figure.

### `fig_scrambling_N512/`
**Figure: `fig_scrambling_N512.png`.** Scrambling time t_scr/N against the SWAP fraction f_A for the
SWAP/WC composite alphabet at N=512.

- `iid_sweep_N512.npz`: i.i.d. (k=0) scrambling times at 8 values of f_A. Keys: `N`, `n_trials`,
  `frac_A_values`, `iid_fA{XX}_scr_times`, `det_fA{XX}_scr_times`. The k=1 and k=2 best-case points are
  set in `scripts/fig_scrambling.py`.

### `fig_4panel_corr/`
**Figure: `fig_4panel_corr.png`.** Operator heatmap, late-time spacetime correlation C(Δx,Δt), and
correlation differences ΔC for the composite and single-layer SWAP/WC alphabets.

- `fig_4panel_corr_data.npz`: N=512, f_A=0.83, λ=−0.15, 5 realizations, seed 42. Keys: `comp_heatmap`
  (panel a), `comp_C_iid` (b), `comp_dC` (c), `sl_dC` (d), plus `N`, `fA`, `lam`, `max_dx`, `max_dt`,
  `comp_heatmap_tscr`.

### `PSWAP_heatmaps_correlations/`
**Figure: `PSWAP_heatmaps_correlations_Vpaper.png`.** Eight panels: operator-density heatmaps (top) and
late-time spacetime correlations (bottom) for PSWAP circuits (θ=π/300) on N=30 qubits, driven by four
sequence rules. Each multiseed file holds 5 seeds (42–46): per-seed arrays such as `corr_seed{i}`, the
seed-averaged correlation `corr_avg`, and the per-seed equilibration times `tws`.

| Panel | File |
|---|---|
| Uniform, E[R]=80 | `uniform_ER80_multiseed_data.npz`; heatmap from `uniform_ER80_data.npz` (single run) |
| k=0 i.i.d., P(A)=0.75, E[R_A]=4, E[R_B]=4/3 | `k0_ER4_multiseed_data.npz` |
| k=1 symmetric, E[R_A]=E[R_B]=4 | `k1_ER4_multiseed_data.npz` |
| k=2 tight/loose, E[R_A]=1.2, E[R_B]=9.0 | `k2_tightloose_multiseed_data.npz` |

`k2_ER4_multiseed_data.npz` (k=2 symmetric, E[R]=4) appeared in the originally submitted version of this
figure and is kept for reference.

### `waterfall_dw_sector/`
**Figure: `waterfall_dw_sector_Vpaper.png`.** Domain-wall sector weight against time for six sequence
rules, showing the frequency content of the weight oscillations about the common late-time plateau.

- `dw_sector_weights_canonical_data.npz`: Z-sector and DW-sector weights against time for all cases
  (θ=π/300, N=30).

## Loading data

```python
import numpy as np
data = np.load("PSWAP_heatmaps_correlations/k0_ER4_multiseed_data.npz")
print(list(data.keys()))
```

## Citation

```bibtex
@article{Akhouri2026kMarkov,
  author  = {Akhouri, Unnati and Huang, Pei-Jun and Rose, Elliott and Shandera, Sarah},
  title   = {Operator dynamics in k-Markov random circuits},
  journal = {Phys. Rev. Research},
  volume  = {8},
  pages   = {043016},
  year    = {2026},
  doi     = {10.1103/p2cj-1kkl},
  eprint  = {2603.18217},
  archivePrefix = {arXiv}
}
```
