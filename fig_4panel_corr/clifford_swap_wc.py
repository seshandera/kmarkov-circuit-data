"""Standalone SWAP/WC Clifford simulator for the self-averaging check.

Extracted from the authors' research code (quantum_batching/examples:
scrambling_study.py and clifford_alphabet_study.py) so that
reproduce_selfaverage_clifford.py runs from this repository alone. The gate
enumeration order and the random draws are unchanged, so a given seed gives
the same heatmap as the original code.

Pauli tracking (Gottesman-Knill): a Pauli string on N qubits is stored as two
binary arrays x[N], z[N], with site i carrying I=(0,0), X=(1,0), Z=(0,1),
Y=(1,1). A two-qubit Clifford gate on sites (a, b) acts as a 4x4 binary
symplectic matrix, an element of Sp(4, F_2), on (x_a, x_b, z_a, z_b).
"""
import numpy as np


def enumerate_sp4_f2():
    """All 720 elements of Sp(4, F_2), in the order used by the paper's code.

    Tests every binary 4x4 matrix M for M^T Omega M = Omega (mod 2), where
    Omega pairs x_a with z_a and x_b with z_b.
    """
    omega = np.zeros((4, 4), dtype=np.uint8)
    omega[0, 2] = omega[1, 3] = omega[2, 0] = omega[3, 1] = 1
    bits = np.arange(1 << 16, dtype=np.uint32)
    all_m = np.zeros((1 << 16, 4, 4), dtype=np.uint8)
    for idx in range(16):
        i, j = divmod(idx, 4)
        all_m[:, i, j] = (bits >> idx) & 1
    mt_o = np.einsum('bji,jk->bik', all_m, omega) % 2
    prod = np.einsum('bij,bjk->bik', mt_o, all_m) % 2
    sp4 = all_m[np.all(prod == omega[None], axis=(1, 2))]
    assert sp4.shape[0] == 720
    return sp4


SP4_ELEMENTS = enumerate_sp4_f2()


def weight_changing_gates():
    """Indices of the 648 weight-changing (WC) gates.

    A gate is weight-preserving if every single-site Pauli (XI, IX, ZI, IZ,
    YI, IY) stays on one site after one application; the 72 such gates
    permute Pauli axes and positions. WC is the complement.
    """
    single_site = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                            [0, 0, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]],
                           dtype=np.uint8)
    wc = []
    for idx, m in enumerate(SP4_ELEMENTS):
        out = (single_site @ m.T) % 2
        weights = (out[:, 0] | out[:, 2]).astype(int) + (out[:, 1] | out[:, 3])
        if np.any(weights != 1):
            wc.append(idx)
    wc = np.array(wc, dtype=int)
    assert len(wc) == 648
    return wc


WC_GATES = weight_changing_gates()


def pair_indices(n_qubits, letter):
    """Sites (a, b) of the even (letter 0) or odd (letter 1) brickwork layer on a ring."""
    if letter == 0:
        a = np.arange(0, n_qubits - 1, 2)
        return a, a + 1
    a = np.arange(1, n_qubits, 2)
    return a, (a + 1) % n_qubits


def apply_swap_layer(x, z, n_qubits, letter):
    a, b = pair_indices(n_qubits, letter)
    x[a], x[b] = x[b].copy(), x[a].copy()
    z[a], z[b] = z[b].copy(), z[a].copy()


def apply_wc_layer(x, z, n_qubits, letter, rng):
    """One brickwork layer of WC gates drawn uniformly at random."""
    a, b = pair_indices(n_qubits, letter)
    v = np.stack([x[a], x[b], z[a], z[b]], axis=1)
    gates = SP4_ELEMENTS[WC_GATES[rng.integers(0, len(WC_GATES), size=len(a))]]
    w = np.einsum('pij,pj->pi', gates, v) % 2
    x[a], x[b], z[a], z[b] = w[:, 0], w[:, 1], w[:, 2], w[:, 3]


def simulate_swap_wc(n_qubits, max_letters, sequence, rng, threshold=0.5,
                     record_heatmap=False):
    """Evolve Z on site 0 under the SWAP/WC alphabet.

    Letter A (0): SWAP on even then odd pairs (a shift by one site).
    Letter B (1): random WC gates on even then odd pairs.

    Returns a dict with the weight at each letter, the scrambling time
    (first letter at which the weight reaches threshold * N), whether it was
    reached, and optionally the (max_letters, N) occupancy heatmap.
    """
    target = max(1, int(np.ceil(threshold * n_qubits)))
    x = np.zeros(n_qubits, dtype=np.uint8)
    z = np.zeros(n_qubits, dtype=np.uint8)
    z[0] = 1
    weights = np.zeros(max_letters, dtype=np.int32)
    heatmap = np.zeros((max_letters, n_qubits), dtype=np.uint8) if record_heatmap else None
    scr_time, hit = max_letters, False
    for t in range(max_letters):
        if int(sequence[t]) == 0:
            apply_swap_layer(x, z, n_qubits, 0)
            apply_swap_layer(x, z, n_qubits, 1)
        else:
            apply_wc_layer(x, z, n_qubits, 0, rng)
            apply_wc_layer(x, z, n_qubits, 1, rng)
        occ = (x | z) > 0
        weights[t] = int(occ.sum())
        if heatmap is not None:
            heatmap[t] = occ
        if not hit and weights[t] >= target:
            scr_time, hit = t + 1, True
    return {'weights': weights, 'scr_time': scr_time, 'hit': hit, 'heatmap': heatmap}
