"""Shared print-size style for the paper figures.

Figures are drawn at the width they occupy on the page (revtex4-2 twocolumn,
\\textwidth ~ 7.06 in), so the font sizes below are the printed sizes.
"""
import os
import matplotlib as mpl

TEXTWIDTH = 7.06  # inches

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = REPO
PSWAP = os.path.join(REPO, "PSWAP_heatmaps_correlations")
OUT = os.environ.get("FIG_OUT", os.path.join(REPO, "figures"))

DPI = 600


def apply():
    mpl.rcParams.update({
        'font.size': 8.5,
        'axes.labelsize': 9.5,
        'axes.titlesize': 9,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8,
        'axes.linewidth': 0.7,
        'xtick.major.width': 0.7,
        'ytick.major.width': 0.7,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.minor.size': 1.8,
        'ytick.minor.size': 1.8,
        'lines.linewidth': 1.2,
        'lines.markersize': 4,
        'legend.handlelength': 1.6,
        'legend.borderpad': 0.35,
        'legend.labelspacing': 0.3,
        'savefig.dpi': DPI,
        'savefig.pad_inches': 0.02,
    })


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    print("Saved", path)
