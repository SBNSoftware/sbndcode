"""Shared matplotlib style for SBND NuIntNN paper figures.

Usage in each notebook (first cell after imports):
    import sys
    sys.path.insert(0, '/exp/sbnd/app/users/svidales/larsoft_develop/srcs/sbndcode/sbndcode/NuIntNNOpReco/2-cnn-training-notebooks/')
    from plot_style import *
"""
import matplotlib as mpl

# ── rcParams ──────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    # Fonts
    'font.family':      'DejaVu Sans',
    'axes.labelsize':   14,
    'axes.titlesize':   14,
    'axes.titleweight': 'bold',
    'xtick.labelsize':  12,
    'ytick.labelsize':  12,
    'legend.fontsize':  12,
    'legend.framealpha': 1.0,
    'legend.edgecolor': '0.7',
    # Lines
    'lines.linewidth':  2.0,
    # Grid (subtle, off by default — enable per-axes with ax.grid(True))
    'axes.grid':       False,
    'grid.alpha':      0.3,
    'grid.linewidth':  0.6,
    # Minor ticks on by default
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top':       True,
    'ytick.right':     True,
    # Figure / saving
    'figure.dpi':      150,
    'savefig.dpi':     300,
    'savefig.bbox':    'tight',
})

# ── Color palette (Paul Tol "vibrant") ───────────────────────────────────────
C_PRIMARY   = '#0077BB'    # Tol vibrant blue
C_SECONDARY = '#EE7733'    # Tol vibrant orange — high contrast on blue backgrounds
C_TERTIARY  = '#009988'    # Tol vibrant teal
C_ALERT     = '#CC3311'    # Tol vibrant red

COLORS = {'LSTM': C_PRIMARY, 'Transformer': C_SECONDARY}   # used in TimeCoord

# ── 2D histogram style ────────────────────────────────────────────────────────
CMAP_2D  = 'PuBu'    # white → light blue → dark blue
SAVE_DPI = 300

# ── Font size constants ───────────────────────────────────────────────────────
FONT_LABEL = 14
FONT_TICK  = 12
FONT_STAT  = 11

# ── Stat box style ────────────────────────────────────────────────────────────
STAT_BOX = dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='0.6', alpha=0.9)

# ── Experiment label helpers ──────────────────────────────────────────────────
EXPERIMENT_LABEL = 'SBND Work In Progress'

def add_experiment_label(fig, label=EXPERIMENT_LABEL, fontsize=11):
    """Add 'SBND Work In Progress' watermark to bottom-right of a figure."""
    fig.text(0.99, 0.01, label, fontsize=fontsize, style='italic',
             va='bottom', ha='right', color='0.45')
