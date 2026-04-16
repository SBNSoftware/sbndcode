"""Shared matplotlib style for SBND NuIntNN paper figures.

Usage in each notebook (first cell after imports):
    import sys
    sys.path.insert(0, '/exp/sbnd/app/users/svidales/larsoft_develop/srcs/sbndcode/sbndcode/NuIntNNOpReco/2-cnn-training-notebooks/')
    from plot_style import *
"""
import matplotlib as mpl

# ── rcParams ──────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    'axes.labelsize':  14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'lines.linewidth': 2.0,
    'figure.dpi':      150,
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

# ── Dataset / experiment labels (available but not added to figures) ──────────
MC_LABEL         = ''       # dataset watermark
EXPERIMENT_LABEL = 'SBND Work In Progress'         # experiment label
# Usage (uncomment to add to a figure):
#   fig.text(0.01, 0.01, MC_LABEL, fontsize=10, style='italic', va='bottom')
#   fig.text(0.99, 0.01, EXPERIMENT_LABEL, fontsize=10, style='italic', va='bottom', ha='right')
