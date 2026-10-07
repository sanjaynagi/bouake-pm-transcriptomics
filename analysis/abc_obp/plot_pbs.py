"""PBS across the ABCH1 region, styled to match Figure 5B.

Reads pbs_2R.csv (written by pbs.py on the cluster) and draws one trace per An. coluzzii
cohort over the same x-range as Figure 5B, with the ABCH1 span marked. The lower axis widens
to the full computed window so the peak can be read against its local background.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
from bouake.paths import REPO, FONT_DIR

STAT_GREY, INK, INK_SOFT = "#5a5e66", "#2b2d31", "#8d99ae"
ABCH1 = (24_825_285, 24_848_323)
XLIM = (24.30e6, 25.40e6)

for f in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(f))
mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.7,
    "axes.edgecolor": "#4a4a4a", "svg.fonttype": "none", "figure.dpi": 150,
})

df = pd.read_csv(f"{REPO}/results/abc_obp/pbs/pbs_2R.csv")
cohorts = [c for c in df.columns if c not in ("position", "n_segregating")]

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.2, 4.2),
                               gridspec_kw=dict(height_ratios=[1, 1], hspace=0.45))

# The Figure 5B view
for c in cohorts:
    ax0.plot(df.position, df[c], lw=0.8, color=STAT_GREY, alpha=0.85)
ax0.axvspan(*ABCH1, color="#e8e9ea", zorder=0)
ax0.annotate("ABCH1", (np.mean(ABCH1), 0.97), xycoords=("data", "axes fraction"),
             ha="center", va="top", fontsize=7, color=INK)
ax0.set_xlim(*XLIM)
sel = df.position.between(*XLIM)
ax0.set_ylim(bottom=min(0, np.nanmin(df.loc[sel, cohorts].values)))
ax0.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v/1e6:.1f}"))
ax0.annotate(f"{len(cohorts)} cohorts", xy=(0.995, 0.95), xycoords="axes fraction",
             ha="right", va="top", fontsize=6.5, color=INK_SOFT)

# The full computed window, for local background
for c in cohorts:
    ax1.plot(df.position, df[c], lw=0.45, color=STAT_GREY, alpha=0.6)
ax1.axvspan(*ABCH1, color="#e8e9ea", zorder=0)
ax1.set_xlim(df.position.min(), df.position.max())
ax1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v/1e6:.1f}"))

for ax in (ax0, ax1):
    ax.set_ylabel("PBS")
    ax.set_xlabel("Chromosome arm 2R position (Mb)", labelpad=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#eceded", lw=0.5)
    ax.set_axisbelow(True)

fig.savefig(f"{REPO}/figures_ms/pbs_abch1.png", dpi=380, bbox_inches="tight")
print("wrote figures_ms/pbs_abch1.png")
