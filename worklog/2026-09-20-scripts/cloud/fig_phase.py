import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BLUE, ORANGE, GREEN = '#0173B2', '#DE8F05', '#029E73'
INK, MUTED = '#1a1a1a', '#6b6b6b'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                     'savefig.dpi': 220, 'savefig.bbox': 'tight'})
U='/mnt/user-data/uploads/iEEGLAB_data/results_sub02'; OUT='/mnt/user-data/outputs'
d=pd.read_csv(f'{U}/itpc_fixed.csv')
g=d.groupby(['cond','freq']).median(numeric_only=True).reset_index()
r=g[g.cond=='real']; s=g[g.cond=='surrogate']
fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.1), sharey=True)
ax=axes[0]
ax.plot(r.freq, r.evoked_post, color=BLUE, lw=2.2, marker='o', ms=4, label='after the pulse')
ax.plot(r.freq, r.evoked_pre, color=MUTED, lw=1.8, ls='--', marker='o', ms=3, label='before the pulse')
ax.plot(s.freq, s.evoked_post, color=ORANGE, lw=1.8, marker='s', ms=3.5, label='surrogate, after')
ax.set_xscale('log'); ax.set_xticks([4,8,16,32,80]); ax.set_xticklabels([4,8,16,32,80])
ax.set_xlabel('frequency (Hz)'); ax.set_ylabel('intertrial phase clustering')
ax.set_title('All trials', fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8)
ax=axes[1]
ax.plot(r.freq, r.induced_post, color=BLUE, lw=2.2, marker='o', ms=4, label='after the pulse')
ax.plot(r.freq, r.induced_pre, color=MUTED, lw=1.8, ls='--', marker='o', ms=3, label='before the pulse')
ax.plot(s.freq, s.induced_post, color=ORANGE, lw=1.8, marker='s', ms=3.5, label='surrogate, after')
ax.set_xscale('log'); ax.set_xticks([4,8,16,32,80]); ax.set_xticklabels([4,8,16,32,80])
ax.set_xlabel('frequency (Hz)')
ax.set_title('Evoked response subtracted', fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8)
ax.set_ylim(0, 1.05)
fig.suptitle('Phase locking is the evoked response: nothing survives its removal  (48 pairs, 10 trials each)',
             fontsize=9.5, x=.02, ha='left')
fig.savefig(f'{OUT}/email_figC_phase.png')
print('ok')
