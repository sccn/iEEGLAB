"""Two minimal figures for the email."""
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
U = '/mnt/user-data/uploads/iEEGLAB_data/results_sub02'
OUT = '/mnt/user-data/outputs'

# ---- Figure A: ICA, reference effect against seed-to-seed variability
m = np.load(f'{U}/ica_match.npz')
fig, ax = plt.subplots(figsize=(5.0, 3.2))
for v, c, lab in [(m['seed'], GREEN, 'same data, two random seeds'),
                  (m['ref'], BLUE, 'before vs after common average')]:
    ax.plot(np.arange(1, len(v) + 1), v, color=c, lw=2.2, label=lab)
ax.set_xlabel('component, ranked by agreement')
ax.set_ylabel('matched |r| of component time courses')
ax.set_ylim(0, 1.03)
ax.set_title('Re-referencing perturbs ICA no more than the random seed',
             fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8.5, loc='lower left')
ax.text(.98, .93, 'HAPwave sub-02, 198 contacts', transform=ax.transAxes,
        ha='right', fontsize=8, color=MUTED)
fig.savefig(f'{OUT}/email_figA_ica.png')

# ---- Figure B: clinician-bad contacts are bimodal, not simply noisy
b = np.load(f'{U}/badchan2.npz', allow_pickle=True)
sd, nb, status, labels = b['sd'], b['nbrcorr'], b['status'], b['labels']
clin = status != 'good'
med_good = np.median(sd[~clin])
quiet = clin & (sd < med_good)
loud = clin & (sd >= med_good)
fig, ax = plt.subplots(figsize=(5.4, 3.2))
ax.scatter(sd[~clin], nb[~clin], s=16, color=BLUE, alpha=.65, label=f'clinician: good (n={int((~clin).sum())})')
ax.scatter(sd[quiet], nb[quiet], s=34, color=ORANGE, alpha=.95, marker='v',
           label=f'clinician: bad, low amplitude (n={int(quiet.sum())})')
ax.scatter(sd[loud], nb[loud], s=34, color=ORANGE, alpha=.95, marker='^',
           label=f'clinician: bad, high amplitude (n={int(loud.sum())})')
ax.axvline(med_good, color=MUTED, lw=.9, ls='--')
ax.text(med_good * 1.07, .12, 'median of good contacts', fontsize=7.5, color=MUTED, rotation=90)
ax.set_xscale('log')
ax.set_xlabel('baseline SD, 1-200 Hz (uV)')
ax.set_ylabel('max correlation within shaft')
ax.set_title('A variance criterion catches 12 of 14 noisy bad contacts and 0 of 25 quiet ones',
             fontsize=8.5, loc='left')
ax.legend(frameon=False, fontsize=7.5, loc='lower left')
fig.savefig(f'{OUT}/email_figB_badchannels.png')
print('quiet', int(quiet.sum()), 'loud', int(loud.sum()),
      'z>5 among quiet', int(((b['z'] > 5) & quiet).sum()),
      'z>5 among loud', int(((b['z'] > 5) & loud).sum()))
