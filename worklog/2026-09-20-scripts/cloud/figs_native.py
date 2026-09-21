"""Figures from the native-rate HAPwave sub-02 results."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, ORANGE, GREEN = '#0173B2', '#DE8F05', '#029E73'
INK, MUTED = '#1a1a1a', '#6b6b6b'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})
U = '/mnt/user-data/uploads/iEEGLAB_data/results_sub02'
OUT = '/mnt/user-data/outputs'

p = pd.read_csv(f'{U}/crp_paired.csv')
resp = p.z > 6

# ---- Figure 1: CRP explained variance, responsive vs not vs surrogate
fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.0), gridspec_kw={'width_ratios': [1.25, 1]})
ax = axes[0]
bins = np.linspace(0, 1, 30)
ax.hist(p.expl_s, bins=bins, color=ORANGE, alpha=.8, label=f'no stimulation (surrogate), n={len(p)}')
ax.hist(p.loc[resp, 'expl'], bins=bins, color=BLUE, alpha=.85,
        label=f'responsive pairs, n={int(resp.sum())}')
ax.set_xlabel('CRP explained variance per trial, 15-1000 ms')
ax.set_ylabel('site x contact pairs')
ax.set_title('Native rate, 2048 Hz, HAPwave sub-02', fontsize=10, loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper left')

ax = axes[1]
grp = [('responsive', p.loc[resp, 'delta']), ('non-responsive', p.loc[~resp, 'delta'])]
for i, (name, v) in enumerate(grp):
    ax.scatter(np.full(len(v), i) + np.random.default_rng(0).normal(0, .07, len(v)), v,
               s=3, alpha=.25, color=[BLUE, MUTED][i])
    ax.plot([i - .3, i + .3], [v.median()] * 2, color=INK, lw=2.5)
    ax.text(i, v.median() + .06, f'{v.median():+.2f}', ha='center', fontsize=9, color=INK)
ax.axhline(0, color=MUTED, lw=.8, ls='--')
ax.set_xticks([0, 1]); ax.set_xticklabels(['responsive', 'non-responsive'])
ax.set_ylabel('explained variance,\nreal minus surrogate (paired)')
ax.set_title('Paired against the same pair without stimulation', fontsize=9, loc='left')
fig.savefig(f'{OUT}/nat_fig1_crp.png')

# ---- Figure 2: ICA, reference effect versus seed-to-seed variability
m = np.load(f'{U}/ica_match.npz')
s = np.load(f'{U}/ica_summary.npz')
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.0))
ax = axes[0]
for v, c, lab in [(m['seed'], GREEN, 'same data, two random seeds'),
                  (m['ref'], BLUE, 'as recorded vs common average')]:
    ax.plot(np.arange(1, len(v) + 1), v, color=c, lw=2, label=f'{lab} (median {np.median(v):.2f})')
ax.set_xlabel('component, ranked'); ax.set_ylabel('matched |r| of time courses')
ax.set_ylim(0, 1.02)
ax.set_title('Reference vs random seed', fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8, loc='lower left')

ax = axes[1]
for i, (v, c) in enumerate([(s['f1'], BLUE), (s['f2'], ORANGE)]):
    ax.scatter(np.full(len(v), i) + np.random.default_rng(1).normal(0, .06, len(v)), v,
               s=8, alpha=.5, color=c)
    ax.plot([i - .25, i + .25], [np.median(v)] * 2, color=INK, lw=2.5)
    ax.text(i, np.median(v) + 1.5, f'{np.median(v):.0f}', ha='center', fontsize=9, color=INK)
ax.set_xticks([0, 1]); ax.set_xticklabels(['as recorded', 'after CAR'])
ax.set_yscale('log')
ax.set_ylabel('contacts >50% of peak (of 198)')
ax.set_title('Maps stay focal either way', fontsize=9.5, loc='left')
fig.savefig(f'{OUT}/nat_fig2_ica.png')

# ---- Figure 3: bad channels and the phase / high-frequency gap
b = np.load(f'{U}/badchan2.npz', allow_pickle=True)
clin = b['status'] != 'good'
fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.0))
ax = axes[0]
ax.scatter(b['sd'][~clin], b['nbrcorr'][~clin], s=14, color=BLUE, alpha=.7, label='clinician: good')
ax.scatter(b['sd'][clin], b['nbrcorr'][clin], s=26, color=ORANGE, alpha=.9, label='clinician: bad')
ax.set_xscale('log')
ax.set_xlabel('baseline SD, 1-200 Hz (uV)')
ax.set_ylabel('max correlation with a contact\non the same shaft')
ax.set_title('Clinician labels are not a noise criterion', fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8, loc='lower left')

ax = axes[1]
ph = pd.read_csv(f'{U}/phase_hfb.csv')
w = .34
for i, col in enumerate(['theta_alpha_itpc_post', 'beta_itpc_post']):
    for k, (cond, c) in enumerate([('real', BLUE), ('surrogate', ORANGE)]):
        v = ph.loc[ph.cond == cond, col]
        ax.bar(i + (k - .5) * w, v.median(), width=w, color=c, label=cond if i == 0 else None)
        ax.text(i + (k - .5) * w, v.median() + .03, f'{v.median():.2f}', ha='center', fontsize=8)
ax.set_xticks([0, 1]); ax.set_xticklabels(['4-12 Hz', '13-30 Hz'])
ax.set_ylim(0, 1.15)
ax.set_ylabel('intertrial phase clustering,\npeak 11-100 ms')
ax.set_title('Phase: 60 responsive pairs', fontsize=9.5, loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper right')
fig.savefig(f'{OUT}/nat_fig3_badchan_phase.png')
print('figures done')
