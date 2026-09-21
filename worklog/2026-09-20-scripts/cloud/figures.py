"""Figures for the Monday meeting (late window, artifact-screened)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from ccep_core import load, epoch, blank_and_baseline, ref_car_persite_good, run_crp, BAD
from ica_test import run_ica, spread

BLUE, ORANGE, GREEN = '#0173B2', '#DE8F05', '#029E73'
INK, MUTED = '#1a1a1a', '#6b6b6b'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})
OUT = '/mnt/user-data/outputs'
WLO, WHI = 0.1, 1.0

df = pd.read_csv('crp_clean.csv')
d = df
real = d[d.cond == 'real'].expl.values
null = d[d.cond == 'surrogate'].expl.values
thr = np.quantile(null, 0.95)

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(5.4, 3.1))
bins = np.linspace(0, 1, 26)
ax.hist(null, bins=bins, color=ORANGE, alpha=.8, label=f'no stimulation (surrogate), n={len(null)}')
ax.hist(real, bins=bins, color=BLUE, alpha=.8, label=f'real stimulation, n={len(real)}')
ax.axvline(thr, color=INK, lw=1, ls='--')
ax.annotate('95th pct of surrogate', xy=(thr, ax.get_ylim()[1] * .55),
            xytext=(thr - .42, ax.get_ylim()[1] * .55), fontsize=8, color=INK, va='center',
            ha='left', arrowprops=dict(arrowstyle='->', color=INK, lw=.8))
ax.set_xlabel(f'CRP explained variance per trial, {WLO*1e3:.0f}-{WHI*1e3:.0f} ms')
ax.set_ylabel('site x contact pairs')
ax.set_title('Fraction of each trial captured by the rank-1 canonical shape',
             fontsize=10, loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper left')
fig.savefig(f'{OUT}/fig1_crp_explained_variance.png')

# ---------------------------------------------------------------- Figure 2
data, labels, srate, ev = load()
onsets, sites_all = ev.onset.values, ev.electrical_stimulation_site.values
ep, t, keep = epoch(data, srate, onsets)
sites = sites_all[keep]
X = ref_car_persite_good(blank_and_baseline(ep, t), labels, sites)
win = (t >= WLO) & (t <= WHI)
tw = t[win]

dr = d[d.cond == 'real'].sort_values('expl', ascending=False)
picks = [dr.iloc[3], dr.iloc[len(dr) // 2]]
fig, axes = plt.subplots(1, 2, figsize=(7.8, 2.9))
for ax, row in zip(axes, picks):
    tr = np.where(sites == row.site)[0]
    ci = labels.index(row.chan)
    V = X[ci][win][:, tr]
    r = run_crp(V, tw)
    ax.plot(tw * 1e3, V, color=MUTED, lw=.6, alpha=.55)
    n = len(r['C'])
    ax.plot(tw[:n] * 1e3, r['C'] * np.median(r['al']), color=BLUE, lw=2.2,
            label='canonical shape')
    ax.set_title(f"{row.site} to {row.chan}, {int(row.n)} trials, expl. var {row.expl:.2f}",
                 fontsize=9, loc='left')
    ax.set_xlabel('time (ms)')
axes[0].set_ylabel('uV')
axes[0].legend(frameon=False, fontsize=8)
fig.suptitle('Single trials and the rank-1 fit: high (left) and median (right) pair',
             fontsize=10, x=.02, ha='left')
fig.savefig(f'{OUT}/fig2_crp_example_fits.png')

# ---------------------------------------------------------------- Figure 3
good = [i for i, l in enumerate(labels) if l not in BAD]
labels = [labels[i] for i in good]
b, a = butter(4, 0.5 / (srate / 2), 'high')
Xc = filtfilt(b, a, data[good], axis=1)
N = len(labels)
S_raw, A_raw, _ = run_ica(Xc, N, 0)
S_car, A_car, _ = run_ica(Xc - Xc.mean(0, keepdims=True), N - 1, 0)
p_raw, p_car = spread(A_raw), spread(A_car)
jit = np.random.default_rng(0).normal(0, .04, N)

fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.1), gridspec_kw={'width_ratios': [1, 1.5]})
ax = axes[0]
for i, (p, c) in enumerate([(p_raw, BLUE), (p_car, ORANGE)]):
    ax.scatter(np.full(len(p), i) + jit[:len(p)], p, color=c, s=20, alpha=.85)
    ax.plot([i - .22, i + .22], [np.median(p)] * 2, color=INK, lw=2)
    ax.text(i, np.median(p) + .35, f'{np.median(p):.1f}', ha='center', fontsize=8, color=INK)
ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1])
ax.set_xticklabels(['as recorded', 'after CAR'])
ax.set_ylabel('participation ratio\n(1 = one contact, 13 = uniform)')
ax.set_title('Apparent spread of IC maps', fontsize=10, loc='left')

ax = axes[1]
order = np.argsort(-np.abs(A_raw).max(0))[:3]
w = .26
for k, ic in enumerate(order):
    v = A_raw[:, ic] / np.abs(A_raw[:, ic]).max()
    ax.barh(np.arange(N) + (k - 1) * w, v, height=w, color=[BLUE, ORANGE, GREEN][k],
            label=f'IC{ic + 1}')
ax.set_yticks(range(N)); ax.set_yticklabels(labels, fontsize=7)
ax.invert_yaxis(); ax.axvline(0, color=MUTED, lw=.8)
ax.set_xlabel('normalised weight')
ax.set_title('Three strongest components, as recorded (bad contacts removed)', fontsize=9, loc='left')
ax.legend(frameon=False, fontsize=8, loc='center right', bbox_to_anchor=(1.0, 0.35))
fig.savefig(f'{OUT}/fig3_ica_maps.png')
print('done')
