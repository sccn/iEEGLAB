"""Re-referencing schemes on ds004696 sub-02: the figures and the statistics of the comparison (effect
sizes with cluster-robust 95% confidence intervals over stimulation sites, areas under the ROC curve with
site-bootstrap intervals), from the output folder of compare_references_ds004696.m. The statistics are
written to ref_numbers.json in the figure folder.
Usage: python figures_rereferencing.py <compare_references output folder> <figure folder>
Needs NumPy, SciPy and Matplotlib."""
import sys, os, csv, json
import numpy as np
import scipy.io as sio
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e2dc'
BLUE, ORANGE, AQUA, LIGHT, GREY = '#2a78d6', '#eb6834', '#1baf7a', '#d9d7d0', '#a3a29c'
plt.rcParams.update({'font.family': ['Segoe UI', 'Arial', 'DejaVu Sans'], 'font.size': 8.5, 'axes.titlesize': 9,
                     'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlelocation': 'left',
                     'axes.titleweight': 'bold', 'legend.frameon': False, 'savefig.dpi': 200})
NAMES = {'none': 'Recording reference', 'car': 'Common average', 'hapwave': 'Lowest-variance 20% (published)',
         'carla': 'CARLA', 'ica_fieldtrip': 'ICA', 'ica': 'ICA, normalized weights'}
# scheme keys are those of compare_references_ds004696: 'ica_fieldtrip' is ICA re-referencing with its default
# (FieldTrip's scaling of the weights, as the authors' script), 'ica' the normalized weights
PUBREF = 'hapwave'   # the paper's reference: common average of the 20% lowest-variance contacts per 64-contact block
schemes = [s for s in NAMES if os.path.isfile(os.path.join(src, f'pairs_{s}.tsv'))]
print('schemes', schemes)


def num(x):
    try: return float(x)
    except ValueError: return np.nan


P = {}
for s in schemes:
    rows = list(csv.DictReader(open(os.path.join(src, f'pairs_{s}.tsv')), delimiter='\t'))
    P[s] = {k: np.array([num(r[k]) for r in rows]) for k in ['tR_ms', 'significant', 'pub_tested', 'pub_significant', 'pub_tR_ms', 'p', 'p_adj']}
    P[s]['site'] = np.array([r['site'] for r in rows]); P[s]['channel'] = np.array([r['channel'] for r in rows])
    P[s]['key'] = np.char.add(np.char.add(P[s]['site'], '|'), P[s]['channel'])

# ---------- summary numbers (printed for the page and the report) ----------
summ = {}
for s in schemes:
    d = P[s]; t = ~np.isnan(d['significant']); m = (d['pub_tested'] == 1) & t
    i, p = d['significant'][m] == 1, d['pub_significant'][m] == 1
    b = m & (d['significant'] == 1) & (d['pub_significant'] == 1)
    dt = np.abs(d['tR_ms'][b] - d['pub_tR_ms'][b])
    summ[s] = dict(tested=int(t.sum()), sig=int((d['significant'] == 1).sum()), pub=int(m.sum()), pub_sig=int(p.sum()),
                   same=float(100 * np.mean(i == p)), both=int((i & p).sum()), only_i=int((i & ~p).sum()),
                   only_p=int((~i & p).sum()), neither=int((~i & ~p).sum()), tr_med=float(np.median(dt)),
                   tr_50=float(100 * np.mean(dt <= 50)), n_tr=int(b.sum()))
    print(s, summ[s])

# ---------- effect sizes (effects_<scheme>.tsv): d of each pair, jackknife 95% CI ----------
E = {}
for s in schemes:
    f = os.path.join(src, f'effects_{s}.tsv')
    if os.path.isfile(f):
        E[s] = {r['site'] + '|' + r['channel']: (num(r['d']), num(r['d_ci_low']), num(r['d_ci_high']),
                                                 num(r['proj_uv']), num(r['proj_sd_uv']))
                for r in csv.DictReader(open(f), delimiter='\t')}
eff = [s for s in schemes if s in E]
print('effect sizes for', eff)


def fmt_d(s, key):
    if s not in E or key not in E[s]: return ''
    return f'd = {E[s][key][0]:.1f}'   # per-pair jackknife CIs are conservative and wide: not shown


# ---------- figure 1: one pair under each scheme, and what each scheme removed ----------
L = {s: sio.loadmat(os.path.join(src, f'averages_{s}.mat'), squeeze_me=True) for s in schemes}
times = L[schemes[0]]['times']; sites = list(L[schemes[0]]['site_names']); labels = list(L[schemes[0]]['labels'])
# example: a published significant pair on which the schemes disagree, with a clear response after the pulse
# artifact (50-500 ms)
w = (times >= 10) & (times <= 300); wr = (times >= 50) & (times <= 500)
d = P['none']; S = np.array([P[s]['significant'] for s in schemes])
# a pair the paper and most referenced schemes find, but the recording reference misses
cand = np.where((d['pub_significant'] == 1) & (d['significant'] == 0) & (np.nanmean(S[1:] == 1, axis=0) > 0.5))[0]
amp = []
for k in cand:
    a = L['none']['A'][sites.index(d['site'][k]), labels.index(d['channel'][k]), wr]
    amp.append(np.nanmax(np.abs(a)) if np.all(np.isfinite(a)) else -1)
k = cand[int(np.argmax(amp))]
site, chan = d['site'][k], d['channel'][k]
gi, ci = sites.index(site), labels.index(chan)
print('example pair', site, chan, {s: P[s]['significant'][k] for s in schemes})
others = [s for s in schemes if s != 'none']
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.2), sharex=False, gridspec_kw={'hspace': 0.55, 'wspace': 0.32})
axs = axs.ravel()
y0 = L['none']['A'][gi, ci]
ys = [y0] + [L[s]['A'][gi, ci] for s in others]
lim = 1.15 * np.nanmax(np.abs(np.concatenate([yy[times >= 30] for yy in ys])))   # the pulse artifact is clipped
for j, s in enumerate(others[:5]):
    ax = axs[j]
    ax.axhline(0, color=GRID, lw=0.8, zorder=0); ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    ax.plot(times, y0, color=GREY, lw=1.2, label='Recording reference')
    ax.plot(times, L[s]['A'][gi, ci], color=BLUE, lw=1.2, label='This scheme')
    ax.set_ylim(-lim, lim); ax.set_xlim(times[0], times[-1])
    ax.set_title(NAMES[s], fontsize=8, pad=13)
    lab = 'responds' if P[s]['significant'][k] == 1 else 'no response'
    if fmt_d(s, d['key'][k]): lab += '; ' + fmt_d(s, d['key'][k])
    ax.text(0, 1.03, lab, transform=ax.transAxes, ha='left', va='bottom', fontsize=7, color=INK2)
    if j % 3 == 0: ax.set_ylabel('µV')
    if j >= 2: ax.set_xlabel('Time (ms)')
fig.legend(*axs[0].get_legend_handles_labels(), ncol=2, loc='lower left', bbox_to_anchor=(0.08, 0.99), fontsize=8, handlelength=1.4)
ax = axs[5]
rem = [float(np.nanmedian(np.sqrt(np.nanmean(L[s]['R'][:, :, w].astype(float) ** 2, axis=2)))) for s in others]
y = np.arange(len(others))
ax.barh(y, rem, color=BLUE, height=0.6)
for yy, v in zip(y, rem): ax.text(v, yy, f' {v:.0f}', va='center', fontsize=7.5, color=INK2)
SHORT = {'car': 'Common average', 'hapwave': 'Published ref.', 'carla': 'CARLA', 'ica_fieldtrip': 'ICA',
         'ica': 'ICA, normalized'}
ax.set_yticks(y, [SHORT[s] for s in others], fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel('µV (median RMS, 10-300 ms)'); ax.set_title('Signal removed', fontsize=8.5)
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0); ax.set_xlim(0, max(rem) * 1.3)
for j in range(len(others), 5): axs[j].set_visible(False)
fig.savefig(os.path.join(out, 'ref1_example.png'), bbox_inches='tight'); plt.close(fig)
print('removed', dict(zip(others, [round(v, 1) for v in rem])))

# ---------- figure 2: responding pairs and agreement with the published CRP ----------
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.3), gridspec_kw={'width_ratios': [1, 1.6, 0.9], 'wspace': 0.12})
y = np.arange(len(schemes)); lab = [NAMES[s] for s in schemes]
ax = axs[0]; v = [summ[s]['sig'] for s in schemes]
ax.barh(y, v, color=BLUE, height=0.6)
for yy, vv in zip(y, v): ax.text(vv, yy, f' {vv}', va='center', fontsize=7.5, color=INK2)
ax.set_yticks(y, lab); ax.invert_yaxis(); ax.set_xlim(0, max(v) * 1.3)
ax.set_xlabel(f"Pairs responding (of {summ[schemes[0]]['tested']})"); ax.set_title('a  All pairs')
ax = axs[1]
cats = [('both', 'Both', BLUE), ('only_i', 'iEEGLAB only', ORANGE), ('only_p', 'Published only', AQUA), ('neither', 'Neither', LIGHT)]
for yy, s in zip(y, schemes):
    left = 0
    for key, _, col in cats:
        vv = summ[s][key]
        ax.barh(yy, vv, left=left, color=col, height=0.6, edgecolor='white', linewidth=1.2)
        if vv >= 30: ax.text(left + vv / 2, yy, str(vv), ha='center', va='center', fontsize=7.2, color='white' if key != 'neither' else INK2)
        left += vv
    ax.text(left, yy, f"  {summ[s]['same']:.0f}%", va='center', fontsize=7.5, color=INK)
ax.set_yticks(y, []); ax.invert_yaxis(); ax.set_xlim(0, summ[schemes[0]]['pub'] * 1.14)
ax.set_xlabel(f"The {summ[schemes[0]]['pub']} published pairs (% same decision)"); ax.set_title('b  Against the published CRP')
ax.legend([Patch(color=c) for _, _, c in cats], [l for _, l, _ in cats], ncol=4, loc='lower left', bbox_to_anchor=(0, 1.1),
          fontsize=7.5, handlelength=1, columnspacing=1)
ax = axs[2]; v = [summ[s]['tr_med'] for s in schemes]
ax.barh(y, v, color=BLUE, height=0.6)
for yy, vv in zip(y, v): ax.text(vv, yy, f' {vv:.0f}', va='center', fontsize=7.5, color=INK2)
ax.set_yticks(y, []); ax.invert_yaxis(); ax.set_xlim(0, max(v) * 1.35)
ax.set_xlabel('Median difference (ms)'); ax.set_title('c  Response duration')
for ax in axs: ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
fig.savefig(os.path.join(out, 'ref2_pairs.png'), bbox_inches='tight'); plt.close(fig)

# ---------- figure 3: the published pairs as matrices ----------
d = P[schemes[0]]; m = d['pub_tested'] == 1
psites = sorted(set(d['site'][m]), key=lambda x: sites.index(x)); pch = sorted(set(d['channel'][m]), key=lambda x: labels.index(x))
def mat(s, published=False):
    M = np.full((len(psites), len(pch)), np.nan); dd = P[s]
    for kk in np.where(dd['pub_tested'] == 1)[0]:
        v = dd['pub_significant'][kk] if published else dd['significant'][kk]
        M[psites.index(dd['site'][kk]), pch.index(dd['channel'][kk])] = v
    return M
panels = [('Published', mat(schemes[0], True))] + [(NAMES[s], mat(s)) for s in schemes]
nc = 4; nr = int(np.ceil(len(panels) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(7.4, 2.15 * nr), gridspec_kw={'hspace': 0.45, 'wspace': 0.08})
cmap = ListedColormap([LIGHT, BLUE]); cmap.set_bad('white')
for ax, (name, M) in zip(axs.ravel(), panels):
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    ax.set_title(name, fontsize=8.5); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(GRID)
for ax in axs.ravel()[len(panels):]: ax.set_visible(False)
axs.ravel()[0].set_ylabel(f'{len(psites)} stimulation sites'); axs.ravel()[0].set_xlabel(f'{len(pch)} recording contacts')
fig.legend([Patch(color=BLUE), Patch(color=LIGHT)], ['Response', 'No response'], ncol=2, loc='lower right', fontsize=7.5,
           bbox_to_anchor=(0.9, 0.02 if nr > 1 else -0.08))
fig.savefig(os.path.join(out, 'ref3_matrices.png'), bbox_inches='tight'); plt.close(fig)
# ---------- figure 4: effect sizes with 95% CIs ----------
# The published pairs that every scheme tests; the mean d over pairs, 95% cluster-robust t-interval with the
# stimulation sites as clusters (pairs of one site share its trials); differences from the paper's reference paired
# per pair. AUC: percentile bootstrap over sites.
effsum = {}
if eff and PUBREF in eff:
    d0 = P[schemes[0]]; m = d0['pub_tested'] == 1
    keys, resp, site_of = d0['key'][m], d0['pub_significant'][m] == 1, d0['site'][m]
    ok = np.array([all(np.isfinite(E[s].get(k, (np.nan,))[0]) for s in eff) for k in keys])
    keys, resp, site_of = keys[ok], resp[ok], site_of[ok]
    D = np.array([[E[s][k][0] for k in keys] for s in eff])                    # schemes x pairs
    us = np.unique(site_of); idx = [np.where(site_of == u)[0] for u in us]
    rng = np.random.default_rng(20261111); B = 2000
    boots = {'auc': np.empty((len(eff), B))}   # AUC: percentile bootstrap over sites

    def auc(a, b):                                                              # P(d responding > d not responding)
        return (np.sum(a[:, None] > b[None, :]) + 0.5 * np.sum(a[:, None] == b[None, :])) / (a.size * b.size)
    for bb in range(B):
        sel = np.concatenate([idx[i] for i in rng.integers(0, len(us), len(us))])
        r, n = sel[resp[sel]], sel[~resp[sel]]
        boots['auc'][:, bb] = [auc(D[j, r], D[j, n]) for j in range(len(eff))]
    j0 = eff.index(PUBREF)
    PM = np.array([[E[s][k][3] for k in keys] for s in eff]); PS = np.array([[E[s][k][4] for k in keys] for s in eff])

    def crt(x, g):
        # 95% cluster-robust t-interval of the mean over pairs, clusters = stimulation sites (CR1, G - 1 df;
        # Cameron & Miller, 2015): a percentile bootstrap of 20 clusters gives intervals that are too narrow
        u = np.unique(g); G = len(u); m = x.mean()
        se = np.sqrt(G / (G - 1) * np.sum([np.sum(x[g == v] - m) ** 2 for v in u])) / x.size
        q = stats.t.ppf(0.975, G - 1)
        return [float(m - q * se), float(m + q * se)]
    for j, s in enumerate(eff):
        ci = lambda x: [float(v) for v in np.percentile(x, [2.5, 97.5])]
        r, n = resp, ~resp
        effsum[s] = dict(resp=float(D[j, r].mean()), resp_ci=crt(D[j, r], site_of[r]),
                         non=float(D[j, n].mean()), non_ci=crt(D[j, n], site_of[n]),
                         diff=float(D[j, r].mean() - D[j0, r].mean()), diff_ci=crt(D[j, r] - D[j0, r], site_of[r]),
                         diff_non=float(D[j, n].mean() - D[j0, n].mean()),
                         diff_non_ci=crt(D[j, n] - D[j0, n], site_of[n]),
                         auc=float(auc(D[j, r], D[j, n])), auc_ci=ci(boots['auc'][j]),
                         proj_ratio=float(np.median(PM[j, r] / PM[j0, r])), sd_ratio=float(np.median(PS[j, r] / PS[j0, r])),
                         n_resp=int(r.sum()), n_non=int(n.sum()), n_sites=int(len(us)),
                         n_sites_resp=int(len(np.unique(site_of[r]))), n_sites_non=int(len(np.unique(site_of[n]))))
        print(s, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in effsum[s].items()})
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 0.5 + 0.33 * len(eff)), sharey=True, gridspec_kw={'wspace': 0.2})
    y = np.arange(len(eff))
    panels = [('resp', 'resp_ci', f"a  Response in the paper ({int(resp.sum())})", 'Mean d'),
              ('non', 'non_ci', f"b  No response in the paper ({int((~resp).sum())})", 'Mean d'),
              ('diff', 'diff_ci', 'c  a, minus the published', 'Difference in d')]
    for ax, (v, c, title, xl) in zip(axs, panels):
        if v != 'resp': ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
        for yy, s in zip(y, eff):
            col = ORANGE if s == PUBREF else BLUE
            lo, hi = effsum[s][c]
            ax.plot([lo, hi], [yy, yy], color=col, lw=2, solid_capstyle='round')
            ax.plot(effsum[s][v], yy, 'o', color=col, ms=5, mec='white', mew=1)
        ax.set_title(title, fontsize=8.5); ax.set_xlabel(xl); ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    axs[0].set_yticks(y, [NAMES[s] for s in eff]); axs[0].invert_yaxis()
    fig.savefig(os.path.join(out, 'ref4_effects.png'), bbox_inches='tight'); plt.close(fig)

    # the same over all pairs: those every scheme finds responding, and those none does
    allk = [set(P[s]['key'][P[s]['significant'] == 1]) for s in schemes]
    nonek = [set(P[s]['key'][P[s]['significant'] == 0]) for s in schemes]
    cons = {'all_respond': set.intersection(*allk), 'none_respond': set.intersection(*nonek)}
    for nm, ks in cons.items():
        ks = [k for k in ks if all(np.isfinite(E[s].get(k, (np.nan,))[0]) for s in eff)]
        effsum[nm] = {s: float(np.mean([E[s][k][0] for k in ks])) for s in eff}; effsum[nm]['n'] = len(ks)
        print(nm, len(ks), {s: round(v, 2) for s, v in effsum[nm].items() if s != 'n'})
fd = os.path.join(src, 'ica_decomposition.mat')
ica_rank = int(sio.loadmat(fd, squeeze_me=True, struct_as_record=False)['info'].rank) if os.path.isfile(fd) else None
json.dump({'example': [site, chan], 'example_d': {s: E[s].get(d['key'][k]) for s in eff} if eff else {},
           'ica_rank': ica_rank,
           'summary': summ, 'removed_uv': dict(zip(others, rem)), 'effects': effsum},
          open(os.path.join(out, 'ref_numbers.json'), 'w'), indent=1)
