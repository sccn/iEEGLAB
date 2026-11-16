"""Figures of the method validations (FieldTrip, CRP, N1, CARLA), from the outputs of the MATLAB
validation scripts (validate_fieldtrip_uci29, validate_crp_scripts, validate_erdetect_ds004080,
validate_carla_huang_ds004977), which run_validations.m writes to one folder.
Usage: python figures_validations.py <validation output folder> <figure folder> [<SubjectUCI29 folder>]
With FieldTrip's SubjectUCI29 folder, also draws the high-frequency test on the participant's pial surface.
Needs NumPy, SciPy and Matplotlib."""
import json, sys, os, csv
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

val, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e2dc'
BLUE, ORANGE, AQUA, LIGHT = '#2a78d6', '#eb6834', '#1baf7a', '#d9d7d0'
DIV = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#f0efec', '#f0a09f', '#b8302f'])
plt.rcParams.update({'font.family': ['Segoe UI', 'Arial', 'DejaVu Sans'], 'font.size': 8.5, 'axes.titlesize': 9,
                     'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlelocation': 'left',
                     'axes.titleweight': 'bold', 'legend.frameon': False, 'savefig.dpi': 220})


def letter(ax, s):
    ax.set_title(s, loc='left', fontweight='bold', color=INK)


def identity(ax, x, y):
    lo, hi = np.nanmin([x, y]), np.nanmax([x, y]); pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=LIGHT, lw=1, zorder=0)
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad); ax.set_aspect('equal')


# ---------- Figure 1: FieldTrip ----------
F = sio.loadmat(os.path.join(val, 'fieldtrip_uci29_figure.mat'), squeeze_me=True)
f, t = F['tf_freqs'], F['tf_times_ms']
A, B = F['tf_mean_ieeglab'], F['tf_mean_fieldtrip']
lim = np.nanmax(np.abs(np.r_[A.ravel(), B.ravel()]))
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.4), gridspec_kw={'width_ratios': [1, 1, 0.75], 'wspace': 0.75})
for ax, M, s in [(axs[0], A, 'a  iEEGLAB'), (axs[1], B, 'b  FieldTrip')]:
    im = ax.pcolormesh(t, f, M, cmap=DIV, vmin=-lim, vmax=lim, shading='nearest', rasterized=True)
    ax.axvline(0, color=INK2, lw=0.6, ls=':')
    ax.set(xlabel='Time (ms)', ylabel='Frequency (Hz)'); letter(ax, s)
cb = fig.colorbar(im, ax=axs[1], fraction=0.05, pad=0.03); cb.set_label('Relative change', color=INK2); cb.outline.set_visible(False)
x, y = F['hfa_fieldtrip'], F['hfa_ieeglab']
ax = axs[2]; identity(ax, x, y)
ax.plot(x, y, 'o', ms=3.2, mfc=BLUE, mec='white', mew=0.4)
ax.set(xlabel='FieldTrip', ylabel='iEEGLAB'); letter(ax, 'c  70-150 Hz, 0-800 ms')
fig.savefig(os.path.join(out, 'fieldtrip_tf_hfa.png'), bbox_inches='tight'); plt.close(fig)
print('fig1: r tf', float(F['tf_r']), 'hfa r', float(F['hfa_r']), 'n', len(x), 'max', float(F['hfa_max_abs_diff']))

# ---------- Figure 2: CRP against the authors' scripts ----------
rows = list(csv.DictReader(open(os.path.join(val, 'crp_vs_crp_scripts.tsv')), delimiter='\t'))
g = lambda k: np.array([float(r[k]) for r in rows])
ta, ti, pa, pi = g('tR_ms_authors'), g('tR_ms_ieeglab'), g('p_tR_authors'), g('p_tR_ieeglab')
fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.5), gridspec_kw={'wspace': 0.5})
ax = axs[0]; identity(ax, ta, ti); ax.plot(ta, ti, 'o', ms=4, mfc=BLUE, mec='white', mew=0.4)
ax.set(xlabel='Authors\' scripts (ms)', ylabel='iEEGLAB (ms)'); letter(ax, 'a  Response duration')
ok = (pa > 0) & (pi > 0)                                     # p below 1e-308 is stored as 0
print('fig2: p zero in authors', int(np.sum(pa == 0)), 'in iEEGLAB', int(np.sum(pi == 0)), 'same files', bool(np.all((pa == 0) == (pi == 0))))
la, li = -np.log10(pa[ok]), -np.log10(pi[ok])
ax = axs[1]; identity(ax, la, li); ax.plot(la, li, 'o', ms=4, mfc=BLUE, mec='white', mew=0.4)
ax.set(xlabel='Authors\' scripts', ylabel='iEEGLAB'); letter(ax, 'b  -log10 p')
fig.savefig(os.path.join(out, 'crp_vs_authors.png'), bbox_inches='tight'); plt.close(fig)
print('fig2: files', len(ta), 'same tR', int(np.sum(ta == ti)), 'max p rel diff', float(np.max(np.abs(pa - pi) / pa)))

# ---------- Figure 4: N1 against ER-detect ----------
E = sio.loadmat(os.path.join(val, 'validate_erdetect_ds004080_out.mat'), squeeze_me=True, struct_as_record=False)['out']
rowsE = [('ER-detect\'s rule', E.B), ('iEEGLAB\'s test', E.C)]
fig, ax = plt.subplots(figsize=(6.0, 1.5))
cats = [('both', 'Both', BLUE), ('ieeglab_only', 'iEEGLAB only', ORANGE), ('erdetect_only', 'ER-detect only', AQUA),
        ('neither', 'Neither', LIGHT)]
for k, (name, S) in enumerate(rowsE):
    left = 0
    for fld, lab, col in cats:
        v = int(getattr(S, fld))
        if v:
            ax.barh(k, v, left=left, color=col, height=0.6, edgecolor='white', linewidth=1.5)
            if v >= 20:
                ax.text(left + v / 2, k, str(v), ha='center', va='center', color='white' if fld != 'neither' else INK2, fontsize=8)
        left += v
ax.set_yticks([0, 1], [r[0] for r in rowsE]); ax.invert_yaxis()
ax.set_xlim(0, 672); ax.set_xlabel('Stimulation-recording pairs (672)')
from matplotlib.patches import Patch
ax.legend([Patch(color=c) for _, _, c in cats], [l for _, l, _ in cats], ncol=4, loc='lower left',
          bbox_to_anchor=(0, 1.02), fontsize=8, handlelength=1)
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
fig.savefig(os.path.join(out, 'n1_vs_erdetect.png'), bbox_inches='tight'); plt.close(fig)
print('fig4:', [(n, int(S.both), int(S.ieeglab_only), int(S.erdetect_only), int(S.neither), round(float(S.kappa), 3)) for n, S in rowsE])

# ---------- Figure 5: CARLA and the common average ----------
J = json.load(open(os.path.join(val, 'carla_fig5.json')))
R, C = J['ref'], J['ieeglab']; fl = lambda v: np.asarray(v, dtype=float)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.4), gridspec_kw={'wspace': 0.35})
ax = axs[0]
vr, vc = fl(R['vars_sorted']), fl(C['vars_sorted'])
ax.axvspan(0.5, R['n'] + 0.5, color='#eef3fa', zorder=0)
ax.plot(np.arange(1, len(vr) + 1), vr, '-', color=INK, lw=1.6, label='Authors\' code')
ax.plot(np.arange(1, len(vc) + 1), vc, '--', color=BLUE, lw=1.2, label='iEEGLAB')
ax.text(R['n'] / 2, np.nanmax(vr) * 0.9, f"{R['n']} of {R['n_total']} channels\nin the average", ha='center', va='top', color=INK2, fontsize=8)
ax.set(xlabel='Channel rank', ylabel='Mean cross-trial covariance'); letter(ax, 'a  CARLA channel selection')
ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.72), fontsize=8)
ax = axs[1]
tt = 1000 * fl(J['tt']); ca, cc = fl(R['CAR_mean']), fl(C['CAR_mean'])
ax.plot(tt, ca, '-', color=INK, lw=1.6, label='Authors\' code'); ax.plot(tt, cc, '--', color=BLUE, lw=1.2, label='iEEGLAB')
ax.axhline(0, color=GRID, lw=0.8, zorder=0)
ax.set(xlabel='Time (ms)', ylabel='µV', xlim=(tt[0], tt[-1])); letter(ax, 'b  Common average')
ax.legend(loc='upper right', fontsize=8)
fig.savefig(os.path.join(out, 'carla_figure5.png'), bbox_inches='tight'); plt.close(fig)
print('fig5: n', R['n'], C['n'], 'CAR r', float(np.corrcoef(ca, cc)[0, 1]), 'max diff', float(np.max(np.abs(ca - cc))), 'peak', float(np.max(np.abs(ca))))

# ---------- FieldTrip data: high-frequency activity on the participant's cortex ----------
# Optional third argument: FieldTrip's SubjectUCI29 folder (SubjectUCI29_elec_acpc_f.mat, freesurfer/surf/lh.pial).
# Left hemisphere seen from the left (anterior to the left), each grid contact coloured by its 70-150 Hz power
# against baseline (fieldtrip_uci29_band_stats.tsv, validate_fieldtrip_uci29); large dots: significant after FDR.
if len(sys.argv) > 3:
    import struct
    from matplotlib.collections import PolyCollection
    uci = sys.argv[3]

    def read_fs_surface(f):
        with open(f, 'rb') as h:
            h.read(3); h.readline(); h.readline()
            nv, nf = struct.unpack('>2i', h.read(8))
            v = np.frombuffer(h.read(nv * 12), '>f4').reshape(nv, 3).astype(float)
            fc = np.frombuffer(h.read(nf * 12), '>i4').reshape(nf, 3)
        return v, fc

    V, Fc = read_fs_surface(os.path.join(uci, 'freesurfer', 'surf', 'lh.pial'))
    el = sio.loadmat(os.path.join(uci, 'SubjectUCI29_elec_acpc_f.mat'), squeeze_me=True, struct_as_record=False)['elec_acpc_f']
    pos = dict(zip([str(x) for x in el.label], np.asarray(el.elecpos, dtype=float)))
    rows = list(csv.DictReader(open(os.path.join(val, 'fieldtrip_uci29_band_stats.tsv')), delimiter='\t'))
    grid = [r for r in rows if r['channel'].startswith(('LPG', 'LTG')) and r['channel'] in pos]
    tri = V[Fc]                                                   # faces x 3 vertices x xyz
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
    light = np.array([-0.85, 0.2, 0.5]); light /= np.linalg.norm(light)   # from the viewer's side, a little above
    shade = 0.35 + 0.6 * np.abs(nrm @ light)
    order = np.argsort(-tri[:, :, 0].mean(1))                     # far faces (larger x) first
    xy = np.stack([-tri[:, :, 1], tri[:, :, 2]], axis=2)[order]  # screen: -y (anterior left), z
    grey = np.clip(shade[order], 0, 1)[:, None] * np.array([0.93, 0.92, 0.9])
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.add_collection(PolyCollection(xy, facecolors=grey, edgecolors='none', antialiaseds=False, rasterized=True))
    db = np.array([float(r['band_db']) for r in grid]); sig = np.array([r['significant'] in ('1', 'true', 'True') for r in grid])
    P = np.array([pos[r['channel']] for r in grid])
    lim = max(np.abs(db).max(), 0.5)
    ax.scatter(-P[~sig, 1], P[~sig, 2], c=db[~sig], cmap=DIV, vmin=-lim, vmax=lim, s=14, edgecolors=MUTED, linewidths=0.4, zorder=3)
    sc = ax.scatter(-P[sig, 1], P[sig, 2], c=db[sig], cmap=DIV, vmin=-lim, vmax=lim, s=46, edgecolors=INK, linewidths=1.0, zorder=4)
    ax.set_xlim(xy[..., 0].min() - 3, xy[..., 0].max() + 3)
    ax.set_ylim(xy[..., 1].min() - 3, xy[..., 1].max() + 3); ax.set_aspect('equal'); ax.axis('off')
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02); cb.set_label('70-150 Hz, 0-800 ms (dB)', color=INK2); cb.outline.set_visible(False)
    fig.savefig(os.path.join(out, 'fieldtrip_surface.png'), bbox_inches='tight', dpi=220); plt.close(fig)
    print('surface:', len(grid), 'grid contacts,', int(sig.sum()), 'significant:', [r['channel'] for r, s in zip(grid, sig) if s])
