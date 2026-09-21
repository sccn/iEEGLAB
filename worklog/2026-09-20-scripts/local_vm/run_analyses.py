"""Analyses behind worklog/2026-09-20-native-rate-numbers.md, as run on 20 Sep 2026.

These steps were run interactively in the Cowork Linux workspace on the local
machine; this file collects them in the order they ran so they can be reproduced.
Prerequisites, in order (each step writes to $HOME, outside the repo):
    python3 extract3.py        # ~/sub02_epochs.npy  (237 ch x 3072 x 695, -0.5..1.0 s, 2048 Hz)
    python3 extract_surr.py    # ~/sub02_surr.npy    (same, random onsets >= 2 s from any pulse)
    python3 extract_long.py    # ~/long_epochs.npz   (-2..+2 s, 48 responsive pairs, 6 sites)
    python3 ica_run.py raw 0 ; python3 ica_run.py car 0 ; python3 ica_run.py raw 1
Data: HAPwave (ds004696) unzipped at ~/Downloads/v1.0.0, sub-02 run-01 (MEF3).
Results were written to Documents/iEEGLAB_data/results_sub02/.
Needs: pymef numpy scipy pandas python-picard mne erdetect.
"""
import os, re, numpy as np, pandas as pd
from scipy.signal import butter, filtfilt, hilbert
H = os.path.expanduser('~')
OUT = os.path.join(H, 'mnt/iEEGLAB_data/results_sub02')


def crp_real_vs_surrogate(n_sites=24, seed=0):
    """CRP per site x contact, real and surrogate epochs -> crp_native.csv, crp_paired.csv."""
    from analysis import X, labels, sites, t, srate, good, BL, run_crp, n1
    S = np.lib.format.open_memmap(H + '/sub02_surr.npy', mode='r')
    rng = np.random.default_rng(seed)
    usites = [s for s in sorted(set(sites)) if (sites == s).sum() >= 10]
    pick = list(rng.choice(usites, size=n_sites, replace=False))
    win = (t >= 0.015) & (t <= 1.0); tw = t[win]
    b, a = butter(4, 0.5 / (srate / 2), 'high')

    def prep(A, idx, site):
        Y = filtfilt(b, a, np.asarray(A[:, :, idx], float), axis=1)
        Y -= Y[:, BL, :].mean(1, keepdims=True)
        stim = set(str(site).split('-'))
        sel = np.array([i for i in range(len(labels)) if good[i] and labels[i] not in stim])
        Y -= Y[sel].mean(0, keepdims=True)
        Y[:, np.abs(t) <= 0.010, :] = np.nan
        return Y, sel

    rows = []
    for s in pick:                       # the original run stopped after 20 sites (time budget)
        tr = np.where(sites == s)[0]
        for tag, A in (('real', X), ('surrogate', S)):
            Y, sel = prep(A, tr, s)
            for ci in sel:
                r = run_crp(Y[ci][win], tw)
                if r is None: continue
                rows.append(dict(cond=tag, site=str(s), chan=str(labels[ci]), K=r['K'],
                                 tR=r['tR'], expl=r['expl'], **n1(Y[ci], t)))
    df = pd.DataFrame(rows); df.to_csv(OUT + '/crp_native.csv', index=False)
    r = df[df.cond == 'real'].set_index(['site', 'chan'])
    s = df[df.cond == 'surrogate'].set_index(['site', 'chan'])
    j = r.join(s, rsuffix='_s', how='inner'); j['delta'] = j.expl - j.expl_s
    j.reset_index().to_csv(OUT + '/crp_paired.csv', index=False)
    return j


def split_half_control(n_sites=8, seed=0):
    """Select responsive pairs on odd trials, evaluate CRP on even trials -> splithalf.csv.
    Result on 20 Sep: 1568 pairs over 8 sites, 416 responsive; responsive delta +0.172
    (68% positive), non-responsive +0.019."""
    from analysis import X, labels, sites, t, srate, good, BL, run_crp, n1
    S = np.lib.format.open_memmap(H + '/sub02_surr.npy', mode='r')
    rng = np.random.default_rng(seed)
    usites = [s for s in sorted(set(sites)) if (sites == s).sum() >= 10]
    pick = list(rng.choice(usites, size=n_sites, replace=False))
    win = (t >= 0.015) & (t <= 1.0); tw = t[win]
    b, a = butter(4, 0.5 / (srate / 2), 'high')

    def prep(A, idx, site):
        Y = filtfilt(b, a, np.asarray(A[:, :, idx], float), axis=1)
        Y -= Y[:, BL, :].mean(1, keepdims=True)
        stim = set(str(site).split('-'))
        sel = np.array([i for i in range(len(labels)) if good[i] and labels[i] not in stim])
        Y -= Y[sel].mean(0, keepdims=True)
        Y[:, np.abs(t) <= 0.010, :] = np.nan
        return Y, sel

    rows = []
    for s in pick:
        tr = np.where(sites == s)[0]
        odd, even = tr[::2], tr[1::2]
        if len(even) < 4: continue
        Yo, sel = prep(X, odd, s); Ye, _ = prep(X, even, s); Se, _ = prep(S, even, s)
        for ci in sel:
            r_e = run_crp(Ye[ci][win], tw); r_s = run_crp(Se[ci][win], tw)
            if r_e is None or r_s is None: continue
            rows.append(dict(site=str(s), chan=str(labels[ci]), z_sel=n1(Yo[ci], t)['z'],
                             expl_even=r_e['expl'], expl_surr=r_s['expl']))
    df = pd.DataFrame(rows); df['delta'] = df.expl_even - df.expl_surr
    df.to_csv(OUT + '/splithalf.csv', index=False)
    return df


def bad_channels():
    """Filtered (1-200 Hz) baseline SD and within-shaft neighbour correlation vs channels.tsv -> badchan2.npz.
    Result: 39/237 clinician-bad, bimodal (25 quiet, 14 loud); z>5 catches 12/14 loud, 0/25 quiet."""
    from analysis import X, labels, status, t, srate
    bl = (t > -0.5) & (t < -0.05); nb = bl.sum(); ntr = min(60, X.shape[2])
    b, a = butter(4, [1 / (srate / 2), 200 / (srate / 2)], 'band')
    B = np.zeros((len(labels), nb * ntr))
    for i in range(len(labels)):
        B[i] = filtfilt(b, a, np.asarray(X[i][bl][:, :ntr], float), axis=0).T.ravel()
    sd = B.std(1); med = np.median(sd); mad = np.median(np.abs(sd - med)) * 1.4826
    z = (sd - med) / mad
    shaft = np.array([re.sub(r'\d+$', '', l) for l in labels]); C = np.corrcoef(B)
    nbr = np.array([np.max(np.abs(C[i, (shaft == shaft[i]) & (np.arange(len(labels)) != i)]))
                    if ((shaft == shaft[i]).sum() > 1) else np.nan for i in range(len(labels))])
    np.savez(OUT + '/badchan2.npz', sd=sd, z=z, nbrcorr=nbr, labels=labels, status=status, shaft=shaft)


def erdetect_comparison(n_sites=20, seed=0):
    """erdetect 2.6.2 on the same preprocessed averages -> erdetect_vs_mine.csv.
    NOTE: erdetect's default baseline (-1, -0.1) s falls outside a -0.5 s epoch and returns
    zero detections silently; it is set to (-0.45, -0.05) here."""
    from erdetect.core.config import set as cset
    from erdetect.core.detection import ieeg_detect_er
    from analysis import X, labels, sites, t, srate, good, prep
    cset((-0.45, -0.05), 'detection', 'std_base', 'baseline_epoch')
    rng = np.random.default_rng(seed)
    usites = [s for s in sorted(set(sites)) if (sites == s).sum() >= 10]
    pick = list(rng.choice(usites, size=n_sites, replace=False)); gi = np.where(good)[0]
    AVG = np.zeros((len(gi), len(pick), X.shape[1]))
    for k, s in enumerate(pick):
        Y, _ = prep(np.where(sites == s)[0], s)
        AVG[:, k, :] = np.nan_to_num(np.nanmean(Y[gi], axis=2))
    onset = int(np.argmin(np.abs(t)))
    neg = ieeg_detect_er(AVG, onset, int(srate), detect_positive=False)
    pos = ieeg_detect_er(AVG, onset, int(srate), detect_positive=True)
    np.savez(H + '/erdetect_out.npz', neg_lat=neg[0], neg_amp=neg[1], pos_lat=pos[0], pos_amp=pos[1],
             chans=labels[gi], sites=np.array(pick), onset=onset, srate=srate, t=t)
    # merged with crp_paired.csv on (site, chan) -> erdetect_vs_mine.csv


def itpc_fixed(seed=3):
    """Morlet ITPC on +-2 s epochs, baseline -1.8..-1.0 s, equal trial counts, evoked and
    evoked-subtracted (induced) -> itpc_fixed.csv. Pre-stim ITPC 0.26-0.29 (= n=10 chance)."""
    from mne.time_frequency import tfr_array_morlet
    Z = np.load(H + '/long_epochs.npz', allow_pickle=True)
    srate = float(Z['srate']); t = Z['times']
    freqs = np.array([4, 6, 8, 12, 16, 24, 32, 48, 80, 120]); ncyc = np.where(freqs < 20, 5, 7).astype(float)
    pre_w = (t >= -1.8) & (t <= -1.0); post_w = (t >= 0.011) & (t <= 0.100)
    b, a = butter(4, 0.5 / (srate / 2), 'high'); rng = np.random.default_rng(seed)

    def prep(V):
        Y = filtfilt(b, a, V.astype(float), axis=0); Y -= Y[pre_w].mean(0, keepdims=True); return Y.T

    def itpc(Y, sub):
        if sub: Y = Y - Y.mean(0, keepdims=True)
        W = tfr_array_morlet(Y[:, None, :], srate, freqs, n_cycles=ncyc, output='complex',
                             zero_mean=True, verbose=False)[:, 0]
        if W.ndim == 4: W = W[:, :, 0, :]
        v = np.abs((W / np.maximum(np.abs(W), 1e-30)).mean(0))
        return v[:, pre_w].mean(1), v[:, post_w].mean(1)

    surr = {k.split('|')[1]: k for k in Z.files if k.startswith('SURR')}
    rows = []
    for k in [k for k in Z.files if '|' in k and not k.startswith('SURR')]:
        site, ch = k.split('|'); Yr = prep(Z[k]); n = Yr.shape[0]
        Ys = prep(Z[surr[ch]]); Ys = Ys[rng.choice(Ys.shape[0], n, replace=False)]
        for tag, Y in (('real', Yr), ('surrogate', Ys)):
            pe, po = itpc(Y, False); pi, poi = itpc(Y, True)
            for fi, f in enumerate(freqs):
                rows.append(dict(cond=tag, site=site, chan=ch, freq=int(f), n=n, evoked_pre=pe[fi],
                                 evoked_post=po[fi], induced_pre=pi[fi], induced_post=poi[fi]))
    pd.DataFrame(rows).to_csv(OUT + '/itpc_fixed.csv', index=False)


def hfb_long(seed=5):
    """70-170 Hz power ratio, 11-100 ms over -1.8..-1.0 s, equal trials -> hfb_long.csv.
    Result: median 1.74 (IQR 1.15-5.17) vs 1.12 surrogate; 24/72 pairs above surrogate 95th."""
    Z = np.load(H + '/long_epochs.npz', allow_pickle=True)
    srate = float(Z['srate']); t = Z['times']
    pre_w = (t >= -1.8) & (t <= -1.0); post_w = (t >= 0.011) & (t <= 0.100)
    bh, ah = butter(4, [70 / (srate / 2), 170 / (srate / 2)], 'band'); b0, a0 = butter(4, 0.5 / (srate / 2), 'high')
    rng = np.random.default_rng(seed); surr = {k.split('|')[1]: k for k in Z.files if k.startswith('SURR')}

    def ratio(V, n=None):
        Y = filtfilt(b0, a0, V.astype(float), axis=0)
        if n is not None and Y.shape[1] > n: Y = Y[:, rng.choice(Y.shape[1], n, replace=False)]
        A = np.abs(hilbert(filtfilt(bh, ah, Y, axis=0), axis=0))
        return A[post_w].mean() / A[pre_w].mean()

    rows = []
    for k in [k for k in Z.files if '|' in k and not k.startswith('SURR')]:
        site, ch = k.split('|'); n = Z[k].shape[1]
        rows.append(dict(site=site, chan=ch, n=n, real=ratio(Z[k]), surrogate=ratio(Z[surr[ch]], n)))
    pd.DataFrame(rows).to_csv(OUT + '/hfb_long.csv', index=False)
