import os, numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
H = os.path.expanduser('~')
m = np.load(H + '/sub02_meta.npz', allow_pickle=True)
X = np.lib.format.open_memmap(H + '/sub02_epochs.npy', mode='r')
labels, status, sites, t = m['labels'], m['status'], m['sites'], m['times']
srate = float(m['srate'])
good = status == 'good'
BL = (t > -0.5) & (t < -0.05)


def ccep_proj(V):
    nrm = np.sqrt((V ** 2).sum(0))
    V0 = np.divide(V, nrm, out=np.zeros_like(V), where=nrm > 0)
    P = V0.T @ V
    np.fill_diagonal(P, np.nan)
    S = P.ravel()
    return S[~np.isnan(S)]


def kt_pca(A):
    w, F = np.linalg.eigh(A.T @ A)
    o = np.argsort(w)[::-1]
    F, s = F[:, o], np.sqrt(np.clip(w[o], 0, None))
    return (A @ F) / np.maximum(s, np.finfo(float).eps)


def run_crp(V, t_win, t_step=20):
    V = np.asarray(V, float)
    ok = np.isfinite(V).all(0)
    V = V[:, ok]
    T, K = V.shape
    if T < 10 or K < 2:
        return None
    sr = 1.0 / np.median(np.diff(t_win))
    tp = np.arange(10, T + 1, t_step)
    mm = np.array([ccep_proj(V[:k, :]).mean() / np.sqrt(sr) for k in tp])
    tR = int(tp[int(np.argmax(mm))])
    VtR = V[:tR, :]
    C = kt_pca(VtR)[:, 0]
    al = C @ VtR
    ep = VtR - np.outer(C, al)
    expl = 1 - (ep ** 2).sum(0) / (VtR ** 2).sum(0)
    return dict(tR=float(t_win[tR - 1]), C=C, al=al, expl=float(np.mean(expl)), K=K)


def prep(trial_idx, site, blank=0.010, hp=0.5):
    """Epochs for one site: baseline-correct, blank, per-site CAR over good non-stim channels."""
    Y = np.asarray(X[:, :, trial_idx], dtype=np.float64)
    b, a = butter(4, hp / (srate / 2), 'high')
    Y = filtfilt(b, a, Y, axis=1)
    Y -= Y[:, BL, :].mean(1, keepdims=True)
    stim = set(str(site).split('-'))
    sel = np.array([i for i in range(len(labels)) if good[i] and labels[i] not in stim])
    Y -= Y[sel].mean(0, keepdims=True)
    Y[:, np.abs(t) <= blank, :] = np.nan
    return Y, sel


def n1(V, t, lo=0.011, hi=0.100, bl=(-0.5, -0.05)):
    """Mean response; returns (neg-peak amp, lat) and (largest-deflection amp, lat, sign) and SNR."""
    mu = np.nanmean(V, 1)
    w = (t >= lo) & (t <= hi)
    base = mu[(t > bl[0]) & (t < bl[1])]
    sd = np.std(base)
    seg, ts = mu[w], t[w]
    i_neg = int(np.argmin(seg))
    i_abs = int(np.argmax(np.abs(seg)))
    return dict(neg_amp=seg[i_neg], neg_lat=ts[i_neg], abs_amp=seg[i_abs], abs_lat=ts[i_abs],
                sign=np.sign(seg[i_abs]), z=abs(seg[i_abs]) / sd if sd > 0 else np.nan)
