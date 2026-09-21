"""CCEP analysis prototype: CRP explained variance + reference/ICA checks.

Port of iEEGLAB's run_CRP.m (Miller et al., 2023) to numpy, plus reference
variants matching ieeglab_car.m conventions.
"""
import numpy as np, pandas as pd, scipy.io as sio

DATA = '/mnt/user-data/uploads/iEEGLAB/tutorial/dataset_seeg'


def load():
    m = sio.loadmat(f'{DATA}/sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set',
                    struct_as_record=False, squeeze_me=True)
    srate = float(m['srate'])
    data = np.asarray(m['data'], dtype=float)            # chan x time
    labels = [str(c.labels) for c in m['chanlocs']]
    ev = pd.read_csv(f'{DATA}/sub-02_ses-ieeg01_task-ccep_run-01_events.tsv', sep='\t')
    ev = ev[ev.status == 'good'].reset_index(drop=True)
    return data, labels, srate, ev


def epoch(data, srate, onsets, tmin=-0.5, tmax=1.0):
    n0, n1 = int(round(tmin * srate)), int(round(tmax * srate))
    t = np.arange(n0, n1) / srate
    idx = np.round(np.asarray(onsets) * srate).astype(int)
    keep = (idx + n0 >= 0) & (idx + n1 < data.shape[1])
    ep = np.stack([data[:, i + n0:i + n1] for i in idx[keep]], axis=2)  # chan x time x trial
    return ep, t, keep


def blank_and_baseline(ep, t, blank=0.015, bl=(-0.5, -0.05)):
    ep = ep.copy()
    b = (t >= bl[0]) & (t < bl[1])
    ep -= np.nanmean(ep[:, b, :], axis=1, keepdims=True)
    ep[:, np.abs(t) <= blank, :] = np.nan          # stimulation artifact blanked
    return ep


# ---------------------------------------------------------------- CRP
def ccep_proj(V):
    """V: T x K. Semi-normalized cross-projections, diagonal removed."""
    nrm = np.sqrt((V ** 2).sum(0))
    V0 = np.divide(V, nrm, out=np.zeros_like(V), where=nrm > 0)
    P = V0.T @ V
    np.fill_diagonal(P, np.nan)
    S = P.ravel()
    return S[~np.isnan(S)]


def kt_pca(X):
    w, F = np.linalg.eigh(X.T @ X)
    order = np.argsort(w)[::-1]
    F, s = F[:, order], np.sqrt(np.clip(w[order], 0, None))
    E = (X @ F) / np.maximum(s, np.finfo(float).eps)
    return E


def run_crp(V, t_win, t_step=5):
    """V: T x K (time x trials). Returns dict mirroring run_CRP.m outputs."""
    V = np.asarray(V, dtype=float)
    ok = np.isfinite(V).all(0)
    V, K = V[:, ok], int(ok.sum())
    T = V.shape[0]
    if T < 10 or K < 2:
        return None
    srate = 1.0 / np.median(np.diff(t_win))
    tpts = np.arange(10, T + 1, t_step)
    m = np.array([ccep_proj(V[:k, :]).mean() / np.sqrt(srate) for k in tpts])
    tt = int(np.argmax(m))
    tR = int(tpts[tt])
    V_tR = V[:tR, :]
    C = kt_pca(V_tR)[:, 0]
    al = C @ V_tR
    ep = V_tR - np.outer(C, al)
    expl = 1 - (ep ** 2).sum(0) / (V_tR ** 2).sum(0)
    S_tR = ccep_proj(V_tR) / np.sqrt(srate)
    return dict(tR_s=float(t_win[tR - 1]), C=C, al=al, expl_var=expl,
                n_trials=K, mean_proj=float(S_tR.mean()),
                t_stat=float(S_tR.mean() / (S_tR.std(ddof=1) / np.sqrt(S_tR.size))))


# ---------------------------------------------------- reference variants
def ref_none(ep, *_):
    return ep


def ref_car_fixed(ep, labels, sites):
    """Plain CAR over all channels, stimulated contacts included (the bug
    iEEGLAB fixed): one fixed linear operator."""
    return ep - np.nanmean(ep, axis=0, keepdims=True)


BAD = ['ROP3', 'ROP5', 'ROP6']   # flagged by quiet-period variance, robust z > 5


def ref_car_persite_good(ep, labels, sites, bad=BAD):
    """As ref_car_persite, with bad channels also kept out of the reference."""
    out = ep.copy()
    badi = [labels.index(b) for b in bad if b in labels]
    for s in np.unique(sites):
        tr = np.where(sites == s)[0]
        excl = [labels.index(c) for c in str(s).split('-') if c in labels] + badi
        sel = [i for i in range(len(labels)) if i not in excl]
        out[:, :, tr] -= np.nanmean(ep[sel][:, :, tr], axis=0, keepdims=True)
    return out


def ref_car_persite(ep, labels, sites):
    """iEEGLAB convention: per stimulation site, average over good channels
    with that site's stimulated pair excluded, applied trial by trial."""
    out = ep.copy()
    for s in np.unique(sites):
        tr = np.where(sites == s)[0]
        excl = [labels.index(c) for c in str(s).split('-') if c in labels]
        sel = [i for i in range(len(labels)) if i not in excl]
        out[:, :, tr] -= np.nanmean(ep[sel][:, :, tr], axis=0, keepdims=True)
    return out
