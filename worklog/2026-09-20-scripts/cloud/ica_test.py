"""Does re-referencing change the ICA decomposition?

A: fixed common average (Scott's case: one fixed linear operator).
B: iEEGLAB's per-site, per-trial average with the stimulated pair excluded.
"""
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.signal import butter, filtfilt
from mne.preprocessing import infomax
from ccep_core import load, epoch, blank_and_baseline, ref_car_persite

rng = 0


def pca_whiten(X, k):
    """X: chan x time. Returns whitened k x time, whitener, mean."""
    mu = X.mean(1, keepdims=True)
    Xc = X - mu
    U, s, _ = np.linalg.svd(Xc @ Xc.T / Xc.shape[1])
    W = (U[:, :k] / np.maximum(s[:k], 1e-30) ** 0.5).T
    return W @ Xc, W, mu, s


def run_ica(X, k, seed=0):
    Z, W, mu, s = pca_whiten(X, k)
    U = infomax(Z.T, extended=True, random_state=seed, max_iter=1000)   # k x k
    unmix = U @ W                       # sources = unmix @ (X - mu)
    S = unmix @ (X - mu)
    A = np.linalg.pinv(unmix)           # chan x k maps
    return S, A, s


def match(S1, S2):
    """Best one-to-one matching of source time courses by |correlation|."""
    C = np.abs(np.corrcoef(S1, S2)[:S1.shape[0], S1.shape[0]:])
    r, c = linear_sum_assignment(-C)
    return np.sort(C[r, c])[::-1]


def spread(A):
    """Participation ratio of each map: 1 = one contact, N = uniform."""
    W = np.abs(A)
    return (W.sum(0) ** 2) / np.maximum((W ** 2).sum(0), 1e-30)


if __name__ == '__main__':
    data, labels, srate, ev = load()
    b, a = butter(4, 0.5 / (srate / 2), 'high')
    X = filtfilt(b, a, data, axis=1)
    N = len(labels)

    print('== eigenvalue spectrum ==')
    _, _, _, s_raw = pca_whiten(X, N)
    Xcar = X - X.mean(0, keepdims=True)
    _, _, _, s_car = pca_whiten(Xcar, N)
    print('as recorded, smallest/largest eigenvalue: %.3e' % (s_raw[-1] / s_raw[0]))
    print('after CAR,   smallest/largest eigenvalue: %.3e' % (s_car[-1] / s_car[0]))

    print('\n== A: fixed CAR, continuous data ==')
    S_raw, A_raw, _ = run_ica(X, N, rng)
    S_car, A_car, _ = run_ica(Xcar, N - 1, rng)
    m = match(S_raw, S_car)
    print('matched |r| of %d component pairs: median %.3f, min %.3f, n>0.9: %d'
          % (len(m), np.median(m), m.min(), (m > 0.9).sum()))

    print('\n== B: per-site, per-trial reference, epoched data ==')
    onsets = ev.onset.values
    sites_all = ev.electrical_stimulation_site.values
    ep, t, keep = epoch(X, srate, onsets)
    sites = sites_all[keep]
    epb = blank_and_baseline(ep, t)
    ep_ps = ref_car_persite(epb, labels, sites)

    def cat(e):
        M = e.reshape(e.shape[0], -1)
        return M[:, np.isfinite(M).all(0)]

    Er, Ep = cat(epb), cat(ep_ps)
    S_er, A_er, _ = run_ica(Er, N, rng)
    S_ep, A_ep, _ = run_ica(Ep, N - 1, rng)
    m2 = match(S_er, S_ep)
    print('matched |r| of %d component pairs: median %.3f, min %.3f, n>0.9: %d'
          % (len(m2), np.median(m2), m2.min(), (m2 > 0.9).sum()))

    print('\n== spatial spread of maps (participation ratio, max = %d) ==' % N)
    for tag, A in [('as recorded', A_raw), ('fixed CAR', A_car), ('per-site/per-trial', A_ep)]:
        p = spread(A)
        print('%-20s median %.2f  range %.2f-%.2f' % (tag, np.median(p), p.min(), p.max()))

    np.savez('ica_out.npz', A_raw=A_raw, A_car=A_car, A_ep=A_ep,
             labels=np.array(labels), m=m, m2=m2)
