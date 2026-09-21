import os, time, numpy as np, pandas as pd
from pymef.mef_session import MefSession
H = os.path.expanduser('~')
V = os.path.expanduser('~/mnt/v1.0.0/sub-02/ses-ieeg01/ieeg')
B = 'sub-02_ses-ieeg01_task-ccep_run-01'
ms = MefSession(os.path.join(V, B + '_ieeg.mefd'), None)
info = {c['name']: c for c in ms.read_ts_channel_basic_info()}
m = np.load(H + '/sub02_meta.npz', allow_pickle=True)
labels, starts, t = m['labels'], m['starts'], m['times']
srate = float(m['srate']); nsamp = int(info[labels[0]]['nsamp'][0])
pre, post = int(round(0.5 * srate)), int(round(1.0 * srate))
rng = np.random.default_rng(7)
cand = []
while len(cand) < len(starts):
    s = int(rng.integers(pre + 1, nsamp - post - 1))
    if np.min(np.abs(starts - s)) > 2 * srate:
        cand.append(s)
cand = np.array(sorted(cand))
path = H + '/sub02_surr.npy'
if not os.path.exists(path):
    Y = np.lib.format.open_memmap(path, mode='w+', dtype=np.float32,
                                  shape=(len(labels), pre + post, len(cand)))
    np.save(H + '/surr_starts.npy', cand)
    open(H + '/progress_s.txt', 'w').write('0')
else:
    Y = np.lib.format.open_memmap(path, mode='r+')
    cand = np.load(H + '/surr_starts.npy')
done = int(open(H + '/progress_s.txt').read().strip())
t0 = time.time()
while done < len(labels) and time.time() - t0 < 140:
    d = np.asarray(ms.read_ts_channels_sample([labels[done]], [0, nsamp]), dtype=np.float32).ravel()
    for i, s in enumerate(cand):
        Y[done, :, i] = d[s - pre:s + post]
    done += 1
    if done % 40 == 0:
        Y.flush(); open(H + '/progress_s.txt', 'w').write(str(done)); print(done, f'{time.time()-t0:.0f}s', flush=True)
Y.flush(); open(H + '/progress_s.txt', 'w').write(str(done))
print('progress', done, '/', len(labels), f'{time.time()-t0:.0f}s', flush=True)
