import os, sys, time, numpy as np, pandas as pd
from pymef.mef_session import MefSession
H = os.path.expanduser('~')
V = os.path.expanduser('~/mnt/v1.0.0/sub-02/ses-ieeg01/ieeg')
B = 'sub-02_ses-ieeg01_task-ccep_run-01'
ms = MefSession(os.path.join(V, B + '_ieeg.mefd'), None)
info = {c['name']: c for c in ms.read_ts_channel_basic_info()}
ch = pd.read_csv(os.path.join(V, B + '_channels.tsv'), sep='\t')
seeg = ch[ch.type.isin(['SEEG', 'ECOG'])]
keep = [c for c in seeg.name.tolist() if c in info]
ev = pd.read_csv(os.path.join(V, B + '_events.tsv'), sep='\t')
ev = ev[(ev.trial_type == 'electrical_stimulation') & (ev.status == 'good')].reset_index(drop=True)
srate = float(info[keep[0]]['fsamp'][0]); nsamp = int(info[keep[0]]['nsamp'][0])
pre, post = int(round(0.5 * srate)), int(round(1.0 * srate))
starts = ev.sample_start.values.astype(int)
ok = (starts - pre >= 0) & (starts + post < nsamp)
ev, starts = ev[ok].reset_index(drop=True), starts[ok]
shape = (len(keep), pre + post, len(starts))
path = os.path.join(H, 'sub02_epochs.npy')
if not os.path.exists(path):
    X = np.lib.format.open_memmap(path, mode='w+', dtype=np.float32, shape=shape)
    np.savez(os.path.join(H, 'sub02_meta.npz'), labels=np.array(keep),
             status=np.array(seeg.set_index('name').loc[keep, 'status'].values, dtype=str),
             srate=srate, sites=ev.electrical_stimulation_site.values.astype(str),
             times=np.arange(-pre, post) / srate, starts=starts)
    open(os.path.join(H, 'progress.txt'), 'w').write('0')
else:
    X = np.lib.format.open_memmap(path, mode='r+')
done = int(open(os.path.join(H, 'progress.txt')).read().strip())
t0 = time.time()
while done < len(keep) and time.time() - t0 < 140:
    d = np.asarray(ms.read_ts_channels_sample([keep[done]], [0, nsamp]), dtype=np.float32).ravel()
    for i, s in enumerate(starts):
        X[done, :, i] = d[s - pre:s + post]
    done += 1
    if done % 20 == 0:
        X.flush(); open(os.path.join(H, 'progress.txt'), 'w').write(str(done))
        print(f'  {done}/{len(keep)} {time.time()-t0:.0f}s', flush=True)
X.flush(); open(os.path.join(H, 'progress.txt'), 'w').write(str(done))
print(f'progress {done}/{len(keep)}  shape {shape}  {time.time()-t0:.0f}s', flush=True)
