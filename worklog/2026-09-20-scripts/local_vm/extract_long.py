import os, sys, time, numpy as np, pandas as pd
from pymef.mef_session import MefSession
H=os.path.expanduser('~'); V=os.path.expanduser('~/mnt/v1.0.0/sub-02/ses-ieeg01/ieeg')
B='sub-02_ses-ieeg01_task-ccep_run-01'
ms=MefSession(os.path.join(V,B+'_ieeg.mefd'), None)
info={c['name']: c for c in ms.read_ts_channel_basic_info()}
meta=np.load(H+'/sub02_meta.npz', allow_pickle=True)
labels, sites, starts = meta['labels'], meta['sites'], meta['starts']
srate=float(meta['srate']); nsamp=int(info[labels[0]]['nsamp'][0])
P=pd.read_csv(H+'/crp_paired.csv'); P=P[P.z>6]
top=P.groupby('site').size().sort_values(ascending=False).head(6).index.tolist()
pairs=P[P.site.isin(top)].groupby('site').apply(lambda d: d.nlargest(12,'z').chan.tolist(), include_groups=False).to_dict()
chans=sorted({c for v in pairs.values() for c in v})
pre=post=int(round(2.0*srate))
rng=np.random.default_rng(11)
surr=[]
while len(surr) < 40:
    s=int(rng.integers(pre+1, nsamp-post-1))
    if np.min(np.abs(starts-s)) > 3*srate: surr.append(s)
surr=np.array(sorted(surr))
site_tr={s: starts[sites==s] for s in top}
out={}
t0=time.time()
for ci,ch in enumerate(chans):
    d=np.asarray(ms.read_ts_channels_sample([ch],[0,nsamp]), dtype=np.float32).ravel()
    for s in top:
        if ch not in pairs[s]: continue
        tr=site_tr[s]
        out[f'{s}|{ch}']=np.stack([d[x-pre:x+post] for x in tr], axis=1)
    out[f'SURR|{ch}']=np.stack([d[x-pre:x+post] for x in surr], axis=1)
    if ci%10==0: print(' ', ci, len(chans), '%.0fs'%(time.time()-t0), flush=True)
np.savez(H+'/long_epochs.npz', srate=srate, pre=pre,
         times=(np.arange(-pre,post)/srate), **out)
print('saved', len(out), 'arrays, %.0fs'%(time.time()-t0), 'chans', len(chans), 'sites', top, flush=True)
