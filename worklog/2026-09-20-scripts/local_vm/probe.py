from pymef.mef_session import MefSession
import os, json
V = os.path.expanduser('~/mnt/v1.0.0/sub-02/ses-ieeg01/ieeg')
p = os.path.join(V, 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd')
ms = MefSession(p, None)
info = ms.read_ts_channel_basic_info()
print('n channels', len(info))
k = info[0]
print({key: k[key] for key in ['name','fsamp','nsamp','ufact','unit','start_time','end_time'] if key in k})
print('names', [c['name'] for c in info[:8]])
