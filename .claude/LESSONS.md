# Lessons (newest last)

- MATLAB `-batch` on the Windows dev machine hangs on any figure or dialog (R2025b/R2026a,
  also with software OpenGL). Tests are compute-only; functions separate compute from
  drawing (`'draw', false`). For figures, export results to JSON and plot in Python.
- The shipped sEEG tutorial is 128 Hz: blanking and N1 cannot be judged on it (N1 locks
  onto residual artifact, >1000 uV). Use native-rate recordings for those checks.
- `run_CRP`'s `p_value_tR` is a t-test at the tau_R chosen from the same data: biased.
  `ieeglab_stats_subject` uses a permutation null that repeats the selection.
- GitHub Actions (Ubuntu) rounds single-precision data differently from Windows: scale
  tolerances with `eps(single(max(abs(data(:)))))`, not fixed absolute limits.
- GitHub job logs need a token; the workflow prints `::error` lines so failing tests are
  readable through the public check-run annotations API.
- `.surf.gii` files: do not assume `darrays[0]` holds vertices; pick the float N x 3 array.
- The tutorial extract was also missing its `channels.tsv` (added 2026-09-21 from ds004696).
  Native-rate HAPwave (MEF3, 2048 Hz, 237 contacts, 695 pulses) is unzipped at
  `C:\Users\ccann\Downloads\v1.0.0`; use it for anything involving N1, blanking, HFB or phase.
- CRP explained variance is ~0.45 on stimulation-free epochs (rank-1 fit of 8-25 short traces).
  Always report it against a matched surrogate. Selecting "responsive" pairs on the same real
  data inflates real-vs-surrogate differences; use a split-half (select odd, evaluate even).
- ITPC/phase: mind the wavelet cone of influence. With -0.5 s epochs a 4 Hz, 5-cycle wavelet
  (0.625 s half-width) reaches across the stimulus and the "baseline" measures the response.
  Use epochs of at least -2..+2 s, baseline -1.8..-1.0 s, equal trial counts vs surrogate, and
  subtract the ERP (induced ITPC) before claiming phase reset. Chance ITPC ~ sqrt(pi)/2/sqrt(n).
- erdetect's rule is 3.4 x max(baseline SD, 50 uV). Its default baseline (-1, -0.1) s returns
  zero detections silently on epochs shorter than 1 s pre-stimulus.
- Clinician-bad contacts are bimodal: quiet (out of brain) and noisy. Variance criteria catch the
  noisy ones only (12/14 vs 0/25 on HAPwave sub-02). channels.tsv stays the source of truth.
- Any "does X change ICA" claim needs a seed-to-seed control: two ICA runs on identical data
  agree only at median |r| ~0.90 with ~200 channels.
- Cowork's local Linux workspace kills background jobs when each call ends and caps calls at
  180 s; chunk long jobs with a memmap plus a progress file. It cannot reach openneuro.org, S3
  or NEMAR, and cannot run MATLAB.
