%% Use FreeSurfer and Fieltrip from MATLAB to generate surface files for 
% visualization of iEEG electrodes (and more) with the iEEGLAB plugin. 
% 
% REQUIREMENTS:
%   - Install XQuartz: https://www.xquartz.org/
%   - Install FreeSurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
%   - Download/clone Fiedtrip and add path to MATLAB
% 
% Author: Cedric Cannard, adapted from Fieldtrip Tutorial. 

clear; close all; clc
eeglab; close
plugin_path = fileparts(which('eegplugin_ieeglab'));
filepath = fullfile(plugin_path, 'tutorial', 'seeg_set');
cd(filepath)

% Fieldtrip defaults
ft_defaults

% Subject folder and ID
subdir = filepath;  % parent folder only
subid  = 'sub-02';

% Create directories
orig_dir = fullfile(subdir, 'freesurfer', 'orig');
if ~exist(orig_dir, 'dir'); mkdir(orig_dir); end

% 1) Import the anatomical MRI into the MATLAB workspace using ft_read_mri
sub_files = {dir(fullfile(filepath)).name}';
mri_filename = contains(lower(sub_files), 't1');
if ~any(mri_filename)
        error("No T1 file detected in subject folder")
end
mri_filename = sub_files{mri_filename};
fprintf("T1 MRI file succesfully detected in subject folder: %s \n", t1_nii)
% mri = ft_read_mri([subjID '_MR_acpc.nii']); % we used the dcm series
mri = ft_read_mri(mri_filename);


% 2) Determine the native orientation of the anatomical MRI's left-right 
% axis using ft_determine_coordsys.
% CRITICAL STEP To correctly fuse the MRI and CT scans at a later step, 
% accuracy in demarcating the right hemisphere landmark in the following 
% step is important for avoiding an otherwise hard to detect flip of the 
% scan's left and right orientation.
mri = ft_determine_coordsys(mri);

% 3) Align the anatomical MRI to the ACPC coordinate system, a preferred 
% convention for the FreeSurfer operation optionally used in a later step. 
% In this coordinate system, the origin (coordinate [0,0,0]) is at the 
% anterior commissure (AC), the Y-axis runs along the line between the 
% anterior and posterior commissure (PC), and the Z-axis lies in the midline
%  dividing the two cerebral hemispheres. Specify the anterior and posterior
%  commissure, an interhemispheric location along the midline at the top of 
% the brain, and a location in the brain's right hemisphere.
% If the scan was found to have a left-to-right orientation in the previous
%  step, the right hemisphere is identified as the hemisphere having larger 
% values along the left-right axis. Vice versa, in a right-to-left system, 
% the right hemisphere has smaller values along that axis than its left
%  counterpart 
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);
% GUIDE: https://wiki.besa.de/index.php?title=Marking_AC-PC_Points_in_BESA_MRI#Marking_the_AC_and_PC_points_in_BESA_MRI


% 4) Write the preprocessed anatomical MRI out to file
cfg           = [];
cfg.filename  = sprintf('%s_MR_acpc',subid);
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);

% 5) Cortical surface extraction with FreeSurfer (Optional)
% Execute FreeSurfer's recon-all functionality from the Linux or 
% MacOS terminal (Windows via VirtualBox), or from the MATLAB command 
% window as below. This set of commands will create a folder named 
% 'freesurfer' in the subject directory, with subdirectories containing a
%  multitude of FreeSurfer-generated files.
% 
% 
fshome      = '/Applications/freesurfer/8.1.0';
% fslicense  = '/Users/cedriccannard/Documents/license.txt';
mrfile      = sprintf('%s.nii', cfg.filename);

setenv('FREESURFER_HOME', fshome);
setenv('SUBJECTS_DIR', subdir);
% setenv('FS_LICENSE', fslicense);
setenv('PATH', [fullfile(FREESURFER_HOME,'bin') ':' fullfile(FREESURFER_HOME,'fsfast','bin') ':' getenv('PATH')]);
nThreads = max(1, feature('numcores')-1);

% FreeSurfer's fully automated segmentation and cortical extraction of a
% raw T1-weighted MRI (CAN TAKE 10 hours!!!!!)
% 
% tic
% system(['export FREESURFER_HOME=' fshome '; ' ...
% 'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
% 'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
% 'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])
% toc; gong

% To do minimal surface reconstruction only (using parallel CPU cores (optional). 
% -autorecon1 = preprocessing (motion correction, Talairach, skull strip, etc.)
% -autorecon2-wm = white matter surface reconstruction (creates surf/lh.white, surf/rh.white)
% -autorecon3 = parcellations, thickness maps/metrics, surfaces on sphere, and stats
% here, ensure the brainmask.mgz file is built after autorecon1
tic
% mridir = fullfile(subdir,'freesurfer','mri');
% system(['export FREESURFER_HOME=' fshome '; ' ...
% 'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
% 'export SUBJECTS_DIR=' subdir '; ' ...
% 'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
% 'recon-all -i ' [subdir '/tmp.nii'] ' -s freesurfer -sd ' subdir ' -autorecon1 -openmp ' num2str(nThreads) '; ' ...
% 'cd "' mridir '"; $FREESURFER_HOME/bin/mri_mask -T 5 nu.mgz brainmask.mgz brain.mgz; cd - >/dev/null; ' ...
% 'recon-all -s freesurfer -sd ' subdir ' -autorecon2 -openmp ' num2str(nThreads)]);

% -autorecon1 → produces nu.mgz, brainmask.mgz; we also create brain.mgz if it's missing (your earlier crash happened right after the BFS/masking step).
% -autorecon2-volonly → creates aseg.presurf.mgz. If this fails, you'll see it immediately in the log.
% -autorecon2-wm/-autorecon2-pial → builds white and pial surfaces. Stats/parcellations are skipped to avoid unrelated failures.
% Final ls shows whether lh.white, rh.white, lh.pial, rh.pial exist in .../freesurfer/surf.
tic
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'export SUBJECTS_DIR="' subdir '"; ' ...
'export PATH="$FREESURFER_HOME/bin:$FREESURFER_HOME/fsfast/bin:$PATH"; ' ...
'mri_convert -c -oc 0 0 0 "' mrfile '" "' [subdir '/tmp.nii'] '"; ' ...
'recon-all -i "' [subdir '/tmp.nii'] '" -s freesurfer -sd "' subdir '" -autorecon1 -openmp ' num2str(nThreads) '; ' ...
'cd "' [subdir '/freesurfer/mri'] '"; ' ...
'if [ ! -f brain.mgz ] && [ -f nu.mgz ] && [ -f brainmask.mgz ]; then mri_mask -T 5 nu.mgz brainmask.mgz brain.mgz; fi; ' ...
'cd - >/dev/null; ' ...
'export OMP_NUM_THREADS=1; export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1; ' ...  % stabilize seg2cc/mri_cc
'recon-all -s freesurfer -sd "' subdir '" -autorecon2-volonly -openmp 1 -noparcstats -noparcstats2 -noparcstats3 -nobalabels; ' ...
'unset OMP_NUM_THREADS; unset ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS; ' ...
'recon-all -s freesurfer -sd "' subdir '" -autorecon2 -openmp ' num2str(nThreads) ' -noparcstats -noparcstats2 -noparcstats3 -nobalabels; ' ...
'ls -l "' [subdir '/freesurfer/surf'] '" ']);
toc; gong



% 6) Import the extracted cortical surfaces into the MATLAB workspace and 
% examine their quality. Repeat the following code using rh.pial to visualize
%  the pial surface of the right hemisphere.
pial_lh = ft_read_headshape('freesurfer/surf/lh.pial');
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud;
camlight


%% Quick FreeSurfer sanity checks

% % Show version and where binaries come from
% system('recon-all --version');
% system('which recon-all');
% system(sprintf('mri_info "%s" | head -n 10', t1_mgz));
% 
% fprintf('\nDone. Subject folder: %s\n', fullfile(subdir,subid));


%% Preprocessing of the anatomical CT


% 9) Import the anatomical CT into the MATLAB workspace. Similar to the MRI,
%  the CT scan comes in the format of a single file with an .img or .nii 
% extension, or a folder containing a series of files with a .dcm or .ima 
% extension
ct = ft_read_mri([subjID '_CT_acpc_f.nii']); % we used the dcm series

% 10) In case this cannot be done on the basis of knowledge of the 
% laterality of electrode implantation, determine the native orientation of
%  the anatomical CT's left- right axis using ft_determine_coordsys, 
% similarly to how it was done with the anatomical MRI in Step 3.
% CRITICAL STEP To correctly fuse the MRI and CT scans at a later step, 
% accuracy in demarcating the right and left preauricular landmark in the 
% following step is important for avoiding an otherwise hard to detect flip 
% of the scan's left and right orientation.
ct = ft_determine_coordsys(ct);


% 11) Align the anatomical CT to the CTF head surface coordinate system by
% specifying the nasion (at the root of the nose), left and right preauricular 
% points (just in front of the ear canals), and an interhemispheric location 
% along the midline at the top of the brain. The CT scan is initially 
% aligned to the CTF head surface coordinate system given that the ACPC 
% coordinate system used for the MRI relies on neuroanatomical landmarks 
% that are not visible in the CT.
cfg           = [];
cfg.method    = 'fiducial'; % 'interactive' 'fiducial'
cfg.coordsys  = 'ctf';
cfg.viewresult= 'yes';
ct_ctf = ft_volumerealign(cfg, ct);


% 12) Automatically convert the CT's coordinate system into an approximation
%  of the ACPC coordinate system, the same system the anatomical MRI was 
% aligned to. The call to ft_hastoolbox ensures that the required SPM 
% toolbox is on the path.
ft_hastoolbox('spm12', 1);
ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');

%% Fusion of the CT with the MRI

% 8) Import the FreeSurfer-processed MRI into the MATLAB workspace for 
% later fusion with the CT scan, and specify the coordinate system to which
% it was aligned.
fsmri_acpc = ft_read_mri('freesurfer/mri/T1.mgz'); % on Windows, use 'SubjectUCI29_MR_acpc.nii'
fsmri_acpc.coordsys = 'acpc';

% 13) Fuse the CT with the MRI, a necessary step to link the electrode 
% locations in the anatomical CT to their corresponding locations in the 
% anatomical MRI. Given that both scans are from the same subject and their 
% common denominator is the skull, a rigid body transformation suffices for 
% their alignment under normal circumstances (the default technique when
%  using the SPM-method in FieldTrip).
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'acpc';
cfg.viewresult  = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);

% 14) Carefully examine the interactive figure that is produced after the 
% coregistration is completed, showing the MRI and fused CT superimposed. 
% A successful fusion will show tight interlocking of CT-positive skull (in blue) and MRI-positive brain and skin tissue (in red).
% CRITICAL STEP Accuracy of the fusion operation is important for correctly 
% placing the electrodes in anatomical context in a following step.


% 15) Write the MRI-fused anatomical CT out to file
cfg           = [];
cfg.filename  = [subjID '_CT_acpc_f'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);


%% Electrode placement

% 16) Import the header information from the recording file, if possible. 
% By giving the electrode labels originating from the header as input to 
% ft_electrodeplacement in the next step, the labels will appear as a to-do
%  list during the interactive electrode placement activity. A second 
% benefit is that the electrode locations can be directly assigned to labels
%  collected from the recording file, obviating the need to sort and rename
%  electrodes to match the electrophysiological data.

load([subjID '_hdr.mat']);
hdr = ft_read_header(hdr);



