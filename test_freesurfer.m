%% FreeSurfer run from MATLAB on macOS (M1)
% Subject: UCI29
% Data folder:
%   /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/dataset5
%
% This script
% 1. Sets environment variables for FreeSurfer inside MATLAB
% 2. Converts T1 NIfTI to 001.mgz in $SUBJECTS_DIR/<SUB>/mri/orig
% 3. Kicks off recon-all (autorecon1 first)

clear; close all; clc
eeglab; close
plugin_path = fileparts(which('eegplugin_ieeglab'));
filepath = fullfile(plugin_path, 'tutorial', 'seeg_set');
cd(filepath)

% Fieldtrip defaults
ft_defaults

SUBJECT_DIR = filepath;  % parent folder only
SUB_ID      = 'sub-02';

% Create directories
orig_dir = fullfile(SUBJECT_DIR, 'freesurfer', 'orig');
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
cfg.filename  = [SUB_ID '_MR_acpc'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);

% 5) Cortical surface extraction with FreeSurfer (Optional)
% Requirements:
%   - Install XQuartz: https://www.xquartz.org/
%   - Install FreeSurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
% 
% Execute FreeSurfer's recon-all functionality from the Linux or 
% MacOS terminal (Windows via VirtualBox), or from the MATLAB command 
% window as below. This set of commands will create a folder named 
% 'freesurfer' in the subject directory, with subdirectories containing a
%  multitude of FreeSurfer-generated files.
% 
% WARNING: FreeSurfer's fully automated segmentation and cortical
%   extraction of a T1-weighted MRI can take 10 hours!!!
% 
% FREESURFER_HOME = '/Applications/freesurfer/8.1.0';
% FS_LICENSE      = '/Users/cedriccannard/Documents/license.txt';
% setenv('FREESURFER_HOME', FREESURFER_HOME);
% setenv('SUBJECTS_DIR', SUBJECT_DIR);
% setenv('FS_LICENSE', FS_LICENSE);
% setenv('FSF_OUTPUT_FORMAT','nii.gz');
% setenv('PATH', [fullfile(FREESURFER_HOME,'bin') ':' fullfile(FREESURFER_HOME,'fsfast','bin') ':' getenv('PATH')]);
fshome     = '/Applications/freesurfer/8.1.0';
subdir     = SUBJECT_DIR;
mrfile     = sprintf('%s.nii', cfg.filename);
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])
gong


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
% fprintf('\nDone. Subject folder: %s\n', fullfile(SUBJECT_DIR,SUB_ID));
