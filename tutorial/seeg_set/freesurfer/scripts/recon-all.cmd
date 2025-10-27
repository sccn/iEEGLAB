\n\n#---------------------------------
# New invocation of recon-all Mon Oct 27 11:55:45 PDT 2025 
\n mri_convert /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/tmp.nii /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig/001.mgz \n
#--------------------------------------------
#@# MotionCor Mon Oct 27 11:55:46 PDT 2025
\n cp /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig/001.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz \n
\n mri_info /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz \n
\n mri_convert /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --conform \n
\n mri_add_xform_to_header -c /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/transforms/talairach.xfm /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz \n
\n mri_info /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz \n
\n mri_synthstrip --threads 9 -i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz -o /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/synthstrip.mgz \n
