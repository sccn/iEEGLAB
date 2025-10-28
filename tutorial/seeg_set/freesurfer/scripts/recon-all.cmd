\n\n#---------------------------------
# New invocation of recon-all Mon Oct 27 15:00:30 PDT 2025 
\n mri_convert /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/tmp.nii /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig/001.mgz \n
#--------------------------------------------
#@# MotionCor Mon Oct 27 15:00:31 PDT 2025
\n cp /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig/001.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz \n
\n mri_info /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz \n
\n mri_convert /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/rawavg.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --conform \n
\n mri_add_xform_to_header -c /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/transforms/talairach.xfm /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz \n
\n mri_info /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz \n
\n mri_synthstrip --threads 9 -i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz -o /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/synthstrip.mgz \n
\n mri_synthseg --i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --o /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/synthseg.rca.mgz --threads 9 --vol /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/stats/synthseg.vol.csv --keepgeom --addctab --cpu \n
\n fs-synthmorph-reg --s freesurfer --threads 9 --i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --test \n
#--------------------------------------------
#@# Nu Intensity Correction Mon Oct 27 15:13:22 PDT 2025
\n mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 --ants-n4 \n
\n mri_add_xform_to_header -c /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/transforms/talairach.xfm nu.mgz nu.mgz \n
#--------------------------------------------
#@# Intensity Normalization Mon Oct 27 15:14:30 PDT 2025
\n mri_normalize -g 1 -seed 1234 -mprage nu.mgz T1.mgz \n
#--------------------------------------
\n#@# MCADura Segmentation Mon Oct 27 15:15:48 PDT 2025
#--------------------------------------
\n#@# VSinus Segmentation Mon Oct 27 15:16:03 PDT 2025
\n mri_mask /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/T1.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/synthstrip.mgz /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/brainmask.mgz \n
#--------------------------------------
\n#@# EntoWM Segmentation Mon Oct 27 15:16:16 PDT 2025
\n\n#---------------------------------
# New invocation of recon-all Mon Oct 27 15:16:25 PDT 2025 
\n fs-synthmorph-reg --s freesurfer --threads 1 --i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --test \n
#--------------------------------------
\n#@# MCADura Segmentation Mon Oct 27 15:16:27 PDT 2025
MCADura Segmentation does not need to be updated
#--------------------------------------
\n#@# VSinus Segmentation Mon Oct 27 15:16:27 PDT 2025
VSinus Segmentation does not need to be updated
#-------------------------------------
#@# EM Registration Mon Oct 27 15:16:27 PDT 2025
\n mri_em_register -uns 3 -mask brainmask.mgz nu.mgz /Applications/freesurfer/8.1.0/average/RB_all_2020-01-02.gca transforms/talairach.lta \n
#--------------------------------------
#@# CA Normalize Mon Oct 27 15:22:40 PDT 2025
\n mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /Applications/freesurfer/8.1.0/average/RB_all_2020-01-02.gca transforms/talairach.lta norm.mgz \n
#--------------------------------------
\n#@# EntoWM Segmentation Mon Oct 27 15:23:25 PDT 2025
EntoWM Segmentation does not need to be updated
#--------------------------------------
#@# CC Seg Mon Oct 27 15:23:25 PDT 2025
\n seg2cc --s freesurfer \n
\n\n#---------------------------------
# New invocation of recon-all Mon Oct 27 15:23:26 PDT 2025 
\n fs-synthmorph-reg --s freesurfer --threads 9 --i /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/seeg_set/freesurfer/mri/orig.mgz --test \n
#--------------------------------------
\n#@# MCADura Segmentation Mon Oct 27 15:23:27 PDT 2025
MCADura Segmentation does not need to be updated
#--------------------------------------
\n#@# VSinus Segmentation Mon Oct 27 15:23:27 PDT 2025
VSinus Segmentation does not need to be updated
#-------------------------------------
#@# EM Registration Mon Oct 27 15:23:27 PDT 2025
\n mri_em_register -uns 3 -mask brainmask.mgz nu.mgz /Applications/freesurfer/8.1.0/average/RB_all_2020-01-02.gca transforms/talairach.lta \n
#--------------------------------------
#@# CA Normalize Mon Oct 27 15:29:14 PDT 2025
\n mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /Applications/freesurfer/8.1.0/average/RB_all_2020-01-02.gca transforms/talairach.lta norm.mgz \n
#--------------------------------------
\n#@# EntoWM Segmentation Mon Oct 27 15:29:59 PDT 2025
EntoWM Segmentation does not need to be updated
#--------------------------------------
#@# CC Seg Mon Oct 27 15:29:59 PDT 2025
\n seg2cc --s freesurfer \n
