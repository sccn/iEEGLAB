# iEEGLAB (in development)

EEGLAB plugin for analyzing intracranial EEG (iEEG) data. 

EEGLAB plugin for analyzing intracranial EEG (iEEG) data. Supports both stereoEEG (sEEG) and eCoG data. The plugin supports continuous iEEG data applications:

- Epilepsy research

- Clinical monitoring


Although it is mainly designed for event-related applications:

- Stimulus-/Response-locked

- Cortico-Cortical Evoked Potentials (CCEP; e.g., single pulse stimulation experiments)


<img src="https://raw.githubusercontent.com/amisepa/iEEGLAB/main/tutorial/images/6_vis_elecs1.png" width="40%" /> <img src="https://raw.githubusercontent.com/amisepa/iEEGLAB/main/tutorial/images/6_vis_elecs2.png" width="40%" />
<img src="https://raw.githubusercontent.com/amisepa/iEEGLAB/main/tutorial/images/6_vis_elecs3.png" width="50%" />
<p align="center">
  <img src="https://raw.githubusercontent.com/amisepa/iEEGLAB/main/tutorial/images/6_vis_elecs3_opt.gif" width="70%">
  <br>
  <em>Interactive 3D rotation of the glass brain and electrodes</em>
</p>


A lot of the code and algorithms implemented in this plugin were adapted from work by Dora Hermes and the Multimodal Neuroimaging Lab (https://github.com/MultimodalNeuroimagingLab). Please cite the following references when using this plugin:

Valencia, G. O.et al., (2023). Signatures of electrical stimulation driven network interactions in the human limbic system. Journal of Neuroscience, 43(39), 6697-6711. https://pubmed.ncbi.nlm.nih.gov/37620159/

Huang, H., Valencia, G. O., Gregg, N. M., Osman, G. M., Montoya, M. N., Worrell, G. A., ... & Hermes, D. (2024). CARLA: Adjusted common average referencing for cortico-cortical evoked potential data. Journal of neuroscience methods, 407, 110153. https://pubmed.ncbi.nlm.nih.gov/38710234/

Miller, K. J., Müller, K. R., Valencia, G. O., Huang, H., Gregg, N. M., Worrell, G. A., & Hermes, D. (2023). Canonical Response Parameterization: Quantifying the structure of responses to single-pulse intracranial electrical brain stimulation. PLoS computational biology, 19(5), e1011105. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011105


## Requirements

- MATLAB installed with license

- EEGLAB Toolbox installed and path added to MATLAB: https://github.com/sccn/eeglab?tab=readme-ov-file#installingcloning

- Vistasoft for 3D glass brain visualizations: clone/download repository (https://github.com/vistalab/vistasoft.git) and add the path in MATLAB manually (see above)

- Data importation plugins to import your iEEG data in EEGLAB (depends on the data format; e.g., .mefd, .edf, .vhdr, etc.; See tutorial for more details).

- Have 3D cartesian (XYZ) electrode locations for each file you wish to analyze (in .tsv/.csv file)

- Have events either directly in the data or in a .tsv/.csv file (for event-related applications)


## Installing the iEEGLAB plugin: https://github.com/amisepa/iEEGLAB/wiki#1-install-the-ieeglab-plugin


## Tutorial

Go to the wiki page: https://github.com/amisepa/iEEGLAB/wiki


