# Instructions for generating figures for the Ten Years Reproducibility Challenge submission

The following instructions assume the presence of data processed by rmouse for at least two experiments (two animals, each of them recorded during control and drug conditions). In the original paper, five animals had been analyzed; for the Reproducibility Challenge, data from only three animals could be retrieved.

### Preparatory steps required for all figures which aggregate data across experiments
1. Edit function set_rv.m and include all quantities of interest in variable rv_commented. set_rv will be default be called by combine_r (see step 3)
2. Call script collect_APDS_all.m, which will create variable RInfo and global variables ANPAR and DSET in the workspace. They are required in the following step.
3. Assemble processed data into one file by calling function combine_r, e.g. like so:
combine_r(RInfo, 'd:\_data\rmouse\ten-years\beta3_wtko\intermediate_data\out_combine_r.mat').

### Figure 3A (power spectra of the raw data and the gamma envelope)
Edit script plot_spectra.m, possibly adjust paths and a few other settings to taste, and run.

### Figure 3B (laminar 'Christmas tree' plot of the power of the gamma envelope at theta frequencies)
Call function rdeal with the corresponding set of input args, particularly rv set to {'gaeThPEMn_auto'}. It is strongly recommended to edit and run script runRdeal (otherwise it is nigh to impossible to figure out the correct input args)

### Figure 3C (laminar 'Christmas tree' plot of theta-gamma peak crosscorrelation)
Call function rdeal with the corresponding set of input args, particularly rv set to {'thgaeCCPeakMn_auto'}. It is strongly recommended to edit and run script runRdeal (otherwise it is nigh to impossible to figure out the correct input args)