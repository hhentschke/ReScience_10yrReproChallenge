# Instructions for generating raw figures for the Ten Years Reproducibility Challenge paper

1. Open scripts run_rmouse_all.m and the scripts it calls (AP_job1, AP_wt_atropine, as well as the experiment-specific scripts in the subfolders 'wt...'). Set analysis parameters according to your needs. Note that data could only be retrieved for three out of five experimental subjects, which is why the 
2. Run script run_rmouse_all.m to generate processed data for each of the listed recordings via the rmouse toolbox. run_rmouse_all_parallel.m does the same in an embarassingly parallel way, using parfeval (Parallel Computing Toolbox needed) and the run_wt00... scripts.
3. Edit function set_rv.m and include all quantities of interest in variable rv_commented. set_rv will be default be called by combine_r (see step 5)
4. Call script collect_APDS_all.m, which will create variable RInfo and global variables ANPAR and DSET in the workspace. They are required in the following step.
5. Assemble processed data into one file by calling function combine_r, e.g. like so:
combine_r(RInfo, 'd:\_data\rmouse\ten-years\beta3_wtko\intermediate_data\out_combine_r.mat').

### Figure 3A (power spectra of the raw data and the gamma envelope)
Edit script plot_spectra.m, possibly adjust paths and a few other settings to taste, and run.

### Figure 3B (laminar 'Christmas tree' plot of the power of the gamma envelope at theta frequencies)
Call function rdeal with the corresponding set of input args, particularly rv set to {'gaeThPEMn_auto'}. It is strongly recommended to edit and run script runRdeal (otherwise it is nigh to impossible to figure out the correct input args)

### Figure 3C (laminar 'Christmas tree' plot of theta-gamma peak crosscorrelation)
Call function rdeal with the corresponding set of input args, particularly rv set to {'thgaeCCPeakMn_auto'}. It is strongly recommended to edit and run script runRdeal (otherwise it is nigh to impossible to figure out the correct input args)

### Any other figure
Inspect the contents of folder 'reservoir': the code for sought-after plot is very likely in there, and will for sure have to be adapted.