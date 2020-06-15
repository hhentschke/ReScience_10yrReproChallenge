# Instructions for running raw data analyses with rmouse for the Ten Years Reproducibility Challenge paper

This directory contains scripts for running rmouse which have been originally set up for the analysis of the data in Hentschke et al. (2007). Scripts specific to the five experimental animals originally analyzed are distributed into the subdirectories starting with 'wt'.

Raw data and behavioral scoring files are not part of this repository, but can be furnished upon request. 

1. Open script run_rmouse_all.m. This is the 'master script'. Also open the scripts it calls (AP_job1, AP_wt_atropine, as well as the experiment-specific scripts in the subfolders 'wt...'). Set analysis parameters according to your needs and adjust the paths.
2. Run script run_rmouse_all.m to generate processed data for each of the recordings listed in it via the rmouse toolbox. run_rmouse_all_parallel.m does the same in an embarassingly parallel way, using parfeval (Parallel Computing Toolbox needed) and the run_wt00... scripts.

For data aggregation, statistics, and specific figures, please see /scripts_paper.
