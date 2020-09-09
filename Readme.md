# Reproduction of phase-amplitude coupling analysis of multi-electrode field potential data from rodent hippocampus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4021389.svg)](https://doi.org/10.5281/zenodo.4021389)

This is a contribution to the ["Ten Years Reproducibility Challenge"]( https://github.com/ReScience/ten-years ) run by [ReScience C](https://rescience.github.io/ ).

Original article: Hentschke, H., Perkins, M.G., Pearce, R.A., & Banks, M.I. (2007) Muscarinic blockade weakens interaction of gamma with theta rhythms in mouse hippocampus. Eur.J.Neurosci., 26, 1642–1656.

## Contents

__paper__ contains the ReScience submission (Latex source and pdf), including figures.

__rmouse__ is the core of the code presented in the ReScience submission. It is a Matlab toolbox for time-and frequency domain analyses of multi-channel local field potential (LFP) data from rodent hippocampus, with a focus on phase-amplitude coupling of theta with gamma rhythms. It has been designed for data recorded from awake animals showing a diverse range of behaviors, and therefore expects as inputs not only the neural recordings but also behavioral scoring data.

The 'main' file rmouse.m digests one data file at a time, and saves the results of the computations into native Matlab data files (*.mat), one per raw data file. It can optionally create and save graphics for most types of analysis, including a condensed summary figure, for each file. For data aggregation across experiments and statistics, a range of post-processing data collection and processing routines are included.

rmouse requires about 70 (seventy) parameters to run the analyses. These are contained in two structs termed DS (DataSet) and AP (Analysis Parameters). They must be global variables, so need not be specified as inputs.

The code has been restructured and cleaned up to a minor degree, and basic functionality with Matlab R2020a has been verfied. To use the code, please consult the following resources:

* rmouse/manual_rmouse.pdf - a general manual explaining the concept and evolution of the toolbox in depth, including the various analyses, post-processing routines, nomenclature of variables and channels, as well as a section on behavioral scoring.
* rmouse/templates - here, you will find templates of scripts which set up above-mentioned global variables DS and AP. All parameters are explained via comments. These files are essential starting points for new analyses.
* the code in subdirectories scripts_data and scripts_paper for concrete use cases
* the original article for an example of the kinds of results produced by rmouse
* the ReScience article for a fresh look at and personal evaluation of the code by its author.

__scripts_data__ contains the scripts for running rmouse which have been originally set up for the analysis of the data in Hentschke et al. (2007). The five subdirectories starting with 'wt' contain scripts specific to the experimental animals originally analyzed. See the Readme file.

__scripts_paper__ contains scripts required for reproducing figures with key findings for the ReScience submission, plus detailed instructions on how to accomplish this.

