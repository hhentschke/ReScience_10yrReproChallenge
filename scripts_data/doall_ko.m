global WP DS AP

% learn data root path from rmouse ini routine
rmouse_ini;

% add path because various commonly used ANPAR files resides there
addpath([WP.rootPath '\beta3_wtko']);

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2445_03214']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
% load artifact time stamps
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2472_02o11']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
a000_r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2483_03303']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2484_02927']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
a002_r1;
% ** gamma (and theta) CC lags get derailed with the standard settings - better use
% restricted intervals 
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko3362_02528']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko3363_02313']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko3382_02801']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all


rmouse_ini;
cd([WP.rootPath '\beta3_wtko']);