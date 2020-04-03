% all wt, including baseline recordings from atropine experiments 
% but excluding atropine recordings 
% ** note there is an overlap with doall_wtAtropineCtrlRedundant **

global WP DS AP
% learn data root path from rmouse ini routine
rmouse_ini;

% add path because various commonly used ANPAR files resides there
addpath([WP.rootPath '\beta3_wtko']);
% -----------------------------------------------------

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% same as for atropine project
% (000 is corrupt)
a001_exc1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0002_04707']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% same as for atropine project
% (alternatively: a001_r1)
a000_r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0003_04730']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% same as for atropine project
% (alternatively: a002_r1 (has very little immobility, though))
a003_r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt2642_03211']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% NOT in atropine project
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt2654_03305']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% NOT in atropine project
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt2009_02319']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% NOT in atropine project
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt2008_02220']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
% NOT in atropine project
r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

% the rejected ones:
if 0
  cd d:\rmouse\wt2141_01d18
  AP=[]; DS=[];
  r1;
  AP_beta3_wtko;
  AP_job1;  
  rmouse;
  % clear; close all
  
  cd d:\rmouse\wt2141_02510
  AP=[]; DS=[];
  r1;
  AP_beta3_wtko;
  AP_job1;  
  rmouse;
  % clear; close all
  
  cd d:\rmouse\wt2001_02424
  AP=[]; DS=[];
  r1;
  AP_beta3_wtko;
  AP_job1;  
  rmouse;
  % clear; close all

  cd d:\rmouse\wt2206_02725
  AP=[]; DS=[];
  r1;
  AP_beta3_wtko;
  AP_job1;  
  rmouse;
  % clear; close all
  
end

rmouse_ini;
cd([WP.rootPath '\beta3_wtko']);