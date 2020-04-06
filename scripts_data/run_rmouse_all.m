% all wt, including baseline recordings from atropine experiments 
% but excluding atropine recordings 
% ** note there is an overlap with doall_wtAtropineCtrlRedundant **

global WP DS AP
% learn data root path from rmouse ini routine - has to be pulled into
% separate variable because WP is set to [] with rmouse
rmouse_ini;
dataRootPath = WP.rootPath;

% script root path 
script_root_path = "d:\hh\projects_programming\ten-years\ReScience_10yrReproChallenge\scripts_data";

tic
%% ----------------------- wt0001_04708 ------------------------------
expName = "wt0001_04708";
cd(fullfile(dataRootPath, 'beta3_wtko', expName));

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
% (000 is corrupt)
run(fullfile(script_root_path, expName, "a001_exc1.m"));
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a002_exc1.m"));
rmouse;

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a003_r1.m"));
rmouse;


%% --------------------- wt0002_04707 --------------------------------
expName = "wt0002_04707";
cd(fullfile(dataRootPath, 'beta3_wtko', expName));

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a001_r1.m"));
rmouse;

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a002_r1.m"));
rmouse;




%% --------------------- wt0003_04730 --------------------------------
expName = "wt0003_04730";
cd(fullfile(dataRootPath, 'beta3_wtko', expName));

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a003_r1.m"));
rmouse;

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a005_r1.m"));
rmouse;


% rmouse_ini;
% cd([dataRootPath '\beta3_wtko\wt2642_03211']);
% AP=[]; DS=[];
% AP_wt_atropine;
% AP_job1;  
% % NOT in atropine project
% r1;
% load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
% rmouse('af_tsl_ext',evt.tsl{1});
% % clear; close all
% 
% rmouse_ini;
% cd([dataRootPath '\beta3_wtko\wt2654_03305']);
% AP=[]; DS=[];
% AP_wt_atropine;
% AP_job1;  
% % NOT in atropine project
% r1;
% load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
% rmouse('af_tsl_ext',evt.tsl{1});
% % clear; close all
% 
% rmouse_ini;
% cd([dataRootPath '\beta3_wtko\wt2009_02319']);
% AP=[]; DS=[];
% AP_wt_atropine;
% AP_job1;  
% % NOT in atropine project
% r1;
% load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
% rmouse('af_tsl_ext',evt.tsl{1});
% % clear; close all
% 
% rmouse_ini;
% cd([dataRootPath '\beta3_wtko\wt2008_02220']);
% AP=[]; DS=[];
% AP_wt_atropine;
% AP_job1;  
% % NOT in atropine project
% r1;
% load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
% rmouse('af_tsl_ext',evt.tsl{1});
% % clear; close all

toc