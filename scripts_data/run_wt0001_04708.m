function run_wt0001_04708
% Processing of data files of wt0001 with rmouse. Use for batch().

global WP DS AP
% learn data root path from rmouse ini routine - has to be pulled into
% separate variable because WP is set to [] with rmouse
rmouse_ini;
dataRootPath = WP.rootPath;

% script root path 
script_root_path = "d:\hh\projects_programming\ten-years\ReScience_10yrReproChallenge\scripts_data";


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

% AP=[]; 
% DS=[];
% AP_wt_atropine;
% AP_job1;
% run(fullfile(script_root_path, expName, "a002_exc1.m"));
% rmouse;

AP=[]; 
DS=[];
AP_wt_atropine;
AP_job1;
run(fullfile(script_root_path, expName, "a003_r1.m"));
rmouse;

