function run_wt0003_04730
% Processing of data files of wt0003 with rmouse. Use for batch().

global WP DS AP
% learn data root path from rmouse ini routine - has to be pulled into
% separate variable because WP is set to [] with rmouse
rmouse_ini;
dataRootPath = WP.rootPath;

% script root path 
script_root_path = "d:\hh\projects_programming\ten-years\ReScience_10yrReproChallenge\scripts_data";



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
