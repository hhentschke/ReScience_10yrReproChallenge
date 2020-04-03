global WP AP DS
% learn data root path from rmouse ini routine
rmouse_ini;

% add path because various commonly used ANPAR files resides there
addpath([WP.rootPath '\beta3_wtko']);

% -----------------------------------------------------
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
for fi=[1 2]
  AP=[]; DS=[];
  eval(['a0' sprintf('%.2i',fi) '_exc1;']);
  AP_beta3_wtko;
  AP_job2;  
  rmouse('callJob','seg_thetaTrigplot'); 
  clear; close all; pack
end
% file(s) designated for comparison control vs. drug
% control: 001
% drug: 002 is good, 003 even better



rmouse_ini;
cd([WP.rootPath '\beta3_wtko']);