global DS AP WP 


% learn data root path from rmouse ini routine
rmouse_ini;

% add path because various commonly used ANPAR files resides there
addpath([WP.rootPath '\beta3_wtko']);


rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
AP_beta3_wtko;
AP_job1;
a002_tmp;
rmouse;


return
return
return
return


rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2484_02927']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job2;  
a002_r1;
load(['streams\' DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all

return
return
return
return


rmouse_ini;
cd([WP.rootPath '\beta3_wtko\ko2483_03303']);
AP=[]; DS=[];
AP_beta3_wtko;
AP_job1;  
r1;
load([DS.abfFn '_diff_d_fpspike_res.mat']);
rmouse('af_tsl_ext',evt.tsl{1});
% clear; close all


return
return
return
return








try
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko\ko3363_02313']);
  AP=[]; DS=[];
  AP_beta3_wtko;
  AP_job1;
  r1;
  rmouse;
  clear; close all; pack
catch
  lasterr
  sendmail('harald.hentschke@uni-tuebingen.de','rmouse problem',lasterr);
end


% -----------------------------------------------------
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
AP_beta3_wtko;
AP_job1;
a001_r1;
rmouse('callJob','seg_thetaTrigStream');
clear; close all; pack



rmouse_ini;
cd([WP.rootPath '\beta3_wtko']);

% rmouse('callJob','rmouse_p_tcthpeaks_greenorange(p_etsl,boe,eoe)');
% rmouse('callJob','rmouse_genDiffTrace');
% rmouse('af_tsl_ext',evt.tsl{1});