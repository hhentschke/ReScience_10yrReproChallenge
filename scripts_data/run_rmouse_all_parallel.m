% Script for running all available analyses through rmouse in parallel
% fashion. Requires the Parallel Computing Toolbox and assumes at least 
% four available workers.
p = gcp;
start = datetime(now, 'ConvertFrom', 'datenum')

job1 = parfeval(@run_wt0001_04708, 0);
job2 = parfeval(@run_wt0002_04707, 0);
job3 = parfeval(@run_wt0003_04730, 0);
disp('dispatched three jobs')

fetchOutputs(job1);
disp('finished')
stop = datetime(now, 'ConvertFrom', 'datenum')