function ts=findBehavTransit(etsl,bcode1,bcode2)
% ** function findBehavTransit(tsl,bcode1,bcode2)
% spits out and displays time points (in seconds) of transitions from 
% one behavior (bcode1) to another (bcode2). Needs an extended time stamp 
% list (produced by rmouse.m, job 'gen_btsl') and codes for the behaviors 
% (classically, 3=immobile and 4=exploring). This function is useful for finding
% interesting raw data segments for plots

imm=etsl(:,4)==bcode1;
ex=etsl(:,4)==bcode2;
% transitions I->E in seconds
ts=etsl(find(imm(1:end-1)&ex(2:end))+1,1)*.001;
disp(num2str(ts,'%5.1f'));