% load file distrib_npa.mat, compare variances of distributions of peak amplitude 
load c:\projects\hippocampus\rmouse\paper_atropine\rawFig\distrib_npa

% subtract medians 
mn_imm_ctrl=mean(imm_ctrl);
mn_imm_atr=mean(imm_atr);
mn_exp_ctrl=mean(exp_ctrl);
mn_exp_atr=mean(exp_atr);



[h,p]=ansaribradley(imm_ctrl-mn_imm_ctrl,imm_atr-mn_imm_atr)
[h,p]=ansaribradley(exp_ctrl-mn_exp_ctrl,exp_atr-mn_exp_atr)


% [h,p]=vartest2(imm_ctrl,exp_ctrl)
% [h,p]=vartest2(imm_ctrl,imm_atr)


% are they really different? make a test with simulated data
if 0
  d1=mean(imm_ctrl)+randn(size(imm_ctrl))*std(imm_ctrl);
  d2=mean(imm_atr)+randn(size(imm_atr))*std(imm_atr);
else
  d1=-0.3117+randn(4497,1)*0.2041;
  d2=-0.3331+randn(3592,1)*0.2357;
end
[h1,x]=hist(d1,100);
[h2,x]=hist(d2,x);

plot(x,h1,'b',x,h2,'r');
[h,p]=ansaribradley(d1,d2)

% repeat test with 1/40th of the data: now, with only ~100 pts of each
% distribution being compared the difference is not significant anymore
[h,p]=ansaribradley(d1(1:40:end),d2(1:40:end))
