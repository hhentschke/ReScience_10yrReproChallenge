% r=orderfields(r);
subplot(2,2,1)
plot(r(4).rawGaPE{11,11},r(4).rawThPE{11,11},'.');
subplot(2,2,2)
plot(r(4).rawGaPE{11,11},r(4).rawPPeak{11,11},'.');

subplot(2,2,3)
plot(r(4).rawThPE{11,11},r(4).thgaeCCPeak(:,11),'.');
subplot(2,2,4)
plot(r(4).rawGaPE{11,11},r(4).thgaeCCPeak(:,11),'.');


plot(r(1).F,r(4).rawPMn{11,11})
set(gca,'xlim',[20 120])
set(gca,'yscale','log')