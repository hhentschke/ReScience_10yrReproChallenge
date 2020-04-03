
nnA=diag(r(3).thCCPeak,1);
nnT=diag(r(3).thCCPeakT,1);
nnA=cat(1,nnA{:});
nnT=cat(1,nnT{:});

for i=1:size(nnA,2)
  plot(nnA(:,i),'-o');
  set(gca,'ylim',[0 1.1]);
  drawnow
  pause(.2)
end

for i=1:size(nnA,2)
  plot(nnT(:,i),nnA(:,i),'-o');
  set(gca,'xlim',[-10 30],'ylim',[0 1.1]);
  drawnow
  pause(.2)
end

break
