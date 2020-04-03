% produces bar plots of the memory taken up by all fields of results variable r
rix=4;
r=orderfields(r);
s=fieldnames(r);
for i=1:length(s)
  eval([s{i} '=r(rix).' s{i} ';']);
  eval(['fs(i)=whos(''' s{i} ''');']);
  eval(['clear ' s{i}]);
end
bytes=cat(1,fs.bytes);
% sort according to bytes occupied and put out 25 largest
n=25;
[y,ix]=sort(bytes,1);
tickl=strvcat(fs(ix(end-n+1:end)).name);
mb=2^20;

% overview
subplot(1,2,1);
barh(bytes/mb);
xlabel('Size of field (Mb)');
title(['all (sum: ' int2str(round(sum(bytes)/mb)) ' Mb)']);

subplot(1,2,2);
barh(bytes(ix(end-n+1:end))/mb);
set(gca,'ytick',1:n,'ydir','reverse','yticklabel',tickl);
nicexax;
xlabel('Size of field (MB)');
title(['largest 20 (sum: ' int2str(round(sum(bytes(ix(end-n+1:end)))/mb)) ' Mb)']);

