% this file does the following conversion of behavioral episode information:
% [episode #] -> [min.sec]
% IT IS MANDATORY TO LIST ALL EPISODES CONTAINED IN FILE 



% 02220001a.abf, wt2008
expl=[1:4 7:9 13 14 18:31 33:37 39 40 51 53:60 69 73:75 80:89 91 95:99 108 109]';
imm=[6 10:12 15 17 32 41:50 52 61:65 68 70:72 76:79 90 92:94 100:107]';
bad=[5 16 38 66 67]';

% 02424001.abf, wt2001
expl=[1:3 6:15 17:23 26 27 32:34 38:44 47 49:53 63:67 69:77 79:82 88 97 105:107]';
imm=[4 5 24 25 28:31 35:37 45 46 59:62 68 78 83:87 89:96 98:104 108 109]';
bad=[16 48 54:58]';

% 02510000.abf, wt2141
expl=[1:4 6:9 12:26 29:35 37 38 40 42 46 56 57 62 64:67 72 73 76 78:82 93]';
imm=[48 49 51:53 70 71 88:91 94 95 97 100:108]';
bad=[5 10 11 27 28 36 39 41 43:45 47 50 54 55 58:61 63 68 69 74 75 77 83:87 92 96 98 99]';


allev=sort([expl;imm;bad]);
% check for multiply assigned episodes
if ~isempty(intersect(expl,imm)) | ~isempty(intersect(expl,bad)) | ~isempty(intersect(imm,bad))
  error('there are multiply assigned episodes'); 
end
% check for omitted episodes
if ~isequal((1:allev(end))',allev)
  setdiff(1:allev(end),allev)
  error('displayed episodes are unassigned'); 
end

% put together & sort
btl=sortrows(cat(1,[expl repmat(7,length(expl),1)],[imm repmat(3,length(imm),1)],[bad repmat(-2,length(bad),1)]));
% now convert from ms to format min.sec, taking into account that time  correspond to ENDS of episodes 
% (unit: episode number; 1 episode = 16384 ms)
btl(:,1)=(btl(:,1)-1)*16384;
btl(:,1)=floor(btl(:,1)/60000)+(btl(:,1)/60000-floor(btl(:,1)/60000))*.6;
% display in appealing way
sprintf('%3.5f   %1.1i\n',btl')