% rmouse('AP_beta3_wtko_iStateTheta');
global ANPAR DSET AP DS WP R D ND

% collect DS and AP
rmouse_ini;
ci=0; i3=0;

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
AP=[]; DS=[];
a001_iState;
AP_beta3_wtko_iStateTheta;
rmouse_APcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

% quick & dirty movie
% load(AP.resFn);
% qdccmovie(r,{'theta'},'beh',{'expl'},'chanComb','all');
% return

% pars={'thCCPeak','thCCPeakT','gaeCCPeak','gaeCCPeakT'};
% pars={'thCCPeakT'};
% pars={'thCCPeak'};
pars={'thCCPeak','thCCPeakT'};

% concatenate results
% rcat({'immobile','exploring'},pars);
rcat({'exploring'},pars);
% rcat({'immobile'},pars);

% here's the place to reduce # of channels
R.AP.rawChAnNm=R.AP.rawChAnNm(4:12); % wt0001
% R.AP.rawChAnNm=R.AP.rawChAnNm(5:10); % ko3363
% R.AP.rawChAnNm=R.AP.rawChAnNm(1:9); % ko2483


% R contains all fields listed in variable pars above; now we want the PC calcs 
% be done with one parameter only 
pars={'thCCPeakT'};
% rcatfish(pars,'chanComb','neigh');
rcatfish(pars,'chanComb','princ');

% return
% normalize:
% - for each par, pick the channel with the highest standard deviation and
% divide all by it, unless we're dealing with only one par
varPerPar=size(D,2)/length(pars);
if length(pars)>1
  for g=1:length(pars)
    [s,ix]=max(std(D(:,(g-1)*varPerPar+1:g*varPerPar),0,1));
    D(:,(g-1)*varPerPar+1:g*varPerPar)=D(:,(g-1)*varPerPar+1:g*varPerPar)/s;
  end
end
% return

% {'line','scatter','density','vector'}

% PCexplore('plotType',{'scatter','density'},'t',R.tb(:,1),'tag',R.tb(:,2),'nPC',5,'normalize',0);
% return

oD=D;
D=D(1:2000,:);
PCexplore('plotType',{'line','scatter','density'},'t',R.tb(1:2000,1),'tag',R.tb(1:2000,2),'nPC',5,'normalize',0);

% return
% invoke single AP here
helper01;



%ccmovie(R,'rv',{'thCCPeak','thCCPeakT'},'mode','all');


% --------- crap ------------ crap -----------
