function fComodSummary
% summary of frequency comodulograms: integral of pixels corresponding to 
% freq combinations listed in combine_fcomod plotted as bars (ctrl and drug
% in one figure;  principal channel and exploring only)

global D
D=[];

behav={'immobile','exploring'};
drug={'ctrl','atrop'};

curBehav=2;

% results variable R:
% - a struct array (as many elements as behaviors)
% - the field(s) are 1D cell arrays (one cell per channel); 
%   each data set will be embedded in the correct position
% - each cell is a 3D array:
%   - freq combo down columns, 
%   - column order control | drug | recovery (optional)
%   - different experiments in slices


% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.7,'lineW',.8,'markSz',8); 
ornt='portrait';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
figdir='';
printas='-dpsc';[];
printas=[];

% --- load 
load fComod_all;

curChInd=hp.princChInd;

nFBand=size(fBand,1);

cnt=0;
[combo,ncombin]=combin(nFBand,'autoC',1);
% --- plot
for bi=curBehav
  nAnimal=size(R(bi).r{curChInd},3);
  mn=mean(R(bi).r{curChInd},3);
  stdd=std(R(bi).r{curChInd},0,3);
  % ---------------------------------------------
  % global input var for PC: observations (=animals under different
  % drug conditions) in rows, variables (the different comodulation types)
  % in columns
  D=permute(R(bi).r{curChInd},[3 1 2]);
  D=[D(:,:,1);D(:,:,2)];
  % tag: control vs. drug
  tag=[ones(nAnimal,1); ones(nAnimal,1)*2];
  % ---------------------------------------------  
  figName=[mfilename '_' behav{bi}];
  orient(ornt);
  
  % bar(mean(R(bi).r{curChInd},3));
  %   mnMat=repmat(nan,4,4);
  %   stdMat=mnMat;

  for g=1:ncombin
    %     mnMat(combo(g,1),combo(g,2))=mn(g,drugi);
    %     stdMat(combo(g,1),combo(g,2))=stdd(g,drugi);
    pupu=permute(R(bi).r{curChInd}(g,:,:),[3 2 1]);
    [h,p]=ttest(pupu(:,1),pupu(:,2));

    subplot(nFBand,nFBand,(combo(g,1)-1)*nFBand+combo(g,2))
    axis off
%    set(gca,
    bh=bar([1 2],mn(g,:),.5);
    set(bh(1),'facecolor','k');
    errorCross([[1; 2] mn(g,:)'],[[0; 0]  stdd(g,:)'],'color','k');
    box on
    set(gca,'xtick',[]);
    nicexyax;
    set(gca,'xlim',[.4 2.6]);
    lh=line(get(gca,'xlim'),[0 0],'color','k','linestyle',':');
    %    set(lh,'linestyle','-','linewidth',
    % rexy('ax',gca,'xfac',.8,'yfac',1.1);
    
    th=ultext(num2str(p,'%1.3f'));
    set(th,'fontsize',6,'color','r');
    th=title([fbLabel{combo(g,1)} ' - ' fbLabel{combo(g,2)}]);
    set(th,'fontsize',9);    

  end

  if ~isempty(printas),
    print(printas,[figdir figName]);
  end
end % for:behaviors

% PCexplore('tag',tag,'normalize',1);
% 
% disp('feddich');


