% extracts selected parameters from results structure R. R is the combination of results
% from several experiments, won by running combine3.m
% !! **** the number and order of parameters listed here must correspond to
%         that listed in combine3.m!! ***
% first index = behavior (1=immobile, 2=exploring)
% second index = parameter 
% 1 (rawPMnPeak) = peak power (theta)
% 2 (rawPMnPeakT) = peak frequency (theta)
% 3 (rawThPEMn) = power in [4 12] Hz
% 4 (rawGaPEMn) = power in [30 90] Hz
% 5 (thgaeCCPeakMn) = peak crosscorrelation theta vs gamma envelope
% example:
% rr=R.d{2,3};  exploring, theta power
pas={'peak power (theta,non-segmental)','peak frequency (theta,non-segmental)',...
     'peak power (theta,segmental)','peak frequency (theta,segmental)',...    
     'power in [4 12] Hz','power in [30 90] Hz',...
     'peak amp CC theta vs gamma envelope','peak freq CC theta vs gamma envelope'};

bs={'immobile','exploring'};

% if nonzero, wash values will be displayed too
wd=1;

for bi=1:2
  disp('*****************************************************');
  disp(['******************** ' bs{bi} ' **************'])
  disp('*****************************************************');  
  for pai=1:length(pas)
    disp(pas{pai})
    rr=R.d{bi,pai};
    indv=R.indv{bi,pai};
    % ------- don't change below this line --------
    % princ channel
    ix=rr(:,1)==0;
    % new column: normalized values
    rr(:,6)=rr(:,3)./rr(:,2);
    if wd
      % new column: normalized values
      rr(:,7)=rr(:,4)./rr(:,2);
      % order of display: conc, control, applic, wash, applic_normalized, wash_normalized, name
      disp([num2str(rr(ix,[5 2 3 4 6 7]),'%2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f')  repmat('  ',length(find(ix)),1)  char(indv(ix,:))]);
    else
      % order of display: conc, control, applic, applic_normalized, name
      disp([num2str(rr(ix,[5 2 3 6]),'%2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f')  repmat('  ',length(find(ix)),1)  char(indv(ix,:))]);
    end
  end
end