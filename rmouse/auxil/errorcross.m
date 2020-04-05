function ebh=errorCross(x,y,avType,varargin)
% ** function ebh=errorCross(x,y,avType,varargin)
% plots error crosses using function errorbar (** new version of errorcross
% written April 2017, not fully documented yet)
%                         >>> INPUT VARIABLES >>>
% NAME      TYPE/DEFAULT          DESCRIPTION
% 
% x        column array        x values
% y        column array        y values
% avType   char, 'mn'          type of average, 'mn' (mean) or 'md'
%                              (median) ** NOTE: if 'md' error bars will
%                              automatically be set to 'qrt'
%        -- DISABLED --
% varargin                     any input into errorbar beyond x,y,yneg,ypos,xneg,xpos

avType = 'mn';

switch avType
  case 'md'
    prc=prctile(x,[25 50 75],1);
    avX=prc(2,:);
    errXNeg=prc(2,:)-prc(1,:);
    errXPos=prc(3,:)-prc(2,:);    
    
    prc=prctile(y,[25 50 75],1);
    avY=prc(2,:);
    errYNeg=prc(2,:)-prc(1,:);
    errYPos=prc(3,:)-prc(2,:);    

  case 'mn'
    avX=nanmean(x,1);
    avY=nanmean(y,1);
    errXNeg=nanstd(x,0,1);
    errXPos=errXNeg;
    errYNeg=nanstd(y,0,1);
    errYPos=errYNeg;

  otherwise
    error('illegal value for input arg avType')
end
ebh=errorbar(avX,avY,errYNeg,errYPos,errXNeg,errXPos,varargin{:});
