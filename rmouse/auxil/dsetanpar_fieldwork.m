job=2;

% generate list of file names including full path of the files to be modified
switch job
  case 1
    % !dir dset*.m /s /b > d:\flist.txt
    !dir r*.m /s /b > d:\flist.txt
    !dir a*.m /s /b >> d:\flist.txt
  case 2
    !dir *.m /s /b > d:\flist.txt
%     !dir r*.m /s /b > d:\flist.txt
%     !dir a*.m /s /b >> d:\flist.txt
end

% read that file, each file listed being assigned to a cell element
fili=textread('d:\flist.txt','%s','delimiter','\n');
disp(strvcat(fili));
bn=questdlg('*.m parameter files as listed in command window will be changed (=potentially messed up) - are you absolutely, positively sure you want to continue?','A little intimidation');
showtag=0;
if strcmpi(bn,'yes')
  for ii=1:length(fili)
    problem='you chickened out';
    % §§§ the problem with code files is that indentation will be gone...
    ftext=textread(fili{ii},'%c','delimiter','\n');
    switch job
      case 1
        % **** insert a line in a specific place within a dset/anpar:
        
        % 1. the new line(s)
        sNew=[];
        sNew{1}='% some recording hardware (notably the Neuralynx) inverts signals by default.';
        sNew{2}='% set this parameter to a nonzero value if that is the case';        
        sNew{3}='DS.rawSignalInverted=          0;';
        % 2. the first (or all) characters of the line after which the new line(s)
        % should be inserted
        s1='DS.abfFn=';

        % 1. the new line(s)
        sNew=[];
        sNew{1}='AP.thetaLo=                    [4 6];';
        sNew{2}='AP.thetaHi=                    [6 12];';
        sNew{3}='AP.beta=                       [15 30];';
        % 2. the first (or all) characters of the line after which the new line(s)
        % should be inserted
        s1='AP.theta=';

        
        % 1. the new line(s)
        sNew=[];
        sNew{1}='AP.thetaLoCFreq=               [4 6];';
        sNew{2}='AP.thetaHiCFreq=               [6 12];';
        sNew{3}='AP.betaCFreq=                  [15 30];';
        % 2. the first (or all) characters of the line after which the new line(s)
        % should be inserted
        s1='AP.thetaCFreq=';
        
        % now check whether the new stuff (the first line of it) is already
        % present 
        sNewIx=strmatch(sNew{1},ftext);
        if isempty(sNewIx)
          % identify the line after which we want to insert
          s1Ix=strmatch(s1,ftext);
          if length(s1Ix)==1
            % raster shift
            nShift=length(sNew);
            ftext(s1Ix+1+nShift:end+nShift)=ftext(s1Ix+1:end);
            ftext(s1Ix+1:s1Ix+1+nShift-1)=sNew;
            % show new contents for first file
            if ~showtag
              disp(strvcat(ftext));
              bn=questdlg('Take a look at the new contents of 1st file in list in the command window (you will not be asked for the others). Continue?','last chance');
              showtag=1;
            end
            if strcmpi(bn,'yes')
              problem='no problem';
            end
          elseif length(s1Ix)>1
            problem=['more than one line in ' fili{ii} ' start with ' s1 ' as the alignment string'];
          else
            problem=['no line in ' fili{ii} ' starting with ' s1 ' as the alignment string'];
          end
        else
          problem=['new string already present in ' fili{ii}];
        end

      case 2
        % **** replace parts of a line:
        % 1. the new part
        sNew='AP.afThresh= -1*';
        % 2. the part to be replaced
        s1='AP.afThresh=';
        
        % 1. the new part
        sNew='AP.gammaCFreq=                 [40 90];';
        % 2. the part to be replaced
        s1='AP.gammaCFreq=                 [30 90];';
        
        % 1. the new part
        sNew='rawCh=rmouse_chan;';
        % 2. the part to be replaced
        s1='rmouse_chan;';

        % now check whether the new stuff (the first line of it) is already
        % present 
        sNewIx=strmatch(sNew,ftext);
        if isempty(sNewIx)
          % identify the line which we want to modify
          s1Ix=strmatch(s1,ftext);
          if length(s1Ix)==1
            % replace
            ftext{s1Ix}=[sNew ftext{s1Ix}(length(s1)+1:end)];
            % show new contents for first file
            if ~showtag
              disp(strvcat(ftext));
              showtag=1;
              bn=questdlg('Take a look at the new contents of 1st file in list in the command window (you will not be asked for the others). Continue?','last chance');
            end
            if strcmpi(bn,'yes')
              problem='no problem';
            end
          elseif length(s1Ix)>1
            problem=['more than one line in ' fili{ii} ' start with ' s1 ' as the alignment string'];
          else
            problem=['no line in ' fili{ii} ' starting with ' s1 ' as the alignment string'];
          end
        else
          problem=['new string already present in ' fili{ii}];
        end
      otherwise
        error('unknown job number');
    end % switch:job

    if strcmpi(problem,'no problem')
      % now write back
      fid=fopen(fili{ii},'w');
      for li=1:length(ftext)
        fprintf(fid,'%s\n',ftext{li});
      end
      fclose(fid);
      disp(['modified ' fili{ii}]);
    else
      warning(problem);
      if strfind(problem,'chicken')
        break
      end
    end
  end
end

! del d:\flist.txt