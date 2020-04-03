try
  doall_ko;
catch
  disp(lasterr)
  sendmail('harald.hentschke@uni-tuebingen.de','rmouse problem: ko',lasterr);
  relax;
  global WP  
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end

try
  doall_wt;  
catch
  disp(lasterr)
  sendmail('harald.hentschke@uni-tuebingen.de','rmouse problem: wt',lasterr);
  relax;
  global WP
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end



return
return



try
  doall_wtAtropine_ctrlRedundant;
catch
  disp(lasterr)
  relax;
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end


try
  doall_wtAtropineIso;
catch
  disp(lasterr)
  relax;
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end

try
  doall_wtProspect;
catch
  disp(lasterr)
  relax;
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end



return

try
  dosome;
catch
  disp(lasterr)
  relax;
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko']);
end




