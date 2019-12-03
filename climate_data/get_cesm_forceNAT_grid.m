function [ clim_f ] = get_cesm_forceNAT_grid( cesm_path,lone,lonw,lats,latn)
% Get cesm forcing data for 1800 year natural run

clim_mat = zeros(1800*12,3); 
vecs = 1:1200:(1800*12); 
vece = 1200:1200:(1800*12);
clim_vars = {'TS', 'PRECC', 'PRECL'} ; 
for i = 1:3
    climv = clim_vars(i);
    climv = climv{1};
    files = dir([cesm_path,climv,'/*.nc']); % all files of CESM LE  
    for ii = 1:length(files)  % should always be 18
      f_in = ncread([cesm_path,climv,'/',files(ii).name],climv) ;
      ta = f_in(lone:lonw,lats:latn,1:1200);
      clim_mat(vecs(ii):vece(ii),i) = squeeze(mean2d(ta));
      f_in = []; % release mem
    end
end

clim_f = [clim_mat(:,1)-273.15, clim_mat(:,2)+clim_mat(:,3)] ;  % convert t, add p
clim_f(:,2) = clim_f(:,2) .* (86400*30) ; % convert precip 1 kg/m2/s = 86400 mm/day 
end