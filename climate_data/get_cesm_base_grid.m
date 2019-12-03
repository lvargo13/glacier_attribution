function [ clim_b ] = get_cesm_base_grid( cesm_path,ind,lone,lonw,lats,latn,styr,edyr )
% Get CESM historic data (t and p) to compare forcing with
%   Read in data, get the grid box that includes the center of the glacier
% hardcoded: right now f_in depends on 20thcent runs being listed 1:34

hist_start_want = styr ; % start is 1850 first run, 1920 rest
hist_end_want = edyr; % end of all runs
hist_length = (hist_end_want-hist_start_want+1)*12 ; 

clim_vars = {'TS', 'PRECC', 'PRECL'} ; 
clim_base = zeros(hist_length,3); 
for i = 1:3
    climv = clim_vars(i);
    climv = climv{1};  % get string from structure
    files = dir([cesm_path,climv,'/','*.nc']); 
    f_in = ncread([cesm_path,climv,'/',files(ind).name],climv) ; % ind should only go to 34, only PI, not RCP
    ta = f_in(lone:lonw,lats:latn,end-(hist_length-1):end);
    clim_base(:,i) =squeeze(mean2d(ta)) ;  
end

clim_b = [clim_base(:,1)-273.15, clim_base(:,2)+clim_base(:,3)] ;  % convert t, add p
clim_b(:,2) = clim_b(:,2) .* (86400*30) ; % convert precip 1 kg/m2/s = 86400 mm/day 
end