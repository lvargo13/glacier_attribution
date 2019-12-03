function [ clim_f ] = get_cesm_force( cesm_path,ind,clon,clat,wantyr )
% Get cesm forcing data
%   in: cesm_runs- number of CESM runs in LE (34 now, could be 40)
%   getting 20th century vc rcp85 is a bit ugly

clim_mat = []; 

clim_vars = {'TS', 'PRECC', 'PRECL'} ; 
for i = 1:3
    climv = clim_vars(i);
    climv = climv{1};
    ind_str = sprintf('%03d', ind); % turn ind into str with padded zeros
    files = dir([cesm_path,climv,'/*',ind_str,'.cam*.nc']); % all files of CESM LE ind 
    temp_array = []; 
    for ii = 1:length(files)
      f_in = ncread([cesm_path,climv,'/',files(ii).name],climv) ;
      temp_array = [temp_array; squeeze(f_in(clon,clat,:))]; % full time (1920 - 2100)
      f_in = []; 
    end
    clim_mat = [clim_mat temp_array]; 
end

% get timing
start_ind = (wantyr(1) - str2num(files(1).name(end-15:end-12))) * 12 +1  ;
clim_length = (wantyr(2) - wantyr(1) +1) *12; % # of months in record
clim_force = clim_mat(start_ind:start_ind+(clim_length-1),:);

clim_f = [clim_force(:,1)-273.15, clim_force(:,2)+clim_force(:,3)] ;  % convert t, add p
clim_f(:,2) = clim_f(:,2) .* (86400*30) ; % convert precip 1 kg/m2/s = 86400 mm/day 
end