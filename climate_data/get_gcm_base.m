function [ clim_base ] = get_gcm_base(data_path,lon_ind,lat_ind,base_run,styr,edyr)
% Get GCM historic data (t and p) to compare forcing with
%   Read in data, get the grid box that includes the center of the glacier
% base_run: 'historical_' or 'historicalNat_'

% read in historical base files, generate .mat files 
clim_vars = {'tas_', 'pr_'} ;  
for i = 1:2
    climv = clim_vars(i);
    climv = climv{1};
    files_str = strcat(data_path,climv,'*',base_run,'*.nc'); 
    files = dir(files_str) ;
    ncfin = [] ;
    for nc = 1:size(files,1)
        a_all = ncread(strcat(data_path,files(nc).name),climv(1:end-1)) ; %end-1 gets rid of underscore 
        ncout = squeeze(a_all(lon_ind,lat_ind,:));
        ncfin = vertcat(ncfin,ncout);
    end
    clim_mat_hist(:,i) = ncfin ;  % this is full record 
end
hist_start_want = styr ; % base could be difference
hist_end_want = edyr; 
hist_length = (hist_end_want-hist_start_want+1)*12 ; 
hist_start = str2num(files(1).name(end-15:end-12));  % get start year of record
start_ind_hist = (hist_start_want - hist_start) * 12 +1 ; 
clim_base = clim_mat_hist(start_ind_hist:start_ind_hist+(hist_length-1),:);
clim_base(:,1) = clim_base(:,1) - 273.15;  % convert t
clim_base(:,2) = clim_base(:,2) .* (86400*30) ; % convert precip 1 kg/m2/s = 86400 mm/day 
end

