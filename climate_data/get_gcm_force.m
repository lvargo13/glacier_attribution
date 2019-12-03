function [ clim_force ] = get_gcm_force( data_path,lon_ind,lat_ind,clim_scen,yr_in )
% Get GCM forcing data
%   read in netcdf files for climate scenario

clim_vars = {'tas_', 'pr_'} ;  
for i = 1:2
    climv = clim_vars(i);
    climv = climv{1};
    files_str = strcat(data_path,climv,'*',clim_scen,'_','*'); 
    files = dir(files_str) ;
    ncfin = [] ;
    for nc = 1:size(files,1)
        a_all = ncread(strcat(data_path,files(nc).name),climv(1:end-1)) ; 
        ncout = squeeze(a_all(lon_ind,lat_ind,:));
        ncfin = vertcat(ncfin,ncout);
    end
    clim_mat(:,i) = ncfin ; 
end
startyr = str2num(files(1).name(end-15:end-12));  % get start year of record
start_ind = (yr_in(1) - startyr) * 12 +1 ; 
climscen_length = (yr_in(2) - yr_in(1) +1) *12; % # of months in record
clim_force = clim_mat(start_ind:start_ind+(climscen_length-1),:);
clim_force(:,1) = clim_force(:,1) - 273.15;  % convert t
clim_force(:,2) = clim_force(:,2) .* (86400*30) ;  % convert precip
end

