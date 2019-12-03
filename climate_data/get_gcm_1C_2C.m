function [ c1,c2,c1_decind,c2_decind,diff_t ] = get_gcm_1C_2C( data_path,lon_ind,lat_ind,clim_scen,clim_base )
% Get GCM forcing data,  read in netcdf files for climate scenario
% INPUTS: clim_base: climate calculated in get_gcm_base
% does not include output time (to know what years are included

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
clim_force(:,1) = clim_mat(:,1) - 273.15;  % convert t
clim_force(:,2) = clim_mat(:,2) .* (86400*30) ;  % convert precip

% up to here, have all time, want to make sure 2006 - 2100 
if length(clim_force) > (95*12) % months in 95 years
   clim_force = clim_force(1:1140,:);  
end
meanbase = mean(clim_base(:,1)); 
vec1 = 1:120:length(clim_force); 
vec2 = 120:120:length(clim_force);
decmean_force = length(vec2); 
for i = 1:length(vec2)  
    decmean_force(i) = mean(clim_force(vec1(i):vec2(i),1));  % get annual mean
end

diff_t = decmean_force - meanbase ; % difference
c1_decind = find(diff_t > 1.3 &  diff_t< 1.7);
c2_decind = find(diff_t > 1.8 & diff_t < 2.2);
c1 = clim_force(vec1(min(c1_decind)):vec2(max(c1_decind)),:); 
c2 = clim_force(vec1(min(c2_decind)):vec2(max(c2_decind)),:);

end