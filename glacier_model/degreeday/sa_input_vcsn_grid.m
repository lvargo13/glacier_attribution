function [edate,erain,tmin,tmax,swr,b_elev,lon_vcsn,lat_vcsn]=sa_input_vcsn_grid(run_years,lat_glac,lon_glac)

% VCSN netcdf lat/long to read in from gridded data
% the output is influenced by the run years

vcsn_dir='/Volumes/arc_03/vargola/glacier_attribution/climate_data/';
lat_vcsn = ncread([vcsn_dir,'tmin_N2_1980010100_2017073100_south-island_p05_daily.nc'],'latitude') ;
lon_vcsn = ncread([vcsn_dir,'tmin_N2_1980010100_2017073100_south-island_p05_daily.nc'],'longitude') ;
[~, lat_ind] = min(abs(lat_glac-lat_vcsn));
[~, lon_ind] = min(abs(lon_glac-lon_vcsn));

tmax=ncread([vcsn_dir,'tmax_N2_1980010100_2017073100_south-island_p05_daily.nc'],'tmax');
tmax=squeeze(tmax(lon_ind,lat_ind,1:end)); %tmin>tmax after that
%tmax_t=squeeze(ncread([vcsn_dir,'tmax_N2_1980010100_2017073100_south-island_p05_daily.nc'],'time'));

tmin=ncread([vcsn_dir,'tmin_N2_1980010100_2017073100_south-island_p05_daily.nc'],'tmin');
tmin=squeeze(tmin(lon_ind,lat_ind,1:end));
tmin_t=squeeze(ncread([vcsn_dir,'tmin_N2_1980010100_2017073100_south-island_p05_daily.nc'],'time'));

elev = ncread([vcsn_dir,'tmax_N2_1980010100_2017073100_south-island_p05_daily.nc'],'elevation') ;
b_elev=squeeze(elev(lon_ind,lat_ind));

% this gets the 1 grid box the glacier center is in, AND surrounding for interpolation (3x3xtime)
rain=ncread([vcsn_dir,'rain_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc'],'rain');
rain=rain(lon_ind-1:lon_ind+1,lat_ind-1:lat_ind+1,:); 
rain_t=squeeze(ncread([vcsn_dir,'rain_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc'],'time'));
lat_vcsn = lat_vcsn(lat_ind-1:lat_ind+1); 
lon_vcsn = lon_vcsn(lon_ind-1:lon_ind+1); 

% swr
swr=ncread([vcsn_dir,'srad_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc'],'srad');
swr=swr(lon_ind-1:lon_ind+1,lat_ind-1:lat_ind+1,:); 
swr_t=squeeze(ncread([vcsn_dir,'srad_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc'],'time'));

% time
start_date=datenum(run_years(1),4,1);
end_date=datenum(run_years(end)+1,3,31);

% generate daily timestep
edate=start_date:end_date;
rain_dn=datenum('1959-12-31 00:00:00')+double(rain_t);  % rain_t is just rain time series
temp_dn=datenum('1959-12-31 00:00:00')+double(tmin_t);
swr_dn=datenum('1959-12-31 00:00:00')+double(swr_t);
[~, rain0] = min(abs(start_date-rain_dn)) ;
[~, rain1] = min(abs(end_date-rain_dn)) ;
[~, temp0] = min(abs(start_date-temp_dn)) ;
[~, temp1] = min(abs(end_date-temp_dn)) ;
[~, swr0] = min(abs(start_date-swr_dn)) ;
[~, swr1] = min(abs(end_date-swr_dn)) ;
    
erain=rain(:,:,rain0:rain1);
tmax = tmax(temp0:temp1); 
tmin = tmin(temp0:temp1); 
swr = swr(:,:,swr0:swr1); 

return