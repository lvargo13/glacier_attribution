clear

% get climate data ready for input into PDD model
% lauren vargo,  july 2018

% run this from directory X, also in directory X are subdir for each GCM
% idea: run 1x for each glacier (will loop all clim scen and GCMs and CESM)

% details: reads in t, p, swr (swr isn't adjusted, only set to same
% repeating cycle). since swr is from vcsn time, use that time for the day
% input into PDD model calculations for SWR (so date is set at same
% repeating cycle as swr)

% places where other variables can be adjusted: 1) get_gcm_base: the past
% period calculated from historical GCM runs, 2) calc_gcm_anom: the number of year period
% of repeating times 

%--- User input -------- 
glac = 'brewster/';  % name of glacier, e.g. 'rolleston/'
lat_glac = -44.0728; 
lon_glac = 169.436 ;

GCM = 1;  % 1 to run, 0 to not
CESM = 1; 
basey = [1961 1990]; 


%--- input and output -------- 
gcm_path = '/Volumes/arc_03/vargola/glacier_attribution/climate_data/GCM_output/' ;
cesm_path = '/Volumes/arc_03/vargola/glacier_attribution/climate_data/CESM_LE/' ;
out_gcm = ['/Volumes/arc_03/vargola/glacier_attribution/climate_data/adjusted_GCM/' glac]; % where to save
out_cesm = ['/Volumes/arc_03/vargola/glacier_attribution/climate_data/adjusted_CESM/' glac];
if exist(out_gcm,'dir') ~= 7   % if directory does not exist, create
  mkdir(out_gcm)
end
if exist(out_cesm,'dir') ~= 7   % if directory does not exist, create
  mkdir(out_cesm)
end

% --- VCSN 
tmx = ncread('tmax_N2_1980010100_2017073100_south-island_p05_daily.nc','tmax') ;
tmn = ncread('tmin_N2_1980010100_2017073100_south-island_p05_daily.nc','tmin') ;
tmin_t=squeeze(ncread('tmin_N2_1980010100_2017073100_south-island_p05_daily.nc','time'));
p = ncread('rain_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc','rain') ; 
rain_t=squeeze(ncread('rain_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc','time'));
swr = ncread('srad_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc','srad'); 
swr_t=squeeze(ncread('srad_vclim_clidb_1972010100_2018071400_south-island_p05_daily.nc','time'));
lat_vcsn = ncread('tmin_N2_1980010100_2017073100_south-island_p05_daily.nc','latitude') ;
lon_vcsn = ncread('tmin_N2_1980010100_2017073100_south-island_p05_daily.nc','longitude') ;
elev = ncread('tmax_N2_1980010100_2017073100_south-island_p05_daily.nc','elevation') ;

[~, vlat] = min(abs(lat_glac-lat_vcsn));  % find index of glacier of interest
[~, vlon] = min(abs(lon_glac-lon_vcsn));

% time - making sure all start and end at the same times
start_date=datenum(1980,1,1);  % start year
end_date=datenum(2016,12,31);  % end
rain_dn=datenum('1959-12-31 09:00:00')+double(rain_t);  % rain_t is just rain time series
temp_dn=datenum('1959-12-31 09:00:00')+double(tmin_t);
swr_dn=datenum('1959-12-31 09:00:00')+double(swr_t);
[~, rain0] = min(abs(start_date-rain_dn)) ;
[~, rain1] = min(abs(end_date-rain_dn)) ;
[~, temp0] = min(abs(start_date-temp_dn)) ;
[~, temp1] = min(abs(end_date-temp_dn)) ;
[~, swr0] = min(abs(start_date-swr_dn)) ;
[~, swr1] = min(abs(end_date-swr_dn)) ;

tmx = squeeze(tmx(vlon,vlat,temp0:temp1)) ;  
tmn = squeeze(tmn(vlon,vlat,temp0:temp1)) ;  
p = squeeze(p(vlon,vlat,rain0:rain1)) ;   
s = squeeze(swr(vlon,vlat,swr0:swr1)) ;
t = ((tmx-273.15) + (tmn-273.15)) ./2 ;  % get mean temp
elev=squeeze(elev(vlon,vlat));
climd_vcsn = [t,p,s]; 

% --- GCMs
if GCM == 1 
gcm_dir = dir(gcm_path);
gcm_dir = gcm_dir(4:end);  % first 3 are non-dir

clim_scen = {'rcp85'; 'historicalNat'}; 
clim_out = {'present'; 'past'}; 
yr_want = [2006 2026; 1901 2005] ; 

for i = 1:length(gcm_dir)

data_path = [gcm_path,gcm_dir(i).name,'/']; 

% read in lat lon
files_meta = strcat(data_path,'*.nc');
f_meta = dir(files_meta) ;
lat = ncread(strcat(data_path,f_meta(1).name), 'lat') ;  % 1 doesn't matter, all have lat lon
lon = ncread(strcat(data_path,f_meta(1).name), 'lon') ;

% find index of GCM for glacier of interest 
[~, glat] = min(abs(lat_glac-lat));
[~, glon] = min(abs(lon_glac-lon));

% get gcm historic (always same time period)
clim_base = get_gcm_base(data_path,glon,glat,'historical_',basey(1),basey(2));  

for ii = 1:length(clim_scen)  % number of climate scenarios
    cs = clim_scen(ii);
    co = clim_out(ii);
    disp(['generating data for ', gcm_dir(i).name, ' ', co{1}])

    % get gcm forcing
    wantyr = yr_want(ii,:);
    clim_force = get_gcm_force(data_path,glon,glat,cs{1},wantyr); 

    % --- calculate GCM anomalies
    [clim_adj, climd_adj, clim_vcsn, anom] = calc_gcm_anom(t,p,s,clim_force,clim_base); 
     
    sv = strcat(out_gcm,gcm_dir(i).name,'_',co{1},'.mat'); 
    save(sv,'clim_vcsn','clim_base', 'clim_force','climd_vcsn','clim_adj','climd_adj','wantyr','elev','anom');

end
end
end

% --- CESM
% get file path
if CESM == 1
cesm_ts = dir([cesm_path,'TS/','*.nc']);
cesm_runs = 34;   %right now using 34 CESM, change if use more/less of LE

% lat long for glacier of interest (same for all of ensemble)
lat = ncread([cesm_path,'TS/',cesm_ts(1).name], 'lat') ;  % 1 doesn't matter
lon = ncread([cesm_path,'TS/',cesm_ts(1).name], 'lon') ; 
[~, clat] = min(abs(lat_glac-lat));  % find index of GCM for glacier of interest 
[~, clon] = min(abs(lon_glac-lon));

% get past data
clim_base = get_cesm_base(cesm_path,1,clon,clat,basey(1),basey(2));  % get base, 1 refers to 001 LE run for anom.
clim_force = get_cesm_forceNAT([cesm_path 'PIControl/'],clon,clat);
[clim_adj, climd_adj, clim_vcsn, anom] = calc_gcm_anom(t,p,s,clim_force,clim_base); 
sv = strcat(out_cesm,'cesm_NAT.mat'); 
save(sv,'clim_vcsn','clim_base', 'clim_force','climd_vcsn','clim_adj','climd_adj','elev','anom');

% define climate scenario
clim_out = {'present'}; 
yr_want = [2006 2026] ; 

for ind = 1:cesm_runs
    
% get gcm base for anomalies
clim_base = get_cesm_base(cesm_path,ind,clon,clat,basey(1),basey(2));   

for ii = 1:length(clim_out)  % # climate scenarios
    
    co = clim_out(ii);
    disp(['generating data for cesm run ', num2str(ind), ' ', co{1}])

    % get gcm forcing
    wantyr = yr_want(ii,:);
    clim_force = get_cesm_force(cesm_path,ind,clon,clat,wantyr); 

    % --- calculate GCM anomalies
    [clim_adj, climd_adj, clim_vcsn, anom] = calc_gcm_anom(t,p,s,clim_force,clim_base);  
    
    sv = strcat(out_cesm,'cesm_le_',num2str(ind),'_',co{1},'.mat'); 
    save(sv,'clim_vcsn','clim_base', 'clim_force','climd_vcsn','clim_adj','climd_adj','wantyr','elev','anom');
end
end
end
