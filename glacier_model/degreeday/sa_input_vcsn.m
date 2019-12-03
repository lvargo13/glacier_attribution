function [edate,erain,etemp]=sa_input_vcsn(run_years,hourly)

% VCSN netcdf lat/long 169.425,44.075, i=62,j=75
%ncks -d latitude,-44.05 -d longitude,169.43 tmax_N2_1972010200_2015080100_south-island_p05_daily.nc brewster_tmax_N2_1972010200_2015080100_south-island_p05_daily.nc
%ncks -d latitude,-44.05 -d longitude,169.43 tmin_N2_1972010200_2015080100_south-island_p05_daily.nc brewster_tmin_N2_1972010200_2015080100_south-island_p05_daily.nc
%ncks -d latitude,-44.05 -d longitude,169.43 rain_vclim_clidb_1972010100_current_south-island_p05_daily.nc brewster_rain_vclim_clidb_1972010100_current_south-island_p05_daily.nc
% short time(time) ;
%                 time:standard_name = "time" ;
%                 time:long_name = "time (end of reporting interval)" ;
%                 time:units = "days since 1959-12-31 09:00:00.0 +12:00" ;
%                 time:bounds = "time_bounds" ;
%                 time:calendar = "gregorian" ;
%                 time:axis = "T" ;

if nargin==1
    hourly=0;
end

vcsn_dir='/Volumes/arc_03/vargola/glacier_attribution/climate_data/';
tmax=squeeze(ncread([vcsn_dir,'brewster_tmax_N2_1972010200_2015080100_south-island_p05_daily.nc'],'tmax'));
tmax=tmax(1:15615); %tmin>tmax after that
tmax_t=squeeze(ncread([vcsn_dir,'brewster_tmax_N2_1972010200_2015080100_south-island_p05_daily.nc'],'time'));
tmax_t=tmax_t(1:15615); %tmin>tmax after that
tmin=squeeze(ncread([vcsn_dir,'brewster_tmin_N2_1972010200_2015080100_south-island_p05_daily.nc'],'tmin'));
tmin=tmin(1:15615);
tmin_t=squeeze(ncread([vcsn_dir,'brewster_tmin_N2_1972010200_2015080100_south-island_p05_daily.nc'],'time'));
tmin_t=tmin_t(1:15615);
rain=squeeze(ncread([vcsn_dir,'brewster_rain_vclim_clidb_1972010100_current_south-island_p05_daily.nc'],'rain'));
rain=rain(1:15859); %Nan at end
rain_t=squeeze(ncread([vcsn_dir,'brewster_rain_vclim_clidb_1972010100_current_south-island_p05_daily.nc'],'time'));
rain_t=rain_t(1:15859); %Nan at end



% [tmax_t-tmin_t] shows that times are indentical
% find shortest time dimension
vcsn_time=max([min(tmax_t) min(tmin_t) min(rain_t)]):min([max(tmax_t) max(tmin_t) max(rain_t)]);
vcsn_dn=datenum('1959-12-31 09:00:00')+double(vcsn_time);

start_date=datenum(run_years(1),4,1);
end_date=datenum(run_years(end)+1,3,31);

if hourly
    rain0=findclose(vcsn_time(1),rain_t,1);
    rain1=findclose(vcsn_time(end),rain_t,1);
    temp0=findclose(vcsn_time(1),tmax_t,1);
    temp1=findclose(vcsn_time(end),tmax_t,1);
    % generate hourly timestep
    edate=start_date:1/24:end_date;
    % quick interpolation using nearest neighbour
    % think about correct offset for 9am data (none at the moment)
    erain=interp1(vcsn_dn,rain(rain0:rain1),edate,'nearest')/24;
    
    % interleave tmin and tmax and interpolate
    tminmax=[tmin(temp0:temp1) tmax(temp0:temp1)];
    tminmax_t=[vcsn_dn'-3/24,vcsn_dn'+9/24];
    %reshape into vectors
    etemp=interp1(reshape(tminmax_t,1,[]),reshape(tminmax,1,[]),edate,'spline')-273.15;
else
    % generate daily timestep
    edate=start_date:end_date;
    rain_dn=datenum('1959-12-31 09:00:00')+double(rain_t);  % rain_t is just rain time series
    temp_dn=datenum('1959-12-31 09:00:00')+double(tmin_t);
    [~, rain0] = min(abs(start_date-rain_dn)) ;
    [~, rain1] = min(abs(end_date-rain_dn)) ;
    [~, temp0] = min(abs(start_date-temp_dn)) ;
    [~, temp1] = min(abs(end_date-temp_dn)) ;
    %rain0=findclose(start_date,rain_dn,1);
    %rain1=findclose(end_date,rain_dn,1);
    %temp0=findclose(start_date,temp_dn,1);
    %temp1=findclose(end_date,temp_dn,1);
    
    erain=rain(rain0:rain1);
    tminmax=[tmin(temp0:temp1) tmax(temp0:temp1)];
    etemp=mean(tminmax,2)-273.15;
end
return
