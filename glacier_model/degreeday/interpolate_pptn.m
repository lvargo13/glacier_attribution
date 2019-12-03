function [ interp_precip ] = interpolate_pptn( p,dem_file,lon_vcsn,lat_vcsn )
%function [ p,lonn,latn,lon_want,lat_want ] = interpolate_pptn( p,dem_file,lon_vcsn,lat_vcsn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lat_dem = ncread(dem_file,'northing'); 
lon_dem = ncread(dem_file,'easting'); 
want_res = 10 ; % output resolution we want (in m)

lat_v = lat_vcsn' ;
lon_v = lon_vcsn' ;
lat_v = round(lat_v,7); 
lon_v = round(lon_v,7); 

[lonn,latn] = geod2nztm(lat_v,lon_v); 
lon_want = lon_dem(1,1):want_res:lon_dem(end,end); 
lat_want = lat_dem(1):want_res:lat_dem(end,end);

p = p';  % now is in map view correctly
interp_precip = interp2(lonn,latn',p,lon_want,lat_want','linear');
interp_precip = flipud(interp_precip); 

end

