function [x y] = geod2nztm(lat,lon)
% function [x y] = geod2nztm(lat, lon)
% returns NZTM eastings (x) and northings (y) given latitude and 
% longitude in WGS84.

proj=load('proj_nztm.mat');

lat=lat(:);
lon=lon(:);
[x y]= projfwd(proj,lat,lon);
