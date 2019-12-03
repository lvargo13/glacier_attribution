function [ zenith, azimuth, S] = inso_sunpos(latitude,daynum)

% function [D, zenith, azimuth, S] = inso_sunpos(kyr,latitude,daynum)
%
% use inso to return sun angle information (azimuth and zenith) as well as 
% daily mean insolation (D, W/m^2) and instantanous insolation at the top of the
% atmosphere (S, W/m^2). kyr is thousands of years before present, and daynum is
% the day of the year (with a fractional part to indicate the hour)

So = 1365; % w/m2: solar constant

% calculate hour angle from fractional part of daynum:
% H = 0 is midday, and Ho is +/- hour angle for sunset/sunrise
hour_angle = (daynum - floor(daynum)) * 2 * pi - pi;        

[ecc,long_perh,delta,lambda] = inso(latitude,daynum);
%[ecc,long_perh,delta,lambda] = daily_insolation(0,latitude,daynum,1)
% added lambda as return argument to inso

omega = long_perh * pi/180;

% delta = declination of the sun (= angle between sun ray and
% equatorial plane)
% Ho = sun angle at sunrise/sunset
% calculate sun altitude, zenith, azimuth
% http://www.usc.edu/dept/architecture/mbs/tools/vrsolar/Help/solar_concepts.html
lat_r=latitude*pi/180;
sin_alt_r = cos(lat_r)*cos(delta)*cos(hour_angle)+sin(lat_r)*sin(delta);
sin2alt = sin_alt_r * sin_alt_r;
cos_alt_r = sqrt(1 - sin2alt);
sunalt = atan(sin_alt_r / cos_alt_r);
zenith = pi/2 - sunalt;

% original
% x_azm = sin(hour_angle) * cos(delta);  % adding *-1 makes it work
% y_azm = (-(cos(hour_angle))*cos(delta)*sin(lat_r))+(cos(lat_r)* sin(delta));
% azimuth = (-1 * atan2(x_azm, y_azm));

% % works, idk why
x_azm = sin(hour_angle) * cos(delta);
y_azm = (-(cos(hour_angle))*cos(delta)*sin(lat_r))+(cos(lat_r)*sin(delta));
azimuth = atan2(x_azm, y_azm);

% corripio - doesn't make sense- need to use atan2 without the +pi/2
% x_azm = -1 * sin(hour_angle) * cos(delta);
% y_azm = (sin(lat_r)*cos(hour_angle)*cos(delta)) -(cos(lat_r)*sin(delta));
% azimuth = atan( (y_azm/x_azm))+ pi/2;

% calculate instantaneous insolation
% Get eccentricity factor
% this is the reciprocal of the normalised earth's sun distance
% berger 1978 eqn (13)
% lambda (or solar longitude) is the angular distance along Earth's orbit
% measured from spring equinox (21 March)
% omega is longitude of the perihelion relative to the moving vernal equinox

eccf  = (( 1. + ecc*cos(lambda-omega) ) / (1. - ecc ^ 2))^2; % berger 1978 eq (13)

S = So * eccf;

return