function [ecc,long_perh,delta,lambda] = inso(lat,day)

% Usage:
%   Fsw = daily_insolation(kyear,lat,day)
%
% Optional inputs/outputs:
%   [Fsw, ecc, obliquity, long_perh, delta, Ho] = daily_insolation(kyear,lat,day,day_type)
%
% Description:
%   Computes daily average insolation as a function of day and latitude at
%   any point during the past 5 million years.
%
% Inputs:
%   kyear:    Thousands of years before present (0 to 5000).
%   lat:      Latitude in degrees (-90 to 90).
%   day:      Indicator of time of year, by default day 1 is Jan 1.
%
% Output:
%   Fsw = Daily average solar radiation in W/m^2.
%   Can also output orbital parameters.
%   delta = declination of the sun
%   Ho = hour angle at sunrise/sunset
%
% Required file: orbital_parameter_data.mat
%
% Detailed description of calculation:
%   Values for eccentricity, obliquity, and longitude of perihelion for the
%   past 5 Myr are taken from Berger and Loutre 1991 (data from
%   ncdc.noaa.gov). If using calendar days, solar longitude is found using an
%   approximate solution to the differential equation representing conservation
%   of angular momentum (Kepler's Second Law).  Given the orbital parameters
%   and solar longitude, daily average insolation is calculated exactly
%   following Berger 1978.
%
% References: 
%   Berger A. and Loutre M.F. (1991). Insolation values for the climate of
%     the last 10 million years. Quaternary Science Reviews, 10(4), 297-317.
%   Berger A. (1978). Long-term variations of daily insolation and
%     Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
%     2362-2367.
%
% Authors:
%   Ian Eisenman and Peter Huybers, Harvard University, August 2006
%   eisenman@fas.harvard.edu

% === Get orbital parameters ===
ecc = 0.017236 ;
epsilon = 0.409209896422591 ;
omega = 4.910832916336445 ;

% For output of orbital parameters
long_perh=omega*180/pi;

% === Calculate insolation ===
% lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
% estimate lambda from calendar day using an approximation from Berger 1978 section 3

delta_lambda_m=(day-80)*2*pi/365.2422;  % original
beta=(1-ecc.^2).^(1/2);
lambda_m0=-2*( (1/2*ecc+1/8*ecc.^3).*(1+beta).*sin(-omega)-...
      1/4*ecc.^2.*(1/2+beta).*sin(-2*omega)+1/8*ecc.^3.*(1/3+beta).*(sin(-3*omega)) );
lambda_m=lambda_m0+delta_lambda_m;
lambda=lambda_m+(2*ecc-1/4*ecc.^3).*sin(lambda_m-omega)+...
  (5/4)*ecc.^2.*sin(2*(lambda_m-omega))+(13/12)*ecc.^3.*sin(3*(lambda_m-omega));

delta=asin(sin(epsilon).*sin(lambda)); % declination of the sun


