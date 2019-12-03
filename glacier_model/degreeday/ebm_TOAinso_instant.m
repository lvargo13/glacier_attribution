function [Sa,zenith,azimuth] = ebm_TOAinso_instant(sdate,location,dem)

% returns instantaneous TOA clear sky radiation:
% take from ebm_inso_instant


daynum=dayofyear(sdate);

% daily average insolation and info on sun angle
[zenith, azimuth, S] = inso_sunpos(location,daynum);
sunalt = pi/2-zenith;  
ta = (0.79 + 0.000024 * dem) * (1 - 0.08 * (pi/2 - sunalt)/pi/2); %transmissivity atm

if zenith<pi/2  % sun is above horizon
   
   Sa = S .* cos(zenith) .* ta;

else
    Sa=zeros(size(dem));   
end
  
return