function [I,shade] = ebm_inso_inst_calcshade(dem,rad_dem,R,Xq,Yq,csize,slope,aspect,cloud,zenith,azimuth,S)

% returns instantaneous insolation:
%   I = at surface including diffuse and direct component, depending on cloudiness (W/m^2)
%
% given:
%   sdate = matlab date number - just used for day of year with inso
%       if sdate<66 its interpreted as day of year
%   location = location structure as required by sun_position
%   dem = digital elevation model matrix
%   csize = grid cell size of dem matrix
%   slope = slope matrix corresponding to dem matrix
%   aspect = aspect matrix corresponding to dem matrix
        
% daily average insolation and info on sun angle, valid for last 5 MA
% 0 refers to kyr (thousands of years before present)
%[zenith, azimuth, S] = inso_sunpos(location.latitude,daynum);
sunalt = pi/2-zenith;  

if zenith<pi/2  % sun is above horizon
    
   %azimuth=azimuth+demangle*pi/180; % adjust if dem angle != 0
   %aspect=aspect+demangle*pi/180; % adjust if dem angle != 0
   costheta=cos(slope)*cos(zenith)+sin(slope)*sin(zenith).*cos(azimuth-aspect);
   full_shade=grd_shade(rad_dem,csize,[azimuth,sunalt]); 
   
   X = R.XWorldLimits(1):csize:R.XWorldLimits(2);
   X = X(1:end-1)+csize/2;
   Y = R.YWorldLimits(1):csize:R.YWorldLimits(2);
   Y = Y(1:end-1)+csize/2;
   shade = interp2(X,Y,full_shade,Xq,Yq); 
   
   % direct radiation
   Idir = (0.2 + 0.65 * (1-cloud)) * S .* costheta .* shade;
   % if incidence angle costheta<0, then so is Qi- not possible so set to 0
   Idir(costheta<0)=0;
   
   % diffuse radiation 
   Idif = (0.8 - 0.65 * (1-cloud)) * S * sin(pi/2 - zenith); % Oerlemans (1992)
   if cloud > 0.3
       Idif = Idif .* sin(pi/2 - zenith) .* (0.9-cloud); 
   end
   
   ta = (0.79 + 0.000024 * dem) * (1 - 0.08 * (pi/2 - sunalt)/pi/2); %transmissivity atm 
   tc = 1 - (0.41 - 0.000065 * dem) .* cloud - 0.37 * cloud .^ 2;
   
   I = ta .* tc .* (Idif + Idir);

else
    I=zeros(size(dem));  
    shade = I; 
end
  
return