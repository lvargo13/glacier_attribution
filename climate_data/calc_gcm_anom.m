function [ clim_adj,clim_adj_d,clim_vcsn,anom] = calc_gcm_anom( t,p,s,clim_force,clim_base )
% Calculate monthly and daily GCM anomalies for gridded climate data

%   Monthly anomaly calculations are similar to methods of Huss & Hock
%   2015 and Clarke et al 2015.

% Inputs: t: VCSN temp, p: VCSN precip, clim_force: GCM forcing T and P,
% clim_base: GCM historic to get anomalies from

% Outputs: clim_adj (length of months x 2; temp and precip monthly data)
% clim_adj_d (length of days x 2; temp and precip daily data)


% ---- define monthly vectors
mos = [31 28 31 30 31 30 31 31 30 31 30 31] ; 
mo_vcsn = zeros(floor(length(p)/(365/12)),1);  % get # months in days vcsn
mo_force = zeros(length(clim_force),1);  % Jan 1980 - July 2017
for i = 1:12   
 vect = i:12:(floor(length(p)/(365/12)));  
 vect_f = i:12:(length(clim_force));
 mo_vcsn(vect) = mos(i); 
 mo_force(vect_f) = mos(i);
end

moadd_vcsn = mo_vcsn;  % monthly vector for length of vcsn
for i = 2:length(mo_vcsn)
    moadd_vcsn(i) = moadd_vcsn(i)+moadd_vcsn(i-1) ; 
end
moadd_vcsn = moadd_vcsn+1; 
mo = vertcat(1,moadd_vcsn);  

moadd_force = mo_force; % monthly vector for length of forcing
for i = 2:length(mo_force)
    moadd_force(i) = moadd_force(i)+moadd_force(i-1) ; 
end
moadd_force = moadd_force+1; 
mo_force = vertcat(1,moadd_force);  

% ---- get monthly data from daily
t_mo = zeros(length(mo)-1,1);  
p_mo = zeros(length(mo)-1,1);
for i = 1:length(mo)-1
  t_mo(i) = mean(t(mo(i):mo(i+1)-1)) ; 
  p_mo(i) = mean(p(mo(i):mo(i+1)-1)) ; % mean daily precip
end

% p_mo is kg m-2 day-1 (daily mean precip for each month)
% get clim_vcsn(:,2) -> kg m-2 month-1 (mm of precip for a month)
clim_vcsn = [t_mo, p_mo]; 
for i = 1:12
  vect = i:12:length(p_mo);
  clim_vcsn(vect,2) = clim_vcsn(vect,2) .* mos(i); 
end
cvp = clim_vcsn(:,2); 
cvp(cvp==0)=min(cvp(cvp>0)); % if monthly precip = 0, make min value gt 0
clim_vcsn(:,2) = cvp; 

% ----- GCM anomalies -------
% get mean monthly values for all climate
vcsn_mean = zeros(12,2) ;
base_mean = zeros(12,2) ; 
for i = 1:12 
 vect_vcsn = i:12:length(clim_vcsn) ;
 vect_base = i:12:length(clim_base) ;
 for j = 1:2   
  vcsn_mean(i,j) = mean(clim_vcsn(vect_vcsn,j)); 
  base_mean(i,j) = mean(clim_base(vect_base,j)); 
 end 
end

% calculate climate adjustment without daily variability
clim_adj = zeros(length(clim_force),2);  
for i = 1:12
 vect = i:12:length(clim_force) ;
 clim_adj(vect,1) = vcsn_mean(i,1) + (clim_force(vect,1) - base_mean(i,1));
 clim_adj(vect,2) = vcsn_mean(i,2) .* (clim_force(vect,2) ./ base_mean(i,2));  % bc ratio, units don't matter as long as same
end


% ---- Apply gcm monthly anom to daily
% define cycle for daily variability 
period = 36 ; % number of years to repeat daily var cycle, must be divis by yrs
if (length(clim_force)/12)<period  %climscen shorter than period
    period = length(clim_force)/12; 
%elseif floor((length(clim_force)/12)/period) ~= (length(clim_force)/12)/period  
end
t_daily_cut = t(1:(period*365)); % get daily variability for length of period
p_daily_cut = p(1:(period*365));
s_daily_cut = s(1:(period*365)); % not adjusting S, but want same cycle
  
fact = ceil((length(clim_force)/12) / period) ;  % yrs in clim_force : yrs in period
t_period = repmat(t_daily_cut,fact,1);  % repeat number of times in fact
p_period = repmat(p_daily_cut,fact,1);
s_period = repmat(s_daily_cut,fact,1);

% get monthly anomalies/ratios
anom = zeros(length(clim_adj),2); 
% fact may result in vectors that are too long
vcsn_rep_t=repmat(clim_vcsn(1:(period*12),1),fact,1);
vcsn_rep_p=repmat(clim_vcsn(1:(period*12),2),fact,1);
% clip to length of required series
anom(:,1) = clim_adj(:,1) - vcsn_rep_t(1:length(clim_adj(:,1))); 
anom(:,2) = clim_adj(:,2) ./ vcsn_rep_p(1:length(clim_adj(:,2))); 

% apply anoms to t_period and p_period
clim_adj_d = zeros(length(clim_force) * 365/12, 3) ;  % gives # days in length of clim_force
for i = 1:length(clim_force)
  clim_adj_d(mo_force(i):mo_force(i+1)-1,1) = anom(i,1) + t_period(mo_force(i):mo_force(i+1)-1) ;
  clim_adj_d(mo_force(i):mo_force(i+1)-1,2) = anom(i,2) .* p_period(mo_force(i):mo_force(i+1)-1) ;
end

% clip to length of required series
clim_adj_d(:,3) = s_period(1:size(clim_adj_d,1)); 

% set up date
clim_adj_d = clim_adj_d(91:end-275,:); % 91 for april 1, -275 to march 31

end
