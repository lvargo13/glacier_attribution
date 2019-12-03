function [ c_pres,c_past,c_pres_sl,c_past_sl ] = read_pdd_gcm( in_files,files )
% read_pdd_gcm Read output of PDD model for all GCMs 

n = 2; % number of climate scenarios, if this changes, other stuff needs to

% to get sizes
matObj = matfile([in_files,files(1).name]);
info = whos(matObj,'mb');
past_yrs = info.size; 
matObj = matfile([in_files,files(2).name]);
info = whos(matObj,'mb');
pres_yrs = info.size; 

past_mb = zeros(length(files)/n,past_yrs(1)); 
pres_mb = zeros(length(files)/n,pres_yrs(1)); 
past_sl = zeros(length(files)/n,past_yrs(1)); 
pres_sl = zeros(length(files)/n,pres_yrs(1));

% past
vec = 1:n:length(files) ;
for i = 1:length(vec)   
 load([in_files,files(vec(i)).name])
 past_mb(i,:) = mb; 
 %past_sl(i,:) = sl(:,2); 
  past_sl(i,:) = ela;
end

% pres
vec = 2:n:length(files) ;
for i = 1:length(vec)   
 load([in_files,files(vec(i)).name]) 
 pres_mb(i,:) = mb;
 %pres_sl(i,:) = sl(:,2);
 pres_sl(i,:) = ela;
end

% mean of all simulations
c_past = reshape(past_mb,[past_yrs(1)*(length(files)/2),1]);
c_pres = reshape(pres_mb,[pres_yrs(1)*(length(files)/2),1]);
c_past_sl = reshape(past_sl,[past_yrs(1)*(length(files)/2),1]);
c_pres_sl = reshape(pres_sl,[pres_yrs(1)*(length(files)/2),1]);

end

