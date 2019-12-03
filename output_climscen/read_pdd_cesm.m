function [ c_pres,c_pres_sl ] = read_pdd_cesm( in_files, files )
% read_pdd_cesm Read output of PDD model for all CESM LE (not NAT)

% to get sizes
matObj = matfile([in_files,files(1).name]);
info = whos(matObj,'mb');
pres_yrs = info.size; 

pres_mb = zeros(length(files),pres_yrs(1)); 
pres_sl = zeros(length(files),pres_yrs(1));

% pres
vec = 1:length(files) ;
for i = 1:length(vec)   
 load([in_files,files(vec(i)).name]) 
 pres_mb(i,:) = mb;
 pres_sl(i,:) = ela; %sl(:,2);
end

% mean of all simulations
c_pres = reshape(pres_mb,[pres_yrs(1)*(length(files)),1]);
c_pres_sl = reshape(pres_sl,[pres_yrs(1)*(length(files)),1]);

end

