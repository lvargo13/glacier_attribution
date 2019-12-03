function [mb,acc,abl,snow_th]=pdd_run_calcshade(dc,temp,pptn,swr,dem,rad_dem,R,Xq,Yq,time,csize,lat_glac)

% takes parameters in CONFIG.DegreeDay and runs a positive degree-day model.
% temp and pptn are m x n x t matices where m and n are spatial dimensions
% and t is the number of timesteps.
% time is absolute time, and dH is change in ice thickness from that in topo_clim_file
% -- This will calculate shade tha toa S for the year, and output to sa_pdd

persistent slope aspect 

snow_depth=(dem-dc.SnowLineElevation)*dc.SnowGradient;
snow_depth(snow_depth<0)=0;
snow_depth(snow_depth>dc.SnowMaxThickness)=dc.SnowMaxThickness;

clear ebm_inso_instant

num_timesteps=size(temp,ndims(temp));

if isempty(slope)
    slope=demslope(dem,csize);
    aspect=demaspect(dem);
end

% set up matrices the correct size
acc=zeros(size(temp));
abl=zeros(size(temp));
mb=zeros(size(temp));
snow_th=zeros(size(temp));
ptt=zeros(size(temp,1),size(temp,2)); 
Ta = zeros(size(temp,1),size(temp,2));  % daily temp > 0, accumulated since snowfall, reset w snow

temp=temp+dc.TempOffset; % apply uniform temperature offset

for ii=1:num_timesteps
   
    temp_step=temp(:,:,ii);
    pptn_step=pptn(:,:,ii);
    swr_step=swr(:,:,ii); 
    
    % calc accumulation (temp_step<dc.SnowTempThreshold is 1 yes, 0 not)
    accum_raw = pptn_step .* (temp_step<dc.SnowTempThreshold) * dc.PptnFactor ;
    acc(:,:,ii) = accum_raw;  
    snow_depth=snow_depth+(accum_raw); % keep snowthickness up to date
    
    % calculate positive temperature sums
    ptt=temp_step.*(temp_step>0); 
    
    % calculate albedo (brock 2000 and pellicciotti 2005)
    Ta(accum_raw==0)=Ta(accum_raw==0)+ptt(accum_raw==0);  % accumulate Ta
    Ta(accum_raw>0)=0;  % if accum, reset Ta to 0
    albedo = dc.AlbedoSnow - dc.AlbedoCoef * log10(Ta+1) ; % give matrix of albedo
    albedo(snow_depth==0) = dc.AlbedoIce ;  
    
    % --- radiation
    % TOA clear sky radiation (hourly temporal res) (for cloud factor)
    ts = 24;
    Sa = zeros(size(temp,1),size(temp,2), ts);
    za = zeros(ts,1); 
    aa = zeros(ts,1); 
    vec = round((0:(24/24):23)/24,5); 
    for i = 1:ts
        [Sa(:,:,i),za(i),aa(i)] = ebm_TOAinso_instant(time(ii)+(vec(i)),lat_glac,dem);
    end
    swr_cf = mean(Sa,3);  % mean daily modeled SWR (not TOA, for cloud factor)
    
    % Cloud factor ( can only be daily temporal res- vcsn swr is daily)
    cloudf = swr_step ./ swr_cf ; 
    cloudf = mean(cloudf(:));  % not spatially variable
    if cloudf > 1 
        cloudf = 0.99;
    end
    cloud = 1-cloudf;  % now low cloud = clear sky
    
    % Surface radiation
    radd = zeros(size(temp,1),size(temp,2), ts);
    %sh = radd;
    for i = 1:ts
      [radd(:,:,i),~] = ebm_inso_inst_calcshade(dem,rad_dem,R,Xq,Yq,csize,slope,aspect,cloud,za(i),aa(i),Sa(:,:,i));
    end
    rad = mean(radd,3);

    ablation=(ptt .* dc.DDF) + (rad .* (1-albedo) * dc.RadiationFactor);
    snow_depth=snow_depth-ablation;
    snow_depth(snow_depth<0)=0;  % make negative thickness 0
    snow_th(:,:,ii) = snow_depth; 
    
    abl(:,:,ii)=ablation;
    mb(:,:,ii)=accum_raw-ablation;
end

return