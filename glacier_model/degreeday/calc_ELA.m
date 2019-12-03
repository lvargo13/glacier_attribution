function [ ela ] = calc_ELA(dem,ice,snowth,run_years,min_lim,max_lim,sens_val)
% calc_ELA Calculate ELA from annual mass balance
 
% in: mb grid, dem, sens val (the number of grid boxes for a chunk of ELA
% to have to be considered (higher value like 50 means only the long lines
% are ELA, lower value like 5 mean the small chunks of - to + MB are
% considered with the larger)
dem(dem==0)=min(dem(dem>0));  % if 0 in DEM, make min DEM value

% create a snow thickness binary image, where only snow thickness within
% wanted range is =1, everything else (lower, higher, NaN) is 0
 
 snowth_binary = zeros(size(snowth)); 
 sl_below = zeros(length(run_years),1); 
 sl_above = zeros(length(run_years),1);
for ii = 1:length(run_years)
    sb = snowth(:,:,ii); 
    si = snowth(:,:,ii);
    si(ice<1)=0; 
    snowth_inds = find( (si>min_lim) & (si<max_lim));
    if length(snowth_inds)>0
        sb(si<min_lim)=0;
        sb(si>max_lim)=0;
        sb(sb>0)=1;
    elseif min(si(ice>0)) > max_lim   % if min snowthickness is gt. max, snow everywhere
        sl_below(ii) = 1; 
    elseif max(si(:)) < min_lim
        sl_above(ii) =1; 
    end
    snowth_binary(:,:,ii) = sb;
    %figure; imagesc(snowth_binary(:,:,ii))
end

% loop through each year
% for each year, comps is groups of pixels (snow thickness in wanted range)
% if size of group is larger than sens_val, then the ELA is calculated with
% those values
ela = zeros(length(run_years),1); 
for ii = 1:length(run_years)
 comps = bwconncomp(snowth_binary(:,:,ii)); % get the pixel bins
 temp_ela = []; 
 
 % check that there are groups with enough pixels, if none, sens_val is now
 % largest group
 check_len = max(cellfun('length',comps.PixelIdxList));
 if check_len < sens_val
     sens_val = check_len-1; 
 end
 
for i = 1:length(comps.PixelIdxList)
 v = comps.PixelIdxList(i); 
 vec = v{1}; 
 if length(vec) >= sens_val
   temp_ela = [temp_ela; dem(vec)]; 
 end
end

ela(ii) = mean(temp_ela); 
sens_val = 50; 
end

% use data from first loop, if ELA is all above or below
for ii = 1:length(run_years)
    if sl_below(ii) == 1
        ela(ii) = min(dem(ice>0)); 
    elseif sl_above(ii) == 1
        ela(ii) = max(dem(ice>0)); 
    end
end

end

