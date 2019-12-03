function [ out_ela, new_ela ] = calc_ELA_slct( mb, dem )
% calc_ELA Calculate ELA from annual mass balance

% before part 2 is the same as calc_ELA
% part 2 gets the number of points in each snowline chunk
% if number is less than 50 (for now), those are ignored

[a,b,c] = size(mb); 
out_ela = zeros(c,3); 
ela = zeros(a,b,c); 

for h = 1:c
mbt = squeeze(mb(:,:,h)); 
el = zeros(a,b);
for i = 1:a-1
    for j = 1:b-1
        if (mbt(i,j) * mbt(i+1,j)) < 0
            el(i,j) = dem(i,j); 
            el(i+1,j) = dem(i+1,j); 
        elseif (mbt(i,j) * mbt(i,j+1)) < 0
            el(i,j) = dem(i,j); 
            el(i,j+1) = dem(i,j+1); 
        end
    end
end

ela(:,:,h) = (el);
end 

% part 2: 

ela_bin = ela; 
ela_bin(ela_bin > 0)=1; % in new matrix, binary is spot w ELA
new_ela = zeros(a,b,c); % will be new output ELA

for h =1:c  % for each year
comps = bwconncomp(ela_bin(:,:,h)); % get the pixel bins
temp_ela = zeros(a,b);

for i = 1:length(comps.PixelIdxList)
 v = comps.PixelIdxList(i); 
 vec = v{1} ;
 if length(vec) > 10
   temp_ela(vec) = dem(vec); 
 end
end

new_ela(:,:,h) = temp_ela; % new matrix of ELAs where area > 50
try
 out_ela(h,1) = min(nonzeros(temp_ela)) ;  
 out_ela(h,2) = mean(nonzeros(temp_ela)) ;  % mean for each year
 out_ela(h,3) = max(nonzeros(temp_ela)) ;  
catch
  out_ela(h,1) =  1680;
  out_ela(h,2) =  1680;
  out_ela(h,3) =  1680;
end

end

end

