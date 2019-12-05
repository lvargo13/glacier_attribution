function [ outT ] = lapseTemp( tmx,tmn,dem,base_elev,yr )
% Calculate mean temp from min and max temps with seasonal lapse rates
%   tmx and tmx both scalars
%   dem is matrix
%   base elev is reference elevation for tmx and tmn
%   i is the day of the year

% lapse rates from Tait & Maraca and Norton
lr_mx_djf = -0.0061 ; % all in C m-1
lr_mx_mam = -0.0063 ;
lr_mx_jja = -0.0064 ;
lr_mx_son = -0.0066 ;
lr_mn_djf = -0.0041 ; 
lr_mn_mam = -0.0032 ;
lr_mn_jja = -0.0030 ;
lr_mn_son = -0.0042 ;

vec1 = [1:61,336:length(yr)];  % april may, next march
vec2 = [62:153];  % june july aug
vec3 = [154:244]; % sept oct nov
vec4 = [245:335]; % dec jan feb 

tmx_adj = zeros(size(dem,1), size(dem,2), length(yr)); 
tmn_adj = zeros(size(dem,1), size(dem,2), length(yr)); 

for i = 1:length(yr)
    if ismember(i,vec1)
           tmx_adj(:,:,i) = tmx(i) + (dem-base_elev) * lr_mx_mam ;
           tmn_adj(:,:,i) = tmn(i) + (dem-base_elev) * lr_mn_mam ;
    elseif ismember(i,vec2)
           tmx_adj(:,:,i) = tmx(i) + (dem-base_elev) * lr_mx_jja ;
           tmn_adj(:,:,i) = tmn(i) + (dem-base_elev) * lr_mn_jja ;
    elseif ismember(i,vec3)
           tmx_adj(:,:,i) = tmx(i) + (dem-base_elev) * lr_mx_son ;
           tmn_adj(:,:,i) = tmn(i) + (dem-base_elev) * lr_mn_son ;
    elseif ismember(i,vec4)
           tmx_adj(:,:,i) = tmx(i) + (dem-base_elev) * lr_mx_djf ;
           tmn_adj(:,:,i) = tmn(i) + (dem-base_elev) * lr_mn_djf ;
    end
end

outT = zeros(size(dem,1), size(dem,2), length(yr)); 
for i = 1:length(yr)
 outT(:,:,i) = ((tmx_adj(:,:,i)-273.15) + (tmn_adj(:,:,i)-273.15)) ./2 ;
end

end

