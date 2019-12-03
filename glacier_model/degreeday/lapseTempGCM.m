function [ t_adj ] = lapseTempGCM( t,dem,base_elev )
% Lapse temp to elevation 
%   t scalar
%   dem is matrix
%   base elev is reference elevation for tmx and tmn

% mean of min and max lapse rates from Tait & Maraca and Norton
lr_djf = (-0.0061+-0.0041)/2 ; % all in C m-1
lr_mam = (-0.0063+-0.0032)/2 ;
lr_jja = (-0.0064+-0.0030)/2 ;
lr_son = (-0.0066+-0.0042)/2 ;

% time starts april 1 and ends march 31
vec1 = [1:61,335:365];  % april may, next march
vec2 = [62:153];  % june july aug
vec3 = [154:244]; % sept oct nov
vec4 = [245:334]; % dec jan feb


t_adj = zeros(size(dem,1), size(dem,2), 365);
for i = 1:365
    if ismember(i,vec1)
           t_adj(:,:,i) = t(i) + (dem-base_elev) * lr_mam ;
    elseif ismember(i,vec2)
           t_adj(:,:,i) = t(i) + (dem-base_elev) * lr_jja ;
    elseif ismember(i,vec3)
           t_adj(:,:,i) = t(i) + (dem-base_elev) * lr_son ;
    elseif ismember(i,vec4)
           t_adj(:,:,i) = t(i) + (dem-base_elev) * lr_djf ;
    end
end

end

