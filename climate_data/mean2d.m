function [ y ] = mean2d( x )
% mean2d Get mean of 1st 2 dimensions

% can be 2 or 3 dim

y = zeros(size(x,3),1); 
for i = 1:size(x,3)
    temp = squeeze(x(:,:,i));
    y(i) = nanmean(temp(:)); 
end

end

