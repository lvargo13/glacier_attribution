function [ind,d]=findclose(num,arr,sgl,tol)

% function [ind,d]=findclose(num,arr,sgl,tol)
%
% find the closest values in an array
%
% num is the number you're loooking for
% arr is the array
% tol is how close you want the values to be
% sgl is 1 if you want just one index for the result
%
% ind is the ind in arr of the nearest value to num
% d is the distance from the nearest value

if nargin<3 sgl=0; end
if nargin<4 tol=0; end
if isnan(num) ind=[]; return; end
if all(isnan(arr)) ind=[]; return; end
diff=abs(arr-num);
md=min(diff);
if tol==0
    ind=find(diff==md);
else
    ind=find(diff<tol);
end
d=arr(ind)-num;

if sgl && ~isempty(ind)
    ind=ind(1);
    d=d(1);
end

return