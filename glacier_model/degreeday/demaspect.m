function asp=demaspect(dem)

% demaspect
%
% calculate dem aspect
% usage aspect = demaspect(dem)


	sd=size(dem);
	asp=zeros(sd)*NaN;
   for row=2:sd(1)-1
      for col=2:sd(2)-1;
         tzx=mean([dem(row-1:row+1,col+1) - dem(row-1:row+1,col-1)]);
         tzy=mean([dem(row+1,col-1:col+1) - dem(row-1,col-1:col+1)]);
         asp(row,col)=atan(tzx/tzy);
		end
   end
   
         
return