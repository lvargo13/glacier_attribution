function slp=demslope(dem,csize)

% demslope
%
% calculate slope from a matrix dem
% usage slope=demslope(dem,csize)
% where csize is an optional cell size

	if nargin==1 csize=1; end
	sd=size(dem);
	slp=zeros(sd)*NaN;
   for row=2:sd(1)-1
      for col=2:sd(2)-1;
%         tzx=mean([dem(row-1:row+1,col+1) - dem(row-1:row+1,col-1)])/csize;
%         tzy=mean([dem(row+1,col-1:col+1) - dem(row-1,col-1:col+1)])/csize;
%         slp(row,col)=atan(sqrt(tzx^2+tzy^2)/2);
         tzx1=(dem(row-1,col+1)+2*dem(row,col+1)+dem(row+1,col+1))/(8*csize);
         tzx2=(dem(row-1,col-1)+2*dem(row,col-1)+dem(row+1,col-1))/(8*csize);
         tzy1=(dem(row+1,col-1)+2*dem(row+1,col)+dem(row+1,col+1))/(8*csize);
         tzy2=(dem(row-1,col-1)+2*dem(row-1,col)+dem(row-1,col+1))/(8*csize);
         slp(row,col)=atan(sqrt((tzx1-tzx2)^2+(tzy1-tzy2)^2));
   
      end
   end
         
return