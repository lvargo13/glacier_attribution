function sombra=grd_shade(dem, csize, SunPos)

%-------------------------------------------------------------------------------
% function sombra=grd_shade(dem, csize, SunPos)
%
% original credit to Javier Corripio, documentation at end of file
%
% calculates shading on a digital elevation model
%
% output argument:
%   shade = an array the same size as the digital elevation model, with
%           value of 1 where the element is in the sun, and zero in the shade
%
% input arguments:
%   dem =   a digital elevation model (m)
%   csize = grid cell size of the dem (m)
%   sunpos = a two-element vector, [azimuth, sun_elevation] (radians)
%
% SunVector is a unit vector pointing from the surface of the Earth towards
% the sun
% coordinate system on surface of earth:
% x axis W-E
% y axis N-S
%-------------------------------------------------------------------------------
 
mn = size(dem);
m=mn(1); n=mn(2);
if m<4 || n<4
    % shading routine won't work on a grid this small so return full sun
    sombra=ones(m,n);
    return
end

% sunpos(1) = azimuth % radians
% sunpos(2) = elevation % radians
% SunVector is a unit vector pointing from the surface of the Earth towards the sun
% coordinate system on surface of earth: x axis W-E, y axis N-S, Z axis up

SunVector(1) = -sin(SunPos(1))*cos(SunPos(2));    %positive Eastwards
SunVector(2) = cos(SunPos(1))*cos(SunPos(2));   %positive Southwards
SunVector(3) = sin(SunPos(2));

% If operated with a different reference system, eg: Y is positive northwards,
%   simply make SunVector(2) = -SunVector(2)


% solvector is opposite to sunvector and used for integer steps scanning
%   along the solar direction
% this makes the largest coordinate of solvector =  1

SolVector = -SunVector/max(abs(SunVector(1:2)));

% SunVectorPXY is the horizontal projection of the sunvector

SunVectorPXY = SunVector;    % horizontal pojection sunvector, step1
SunVectorPXY(3) = 0;         % horizontal pojection sunvector, step2

% vector normal to sun in the horizontal plane
SunVectorH = cross(SunVectorPXY,SunVector);
SunVectorH = SunVectorH/sqrt(sum(SunVectorH.^2)); %** Unit SunVectorH

% unit vector normal to sun upwards
SunVectorU = cross(SunVectorH,SunVector);

%disp(['SunVector: ',num2str(SunVector)])
%disp(['SunVectorU: ',num2str(SunVectorU)])
%disp(['SolVector: ',num2str(SolVector)])


% Vector to Origin of coordinates, for projection onto plane perpendicular to sun
% It doesn't need to be (0,0,0) any arbitrary point is fine (0,0,0) is simpler

[xvalues, yvalues] = meshgrid(1:n,1:m);
xvalues = xvalues * csize;
yvalues = yvalues * csize;
zvalues = dem;

% VectortoOrigin[i] is [xvalues[i],yvalues[i],zvalues[i]

% zprojection is the projection of the DEM onto the solar plane
% which is perpendicular to SolarVector, defined by 
% SolarVectorH and SolarVectorU
% dot product [x,y,z] * SunVectorU
zprojection = xvalues * SunVectorU(1) +...
       yvalues * SunVectorU(2) +...
       zvalues * SunVectorU(3);
zprojection=zprojection';
% Determine origin of scanning lines in the direction of the sun

if (SunVector(1)<=0) Xend=n-1; end 		% sun is on the West
if (SunVector(1)>0) Xend=0; end   		% sun is on the East
if (SunVector(2)<=0) Yend = m-1; end   	% sun is on the North
if (SunVector(2)>0) Yend = 0; end       % sun is on the South
% Xend
% Yend
% SunVector
% SolVector
% sombra is the array to store binary shadow values
sombra = ones(m,n);

%------------- SCANNING CASE SUN at  0, 90, 180, 270 degrees -----------
lastloop=0;
scanindex=0;
if (abs(SunVector(1)) <= 1e-6 | abs(SunVector(2)) <= 1e-6)
	% 0 or 180 degrees
	if (abs(SunVector(1)) <= 1e-6) 
       % SolVector(2) is 1 or -1 depending on whether the sun is 0 or 180
        lastloop = n-1;
        lengthS = m-1;
        firstI = 0;
        firstJ = m-Yend-1;
        scanindex=[1,2];
        ScanVector = zeros(2,lengthS);
        ScanVector(2,:) = firstJ + [0:lengthS-1] * SolVector(2);  
	end
	
	% 90 or 270 degrees
	if (abs(SunVector(2)) <= 1e-6) 
       % SolVector(1) is 1 or -1 depending on whether the sun is 90 or 270
        lastloop = m-1;
        lengthS = n-1;
        firstI = n-Xend-1;
        firstJ = 0;
        scanindex=[2,1];
        ScanVector = zeros(2,lengthS);
        ScanVector(1,:) = firstI + [0:lengthS-1] * SolVector(1);
    end
    
   for i = 0:lastloop
        ScanVector(scanindex(1),:) = i;
        % compare z projection
	%______________________________________________________________
	% Array of Z_Projections along scan line
    zpcompare =-1e30;
    ZPdxdy=zeros(length(ScanVector));
    for k=1:length(ScanVector)
        ZPdxdy(k) = zprojection(ScanVector(1,k)+1,ScanVector(2,k)+1);
        % Some antialiasing device is required.  A temporary improvement is
        %    achieved by averaging the cell projections with the previous value
        %    along the scan line
        if k>1
            ZPdxdy(k) = (ZPdxdy(k) + ZPdxdy(k-1))/2;
        end
        if (ZPdxdy(k)<zpcompare)
            sombra(ScanVector(2,k)+1,ScanVector(1,k)+1) = 0;  % this line sucks up 80% of time
        else
            zpcompare = ZPdxdy(k);
        end
    end
   
   %_____________________________________________________________
   end
else
%-------------------------  SCANNING BLOCK I ----------------------------------
	
	%***************** scanning along X-AXIS  ****************
    for i = 0:n-1
        % length of scanning line (vector)
        lengthX1 = abs(Xend - i) / SolVector(1); 	% solve for: l*Solvector=Xmax
        lengthX2 = (m-1) / SolVector(2);     		% solve for: l*Sv=Ymax
        lengthS = floor(min(abs([lengthX1,lengthX2])));

        if (lengthS==1) continue; end % this means a corner = sunny

        vectorlength = (0:lengthS-1);

        % X-origin of scan line
        firstI = i;
        % Y-origin of scan line
        firstJ = m-Yend-1;  %bma

        ScanVector=[];
        % ScanVectorX: arrays of xvalues along scan line  ScanVectorY: Yvalues
        ScanVector(1,:) = floor(firstI + vectorlength * SolVector(1));  
        % ScanVectorY: Yvalues
        ScanVector(2,:) = floor(firstJ + vectorlength * SolVector(2));
        % compare z projection
        %______________________________________________________________

        % Array of Z_Projections along scan line
        zpcompare =-1e30;
        ZPdxdy=zeros(length(ScanVector));
        for k=1:length(ScanVector)
            ZPdxdy(k) = zprojection(ScanVector(1,k)+1,ScanVector(2,k)+1);
            % Some antialiasing device is required.  A temporary improvement is
            %    achieved by averaging the cell projections with the previous value
            %    along the scan line
            if k>1
                ZPdxdy(k) = (ZPdxdy(k) + ZPdxdy(k-1))/2;
            end
            if (ZPdxdy(k)<zpcompare)
                sombra(ScanVector(2,k)+1,ScanVector(1,k)+1) = 0;  % this line sucks up 80% of time
            else
                zpcompare = ZPdxdy(k);
            end
        end

        %_____________________________________________________________
    end
	
   Xend = abs(Xend-i);
	
	%==============================================================================
	%-------------------------  SCANNING BLOCK J ----------------------------------
	
	%*** scanning   Along Y-AXIS
	
    for j = 0:m-1
        %length of scanning line (vector)
        lengthY1 = abs((Yend - j) / SolVector(2));    	% solve for: l*Sv=Ymax 
        lengthY2 = (n-1) / SolVector(1); 				% solve for: l*Sv=Xmax
        lengthS = floor(min(abs([lengthY1,lengthY2])));
%        disp([j lengthY1 lengthY2 lengthS])
        if (lengthS==1) continue; end % this means a corner = sunny

        vectorlength = (0:lengthS-1);
        % X-origin of scan line
        firstI = Xend;
        % Y-origin of scan line
        firstJ =  j;
        ScanVector=[];
        % ScanVectorX: arrays of xvalues along scan line
        ScanVector(1,:) = floor(Xend + vectorlength * SolVector(1));  

        % ScanVectorY: Yvalues
        ScanVector(2,:) = floor(firstJ + vectorlength * SolVector(2));
        % compare z projection
        %______________________________________________________________

        % Array of Z_Projections along scan line
        zpcompare =-1e30;
        for k=1:length(ScanVector)
            ZPdxdy(k) = zprojection(ScanVector(1,k)+1,ScanVector(2,k)+1);
            % Some antialiasing device is required.  A temporary improvement is
            %    achieved by averaging the cell projections with the previous value
            %    along the scan line
            if k>1
                ZPdxdy(k) = (ZPdxdy(k) + ZPdxdy(k-1))/2;
            end
            if (ZPdxdy(k)<zpcompare)
                sombra(ScanVector(2,k)+1,ScanVector(1,k)+1) = 0;  % this line sucks up 80% of time
            else
                zpcompare = ZPdxdy(k);
            end
        end

        %_____________________________________________________________
    end
end

% border columns and rows are set to 1st (or 2nd) inner column and row
sombra(1,:) = sombra(2,:);
sombra(m,:) = sombra(m-1,:);
sombra(:,1) = sombra(:,2);
sombra(:,n-1) = sombra(:,n-1);

return

%   Corripio, J.G. 2003: "Vectorial algebra algorithms for calculating...
%   Written by:  Javier G. Corripio. March, 2002
