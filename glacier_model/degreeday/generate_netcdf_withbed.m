
% generate .nc file for PDD model
% for a glacier with 1) surface dem AND 2) ice thickness

nc_out = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/brewster_10m.nc';

% read data
[dem,R] = geotiffread('/Volumes/arc_02/ELMER/lauren_elmer/brewster_dems/2017_surf_DEM_10m_clip.tif'); 
[ice_thickness,R] = geotiffread('/Volumes/arc_01/vargo/Brewster/2017/2017_ice_thickness_10m_cut.tif');

r= R.RasterSize(1);
c= R.RasterSize(2);
x1 = R.XLimWorld(1); 
x2 = R.XLimWorld(2);
y1 = R.YLimWorld(1);
y2 = R.YLimWorld(2);

x = linspace(x1,x2,c); 
y = linspace(y1,y2,r); 
[easting,northing] = meshgrid(x,y); 
northing = flipud(northing); 

% make -9999 in dem 1950 (just to see if it changes shade calc
% get rid of not values and make all thicknesses bt -1 and 2 = 2
nv = ice_thickness(1,1) ;  % make sure this is not a value
min_th = 20 ;
for i = 1:r
    for j = 1:c
       if dem(i,j) < 0
          dem(i,j) = 1950;  
       end
        
       if ice_thickness(i,j) == nv ; 
           ice_thickness(i,j) = 0 ;
       elseif (ice_thickness(i,j) > -1) && (ice_thickness(i,j) < min_th)
           ice_thickness(i,j) = min_th ;
       end 
    end 
end



dem = dem(21:235,55:239); 
ice_thickness = ice_thickness(21:235,55:239);
easting = easting(21:235,55:239);
northing = northing(21:235,55:239);
[a,b] = size(dem);

% write netcdf file
ncw = 1; 
if ncw == 1
    delete /Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/brewster_10m.nc
    nccreate(nc_out,'dem','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'dem',flipud(dem)')
    nccreate(nc_out,'ice_thickness','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out','ice_thickness',flipud(ice_thickness)')
    nccreate(nc_out,'easting','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'easting',flipud(easting)')
    nccreate(nc_out,'northing','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'northing',flipud(northing)')
end
