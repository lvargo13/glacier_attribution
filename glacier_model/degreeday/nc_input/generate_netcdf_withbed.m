
% generate .nc file for PDD model
% for a glacier with 1) surface dem AND 2) ice thickness

nc_out = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/nc_input/brewster_2011.nc';

% read data
%[dem,R] = geotiffread('/Volumes/arc_02/ELMER/lauren_elmer/brewster_dems/2017_surf_DEM_10m_clip.tif'); 
%[ice_thickness,R] = geotiffread('/Volumes/arc_01/vargo/Brewster/2017/2017_ice_thickness_10m_cut.tif');
[dem,~] = geotiffread('/Volumes/arc_03/vargola/eoss_images/BREWSTER/20110312_DEM_10m.tif');
%[dem,~] = geotiffread('/Volumes/arc_01/vargo/Brewster/2018/2018_d800_dem_10m.tif');
[ice_thickness,R] = geotiffread('/Volumes/arc_03/vargola/eoss_images/BREWSTER/brewster_2011_outline.tif');

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

% make -9999 in dem 0 
adjust_input = 1; % if needing to adjust values of input
[a,b] = size(dem); 
if adjust_input == 1
    for i = 1:a
        for j = 1:b
           if dem(i,j) < 0
              dem(i,j) = 0;  
           end
        end 
    end
end

% write netcdf file
ncw = 1; 
if ncw == 1
    %delete /Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/brewster_10m.nc
    %delete nc_out
    nccreate(nc_out,'dem','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'dem',flipud(dem)')
    nccreate(nc_out,'ice_thickness','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out','ice_thickness',flipud(ice_thickness)')
    nccreate(nc_out,'easting','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'easting',flipud(easting)')
    nccreate(nc_out,'northing','Dimensions',{'columns',b,'rows',a},'Format','classic')
    ncwrite(nc_out,'northing',flipud(northing)')
end
