
% script to loop through and run pdd_gcm for multiple glaciers 
% need to do a manual clear all before running
 
% glacs_in = {'rolleston/'; 'ridge/'; 'glenmary/'; 'thurneyson/'; 'southcameron/'}; 
% glacs_out = glacs_in; % this could be different, if needing to run multi simulations for each glacier without overwriting
% dem_fs = {'rolleston_10m.nc'; 'ridge_10m.nc'; 'glenmary_10m.nc'; 'thurneyson_10m.nc'; 'southcameron_10m.nc'}; 
% rad_dem_f = {'rolleston_pdd.tif'; 'ridge_pdd.tif'; 'glenmary_pdd.tif'; 'thurneyson_pdd.tif'; 'southcameron_pdd.tif'}; 
% lat = [-42.889; -43.620; -43.992; -44.165; -43.355]; 
% padj = [1.1; 1.16; 1.82; 1.81; 0.81] ; 
% tfadj = zeros(length(lat),1) + 0.75; 
% rfadj = zeros(length(lat),1) + 0.22;
% sradj = zeros(length(lat),1) + 1;
% tadj = zeros(length(lat),1) + -1.25;

% glacs_in = {'salisbury/'; 'chancellor/'; 'vertebrae12/'; 'vertebrae25/'; 'parkpass/'}; 
% glacs_out = glacs_in; % this could be different, if needing to run multi simulations for each glacier without overwriting
% dem_fs = {'salisbury_10m.nc'; 'chancellor_10m.nc'; 'vertebrae12_10m.nc'; 'vertebrae25_10m.nc'; 'parkpass_10m.nc'}; 
% rad_dem_f = {'salisbury_pdd.tif'; 'chancellor_pdd.tif'; 'vert_pdd.tif'; 'vert_pdd.tif'; 'parkpass_pdd.tif'}; 
% lat = [-43.47; -43.511; -43.32; -43.32; -44.586]; 
% padj = zeros(length(lat),1) + 0.8;
% tfadj = [0.85; 0.75; 0.85; 0.8; 0.76];
% rfadj = (tfadj*.22)/0.75 ;
% %sradj = zeros(length(lat),1) + 1;
% tadj = [-0.95; -1.25; -1.18; -1.25; -1.25];

glacs_in = {'brewster/'}; 
glacs_out = {'brewster_2011_a/'}; % this could be different, if needing to run multi simulations for each glacier without overwriting
dem_fs = {'brewster18_10m.nc'}; 
rad_dem_f = {'brewster_pdd.tif'}; 
lat = -44.072893 ; 
padj = 1.3;
tadj = -1.25;
tf = 0.88;
rf = 0.2;

type = {'gcm'; 'cesm'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(glacs_in)
for k = 1:length(type)
    
    dem_file=dem_fs{j}; % needs to be in /in_pre
    rad_dem_file = rad_dem_f{j}; % full dem for shade calc
    lat_glac = lat(j) ; 
    glac_in = glacs_in{j}; % name of glacier, e.g. 'rolleston/'
    glac_out = glacs_out{j};
    runtype = type{k};  % either 'gcm' or 'cesm'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch runtype
    case 'gcm'
        in_files = ['/Volumes/arc_03/vargola/glacier_attribution/climate_data/adjusted_GCM/' glac_in]; 
        out_files = ['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/gcm/' glac_out]; 
    case 'cesm'
        in_files = ['/Volumes/arc_03/vargola/glacier_attribution/climate_data/adjusted_CESM/' glac_in]; 
        out_files = ['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/cesm/' glac_out];          
    otherwise
        warning('Incorrect run type entered (line 11)')
end

if exist(out_files,'dir') ~= 7   % if directory does not exist, create
  mkdir(out_files)
end

% function to set up parameters
sa_pdd_parameters  
dc=CONFIG.DegreeDay;
dc.PptnFactor = padj(j);
dc.DDF = tf(j);
dc.RadiationFactor = rf(j);
%dc.SnowTempThreshold = sradj(j);
dc.TempOffset = tadj(j);

% read DEM to get size
in_pre = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/nc_input/';
topo_clim_file=[in_pre dem_file];
dem=ncread(topo_clim_file,'dem')';
dem(dem==0)=min(dem(dem>0));  % if 0 in DEM, make min DEM value
dem(dem==-9999)=min(dem(dem>0));
ice=ncread(topo_clim_file,'ice_thickness')';
easting=ncread(topo_clim_file,'easting');
northing=ncread(topo_clim_file,'northing');
dcsize=median(round(sqrt(mean(diff(easting)).^2+mean(diff(northing)).^2)));
Xq = min(easting(:)):dcsize:max(easting(:)); 
Yq = min(northing(:)):dcsize:max(northing(:));
[a,b] = size(dem);
if length(Xq) == b+1
    Xq = (min(easting(:)):dcsize:(max(easting(:))-dcsize) +dcsize/2);  
end
if length(Yq) == a+1
    Yq = (min(northing(:)):dcsize:(max(northing(:))-dcsize) +dcsize/2);  
end
Yq = Yq';  % need this for interp2

rad_dem_name = ['/Volumes/arc_03/vargola/glacier_attribution/linz_nzdem/' rad_dem_file];
[rad_dem,R] = geotiffread(rad_dem_name);
csize=R.CellExtentInWorldX;
rad_dem = flipud(rad_dem);

all_runs = dir([in_files,'*.mat']);
edate = 736421:736785; % just 365 days to define something

for ind = 1:length(all_runs)
    
% read in output of gcm_vcsn_adj
load([in_files,all_runs(ind).name]);  % load file 
disp(['running ', all_runs(ind).name ])  

run_years = 1:length(climd_adj)/365; 

mb = zeros(size(dem,1),size(dem,2),length(run_years)); 
snowth = zeros(size(dem,1),size(dem,2),length(run_years));

% generate inputdata a year at a time
st_ind = 1:365:length(climd_adj); 
end_ind = 365:365:length(climd_adj);

parfor i = 1:length(run_years)

    year_ind = st_ind(i):end_ind(i); % get vector of day of year we want (ex.1:365 yr 1, 366:730 yr 2)
    
    % interpolating precip & swr
    pptn_grid=single(zeros(size(dem,1),size(dem,2),length(year_ind)));
    swr_grid = pptn_grid; 
    for ii = 1:length(year_ind)
       pptn_grid(:,:,ii) = pptn_grid(:,:,ii) + climd_adj(year_ind(ii),2) ;
       swr_grid(:,:,ii) = swr_grid(:,:,ii) + climd_adj(year_ind(ii),3) ;
    end
    
    % apply lapse rate for temperature 
    temp_grid=lapseTempGCM(climd_adj(year_ind,1),dem,elev); 
    
    % run model
    [mb_grid,~,~,snow]=pdd_run_calcshade(dc,temp_grid,pptn_grid,swr_grid,dem,rad_dem,R,Xq,Yq,edate,csize,lat_glac);

    temp_grid=[]; % release memory
    pptn_grid=[];
    
    mb(:,:,i) = sum(mb_grid,3);  % saves annual mb (not daily tho)
    snowth(:,:,i) = snow(:,:,end);
end

ice3 = repmat(ice, [1 1 length(run_years)]);
mb(ice3<1) = NaN;
mb = mean2d(mb); 

ela = calc_ELA( dem, ice, snowth, run_years, dc.SLThreshMin, dc.SLThreshMax, 30); 

sv = [out_files,all_runs(ind).name]; 

save(sv,'mb','ela','snowth','dc','CONFIG');
end
clear pdd_run_calcshade
end
end