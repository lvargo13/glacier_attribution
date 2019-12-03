% script to loop through and run pdd_gcm for multiple param combos (1 glac)
% need to do a manual clear all before running (whos global )

% input: 
gla = 'ridge'; 
lat = -43.620 ;  % ridge


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indat = [gla '_psuiteSL.mat'];
load(['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test_postreview/',indat])
glacs_in = [gla '/'];
dem_fs = [gla '_10m.nc'];  % needs to be in /nc_input/
rad_dem_f = [gla '_pdd.tif'] ; 

type = {'gcm'; 'cesm'}; 

for j = 1:19; %1:length(glacs_out)
for k = 1:length(type)
    
    dem_file=dem_fs; % needs to be in /in_pre
    rad_dem_file = rad_dem_f; % full dem for shade calc
    lat_glac = lat ; 
    glac_in = glacs_in; % name of glacier, e.g. 'rolleston/'
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
dc.PptnFactor = pad(j);
dc.TempOffset = tad(j);
dc.DDF = tf(j);
dc.RadiationFactor = rf(j);

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