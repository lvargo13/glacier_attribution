%clear all
% PDD model- altered by Lauren from Brian's sample

run_years=1980:2016;

dem_file='brewster18_10m.nc';  % needs to be in /nc_input/
%dem_file='brewster_2011.nc';
rad_dem_file = 'brewster_pdd.tif' ; 

lat_glac = -44.072893 ; 
lon_glac = 169.436038 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[edate,erain,tmx,tmn,swr,elev,lon_vcsn,lat_vcsn]=sa_input_vcsn_grid(run_years,lat_glac,lon_glac);

sa_pdd_parameters  % function to set up parameters

dc=CONFIG.DegreeDay;

in_pre = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/nc_input/';
topo_clim_file=[in_pre dem_file];
dem=ncread(topo_clim_file,'dem')';
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

mb = zeros(size(dem,1),size(dem,2),length(run_years)); 
mb_day = zeros(366,length(run_years)); 
acc_day = zeros(366,length(run_years)); 
abl_day = zeros(366,length(run_years)); 

%sl_day = 31 - [14 4 9 14 3 6 30 20 12 11 14];  % mb time (2004/5 - 14/15)
sl_day = 31-[20 4 22 20 15 31 6 25 31 31 31 31 15 10 4 5 12 12 16 8 6 11 7 17 14 4 9 14 3 6 30 20 12 11 14 30 9];  % full (1980/1 - 14/15)
snowth = zeros(size(dem,1),size(dem,2),length(run_years));

% generate inputdata a year at a time
eyear=year(edate(1)):year(edate(end))-1;
parfor i = 1:length(eyear)
    % get data between April 1st and Next March 31
    year_ind=find((year(edate)==eyear(i) & month(edate)>3) | (year(edate)==eyear(i)+1 & month(edate)<=3));
    
    % interpolating precip & swr
    pptn_grid=single(zeros(size(dem,1),size(dem,2),length(year_ind)));
    swr_grid = pptn_grid; 
    for ii = 1:length(year_ind)
       pptn_grid(:,:,ii) = pptn_grid(:,:,ii) + erain(2,2,year_ind(ii)) ;
       swr_grid(:,:,ii) = swr_grid(:,:,ii) + swr(2,2,year_ind(ii)) ;
    end
    
    % apply lapse rate for temperature 
    temp_grid=lapseTemp(tmx(year_ind),tmn(year_ind),dem,elev,year_ind); 
    
    [mb_grid,acc_grid,abl_grid,snow]=pdd_run_calcshade(dc,temp_grid,pptn_grid,swr_grid,dem,rad_dem,R,Xq,Yq,edate(year_ind),csize,lat_glac);
    %[mb_grid,acc_grid,abl_grid,snow]=pdd_run(dc,temp_grid,pptn_grid,swr_grid,dem,rad_dem,R,Xq,Yq,edate(year_ind),csize,lat_glac);
    
    temp_grid=[]; % release memory
    pptn_grid=[];
     
    mb(:,:,i) = sum(mb_grid,3);  % saves annual mb  
    snowth(:,:,i) = snow(:,:,end-sl_day(i)); 
    
    for iii = 1:length(year_ind)  
        m = mb_grid(:,:,iii);
        m(ice<1) = NaN;
        mb_grid(:,:,iii) = m; 
        ac = acc_grid(:,:,iii);
        ac(ice<1) = NaN;
        acc_grid(:,:,iii) = ac;
        ab = abl_grid(:,:,iii);
        ab(ice<1) = NaN;
        abl_grid(:,:,iii) = ab;            
    end
    
    try
     mb_day(:,i) = mean2d(mb_grid);  % leap year
     acc_day(:,i) = mean2d(acc_grid);  
     abl_day(:,i) = mean2d(abl_grid);
    catch
     mb_day(:,i) = [mean2d(mb_grid); NaN];    % not leap year, need to add extra value
     acc_day(:,i) = [mean2d(acc_grid); NaN];
     abl_day(:,i) = [mean2d(abl_grid); NaN];
    end
    
end

mbday_cum = cumsum(mb_day); 
acccum = cumsum(acc_day); 
ablcum = cumsum(abl_day);
[~, st_mbday] = min(mbday_cum(1:90,:));

% mb NaN where ice < 10
ice3 = repmat(ice, [1 1 length(run_years)]);
mb(ice3<1) = NaN; 
mmb = mean2d(mb);

% % measured mb data
% measmb = [1376 691 692 -1698 -702 -74 -1728 -565 201 470 215 -1193 553 ]; % -2217];
% mberr = [247 283 403 228 211 241 187 281 418 373 314 305 215]; % 323];
% end_mbday= [11 47 9 -19 14 3 20 11 11 -7 15 6 16]; % number of days (april 1 - day mb meas)
% meas_yr = 2005:2017;
% 
% %adjust modeled date range
% mmb_dateadj = zeros(length(run_years),1); 
% if length(run_years) < 15 
%     for i = 1:length(run_years)
%      lpyr = isnan(mb_day(end,i));     
%      if lpyr == 1  % if not leapyear, end is nan
%         try
%          mmb_dateadj(i) = sum(mb_day(st_mbday(i):end-(end_mbday(i)+1),i));
%         catch 
%             % mb_day(st_mbday(i):end-1,i (end-1 because last value is nan)
%          mmb_dateadj(i) = sum([mb_day(st_mbday(i):end-1,i); mb_day(1:1+((end_mbday(i)*-1)),i+1)]);  % when end is -
%         end
%      else
%         try   
%          mmb_dateadj(i) = sum(mb_day(st_mbday(i):end-end_mbday(i),i));  % yep leap year
%         catch
%          mmb_dateadj(i) = sum([mb_day(st_mbday(i):end,i); mb_day(1:1+((end_mbday(i)*-1)),i+1)]); 
%         end
%      end
%     end
% 
%     rmse = sqrt(nanmean((mmb_dateadj' - measmb).^2));
% 
% %     figure; errorbar(meas_yr, measmb, mberr,'.','LineWidth',2.5) ; hold on
% %     plot(meas_yr, measmb, 'ko--')
% %     plot(run_years+1,mmb_dateadj,'.--','MarkerSize',10)
% %     xlabel('year'); ylabel('mass balance (mm w.e.) (date adj)')
% end

% figure; imagesc(flipud(mb)); colorbar; colormap(brewermap([],'Spectral'))
% figure; plot(mb(:,:,1),dem,'.'); xlabel('mass balance (mm w.e.)'); ylabel('elevation (m)')

% ELAs    
load('compare_clim_sl/brew_sldata')
snowth(ice3<1) = NaN;
ela = calc_ELA( dem, ice, snowth, run_years, dc.SLThreshMin, dc.SLThreshMax, 30);
%ela_rmse = sqrt(nanmean((sl_mean(end-12:end) - ela').^2));
ela_rmse = sqrt(nanmean((sl_mean - ela').^2));
% figure; plot(run_years+1,ela,'o--'); hold on; plot(sl_years,sl_mean,'o--')
% xlabel('year'); ylabel('ELA (m a.s.l.)')


