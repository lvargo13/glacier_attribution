
% PDD model- altered by Lauren from Brian's sample
% need to clear all if run for different glacier- global geometry variables
% inputs: glacier name, run years
% outputs: mb, snowlines for run years
% this script isnt used for attribution calculations- just model calibration

run_years=1980:2016;
glac = 'vertebrae25';
pltsl = 0; % 1 to plot snowlines 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dem_file = [glac '_10m.nc'];  % needs to be in /nc_input/
rad_dem_file = [glac '_pdd.tif'] ; 

switch glac
    case 'rolleston'
        lat_glac = -42.889 ;
        lon_glac = 171.527 ;
        SL = [1760 1773 1640 1755 1755 1776 1778 1762 1770 NaN NaN NaN 1620 1745 1620 1750 1640 1803 1850 1868 1755 1860 1669 1762 1740 1800 1795 1815 1814 1768 1865 1813 1765 1800 1765 1865 1765]; % rolleston starts 81
    case 'thurneyson'
        lat_glac = -44.165 ;  % thurneyson
        lon_glac = 169.598 ;
        SL = [1943 NaN 1865 1882 1905 1938 1918 1950 1970 NaN NaN 1930 1910 1904 1868 1938 1900 1965 2112 2132 1878 2105 1873 1935 1900 1980 1963 2110 2100 1975 2150 2110 1938 1915 2080 2120 1938]; %thurneyson starts '81
    case 'kahutea'
        lat_glac = -43.01 ;  % kahutea
        lon_glac = 171.376 ;
    case 'southcameron'
        lat_glac = -43.355 ;  % southcameron
        lon_glac = 170.988 ;
        SL = [NaN 2290 2150 2130 2150 2250 NaN 2275 2340 NaN NaN NaN 2150 2250 2130 2250 2135 2340 2400 2365 2220 2300 NaN 2275 2220 2340 2260 2400 2290 2365 2550 2300 2270 2250 2280 2450 2250] ; %S Cameron
    case 'vertebrae12'
        lat_glac = -43.32 ;  % vertebrae
        lon_glac = 170.615 ;
        SL = [1813 1850 1800 NaN 1808 1820 NaN 1871 1825 NaN NaN 1768 1778 1807 1778 1807 1791 1855 2085 1985 1805 1995 1816 1807 1790 1830 1856 2030 1924 1840 2090 2020 1840 1830 1870 2030 1830]; %vert12
    case 'vertebrae25'
        lat_glac = -43.32 ;  % vertebrae
        lon_glac = 170.615 ;
        SL = [1802 1840 1778 NaN 1790 1813 NaN 1843 1820 NaN NaN 1746 1756 1786 1756 1786 1770 1835 1965 1910 1807 1920 1789 1795 1765 1834 1824 1950 1859 1834 1995 1940 1834 1834 1860 1950 1834]; % vert25
    case 'salisbury'    
        lat_glac = -43.47 ;  % salisbury
        lon_glac = 170.22 ;
        SL = [1752 1827 1718 1760 1734 1809 1775 1729 1744 NaN NaN 1681 1710 1772 1645 1752 1726 1852 2030 1982 1715 1860 1715 1732 1715 1850 1810 1950 1950 1780 2095 1980 1780 1780 1780 1980 1780]; % salisbury
    case 'parkpass'
        lat_glac = -44.586 ;  % parkpass
        lon_glac = 168.238 ;
        SL = [1778 1858 1762 1765 1702 1863 1843 NaN 1794 NaN NaN NaN NaN 1783 1635 1808 1745 1880 1955 1943 1748 1910 1661 1665 1670 1900 1850 1910 1910 1880 2005 1940 1850 1843 1910 1943 1843]; % parkpass, start 81
    case 'glenmary'
        lat_glac = -43.992 ;  % glenmary
        lon_glac = 169.883 ;
        SL = [2130 2181 2135 2020 2108 NaN NaN 2180 2180 NaN NaN 2160 2107 2145 2140 2130 2045 2195 2290 2245 2138 2210 2115 2145 2110 2205 2190 2280 2210 2180 2305 2245 2170 2150 2210 2305 2170]; % glenmary start 81
    case 'ridge'
        lat_glac = -43.620 ;  % ridge
        lon_glac = 170.360 ;
        SL = [2228 2236 2211 2085 2194 2217 NaN 2285 2277 NaN NaN NaN 2100 2170 2132 2138 2090 2235 2335 2305 2110 2300 2163 2228 2165 2310 2256 2325 2305 2280 2490 2325 2235 2230 2235 2490 2230]; %ridge
    case 'chancellor'
        lat_glac = -43.511 ;  % chancellor
        lon_glac = 170.125 ;
        SL = [1663 1848 1545 NaN 1609 1678 1728 1808 1678 NaN NaN NaN 1550 1609 1545 1720 1580 1848 1965 1960 1570 1860 1570 1605 1575 1850 1735 1865 1865 1735 2000 1865 1609 1735 1750 1965 1735]; % chancellor
end

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
snowth = zeros(size(dem,1),size(dem,2),length(run_years));
tgrid = snowth; 

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
    
    tgrid(:,:,i) = mean(temp_grid,3); 
    temp_grid=[]; % release memory
    pptn_grid=[];
     
    mb(:,:,i) = sum(mb_grid,3);  % saves annual mb  
    snowth(:,:,i) = snow(:,:,end);
    
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
snowth(ice3<1) = NaN; 
mmb = mean2d(mb); 

% measured rolleston mb data
% measmb = [-2039 -429 740 -38 657 -1006]; 
% end_mbday= [7 3 17 9 11 9]; % number of days (april 1 - day mb meas) (email from heather 9/9/19)

% adjust modeled date range
mmb_dateadj = zeros(length(run_years),1); 
if length(run_years) < 15 
    for i = 1:length(run_years)
     lpyr = isnan(mb_day(end,i));     
     if lpyr == 1  % if not leapyear, end is nan
        try
         mmb_dateadj(i) = sum(mb_day(st_mbday(i):end-(end_mbday(i)+1),i));
        catch 
            % mb_day(st_mbday(i):end-1,i (end-1 because last value is nan)
         mmb_dateadj(i) = sum([mb_day(st_mbday(i):end-1,i); mb_day(1:1+((end_mbday(i)*-1)),i+1)]);  % when end is -
        end
     else
        try   
         mmb_dateadj(i) = sum(mb_day(st_mbday(i):end-end_mbday(i),i));  % yep leap year
        catch
         mmb_dateadj(i) = sum([mb_day(st_mbday(i):end,i); mb_day(1:1+((end_mbday(i)*-1)),i+1)]); 
        end
     end
    end

figure; plot(run_years+1,mmb_dateadj,'o--');hold on  % date adjusted MB
xlabel('year'); ylabel('Mass Balance (mm w.e.)')
plot(run_years+1,measmb,'o--')
rmse = sqrt(nanmean((mmb_dateadj' - measmb).^2));
end

% calculate and plot snowline
SL(SL<min(dem(ice>0))) = min(dem(ice>0));
SL(SL>max(dem(ice>0))) = max(dem(ice>0));
ela = calc_ELA( dem, ice, snowth, run_years, dc.SLThreshMin, dc.SLThreshMax, 30);

if pltsl == 1
    figure;plot(run_years+1,ela,'o--'); hold on; plot(run_years+1,SL,'o--')
    ylabel('ELA (m a.s.l.)')
    legend('modeled ELA', 'measured EoSS')
    rmse = sqrt(nanmean((SL - ela').^2));
end

