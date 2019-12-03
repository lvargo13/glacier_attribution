
% ---- read in data

ELAeq = 1935;  % mean ELA
ELAi = 1740;  % ELA year i 
acc_grad = 7.1;  % mm w.e. m-1
abl_grad = 14.9;  % mm w.e. m-1

% dem
dem_file='brewster_10m_t2.nc';  % needs to be in /nc_input/
in_pre = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/nc_input/';
topo_clim_file=[in_pre dem_file];
dem=ncread(topo_clim_file,'dem')';
ice=ncread(topo_clim_file,'ice_thickness')';
max_elev = max(dem(ice>1)); 
min_elev = min(dem(ice>1)); 

bins = 37;
figure; h = histogram(dem(ice>1),bins); 

% --- calculate mb
mb_elev = zeros(1,h.NumBins);
m_elev = zeros(1,h.NumBins); 
for i = 1:h.NumBins  % for each elev band
    m_elev(i) = (h.BinEdges(i) + h.BinEdges(i+1))/2 ; % mean elev in band
     
end
mb_abl = (m_elev - ELAi) .* abl_grad; 
mb_acc = (m_elev - ELAi) .* acc_grad; 

%mb_bins(i) = mb_elev(i) * h.Values(i);  % sum of mb at that elevation band

%figure; plot(mb_elev,m_elev,'.')

figure; plot(mb_abl,m_elev,'.'); hold on
plot(mb_acc,m_elev,'.')

