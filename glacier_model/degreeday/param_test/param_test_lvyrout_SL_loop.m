% after running run_param_test
% for leave-one-year-out calibration
% for years left out (yo), calculates and prints RMSE for that year for new calibration
% prints out mean annual and annual snowline RMSE, and RMSE for year i, when there are mulitple param
% options that are good, use lowest annual RMSE 

% inputs: 1) glacier (ex. gla = 'rolleston';)
% 2) index of best run in parameter sweep (param_test_SL.m)


gla = 'glenmary';
ind = 1457; % index of best run (from param_test_SL.m)


%%%%%%%%%%%%%%%%%%%%%%%

% load input
fils_path = ['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/' gla 'SL/']; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);
load([gla '_SL.mat']); % measured snowline data

yo = 1:37;  % years left out
yr_var = zeros(37,3); 

for j = 1:length(yo)
    
yrout = yo(j);     
meassl = SL;
mw_m = meassl(yrout); 
meassl(yrout) = [];  % delete year out
m = nanmean(meassl);

mmb_all = zeros(length(fils),1);
mmb_rms = zeros(length(fils),1);
rmserr = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
mw = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils(i).name]) 
   mw(i) = ela(yrout);
   ela(yrout)=[]; 
   mmb_all(i) = mean(ela);
   mmb_rms(i) = sqrt( ((nanmean(ela) -m).^2));
   rmserr(i) = sqrt(nanmean((ela' - meassl).^2));
end

l_mb = min(rmserr);
lf_mb = find(rmserr==l_mb);
l_mmb = min(mmb_rms); 
lf_mmb = find(mmb_rms==l_mmb);

if lf_mb == lf_mmb   
    C=lf_mb;
else 
   lf_mb = find(rmserr<l_mb+1); 
   lf_mmb = find(mmb_rms<l_mmb+1);
   C = intersect(lf_mb,lf_mmb);
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+2); 
       lf_mmb = find(mmb_rms<l_mmb+2);
       C = intersect(lf_mb,lf_mmb);
   end
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+3); 
       lf_mmb = find(mmb_rms<l_mmb+3);
       C = intersect(lf_mb,lf_mmb);
   end
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+4); 
       lf_mmb = find(mmb_rms<l_mmb+4);
       C = intersect(lf_mb,lf_mmb);
   end
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+5); 
       lf_mmb = find(mmb_rms<l_mmb+5);
       C = intersect(lf_mb,lf_mmb);
   end   
end

if length(C) > 1  % get C to 1
    if C == 0
        'bad'
    else
    [n,m] = min(mmb_rms(C)); % find min 
    C = C(m);  
    end
end

yr_var(j,2) = rmserr(C); 
yr_var(j,3) = mmb_rms(C);
yr_var(j,1) = sqrt((mw_m-mw(C)).^2); 
end


% load file of best model run parameters to compare
fils_path = ['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/' gla 'SL/']; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);
load([fils_path fils(ind).name])
mod_SL = sqrt((SL' - ela).^2); 

% now get variation for each year, for best calib vs best calib when year
% left out
outp = sqrt((mod_SL - yr_var(:,1)).^2);
nanmean(outp)
nanstd(outp)
min(outp)
max(outp)
