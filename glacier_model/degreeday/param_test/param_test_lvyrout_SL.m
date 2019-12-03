%close all
% after running run_param_test
% for leave-one-year-out calibration

yrout = 3; 

%fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test_postreview/brew_highrange/'; 
fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/southcameronSL/'; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

%load('/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/compare_clim_sl/brew_sldata')
load southcameron_SL.mat
%meassl = sl_mean(end-12:end);
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
       lf_mb = find(rmserr<l_mb+3); 
       lf_mmb = find(mmb_rms<l_mmb+3);
       C = intersect(lf_mb,lf_mmb);
   end
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+6); 
       lf_mmb = find(mmb_rms<l_mmb+6);
       C = intersect(lf_mb,lf_mmb);
   end   
end


rmserr(C)
mmb_rms(C)
% tf_all(C)
% rf_all(C)
% sv = strcat('calib_yrout_',num2str(yrout),'.mat'); 
% save([fils_path sv],'rmserr','mmb_rms','tf_all','rf_all','C');

for ii = 1:length(C)
 sqrt((mw_m-mw(C(ii))).^2)
end
