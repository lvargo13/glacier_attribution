close all
% after running run_param_test
% for leave-one-year-out calibration

yrout = 6; 

%fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test_postreview/brew_highrange/'; 
fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/rolleston4full/'; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

%measmb = [1376 691 692 -1698 -702 -74 -1728 -565 201 470 215 -1193 553];
measmb = [-2039 -429 740 -38 657 -1006];
mw_m = measmb(yrout); 
measmb(yrout) = [];  % delete year out
m = mean(measmb);

mmb_all = zeros(length(fils),1);
mmb_rms = zeros(length(fils),1);
rmserr = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
mw = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils(i).name]) 
   mw(i) = mmb_dateadj(yrout);
   mmb_dateadj(yrout)=[]; 
   mmb_all(i) = mean(mmb_dateadj);
   mmb_rms(i) = sqrt( ((mean(mmb_dateadj) -m).^2));
   rmserr(i) = sqrt(mean((mmb_dateadj' - measmb).^2));
end

l_mb = min(rmserr);
lf_mb = find(rmserr==l_mb);
l_mmb = min(mmb_rms); 
lf_mmb = find(mmb_rms==l_mmb);

if lf_mb == lf_mmb   % brewster
    'same'
else
   lf_mb = find(rmserr<l_mb+3); 
   lf_mmb = find(mmb_rms<l_mmb+40);
   C = intersect(lf_mb,lf_mmb);
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+6); 
       lf_mmb = find(mmb_rms<l_mmb+50);
       C = intersect(lf_mb,lf_mmb);
   end
   if isempty(C) == 1 ; % is empty
       lf_mb = find(rmserr<l_mb+10); 
       lf_mmb = find(mmb_rms<l_mmb+60);
       C = intersect(lf_mb,lf_mmb);
   end   
end

% if lf_mb == lf_mmb   % brewster
%     'same'
% else
%    lf_mb = find(rmserr<l_mb+1); 
%    lf_mmb = find(mmb_rms<l_mmb+11);
%    C = intersect(lf_mb,lf_mmb);
%    if isempty(C) == 1 ; % is empty
%        lf_mb = find(rmserr<l_mb+3); 
%        lf_mmb = find(mmb_rms<l_mmb+15);
%        C = intersect(lf_mb,lf_mmb);
%    end
%    if isempty(C) == 1 ; % is empty
%        lf_mb = find(rmserr<l_mb+10); 
%        lf_mmb = find(mmb_rms<l_mmb+15);
%        C = intersect(lf_mb,lf_mmb);
%    end   
% end


rmserr(C)
mmb_rms(C)
% tf_all(C)
% rf_all(C)
% sv = strcat('calib_yrout_',num2str(yrout),'.mat'); 
% save([fils_path sv],'rmserr','mmb_rms','tf_all','rf_all','C');

for ii = 1:length(C)
 sqrt((mw_m-mw(C(ii))).^2)
end
