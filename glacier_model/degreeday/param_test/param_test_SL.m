% after running run_param_test
% to find the optimum SNOWLINE parameters for all years
% for snowlines, for all glaciers except brewster

% input data:
gla = 'chancellor';  % glacier name
r1= 1650; % min and max SL range for plotting scatter plot  
r2 = 1850;
ann_thresh = 56.5; % threshold for annual RMSE for parameters used
mean_thresh = 8; % threshold for mean RMSE for parameters used
figsave = 0; % 1 to save figs and write text output, 0 not to (generally check first then save)

%%%%%%%%%%%%%%%%

savedat = [gla '_psuiteSL.mat'] ; %file to save for input to run att
savefil = [gla '_psuiteSL']; % names of each file 

fils_path = ['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/' gla 'SL/']; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

load([gla '_SL.mat']); % measured snowline data
em = nanmean(SL);
na = find(isnan(SL)>0);


ela_m = zeros(length(fils),1);
ela_all = zeros(length(fils),37);
ela_rms = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
ta_all = zeros(length(fils),1);
pa_all = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils(i).name]) 
   ela_all(i,:) = ela; 
   ela_rms(i) = sqrt(nanmean((SL - ela').^2));
   ela(na) = [];  % beacuse na ind are nan, dont want to comp those
   ela_m(i) =  sqrt((nanmean(ela) - em) .^2);
   tf_all(i) = CONFIG.DegreeDay.DDF;
   rf_all(i) = CONFIG.DegreeDay.RadiationFactor;
   ta_all(i) = CONFIG.DegreeDay.TempOffset;
   pa_all(i) = CONFIG.DegreeDay.PptnFactor;
end

% for writing text output
lf_mb = find(ela_rms<ann_thresh);
lf_mmb = find(ela_m<mean_thresh);
D = intersect(lf_mb,lf_mmb);
% tf = tf_all(D); 
% rf = rf_all(D); 
% tad = ta_all(D); 
% pad = pa_all(D);
% r = ela_rms(D); 
% rm = ela_m(D); 
% glacs_out = cell(length(D),1);
% for i=1:length(D)
%     tmp = [savefil,num2str(i),'/'];
%     glacs_out(i) = {tmp};
% end
% if figsave == 1 
%     save(savedat,'tf','rf','tad','pad','glacs_out');
% end
% 
% % to get one best run, for leave-one-year-out test
% l_mb = min(r);
% lf_mb = find(r==l_mb);
% l_mmb = min(rm); 
% lf_mmb = find(rm==l_mmb);
% C = intersect(lf_mb,lf_mmb);
% stps = 0.1; % can change 
% % need an 'if c is empty quit', here
% while length(C) < 1
%    lf_mb = find(r<l_mb+stps);
%    lf_mmb = find(rm<l_mmb+stps);
%    stps = stps+0.1;
%    C = intersect(lf_mb,lf_mmb);
% end
% D(C)
% ela_m(D(C))
% ela_rms(D(C))
% 
% figure; plot(1981:2017, SL, 'bo--');hold on
% plot(1981:2017,ela_all(D,:),'ko--')
% xlabel('year'); ylabel('snowline elevation (m a.s.l.)')
% xlim([1980 2018])
% if figsave == 1
%     saveas(gcf,['/Volumes/arc_03/vargola/' gla '_SL_caliball.pdf'])
% end
% 
% figure; plot(SL,ela_all(D(:),:),'.k'); hold on
% plot([r1, r2], [r1, r2],'--k'); 
% xlabel('measured snowline'); ylabel('modeled snowline')
% xlim([r1 r2]); ylim([r1 r2])
% if figsave == 1
%     saveas(gcf,['/Volumes/arc_03/vargola/' gla '_SL_scatter.pdf'])
% end
% 
% figure
% subplot(2,2,1); plot(pad,'o'); ylim([0.8 1.8]); ylabel('precip adj')
% subplot(2,2,2); plot(tad,'o'); ylim([-2 -0.6]); ylabel('temp adj')
% subplot(2,2,3); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
% subplot(2,2,4); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')

% figure; plot(ela_rms,'o')
% figure; plot(ela_m,'o')
