% after running run_param_test
% to find the optimum mass balance parameters for all years
% just for rolleston mass balance

% input data:
r1= -2200; % min and max SL range for plotting scatter plot  
r2 = 800;
ann_thresh = 600; % threshold for annual RMSE for parameters used
mean_thresh = 60; % threshold for mean RMSE for parameters used
figsave = 1; % 1 to save figs and write text output, 0 not to (generally check first then save)


%%%%%%%%%%%%%%%%

savedat = 'rolleston_psuite.mat' ; %file to save for input to run att
savefil = 'rolleston_p'; % names of each file 

fils_dir = 'rolleston4full/';
fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/'; 
fils = dir([fils_path fils_dir '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

measmb = [-2039 -429 740 -38 657 -1006]; 
m = mean(measmb);

mmb_all = zeros(length(fils),6); % 6 years of meausred data
mmb_rms = zeros(length(fils),1);
rmserr = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
ta_all = zeros(length(fils),1);
pa_all = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils_dir fils(i).name]) 
   mmb_all(i,:) = mmb_dateadj;
   mmb_rms(i) = sqrt( ((mean(mmb_dateadj) -m).^2));
   rmserr(i) = sqrt(mean((mmb_dateadj' - measmb).^2));
   tf_all(i) = CONFIG.DegreeDay.DDF;
   rf_all(i) = CONFIG.DegreeDay.RadiationFactor;
   ta_all(i) = CONFIG.DegreeDay.TempOffset;
   pa_all(i) = CONFIG.DegreeDay.PptnFactor;
end

% % Get an idea of best parameters
% l_mb = min(rmserr);
% lf_mb = find(rmserr==l_mb);
% l_mmb = min(mmb_rms); 
% lf_mmb = find(mmb_rms==l_mmb);
% 
% if lf_mb == lf_mmb
%     'same'
% else
%    lf_mb = find(rmserr<l_mb+5); 
%    lf_mmb = find(mmb_rms<l_mmb+11);
%    C = intersect(lf_mb,lf_mmb);
% end
% tf_all(C)
% rf_all(C)

% for writing text output
lf_mb = find(rmserr<ann_thresh);
lf_mmb = find(mmb_rms<mean_thresh);
D = intersect(lf_mb,lf_mmb);
tf = tf_all(D); 
rf = rf_all(D); 
tad = ta_all(D); 
pad = pa_all(D);
r = rmserr(D); 
rm = mmb_rms(D); 
glacs_out = cell(length(D),1);
for i=1:length(D)
    tmp = [savefil,num2str(i),'/'];
    glacs_out(i) = {tmp};
end
if figsave == 1 
    save(savedat,'tf','rf','tad','pad','glacs_out');
end

figure; plot(2011:2016,measmb,'o--'); hold on 
plot(2011:2016,mmb_all(D,:),'ko--')
xlabel('year'); ylabel('mass balance (mm w.e.)')
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/rolleston_MB_caliball.pdf'])
end

figure; plot(measmb,mmb_all(D(:),:),'.k'); hold on
plot([r1, r2], [r1, r2],'--k'); 
xlabel('measured mass balance'); ylabel('modeled mass balance')
xlim([r1 r2]); ylim([r1 r2])
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/rolleston_MB_scatter.pdf')
end

figure
subplot(2,2,1); plot(pad,'o'); ylim([0.8 1.8]); ylabel('precip adj')
subplot(2,2,2); plot(tad,'o'); ylim([-2 -0.6]); ylabel('temp adj')
subplot(2,2,3); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
subplot(2,2,4); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')


