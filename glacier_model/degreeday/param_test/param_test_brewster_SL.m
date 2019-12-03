% after running run_param_test
% to find the optimum SNOWLINE parameters for all years for brewster only


% input data:
r1= 1750; % min and max SL range for plotting scatter plot  
r2 = 2330;
ann_thresh = 95; % threshold for annual RMSE for parameters used
mean_thresh = 18; % threshold for mean RMSE for parameters used
figsave = 1; % 1 to save figs and write text output, 0 not to (generally check first then save)


%%%%%%%%%%%%%%%%%
savedat = 'brewster_psuiteSL40.mat' ; %file to save for input to run att
savefil = 'brewster_psuiteSL40'; % names of each file 

fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/brew_highrange/'; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

load('/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/compare_clim_sl/brew_sldata')

em = mean(sl_mean(end-12:end)); 

ela_m = zeros(length(fils),1);
ela_all = zeros(length(fils),13);
ela_rms = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils(i).name]) 
   ela_m(i) =  sqrt((nanmean(ela(end-12:end)) - em) .^2);
   ela_all(i,:) = ela(end-12:end); 
   ela_rms(i) = sqrt(nanmean((sl_mean(end-12:end) - ela(end-12:end)').^2));
   tf_all(i) = CONFIG.DegreeDay.DDF;
   rf_all(i) = CONFIG.DegreeDay.RadiationFactor;
end

ddf = 0.5:0.04:1.7; %original
radf = 0.13:0.01:0.25 ;
% ddf = 0.5:0.06:1.7; %original
% radf = 0.15:0.015:0.25 ;

e_rmse = reshape(ela_rms,length(radf),length(ddf));
figure; 
subplot(2,2,1); imagesc(radf,ddf,e_rmse')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('annual ela rmse (m)')

subplot(2,2,3); imagesc(radf,ddf,e_rmse')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('annual ela rmse (m)'); caxis([75 85])

em_rms = reshape(ela_m,length(radf),length(ddf));
subplot(2,2,2); imagesc(radf,ddf,em_rms')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('mean ela rmse (m)')

subplot(2,2,4); imagesc(radf,ddf,em_rms')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('mean ela rmse (m)'); caxis([0 18])

% l_mb = min(ela_rms);
% lf_mb = find(ela_rms==l_mb);
% l_mmb = min(ela_m); 
% lf_mmb = find(ela_m==l_mmb);


% for writing text output
lf_mb = find(e_rmse<ann_thresh);
lf_mmb = find(em_rms<mean_thresh);
D = intersect(lf_mb,lf_mmb);
tf = tf_all(D); 
rf = rf_all(D); 
r = ela_rms(D); 
rm = ela_m(D); 
glacs_out = cell(length(D),1);
for i=1:length(D)
    tmp = [savefil,num2str(i),'/'];
    glacs_out(i) = {tmp};
end
if figsave == 1 
    save(savedat,'tf','rf','glacs_out');
end

figure; plot(2005:2017, sl_mean(end-12:end), 'bo--');hold on
plot(2005:2017,ela_all(D,:),'ko--')
xlabel('year'); ylabel('snowline elevation (m a.s.l.)')
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/brewster_SL_caliball.pdf'])
end

figure; plot(sl_mean(end-12:end),ela_all(D(:),:),'.k'); hold on
plot([r1, r2], [r1, r2],'--k'); 
xlabel('measured snowline'); ylabel('modeled snowline')
xlim([r1 r2]); ylim([r1 r2])
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/brewster_SL_scatter.pdf']) 
end

% figure; subplot(1,2,1); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
% subplot(1,2,2); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')


