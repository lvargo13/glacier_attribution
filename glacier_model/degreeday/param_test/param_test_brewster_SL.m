% after running run_param_test
% to find the optimum SNOWLINE parameters for all years for brewster only


% input data:
r1= 1750; % min and max SL range for plotting scatter plot  
r2 = 2330;
ann_thresh = 78+(78/2); % threshold for annual RMSE for parameters used
mean_thresh = (78/2); % threshold for mean RMSE for parameters used
figsave = 0; % 1 to save figs and write text output, 0 not to (generally check first then save)


%%%%%%%%%%%%%%%%%
savedat = 'brewster_27feb_SL.mat' ; %file to save for input to run att
savefil = 'brewster_27feb_SL'; % names of each file 

fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/brew_highrange/'; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

load('/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/compare_clim_sl/brew_sldata')

em = mean(sl_mean(end-12:end)); 

ddf = 0.5:0.08:1.7;
radf = 0.13:0.02:0.25 ;
fi = length(ddf)*length(radf); 

% ela_m = zeros(fi,1);
% ela_all = zeros(fi,13);
% ela_rms = zeros(fi,1);
% tf_all = zeros(fi,1);
% rf_all = zeros(fi,1);
ela_m = [];
ela_all = [];
ela_rms = [];
tf_all = [];
rf_all = [];
for i = 1:length(fils)
   load([fils_path fils(i).name])  
   dt = find((abs(ddf-CONFIG.DegreeDay.DDF)<0.001)) ;
   dr = find((abs(radf-CONFIG.DegreeDay.RadiationFactor)<0.001)) ; 
   if (~isempty(dt) && ~isempty(dr))
       tf_all = [tf_all CONFIG.DegreeDay.DDF];
       rf_all = [rf_all CONFIG.DegreeDay.RadiationFactor]; 
       ela_m =  [ela_m sqrt((nanmean(ela(end-12:end)) - em) .^2)];
       ela_all = [ela_all ela(end-12:end)]; 
       ela_rms = [ela_rms sqrt(nanmean((sl_mean(end-12:end) - ela(end-12:end)').^2))];
   end
end

%e_rmse = reshape(ela_rms,length(radf),length(ddf));
% figure; 
% subplot(2,2,1); imagesc(radf,ddf,e_rmse')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('annual ela rmse (m)')
% 
% subplot(2,2,3); imagesc(radf,ddf,e_rmse')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('annual ela rmse (m)'); caxis([75 85])
% 
% em_rms = reshape(ela_m,length(radf),length(ddf));
% subplot(2,2,2); imagesc(radf,ddf,em_rms')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('mean ela rmse (m)')
% 
% subplot(2,2,4); imagesc(radf,ddf,em_rms')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('mean ela rmse (m)'); caxis([0 18])

% l_mb = min(ela_rms);
% lf_mb = find(ela_rms==l_mb);
% l_mmb = min(ela_m); 
% lf_mmb = find(ela_m==l_mmb);


% for writing text output
lf_mb = find(ela_rms<ann_thresh);
lf_mmb = find(ela_m<mean_thresh);
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

n = zeros(length(D),1);
for i = 1:length(D)
   n(i) = (std(ela_all(:,D(i))))/sqrt(13);
end
save stderr_brewster_27feb_SL n 

l_mb = min(r);
lf_mb = find(r==l_mb);
l_mmb = min(rm); 
lf_mmb = find(rm==l_mmb);
C = intersect(lf_mb,lf_mmb);
stps = 0.1; % can change 
% need an 'if c is empty quit', here
while length(C) < 1
   lf_mb = find(r<l_mb+stps);
   lf_mmb = find(rm<l_mmb+stps);
   stps = stps+0.1;
   C = intersect(lf_mb,lf_mmb);
end
D(C)

figure; plot(2005:2017, sl_mean(end-12:end), 'bo--');hold on
plot(2005:2017,ela_all(:,D),'ko--')
xlabel('year'); ylabel('snowline elevation (m a.s.l.)')
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/brewster_SL_caliball.pdf')
end

figure; plot(sl_mean(end-12:end),ela_all(:,D(:)),'.k'); hold on
plot([r1, r2], [r1, r2],'--k'); 
xlabel('measured snowline'); ylabel('modeled snowline')
xlim([r1 r2]); ylim([r1 r2])
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/brewster_SL_scatter.pdf') 
end

%ed = ela_all(:,D);

% figure; subplot(1,2,1); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
% subplot(1,2,2); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')


