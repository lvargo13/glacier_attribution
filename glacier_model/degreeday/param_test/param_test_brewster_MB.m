% after running run_param_test
% to find the optimum mass balance parameters for all years for brewster mass balance

% input data:
% r1= 1600; % min and max SL range for plotting scatter plot  
% r2 = 2300;
minrmse = 560
ann_thresh = minrmse+(minrmse/2); % threshold for annual RMSE for parameters used
mean_thresh = minrmse/2; % threshold for mean RMSE for parameters used
figsave = 0; % 1 to save figs and write text output, 0 not to (generally check first then save)

%%%%%%%%%%%%%%%%
savedat = 'brewster_psuite14feb.mat' ; %file to save for input to run att
savefil = 'brewster_14feb'; % names of each file 

fils_dir = 'brewster_14feb/';
fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/'; 
fils = dir([fils_path fils_dir '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

measmb = [1376 691 692 -1698 -702 -74 -1728 -565 201 470 215 -1193 553 ]; 
mberr = [247 283 403 228 211 241 187 281 418 373 314 305 215]; 
m = mean(measmb);

mmb_all = zeros(length(fils),13);
mmb_rms = zeros(length(fils),1);
rmserr = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
rf_all = zeros(length(fils),1);
std_all = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils_dir fils(i).name]) 
   mmb_all(i,:) = mmb_dateadj;
   mmb_rms(i) = sqrt( ((mean(mmb_dateadj) -m).^2));
   rmserr(i) = sqrt(mean((mmb_dateadj' - measmb).^2));
   tf_all(i) = CONFIG.DegreeDay.DDF;
   rf_all(i) = CONFIG.DegreeDay.RadiationFactor;
   std_all(i) = std(mmb_dateadj); 
end
% ddf = 0.5:0.04:1.7; %original
% radf = 0.13:0.01:0.25 ;

ddf = 0.5:0.08:1.7; %original
radf = 0.13:0.02:0.25 ;

rmserr_rs = reshape(rmserr,length(radf),length(ddf));
std_rs = reshape(std_all,length(radf),length(ddf));

figure; 
% subplot(2,2,1); imagesc(radf,ddf,rmserr_rs')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('annual mb rmse (mm w.e.)')

subplot(1,2,1); imagesc(radf,ddf,rmserr_rs')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('annual mb rmse (mm w.e.)'); caxis([560 ann_thresh])

mmb_rmse = reshape(mmb_rms,length(radf),length(ddf));
% subplot(2,2,2); imagesc(radf,ddf,mmb_rmse')
% xlabel('radiation factor')
% ylabel('temperature factor')
% colorbar; title('mean mb rmse (mm w.e.)')

subplot(1,2,2); imagesc(radf,ddf,mmb_rmse')
xlabel('radiation factor')
ylabel('temperature factor')
colorbar; title('mean mb rmse (mm w.e.)'); caxis([0 300])

% l_mb = min(rmserr);
% lf_mb = find(rmserr==l_mb);
% l_mmb = min(mmb_rms); 
% lf_mmb = find(mmb_rms==l_mmb);

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
r = mmb_rms(D);
rm = rmserr(D);
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
   n(i) = (std(mmb_all(D(i),:)))/sqrt(13);
end

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

figure; 
errorbar(2005:2017, measmb, mberr,'.','LineWidth',2.5) ; hold on
plot(2005:2017, measmb, 'bo--')
plot(2005:2017,mmb_all(D,:),'ko--')
xlabel('year'); ylabel('mass balance (mm w.e.)')
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/brewster_MB_caliball.pdf')
end

mt = mmb_all(D,:); 
ts = zeros(40, 13); 
for i = 1:length(D)
    ts(:,i) = mt(i,:)-measmb; 
end

vh = zeros(40,1);
for ii = 1:13
   figure; h = histogram(ts(:,ii),'BinEdges',-2000:100:2000);  
   vh(:,ii) = h.Values; 
end
ts = sum(vh,1); 
figure; hold on
for ii = 1:6
    stairs(vh(:,ii))
end
% figure, plot(n,rm,'o')
% xlabel('annual modeled standard deviation (mm w.e.)'); 
% ylabel('annual RMSE (measured & modeled) (mm w.e.)'); 

figure; plot(measmb,mmb_all(D(:),:),'.k'); hold on
plot([-1700, 1500], [-2000, 1200],'--k'); 
plot([-2000, 1200], [-1700, 1500],'--k'); 
plot([-2000, 1500], [-2000, 1500],'--k','LineWidth',2);

% plot([r1, r2], [r1, r2],'--k'); 
% xlabel('measured snowline'); ylabel('modeled snowline')
% xlim([r1 r2]); ylim([r1 r2])
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/brewster_MB_scatter.pdf')
end

