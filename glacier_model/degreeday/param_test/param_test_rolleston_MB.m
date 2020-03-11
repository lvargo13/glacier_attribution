% after running run_param_test
% to find the optimum mass balance parameters for all years
% just for rolleston mass balance

% input data:
r1= -2300; % min and max SL range for plotting scatter plot  
r2 = 800;
ann_thresh = 580+(580/2); % threshold for annual RMSE for parameters used
mean_thresh = (580/2); % threshold for mean RMSE for parameters used
figsave = 0; % 1 to save figs and write text output, 0 not to (generally check first then save)


%%%%%%%%%%%%%%%%

savedat = 'rolleston_psuite27feb.mat' ; %file to save for input to run att
savefil = 'rolleston_27feb'; % names of each file 

fils_dir = 'rolleston4full/';
fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/'; 
fils = dir([fils_path fils_dir '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

measmb = [-2039 -429 740 -38 657 -1006]; 
m = mean(measmb);

ddf = 0.5:0.3:1.7;
radf = 0.13:0.03:0.25 ;
ta = -2:0.45:-0.5;
pa = 0.8:0.45:1.8 ;
fi = length(ddf)*length(radf)*length(ta)*length(pa); 

tf_all = [];
rf_all = [];
ta_all = [];
pa_all = [];
mmb_all = []; %zeros(length(fils),6); % 6 years of meausred data
mmb_rms = []; %zeros(length(fils),1);
rmserr = []; %zeros(length(fils),1);

for i = 1:length(fils)
   load([fils_path fils_dir fils(i).name])
   dt = find(ddf==CONFIG.DegreeDay.DDF) ;
   dr = find(radf==CONFIG.DegreeDay.RadiationFactor) ;    
   dta = find(ta==CONFIG.DegreeDay.TempOffset) ;
   dpa = find(pa==CONFIG.DegreeDay.PptnFactor) ;
   if (~isempty(dt) && ~isempty(dr) && ~isempty(dta) && ~isempty(dpa)) 
       load([fils_path fils_dir fils(i).name]) 
       mmb_all = [mmb_all mmb_dateadj];
       mmb_rms = [mmb_rms sqrt( ((mean(mmb_dateadj) -m).^2))];
       rmserr = [rmserr sqrt(mean((mmb_dateadj' - measmb).^2))];
       tf_all = [tf_all CONFIG.DegreeDay.DDF];
       rf_all = [rf_all CONFIG.DegreeDay.RadiationFactor];
       ta_all = [ta_all CONFIG.DegreeDay.TempOffset];
       pa_all = [pa_all CONFIG.DegreeDay.PptnFactor];
   end
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
param_usd = [min(tf), max(tf); min(rf), max(rf); min(tad), max(tad); min(pad), max(pad)];
disp(param_usd)

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

n = zeros(length(D),1);
for i = 1:length(D)
    %n(i) = std(mmb_all(D(i),:));
   n(i) = (std(mmb_all(:,D(i))))/sqrt(6);
end
if figsave == 1
    save(['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/stderr_' savefil '.mat'],'n'); 
end

mberr = zeros(6,1)+300;
figure; 
errorbar(2011:2016, measmb, mberr,'.','LineWidth',2.5) ; hold on
plot(2011:2016,measmb,'bo--'); hold on 
plot(2011:2016,mmb_all(:,D),'ko--')
xlabel('year'); ylabel('mass balance (mm w.e.)')
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/rolleston_MB_caliball.pdf'])
end

figure;
plot(measmb,mmb_all(:,D(:)),'.k'); hold on
plot([r1, r2], [r1, r2],'--k','LineWidth',2); 
xlabel('measured mass balance'); ylabel('modeled mass balance')
xlim([r1 r2]); ylim([r1 r2])
plot([-2000, 800], [-2300, 500],'--k'); 
plot([-2300, 500], [-2000, 800],'--k'); 
%plot([-2300, 800], [-2600, 800],'--k','LineWidth',2);
if figsave == 1
    saveas(gcf,'/Volumes/arc_03/vargola/rolleston_MB_scatter.pdf')
end

figure
subplot(2,2,1); plot(pad,'o'); ylim([0.8 1.8]); ylabel('precip adj')
subplot(2,2,2); plot(tad,'o'); ylim([-2 -0.6]); ylabel('temp adj')
subplot(2,2,3); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
subplot(2,2,4); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')

whos D
min_val = [min(r), max(r); min(rm), max(rm)];
disp(min_val)
