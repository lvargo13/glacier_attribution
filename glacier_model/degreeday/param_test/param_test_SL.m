% after running run_param_test
% to find the optimum SNOWLINE parameters for all years
% for snowlines, for all glaciers except brewster

% input data:
gla = 'salisbury';  % glacier name
r1= 1700; % min and max SL range for plotting scatter plot  
r2 = 1900;
th=106;  % this should be the minimum RMSE
ann_thresh = th+(th/2); % threshold for annual RMSE for parameters used
mean_thresh = th/2; % threshold for mean RMSE for parameters used
figsave = 0; % 1 to save figs and write text output, 0 not to (generally check first then save)

%%%%%%%%%%%%%%%%

savedat = [gla '_psuite28feb_SL.mat'] ; %file to save for input to run att
savefil = [gla '_28feb_SL']; % names of each file 

fils_path = ['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/' gla 'SL/']; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

ddf = 0.5:0.4:1.7;
radf = 0.13:0.04:0.21 ;
ta = -2:0.4:-0.6;
pa = 0.8:0.4:1.8 ;
fi = length(ddf)*length(radf)*length(ta)*length(pa);

load([gla '_SL.mat']); % measured snowline data
em = nanmean(SL);
na = find(isnan(SL)>0);

tf_all = [];
rf_all = [];
ta_all = [];
pa_all = [];
ela_all = []; %zeros(length(fils),6); % 6 years of meausred data
ela_rms = []; %zeros(length(fils),1);
ela_m = []; %zeros(length(fils),1);
compot = [];
for i = 1:length(fils)
   load([fils_path fils(i).name])
   dt = find((abs(ddf-CONFIG.DegreeDay.DDF)<0.001)) ;
   dr = find((abs(radf-CONFIG.DegreeDay.RadiationFactor)<0.001)) ;    
   dta = find((abs(ta-CONFIG.DegreeDay.TempOffset)<0.001)) ;
   dpa = find((abs(pa-CONFIG.DegreeDay.PptnFactor)<0.001)) ;
   if (~isempty(dt) && ~isempty(dr) && ~isempty(dta) && ~isempty(dpa))
       %load([fils_path fils(i).name]) 
       ela_all = [ela_all ela]; 
       ela_rms = [ela_rms sqrt(nanmean((SL - ela').^2))];
       ela(na) = [];  % beacuse na ind are nan, dont want to comp those
       ela_m =  [ela_m sqrt((nanmean(ela) - em) .^2)];
       tf_all = [tf_all CONFIG.DegreeDay.DDF];
       rf_all = [rf_all CONFIG.DegreeDay.RadiationFactor];
       ta_all = [ta_all CONFIG.DegreeDay.TempOffset];
       pa_all = [pa_all CONFIG.DegreeDay.PptnFactor];
   else
       compot = [compot i];
   end
end

% for writing text output
lf_mb = find(ela_rms<ann_thresh);
lf_mmb = find(ela_m<mean_thresh);
D = intersect(lf_mb,lf_mmb);
tf = tf_all(D); 
rf = rf_all(D); 
tad = ta_all(D); 
pad = pa_all(D);
r = ela_rms(D); 
rm = ela_m(D); 
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

% to get one best run, for leave-one-year-out test
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
ela_m(D(C));
ela_rms(D(C));

n = zeros(length(D),1);
for i = 1:length(D)
   n(i) = (std(ela_all(:,D(i))))/sqrt(13);
end
safil =  ['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/stderr_' gla '_28feb_SL.mat'];
save(safil,'n'); 

figure; plot(1981:2017, SL, 'bo--');hold on
plot(1981:2017,ela_all(:,D),'ko--')
xlabel('year'); ylabel('snowline elevation (m a.s.l.)')
xlim([1980 2018])
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/' gla '_SL_caliball.pdf'])
end

figure; plot(SL,ela_all(:,D(:)),'.k'); hold on
plot([r1, r2], [r1, r2],'--k'); 
xlabel('measured snowline'); ylabel('modeled snowline')
xlim([r1 r2]); ylim([r1 r2])
if figsave == 1
    saveas(gcf,['/Volumes/arc_03/vargola/' gla '_SL_scatter.pdf'])
end

ed = ela_all(:,D); 
mx = max(ed'); 
mn = min(ed'); 
et = (SL<=mx & SL>=mn); 
sum(et)/(37-sum(isnan(SL)))

figure
subplot(2,2,1); plot(pad,'o'); ylim([0.8 1.8]); ylabel('precip adj')
subplot(2,2,2); plot(tad,'o'); ylim([-2 -0.6]); ylabel('temp adj')
subplot(2,2,3); plot(tf,'o'); ylim([0.5 1.7]); ylabel('temp fact')
subplot(2,2,4); plot(rf,'o'); ylim([0.13 0.21]); ylabel('rad fact')

figure; plot(ela_rms,'o')
figure; plot(ela_m,'o')

% whos D
% min_val = [min(r), max(r); min(rm), max(rm)];
% disp(min_val)

