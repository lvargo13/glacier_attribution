
% after running run_param_test
% for pdd with TF only and no SWR

fils_path = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test/brewster_2018_noSWR/'; 
fils = dir([fils_path '*.mat']);
fils = fils(~[fils.isdir]);
[~,idx] = sort([fils.datenum]);
fils = fils(idx);

measmb = [1376 691 692 -1698 -702 -74 -1728 -565 201 470 215 -1193 553];
%mberr = [247 283 403 228 211 241 187 281 418 373 314 305 215]; 

ela_rms = zeros(length(fils),1);
mmb_rms = zeros(length(fils),1);
tf_all = zeros(length(fils),1);
%rf_all = zeros(length(fils),1);
for i = 1:length(fils)
   load([fils_path fils(i).name]) 
   ela_rms(i) = ela_rmse;
   mmb_rms(i) = rmse;
   tf_all(i) = CONFIG.DegreeDay.DDF;
   %rf_all(i) = CONFIG.DegreeDay.RadiationFactor;
end

yyaxis left
plot(tf_all,mmb_rms,'--o'); ylabel('annual mass balance RMSE (mm)')
yyaxis right
plot(tf_all,ela_rms,'--o'); ylabel('annual snowline RMSE (m)')
xlabel('Temperature factor'); 

