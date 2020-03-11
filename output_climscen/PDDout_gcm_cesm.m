% read in output files of sa_pdd_gcm_loop_param
% similar to PDDout_loop_hist, but with the goal of getting individual
% probabilities for GCM and CESM for the Supplemental table

glacrun = {'rolleston','salisbury','parkpass','southcameron','thurneyson','vertebrae12','ridge'};

for y =1:length(glacrun)
glac = [glacrun{y} '_28feb_SL*/'];

gcmflg = 0; % flag, 1 for GCM, 0 for CESM
sl=1;  % flag, 0 for mb, 1 for sl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch glac
    case 'brewster_27feb_SL*/'
        glac_11_sl = 2303 ; 
        glac_18_sl = 2311 ;  
        glac_11_mb = -1728; 
        glac_18_mb = -2217;
        glac_18_mbsig = 323;
        glac_11_mbsig = 187;
        sl_bins = 1600:50:2450;
    case 'rolleston_27feb*/'
        glac_11_sl = 1857 ;
        glac_18_sl = 1847 ; 
        glac_11_mb = -2039-300;
        glac_18_mb = -2456-300;
        sl_bins = 1650:25:2000;
    case 'parkpass_28feb_SL*/'
        glac_18_sl = 1966 ; 
        glac_11_sl = 2040 ; 
        sl_bins = 1350:50:2450;
    case 'southcameron_28feb_SL*/'
        glac_18_sl = 2426 ; 
        glac_11_sl = 2395 ; 
        sl_bins = 1900:35:2800; 
    case 'salisbury_28feb_SL*/'
        glac_18_sl = 1967 ; 
        glac_11_sl = 1967 ; % no 2011 SL, so just to run
        sl_bins = 1300:50:2500;
    case 'thurneyson_28feb_SL*/'
        glac_18_sl = 2232 ; 
        glac_11_sl = 2125 ; 
        sl_bins = 1650:40:2450;
    case 'ridge_28feb_SL*/'
        glac_11_sl = 2436 ;
        glac_18_sl = 2434 ;
        sl_bins = 1950:40:2700;
    case 'vertebrae12_28feb_SL*/'
        glac_18_sl = 2086 ; 
        glac_11_sl = 2094 ; 
        sl_bins = 1600:40:2300;
    case 'glenmary_28feb_SL*/' 
        glac_11_sl = 2297 ;
        glac_18_sl = 2327 ; 
        sl_bins = 1900:30:2600; 
    case 'vertebrae25_28feb_SL*/'
        glac_18_sl = 2001.8 ; % because SL is over the top
        glac_11_sl = 2001.8 ; 
        sl_bins = 1650:30:2150;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---- read in data
a = 2;  % scenarios
mb_bins = -5500:300:4100;

in_dir_g = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/gcm/';
fi_g = dir([in_dir_g glac]);
[~,idx] = sort([fi_g.datenum]); fi_g = fi_g(idx); % sort for using std err

in_dir_c = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/cesm/';
fi_c = dir([in_dir_c glac]);
[~,idx] = sort([fi_c.datenum]); fi_c = fi_c(idx); % sort for using std err
glacs_out = zeros(length(fi_c),1); %how many param. combos were run

mbp = zeros(length(mb_bins),a,length(glacs_out)); % bins x 4 scenarios x number of pruns
slp = zeros(length(sl_bins),a,length(glacs_out)); % bins x 4 scenarios x number of pruns

mb_prob = zeros(length(mb_bins),a); 
sl_prob = zeros(length(sl_bins),a); 
load(['stderr/stderr_' glac(1:end-2) '.mat']);  % loading the std err for each param.

hn = 100; % number of points in each histogram
if gcmflg == 1
    pastallhs = zeros(1664*hn,length(glacs_out));
    presallhs = zeros(320*hn,length(glacs_out));
else
    pastallhs = zeros(1799*hn,length(glacs_out));
    presallhs = zeros(680*hn,length(glacs_out));
end

for ii = 1:length(glacs_out)  
 
% read in gcm data
if gcmflg == 1
    files = dir([in_dir_g fi_g(ii).name,'/*.mat']);
    [g_pres,g_past,g_pres_sl,g_past_sl] = read_pdd_gcm([in_dir_g fi_g(ii).name, '/'],files);
else
    % read in cesm data
    files_c = [in_dir_c fi_c(ii).name]; 
    cfiles = dir([files_c, '/', '*present.mat']);  % dont want NAT run
    [g_pres,g_pres_sl] = read_pdd_cesm([in_dir_c fi_c(ii).name, '/'],cfiles); 
    cfilesNAT = dir([files_c, '/', '*NAT.mat']);
    load([files_c,'/', cfilesNAT.name])
    g_past = mb; 
    g_past_sl = ela; 
end

if sl == 0
    presall = g_pres;
    pastall = g_past;
else
    presall = g_pres_sl;
    pastall = g_past_sl; 
end

% implement incorporating inherant model error
presallh = zeros(length(presall), hn); 
for k = 1:length(presall)
   presallh(k,:) = normrnd(presall(k),n(ii),[1,hn]);  
end
presall_h = presallh(:); 

pastallh = zeros(length(pastall), hn);
for k = 1:length(pastall)
   pastallh(k,:) = normrnd(pastall(k),n(ii),[1,hn]);  
end
pastall_h = pastallh(:);
pastallhs(:,ii) = pastall_h;
presallhs(:,ii) = presall_h;
climscen = {presall_h, pastall_h}; 

% probability
if sl == 0
    for cs = 1:a  % mass balance probability
      cs_var = climscen{cs};  
        for i = 1:length(cs_var)
              [~, ind] = min(abs(mb_bins-cs_var(i)));  % get ind of closest bin
              mb_prob(ind,cs) = mb_prob(ind,cs)+1;
        end
       mb_prob(:,cs) = mb_prob(:,cs) / length(cs_var); 
    end
else
    for cl = 1:a   
      cl_var = climscen{cl};
      cl_var = cl_var(~isnan(cl_var));  % get rid of nans
        for i = 1:length(cl_var)
              [~, indl] = min(abs(sl_bins-cl_var(i)));  % get ind of closest bin
              sl_prob(indl,cl) = sl_prob(indl,cl)+1;
        end
       sl_prob(:,cl) = sl_prob(:,cl) / length(cl_var); 
    end
end

if sl == 0
    mbp(:,:,ii) = mb_prob;
else
    slp(:,:,ii) = sl_prob; 
end
end

% percent chances
pres_slp = zeros(length(glacs_out),1);
past_slp = zeros(length(glacs_out),1);

if sl == 0  % if MB
    past_11 = length(find(pastallhs(:)<=glac_11_mb)) ./ length(pastallhs(:));
    pres_11 = length(find(presallhs(:)<=glac_11_mb)) ./ length(presallhs(:));
    past_18 = length(find(pastallhs(:)<=glac_18_mb)) ./ length(pastallhs(:));
    pres_18 = length(find(presallhs(:)<=glac_18_mb)) ./ length(presallhs(:));
else
    past_11 = length(find(pastallhs(:)>=glac_11_sl)) ./ length(pastallhs(:));
    pres_11 = length(find(presallhs(:)>=glac_11_sl)) ./ length(presallhs(:));
    past_18 = length(find(pastallhs(:)>=glac_18_sl)) ./ length(pastallhs(:));
    pres_18 = length(find(presallhs(:)>=glac_18_sl)) ./ length(presallhs(:));
end
meanprob = [past_11*100 pres_11*100 past_18*100 pres_18*100 pres_11/past_11 pres_18/past_18];

if gcmflg == 1
    fsv = ['gcm_att_out_' glac(1:end-2)]; 
else
    fsv = ['cesm_att_out_' glac(1:end-2)]; 
end
save(fsv,'meanprob','glac'); 

 clearvars -except glacrun y
end
