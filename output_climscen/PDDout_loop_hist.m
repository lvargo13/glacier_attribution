% read in output files of sa_pdd_gcm_loop_param
% calculate and plot probability
% outputs: plots of natural and present probabilities 
% can also see probabilities without bootstrapping here

% added step to include model uncertainty as of march 2020

glacrun = {'rolleston','salisbury','parkpass','southcameron','thurneyson','vertebrae12','ridge'};

for y =1:length(glacrun)
glac = [glacrun{y} '_28feb_SL*/'];
sl=1;  % flag, 0 for mb, 1 for sl
svdat= 0; % flag to save fig and .mat if 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch glac
    case 'brewster_27feb_SL*/'
        glac_11_sl = 2303 ; 
        glac_18_sl = 2311 ;  
        glac_11_mb = -1728; 
        glac_18_mb = -2217;
        glac_18_mbsig = 323;
        glac_11_mbsig = 187;
        sl_bins = 1300:50:2650;
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
        sl_bins = 1700:30:2800; 
    case 'vertebrae25_28feb_SL*/'
        glac_18_sl = 2001.8 ; % because SL is over the top
        glac_11_sl = 2001.8 ; 
        sl_bins = 1600:30:2200;
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

% store for bootstrapping
hn = 100; % number of points in each histogram
pastallhs = zeros(3463*hn,length(glacs_out));
presallhs = zeros(1000*hn,length(glacs_out));

mb_prob = zeros(length(mb_bins),a); 
sl_prob = zeros(length(sl_bins),a); 
load(['stderr/stderr_' glac(1:end-2) '.mat']);  % loading the std err for each param.

for ii = 1:length(glacs_out)  
 
% read in gcm data
files = dir([in_dir_g fi_g(ii).name,'/*.mat']);
[g_pres,g_past,g_pres_sl,g_past_sl] = read_pdd_gcm([in_dir_g fi_g(ii).name, '/'],files);

% read in cesm data
files_c = [in_dir_c fi_c(ii).name]; 
cfiles = dir([files_c, '/', '*present.mat']);  % dont want NAT run
[c_pres,c_pres_sl] = read_pdd_cesm([in_dir_c fi_c(ii).name, '/'],cfiles); 
cfilesNAT = dir([files_c, '/', '*NAT.mat']);
load([files_c,'/', cfilesNAT.name])
c_past = mb; 
c_past_sl = ela; 

if sl == 0
    presall = [g_pres;c_pres];
    pastall = [g_past;c_past];
else
    presall = [g_pres_sl;c_pres_sl];
    pastall = [g_past_sl;c_past_sl]; 
end

% implement incorporating inherant model error
presallh = zeros(length(presall), hn); 
for k = 1:length(presall)
   presallh(k,:) = normrnd(presall(k),n(ii),[1,hn]);  
end
presall_h = presallh(:); 
%presall_h = normrnd(mean(presall_h),std(presall_h)*1.2,[1,length(presall_h)]);

pastallh = zeros(length(pastall), hn);
for k = 1:length(pastall)
   pastallh(k,:) = normrnd(pastall(k),n(ii),[1,hn]);  
end
pastall_h = pastallh(:);
%pastall_h = normrnd(mean(pastall_h),std(pastall_h)*1.2,[1,length(pastall_h)]);

pastallhs(:,ii) = pastall_h;
presallhs(:,ii) = presall_h;
climscen = {presall_h, pastall_h};  % this includes model error

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

%---- plot mass balance
if sl == 0
    bns = zeros(length(mb_bins),2, 3); % 2 past and pres
    for j = 1:length(mb_bins)
        bns(j,1,1) = min(mbp(j,1,:)); % past
        bns(j,1,2) = mean(mbp(j,1,:)); 
        bns(j,1,3) = max(mbp(j,1,:)); 
        bns(j,2,1) = min(mbp(j,2,:)); %pres
        bns(j,2,2) = mean(mbp(j,2,:)); 
        bns(j,2,3) = max(mbp(j,2,:)); 
    end
    
    mb_bins = mb_bins -(diff(mb_bins(1:2))/2);
    figure;
    stairs(mb_bins,bns(:,1,2),'r','LineWidth',2); hold on
    stairs(mb_bins,bns(:,2,2),'k','LineWidth',2)
    legend('present','natural')
    stairs(mb_bins,bns(:,1,1),'r'); hold on
    stairs(mb_bins,bns(:,1,3),'r')
    stairs(mb_bins,bns(:,2,1),'k')
    stairs(mb_bins,bns(:,2,3),'k')
    for ii = 1:length(mb_bins)-1
        x2 = [mb_bins(ii) mb_bins(ii+1) mb_bins(ii+1) mb_bins(ii)];
        fill_past= [bns(ii,1,1) bns(ii,1,1) bns(ii,1,3) bns(ii,1,3)];
        fill_pres= [bns(ii,2,1) bns(ii,2,1) bns(ii,2,3) bns(ii,2,3)];
        patch(x2, fill_past, 'r','FaceAlpha',.3,'LineStyle','none'); hold on
        fill(x2, fill_pres, 'k','FaceAlpha',.3,'LineStyle','none')
    end

    plot([glac_11_mb, glac_11_mb], [0, 0.25],'--k');
    plot([glac_18_mb, glac_18_mb], [0, 0.25],'--k');
    xlabel('Mass Balance (mm w.e.)')
    ylabel('Probability')
    xlim([-5500 4000])
    ylim([0 0.21])
end

%---- plot snowline probability
if sl == 1
    bs = zeros(length(sl_bins),2, 3); % 2 past and pres
    for j = 1:length(sl_bins)
        bs(j,1,1) = min(slp(j,1,:)); % past
        bs(j,1,2) = mean(slp(j,1,:)); 
        bs(j,1,3) = max(slp(j,1,:)); 
        bs(j,2,1) = min(slp(j,2,:)); %pres
        bs(j,2,2) = mean(slp(j,2,:));
        bs(j,2,3) = max(slp(j,2,:));
    end

    sl_bins = sl_bins -(diff(sl_bins(1:2))/2); 
    figure;
    stairs(sl_bins,bs(:,1,2),'r','LineWidth',2); hold on
    stairs(sl_bins,bs(:,2,2),'k','LineWidth',2)
    legend('natural','present')
    stairs(sl_bins,bs(:,1,1),'r'); hold on
    stairs(sl_bins,bs(:,1,3),'r')
    stairs(sl_bins,bs(:,2,1),'k')
    stairs(sl_bins,bs(:,2,3),'k')
    for ii = 1:length(sl_bins)-1
        x3 = [sl_bins(ii) sl_bins(ii+1) sl_bins(ii+1) sl_bins(ii)];
        fill_past= [bs(ii,1,1) bs(ii,1,1) bs(ii,1,3) bs(ii,1,3)];
        fill_pres= [bs(ii,2,1) bs(ii,2,1) bs(ii,2,3) bs(ii,2,3)];
        patch(x3, fill_past, 'r','FaceAlpha',.3,'LineStyle','none'); hold on
        fill(x3, fill_pres, 'k','FaceAlpha',.3,'LineStyle','none')
    end
    plot([glac_18_sl, glac_18_sl], [0,0.19],'--k');
    plot([glac_11_sl, glac_11_sl], [0,0.19],'--k');
    xlabel('snowline elevation (m)')
    ylabel('Probability')
    if svdat == 1
        savefig([glac(1:end-7) 'att_SL.fig'])
    end
end

% percent chances
pres_slp = zeros(length(glacs_out),1);
past_slp = zeros(length(glacs_out),1);

% getitng mean of all 1) modeled mb/sl years, 2) param suites, and 3) added
% points from inherent model error histograms
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

% boostrapping to get the 5-95th confidence levels
po = zeros(3,6,length(glacs_out));
if sl == 0
    for i=1:length(glacs_out)
        po(:,:,i) = bootstrp_mb_fn(pastallhs(:,i),presallhs(:,i),glac_11_mb,glac_18_mb,10000,0);
    end
else
    for i=1:length(glacs_out)
        po(:,:,i) = bootstrp_mb_fn(pastallhs(:,i),presallhs(:,i),glac_11_sl,glac_18_sl,10000,1);
    end
end

pmean = [mean(po(1,1,:)),mean(po(1,2,:)),mean(po(1,3,:)),mean(po(1,4,:)),nanmean(po(1,5,:)),mean(po(1,6,:))];
pmin = [min(po(2,1,:)),min(po(2,2,:)),min(po(2,3,:)),min(po(2,4,:)),min(po(2,5,:)),min(po(2,6,:))];
pmax = [max(po(3,1,:)),max(po(3,2,:)),max(po(3,3,:)),max(po(3,4,:)),max(po(3,5,:)),max(po(3,6,:))];
if svdat == 1
    fsv = ['att_out_' glac(1:end-2)]; 
    save(fsv,'meanprob','pmin','pmax'); 
end

clearvars -except glacrun y
end
