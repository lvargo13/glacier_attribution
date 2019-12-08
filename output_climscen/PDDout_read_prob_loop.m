% read in output files of sa_pdd_gcm_loop_param
% calculate and plot probability
% outputs: plots of natural and present probabilities 
% can also see probabilities without bootstrapping here

glac = 'rolleston_p*/';
mbplt = 0;  % 0 if no MB, 1 if yes MB
yrprob = 2; % 1 for 2011 probabilities; not 1 for 2018 prob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch glac
    case 'brewster_p*/'
        glac_11_sl = 2303 ; 
        glac_18_sl = 2311 ;  
        glac_11_mb = -1728 ; 
        glac_18_mb = -2217;
        glac_18_mbsig = 323;
        glac_11_mbsig = 187;
        sl_bins = 1600:50:2400;
    case 'rolleston_p*/'
        glac_11_sl = 1857 ;
        glac_18_sl = 1847 ; 
        glac_11_mb = -2039;
        glac_18_mb = -2456;
        sl_bins = 1725:20:1905;
    case 'parkpass_p*/'
        glac_18_sl = 1966 ; 
        glac_11_sl = 2040 ; 
        sl_bins = 1450:50:2250;
    case 'southcameron_p*/'
        glac_18_sl = 2426 ; 
        glac_11_sl = 2395 ; 
        sl_bins = 2000:35:2700; 
    case 'salisbury_p*/'
        glac_18_sl = 1967 ; 
        glac_11_sl = 1967 ; % no 2011 SL, so just to run
        sl_bins = 1300:50:2500;
    case 'thurneyson_p*/'
        glac_18_sl = 2232 ; 
        glac_11_sl = 2125 ; 
        sl_bins = 1750:30:2350;
    case 'ridge_p*/'
        glac_11_sl = 2436 ;
        glac_18_sl = 2434 ;
        sl_bins = 2050:40:2600;
    case 'vertebrae12_p*/'
        glac_18_sl = 2086 ; 
        glac_11_sl = 2094 ; 
        sl_bins = 1700:40:2200;
    case 'glenmary_p*/' 
        glac_11_sl = 2297 ;
        glac_18_sl = 2327 ; 
        sl_bins = 2020:30:2490; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---- read in data
a = 4;  % scenarios
mb_bins = -3450:300:3450;

in_dir_g = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/gcm/';
fi_g = dir([in_dir_g glac]);
in_dir_c = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/cesm/';
fi_c = dir([in_dir_c glac]);
glacs_out = zeros(length(fi_g),1); %how many param. combos were run

mbp = zeros(length(mb_bins),a,length(glacs_out)); % bins x 4 scenarios x number of pruns
slp = zeros(length(sl_bins),a,length(glacs_out)); % bins x 4 scenarios x number of pruns

g_pres_a = zeros(320,length(glacs_out)); 
g_past_a = zeros(1664,length(glacs_out)); 
c_pres_a = zeros(680,length(glacs_out)); 
c_past_a = zeros(1799,length(glacs_out));
g_pres_sla = zeros(320,length(glacs_out)); 
g_past_sla = zeros(1664,length(glacs_out)); 
c_pres_sla = zeros(680,length(glacs_out)); 
c_past_sla = zeros(1799,length(glacs_out));

mb_prob = zeros(length(mb_bins),a); 
sl_prob = zeros(length(sl_bins),a); 
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
c_past_sl = ela; %sl(:,2); 

climscen = {g_past, g_pres, c_past, c_pres};
% mass balance probability
for cs = 1:a  
  cs_var = climscen{cs};  
    for i = 1:length(cs_var)
          [~, ind] = min(abs(mb_bins-cs_var(i)));  % get ind of closest bin
          mb_prob(ind,cs) = mb_prob(ind,cs)+1;
    end
   mb_prob(:,cs) = mb_prob(:,cs) / length(cs_var); 
end

% snowline probability
climscens = {g_past_sl, g_pres_sl, c_past_sl, c_pres_sl};
for cl = 1:a   
  cl_var = climscens{cl};
  cl_var = cl_var(~isnan(cl_var));  % get rid of nans
    for i = 1:length(cl_var)
          [~, indl] = min(abs(sl_bins-cl_var(i)));  % get ind of closest bin
          sl_prob(indl,cl) = sl_prob(indl,cl)+1;
    end
   sl_prob(:,cl) = sl_prob(:,cl) / length(cl_var); 
end

mbp(:,:,ii) = mb_prob;
slp(:,:,ii) = sl_prob; 
g_pres_a(:,ii) = g_pres; 
g_past_a(:,ii) = g_past; 
c_pres_a(:,ii) = c_pres; 
c_past_a(:,ii) = c_past;
g_pres_sla(:,ii) = g_pres_sl; 
g_past_sla(:,ii) = g_past_sl; 
c_pres_sla(:,ii) = c_pres_sl; 
c_past_sla(:,ii) = c_past_sl;
end

%---- plot mass balance
if mbplt == 1
    bns = zeros(length(mb_bins),2, 3); % 2 past and pres
    for j = 1:length(mb_bins)
        bns(j,1,1) = min([squeeze(mbp(j,1,:)); squeeze(mbp(j,3,:))]); % past
        bns(j,1,2) = mean([squeeze(mbp(j,1,:)); squeeze(mbp(j,3,:))]); 
        bns(j,1,3) = max([squeeze(mbp(j,1,:)); squeeze(mbp(j,3,:))]); 
        bns(j,2,1) = min([squeeze(mbp(j,2,:)); squeeze(mbp(j,4,:))]); %pres
        bns(j,2,2) = mean([squeeze(mbp(j,2,:)); squeeze(mbp(j,4,:))]); 
        bns(j,2,3) = max([squeeze(mbp(j,2,:)); squeeze(mbp(j,4,:))]); 
    end

    figure;
    stairs(mb_bins,bns(:,1,2),'k','LineWidth',2); hold on
    stairs(mb_bins,bns(:,2,2),'r','LineWidth',2)
    legend('natural','present')
    stairs(mb_bins,bns(:,1,1),'k'); hold on
    stairs(mb_bins,bns(:,1,3),'k')
    stairs(mb_bins,bns(:,2,1),'r')
    stairs(mb_bins,bns(:,2,3),'r')
    for ii = 1:length(mb_bins)-1
        x2 = [mb_bins(ii) mb_bins(ii+1) mb_bins(ii+1) mb_bins(ii)];
        fill_past= [bns(ii,1,1) bns(ii,1,1) bns(ii,1,3) bns(ii,1,3)];
        fill_pres= [bns(ii,2,1) bns(ii,2,1) bns(ii,2,3) bns(ii,2,3)];
        patch(x2, fill_past, 'k','FaceAlpha',.3,'LineStyle','none'); hold on
        fill(x2, fill_pres, 'r','FaceAlpha',.3,'LineStyle','none')
    end

    plot([glac_11_mb, glac_11_mb], [0, 0.25],'--k');
    plot([glac_18_mb, glac_18_mb], [0, 0.25],'--k');
    xlabel('Mass Balance (mm w.e.)')
    ylabel('Probability')

    % percent chances 
    % (this is also done in bootstrp_sl_prob* including bootstrapping)
    c_pres_p = zeros(length(glacs_out),1); 
    c_past_p = zeros(length(glacs_out),1);
    g_pres_p = zeros(length(glacs_out),1);
    g_past_p = zeros(length(glacs_out),1);
    if yrprob == 1
        f = glac_11_mb;
    else
        f = glac_18_mb;
    end
    for ii = 1:length(glacs_out)
        c_pres_p(ii) = (length(find(c_pres_a(:,ii)<=f)) / length(c_pres_a)) *100;
        c_past_p(ii) = (length(find(c_past_a(:,ii)<=f)) / length(c_past_a)) *100;
        g_pres_p(ii) = (length(find(g_pres_a(:,ii)<=f)) / length(g_pres_a)) *100;
        g_past_p(ii) = (length(find(g_past_a(:,ii)<=f)) / length(g_past_a)) *100;  
    end
    pres = [g_pres_p; c_pres_p];
    past = [g_past_p; c_past_p];
    [min(g_past_p) mean(g_past_p) max(g_past_p)]
    [min(c_past_p) mean(c_past_p) max(c_past_p)]
    [min(g_pres_p) mean(g_pres_p) max(g_pres_p)]
    [min(c_pres_p) mean(c_pres_p) max(c_pres_p)]
    [min(g_pres_p)/max(g_past_p) mean(g_pres_p)/mean(g_past_p) max(g_pres_p)/min(g_past_p)]
    [min(c_pres_p)/max(c_past_p) mean(c_pres_p)/mean(c_past_p) max(c_pres_p)/min(c_past_p)]
end


%---- plot snowline probability
bs = zeros(length(sl_bins),2, 3); % 2 past and pres
for j = 1:length(sl_bins)
    bs(j,1,1) = min([squeeze(slp(j,1,:)); squeeze(slp(j,3,:))]); % past
    bs(j,1,2) = mean([squeeze(slp(j,1,:)); squeeze(slp(j,3,:))]); 
    bs(j,1,3) = max([squeeze(slp(j,1,:)); squeeze(slp(j,3,:))]); 
    bs(j,2,1) = min([squeeze(slp(j,2,:)); squeeze(slp(j,4,:))]); %pres
    bs(j,2,2) = mean([squeeze(slp(j,2,:)); squeeze(slp(j,4,:))]); 
    bs(j,2,3) = max([squeeze(slp(j,2,:)); squeeze(slp(j,4,:))]); 
end

figure;
stairs(sl_bins,bs(:,1,2),'k','LineWidth',2); hold on
stairs(sl_bins,bs(:,2,2),'r','LineWidth',2)
legend('natural','present')
stairs(sl_bins,bs(:,1,1),'k'); hold on
stairs(sl_bins,bs(:,1,3),'k')
stairs(sl_bins,bs(:,2,1),'r')
stairs(sl_bins,bs(:,2,3),'r')
for ii = 1:length(sl_bins)-1
    x3 = [sl_bins(ii) sl_bins(ii+1) sl_bins(ii+1) sl_bins(ii)];
    fill_past= [bs(ii,1,1) bs(ii,1,1) bs(ii,1,3) bs(ii,1,3)];
    fill_pres= [bs(ii,2,1) bs(ii,2,1) bs(ii,2,3) bs(ii,2,3)];
    patch(x3, fill_past, 'k','FaceAlpha',.3,'LineStyle','none'); hold on
    fill(x3, fill_pres, 'r','FaceAlpha',.3,'LineStyle','none')
end
plot([glac_18_sl, glac_18_sl], [0,0.19],'--k');
plot([glac_11_sl, glac_11_sl], [0,0.19],'--k');
%xlim([1620 2420])
%ylim([0 0.19])
xlabel('snowline elevation (m)')
ylabel('Probability')

% percent chances
% (this is also done in bootstrp_sl_prob* including bootstrapping)
c_pres_slp = zeros(length(glacs_out),1); 
c_past_slp = zeros(length(glacs_out),1);
g_pres_slp = zeros(length(glacs_out),1);
g_past_slp = zeros(length(glacs_out),1);
if yrprob == 1
    f = glac_11_sl;
else
    f = glac_18_sl;
end
for ii = 1:length(glacs_out)
    c_pres_slp(ii) = (length(find(c_pres_sla(:,ii)>=f)) / length(c_pres_sla)) *100;
    c_past_slp(ii) = (length(find(c_past_sla(:,ii)>=f)) / length(c_past_sla)) *100;
    g_pres_slp(ii) = (length(find(g_pres_sla(:,ii)>=f)) / length(g_pres_sla)) *100;
    g_past_slp(ii) = (length(find(g_past_sla(:,ii)>=f)) / length(g_past_sla)) *100;   
end
pressl = [g_pres_slp; c_pres_slp];
pastsl = [g_past_slp; c_past_slp];
[min(g_past_slp) mean(g_past_slp) max(g_past_slp)]
[min(c_past_slp) mean(c_past_slp) max(c_past_slp)]
[min(g_pres_slp) mean(g_pres_slp) max(g_pres_slp)]
[min(c_pres_slp) mean(c_pres_slp) max(c_pres_slp)]
[min(g_pres_slp)/max(g_past_slp) mean(g_pres_slp)/mean(g_past_slp) max(g_pres_slp)/min(g_past_slp)]
[min(c_pres_slp)/max(c_past_slp) mean(c_pres_slp)/mean(c_past_slp) max(c_pres_slp)/min(c_past_slp)]

 
