% edited from bootstrp_sl_prob for parameter suite

% 1) read in all gcm past and pres
% 2) use bootstrapping to get uncertainties
% 3) print and save(if wanted) mean, and 95% confidence level

clear

glac = 'brewster_p*/';
sl=0;  % flag, 0 for mb, 1 for sl
plotflag=1;  % 1 to plot
savedat=0;  % 1 to save 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch glac
    case 'rolleston_p*/'
        glac_11_sl = 1857 ;
        glac_18_sl = 1847 ; 
        glac_11_mb = -2039;
        glac_18_mb = -2456;
    case 'brewster_p*/'
        glac_11_mbsig = 187;
        glac_18_mbsig = 323;
        glac_11_sl = 2303 ; 
        glac_18_sl = 2311 ;  
        glac_11_mb = -1728; 
        glac_18_mb = -2217+glac_18_mbsig;
     case 'glenmary_p*/'
        glac_11_sl = 2297 ;
        glac_18_sl = 2327 ;   
     case 'ridge_p*/'
        glac_11_sl = 2436 ;
        glac_18_sl = 2434 ;   
     case 'thurneyson_p*/'
        glac_18_sl = 2232 ; 
        glac_11_sl = 2125 ;          
     case 'southcameron_p*/'
        glac_18_sl = 2426 ; 
        glac_11_sl = 2395 ;    
     case 'salisbury_p*/'
        glac_18_sl = 1967 ; 
        glac_11_sl = 1967 ; % no 2011 SL, so just to run
     case 'chancellor/'
        glac_18_sl = 1834 ; % both over top
        glac_11_sl = 1834 ;         
     case 'vertebrae12_p*/'
        glac_18_sl = 2086 ; 
        glac_11_sl = 2094 ; 
     case 'vertebrae25/'
        glac_18_sl = 2001 ; 
        glac_11_sl = 2001 ; 
     case 'parkpass_p*/'
        glac_18_sl = 1966 ; 
        glac_11_sl = 2040 ; 
end
        
% --------- read in data 
% gcm 
in_dir_gcm = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/gcm/';
fi_g = dir([in_dir_gcm glac]);
g_pres = zeros(320,length(fi_g));
g_past = zeros(1664,length(fi_g));
g_pres_sl = zeros(320,length(fi_g));
g_past_sl = zeros(1664,length(fi_g));
for ii = 1:length(fi_g)
    files = dir([in_dir_gcm '/' fi_g(ii).name '/','*.mat']);
    [g_pres(:,ii),g_past(:,ii),g_pres_sl(:,ii),g_past_sl(:,ii)] = read_pdd_gcm([in_dir_gcm '/' fi_g(ii).name '/'],files); 
end

% cesm
in_dir_cesm = '/Volumes/arc_03/vargola/glacier_attribution/output_climscen/cesm/';
fi_c = dir([in_dir_cesm glac]);
c_pres = zeros(680,length(fi_c));
c_pres_sl = zeros(680,length(fi_c));
c_past = zeros(1799,length(fi_c));
c_past_sl = zeros(1799,length(fi_c));
for ii = 1:length(fi_g)
    cfiles = dir([in_dir_cesm '/' fi_c(ii).name '/' '*present.mat']);
    [c_pres(:,ii),c_pres_sl(:,ii)] = read_pdd_cesm([in_dir_cesm '/' fi_c(ii).name '/'],cfiles);
    cfilesNAT = dir([in_dir_cesm '/' fi_c(ii).name '/' '*NAT.mat']);
    load([in_dir_cesm fi_c(ii).name '/cesm_NAT.mat' ])
    c_past(:,ii) = mb; 
    c_past_sl(:,ii) = ela;
end

if sl == 0
    past = [c_past(:)' g_past(:)'];
    pres = [c_pres(:)' g_pres(:)'];
    past_11p = (length(find(past<=glac_11_mb)) / length(past)) *100;
    pres_11p = (length(find(pres<=glac_11_mb)) / length(pres)) *100;
    past_18p = (length(find(past<=glac_18_mb)) / length(past)) *100;
    pres_18p = (length(find(pres<=glac_18_mb)) / length(pres)) *100;
else
    past = [c_past_sl(:)' g_past_sl(:)'];
    pres = [c_pres_sl(:)' g_pres_sl(:)'];
    past_11p = (length(find(past>=glac_11_sl)) / length(past)) *100;
    pres_11p = (length(find(pres>=glac_11_sl)) / length(pres)) *100;
    past_18p = (length(find(past>=glac_18_sl)) / length(past)) *100;
    pres_18p = (length(find(pres>=glac_18_sl)) / length(pres)) *100;
end
meanlike11 = pres_11p/past_11p;
meanlike18 = pres_18p/past_18p;

% --------- bootstrapping
num_bs = 10000;
like11 = zeros(num_bs,1);
like18 = zeros(num_bs,1);
past_11 = zeros(num_bs,1);
pres_11 = zeros(num_bs,1);
past_18 = zeros(num_bs,1);
pres_18 = zeros(num_bs,1);

le_past = round(length(past)/2);  % get 50% length
le_pres = round(length(pres)/2);

% this is what get repeated 
% need to get random different 50% each time, why this is in own loop
for i = 1:num_bs
  r = round( (length(past)-1).*rand(le_past,1) + 1);  % generate indicies for sampling
  [~,boot_ind]=bootstrp(1,[],r);
  tmp = r(boot_ind);
  out_past = past(tmp);
  
  r1 = round( (length(pres)-1).*rand(le_pres,1) + 1);  % generate indicies for sampling
  [~,boot_ind]=bootstrp(1,[],r1);
  tmp = r1(boot_ind);
  out_pres = pres(tmp);  

  if sl == 0
      past_11(i) = length(find(out_past<=glac_11_mb)) ./ length(out_past);
      pres_11(i) = length(find(out_pres<=glac_11_mb)) ./ length(out_pres);
      past_18(i) = length(find(out_past<=glac_18_mb)) ./ length(out_past);
      pres_18(i) = length(find(out_pres<=glac_18_mb)) ./ length(out_pres);
  else
      past_11(i) = length(find(out_past>=glac_11_sl)) ./ length(out_past);
      pres_11(i) = length(find(out_pres>=glac_11_sl)) ./ length(out_pres);
      past_18(i) = length(find(out_past>=glac_18_sl)) ./ length(out_past);
      pres_18(i) = length(find(out_pres>=glac_18_sl)) ./ length(out_pres);
  end
  
  like11(i) = (pres_11(i) ./ past_11(i)); 
  like18(i) = (pres_18(i) ./ past_18(i));
end
% like11((like11)==Inf)=NaN;
% like11 = like11(~isnan(like11));
% like18((like18)==Inf)=NaN; 
% like18 = like18(~isnan(like18));

% --------- plot
if plotflag == 1
    figure; histogram(like11,20); xlabel('Likelihood'); ylabel('Probability'); title('2011'); hold on
    plot([mean(like11), mean(like11)], [0,150],'--k','LineWidth',2);
    std11 = 2*std(like11);
    plot([mean(like11)-std11, mean(like11)-std11], [0,150],'--k');
    plot([mean(like11+std11), mean(like11)+std11], [0,150],'--k');

    figure; histogram(like18,20); xlabel('Likelihood'); ylabel('Probability'); title('2018'); hold on
    plot([mean(like18), mean(like18)], [0,160],'--k','LineWidth',2);
    std18 = 2*std(like18);
    plot([mean(like18)-std18, mean(like18)-std18], [0,160],'--k');
    plot([mean(like18+std18), mean(like18)+std18], [0,160],'--k');
end

% --------- uncertainties
varall = {past_11, pres_11, past_18, pres_18, like11, like18}; 
po = zeros(3, length(varall)); 

% calculate uncertainies using rank
% this po gives absolute min and max, not difference from mean
for i = 1:length(varall)
    tp = varall{i};
    tp((tp)==Inf)=NaN;
    po(1,i) = nanmean(tp); % mean w/o nan or inf
    
    tmp = varall{i}; 
    tmp(isnan(tmp)) = 0;
    st = sort(tmp); 
    po(2,i) = st(250);
    po(3,i) = st(9750);
end
po = [(po(:,1:4).*100) po(:,5:6)]; 
disp('past11, pres11, past_18, pres_18, like11, like18')
po  % print in command line

% --------- save
if savedat == 1
    if sl == 0
        sv = [ glac(1:end-1) '_mb_prob'];
    else
        sv = [ glac(1:end-1) '_sl_prob']; 
    end
    save(sv, 'po');
end
