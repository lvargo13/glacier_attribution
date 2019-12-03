clear 

% read in output files of sa_pdd_gcm, calculate and plot probability

glac = 'brewster_SL_p3/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch glac
    case 'rolleston/'
        glac_11_sl = 1857 ;
        glac_18_sl = 1847 ; 
        glac_11_mb = -2039;
        glac_18_mb = -2456;
        mb_bins_step = 300;
        sl_bins_step = 19;
    case 'brewster_SL_p3/'
        glac_11_sl = 2303 ; 
        glac_18_sl = 2311 ;  
        glac_11_mb = -1728 ; 
        glac_11_mbsig = 187;
        glac_18_mb = -2217;
        glac_18_mbsig = 323;
        mb_bins_step = 300;
        sl_bins_step = 50;
     case 'glenmary/'
        glac_11_sl = 2297 ;
        glac_18_sl = 2327 ; 
        mb_bins_step = 250;
        sl_bins_step = 20;    
     case 'stocking/'
        glac_18_sl = 2312 ; % sl over top
        glac_11_sl = 2130 ;
        mb_bins_step = 250;
        sl_bins_step = 35;   
     case 'ridge/'
        glac_11_sl = 2436 ;
        glac_18_sl = 2434 ;
        mb_bins_step = 250;
        sl_bins_step = 20;  
     case 'dainty/'
        glac_18_sl = 2139 ; % sl over top
        mb_bins_step = 250;
        sl_bins_step = 20;   
     case 'thurneyson/'
        glac_18_sl = 2232 ; 
        glac_11_sl = 2125 ; 
        mb_bins_step = 250;
        sl_bins_step = 40;         
     case 'southcameron/'
        glac_18_sl = 2426 ; 
        glac_11_sl = 2395 ; 
        mb_bins_step = 250;
        sl_bins_step = 20;    
     case 'salisbury/'
        glac_18_sl = 1967 ; 
        glac_11_sl = 1967 ; % no 2011 SL, so just to run
        mb_bins_step = 250;
        sl_bins_step = 50; 
     case 'chancellor/'
        glac_18_sl = 1834 ; % both over top
        glac_11_sl = 1834 ; 
        mb_bins_step = 250;
        sl_bins_step = 35;         
     case 'vertebrae12/'
        glac_18_sl = 2086 ; 
        glac_11_sl = 2094 ; 
        mb_bins_step = 250;
        sl_bins_step = 35; 
     case 'vertebrae25/'
        glac_18_sl = 2001 ; 
        glac_11_sl = 2001 ; 
        mb_bins_step = 250;
        sl_bins_step = 20; 
     case 'parkpass/'
        glac_18_sl = 1966 ; 
        glac_11_sl = 2040 ; 
        mb_bins_step = 250;
        sl_bins_step = 50;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in data
a = 4; 
% gcm 
in_files_gcm = ['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/gcm/' glac];  
files = dir([in_files_gcm,'*.mat']);
[g_pres,g_past,g_pres_sl,g_past_sl] = read_pdd_gcm(in_files_gcm,files); 

% cesm
in_files_cesm = ['/Volumes/arc_03/vargola/glacier_attribution/output_climscen/cesm/' glac]; 
cfiles = dir([in_files_cesm,'*present.mat']);  % dont want NAT run
[c_pres,c_pres_sl] = read_pdd_cesm(in_files_cesm,cfiles); 
cfilesNAT = dir([in_files_cesm,'*NAT.mat']);
load([in_files_cesm cfilesNAT.name])
c_past = mb; 
c_past_sl = ela; %sl(:,2); 


%----- mass balance probability
mxg = round(max(g_past),-2);
mxc = round(max(c_past),-2);
maxmb = max(mxg,mxc); 
mng = round(min(g_pres),-2); 
mnc = round(min(c_pres),-2); 
minmb = min(mng,mnc);  
mb_bins = minmb-mb_bins_step:mb_bins_step:maxmb+(mb_bins_step*2); 
mb_prob = zeros(length(mb_bins),a); 
climscen = {g_past, g_pres, c_past, c_pres};

for cs = 1:a
  cs_var = climscen{cs};  
for i = 1:length(cs_var)
      [~, ind] = min(abs(mb_bins-cs_var(i)));  % get ind of closest bin
      mb_prob(ind,cs) = mb_prob(ind,cs)+1;
end
   mb_prob(:,cs) = mb_prob(:,cs) / length(cs_var); 
end

% plot
figure; 
mb_bins = mb_bins -(mb_bins_step/2); 
%figure
stairs(mb_bins,mb_prob(:,1),'r','LineWidth',2); hold on 
stairs(mb_bins,mb_prob(:,2),'k','LineWidth',2)
stairs(mb_bins,mb_prob(:,3),'r:','LineWidth',2)
stairs(mb_bins,mb_prob(:,4),'k:','LineWidth',2)

switch glac
    case 'rolleston/'
        plot([-2039, -2039], [0, 0.23],'--k'); 
        plot([-2456, -2456], [0, 0.23],'--k'); 
        xlim([-3500 3500])
        ylim([0 0.23])
    case 'brewster_2018/'
        plot([glac_11_mb, glac_11_mb], [0, 0.23],'--k');
%         plot([glac_11_mb-glac_11_mbsig, glac_11_mb-glac_11_mbsig], [0, 0.23],'k');
%         plot([glac_11_mb+glac_11_mbsig, glac_11_mb+glac_11_mbsig], [0, 0.23],'k');
        plot([glac_18_mb, glac_18_mb], [0, 0.23],'--k');
%         plot([glac_18_mb-glac_18_mbsig, glac_18_mb-glac_18_mbsig], [0, 0.23],'k');
%         plot([glac_18_mb+glac_18_mbsig, glac_18_mb+glac_18_mbsig], [0, 0.23],'k');
        xlim([-3500 3500])
        ylim([0 0.23])
    otherwise
        disp('no measured data')
end
xlabel('Mass Balance (mm w.e.)')
ylabel('Probability')
legend('GCM NAT','GCM present','CESM NAT','CESM present')

% percent chances
v = exist('glac_18_mb');

if v ~= 0
    c_pres_p = (length(find(c_pres<=glac_11_mb)) / length(c_pres)) *100;
    c_past_p = (length(find(c_past<=glac_11_mb)) / length(c_past)) *100;
    g_pres_p = (length(find(g_pres<=glac_11_mb)) / length(g_pres)) *100;
    g_past_p = (length(find(g_past<=glac_11_mb)) / length(g_past)) *100;   
end

%----- snowline probability
% c_past_sl(c_past_sl<1515)=1515; % only for park pass
% g_past_sl(g_past_sl<1515)=1515;

maxsl = max([max(c_pres_sl),max(g_pres_sl),max(c_past_sl),max(g_past_sl)])+(sl_bins_step*3); 
minsl = min([min(c_pres_sl),min(g_pres_sl),min(c_past_sl),min(g_past_sl)])-(sl_bins_step*3); 
mb_bins = minsl:sl_bins_step:maxsl;
mb_prob = zeros(length(mb_bins),a); 
climscen = {g_past_sl, g_pres_sl, c_past_sl, c_pres_sl};

for cs = 1:a
  cs_var = climscen{cs};
  cs_var = cs_var(~isnan(cs_var));  % get rid of nans
for i = 1:length(cs_var)
      [~, ind] = min(abs(mb_bins-cs_var(i)));  % get ind of closest bin
      mb_prob(ind,cs) = mb_prob(ind,cs)+1;
end
   mb_prob(:,cs) = mb_prob(:,cs) / length(cs_var); 
end

% plot
mb_bins = mb_bins -(sl_bins_step/2); 

figure
stairs(mb_bins,mb_prob(:,1),'r','LineWidth',2); hold on 
stairs(mb_bins,mb_prob(:,2),'k','LineWidth',2)
stairs(mb_bins,mb_prob(:,3),'r:','LineWidth',2)
stairs(mb_bins,mb_prob(:,4),'k:','LineWidth',2)

switch glac
    case 'rolleston/'
        plot([glac_18_sl, glac_18_sl], [0, max(mb_prob(:))],'--k');
        plot([glac_11_sl, glac_11_sl], [0, max(mb_prob(:))],'--k');
        xlim([1720 1890])
        ylim([0 0.55])
    case 'brewster_SL_p3/'
        plot([glac_18_sl, glac_18_sl], [0,0.19],'--k');
        plot([glac_11_sl, glac_11_sl], [0,0.19],'--k');
        xlim([1620 2420])
        ylim([0 0.19])
    case 'ridge/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k'); 
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');
        xlim([2080 2500])
     case 'glenmary/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');
        xlim([2020 2400])       
     case 'stocking/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1570 2370])  
        ylim([0 0.21])
     case 'thurneyson/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1760 2330])  
        ylim([0 0.27])   
     case 'southcameron/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([2020 2685])  
        %ylim([0 0.21])        
     case 'dainty/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1850 2280])   
        ylim([0 0.56])  
     case 'essex/'
        %plot([2243, 2243], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1450 2500])     
      case 'salisbury/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1300 2350])        
      case 'chancellor/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1650 1860])  
        ylim([0 0.62])
      case 'vertebrae12/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1700 2150])  
        ylim([0 0.62])
      case 'vertebrae25/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1710 2040])  
        ylim([0 0.65])
      case 'parkpass/'
        plot([glac_18_sl, glac_18_sl], [0,max(mb_prob(:))],'--k');  %top elev
        plot([glac_11_sl, glac_11_sl], [0,max(mb_prob(:))],'--k');  %top elev
        xlim([1420 2250])  
        ylim([0 0.17])   
     otherwise
        plot([meas_sl, meas_sl], [0,max(mb_prob(:))],'--k');
        xlim([xmin xmax])
end
xlabel('snowline elevation (m)')
ylabel('Probability')
legend('GCM NAT','GCM present','CESM NAT','CESM present')

% --- snowline probibilities
c_pres_sl_11_p = (length(find(c_pres_sl>=glac_11_sl)) / length(c_pres_sl)) *100;
c_past_sl_11_p = (length(find(c_past_sl>=glac_11_sl)) / length(c_past_sl)) *100;
g_pres_sl_11_p = (length(find(g_pres_sl>=glac_11_sl)) / length(g_pres_sl)) *100;
g_past_sl_11_p = (length(find(g_past_sl>=glac_11_sl)) / length(g_past_sl)) *100;        
c_pres_sl_18_p = (length(find(c_pres_sl>=glac_18_sl)) / length(c_pres_sl)) *100;
c_past_sl_18_p = (length(find(c_past_sl>=glac_18_sl)) / length(c_past_sl)) *100;
g_pres_sl_18_p = (length(find(g_pres_sl>=glac_18_sl)) / length(g_pres_sl)) *100;
g_past_sl_18_p = (length(find(g_past_sl>=glac_18_sl)) / length(g_past_sl)) *100; 

sv = [ glac(1:end-1) '_sl_prob']; 
save(sv, 'c_pres_sl_11_p','c_past_sl_11_p','g_pres_sl_11_p','g_past_sl_11_p','c_pres_sl_18_p','c_past_sl_18_p','g_pres_sl_18_p','g_past_sl_18_p');

% ---- normalize --------------------

norm_bins = linspace(0,1,20);
glac_norm_prob = zeros(length(norm_bins),a); 

for cs = 1:4
  cs_var = climscen{cs};
  cs_var = cs_var(~isnan(cs_var));  % get rid of nans
  glac_norm = (cs_var - min(cs_var)) / ( max(cs_var) - min(cs_var) );
  for i = 1:length(cs_var)
      [~, ind] = min(abs(norm_bins-glac_norm(i)));  % get ind of closest bin
      glac_norm_prob(ind,cs) = glac_norm_prob(ind,cs)+1;
  end
   glac_norm_prob(:,cs) = glac_norm_prob(:,cs) / length(cs_var); 
end

sv = [ glac(1:end-1) '_norm']; 
save(sv, 'glac_norm_prob');
