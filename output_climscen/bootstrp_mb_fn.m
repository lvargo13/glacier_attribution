function [ po ] = bootstrp_mb_fn( past, pres, glac_11_mb, glac_18_mb, num_bs,sl_flg )
%  bootstrp_fn
% goal: do the same as bootstrp_sl_prob_loop but as a function to run
% within PDDout_loop_hist.m

% input: 

% --------- bootstrapping
like11 = zeros(num_bs,1);
like18 = zeros(num_bs,1);
past_11 = zeros(num_bs,1);
pres_11 = zeros(num_bs,1);
past_18 = zeros(num_bs,1);
pres_18 = zeros(num_bs,1);

le_past = round(length(past)*0.5);  % get 50% length
le_pres = round(length(pres)*0.5);

% this is what get repeated 
% need to get random different 50% each time, why this is in own loop
parfor i = 1:num_bs
  r = round( (length(past)-1).*rand(le_past,1) + 1);  % generate indicies for sampling
  [~,boot_ind]=bootstrp(1,[],r);
  tmp = r(boot_ind);
  out_past = past(tmp);
  
  r1 = round( (length(pres)-1).*rand(le_pres,1) + 1);  % generate indicies for sampling
  [~,boot_ind]=bootstrp(1,[],r1);
  tmp = r1(boot_ind);
  out_pres = pres(tmp);  

  if sl_flg == 0 
    past_11(i) = length(find(out_past<=glac_11_mb)) ./ length(out_past);
    pres_11(i) = length(find(out_pres<=glac_11_mb)) ./ length(out_pres);
    past_18(i) = length(find(out_past<=glac_18_mb)) ./ length(out_past);
    pres_18(i) = length(find(out_pres<=glac_18_mb)) ./ length(out_pres);
  else 
    past_11(i) = length(find(out_past>=glac_11_mb)) ./ length(out_past);
    pres_11(i) = length(find(out_pres>=glac_11_mb)) ./ length(out_pres);
    past_18(i) = length(find(out_past>=glac_18_mb)) ./ length(out_past);
    pres_18(i) = length(find(out_pres>=glac_18_mb)) ./ length(out_pres);   
  end
  
   like11(i) = (pres_11(i) ./ past_11(i)); 
   like18(i) = (pres_18(i) ./ past_18(i));
end

% % --------- plot
% figure; histogram(like11,20); xlabel('Likelihood'); ylabel('Probability'); title('2011'); hold on
% plot([mean(like11), mean(like11)], [0,150],'--k','LineWidth',2);
% std11 = 2*std(like11);
% plot([mean(like11)-std11, mean(like11)-std11], [0,150],'--k');
% plot([mean(like11+std11), mean(like11)+std11], [0,150],'--k');
% 
% figure; histogram(like18,20); xlabel('Likelihood'); ylabel('Probability'); title('2018'); hold on
% plot([mean(like18), mean(like18)], [0,160],'--k','LineWidth',2);
% std18 = 2*std(like18);
% plot([mean(like18)-std18, mean(like18)-std18], [0,160],'--k');
% plot([mean(like18+std18), mean(like18)+std18], [0,160],'--k');

% --------- uncertainties
varall = {past_11, pres_11, past_18, pres_18, like11, like18}; 
po = zeros(3, length(varall)); 

% calculate uncertainies using rank
% this po gives absolute min and max, not difference from mean
for i = 1:length(varall)
    tp = varall{i};
    po(1,i) = mean(tp);
    tmp = varall{i}; 
    tmp(isnan(tmp)) = 0;
    st = sort(tmp); 
    po(2,i) = st(num_bs*0.05);
    po(3,i) = st(num_bs - (num_bs*0.05));
end
po = [(po(:,1:4).*100) po(:,5:6)]; 

end