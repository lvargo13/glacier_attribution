function [ po ] = uncert_distrib( varall )
% uncert_distrib- calculate uncertainties in using distributions of data
%   Calculate uncertainties in data that's been bootstrapped, using the
%   distribution of data

% input: 1) varall: 
% varall = {past_11, pres_11, past_18, pres_18, like11, like18}; 
% output: po (uncertainties)

dist={'exponential','normal','lognormal','rayleigh'};
fn_want = {'','','','','',''};
pdcert = zeros(3,length(varall)); 

for j = 1:length(varall)
    try
        f = fitmethis(varall{j},'figure','off','output','off');

        % get matching distribution
        a = 1;
        for ii = 1:10
            a = ~any(strcmp(dist,f(ii).name));
            if a == 0
                fn_want{j} = f(ii).name;
                break
            end
        end

        pd = fitdist(varall{j},fn_want{j});
        try 
          pdcert(1,j) = pd.mu;
          tf = strcmp(fn_want{j},'lognormal');
           if tf == 1  % if lognormal, need different mean than calculated
               pdcert(1,j) = mean(varall{j});
           end
        catch
          pdcert(1,j) = pd.B;  
        end

        tf = strcmp(fn_want{j},'exponential');
        rf = strcmp(fn_want{j},'rayleigh');
        ci = paramci(pd,'Alpha',0.05);
        if tf+rf == 1 % if expon. need to calculate sigma from range
            try
                pdcert(2,j)= (pd.mu - ci(1))*2;  % min  (exp uses mu)
                pdcert(3,j)= (ci(2) - pd.mu)*2;  % max 
            catch
                pdcert(2,j)= (pd.B - ci(1))*2;  % min  (rayleigh uses B)
                pdcert(3,j)= (ci(2) - pd.B)*2;  % max
            end
        else
            pdcert(2:3,j) = ci(:,end)*2;  %col2 in normal, but exp only has 1col
            %pdcert(2:3,j) = ci(:,1);
        end
    catch
        pdcert(:,j) = NaN;
    end
end
po = [pdcert(:,1:4)*100 pdcert(:,5:6)];

end

