%clear all
% run through a suite of parameters
% need to comment out the parameters being run through in sa_pdd_parameters

% parameter_tests/param_test.m will compare output with mb data 

global CONFIG

ddf = 0.5:0.06:1.7; %original
radf = 0.15:0.015:0.25 ;
% ta = -1.65:0.1:-0.85; 
% pa = 0.8:0.1:1.7 ;

%runs = length(ddf) * length(radf); 

test_dir = '/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test_postreview/brewster_SL40/';

for ind1 = 1:length(ddf)
    for ind2 = 1:length(radf)
                      
        CONFIG.DegreeDay.DDF = ddf(ind1) ;
        CONFIG.DegreeDay.RadiationFactor = radf(ind2);
        CONFIG.DegreeDay.TempOffset = -1.25 ;
        CONFIG.DegreeDay.PptnFactor = 1.3;
        
        run('sa_pdd_brewster.m')
        
        sv = strcat('run_',num2str(ind1),'_',num2str(ind2),'.mat'); 
        %save([test_dir sv],'ela','snowth','ela_rmse','mb','mmb','mmb_dateadj','rmse','CONFIG');
        save([test_dir sv],'ela','snowth','ela_rmse','CONFIG');
    
    end
end