% run through a suite of parameters
% parameter_tests/param_test.m will compare output with mb data 

gla = 'vertebrae25' ;


%%%%%%%%%%%%%%%%%%%%%%
global CONFIG

ddf = 0.5:0.2:1.7;
radf = 0.13:0.02:0.225 ;
ta = -2:0.2:-0.5 ;
pa = 0.8:0.2:1.8; 

runs = length(ddf) * length(radf) * length(ta) * length(pa); 

test_dir = ['/Volumes/arc_03/vargola/glacier_attribution/glacier_model/degreeday/param_test_postreview/' gla 'SL/'];

for ind1 = 1:length(ddf)
    
    for ind2 = 1:length(radf)
        
        for ind3 = 1:length(ta)
            
            for ind4 = 1:length(pa)
                      
        CONFIG.DegreeDay.DDF = ddf(ind1) ;
        CONFIG.DegreeDay.RadiationFactor = radf(ind2);
        CONFIG.DegreeDay.TempOffset = ta(ind3);
        CONFIG.DegreeDay.PptnFactor = pa(ind4);
        
        run('sa_pdd.m')
        
        sv = strcat('run_',num2str(ind1),'_',num2str(ind2),'_',num2str(ind3),'_',num2str(ind4),'.mat'); 
        save([test_dir sv],'ela','snowth','mmb','CONFIG');
        %save([test_dir sv],'mmb','mmb_dateadj','rmse','CONFIG');
            
            end
        end
    end
end