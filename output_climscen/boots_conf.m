function [ b, pchance ] = boots_conf( in_data, n, measmb )
%boots_conf 
%   1) get bootstrap sample for input data
%   2) calculate confidance levels 
% inputs: in_data- vector of mass balance of snowline data
% n- number of bootstrap samples

% run ex. [bo,pc] = boots_conf(c_pres,100,glac_11_mb);


[~,bootsamp]=bootstrp(n,[],in_data);

boot_out  = zeros(length(in_data),n); 
for i = 1:n
    boot_out(:,i) = in_data(bootsamp(:,i)); 
end

b = reshape(boot_out,(length(in_data)*n),1);  

pchance = (length(find(b<=measmb)) / length(b)) *100;

end

