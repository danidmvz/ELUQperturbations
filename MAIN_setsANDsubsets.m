%% EMPIRICALLY MOMENTS OF SETS VS TOTAL SET ===============================

close all
clear 
clc

Mp       = 6;
ns       = 1e5;
b1       = betarnd(2,5,ns,1);
mean_b1  = mean(b1);
prime_b1 = b1-mean_b1;
b1b1     = mean(prime_b1.^2);
b1b1b1   = mean(prime_b1.^3);
b1b1b1b1 = mean(prime_b1.^4);
min_b1   = min(b1);
max_b1   = max(b1);
b1Lims   = linspace(min_b1,max_b1,Mp+1);
b1mins   = b1Lims(1:end-1);
b1maxs   = b1Lims(2:end);
for k=1:Mp
    cont = 1;
    for i=1:ns
        if b1(i) >= b1mins(k) && b1(i) <= b1maxs(k)
            B(cont) = b1(i);
            cont    = cont + 1;
        end
    end
    B1{k} = B; clear B
    w(k)  = length(B1{k})/ns;
    % moments
    mean_B1(k)  = mean(B1{k});
    prime_B1    = B1{k}-mean_B1(k); 
    B1B1(k)     = mean(prime_B1.^2);  
    B1B1B1(k)   = mean(prime_B1.^3);  
    B1B1B1B1(k) = mean(prime_B1.^4);  
    % plot
    plot(B1{k},'.'); hold on
end

disp('Mean')
group_mean_b1 = sum(w.*mean_B1);
[mean_b1 group_mean_b1]'

disp('Var')
group_b1b1 = sum(w.*B1B1) + sum(w.*(mean_B1-mean_b1).^2);
[b1b1 group_b1b1]'

disp('3rd central')
group_b1b1b1 = sum(w.*B1B1B1) + 3*sum(w.*(mean_B1-mean_b1).*B1B1) ...
    + sum(w.*(mean_B1-mean_b1).^3);
[b1b1b1 group_b1b1b1]'

disp('4rd central')
group_b1b1b1b1 = sum(w.*B1B1B1B1) + 4*sum(w.*(mean_B1-mean_b1).*B1B1B1) ...
    + 6*sum(w.*(mean_B1-mean_b1).^2.*B1B1) + sum(w.*(mean_B1-mean_b1).^4);
[b1b1b1b1 group_b1b1b1b1]'


