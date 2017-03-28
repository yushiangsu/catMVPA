function [MC_corrvec] = MonteCarlo_ras_sim(MC_n, run_n, cond_class, sp_size)

%% Preparation
cond_n = length(cond_class);

%% define tager matrix
cond_class_elem = unique(cond_class);
dm_target = zeros(cond_n, cond_n);
for class_it = 1:length(cond_class_elem)
    dm_target = dm_target + ...
        double(cond_class == cond_class_elem(class_it))' * double(cond_class == cond_class_elem(class_it));
end
dm_target = dm_target - eye(cond_n);

%% simulation for Monte Carlo method
MC_corrvec = nan(MC_n, 1);
fprintf('%s','Start Simulation for Monte Carlo method: 00.00%');
for MC_it = 1: MC_n
    simdm = nan((cond_n*(cond_n -1)/2), 1);
    for simdm_it = 1: (cond_n*(cond_n -1)/2)
        simdm(simdm_it) = corr(rand(run_n*sp_size,1),rand(run_n*sp_size,1));
    end
    MC_corrvec(MC_it) = corr(simdm, squareform(dm_target, 'tovector')');
    fprintf('\b\b\b\b\b\b%5.2f%s', MC_it/MC_n*100, '%')
    
end
fprintf('%s\n', ' Done!');
MC_crit_p = prctile(MC_corrvec, 95);
fprintf('%s%6.4f\n', 'Critical r from Monte-Carlo simulation is: ', MC_crit_p);



