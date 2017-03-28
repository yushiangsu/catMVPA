function [corr_z] = MVPA_dsm_MCcorr(dsm_mat, dsm_target, mcper_n, method)

if nargin < 4
    method = 'Pearson';
end
X_obn = size(dsm_mat,2);
X_prn = length(dsm_target);
sp_n = size(dsm_mat,1);
             
X = nan(X_obn, X_prn);
for X_it = 1:X_prn
    X(:, X_it) = squareform(dsm_target{X_it},'tovector');
end

fprintf('%s', 'Calculating and permutating ...');

corr_z = nan(sp_n, X_prn);
fprintf('%s', ' 00.000%');
for sp_it = 1:sp_n
    % correlation of dissimilarity pattern and target matrix
    dsm_pat = dsm_mat(sp_it,:)';
    corr_r = corr(dsm_pat, X, 'type', method);

    % monte carlo permutation
    X_per = nan(X_obn, X_prn, mcper_n);
    for mcper_it = 1: mcper_n
        X_per(:, :, mcper_it) = X(randperm(X_obn),:);
    end
    
    % estimate correlation p-value by monte carlo's results
    corr_p = nan(X_prn,1);
    for X_it = 1:X_prn
        mcper_corr = corr(dsm_pat, squeeze(X_per(:,X_it,:)));
        corr_p(X_it) = sum(corr_r(X_it) > mcper_corr)/mcper_n;
    end

    % transform p-value to z-value
    corr_z(sp_it, :) = (corr_p) -sqrt(2) * erfcinv(corr_p*2);

    fprintf('\b\b\b\b\b\b\b%06.3f%s', sp_it/sp_n*100, '%');
end
fprintf('%s\n', ' Done!');

