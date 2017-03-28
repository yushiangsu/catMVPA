function [corr_r, corr_Z] = MVPA_dsm_corrZ(dsm_mat, dsm_target, method)

if nargin < 4
    method = 'Spearman';
end
X_obn = size(dsm_mat,2);
X_prn = length(dsm_target);
sp_n = size(dsm_mat,1);
             
X = nan(X_obn, X_prn);
for X_it = 1:X_prn
    X(:, X_it) = squareform(dsm_target{X_it},'tovector');
end

fprintf('%s', 'Calculating dsm correlation for each searchlight...');

corr_r = nan(sp_n, X_prn);
corr_Z = nan(sp_n, X_prn);
fprintf('%s', ' 00.000%');
for sp_it = 1:sp_n
    
    % correlation of dissimilarity pattern and target matrix
    dsm_pat = dsm_mat(sp_it,:)';
    r = corr(dsm_pat, X, 'type', method);
    corr_r(sp_it, :) = r;
    
    % fisher r-to-z transformation
    corr_Z(sp_it, :)=.5.*log((1+r)./(1-r));
    
    fprintf('\b\b\b\b\b\b\b%06.3f%s', sp_it/sp_n*100, '%');
end
fprintf('%s\n', ' Done!');

