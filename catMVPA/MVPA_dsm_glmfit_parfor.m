function [dsm_glm_b] = MVPA_dsm_glmfit_parfor(dsm_mat, dsm_target)

glmX_n = size(dsm_mat,2);
glmX_p = length(dsm_target);
sp_n = size(dsm_mat,1);
             
glmX = nan(glmX_n, glmX_p);
for glmX_it = 1:glmX_p
    glmX(:, glmX_it) = squareform(dsm_target{glmX_it},'tovector');
end

fprintf('%s', 'Calculating glm for each searchlight ...');
dsm_glm_b = nan(sp_n, (glmX_p + 1));

% fprintf('%s', ' 00.000%');
parfor sp_it = 1:sp_n
    dsm_glm_b(sp_it, :) = glmfit([glmX ones(glmX_n, 1)], ... % GLM X
        dsm_mat(sp_it,:), ... % GLM Y
        'normal',... % GLM distribution
        'constant', 'off'); % because we input constant manually
    % fprintf('\b\b\b\b\b\b\b%06.3f%s', sp_it/sp_n*100, '%');
end
fprintf('%s\n', ' Done!');

