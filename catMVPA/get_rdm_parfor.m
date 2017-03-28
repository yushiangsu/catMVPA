function dsm_mat = get_rdm_parfor(imgsample, Searchlight_Mat, corr_measure)

if nargin < 3
    fprintf('%s\n%s\n', 'You did not specify which correlation method should be used.', 'Spearman rank will be use as default.')
    corr_measure = 'Spearman';
else
    if isempty(find(strcmp(corr_measure, {'Pearson', 'Spearman'}), 1))
        fprintf('%s\n', 'Unable to recognize your correlation method, Spearman rank will be use as default.')
        corr_measure = 'Spearman';
    end
end

%% Check Data structure
if ~isfield(imgsample, 'img_vol_demean')
    fprintf('%s\n','Can''t find img_vol_demean!');
end

%% Some information about imgsample
smp_n = length(imgsample);
sp_n = size(Searchlight_Mat,1);
sp_size = size(Searchlight_Mat,2);
%% RSA
fprintf('%s', 'Calculating correlation ...');
dsm_mat = nan(sp_n, (smp_n*(smp_n-1)/2));
% tic
% fprintf('%s', ' 00.000%');

parfor sp_iteration = 1:sp_n
    
    %% Get DM in sphere
    sp_voxsmp = nan(sp_size, smp_n);
    for smp_it = 1:smp_n
        sp_voxsmp(:, smp_it) = imgsample(smp_it).img_vol_demean(Searchlight_Mat(sp_iteration,:));
    end
    
    % remove zero or NAN
    sp_voxsmp(sum(sp_voxsmp == 0 | isnan(sp_voxsmp), 2) == smp_n, :) = [];
    
    % correlation by column
    switch corr_measure
        case 'Pearson'
            rdm_rmat = 1 - corr(sp_voxsmp, 'type', 'Pearson');
        case 'Spearman'
            rdm_rmat = 1 - corr(sp_voxsmp, 'type', 'Spearman');
    end
    rdm_vec = squareform(rdm_rmat,'tovector');
    
    %% Correlation for that DM and get p-value
    dsm_mat(sp_iteration, :) = rdm_vec;
%     fprintf('\b\b\b\b\b\b\b%06.3f%s', sp_iteration/sp_n*100, '%');
end
% toc
fprintf(' %s\n', 'Done!');
