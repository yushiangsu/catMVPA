function dsm_mat = get_rdm(imgsample, Searchlight_Mat, corr_measure)

% Get representational dissimilarity matrix from image between conditions.

% [dsm_mat] = get_rdm(imgsample, Searchlight_Mat, corr_measure)
%
%
% Created by Yu-Shiang Su (2016/10/20)


if nargin < 3
    fprintf('%s\n%s\n', 'You did not specify which distance method should be used.', 'Pearson correlation will be use as default.')
    corr_measure = 'Pearson';
else
    if isempty(find(strcmp(corr_measure, {'Pearson', 'Spearman', 'Euclidean'}), 1))
        fprintf('%s\n', 'Unable to recognize your distance method, Pearson will be use as default.')
        corr_measure = 'Pearson';
    end
end

%% Check Data structure
if ~isfield(imgsample, 'img_vol_demean')
    fprintf('%s\n','Can''t find img_vol_demean! Will use raw value (img_vol)');
    rawvalue = 1;
else
    rawvalue = 0;
end

%% Some information about imgsample
smp_n = length(imgsample);
sp_n = size(Searchlight_Mat,1);
sp_size = size(Searchlight_Mat,2);
%% RSA
fprintf('%s', 'Calculating correlation ...');
dsm_mat = nan(sp_n, (smp_n*(smp_n-1)/2));
% tic
fprintf('%s', ' 00.000%');
for sp_iteration = 1:sp_n
    
    %% Get DM in sphere
%     rdm_rmat = zeros(smp_n, smp_n);
    
    sl_idx = Searchlight_Mat(sp_iteration,:);
    sl_idx = sl_idx(sl_idx ~= 0);
    sp_voxsmp = nan(length(sl_idx), smp_n);    
    if rawvalue
        for smp_it = 1:smp_n
            sp_voxsmp(:, smp_it) = imgsample(smp_it).img_vol(sl_idx);
        end
    else
        for smp_it = 1:smp_n
            sp_voxsmp(:, smp_it) = imgsample(smp_it).img_vol_demean(sl_idx);
        end
    end
    
    % remove zero or NAN
    sp_voxsmp(sum(sp_voxsmp == 0 | isnan(sp_voxsmp), 2) == smp_n, :) = [];
    
    % measure distance between smp_n
    switch corr_measure
        case 'Pearson'
            rdm_mat = 1 - corr(sp_voxsmp, 'type', 'Pearson');
            rdm_mat = squareform(rdm_mat,'tovector');
        case 'Spearman'
            rdm_mat = 1 - corr(sp_voxsmp, 'type', 'Spearman');
            rdm_mat = squareform(rdm_mat,'tovector');
        case 'Euclidean'
            rdm_mat = pdist(sp_voxsmp', 'Euclidean');
            
    end
    %% Correlation for that DM and get p-value
    dsm_mat(sp_iteration, :) = rdm_mat;
    fprintf('\b\b\b\b\b\b\b%06.3f%s', sp_iteration/sp_n*100, '%');
end
% toc
fprintf(' %s\n', 'Done!');
