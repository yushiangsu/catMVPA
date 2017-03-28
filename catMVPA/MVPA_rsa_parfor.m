function [img_rsaz] = MVPA_rsa_parfor(imgsample, Searchlight_Mat, MC_corrvec)

%% Check Data structure
if ~isfield(imgsample.run{1}.cond{1}, 'cond_class')
    fprintf('%s\n','Can''t find cond_class! Make sure you enter the cond_class in your structure!');
    fprintf('%s\n','See "help MVPA_nfoldcrossval" for more information');
end
if ~isfield(imgsample.run{1}.cond{1}, 'img_vol_z')
    fprintf('%s\n','Can''t find img_vol_z! you have not z-score yet!');
end

%% Some information about imgsample
run_n = length(imgsample.run);
cond_n = length(imgsample.run{1}.cond);
sp_n = size(Searchlight_Mat,1);
sp_size =  size(Searchlight_Mat,2);
MC_n = length(MC_corrvec);

%% define tager matrix
tmp_mat = cell2mat(imgsample.run{1}.cond);
cond_vec = [tmp_mat(:).cond_class];
cond_class = unique(cond_vec);
dm_target = zeros(cond_n, cond_n);
for class_it = 1:length(cond_class)
    dm_target = dm_target + ...
        double(cond_vec == cond_class(class_it))' * double(cond_vec == cond_class(class_it));
end
dm_target = (dm_target*-1 + 1);

%% RSA
fprintf('%s', 'Calculating correlation for each Representational Dissimilarity Matrix ...');
rsa_z = nan(sp_n, 1);
% tic
% fprintf('%s', ' 00.00%');
parfor sp_iteration = 1:sp_n
    
    %% Get DM in sphere
    dm_rmat = zeros(cond_n, cond_n);
    %dm_pmat = zeros(cond_n, cond_n);
    for dm_x = 1:cond_n
        for dm_y = (dm_x+1):cond_n
            % Calculate coorelation for each pair (from dm_x and dm_y)
            pair_voxinx = nan(run_n, sp_size);
            pair_voxiny = nan(run_n, sp_size);
            
            for run_it = 1:run_n
                pair_voxinx(run_it,:) = ...
                    imgsample.run{run_it}.cond{dm_x}.img_vol_z(Searchlight_Mat(sp_iteration,:));
                pair_voxiny(run_it,:) = ...
                    imgsample.run{run_it}.cond{dm_y}.img_vol_z(Searchlight_Mat(sp_iteration,:));
            end
            [pair_r, pair_p] = corr(pair_voxinx(~isnan(pair_voxinx)), pair_voxiny(~isnan(pair_voxiny)));
            dm_rmat(dm_x, dm_y) = 1 - pair_r;
            dm_rmat(dm_y, dm_x) = 1 - pair_r;
        end
    end
    
    %% Correlation for that DM and get p-value
    sp_r = corr(squareform(dm_rmat,'tovector')',squareform(dm_target, 'tovector')');
    sp_p = sum(sp_r >= MC_corrvec) ./ MC_n;
    sp_z = (sp_p) -sqrt(2) * erfcinv(sp_p*2);
    
    rsa_z(sp_iteration) = sp_z;
    % fprintf('\b\b\b\b\b\b%5.2f%s', sp_iteration/sp_n*100, '%');
end
% toc
img_rsaz = nan(size(imgsample.run{1}.cond{1}.img_vol_z));
for sp_iteration = 1:sp_n
    img_rsaz(Searchlight_Mat(sp_iteration,1)) = rsa_z(sp_iteration);
end
fprintf(' %s\n', 'Done!');
