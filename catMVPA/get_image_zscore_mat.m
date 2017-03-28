
function image_Vol_z = get_image_zscore_mat(Searchlight_Mat, image_Vol)

% Z-scoring all value in the image
%
% [image_Vol_z] = get_image_zscore_mat(Searchlight_Mat, image_Vol)
%
% Use the Searchlight_Mat from get_searchlight_mat.m to get all possible
% voxel inside given image. And z-score all these voxels, return a matrix 
% with z-scored 3D image matrix.
%
% Note: This z-score method will exclude any NaN generated from SPM GLM
% estimation. (The NaN were produced because SPM generate any NaN during 
% pre-processing procedure). And therefore, this function will generate the
% same NaN in resulted matrix.
% 
% Created by Yu-Shiang Su (2016/08/29)

%% found all voxel which will include in searchlight with index
allspvox_ind = unique(Searchlight_Mat);

%% remove NaNs
allspvox_ind = allspvox_ind(~isnan(image_Vol(allspvox_ind)));

%% get the value from these voxels and z-score.
image_allspvox = image_Vol(allspvox_ind);

%% plot distribution, report mean and std.
% mean(image_allspvox)
% std(image_allspvox)
% histfit(image_allspvox);

%% z-score
image_allspvox_z = zscore(image_allspvox);

%% report distribution, mean and std again.
% mean(image_allspvox_z)
% std(image_allspvox_z)
% histfit(image_allspvox_z);

%% putting it back to matrix
image_Vol_z = nan(size(image_Vol));
for n_vox = 1:length(image_allspvox_z)
    [x_corr, y_corr, z_corr] = ind2sub(size(image_Vol), allspvox_ind(n_vox));
    image_Vol_z(x_corr, y_corr, z_corr) = image_allspvox_z(n_vox);
end