
function imgsample = get_image_demeanvol(Searchlight_Mat, imgsample)

% De-mean value for given images
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
allspvox_ind = allspvox_ind(~isnan(imgsample(1).img_vol(allspvox_ind)));

%% creat new structure for output
smp_size = size(imgsample(1).img_vol);
for smp_it = 1:length(imgsample)
    imgsample(smp_it).img_vol_demean = zeros(smp_size);
end

%% subtracting mean
% for each voxel, calculate mean across length of imgsample and 
% subtracting mean.
fprintf('%s', 'De-mean: 00.000%');
for allspvox_it = 1:length(allspvox_ind)
    vox_idx = allspvox_ind(allspvox_it);
    val_vec = nan(length(imgsample), 1);
    for smp_it = 1:length(imgsample)
        val_vec(smp_it) = imgsample(smp_it).img_vol(vox_idx);
    end
    dm_val_vec = val_vec - mean(val_vec);
    for smp_it = 1:length(imgsample)
        imgsample(smp_it).img_vol_demean(vox_idx) = dm_val_vec(smp_it);
    end
    fprintf('\b\b\b\b\b\b\b%06.3f%s', (allspvox_it/length(allspvox_ind))*100, '%');
end
fprintf('%s\n', ' Done!');