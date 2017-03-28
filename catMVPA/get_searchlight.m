function Searchlight_Mat = get_searchlight(radius, mask_file, ref_file)

% Get a matrix for searchlight analysis
%
% [Searchlight_Mat] = get_searchlight_mat(radius, mask_file, ref_file)
%
% Create a searchlight matrix that contain "numbers of voxels in the mask"
% x "numbers of voxels in sphere". Each row present the coordinates of
% voxel in searchlight. Specifically, the first column will be the center
% of the searchlgiht. So the corresponding function could easily index the
% voxels.
%
% Note: the mask should better have the same matrix size and voxel size 
% as your images for MVPA analysis. Otherwise, it required ref_file and 
% call spm_mask to get a new and temporary mask.
%
% Dependencies: SPM (spm_vol, spm_read_vols, spm_mask)
%
% Created by Yu-Shiang Su (2016/08/26)


%% Configuration
if nargin == 3
    ref_file_V = spm_vol(ref_file);
    mask_file_V = spm_vol(mask_file);
    if ref_file_V.dim == mask_file_V.dim
        fprintf('%s\n', 'mask and ref have identical space.')
        mask_Vol = spm_read_vols(mask_file_V);
    else
        fprintf('%s\n', 'Reference image us in different space from mask file, generated a temporary mask file.')
        [ref_path, ref_name, ref_ext] = fileparts(ref_file);
        spm_mask(mask_file, ref_file, 0.5);
        tempmask_file = [ref_path '/m' ref_name ref_ext];
        mask_Vol = spm_read_vols(spm_vol(tempmask_file));
        % transform the mask to [1, 0], not [any number, NaN]
        if any(isnan(mask_Vol))
            fprintf('%s\n', 'Found your temporary mask contain NaN. Transform to [1, 0] mask.');
            mask_Vol(~isnan(mask_Vol)) = 1;
            mask_Vol(isnan(mask_Vol)) = 0;
        end
    end
% elseif nargin == 2
%     % Load the mask as a matrix
%     mask_Vol = spm_read_vols(spm_vol(mask_file));
% elseif nargin == 1
%     mask_file = spm_select(1,'image','Select mask file',[],pwd);
%     mask_Vol = spm_read_vols(spm_vol(mask_file));
end


%% define sphere by radius and create a matrix present all sphere location
Sphere_Mat = [];
% the centre of the sphere (i.e. [0 0 0]) should be the first element.
for x = [0, -radius:1:-1 1:1:radius]
    for y = [0, -radius:1:-1 1:1:radius]
        for z = [0, -radius:1:-1 1:1:radius]
            if sqrt(x.^2 + y.^2 + z.^2) <= radius
                Sphere_Mat(end +1, :) = [x, y, z];
            end
        end
    end
end

%% Load the mask and create matrix 
% find the index for each voxel
% Sometime the mask contrain 1.000 rather than 1. To avoid confussing and
% some special case that the mask is generated from probability map. I use
% threshold = 0.5 to index the true voxels from masked voxels.
% Exclude zero voxel from ref_file, these voxel are outside of bounding
% box of original space (if you do realign, coregistrate, normalization)
ref_file_vol = spm_read_vols(ref_file_V);
mask_Idx = find(mask_Vol > 0.5 & ref_file_vol ~= 0);

% a for loop to get searchlight for each voxel
Searchlight_Mat = nan(length(mask_Idx), size(Sphere_Mat, 1));
fprintf('%s%d\n', 'Number of voxels inside the mask: ', length(mask_Idx));
fprintf('%s', 'Generating searchlight matrix: ...');
for i = 1:length(mask_Idx)
    [x_corr, y_corr, z_corr] = ind2sub(size(mask_Vol), mask_Idx(i));
    Searchlight_corr = repmat([x_corr y_corr z_corr], ...
        size(Sphere_Mat, 1), 1) + Sphere_Mat;
    
    % For each serachlight, if any voxel in searchlight is outside of
    % bounding box for given mask (the mask_file you input or temporary
    % mask from ref_file), the searchlight don't count and ignore. If all
    % voxels in searchlight is inside the bounding box, transform the
    % subscript to index and store in a matrix to reduce memory loading
    % for output matrix.
    % Please note that, for the searchlight, which have voxel is outside of
    % the mask, is still count!
    if find(Searchlight_corr(:,1) < 1 | Searchlight_corr(:,1) > size(mask_Vol, 1) | ...
            Searchlight_corr(:,2) < 1 | Searchlight_corr(:,2) > size(mask_Vol, 2) | ...
            Searchlight_corr(:,3) < 1 | Searchlight_corr(:,3) > size(mask_Vol, 3))
        Searchlight_Mat(i, :) = NaN;
    else
        sl_idx = sub2ind(size(mask_Vol), ...
            Searchlight_corr(:,1), Searchlight_corr(:,2), Searchlight_corr(:,3));
        % exclude the searchlight contain too many zero in ref_file (these
        % voxel are close to bounding box of original image space.)
        if sum(ref_file_vol(sl_idx) == 0) >= (size(Sphere_Mat, 1)/2)
            Searchlight_Mat(i, :) = NaN;
        else
            Searchlight_Mat(i, :) = sl_idx;
        end
    end
end

% remove NaN (the invalid voxel) from matrix
Searchlight_Mat(isnan(Searchlight_Mat(:,1)), :) = [];

if ref_file_V.dim ~= mask_file_V.dim
    if exist(tempmask_file, 'file')
        delete(tempmask_file);
        if strcmp(ref_ext, 'img')
            delete([pwd '/m' ref_name '.hdr'])
        end
    end
end

fprintf('%s\n', 'Done!');
fprintf('%d%s\n', size(Searchlight_Mat,1), ' voxels include in searchlight.')