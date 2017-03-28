function Searchlight_Mat = get_searchlight_new(radius, mask_file, imgsample, resimg)

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

%%
if nargin < 4
    resimg_check = 0;
else
    resimg_check = 1;
end


%% Check 
fprintf('%s\n', 'Checking ...')
ref_file = imgsample(1).img_file;
ref_file_V = spm_vol(ref_file);
mask_file_V = spm_vol(mask_file);
if ref_file_V.dim == mask_file_V.dim
    fprintf('%s\n', 'Good! mask and data have identical space.')
    mask_Vol = spm_read_vols(mask_file_V);
else
    fprintf('%s\n', 'data images have different space from mask file, generated a temporary mask file.')
    [ref_path, ref_name, ref_ext] = fileparts(ref_file_V.fname);
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

%% Information
img_size = size(mask_Vol);
sl_t = 0.5;
fprintf('%s\n', 'Generating Searchlights ...')

%% find out good voxels and bad voxels
if resimg_check == 0
    vox_invalid = unique([find(mask_Vol < 0.5); ...
        find(isnan(sum(cat(4,imgsample.img_vol),4)))]);
    vox_valid = intersect(find(mask_Vol > 0.5), ...
        find(~isnan(sum(cat(4,imgsample.img_vol),4))));
    
else
    vox_invalid = unique([find(mask_Vol < 0.5); ...
        find(isnan(sum(cat(4,imgsample.img_vol),4))); ...
        find(sum(cat(4, resimg.mask_vol), 4) < 5)]);
    vox_valid = intersect(intersect(find(mask_Vol > 0.5), ...
        find(~isnan(sum(cat(4,imgsample.img_vol),4)))), ...
        find(sum(cat(4, resimg.mask_vol), 4) == 5));
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
sl_voxn = size(Sphere_Mat, 1);

%%
Searchlight_Mat = nan(length(vox_valid), sl_voxn);
for i = 1:length(vox_valid)
    [vox_subx, vox_suby, vox_subz] = ind2sub(img_size, vox_valid(i));
    sl_corr = repmat([vox_subx, vox_suby, vox_subz], sl_voxn, 1) ...
        + Sphere_Mat;
    
    % remove voxels which is outside bpundary
    sl_corr((sl_corr(:,1) < 1 | sl_corr(:,1) > img_size(1) | ...
        sl_corr(:,2) < 1 | sl_corr(:,2) > img_size(2) | ...
        sl_corr(:,3) < 1 | sl_corr(:,3) > img_size(3)), :) = NaN;
    
    % remove voxels which is nan in your sample
    sl_idx = sub2ind(size(mask_Vol), ...
        sl_corr(:,1), sl_corr(:,2), sl_corr(:,3));
    sl_idx(ismember(sl_idx, vox_invalid)) = NaN;
    
    % remove NaN (invalid voxels)
    sl_idx(isnan(sl_idx(:,1))) = [];
    
    % If there are too many invalid voxel (sl_t% of voxels), the searchlight
    % doesn't count. If not, zero padded and store the vector.
    if length(sl_idx)/sl_voxn < sl_t
        Searchlight_Mat(i, :) = NaN;
    else
        Searchlight_Mat(i, :) = [sl_idx; zeros(sl_voxn - length(sl_idx),1)];
    end
end

% remove invalid searchlight (too many NaN inside the searchlight)
Searchlight_Mat(isnan(Searchlight_Mat(:,1)), :) = [];

%
if ref_file_V.dim ~= mask_file_V.dim
    if exist(tempmask_file, 'file')
        delete(tempmask_file);
        if strcmp(ref_ext, 'img')
            delete([ref_path '/m' ref_name '.hdr'])
        end
    end
end

fprintf('%s\n', 'Done!');
fprintf('%d%s\n', size(Searchlight_Mat,1), ' voxels included in searchlight.')