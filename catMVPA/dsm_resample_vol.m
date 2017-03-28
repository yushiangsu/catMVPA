function dsmimg_vol = dsm_resample_vol(imgsample, dsm_glm_b, Searchlight_Mat)

img_size = size(imgsample(1).img_vol);
glmX_p = size(dsm_glm_b, 2) - 1;
sp_n = size(Searchlight_Mat, 1);

dsmimg_vol = cell(glmX_p,1);
for glmX_it = 1:glmX_p
    dsmimg_vol{glmX_it} = nan(img_size);
end

for sp_it = 1:sp_n
    [idx_x, idx_y, idx_z] = ind2sub(img_size, Searchlight_Mat(sp_it,1));
    for glmX_it = 1:glmX_p
        dsmimg_vol{glmX_it}(idx_x, idx_y, idx_z) = dsm_glm_b(sp_it, glmX_it);
    end
end