function [img_corr, img_acc] = MVPA_nfoldcrossval_parfor(imgsample, Searchlight_Mat, svmcost)
% N fole crossvalidation to z-scored image.
%
% [img_corr, img_acc] = MVPA_nfoldcrossval(imgsample, Searchlight_Mat, svmcost)
%
% This function will get z-scored image voxel matrix from "imgsample" and
% use the Searchlight_Mat to do SVM with linear kernal to classcify each 
% searchlight. The class of each condition should be specified in
% "imgsample". Then the SVM will do for run_n times for crossvalidation, 
% resulting a 3D matrix with mean correlation and accuracy for each voxel.
% A temporary text file will be generated to do svmtrain, it will be
% removed by the end of this function.
% 
% Note: make sure you have the structure with field like below:
% imgsample.run{run_n}.cond{cond_n}.img_vol_z
% imgsample.run{run_n}.cond{cond_n}.cond_class
% 
% Dependicies: libsvmread, svmtrain, svmpredict from "libsvm"
%
% Created by Yu-Shiang Su (2016/09/02)

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

%% N-fold crossvalidation
fprintf('%s', 'SVM classifying...  ');
sp_n = size(Searchlight_Mat,1);
svm_corr = nan(sp_n, 1);
svm_acc = nan(sp_n, 1);
parfor sp_iteration = 1:sp_n
    fold_corr = nan(run_n,1);
    fold_acc = nan(run_n,1);
    
    svm_lab = nan((run_n .* cond_n), 1);
    svm_ins = nan((run_n .* cond_n), size(Searchlight_Mat,2));
    svm_run = nan((run_n .* cond_n), 1);
    
    for run_it = 1:run_n
        for cond_it = 1:cond_n
            svm_lab(((run_it-1)*cond_n) + cond_it) = imgsample.run{run_it}.cond{cond_it}.cond_class;
            svm_ins(((run_it-1)*cond_n) + cond_it, :) = ...
                imgsample.run{run_it}.cond{cond_it}.img_vol_z(Searchlight_Mat(sp_iteration,:));
            svm_run(((run_it-1)*cond_n) + cond_it) = run_it;
        end
    end
    
    % Transform nan to zero, which will be squzeed out through following
    % sparse matrix transformation.
    svm_ins(isnan(svm_ins)) = 0;
    
    for fold_n = 1:length(imgsample.run)
        svm_lab_tr = svm_lab(svm_run ~= fold_n);
        svm_lab_te = svm_lab(svm_run == fold_n);
        svm_ins_tr = sparse(svm_ins(svm_run ~= fold_n,:));
        svm_ins_te = sparse(svm_ins(svm_run == fold_n,:));
        
        % libsvm
        svm_model = svmtrain(svm_lab_tr, svm_ins_tr, ['-c ' num2str(svmcost) ' -t 0 -q']);
        [svm_lab_pred, svm_acc_pred, ~] = svmpredict(svm_lab_te, svm_ins_te, svm_model, '-q');
        
        fold_corr(fold_n) = corr(svm_lab_te, svm_lab_pred);
        fold_acc(fold_n) = svm_acc_pred(1);
        
    end
    svm_corr(sp_iteration) = .5.*log((1+mean(fold_corr))./(1-mean(fold_corr)));
    svm_acc(sp_iteration) = mean(fold_acc); 
end

img_corr = nan(size(imgsample.run{1}.cond{1}.img_vol_z));
img_acc = nan(size(imgsample.run{1}.cond{1}.img_vol_z));
for sp_iteration = 1:sp_n
    img_corr(Searchlight_Mat(sp_iteration,1)) = svm_corr(sp_iteration);
    img_acc(Searchlight_Mat(sp_iteration,1)) = svm_acc(sp_iteration);
end
fprintf(' %s\n', 'Done!');
