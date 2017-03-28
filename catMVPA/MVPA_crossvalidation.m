
%% Configuration
SAMPLE_DIR = '/Users/yu-shiang/Research/MVPA_Try/FIRSTLEVEL_MVPA_ca15Cond/';
mask_file = '/Users/yu-shiang/Research/MVPA_Try/PROCESSED_SPM_Task/BrainMask/BrainMask.nii';

run_name = {'DM1','DM2','DM3','DM4','DM5'};
cond_name = {'HH_H','HH_M','HH_L',...
    'MH_H','MH_M','MH_L',...
    'MM_H','MM_M','MM_L',...
    'ML_H','ML_M','ML_L',...
    'LL_H','LL_M','LL_L'};
cond_class = {5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1};

radius = 2;
svmcost = 1;

ref_file = [SAMPLE_DIR '/' run_name{1} '/con_' sprintf('%04d', 1) '.nii'];

%% Prepare data structure
imgsample = [];
for run_n = 1:length(run_name)
    for cond_n = 1:length(cond_name)
        imgsample.run{run_n}.cond{cond_n}.img_file = ...
            [SAMPLE_DIR '/' run_name{run_n} '/con_' sprintf('%04d', cond_n) '.nii'];
        imgsample.run{run_n}.cond{cond_n}.cond_name = ...
            cond_name{cond_n};
        imgsample.run{run_n}.cond{cond_n}.cond_class = ...
            cond_class;
        imgsample.run{run_n}.cond{cond_n}.img_vol = ...
            spm_read_vols(spm_vol(imgsample.run{run_n}.cond{cond_n}.img_file));
    end
end

%% Get searchlight_mat
Searchlight_Mat = get_searchlight_mat(radius, mask_file, ref_file);

%% Z-score sample
for run_n = 1:length(run_name)
    fprintf('%s%d%s', 'Z-score run ', run_n, ' :');
    for cond_n = 1:length(cond_name)
        fprintf('%s%d', ' cond', cond_n);
        imgsample.run{run_n}.cond{cond_n}.img_vol_z = ...
            get_image_zscore_mat(Searchlight_Mat, imgsample.run{run_n}.cond{cond_n}.img_vol);
    end
    fprintf('%s\n', ' Done!');
end

%% N-fold crossvalidation
img_corr = nan(size(imgsample.run{1}.cond{1}.img_vol));
img_acc = nan(size(imgsample.run{1}.cond{1}.img_vol));
fprintf('%s', 'SVM classifying: 00.00%');
format long
for sp_iteration = 1:size(Searchlight_Mat,1)
    fold_corr = nan(length(run_name),1);
    fold_acc = nan(length(run_name),1);
    for fold_n = 1:length(run_name)
        txt_train = fopen(['./MVPA_train', num2str(fold_n)], 'w+');
        txt_test = fopen(['./MVPA_test', num2str(fold_n)], 'w+');
        spfeature = [];
        for run_n = 1:length(run_name)
            if run_n == fold_n
                for cond_n = 1:length(cond_name)
                    fprintf(txt_test, '%d ', cond_class{cond_n});
                    sp_vol = imgsample.run{run_n}.cond{cond_n}.img_vol_z(Searchlight_Mat(sp_iteration,:));
                    for vox_n = 1:length(sp_vol)
                        if ~isnan(sp_vol(vox_n))
                            fprintf(txt_test,'%d%s%f ', vox_n, ':', sp_vol(vox_n));
                        end
                    end
                    fprintf(txt_test, '\n');
                end
            else
                for cond_n = 1:length(cond_name)
                    fprintf(txt_train, '%d ', cond_class{cond_n});
                    sp_vol = imgsample.run{run_n}.cond{cond_n}.img_vol_z(Searchlight_Mat(sp_iteration,:));
                    for vox_n = 1:length(sp_vol)
                        if ~isnan(sp_vol(vox_n))
                            fprintf(txt_train, '%d%s%f ', vox_n, ':', sp_vol(vox_n));
                        end
                    end
                    fprintf(txt_train, '\n');
                end
            end
        end
        fclose(txt_train); fclose(txt_test);
        
        [MVPA_label_train, MVPA_inst_train] = libsvmread(['./MVPA_train', num2str(fold_n)]);
        MVPA_model = svmtrain(MVPA_label_train, MVPA_inst_train, ['-c ' num2str(svmcost) ' -t 0 -q']);
        [MVPA_label_test, MVPA_inst_test] = libsvmread(['./MVPA_test', num2str(fold_n)]);
        [predict_label, MVPA_accuracy, dec_values] = svmpredict(MVPA_label_test, MVPA_inst_test, MVPA_model, '-q');
        fold_corr(fold_n) = corr(MVPA_label_test, predict_label);
        fold_acc(fold_n) = MVPA_accuracy(1);
        
    end
    img_corr(Searchlight_Mat(sp_iteration,1)) = mean(fold_corr);
    img_acc(Searchlight_Mat(sp_iteration,1)) = mean(fold_acc);
    fprintf('\b\b\b\b\b\b%6.3f%s', ((sp_iteration/size(Searchlight_Mat,1))*100), '%');
end
fprintf(' %s\n', 'Done!');

ref_V = spm_vol(ref_file);
ref_V.fname = './img_corr.nii';
spm_write_vol(ref_V, img_corr);
