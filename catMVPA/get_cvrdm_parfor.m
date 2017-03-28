function cvdsm_mat = get_cvrdm_parfor(imgsample, resimg, Searchlight_Mat)

% Get representational dissimilarity matrix from image between conditions.

% [dsm_mat] = get_rdm(imgsample, Searchlight_Mat, corr_measure)
%
%
% Created by Yu-Shiang Su (2016/10/20)



%% Check Data structure


%% Some information about imgsample
sl_n = size(Searchlight_Mat, 1);
run_n = length(resimg);
cond_n = max([imgsample.cond]);
resimg_size = size(resimg(1).res_vol);

if length(imgsample) ~= run_n*cond_n
    fprintf('%s\n', 'the length of imgsample is not consistent to design!')
end

%% Contrast matrix
C = rdm_contrastmatrix(cond_n);

%% Get RDM
fprintf('%s', 'Calculating cv Mahalanobis dist ...');
cvdsm_mat = nan(sl_n, (cond_n*(cond_n-1)/2));
 tic
% fprintf('%s', ' 00.000%');
parfor sl_it = 1:sl_n

    % read out the index for this searchlight
    sl_idx = Searchlight_Mat(sl_it,:);
    sl_idx(sl_idx == 0) = [];
    sl_beta = nan(cond_n, length(sl_idx), run_n);
    
    % get error matrix and beta matrix
    [idxx, idxy, idxz] = ind2sub(resimg_size(1:3),sl_idx);
    Ereg = cell(run_n, 1);
    for run_it = 1:run_n
        E = nan(resimg_size(4), length(sl_idx));
        for vox_it = 1:length(sl_idx)
            E(:,vox_it) = resimg(run_it).res_vol(idxx(vox_it), idxy(vox_it), idxz(vox_it), :);
        end
        
        Ereg{run_it} = covdiag(E);
        %         [E_t, E_n]=size(E{run_it});
        %         E = E - repmat(mean(E{run_it}), E_t, 1);
        %         E_sample = (1/E_t).*(E{run_it}'*E{run_it});
        %
        %         E_prior=diag(diag(E_sample));
        %
        %         E_d = 1/E_n*norm(E_sample-E_prior,'fro')^2;
        %         E_y = E.^2;
        %         E_r2 = 1/E_n/E_t^2*sum(sum(E_y'*E_y))-1/E_n/E_t*sum(sum(E_sample.^2));
        %         E_shrinkage = max(0, min(1, E_r2/E_d));
        %         fprintf('%6.5f\n', E_shrinkage)
        %         Ereg{run_it} = E_shrinkage * E_prior + (1 - E_shrinkage) * E_sample;
        
        for cond_it = 1:cond_n
            counter = (run_it-1)*cond_n + cond_it;
            sl_beta(cond_it, :, run_it) = imgsample(counter).img_vol(sl_idx);
        end
    end
    
    % cross validate with any possible pair of run
    cvdsm = nan(cond_n*(cond_n-1)/2, run_n*(run_n-1)/2);
    counter = 0;
    for run_A = 1:(run_n-1)
        for run_B = (1+run_A):run_n
            counter = counter + 1;
            cvdsm(:, counter) = diag(C*sl_beta(:,:,run_A)*(Ereg{run_A}^-0.5)*(Ereg{run_B}^-0.5)*sl_beta(:,:,run_B)'*C');
        end
    end
    cvdsm = mean(cvdsm, 2)/length(sl_idx);
    
    %% Correlation for that DM and get p-value
    cvdsm_mat(sl_it, :) = cvdsm;
    % fprintf('\b\b\b\b\b\b\b%06.3f%s', sl_it/sl_n*100, '%');
end
 toc
fprintf(' %s\n', 'Done!');
