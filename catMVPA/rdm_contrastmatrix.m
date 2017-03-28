function C = rdm_contrastmatrix(cond_n)

% generate matrix C with deminsion K*(K-1)/2 x K

C = zeros(cond_n*(cond_n-1)/2, cond_n);
counter = 0;
for i = 1:cond_n
    for j = i+1:cond_n
        counter = counter + 1;
        C(counter, i) = 1;
        C(counter, j) = -1;
    end
end