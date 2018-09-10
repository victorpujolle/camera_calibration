function [val, idx] = argmin_J(mat)
%ARGMAX 이 함수의 요약 설명 위치
%   자세한 설명 위치
    [M, I] = min(mat(:));
    size_mat = size(mat);
    dim_mat = length(size_mat);
    dim_size = 1;
    idx = int64(zeros([1,dim_mat]));
    for i=1:dim_mat
        idx(i) = mod(floor((I-1)/dim_size),size_mat(i))+1;
        dim_size = dim_size*size_mat(i);
    end
    val = M;
end