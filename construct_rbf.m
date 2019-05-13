function [K, L] = construct_rbf(data,sigma,flag_decomp)
% This function constructs a gaussian kernel K for the input matrix data 
% the input vectors are rows of the matrix data
% and if flag_decompose is true, it outputs its Cholesky decomposition in L
n=size(data,1);
K=data*data'/sigma^2;
d=diag(K);
K=K-ones(n,1)*d'/2;
K=K-d*ones(1,n)/2;
K=exp(K);

if flag_decomp
    [L,p] = chol(K); %decompose_kernel(K);
    L = L';
else
    L = 0;
end
