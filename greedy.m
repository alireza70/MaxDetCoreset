function [ x , OptI] = greedy( A, k , sigma)
% This function computes a coreset of size k over the set of rows of matrix
% A using the Greedy algorithm
% sigma is the input parameter to be used for gaussian kernel
% The function outputs OptI which contains indices of the rows of A that 
% are chosen in the core-set and x is their corresponding determinant.

n = size(A,1);
OptI = [];
x = 0; 
for i=1:k
    d = -1; %value of optimal determinant after the ith iteration
    for j=1:n %checking all candidate vectors
        I = OptI; 
        I = [I , j]; %candidate indices if j were to be in the core-set
        [B,C] = construct_rbf( A(I,:) , sigma , false);
        if (det(B) > d) %if j is the best candidate so far
            d = det(B);
            ind = j; %ind contains the optimal index
        end
    end
    OptI = [OptI , ind]; %add the optimal index to core-set
    x = d;
end
end

