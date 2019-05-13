function [ x , OptI] = LS(A, k , sigma)
% This function computes a coreset of size k over the set of rows of matrix
% A using the Local Search algoirthm.
% sigma is the input parameter to be used for gaussian kernel
% The function outputs OptI which contains indices of the rows of A that 
% are chosen in the core-set and x is their corresponding determinant.

n = size(A,1);
[x,OptI] = greedy(A,k,sigma); %start with the solution of Greedy
[OptB,C] = construct_rbf( A(OptI,:) , sigma , false);

ep = 1e-5; %value of epsilon
improved = 1; %shows if we have made any progress in each iteration
while improved;
    improved = 0;
    for i=1:n  % iterate over all the points to be swapped in
        for j=1:k  % iterate over all points in the core-set to be swapped out
            I = OptI; % I is the set of new indices after the swap
            I(j) = i;
            [B,C] = construct_rbf( A(I,:) , sigma , false);
            if (det(B) > det(OptB)+ep) % if the swap improves the determinant
                OptI = I;
                OptB = B;
                x = det(OptB);
                improved = 1;
                break;
            end 
        end
    end
end
end

