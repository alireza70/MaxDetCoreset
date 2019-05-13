function [M] = coreset_evaluation( A, k , parts, alpha, runs, sigma)
% This function compares the three methods for computing coresets: Greedy,
% Local Search and the LP-based algorithm.
% The input vectors are the rows of matrix A and the algorithm randomly
% partitions them into a few parts. The input parameter "parts" shows the
% number of parts.
% "runs" shows how many times to repeat each run of the algorithm, as the
% algrithm is randomized.
% sigma is the parameter to be used for the gaussian kernel.
% alpha is the parameter to be used for the LP-based approach.
% M will contain the a detailed result about the experiment: for each run,
% it keeps the value of the value of the determinant using different
% methods and their ratios along with the time it takes to compute them.
M = [k];
for run=1:runs
    n = size(A,1);
    for i=1:parts
        I{i} = []; %I{i} will contain the ith data set
    end
    for i=1:n % puts every point in a random part
        j = randi(parts);
        I{j}= [I{j},i];
    end
    CS_g=[];%union of greedy coresets
    CS_ls=[];%union of local search coresets
    CS_lp=[];%union of lp-based coresets

    gtime=0; % time to construct greedy coresets
    lstime = 0; % time to construct LS coresets
    lptime = 0; % time to construct LP-based coresets
    for i=1:parts % for every part we construct coresets using three methods
        % construct coresets using Greedy
        tic;
        [~,cs_temp_g] = greedy(A(I{i},:),k,sigma);
        gtime = gtime + toc;
        for j=1:k
            CS_g = [CS_g , I{i}(cs_temp_g(j))];
        end
        
        % construct coresets using LS
        tic;
        [~,cs_temp_ls] = LS(A(I{i},:),k,sigma);
        lstime = lstime + toc;
        for j=1:k
            CS_ls = [CS_ls , I{i}(cs_temp_ls(j))];
        end
        
        % construct coresets using the LP-based approach
        tic;
        cs_temp_lp = LP_Based(A(I{i},:),k, alpha, sigma);
        lptime = lptime + toc;
        for j=1:size(cs_temp_lp,2)
            CS_lp = [CS_lp , I{i}(cs_temp_lp(j))];
        end 
    end
    [ggc , ~] = greedy(A(CS_g,:),k, sigma); %GD/GD: greedy on union of greedy coresets

    [glsc , ~] = greedy(A(CS_ls,:),k, sigma); %GD/LS: greedy on union of local search coresets

    [glpc , ~] = greedy(A(CS_lp,:),k, sigma); %GD/LP: greedy on union of LP-based coresets
    
    [lslsc , ~] = LS(A(CS_ls,:),k, sigma); %LS/LS: LS on union of local search coresets

    disp('GD/GD GD/LS GD/LP LS/LS gtime lstime lptime');
    disp([ggc, glsc, glpc, lslsc, gtime/parts, lstime/parts, lptime/parts]);
    M = [M , run, ggc, glsc, glpc, lslsc, gtime/parts , lstime/parts, lptime/parts]; %generating the output
    M = [M , glsc/ggc , lslsc/ggc, glsc/glpc, lstime/gtime, lptime/lstime];
end
end

