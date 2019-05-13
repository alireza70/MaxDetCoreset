

function [OptI] = LP_Based( data, k, alpha, sigma)
%Computes the core-sets by the LP-based method on point set data (after passes the data thorough the RBF kernel. 
%Parameters: 1. data: point set 2. k: number of points (k parameter in the k-dpp)3. alpha: the constant in the LP which is the approximation 
% of the selected core-set in the projected space. 4. sigma: standard deviation for the RBF kernel that filters the data.


[~, A]=construct_rbf(data, sigma, 1); %Constructing the kernel by passing the data through the RBF kernel

[x,OptI] = greedy(data,2*k,sigma); % The algorithm first picks 2k vector by greedy


n = size(A,1);
B = orth(A(OptI,: )'); % Then projects all the data into the space of the these selected vectors.
S = (B')*(A');
start = 1;

%Next we repeat this process for this projected vectors. At the beginning S=Greedy output. Every time we find a direction "dir" so that 
%for some vector v which is not selected, its projection v' has this property <dir,v'> > alpha*|<dir,s>| for every s \in S 
while 1;
  [index, found , dir] = find_direction(S , OptI, alpha,start); %calling the module to solve the LP for already selected vector OptI and the set of all projected vectors S.
  if found==0;
      break;
  end
  start = index;
  candidate_index = 1;
  % Then we add the vector with largest dot product in the direction of dir to "S".
  for i=2:n
      if abs(dot(dir,S(:,i))) > abs(dot(dir,S(:,candidate_index)))
          candidate_index = i;
      end
  end
  OptI = [OptI,candidate_index];
end
end


function [index, found , dir] = find_direction(S,C,alpha,start)
% The function which solves multiple LPs, at most one for each unselected vector to find  a possible bad direction.
% S is the set of all vectors projected onto the output of greedy and C is the currently selected vectors.
[dim,n] = size(S);
found = 0;
dir = zeros(dim);
index= start;

m = size(C,2); % m is the number of vectors already choosen

% In this loop, we first fix the vector v and solve the LP <v,dir> >= \alpha*<dir,u> for all u already selected indexed by C in the code.
for i=start:n
    b = zeros([2*m+2,1]);
    A = zeros([2*m+2,dim+1]);
    for j=1:m
        A(j,:) = [S(:,C(j))',-1];
        A(j+m,:) = [-S(:,C(j))',-1];
    end
	
% It is a LP with 2m+1 constraints and dim+1 variables. The first dim variables are gonna be coordinates of "dir". The last one is going to be a 
%number l so that |<dir,u>|\leq l for all selected vectors u, and <dir,v_i> >alpha*l for vector i which can be accessed by S(:,i). 
	A(2*m+1,:) = [-S(:,i)', alpha];
    A(2*m+2,:) = [-S(:,i)',0];
    b(2*m+2,1)= -1;
    
    options = optimset('Display','none');
    
	
    [x,fval,exitflag,output] = linprog([],A,b,[],[],[],[],options);
    
    if exitflag == 1
        found = 1;
        dir = x(1:dim,1);
        index=i ;
        return;
	%%% if found ==1: a bad direction is found and returns as "dir"
    end
end
end
