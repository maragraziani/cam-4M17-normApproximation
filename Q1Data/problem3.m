%% Basis pursuit denoising with L1 norm regularisation
% We want to solve 
% minimize f(x)+ lambda L1(x) 
% where L1 is the L1 norm of x
% and f(x) is the L2 norm of Ax-b

% Polishing:
% Use l1 heuristic to find x^ with the required sparsity pattern
% fix the pattern
% rerun the convex optimisation problem so to obtain the final (heuristic) solution

%% Cleaning the workspace
    clear all
    % Allocating spacefor saving the output
    n = zeros(5,1); % dimension n of each problem
    timings = zeros(5,1); % running time for each problem 
    results = zeros(5,2); % evaluation of the objective function and # iterations
    
%% Creating global variables 
% Note: in this way we can compute the function without passing arguments
    global A b 
    global eps
    global Gamma beta alpha
    
%% BATCH PROCESSING ALL PROBLEMS
%for problem_number=1:5 
problem_number=3
    %% Data loading
    if problem_number==1
    load('A1.mat')
    load('b1.mat')
    A=A1;
    b=b1;
    end
    if problem_number==2
    load('A2.mat')
    load('b2.mat')
    A=A2;
    b=b2;
    end
    if problem_number==3
    load('A3.mat')
    load('b3.mat')
    A=A3;
    b=b3;
    end
    if problem_number==4
    load('A4.mat')
    load('b4.mat')
    A=A4;
    b=b4;
    end
    if problem_number==5
    load('A5.mat')
    load('b5.mat') 
    A=A5;
    b=b5;
    end
   
    dim = size(A); % computing the size of the matrix A
    eps=0.001;
    alpha=0.49;
    beta=0.9;
%% Tuning gamma
Gamma = 0.1;
cardinality=64
while cardinality>8
    Gamma=Gamma+0.001
    % Performing gradient descent with Armjio backtracking linesearch
    x_start=randn(dim(2),1); %x_start
    %x_start=zeros(dim(2),1);
    stopping_crit=1e-5;
    tic
    [x, iterations]=backtracking_gradient_descent(x_start, @func3, @grad3 ,stopping_crit);
    time=toc
    disp('Done')
    f_value=func3(x)
%%
     min(abs(x))
    x(find(abs(x) < eps))=0;
     %%
    cardinality=card(x)
end
%% Finding the sparsity pattern
    sparsity_pattern=find(abs(x)>0);

%% Minimisation of the 2-norm
   
    A_small=A(:,sparsity_pattern);
    tic
    x_norm2=linsolve(A_small,b);
    time_2=toc;
    x_star=zeros(dim(2),1);
    x_star(sparsity_pattern)=x_norm2;
    value_2= sqrt(sum(power((A*x_star-b),2))) %2-norm computation
    sparsity_pattern
    x_norm2
    
%% Optional comparison with LASSO

    [B,FitInfo] = lasso(A,b);
    index = find(FitInfo.DF == 8);
    lambda = FitInfo.Lambda(index(1))
    x_lasso=B(:,index(1));
    sparsity_pattern_lasso=find(x_lasso~=0)
    x_nonzero=x_lasso(sparsity_pattern_lasso)
    value_2_lasso= sqrt(sum(power((A*x_lasso-b),2))) %2-norm computation
    
%% Optional comparison with L2 heuristic

n = 64;       % signal length
% reconstruct signal using l1 regularization
cvx_begin
variable x1(n);
minimize(norm(b-A*x1,2)+gamma*norm(x1,1) )
cvx_end
% reconstruct signal using l2 regularization
cvx_begin
variable x2(n);
minimize(norm(b-A*x2,2)+gamma*norm(x2,2) )
cvx_end
% plot signals
figure; subplot(3,1,1);
bar(b); axis([1 n -1.5 1.5]);
subplot(3,1,2);
bar(abs(A*x1)); axis([1 n -1.5 1.5]);
subplot(3,1,3);
bar(abs(A*x2)); axis([1 n -1.5 1.5]);

    