%% Gradient descent with Armijo backtracking linesearch
% See also backtracking_gradient_descent.m
% Pseudocode:
% 1. delta_x=-grad(f(x))
% 
% 2. until the stopping condition do repeat:
%    linesearch to compute best t:
%           a.t=1
%           b. until f(x+t * delta_x)<f(x)+alpha*t * grad(f)' * delta_x do
%              repeat:
%              t=beta * t
% 3. update at each it x=x+t delta_x

%% Preparing the workspace
   
    clear all
    close all
    
    % Initialisation
    n = zeros(5,1); % dimension n of each problem
    timings_smoothed = zeros(5,1); % running time for each problem 
    results = zeros(5,2); % evaluation of the objective function and # iterations
    
%% Creating global variables 
% Note: in this way we can compute the function without passing arguments

    global A b 
    global eps
    global alpha beta
    
%% BATCH PROCESSING ALL PROBLEMS

for problem_number=1:5
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
   
    %global and local constants
    eps=0.01;
    alpha=0.3;
    beta=0.1;
    stopping_crit=1e-5;
    
%% Performing gradient descent with Armjio backtracking linesearch

for i=1:5
    x_start=zeros(dim(2),1);
    tic
    [x, iterations]=backtracking_gradient_descent(x_start, @func, @grad, stopping_crit)
    time=toc
    disp('Done')
    f_value=func(x)

%% Saving
    n(problem_number)=dim(2);
    timings_smoothed(problem_number) = time;
    results(problem_number, :) = [f_value, iterations];
%     n(i)=dim(2);
%     timings(i) = time;
%     results(i, :) = [f_value, iterations];
%     plot(timings)
end
end
%% Plotting the running time of the algorithms against the dimension n of x
% and comparing results with those from problem1

problem1 %obtaining results from problem1
close all %cleaning graphs
figure
hold on 
plot(n, timings(:,1), 'LineWidth',2)
plot(n, timings_smoothed, 'LineWidth',2)

title('Running time against dimensionality');
hold off

difference_timings=timings(:,1)-timings_smoothed(:,1);

%% Analysis of convergence rate wrt eps

problem_number=1
A=A1;
b=b1;
dim=size(A)
timings_eps=zeros(4,1);
results_eps=zeros(4,2);
%eps_range=[1e-4, 1e-3, 1e-2, 1e-1];
eps_range=linspace(1e-4, 1e-1);
i=0;
for i=1:100
% Performing gradient descent with Armjio backtracking linesearch
    eps=eps_range(i)
    x_start=zeros(dim(2),1); %x_start
    stopping_crit=1e-5;
    beta=0.5;
    tic
    [x, iterations]=backtracking_gradient_descent(x_start, @func, @grad, stopping_crit)
    time=toc
    disp('Done')
    f_value=func(x)

% Saving
    timings_eps(i) = time;
    results_eps(i, :) = [f_value, iterations];
    i=i+1;
end

%% Plotting the convergence rate for different values of epsilon
figure
hold on
title('Convergence rate for different values of \epsilon');
plot(eps_range,timings_eps, 'LineWidth',1.5)
