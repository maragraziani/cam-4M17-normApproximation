function [x,histout,costdata] = probl2(x0,f,tol,ro, maxit)
% steepest descent with backtracking Armijo rule: 
%The iteration process stops if Armijo condition is not satisfied or
% a specifed maximum number of backtrackings =maxbt is reached.
%
% Input: x0 = initial iterate
%        f = objective function,
%            the calling sequence for f should be
%            [fout,gout]=f(x) where fout=f(x) is a scalar
%              and gout = grad f(x) is a COLUMN vector
%        tol = termination criterion norm(grad) < tol
%              e.g., tol= = 1.e-6
%         ro =reduction factor in the backtracking algorithm
%        maxit = maximum iterations 
%        
% Output: x = solution
%         histout = iteration history   
%             Each row of histout is
%            [ f, norm(grad), # of step length reductions j, iteration count]
%         costdata = [num f, num grad, num hess] (for steep, num hess=0)
%
%
%
% bactracking parms
% 
c=1.d-4;
maxbt=10; % maximum of tries in  the backtracking algorithm
%

itc=1; %iteration counter
xc=x0; 
[fc,gc]=f(xc);
numf=1; numg=1; numh=0;
ithist=zeros(maxit,5);
ithist(1,1)=fc; ithist(1,2) = norm(gc); ithist(1,3)=0; ithist(1,4)=1;
ithist(1,5)=itc-1;
while(norm(gc) > tol & itc <= maxit)
%
%      the original step apha=1
%
	 
    alpha=1;
    xt=xc-alpha*gc; ft=f(xt);
    fgoal=fc-c*alpha*(gc'*gc);
        numf=numf+1;
	j=0; itc=itc+1;
        
%
%       backtrack line search
%
      
	while(ft > fgoal & j<=maxbt)
		j=j+1;
        alpha=ro*alpha;
		xt=xc-alpha*gc;
		ft=f(xt); numf = numf+1; 
        fgoal=fc-c*alpha*(gc'*gc);
    end
    %
    %signal that Armijo might not be OK after maxbt backtrackings
    %
    if(ft > fgoal) 
		disp(sprintf('Max BackTrackings=%d atteined in step k=%d Armijo cond. failed.',maxbt, itc-1))
        x=xc; histout=ithist(1:itc-1,:); costdata=[numf, numg, numh];
    return; end
    
	xc=xt; [fc,gc]=f(xc); numf=numf+1; numg=numg+1;
	ithist(itc,1)=fc; ithist(itc,2) = norm(gc); ithist(itc,3)=j; ithist(itc,4)=alpha;
	ithist(itc,5)=itc-1; 
end
x=xc; 
histout=ithist(1:itc,:); costdata=[numf, numg, numh];