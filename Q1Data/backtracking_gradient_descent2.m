function [x, iterations] = backtracking_gradient_descent2(x0,f,stopping_criterion,beta)
alpha=0.1;
max_it=1000;
xc=x0; 
[fc,grad_x]=f(xc);
iterations=0;
while(norm(grad_x) > stopping_criterion && iterations<max_it)
    iterations=iterations+1;
    t=1;
    xt=xc-t*grad_x;
    ft=f(xt);
    fgoal=fc-alpha*t*(grad_x'*grad_x);
    it=0;
    
    t=backtr(t,ft,grad_x, @func3, stopping_criterion, 0.5)
    xt=xc-t*grad_x;
    ft=f(xt)
% 	while(ft > fgoal && it<max_it)
%         it=it+1;
%         t=beta*t;
%         xt=xc-t*grad_x;
% 		ft=f(xt);
%         fgoal=fc-alpha*t*(grad_x'*grad_x);
%         fgoal-ft;
%     end
	xc=xt;
    [fc,grad_x]=f(xc);
end
x=xc;
end

