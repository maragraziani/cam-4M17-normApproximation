function [x, iterations] = backtracking_gradient_descent(x0,f,g,stopping_criterion)
global alpha beta Gamma
xc=x0; 
fc=f(xc);
grad_x=g(xc);
iterations=0;
while(norm(grad_x,2) > stopping_criterion)
   % n=norm(grad_x,2)
    iterations=iterations+1;
    %1. linesearch
    t=1;
    ft=f(xc-t*grad_x);
    fgoal=fc-alpha*t*(grad_x'*grad_x);
	while(ft > fgoal)
        t=beta*t;
		ft=f(xc-t*grad_x);
        %fgoal=fc-alpha*t*(g(xc)'*grad_x);
        fgoal=fc-alpha*t*(grad_x'*grad_x);
        go=ft-fgoal;
    end
	xc=xc-t*grad_x;
    fc=f(xc);
    grad_x=g(xc);
end
x=xc;

end

