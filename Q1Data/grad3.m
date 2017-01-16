function g = grad3(x)
    global A b eps Gamma
    %g=ones(length(x))*Gamma;
    g=(A'*(A*x-b))/norm(A*x-b,2)+Gamma* x./sqrt(power(x,2)+eps^2);    
% t=1;
% u=ones(length(x));
% q1 = 1./(u+x); 
% q2 = 1./(u-x);
% d1 = (q1.^2+q2.^2)/t;   d2 = (q1.^2-q2.^2)/t;
% 
% 
% % calculate gradient
% 
% gradphi = [A'*(A*x-b*2)-(q1-q2)/t; Gamma*ones(length(x),1)-(q1+q2)/t];
% g=gradphi';
end