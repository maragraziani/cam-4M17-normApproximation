function [y, g] = func_opt(x)
   global A b eps Gamma

   y=norm((A*x-b),2)+Gamma*norm(x,2);
end
