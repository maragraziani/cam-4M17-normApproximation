function [y, g] = func3(x)
   global A b eps Gamma

   y=norm((A*x-b),2)+Gamma*sum(sqrt(power(x,2)+eps^2));
end
