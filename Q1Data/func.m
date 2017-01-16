function [y, g] = func(x)
    global A b eps

    y= sum(sqrt(power(A*x-b,2)+eps^2));
%    g=grad(x);
end
