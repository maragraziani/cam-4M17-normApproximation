function g = grad(x)
    global A b eps

    den=sqrt(power(A*x-b,2)+eps^2); 
    num=(A*x-b);
    g=A'* (num./den);
end
