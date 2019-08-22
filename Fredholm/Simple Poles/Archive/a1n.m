function [a1n]=a1n(alpha,sigma,d)
    g=gam(alpha,sigma);
    a1n=exp(-i*g*d(1));
end