function [a2n]=a2n(alpha,sigma,d)
    g=gam(alpha,sigma);
    a2n=exp(-i*g*d(2));
end