function [an]=an(alpha,sigma,d)
    g=gam(alpha,sigma);
    an=exp(-i*g*d);
end